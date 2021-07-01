/*
 *  warp3.c  -  registers two grayscale images, producing a map of
 *                the corresponding points in the two images
 *
 *   Copyright (c) 2006-2009 Pittsburgh Supercomputing Center
 *
 *  This program is distributed in the hope that it will    *
 *  be useful, but WITHOUT ANY WARRANTY; without even the   *
 *  implied warranty of MERCHANTABILITY or FITNESS FOR A    *
 *  PARTICULAR PURPOSE.  Neither Carnegie Mellon University *
 *  nor any of the authors assume any liability for         *
 *  damages, incidental or otherwise, caused by the         *
 *  installation or use of this software.                   *
 *                                                          *
 *  CLINICAL APPLICATIONS ARE NOT RECOMMENDED, AND THIS     *
 *  SOFTWARE HAS NOT BEEN EVALUATED BY THE UNITED STATES    *
 *  FDA FOR ANY CLINICAL USE.                               *
 *
 *  HISTORY
 *    2006  Written by Greg Hood (ghood@psc.edu)
 *    2009    Simplified move generator, and refined
 *             confidence value computation to take into
 *             account image and reference masks (ghood@psc.edu)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <unistd.h>
#include <limits.h>
#include <dirent.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <mpi.h>
#include <sched.h>

#define GRAPHICS 0

#if GRAPHICS
#include <GL/glut.h>
#include <pthread.h>
#endif

#include "imio.h"
#include "dt.h"
#include "par.h"

#define DEBUG_MOVES	0
#define MASKING		1
#define FOLDING		0

#define MAX_LEVELS	32       /* max image size (w or h) is 2^MAX_LEVELS */

#define IMAGE(i,w,ix,iy)		i[(iy)*((size_t) w) + (ix)]
#define MASK(m,mbpl,ix,iy)		(m[(iy)*((size_t) mbpl) + ((ix) >> 3)] & (0x80 >> ((ix) & 7)))
#define SETMASKBIT(m,mbpl,ix,iy)  	m[(iy)*((size_t) mbpl) + ((ix) >> 3)] |= 0x80 >> ((ix) & 7)
#define CLEARMASKBIT(m,mbpl,ix,iy)  	m[(iy)*((size_t) mbpl) + ((ix) >> 3)] &= ~(0x80 >> ((ix) & 7))
#define MAP(map,w,ix,iy)		map[(iy)*((size_t) w) + (ix)]
#define GETMAP(map,w,ix,iy,xv,yv,cv)	{ MapElement *e = &MAP(map,w,ix,iy); *(xv) = e->x; *(yv) = e->y; *(cv) = e->c; }
#define SETMAP(map,w,ix,iy,xv,yv,cv)	{ MapElement *e = &MAP(map,w,ix,iy); e->x = xv; e->y = yv; e->c = cv; }
#define LINE_LENGTH	255


typedef struct Context {
  /* NOTE: any new fields added to this struct should
     be also added to PackContext and UnpackContext */ 
  char type;
  char imageBasename[PATH_MAX];
  char maskBasename[PATH_MAX];
  char discontinuityBasename[PATH_MAX];
  int strictMasking;

  char cptsName[PATH_MAX];
  char inputMapName[PATH_MAX];
  char constrainingMapName[PATH_MAX];

  char outputMapBasename[PATH_MAX];
  char outputImagesName[PATH_MAX];
  char outputCorrelationName[PATH_MAX];

  char logBasename[PATH_MAX];

  int startLevel;
  int outputLevel;
  int minResolution;
  int depth;
  double distortion;
  double correspondence;
  double correspondenceThreshold;
  double constraining;
  double constrainingThreshold;
  double quality;
  double minOverlap;

  int correlationHalfWidth;
  int writeAllMaps;
  int update;
  int partial;
} Context;

typedef struct Pair {
  char *imageName[2];
  int imageMinX[2], imageMaxX[2], imageMinY[2], imageMaxY[2];
  char *pairName;
} Pair;

typedef struct Task {
  /* NOTE: any new fields added to this struct should
     be also added to PackTask and UnpackTask */ 
  Pair pair;
} Task;

typedef struct Result {
  Pair pair;
  int updated;
  double distortion;
  double correlation;
  double correspondence;
  double constraining;
  char *message;
} Result;

typedef struct CPoint {
  float ix, iy;
  float rx, ry;
  double energy;
  double newEnergy;
} CPoint;

/* GLOBAL VARIABLES FOR MASTER */
int nResults = 0;
Result* results = 0;
#define DIR_HASH_SIZE	8192
char *dirHash[DIR_HASH_SIZE];
char summaryName[PATH_MAX] = "";
int vis = 0;

/* GLOBAL VARIABLES FOR MASTER & WORKER */
Context c;
Task t;
Result r;
FILE *logFile = NULL;

/* GLOBAL VARIABLES FOR WORKER */
int nLevels;                         /* number of resolution levels;
					0 = finest resolution */
int imageWidth[2][MAX_LEVELS],       /* image/reference width/height at */
    imageHeight[2][MAX_LEVELS];      /*   each level */
int imageOffsetX[2][MAX_LEVELS],     /* the amount each image is offset from */
    imageOffsetY[2][MAX_LEVELS];     /*   the origin in coordinates at that
		                     /*   level */
float* images[2][MAX_LEVELS];        /* images to be warped (at various */
                                     /*    resolution levels) */
unsigned char *masks[2][MAX_LEVELS]; /* masks of images to be warped */
unsigned char *idisc[2][MAX_LEVELS]; /* image discontinuity masks */

int mapWidth[MAX_LEVELS],            /* map width/height at each level */
    mapHeight[MAX_LEVELS];
int mapOffsetX[MAX_LEVELS],          /* the amount each map is offset from */
    mapOffsetY[MAX_LEVELS];          /*   the origin in coordinates at that */
                                     /*   level */
MapElement* maps[MAX_LEVELS];        /* the maps at each level; the coordinate
                                     /*   values within these maps are */
                                     /*   expressed in absolute coordinates */
                                     /*   at that level */
size_t minMaskCount[MAX_LEVELS];

int inputMapLevel;
int inputMapFactor;
MapElement *inputMap;
int mpwi, mphi;
int offxi, offyi;

int constrainingMapLevel;
int constrainingMapFactor;
MapElement *constrainingMap;
int mpwc, mphc;
int offxc, offyc;

int nCpts = 0;
CPoint *cpts = 0;

int windowWidth = 1024;
int windowHeight = 1024;
int displayLevel = -1;

/* FORWARD DECLARATIONS */
void MasterTask (int argc, char **argv, char **envp);
void MasterResult ();
void WorkerContext ();
void WorkerTask ();
void PackContext ();
void UnpackContext ();
void PackTask ();
void UnpackTask ();
void PackResult ();
void UnpackResult ();
int Init();
void Compute(char *outputName);
void ComputeWarpedImage (float *warped, unsigned char *valid,
			 int w, int h,
			 int imgox, int imgoy,
			 float *image, unsigned char *mask,
			 int iw, int ih,
			 MapElement *map,
			 int mapFactor,
			 int mpw, int mph,
			 int mox, int moy);
void ComputeCorrelation (float *correlation,
			 float *a, float *b,
			 unsigned char *valid,
			 int w, int h,
			 int hw);
void WriteOutputMap (char *outputName, int level, MapElement *map,
		     int mpw, int mph, int mox, int moy);
int WriteOutputImage (char *outputName, int level, float *output,
		      int w, int h, float scale);
int Compare (const void *x, const void *y);
int SortBySlice (const void *x, const void *y);
int SortByEnergy (const void *x, const void *y);
int ParseRange (char *s, int *pos, int *minValue, int *maxValue);
int ParseValue (char *s, int *pos, int *value);
size_t CountBits (unsigned char *p, size_t n);
int CreateDirectories (char *fn);
void Error (char *fmt, ...);
void Log (char *fmt, ...);
double GetCurrentTime ();
void* mainGraphics (void *arg);
void myDisplay ();
void myKeyboard (unsigned char k, int mouseX, int mouseY);
void myMotion (int mouseX, int mouseY);
void myMouse (int button, int state, int mouseX, int mouseY);
void myReshape (int width, int height);
void myIdle ();

int
main (int argc, char **argv, char **envp)
{
#if 0
  // FIX -- TEMPORARY HACK
  cpu_set_t cpumask;
  sleep(1);
  sleep(1);
  sleep(1);
  sleep(1);
  sleep(1);
  sleep(1);
  sleep(1);
  sleep(1);
  sleep(1);
  sleep(1);
  sleep(1);
  sleep(1);
  sleep(1);
  sleep(1);
  CPU_ZERO(&cpumask);
  CPU_SET(0, &cpumask);
  CPU_SET(1, &cpumask);
  CPU_SET(2, &cpumask);
  CPU_SET(3, &cpumask);
  CPU_SET(4, &cpumask);
  CPU_SET(5, &cpumask);
  CPU_SET(6, &cpumask);
  CPU_SET(7, &cpumask);
  sched_setaffinity(0, sizeof(cpu_set_t), &cpumask);
#endif

#if GRAPHICS
  pthread_t mgt;
  pthread_attr_t attr;
  int i;

  for (i = 1; i < argc; ++i)
    if (strcmp(argv[i], "-vis") == 0)
      {
	glutInit(&argc, argv); 
	if (pthread_attr_init(&attr) != 0)
	  Error("pthread_attr_init failed\n");
	if (pthread_create(&mgt, &attr, mainGraphics, NULL) != 0)
	  Error("pthread_create failed\n");
	break;
      }
#endif

  t.pair.imageName[0] = NULL;
  t.pair.imageName[1] = NULL;
  t.pair.pairName = NULL;
  r.pair.imageName[0] = NULL;
  r.pair.imageName[1] = NULL;
  r.pair.pairName = NULL;
  r.message = NULL;

  par_process(argc, argv, envp,
              (void (*)()) MasterTask, MasterResult,
              WorkerContext, WorkerTask, NULL,
              PackContext, UnpackContext,
              PackTask, UnpackTask,
              PackResult, UnpackResult);
  return(0);
}

#if GRAPHICS
/* GLUT FUNCTIONS */

void *
mainGraphics (void *arg)
{
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB); 
  glutInitWindowSize(windowWidth, windowHeight); 
  glutInitWindowPosition(100, 100);
  glutCreateWindow("warp3");
  glutDisplayFunc(myDisplay);
  glutReshapeFunc(myReshape); 
  glutMouseFunc(myMouse); 
  glutMotionFunc(myMotion);
  glutKeyboardFunc(myKeyboard); 
  glutIdleFunc(myIdle);
  glClearColor(1.0,1.0,1.0,0.0); 
  glutMainLoop();
  return(NULL);
}

void myDisplay()
{
  glClear(GL_COLOR_BUFFER_BIT);
  glFlush();
    //  if (!displayMapValid)
    //    return;
  
}

void myKeyboard (unsigned char k, int mouseX, int mouseY) 
{
  switch (k) 
    {
    case 27:  /* Escape key to exit */
      exit(0);
      break;

    default: /* Otherwise keep going */
      return;
    }
}

void myMotion (int mouseX, int mouseY) 
{
  //  xo = mouseX / 1024.0;
  //  yo = mouseY / 1024.0;
  //  myDisplay();
}

void myMouse(int button, int state, int mouseX, int mouseY) 
{
  switch (button)
    {
    case GLUT_LEFT_BUTTON:
      if (state == GLUT_DOWN)
        {
	  //  xo = mouseX / 1024.0;
	  //  yo = mouseY / 1024.0;
          // myDisplay();
        }
      break;

    case GLUT_RIGHT_BUTTON:  /* Right button to exit */
      exit(0);
      break;

    default: /* Otherwise keep going */
      return;
    }
}

void myReshape(int width, int height) 
{
  glViewport(0, 0, (GLsizei) width, (GLsizei) height); 
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  gluOrtho2D(0.0, 1.0, 0.0, 1.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

void
myIdle ()
{
  int mpw, mph;
  MapElement *map;
  int x, y;
  float cx, cy;
  float displayFactor;
  float xv, yv, cv;
  int valid;

  if (displayLevel < 0)
    return;
    
  // display the current map
  map = maps[displayLevel];
  mpw = mapWidth[displayLevel];
  mph = mapHeight[displayLevel];
  glClear(GL_COLOR_BUFFER_BIT);
  displayFactor = 0.5 / imageWidth[1][displayLevel];
  if (0.5 / imageHeight[1][displayLevel] < displayFactor)
    displayFactor = 0.5 / imageHeight[1][displayLevel];
  glBegin(GL_LINE_STRIP);
  glColor3f(1.0, 0.0, 0.0);
  glVertex2f(0.25, 0.75);
  glColor3f(1.0, 0.0, 0.0);
  glVertex2f(0.25 + displayFactor * imageWidth[1][displayLevel], 0.75);
  glColor3f(0.0, 1.0, 0.0);
  glVertex2f(0.25 + displayFactor * imageWidth[1][displayLevel],
	     0.75 - displayFactor * imageHeight[1][displayLevel]);
  glColor3f(1.0, 0.0, 0.0);
  glVertex2f(0.25,
	     0.75 - displayFactor * imageHeight[1][displayLevel]);
  glColor3f(1.0, 0.0, 0.0);
  glVertex2f(0.25, 0.75);
  glEnd();
  for (y = 0; y < mph; ++y)
    {
      valid = 0;
      for (x = 0; x < mpw; ++x)
	{
	  GETMAP(map, mpw, x, y, &xv, &yv, &cv);
	  if (cv == 0.0)
	    {
	      if (valid)
		{
		  glEnd();
		  valid = 0;
		}
	      continue;
	    }
	  if (!valid)
	    {
	      glBegin(GL_LINE_STRIP);
	      valid = 1;
	    }
	  cx = 0.25 + displayFactor * (xv - imageOffsetX[1][displayLevel]);
	  cy = 0.75 - displayFactor * (yv - imageOffsetY[1][displayLevel]);
	  glColor3f(0.0, 0.0, 0.0);
	  glVertex2f(cx, cy);
	}
      if (valid)
	glEnd();
    }
  for (x = 0; x < mpw; ++x)
    {
      valid = 0;
      for (y = 0; y < mph; ++y)
	{
	  GETMAP(map, mpw, x, y, &xv, &yv, &cv);
	  if (cv == 0.0)
	    {
	      if (valid)
		{
		  glEnd();
		  valid = 0;
		}
	      continue;
	    }
	  if (!valid)
	    {
	      glBegin(GL_LINE_STRIP);
	      valid = 1;
	    }
	  cx = 0.25 + displayFactor * (xv - imageOffsetX[1][displayLevel]);
	  cy = 0.75 - displayFactor * (yv - imageOffsetY[1][displayLevel]);
	  glColor3f(0.0, 0.0, 0.0);
	  glVertex2f(cx, cy);
	}
      if (valid)
	glEnd();
    }
  glFlush();
  displayLevel = -1;
}
#endif

/* MASTER PROCEDURES */

void
MasterTask (int argc,
            char **argv,
            char **envp)
{
  int i;
  int n;
  int len;
  char pairsFile[PATH_MAX];
  int nPairs;
  Pair *pairs;
  int minSlice, maxSlice;
  int error;
  int pos;
  DIR *dir;
  struct dirent *de;
  char fn[PATH_MAX];
  FILE *f;
  char inputDirName[PATH_MAX];
  char inputPrefix[PATH_MAX];
  int inputPrefixLen;
  int minN, maxN;
  int digitLen;
  char maskDirName[PATH_MAX];
  char maskPrefix[PATH_MAX];
  int maskPrefixLen;
  int maskMinN, maskMaxN;
  char outputDirName[PATH_MAX];
  char tc;
  unsigned int w, h;
  int m;
  int xRes, yRes;
  int minS, maxS;
  int forward;
  struct stat statBuf;
  int pn;
  char imgn[2][PATH_MAX];
  char pairn[PATH_MAX];
  int imgMinX[2], imgMaxX[2], imgMinY[2], imgMaxY[2];
  unsigned int iw, ih;
  int imageSlice, refSlice;
  int imi;

  error = 0;
  c.type = '\0';
  c.imageBasename[0] = '\0';
  c.maskBasename[0] = '\0';
  c.discontinuityBasename[0] = '\0';
  c.strictMasking = 0;
  c.inputMapName[0] = '\0';
  c.cptsName[0] = '\0';
  c.constrainingMapName[0] = '\0';
  c.outputMapBasename[0] = '\0';
  c.outputImagesName[0] = '\0';
  c.outputCorrelationName[0] = '\0';
  c.logBasename[0] = '\0';
  c.startLevel = -1;
  c.outputLevel = 0;
  c.minResolution = 1;
  c.distortion = 1.0;
  c.correspondence = 1.0;
  c.correspondenceThreshold = 0.0;
  c.constraining = 1.0;
  c.constrainingThreshold = 0.0;
  c.quality = 10.0;
  c.minOverlap = 20.0;
  c.correlationHalfWidth = 7;
  c.writeAllMaps = 0;
  c.depth = 0;
  c.update = 0;
  c.partial = 0;
  r.pair.imageName[0]  = r.pair.imageName[1] = NULL;
  r.pair.pairName = NULL;
  if (r.message == NULL)
    r.message = (char *) malloc(PATH_MAX + 1024);
  forward = 1;
  pairsFile[0] = '\0';
  nPairs = 0;
  pairs = 0;
  minSlice = 0;
  maxSlice = 1000000000;
  memset(dirHash, 0, DIR_HASH_SIZE * sizeof(char*));

  for (i = 1; i < argc; ++i)
    if (strcmp(argv[i], "-depth") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d", &c.depth) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-min_res") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d", &c.minResolution) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-input") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.imageBasename, argv[i]);
      }
    else if (strcmp(argv[i], "-mask") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.maskBasename, argv[i]);
      }
    else if (strcmp(argv[i], "-input_map") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.inputMapName, argv[i]);
      }
    else if (strcmp(argv[i], "-cpts") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.cptsName, argv[i]);
      }
    else if (strcmp(argv[i], "-constraining_map") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.constrainingMapName, argv[i]);
      }
    else if (strcmp(argv[i], "-output") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.outputMapBasename, argv[i]);
      }
    else if (strcmp(argv[i], "-start_level") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d", &c.startLevel) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-output_level") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d", &c.outputLevel) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-output_images") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.outputImagesName, argv[i]);
      }
    else if (strcmp(argv[i], "-output_correlation") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.outputCorrelationName, argv[i]);
      }
    else if (strcmp(argv[i], "-logs") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.logBasename, argv[i]);
      }
    else if (strcmp(argv[i], "-pairs") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(pairsFile, argv[i]);
      }
    else if (strcmp(argv[i], "-summary") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(summaryName, argv[i]);
      }
    else if (strcmp(argv[i], "-distortion") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%lf", &c.distortion) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-correspondence") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%lf", &c.correspondence) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-correspondence_threshold") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%lf", &c.correspondenceThreshold) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-constraining") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%lf", &c.constraining) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-constraining_threshold") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%lf", &c.constrainingThreshold) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-quality") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%lf", &c.quality) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-correlation_half_width") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%d", &c.correlationHalfWidth) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-all_maps") == 0)
      c.writeAllMaps = 1;
    else if (strcmp(argv[i], "-forward") == 0)
      forward = 1;
    else if (strcmp(argv[i], "-reverse") == 0)
      forward = 0;
    else if (strcmp(argv[i], "-update") == 0)
      c.update = 1;
    else if (strcmp(argv[i], "-partial") == 0)
      c.partial = 1;
    else if (strcmp(argv[i], "-strict_masking") == 0)
      c.strictMasking = 1;
    else if (strcmp(argv[i], "-tif") == 0)
      c.type = 't';
    else if (strcmp(argv[i], "-vis") == 0)
      vis = 1;
    else error = 1;

  if (error)
    {
      fprintf(stderr, "Usage: warp3 -input file_prefix -output file_prefix\n");
      fprintf(stderr, "              [-skip skip_file] [-slice slice_range]\n");
      fprintf(stderr, "              [-resolution HORIZxVERT]\n");
      fprintf(stderr, "              [-distortion distortion_weight]\n");
      fprintf(stderr, "              [-correspondence correspondence_weight]\n");
      fprintf(stderr, "              [-correspondence_threshold pixels]\n");
      fprintf(stderr, "              [-quality quality_factor]\n");
      fprintf(stderr, "              [-min_res minimum_resolution_in_pixels]\n");
      fprintf(stderr, "              [-depth delta_depth]\n");
      fprintf(stderr, "              [-all_maps]\n");
      fprintf(stderr, "              [-update]\n");
      fprintf(stderr, "              [-partial]\n");
      fprintf(stderr, "              [-pairs <pair_file>]\n");
      fprintf(stderr, "              [-input_map <input_map_prefix>]\n");
      fprintf(stderr, "              [-constraining_map <constraining_map_prefix>]\n");
      fprintf(stderr, "              [-constraining constraining_weight]\n");
      fprintf(stderr, "              [-constraining_threshold pixels]\n");
      fprintf(stderr, "              [-strict_masking]\n");
      fprintf(stderr, "   where ranges are expressed as: integer\n");
      fprintf(stderr, "                              or: integer-integer\n");
      fprintf(stderr, "                              or: integer-\n");
      fprintf(stderr, "                              or: -integer\n");
      exit(1);
    }

  Log("MASTER starting on node %d\n", par_instance());

  Log("argc = %d\n", argc);
  for (i = 0; i < argc; ++i)
    Log("ARGV[%d] = %s\n", i, argv[i]);

  /* check that at least minimal parameters were supplied */
  if (c.imageBasename[0] == '\0' || c.outputMapBasename[0] == '\0' ||
      pairsFile[0] == '\0')
    Error("-input, -output, and -pairs parameters must be specified.\n");

  f = fopen(pairsFile, "r");
  if (f == NULL)
    Error("Could not open pairs file %s\n", pairsFile);
  while (fscanf(f, "%s %d %d %d %d %s %d %d %d %d %s",
		imgn[0], &imgMinX[0], &imgMaxX[0], &imgMinY[0], &imgMaxY[0],
		imgn[1], &imgMinX[1], &imgMaxX[1], &imgMinY[1], &imgMaxY[1],
		pairn) == 11)
    {
      if ((nPairs & 1023) == 0)
	pairs = (Pair *) realloc(pairs, (nPairs + 1024) * sizeof(Pair));
      for (imi = 0; imi < 2; ++imi)
	{
	  pairs[nPairs].imageName[imi] = (char *) malloc(strlen(imgn[imi])+1);
	  strcpy(pairs[nPairs].imageName[imi], imgn[imi]);
	  pairs[nPairs].imageMinX[imi] = imgMinX[imi];
	  pairs[nPairs].imageMaxX[imi] = imgMaxX[imi];
	  pairs[nPairs].imageMinY[imi] = imgMinY[imi];
	  pairs[nPairs].imageMaxY[imi] = imgMaxY[imi];
	}
      pairs[nPairs].pairName = (char *) malloc(strlen(pairn)+1);
      strcpy(pairs[nPairs].pairName, pairn);
      ++nPairs;
    }
  fclose(f);

  /* check that output directories are writeable */ 
  Log("MASTER checking output directories\n");
  sprintf(fn, "%sTEST.map", c.outputMapBasename);
  f = fopen(fn, "w");
  if (f == NULL)
    {
      for (i = strlen(c.outputMapBasename)-1; i >= 0 && c.outputMapBasename[i] != '/'; --i) ;
      if (i != 0)
	strncpy(outputDirName, c.outputMapBasename, i+1);
      else
	outputDirName[0] = '.';
      outputDirName[i+1] = '\0';
      Error("Could not open test output file %s --\n       does directory %s exist and is it writeable?\n", fn, outputDirName);
    }
  fclose(f);
  unlink(fn);

  if (c.outputImagesName[0] != '\0')
    {
      sprintf(fn, "%sTEST.pgm", c.outputImagesName);
      f = fopen(fn, "w");
      if (f == NULL)
	{
	  for (i = strlen(c.outputImagesName)-1; i >= 0 && c.outputImagesName[i] != '/'; --i) ;
	  if (i != 0)
	    strncpy(outputDirName, c.outputImagesName, i+1);
	  else
	    outputDirName[0] = '.';
	  outputDirName[i+1] = '\0';
	  Error("Could not open test output image file %s --\n      does directory %s exist and is it writeable?\n", fn, outputDirName);
	}
      fclose(f);
      unlink(fn);
    }

  Log("MASTER setting context\n");

  par_set_context();

  Log("MASTER set context\n");

  /* for all slices */
  printf("Processing slices: ");
  fflush(stdout);
  for (imi = 0; imi < 2; ++imi)
    t.pair.imageMinX[imi] = t.pair.imageMaxX[imi] = t.pair.imageMinY[imi] = t.pair.imageMaxY[imi] = -1;
  Log("nPairs = %d\n", nPairs);
  for (pn = 0; pn < nPairs; ++pn)
    {
      memcpy(&(t.pair), &(pairs[pn]), sizeof(Pair));

      // make sure that output directories exist
      sprintf(fn, "%s%s.map", c.outputMapBasename, t.pair.pairName);
      if (!CreateDirectories(fn))
	continue;
      if (c.outputImagesName[0] != '\0')
	{
	  sprintf(fn, "%s%s.pgm", c.outputImagesName, t.pair.pairName);
	  if (!CreateDirectories(fn))
	    continue;
	}
      if (c.outputCorrelationName[0] != '\0')
	{
	  sprintf(fn, "%s%s.pgm", c.outputCorrelationName, t.pair.pairName);
	  if (!CreateDirectories(fn))
	    continue;
	}

      Log("Delegating pair %d\n", pn);
      par_delegate_task();
    }
  par_finish();

  if (summaryName[0] != '\0')
    {
      f = fopen(summaryName, "w");
      if (f == NULL)
	Error("Could not open summary output file %s\n", summaryName);

      qsort(results, nResults, sizeof(Result), SortBySlice);
      fprintf(f, "Sorted by slice:\n");
      fprintf(f, "         IMAGE     REFERENCE   CORRELATION    DISTORTION    CORRESPOND    CONSTRAIN     ENERGY\n");
      for (i = 0; i < nResults; ++i)
	fprintf(f, "%s  %s  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f\n",
		results[i].pair.imageName[0],
		results[i].pair.imageName[1],
		results[i].correlation,
		results[i].distortion,
		results[i].correspondence,
		results[i].constraining,
		results[i].distortion * c.distortion - results[i].correlation +
		results[i].correspondence * c.correspondence + results[i].constraining * c.constraining);
      fprintf(f, "\n\nSorted by energy:\n");
      fprintf(f, "         IMAGE     REFERENCE   CORRELATION    DISTORTION    CORRESPOND    CONSTRAIN     ENERGY\n");
      qsort(results, nResults, sizeof(Result), SortByEnergy);
      for (i = 0; i < nResults; ++i)
	fprintf(f, "%s  %s  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f\n",
		results[i].pair.imageName[0],
		results[i].pair.imageName[1],
		results[i].correlation,
		results[i].distortion,
		results[i].correspondence,
		results[i].constraining,
		results[i].distortion * c.distortion - results[i].correlation +
		results[i].correspondence * c.correspondence + results[i].constraining * c.constraining);
      fclose(f);
    }

  printf(" %d\nAll slices completed.\n", nResults);
}

void
MasterResult ()
{
  if (r.message[0] != '\0')
    Error("\nThe following error was encountered by one of the worker processes:%s\n", r.message);

  results = (Result*) realloc(results, (nResults + 1) * sizeof(Result));
  results[nResults] = r;
  results[nResults].message = NULL;
  r.pair.imageName[0] = NULL;
  r.pair.imageName[1] = NULL;
  r.pair.pairName = NULL;

  if ((nResults % 50) == 0 && nResults != 0)
    printf(" %d \n                   ", nResults);
  printf(".");
  ++nResults;
}

int
Compare (const void *x, const void *y)
{
  if (*((int *) x) < *((int *) y))
    return(-1);
  else if (*((int *) x) == *((int *) y))
    return(0);
  else
    return(1);
}

int
SortBySlice (const void *x, const void *y)
{
  Result *rx, *ry;
  rx = (Result *) x;
  ry = (Result *) y;
  return(strcmp(rx->pair.imageName[0], ry->pair.imageName[0]));
}

int
SortByEnergy (const void *x, const void *y)
{
  Result *rx, *ry;
  double ex, ey;
  rx = (Result *) x;
  ry = (Result *) y;
  ex = rx->distortion * c.distortion - rx->correlation + rx->correspondence * c.correspondence +
    rx->constraining * c.constraining;
  ey = ry->distortion * c.distortion - ry->correlation + ry->correspondence * c.correspondence +
    ry->constraining * c.constraining;
  if (ex > ey)
    return(-1);
  else if (ex == ey)
    return(0);
  else
    return(1);
}

int
ParseRange (char *s, int *pos, int *minValue, int *maxValue)
{
  int i, j;
  int v;
  int minSpecified, maxSpecified;

  i = *pos;
  v = 0;
  while (s[i] != '\0' && s[i] >= '0' && s[i] <= '9')
    v = 10 * v + (s[i++] - '0');
  minSpecified = (i != *pos);
  *minValue = v;
  if (s[i] == '-')
    {
      j = ++i;
      v = 0;
      while (s[i] != '\0' && s[i] >= '0' && s[i] <= '9')
	v = 10 * v + (s[i++] - '0');
      maxSpecified = (i != j);
      if (maxSpecified)
	*maxValue = v;
      else
	*maxValue = 0x7fffffff;
    }
  else
    if (minSpecified)
      *maxValue = v;
    else
      return(0);
  *pos = i;
  return(1);
}

int
ParseValue (char *s, int *pos, int *value)
{
  int i, j;
  int v;

  i = *pos;
  v = 0;
  while (s[i] != '\0' && s[i] >= '0' && s[i] <= '9')
    v = 10 * v + (s[i++] - '0');
  if (i == *pos)
    return(0);
  *pos = i;
  return(1);
}

/* WORKER PROCEDURES */

void
WorkerContext ()
{
  Log("WorkerContext called\n");
  r.pair.imageName[0] = NULL;
  r.pair.imageName[1] = NULL;
  r.pair.pairName = NULL;
  if (r.message == NULL)
    r.message = (char *) malloc(PATH_MAX + 1024);
  r.message[0] = '\0';
  t.pair.imageName[0] = NULL;
  t.pair.imageName[1] = NULL;
  t.pair.pairName = NULL;
}

void
WorkerTask ()
{
  FILE *f;
  unsigned int w, h;
  int i;
  int n;
  char fn[PATH_MAX];
  int m;
  int mx, my;
  int x, y;
  int irx, iry;
  double rx, ry;
  double rrx, rry;
  size_t imbpl;
  float *warp;
  unsigned char *warpout;
  int level;
  int imw, imh;
  int mpw, mph;
  int maskPresent;
  float imx, imy;
  float refx, refy;
  char line[LINE_LENGTH+1];
  struct stat sb;
  int computeMap;
  double outputTime;
  double energy;
  size_t imagePixels;
  size_t ii;
  int mapHeader[3];
  char imageName[2][PATH_MAX], maskName[2][PATH_MAX], discontinuityName[2][PATH_MAX];
  char cptsName[PATH_MAX], inputMapName[PATH_MAX], constrainingMapName[PATH_MAX];
  char outputName[PATH_MAX], outputImagesName[PATH_MAX], outputCorrelationName[PATH_MAX];
  int cw, ch;
  int factor;
  int dx, dy;
  unsigned char *correlation;
  double xv, yv;
  int ixv, iyv;
  double rv, iv;
  double mi, mr;
  char errorMsg[PATH_MAX+256];
  int imi;
  unsigned char *image_in;
  float *img;
  cpu_set_t cpumask;

  Log("WORKER starting on node %d\n", par_instance());
  //  CPU_ZERO(&cpumask);
  //  CPU_SET(getpid() & 7, &cpumask);
  //  sched_setaffinity(0, sizeof(cpu_set_t), &cpumask);

  /* construct filenames */
  for (imi = 0; imi < 2; ++imi)
    {
      sprintf(imageName[imi], "%s%s.%s", c.imageBasename, t.pair.imageName[imi],
	      c.type == 't' ? "tif" : "pgm");
      if (c.maskBasename[0] != '\0')
	sprintf(maskName[imi], "%s%s.pbm", c.maskBasename, t.pair.imageName[imi]);
      else
	maskName[imi][0] = '\0';
      if (c.discontinuityBasename[0] != '\0')
	sprintf(discontinuityName[imi], "%s%s.pbm", c.discontinuityBasename, t.pair.imageName[imi]);
      else
	discontinuityName[imi][0] = '\0';
    }
  if (c.cptsName[0] != '\0')
    sprintf(cptsName, "%s%s.pts", c.cptsName, t.pair.pairName);
  else
    cptsName[0] = '\0';
  if (c.inputMapName[0] != '\0')
    sprintf(inputMapName, "%s%s.map", c.inputMapName, t.pair.pairName);
  else
    inputMapName[0] = '\0';
  if (c.constrainingMapName[0] != '\0')
    sprintf(constrainingMapName, "%s%s.map", c.constrainingMapName, t.pair.pairName);
  else
    constrainingMapName[0] = '\0';
  sprintf(outputName, "%s%s", c.outputMapBasename, t.pair.pairName);
  if (c.outputImagesName[0] != '\0')
    sprintf(outputImagesName, "%s%s.pgm", c.outputImagesName, t.pair.pairName);
  else
    outputImagesName[0] = '\0';
  if (c.outputCorrelationName[0] != '\0')
    sprintf(outputCorrelationName, "%s%s.pgm", c.outputCorrelationName, t.pair.pairName);
  else
    outputCorrelationName[0] = '\0';

  /* check if we can skip this task because of the -update option */
  computeMap = 0;
  if (!c.update)
    computeMap = 1;
  /* check if output map file already exists */
  sprintf(fn, "%s.score", outputName);
  if (stat(fn, &sb) == 0)
    outputTime = (double) sb.st_mtime;
  else
    computeMap = 1;
  for (imi = 0; imi < 2; ++imi)
    {
      if (!computeMap)
	{
	  if (stat(imageName[imi], &sb) == 0 &&
	      (double) (sb.st_mtime) > outputTime)
	    computeMap = 1;
	}
      if (!computeMap && maskName[imi][0] != '\0')
	{
	  if (stat(maskName[imi], &sb) == 0 &&
	      (double) (sb.st_mtime) > outputTime)
	    computeMap = 1;
	}
      if (!computeMap && discontinuityName[imi][0] != '\0')
	{
	  if (stat(discontinuityName[imi], &sb) == 0 &&
	      (double) (sb.st_mtime) > outputTime)
	    computeMap = 1;
	}	  
    }
  if (!computeMap && cptsName[0] != '\0')
    {
      if (stat(cptsName, &sb) == 0 &&
	  (double) (sb.st_mtime) > outputTime)
	computeMap = 1;
    }
  if (!computeMap && inputMapName[0] != '\0')
    {
      if (stat(inputMapName, &sb) == 0 &&
	  (double) (sb.st_mtime) > outputTime)
	computeMap = 1;
    }
  if (!computeMap && constrainingMapName[0] != '\0')
    {
      if (stat(constrainingMapName, &sb) == 0 &&
	  (double) (sb.st_mtime) > outputTime)
	computeMap = 1;
    }
  if (!computeMap)
    {
      sprintf(fn, "%s.score", outputName);
      f = fopen(fn, "r");
      if (f == NULL ||
	  fscanf(f, "%lf%lf%lf%lf%lf",
		 &energy,
		 &(r.correlation),
		 &(r.distortion),
		 &(r.correspondence),
		 &(r.constraining)) != 5)
	computeMap = 1;
      if (f != NULL)
	fclose(f);
    }
  if (!computeMap)
    {
      Log("WORKER TASK skipping %s to %s since up-to-date\n",
	  imageName[0], imageName[1]);
      memcpy(&(r.pair), &(t.pair), sizeof(Pair));
      r.updated = 0;
      r.message[0] = '\0';
      return;
    }

  if (c.partial)
    {
      computeMap = 1;
      for (imi = 0; imi < 2; ++imi)
	{
	  if ((f = fopen(imageName[imi], "r")) != NULL)
	    fclose(f);
	  else
	    {
	      computeMap = 0;
	      break;
	    }
	  if (maskName[imi][0] != '\0')
	    {
	      if ((f = fopen(maskName[imi], "r")) != NULL)
		fclose(f);
	      else
		{
		  computeMap = 0;
		  break;
		}
	    }
	}
      if (!computeMap)
	{
	  Log("WORKER TASK skipping %s to %s since input files missing\n",
	      imageName[0], imageName[1]);
	  memcpy(&(r.pair), &(t.pair), sizeof(Pair));
	  r.updated = 0;
	  r.message[0] = '\0';
	  return;
	}
    }

  Log("STARTING WORKER TASK\n");

  inputMap = NULL;
  constrainingMap = NULL;

  for (imi = 0; imi < 2; ++imi)
    {
      Log("WORKER reading image %s\n", imageName[imi]);

      image_in = NULL;
      if (!ReadImage(imageName[imi], &image_in,
		     &imageWidth[imi][0], &imageHeight[imi][0],
		     t.pair.imageMinX[imi], t.pair.imageMaxX[imi],
		     t.pair.imageMinY[imi], t.pair.imageMaxY[imi],
		     r.message))
	{
	  sprintf(r.message, "Could not read reference image %s\n", imageName[imi]);
	  return;
	}
      Log("image_in = %llx  w = %d h = %d\n", image_in,
	  imageWidth[imi][0], imageHeight[imi][0]);
      imageOffsetX[imi][0] = (t.pair.imageMinX[imi] >= 0) ? t.pair.imageMinX[imi] : 0;
      imageOffsetY[imi][0] = (t.pair.imageMinY[imi] >= 0) ? t.pair.imageMinY[imi] : 0;
      Log("WORKER READ IMAGE: %s\n", imageName[imi]);

      imagePixels = ((size_t) imageWidth[imi][0]) * imageHeight[imi][0];
      images[imi][0] = (float *) malloc(imagePixels * sizeof(float));
      if (images[imi][0] == NULL)
	{
	  sprintf(r.message, "Could not allocate image arrays\n");
	  return;
	}
      Log("images[0] = %llx  imp = %llu\n", (long long) images[imi][0], imagePixels);
      img = images[imi][0];
      for (ii = 0; ii < imagePixels; ++ii)
	img[ii] = image_in[ii];
      free(image_in);
      image_in = NULL;

#if MASKING
      maskPresent = 0;
      imbpl = (imageWidth[imi][0] + 7) >> 3;
      if (maskName[imi][0] != '\0')
	{
	  masks[imi][0] = NULL;
	  if (!ReadBitmap(maskName[imi], &masks[imi][0],
			  &imw, &imh,
			  t.pair.imageMinX[imi], t.pair.imageMaxX[imi],
			  t.pair.imageMinY[imi], t.pair.imageMaxY[imi],
			  r.message))
	    {
	      sprintf(r.message, "Could not read reference mask %s\n", maskName);
	      return;
	    }
	  Log("WORKER READ THE IMAGE MASK\n");

	  Log("imw[%d] = %d imh[%d] = %d imageWidth[%d][0] = %d imageHeight[%d][0] = %d\n",
	      imi, imw, imi, imh, imi, imageWidth[imi][0], imi, imageHeight[imi][0]);

	  if (imw != imageWidth[imi][0] ||
	      imh != imageHeight[imi][0])
	    {
	      sprintf(r.message, "Error: incorrect reference mask size: %s\n", maskName[imi]);
	      return;
	    }
	  
	  img = images[imi][0];
	  for (my = 0; my < imh; ++my)
	    for (mx = 0; mx < imw; ++mx)
	      if (!MASK(masks[imi][0], imbpl, mx, my))
		IMAGE(img, imw, mx, my) = 0;
	  Log("WORKER APPLIED THE IMAGE MASK\n");
	  maskPresent = 1;
	}
      if (!maskPresent)
	{
	  masks[imi][0] = (unsigned char *) malloc(imageHeight[imi][0]*imbpl);
	  if (masks[imi][0] == NULL)
	    {
	      sprintf(r.message, "Could not allocate image arrays\n");
	      return;
	    }
	  memset(masks[imi][0], 0xff, imageHeight[imi][0]*imbpl);
	}
#endif

#if FOLDING
      imbpl = (imageWidth[imi][0] + 7) >> 3;
      if (discontinuityName[imi][0] != '\0')
	{
	  if (!ReadBitmap(discontinuityName, &idisc[imi][0],
			  &imw, &imh,
			  t.pair.imageMinX[imi], t.pair.imageMaxX[imi],
			  t.pair.imageMinY[imi], t.pair.imageMaxY[imi],
			  r.message))
	    {
	      sprintf(r.message, "Could not read reference discontinuity mask %s\n", discontinuityName);
	      return;
	    }
	  Log("WORKER READ THE DISCONTINUITY MASK\n");
	  if (imw != imageWidth[imi][0] || imh != imageHeight[imi][0])
	    {
	      sprintf(r.message, "Error: incorrect reference discontinuity mask size: %s\n", discontinuityName);
	      return;
	    }
	}
#endif
    }

  /* read in the correspondence points */
  if (cpts != NULL)
    {
      nCpts = 0;
      free(cpts);
      cpts = NULL;
    }
  if (cptsName[0] != '\0')
    {
      f = fopen(cptsName, "r");
      if (f != NULL)
	{
	  while (fgets(line, LINE_LENGTH, f) != NULL)
	    {
	      sscanf(line, "%f%f%f%f", &imx, &imy, &refx, &refy);
	      cpts = (CPoint*) realloc(cpts, (nCpts+1)*sizeof(CPoint));
	      cpts[nCpts].ix = imx;
	      cpts[nCpts].iy = imy;
	      cpts[nCpts].rx = refx;
	      cpts[nCpts].ry = refy;
	      ++nCpts;
	    }
	  fclose(f);
	}
    }

  /* get the input map, if present */
  if (inputMapName[0] != '\0')
    {
      if (!ReadMap(inputMapName,
		   &inputMap,
		   &inputMapLevel,
		   &mpwi, &mphi,
		   &offxi, &offyi,
		   NULL, NULL,
		   errorMsg))
	Error("Error reading map:\n %s\n", errorMsg);
      inputMapFactor = 1 << inputMapLevel;
      Log("Read input map %s of size %d x %d;  inputMapFactor = %d\n",
	  inputMapName, mpwi, mphi, inputMapFactor);
    }
  else
    inputMapFactor = 0;

  /* read in the constraining map, if present */
  if (constrainingMapName[0] != '\0')
    {
      if (!ReadMap(constrainingMapName,
		   &constrainingMap,
		   &constrainingMapLevel,
		   &mpwc, &mphc,
		   &offxc, &offyc,
		   NULL, NULL,
		   errorMsg))
	Error("Error reading map:\n %s\n", errorMsg);
      constrainingMapFactor = 1 << constrainingMapLevel;
      Log("Read constraining map %s of size %d x %d;  constrainingMapFactor = %d\n",
	  constrainingMapName, mpwc, mphc, constrainingMapFactor);
    }
  else
    constrainingMapFactor = 0;

  if (!Init())
    {
      Log("Init was unsuccessful.\n");
      return;
    }

  Compute(outputName);

  if (inputMap != NULL)
    {
      free(inputMap);
      inputMap = NULL;
    }
  if (constrainingMap != NULL)
    {
      free(constrainingMap);
      constrainingMap = NULL;
    }
#if FOLDING
  free(idisc);
#endif
  for (level = 0; level < nLevels; ++level)
    {
      for (imi = 0; imi < 2; ++imi)
	{
	  free(images[imi][level]);
#if MASKING
	  free(masks[imi][level]);
#endif
#if FOLDING
	  free(imdisc[imi][level]);
#endif
	}
      free(maps[level]);
    }
}

int
Init ()
{
  int level;
  unsigned int iw, ih;
  float *image;
  size_t imagePixels;
  size_t ii;
  int i;
  int x, y;
  size_t imbpl;
  int b;
  int dx, dy;
  unsigned char *iCount;
  size_t impw, imph;
  size_t impbpl;
  int imi;
  size_t maskCount;
  int minX, maxX, minY, maxY;
  float *src_image;
  int src_iw, src_ih;
  int src_x, src_y;
  unsigned char *src_mask;
  int src_mbpl;
  int delta_x, delta_y;
  int ix, iy;
  unsigned char *m;
  int factor;
  int mapMinX, mapMaxX, mapMinY, mapMaxY;

  for (imi = 0; imi < 2; ++imi)
    {
      imagePixels = ((size_t) imageWidth[imi][0]) * imageHeight[imi][0];
      Log("iw[%d] = %d ih[%d] = %d pixels = %llu\n",
	  imi, imageWidth[imi][0], imi, imageHeight[imi][0], imagePixels);

      iw = imageWidth[imi][0];
      ih = imageHeight[imi][0];
      imbpl = (iw + 7) >> 3;
#if MASKING
      Log("maskcount = %ld\n", CountBits(masks[imi][0], ih * imbpl));
#endif
    }
  mapWidth[0] = imageWidth[0][0];
  mapHeight[0] = imageHeight[0][0];
  mapOffsetX[0] = imageOffsetX[0][0];
  mapOffsetY[0] = imageOffsetY[0][0];

  /* create hierarchical representations of images */
  for (level = 1; mapWidth[level-1] > 2 || mapHeight[level-1] > 2; ++level)
    {
      Log("creating level %d\n", level);
      factor = 1 << level;

      mapMinX = imageOffsetX[0][0] / factor;
      mapMaxX = (imageOffsetX[0][0] + imageWidth[0][0] - 1) / factor + 1;
      mapMinY = imageOffsetY[0][0] / factor;
      mapMaxY = (imageOffsetY[0][0] + imageHeight[0][0] - 1) / factor + 1;
      mapWidth[level] = mapMaxX - mapMinX + 1;
      mapHeight[level] = mapMaxY - mapMinY + 1;
      mapOffsetX[level] = mapMinX;
      mapOffsetY[level] = mapMinY;
      Log("map dims: %d x %d + %d + %d\n", mapWidth[level], mapHeight[level],
	  mapOffsetX[level], mapOffsetY[level]);

      for (imi = 0; imi < 2; ++imi)
	{
	  // calculate the ranges for level-1
	  minX = imageOffsetX[imi][level-1];
	  maxX = minX + imageWidth[imi][level-1] - 1;
	  minY = imageOffsetY[imi][level-1];
	  maxY = minY + imageHeight[imi][level-1] - 1;

	  // now calculate the ranges for level
	  minX = (minX + 1) / 2;
	  maxX = (maxX + 1) / 2 - 1;
	  minY = (minY + 1) / 2;
	  maxY = (maxY + 1) / 2 - 1;

	  if (minX > maxX)
	    minX = maxX = imageOffsetX[imi][level-1] / 2;
	  if (minY > maxY)
	    minY = maxY = imageOffsetY[imi][level-1] / 2;

	  iw = imageWidth[imi][level] = maxX - minX + 1;
	  ih = imageHeight[imi][level] = maxY - minY + 1;
	  imageOffsetX[imi][level] = minX;
	  imageOffsetY[imi][level] = minY;
	  imagePixels = ((size_t) imageWidth[imi][level]) * imageHeight[imi][level];
	  imbpl = (iw + 7) >> 3;
	  impw = iw + 1;
	  imph = ih + 1;
	  impbpl = (impw  + 7) >> 3;

	  Log("imagePixels = %d\n", imagePixels);
	  images[imi][level] = (float*) malloc(imagePixels * sizeof(float));
	  if (images[imi][level] == NULL)
	    {
	      sprintf(r.message, "Could not allocate image arrays\n");
	      return;
	    }

#if MASKING
	  masks[imi][level] = (unsigned char*) malloc(ih*imbpl);
	  if (masks[imi][level] == NULL)
	    {
	      sprintf(r.message, "Could not allocate mask arrays\n");
	      return;
	    }

	  iCount = (unsigned char*) malloc(imagePixels * sizeof(unsigned char));
	  if (iCount == NULL)
	    {
	      sprintf(r.message, "Could not allocate count arrays\n");
	      return;
	    }
#endif

	  image = images[imi][level];
	  memset(image, 0, imagePixels * sizeof(float));
#if MASKING
	  memset(iCount, 0, imagePixels*sizeof(unsigned char));
#endif

	  src_image = images[imi][level-1];
	  src_ih = imageHeight[imi][level-1];
	  src_iw = imageWidth[imi][level-1];
	  src_mask = masks[imi][level-1];
	  src_mbpl = (src_iw + 7) >> 3;
	  delta_x = 2*imageOffsetX[imi][level] - imageOffsetX[imi][level-1];
	  delta_y = 2*imageOffsetY[imi][level] - imageOffsetY[imi][level-1];
	  for (y = 0; y < ih; ++y)
	    {
	      iy = 2*y + delta_y;
	      for (x = 0; x < iw; ++x)
		{
		  ix = 2 * x + delta_x;
		  for (dy = 0; dy < 2; ++dy)
		    {
		      src_y = iy + dy;
		      if (src_y < 0 || src_y >= src_ih)
			continue;
		      for (dx = 0; dx < 2; ++dx)
			{
			  src_x = ix + dx;
			  if (src_x < 0 || src_x >= src_iw)
			    continue;
			  IMAGE(image, iw, x, y) += IMAGE(src_image, src_iw, src_x, src_y);
#if MASKING
			  if (MASK(src_mask, src_mbpl, ix+dx, iy+dy))
			    ++IMAGE(iCount, iw, x, y);
#endif
			}
		    }
#if MASKING
		  if (IMAGE(iCount, iw, x, y) > 0)
		    IMAGE(image, iw, x, y) *= 1.0 / IMAGE(iCount, iw, x, y);
		  else
		    IMAGE(image, iw, x, y) = 0.0;
#else
		  IMAGE(image, iw, x, y) *= 1.0 / 4.0;
#endif		  
		}
	    }

#if MASKING
	  m = masks[imi][level];
	  memset(m, 0, ih*imbpl);
	  if (c.strictMasking)
	    {
	      for (y = 0; y < ih; ++y)
		for (x = 0; x < iw; ++x)
		  if (IMAGE(iCount, iw, x, y) == 4)
		    SETMASKBIT(m, imbpl, x, y);
	    }
	  else
	    {
	      for (y = 0; y < ih; ++y)
		for (x = 0; x < iw; ++x)
		  if (IMAGE(iCount, iw, x, y) > 0)
		    SETMASKBIT(m, imbpl, x, y);
	    }
#endif

#if FOLDING
	  imdisc[imi][level] = (unsigned char*) malloc(imph * impbpl);
	  if (imdisc[imi][level] == NULL)
	    {
	      sprintf(r.message, "Could not allocate discontinuity arrays\n");
	      return;
	    }
	  memset(imdisc[imi][level], 0, imph * impbpl);

	  // allocate a bitmap array for one map square at this level
	  dmadim = 1 << level;
	  dmabpl = (dmadim + 7) >> 3;
	  dmasize = dmadim * dmabpl;
	  dma = (unsigned char *) malloc(dmasize);
	  if (dma == NULL)
	    {
	      sprintf(r.message, "Could not allocate arrays\n");
	      return;
	    }
	  for (y = 0; y < ih; ++y)
	    for (x = 0; x < iw; ++x)
	      {
		ox = x * dmadim + (dmadim >> 1);
		oy = y * dmadim + (dmadim >> 1);
		memset(dma, 0, dmasize);
		first = 0;
		for (i = 0; i < 4; ++i)
		  {
		    switch (i)
		      {
		      case 0:
			ix = 0;
			iy = 0;
			break;
		      case 1:
			ix = dmadim - 1;
			iy = 0;
			break;
		      case 2:
			ix = dmadim - 1;
			iy = dmadim - 1;
			break;
		      case 3:
			ix = 0;
			iy = dmadim - 1;
			break;
		      }
		    if (ox + ix < 0 || ox + ix >= imageWidth[0] ||
			oy + iy < 0 || oy + iy >= imageHeight[0])
		      continue;
		    if (first)
		      {
			FloodFill(ix, iy, dmadim, dmadim,
				  idisc, imageWidth[0], imageHeight[0],
				  ox, oy);
			first = 0;
		      }
		    else if (dma[y*dmabpl + (ix >> 3)] & (0x80 >> (ix & 7)) == 0)
		      {
			imdisc[level][(y+1)*impbpl + ((x+1) >> 3)] |=
			  0x80 >> ((x+1) & 7);
			break;
		      }
		  }
	      }

	  for (y = 0; y < referenceHeight[level-1]; ++y)
	    for (x = 0; x < referenceWidth[level-1]; ++x)
	      if ((rdisc[level-1][y*rmbpl1 + (x >> 3)] & (0x80 >> (x & 7))) != 0)
		++rCount[(y >> 1) * ((size_t) rw) + (x >> 1)];

	  memset(rdisc[level], 0, rh*rmbpl);

	  for (y = 0; y < ih; ++y)
	    for (x = 0; x < iw; ++x)
	      if (iCount[y * ((size_t) iw) + x] == 4)
		idisc[level][y * imbpl + (x >> 3)] |= 0x80 >> (x & 7);
#endif

#if MASKING
	  free(iCount);
#endif
	}
    }
  nLevels = level;
  Log("nLevels = %d\n", nLevels);

  if (c.outputLevel < 0 || c.outputLevel >= nLevels)
    Error("outputLevel (%d) is invalid for image (nLevels = %d)\n",
	  c.outputLevel, nLevels);

#if MASKING
  for (level = 0; level < nLevels; ++level)
    {
      minMaskCount[level] = 1ULL << 62;
      for (imi = 0; imi < 2; ++imi)
	{
	  iw = imageWidth[imi][level];
	  ih = imageHeight[imi][level];
	  imagePixels = ((size_t) iw) * ih;
	  imbpl = (iw + 7) >> 3;

	  maskCount = CountBits(masks[imi][level], ih*imbpl);
	  if (maskCount < minMaskCount[level])
	    minMaskCount[level] = maskCount;
	  Log("maskCount[%d][%d] = %ld (%dx%d = %ld) %f%%\n",
	      imi, level, maskCount,
	      iw, ih, imagePixels, 100.0 * ((float) maskCount) / imagePixels);
	}
    }
#endif

#if 0
      /* make the full resolution image significance array by taking the standard
	 deviation of the pixel values over a circular region centered on each pixel */
      isig = (float *) malloc(imagePixels * sizeof(float));
#endif    

  Log("Constructed the multi-resolution image set.\n");
  return(1);
}


void
Compute (char *outputName)
{
  int level;
  MapElement *prop;
  int i;
  size_t ii;
  int x, y;
  double x00, x01, x10, x11;
  double y00, y01, y10, y11;
  double c00, c01, c10, c11;
  double rx, ry;
  int irx, iry;
  double rrx, rry;
  double r00, r01, r10, r11;
  double xd, yd;
  double sir, si2, sr2;
  double de;
  double dx, dy;
  double radius;
  double offset;
  double sqrtRadius, sqrtOffset;
  double theta;
  int changeMinX, changeMaxX, changeMinY, changeMaxY;
  int icx, icy;
  double cx, cy, cc;
  int xdir, ydir;
  int ix, iy;
  double dxp, dyp;
  double mrd;
  double nx, ny;
  double ode, nde;
  double newDistortion, newCorrelation, newCorrespondence, newConstraining;
  double newEnergy;
  double aCorrelation, aDistortion, aCorrespondence;
  int upd0, upd1, upd2, upd3;
  double iv, rv;
  double distortion, correlation, correspondence, constraining;
  double corr;
  int move;
  int nMoves;
  int accept;
  double prob;
  double ran;
  float *image;
  float *ref;
  float *warp;
  int mpw1, mph1;
  MapElement *map;
  MapElement *map1;
  double dsi, dsi2, dsir, dsr2, dsr;
  long dPoints;
  long newDPoints;
  unsigned int iw, ih;
  unsigned int rw, rh;
  size_t effectivePoints;
  double xi, yi;
  double si, sr;
  double mi, mr;
  double denom;
  double newsi, newsr;
  double newsi2, newsr2;
  double newsir;
  double newmi, newmr;
  double newDenom;
  size_t newPoints;
  size_t newEffectivePoints;
  double newCorr;
  long cPoints;
  size_t n;
  size_t updcnt;
  double x0, x1, x2, x3;
  double y0, y1, y2, y3;
  double c0, c1, c2, c3;
  double area;
  double nx0, nx1, nx2, nx3;
  double ny0, ny1, ny2, ny3;
  double nc0, nc1, nc2, nc3;
  double narea;
  size_t moveCount, acceptedMoveCount;
  double l0, l1, l2, l3;
  double nl0, nl1, nl2, nl3;
  double goalCoverage;
  double actualCoverage;
  size_t imbpl, rmbpl;
  unsigned char *mask;
  int cLevel;
  int nb;
  unsigned char *cimask, *crmask;
  size_t icmbpl, rcmbpl;
  float *cimage;
  float *cref;
  double rx00, rx01, rx10, rx11, ry00, ry01, ry10, ry11;
  size_t nPoints;
  size_t requiredPoints;
  int factor;
  int sx, sy, ex, ey;
  int mpw, mph;
  int mpw_minus_1, mph_minus_1;
  int ixv, iyv;
  double xv, yv;
  double energy;
  double ea;
  double check_sir, check_si2, check_sr2, check_si, check_sr;
  double check_sum;
  size_t check_points;
  size_t check_effectivePoints;
  double check_mi, check_mr, check_denom, check_r;
  double sqrtMaxMoveRadius, sqrtMaxMoveOffset;
  double kFactor;
  double lFactor;
  double distance;
  double ce;
  double ces;
  size_t mSize;
  double cth;
  double thetaRange;
  size_t statLogRadius[21];
  size_t statTheta[21];
  double statDeltaE[21];
  MapElement *mapCons;
  double oldEnergy;
  double logMinRadius, logMaxRadius, logRadiusRange;
  double minDist2;
  double d2;
  int ix0, iy0;
  double rc;
  double rc00, rc01, rc10, rc11;
  int multiplier;
  int startLevel;
  unsigned char *cidisc, *crdisc, *cirdisc;
  char fn[PATH_MAX];
  FILE *f;
  float rnd;
  float logRadius;
  int ncc;
  float scc;
  double udLimit, diagLimit;
  unsigned char *trimmed;
  MapElement *e, *ec, *ep;
  MapElement *mp;
  int imgox, imgoy;
  int refox, refoy;
  int mox1, moy1;
  float *dist;
  size_t nBytes;
  float threshold;
  int border;
  int rpbw, rpbh;
  size_t rpbbpl;
  unsigned char *rpbmask;
  int ddx, ddy;
  int valid;
  int mox, moy;
  int nTrimmed;
  int linePos;
  char line[65536];
  float *warpedArray;
  unsigned char *validArray;
  float *correlationArray;

  kFactor = 1.0 / sqrt((double) (((size_t) imageWidth[1][0]) * imageWidth[1][0] +
	                         ((size_t) imageHeight[1][0]) * imageHeight[1][0]));
  mapCons = NULL;
  Log("Warping started with %d total resolution levels\n", nLevels);

  if (c.startLevel >= 0)
    startLevel = c.startLevel;
  else if (inputMapFactor != 0)
    {
      startLevel = 0;
      while ((1 << startLevel) < inputMapFactor)
	++startLevel;
      if ((1 << startLevel) != inputMapFactor)
	{
	  Log("inputMapFactor not a power of 2 -- ignoring input map\n");
	  inputMapFactor = 0;
	  startLevel = nLevels - 1;
	}
     }
  else
    startLevel = nLevels - 1;

  /* go down hierarchy one level at a time */
  for (level = startLevel; level >= c.outputLevel; --level)
    {
      Log("Considering level %d\n", level);
      mpw = mapWidth[level];
      mph = mapHeight[level];
      mpw_minus_1 = mpw - 1;
      mph_minus_1 = mph - 1;
      mox = mapOffsetX[level];
      moy = mapOffsetY[level];
      lFactor = 1 << level;
#if FOLDING
      irdisc = (unsigned char *) malloc(imbpl * mph);
#endif

      cLevel = level - c.depth;
      if (cLevel < 0)
	cLevel = 0;
      while (cLevel > 0 &&
	     (imageWidth[0][cLevel] < c.minResolution ||
	      imageHeight[0][cLevel] < c.minResolution))
	--cLevel;
      nb = level - cLevel;
      factor = (1 << nb);
      ih = imageHeight[0][cLevel];
      iw = imageWidth[0][cLevel];
      rh = imageHeight[1][cLevel];
      rw = imageWidth[1][cLevel];
      imgox = imageOffsetX[0][cLevel];
      imgoy = imageOffsetY[0][cLevel];
      refox = imageOffsetX[1][cLevel];
      refoy = imageOffsetY[1][cLevel];
      icmbpl = (iw + 7) >> 3;
      rcmbpl = (rw + 7) >> 3;
      cimage = images[0][cLevel];
      cref = images[1][cLevel];
#if MASKING
      cimask = masks[0][cLevel];
      crmask = masks[1][cLevel];
#endif

#if 0
      // TEMPORARY
      if (level == 13)
	{
	  char *img = (unsigned char *) malloc(iw * ih);
	  int x, y;
	  for (y = 0; y < ih; ++y)
	    for (x = 0; x < iw; ++x)
	      img[y*iw+x] = (int) cimage[y*iw+x];
	  FILE *of = fopen("testimg.pgm", "w");
	  fprintf(of, "P5\n%d %d\n255\n", iw, ih);
	  fwrite(img, iw * ih, 1, of);
	  fclose(of);
	  free(img);
	  of= fopen("testref.pgm", "w");
	  fprintf(of, "P5\n%d %d\n255\n", rw, rh);
	  img = (unsigned char *) malloc(rw * rh);
	  for (y = 0; y < rh; ++y)
	    for (x = 0; x < rw; ++x)
	      img[y*rw+x] = (int) cref[y*rw+x];
	  fwrite(img, rw * rh, 1, of);
	  fclose(of);
	  free(img);
	}
#endif

      mSize = ((size_t) mpw) * mph;
      maps[level] = (MapElement*) malloc(mSize * sizeof(MapElement));
      if (maps[level] == NULL)
	{
	  sprintf(r.message, "Could not allocate map arrays\n");
	  return;
	}
      map = maps[level];
#if MASKING
      mask = masks[0][level];
#endif
      if (constrainingMapFactor != 0)
	{
	  if (mapCons != NULL)
	    free(mapCons);
	  mapCons = (MapElement *) malloc(mSize * sizeof(MapElement));
	  if (mapCons == NULL)
	    {
	      sprintf(r.message, "Could not allocate constraining arrays\n");
	      return;
	    }
	}

      if (level == startLevel)
	if (inputMapFactor != 0)
	  {
	    /* use the provided input map as the initial map */
	    for (y = 0; y < mph; ++y)
	      for (x = 0; x < mpw; ++x)
		{
		  /* translate into absolute image coordinates */
		  xv = (x + mox) * lFactor;
		  yv = (y + moy) * lFactor;

		  /* lookup in the map */
		  xv = (xv / inputMapFactor) - offxi;
		  yv = (yv / inputMapFactor) - offyi;
		  ixv = (int) floor(xv);
		  iyv = (int) floor(yv);
		  rrx = xv - ixv;
		  rry = yv - iyv;
#if 0
		  if (ixv < -1 || ixv >= mpwi ||
		      iyv < -1 || iyv >= mphi)
		    {
		      SETMAP(map, mpw, x, y, 0.0, 0.0, 0.0);
		      continue;
		    }
#endif
		  while (ixv < 0)
		    {
		      ++ixv;
		      rrx -= 1.0;
		    }
		  while (ixv >= mpwi-1)
		    {
		      --ixv;
		      rrx += 1.0;
		    }
		  while (iyv < 0)
		    {
		      ++iyv;
		      rry -= 1.0;
		    }
		  while (iyv >= mphi-1)
		    {
		      --iyv;
		      rry += 1.0;
		    }
		  GETMAP(inputMap, mpwi, ixv, iyv, &rx00, &ry00, &rc00);
		  GETMAP(inputMap, mpwi, ixv, iyv+1, &rx01, &ry01, &rc01);
		  GETMAP(inputMap, mpwi, ixv+1, iyv, &rx10, &ry10, &rc10);
		  GETMAP(inputMap, mpwi, ixv+1, iyv+1, &rx11, &ry11, &rc11);
		  rx = rx00 * (rrx - 1.0) * (rry - 1.0)
		    - rx10 * rrx * (rry - 1.0) 
		    - rx01 * (rrx - 1.0) * rry
		    + rx11 * rrx * rry;
		  ry = ry00 * (rrx - 1.0) * (rry - 1.0)
		    - ry10 * rrx * (rry - 1.0) 
		    - ry01 * (rrx - 1.0) * rry
		    + ry11 * rrx * rry;

		  /* translate into absolute image coordinates */
		  rx = inputMapFactor * rx;
		  ry = inputMapFactor * ry;

		  /* translate into current map coordinates */
		  SETMAP(map, mpw, x, y,
			 rx / lFactor, ry / lFactor, 1.0);
		}
	  }
	else
	  {
	    /* use the identity map; note that this displaces the positions
	       so that the upper left of the image subimage maps to the
	       upper left of the reference subimage (which is probably what
	       the user intended) */
	    for (y = 0; y < mph; ++y)
	      for (x = 0; x < mpw; ++x)
		SETMAP(map, mpw, x, y,
		       (float) (x + mox) + ((float) (refox - imgox)) / factor,
		       (float) (y + moy) + ((float) (refoy - imgoy)) / factor,
		       1.0);
	    Log("Setting initial map to the identity map (mpw=%d mph=%d)\n",
		mpw, mph);
	  }
      else
	{
	  /* interpolate/extrapolate to construct the current level map from
	     the next higher level map */
	  map1 = maps[level+1];
	  mpw1 = mapWidth[level+1];
	  mph1 = mapHeight[level+1];
	  mox1 = mapOffsetX[level+1];
	  moy1 = mapOffsetY[level+1];
	  for (y = 0; y < mph; ++y)
	    for (x = 0; x < mpw; ++x)
	      {
		/* translate into absolute image coordinates */
		xv = (x + mox) / 2.0 - mox1;
		yv = (y + moy) / 2.0 - moy1;
		
		/* lookup in the map */
		ixv = (int) floor(xv);
		iyv = (int) floor(yv);
		rrx = xv - ixv;
		rry = yv - iyv;
		while (ixv < 0)
		  {
		    ++ixv;
		    rrx -= 1.0;
		  }
		while (ixv >= mpw1-1)
		  {
		    --ixv;
		    rrx += 1.0;
		  }
		while (iyv < 0)
		  {
		    ++iyv;
		    rry -= 1.0;
		  }
		while (iyv >= mph1-1)
		  {
		    --iyv;
		    rry += 1.0;
		  }
		GETMAP(map1, mpw1, ixv, iyv, &rx00, &ry00, &rc00);
		GETMAP(map1, mpw1, ixv, iyv+1, &rx01, &ry01, &rc01);
		GETMAP(map1, mpw1, ixv+1, iyv, &rx10, &ry10, &rc10);
		GETMAP(map1, mpw1, ixv+1, iyv+1, &rx11, &ry11, &rc11);
		rx = rx00 * (rrx - 1.0) * (rry - 1.0)
		  - rx10 * rrx * (rry - 1.0) 
		  - rx01 * (rrx - 1.0) * rry
		  + rx11 * rrx * rry;
		ry = ry00 * (rrx - 1.0) * (rry - 1.0)
		  - ry10 * rrx * (rry - 1.0) 
		  - ry01 * (rrx - 1.0) * rry
		  + ry11 * rrx * rry;

		if (isnan(rx) || isnan(ry))
		  Error("ISNAN rx or ry\n");

		/* translate into current map coordinates */
		SETMAP(map, mpw, x, y, 2.0 * rx, 2.0 * ry, 1.0);
	      }
	  Log("Interpolated level %d map (mpw=%d mph=%d + boundary) to obtain level %d map (mpw=%d mph=%d)\n",
	      level+1, mpw1, mph1, level, mpw, mph);
	  if (vis)
	    sleep(3);
	}

      /* trim off exterior of map */
      /* remove any map point that is farther than sqrt(2)*2^depth pixels
	 from an unmasked image pixel */
      nTrimmed = 0;
      dist = (float *) malloc(((size_t) ih) * iw * sizeof(float));
      if (dist == NULL)
	{
	  sprintf(r.message, "Could not allocate dist array\n");
	  return;
	}
      computeDistance(EUCLIDEAN_DISTANCE_SQUARED, iw, ih, cimask, dist);

      threshold = 2.0 * factor * factor;
      for (y = 0; y < mph; ++y)
	{
	  iy = (y + moy) * factor - imgoy;
	  if (iy >= ih)
	    iy = ih-1;
	  for (x = 0; x < mpw; ++x)
	    {
	      ix = (x + mox) * factor - imgox;
	      if (ix >= iw)
		ix = iw-1;
	      if (IMAGE(dist, iw, ix, iy) > threshold)
		{
		  ++nTrimmed;
		  if (x == 210 && y == 175)
		    Log("SCB 0: %f %f ix=%d iy=%d mox=%d moy=%d iw=%d ih=%d mpw=%d mph=%d factor=%d imgox=%d imgoy=%d\n",
			IMAGE(dist, iw, ix, iy), threshold,
			ix, iy, mox, moy, iw, ih, mpw, mph,
			factor, imgox, imgoy);
		  MAP(map, mpw, x, y).c = 0.0;
		}
	    }
	}
      free(dist);
      Log("Trimmed %d map points in stage 1\n", nTrimmed);

      /* remove any map point that maps to a point farther than 2*2^depth
	 pixels from an unmasked reference pixel */
      nTrimmed = 0;
      border = 2*factor;
      if ((border & 7) != 0 || border == 0)
	border = (border + 8) & ~7;
      rpbw = rw + 2 * border;
      rpbh = rh + 2 * border;
      rpbbpl = (rpbw + 7) >> 3;
      rpbmask = (unsigned char *) malloc(rpbh * rpbbpl);
      if (rpbmask == NULL)
	{
	  sprintf(r.message, "Could not allocate rpbmask array\n");
	  return;
	}
      memset(rpbmask, 0, rpbh * rpbbpl);
      for (y = 0; y < rh; ++y)
	memcpy(&rpbmask[(y+border)*rpbbpl + (border>>3)],
	       &crmask[y*rcmbpl], rcmbpl);
      nBytes = rpbh * rpbbpl;
#if 0
      Log("factor = %d border = %d rw = %d rh = %d rpbw = %d rpbh = %d\n",
	  factor, border, rw, rh, rpbw, rpbh);
      for (y = 0; y < nBytes; ++y)
      	Log("rpbmask[%d] = %d\n", y, rpbmask[y]);
      for (y = 0; y < rpbh; ++y)
	{
	  linePos = 0;
	  for (x = 0; x < rpbw; ++x)
	    {
	      if (rpbmask[y * rpbbpl + (x >> 3)] & (0x80 >> (x & 7)))
		line[linePos++] = '1';
	      else
		line[linePos++] = '0';
	      line[linePos++] = ' ';
	    }
	  line[linePos++] = '\0';
	  Log("rpbmask: %d: %s\n", y, line);
	}
#endif

      dist = (float *) malloc(((size_t) rpbh) * rpbw * sizeof(float));
      if (dist == NULL)
	{
	  sprintf(r.message, "Could not allocate dist array\n");
	  return;
	}

      computeDistance(EUCLIDEAN_DISTANCE_SQUARED, rpbw, rpbh, rpbmask, dist);
      threshold = 4.0 * factor * factor;

#if 0
      Log("threshold = %f\n", threshold);
      line[0] = '\0';
      for (y = 0; y < rpbh; ++y)
	{
	  for (x = 0; x < rpbw; ++x)
	    {
	      linePos = strlen(line);
	      sprintf(&line[linePos], "%f ", IMAGE(dist, rpbw, x, y));
	    }
	  Log("dist: %s\n", line);
	}
#endif

      for (y = 0; y < mph; ++y)
	for (x = 0; x < mpw; ++x)
	  {
	    GETMAP(map, mpw, x, y, &rx00, &ry00, &rc00);
	    if (rc00 == 0.0)
	      continue;
	    xv = factor * rx00 - refox;
	    yv = factor * ry00 - refoy;
	    ixv = ((int) floor(xv)) + border;
	    iyv = ((int) floor(yv)) + border;
	    if (ixv < 0 || ixv >= rpbw ||
		iyv < 0 || iyv >= rpbh ||
		IMAGE(dist, rpbw, ixv, iyv) > threshold)
	      {
		++nTrimmed;
		if (x == 210 && y == 175)
		  Log("SCB 1: %f %f\n",
		      IMAGE(dist, rpbw, ixv, iyv), threshold);
		MAP(map, mpw, x, y).c = 0.0;
	      }
	  }
      free(dist);
      free(rpbmask);
      Log("Trimmed %d map points in stage 2\n", nTrimmed);

      // remove any valid map points that aren't part of a square
      //   of valid map points
      nTrimmed = 0;
      for (y = 0; y < mph; ++y)
	for (x = 0; x < mpw; ++x)
	  {
	    if (MAP(map, mpw, x, y).c == 0.0)
	      continue;
	    valid = 0;
	    for (dy = 0; !valid && dy < 2; ++dy)
	      for (dx = 0; !valid && dx < 2; ++dx)
		{
		  valid = 1;
		  for (ddy = 0; valid && ddy < 2; ++ddy)
		    {
		      iy = y + dy + ddy - 1;
		      if (iy < 0 || iy >= mph)
			{
			  valid = 0;
			  break;
			}
		    for (ddx = 0; ddx < 2; ++ddx)
		      {
			ix = x + dx + ddx - 1;
			if (ix < 0 || ix >= mpw ||
			    MAP(map, mpw, ix, iy).c == 0.0)
			  {
			    valid = 0;		
			    break;
			  }
		      }
		    }
		}
	    if (!valid)
	      {
		++nTrimmed;
		if (x == 210 && y == 175)
		  Log("SCB 2: %d\n", nTrimmed);
		MAP(map, mpw, x, y).c = 0.0;
	      }
	  }
      Log("Trimmed %d map points in stage 3\n", nTrimmed);

      /* compute the constraining map if one was given */
      ncc = 0;
      if (constrainingMapFactor != 0)
	for (y = 0; y < mph; ++y)
	  for (x = 0; x < mpw; ++x)
	    {
	      rx = lFactor * (x + mox);
	      ry = lFactor * (y + moy);
	      x0 = rx / constrainingMapFactor - offxc;
	      y0 = ry / constrainingMapFactor - offyc;
	      ix0 = (int) floor(x0);
	      iy0 = (int) floor(y0);
	      rrx = x0 - ix0;
	      rry = y0 - iy0;
	      if (ix0 < -1 || ix0 >= mpwc ||
		  iy0 < -1 || iy0 >= mphc)
		{
		  SETMAP(mapCons, mpw, x, y, 0.0, 0.0, 0.0);
		  continue;
		}
	      while (ix0 < 0)
		{
		  ++ix0;
		  rrx -= 1.0;
		}
	      while (ix0 >= mpwc-1)
		{
		  --ix0;
		  rrx += 1.0;
		}
	      while (iy0 < 0)
		{
		  ++iy0;
		  rry -= 1.0;
		}
	      while (iy0 >= mphc-1)
		{
		  --iy0;
		  rry += 1.0;
		}
	      GETMAP(constrainingMap, mpwc, ix0, iy0, &rx00, &ry00, &rc00);
	      GETMAP(constrainingMap, mpwc, ix0, iy0+1, &rx01, &ry01, &rc01);
	      GETMAP(constrainingMap, mpwc, ix0+1, iy0, &rx10, &ry10, &rc10);
	      GETMAP(constrainingMap, mpwc, ix0+1, iy0+1, &rx11, &ry11, &rc11);
	      rx = rx00 * (rrx - 1.0) * (rry - 1.0)
		- rx10 * rrx * (rry - 1.0) 
		- rx01 * (rrx - 1.0) * rry
		+ rx11 * rrx * rry;
	      ry = ry00 * (rrx - 1.0) * (rry - 1.0)
		- ry10 * rrx * (rry - 1.0) 
		- ry01 * (rrx - 1.0) * rry
		+ ry11 * rrx * rry;
	      rc = rc00 * (rrx - 1.0) * (rry - 1.0)
		- rc10 * rrx * (rry - 1.0) 
		- rc01 * (rrx - 1.0) * rry
		+ rc11 * rrx * rry;

	      rx = rx * constrainingMapFactor / lFactor;
	      ry = ry * constrainingMapFactor / lFactor;

	      SETMAP(mapCons, mpw, x, y, rx, ry, rc);
	      ++ncc;
	      //	      Log("At level %d: mapCons(%d, %d) = %f %f %f\n",
	      //		  level, x, y, rx, ry, rc);
	    }
      Log("constraining map at this level has %d points\n", ncc);

      /* calculate energy of mapping */
      distortion = 0.0;
      correlation = 0.0;
      correspondence = 0.0;
      energy = 0.0;

      /* contribution from correlation */
      sir = 0.0;
      si2 = 0.0;
      sr2 = 0.0;
      si = 0.0;
      sr = 0.0;
      nPoints = 0;

      for (y = 0; y < ih; ++y)
	for (x = 0; x < iw; ++x)
	  {
#if MASKING
	    if ((cimask[y*icmbpl + (x >> 3)] & (0x80 >> (x & 7))) == 0)
	      continue;
#endif
	    /* use bilinear interpolation to find value */
	    xv = (x + 0.5 + imgox) / factor - mox;
	    yv = (y + 0.5 + imgoy) / factor - moy;
	    ixv = ((int) (xv + 2.0)) - 2;
	    iyv = ((int) (yv + 2.0)) - 2;
	    rrx = xv - ixv;
	    rry = yv - iyv;
	    if (ixv < 0 || ixv >= mpw_minus_1 ||
		iyv < 0 || iyv >= mph_minus_1)
	      continue;

	    GETMAP(map, mpw, ixv, iyv, &rx00, &ry00, &rc00);
	    GETMAP(map, mpw, ixv, iyv+1, &rx01, &ry01, &rc01);
	    GETMAP(map, mpw, ixv+1, iyv, &rx10, &ry10, &rc10);
	    GETMAP(map, mpw, ixv+1, iyv+1, &rx11, &ry11, &rc11);
	    if (rc00 == 0.0 || rc01 == 0.0 || rc10 == 0.0 || rc11 == 0.0)
	      continue;
	    rx = rx00 * (rrx - 1.0) * (rry - 1.0)
	      - rx10 * rrx * (rry - 1.0) 
	      - rx01 * (rrx - 1.0) * rry
	      + rx11 * rrx * rry;
	    ry = ry00 * (rrx - 1.0) * (rry - 1.0)
	      - ry10 * rrx * (rry - 1.0) 
	      - ry01 * (rrx - 1.0) * rry
	      + ry11 * rrx * rry;
	    rx = factor * rx - 0.5;
	    ry = factor * ry - 0.5;
	    irx = ((int) (rx + 1.0)) - 1;
	    iry = ((int) (ry + 1.0)) - 1;
	    rrx = rx - irx;
	    rry = ry - iry;
	    if (irx < 0 || irx >= rw-1 || iry < 0 || iry >= rh-1)
	      continue;

#if MASKING
	    if ((crmask[iry*rcmbpl + (irx >> 3)] & (0x80 >> (irx & 7))) == 0 ||
		rrx > 0.0 && (crmask[iry*rcmbpl + ((irx+1) >> 3)] & (0x80 >> ((irx+1) & 7))) == 0 ||
		rry > 0.0 && (crmask[(iry+1)*rcmbpl + (irx >> 3)] & (0x80 >> (irx & 7))) == 0 ||
		rrx > 0.0 && rry > 0.0 && (crmask[(iry+1)*rcmbpl +((irx+1) >> 3)] & (0x80 >> ((irx+1) & 7))) == 0)
	      continue;
#endif
	    ++nPoints;
	    r00 = IMAGE(cref, rw, irx, iry);
	    r01 = IMAGE(cref, rw, irx, iry + 1);
	    r10 = IMAGE(cref, rw, irx + 1, iry);
	    r11 = IMAGE(cref, rw, irx + 1, iry + 1);

	    rv = r00 * (rrx - 1.0) * (rry - 1.0)
	      - r10 * rrx * (rry - 1.0) 
	      - r01 * (rrx - 1.0) * rry
	      + r11 * rrx * rry;
	    if (isnan(rv))
	      Error("ISNAN internal error\nrx = %f irx = %d ry = %f iry = %d rrx = %f rry = %f\n",
		    rx, irx, ry, iry, rrx, rry);
	    iv = IMAGE(cimage, iw, x, y);
	    si += iv;
	    si2 += iv * iv;
	    sr += rv;
	    sr2 += rv * rv;
	    sir += iv * rv;
	  }
#if MASKING
      requiredPoints = (size_t) ceil(minMaskCount[cLevel] * c.minOverlap * 0.01);
#else
      requiredPoints = ((size_t) iw) * ih * c.minOverlap * 0.01;
#endif
      if (nPoints < requiredPoints)
	{
	  sr += (requiredPoints - nPoints) * 255.0;
	  sr2 += (requiredPoints - nPoints) * 255.0 * 255.0;
	  effectivePoints = requiredPoints;
	}
      else
	effectivePoints = nPoints;
      if (effectivePoints != 0)
	{
	  mi = si / effectivePoints;
	  mr = sr / effectivePoints;
	  denom = (si2 - 2.0 * mi * si + effectivePoints * mi * mi) *
	    (sr2 - 2.0 * mr * sr + effectivePoints * mr * mr);
	  /*      printf("denom = %f mi = %f mr = %f si = %f si2 = %f sr = %f h = %d w = %d\n",
		  denom, mi, mr, si, si2, sr, h, w); */
	  if (denom < 0.001)
	    corr = -1000000.0;
	  else
	    corr = (sir - mi * sr - mr * si + effectivePoints * mi * mr) / sqrt(denom);
	}
      else
	corr = -1000000.0;
      correlation = corr;
      Log("Level %d map has initial correlation energy of %f\n",
	     level, correlation);
      Log("mi %f mr %f si %f sr %f ih %d iw %d rh %d rw %d si2 %f sr2 %f denom %f sir %f r %f\n",
	     mi, mr, si, sr, ih, iw, rh, rw, si2, sr2, denom, sir, r);
      Log("nPoints = %d requiredPoints = %d effectivePoints = %d\n",
	     nPoints, requiredPoints, effectivePoints);

      /* contribution from distortion */
#if FOLDING
      memset(irdisc, 0xff, mph * imbpl);
#endif
      dPoints = 0;
      de = 0.0;
      for (y = 0; y < mph_minus_1; ++y)
	for (x = 0; x < mpw_minus_1; ++x)
	  {
#if FOLDING
	    if ((idisc[y*imbpl + (x >> 3)] & (0x80 >> (x & 7))) == 0)
	      continue;
#endif

	    GETMAP(map, mpw, x, y, &x0, &y0, &c0);
	    GETMAP(map, mpw, x, y+1, &x1, &y1, &c1);
	    GETMAP(map, mpw, x+1, y+1, &x2, &y2, &c2);
	    GETMAP(map, mpw, x+1, y, &x3, &y3, &c3);

	    if (c0 == 0.0 || c1 == 0.0 || c2 == 0.0 || c3 == 0.0)
	      continue;

#if FOLDING
	    rx = 0.25 * (x0 + x1 + x2 + x3);
	    ry = 0.25 * (y0 + y1 + y2 + y3);
	    ix = floor(rx);
	    iy = floor(ry);
	    if (ix >= 0 && ix < refw && iy >= 0 && iy < refh &&
		(rdisc[iy*rmbpl + (ix >> 3)] & (0x80 >> (ix & 7))) == 0)
	      {
		irdisc[y*imbpl + (x >> 3)] &= ~(0x80 >> (x & 7));
		continue;
	      }
#endif

	    area = - 0.5 *((x2 - x0) * (y3 - y1) - (x3 - x1) * (y2 - y0));
	    l0 = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
	    l1 = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
	    l2 = sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2));
	    l3 = sqrt((x0 - x3) * (x0 - x3) + (y0 - y3) * (y0 - y3));
	    
	    de += (area - 1.0) * (area - 1.0) +
	      (l0 - 1.0) * (l0 - 1.0) +
	      (l1 - 1.0) * (l1 - 1.0) +
	      (l2 - 1.0) * (l2 - 1.0) +
	      (l3 - 1.0) * (l3 - 1.0);
	    ++dPoints;
	  }
      if (dPoints > 0)
	distortion = de / dPoints;
      else
	distortion = 1000000.0;

      /* contribution from correspondence points */
      for (i = 0; i < nCpts; ++i)
	{
	  xv = cpts[i].ix / lFactor - mox;
	  yv = cpts[i].iy / lFactor - moy;
	  ixv = ((int) (xv + 2.0)) - 2;
	  iyv = ((int) (yv + 2.0)) - 2;
	  rrx = xv - ixv;
	  rry = yv - iyv;
	  if (ixv < 0 || ixv >= mpw-1 || iyv < 0 || iyv >= mph-1)
	    {
	      cpts[i].energy = 1000000.0;
	      correspondence += cpts[i].energy;
	      continue;
	    }
	  GETMAP(map, mpw, ixv, iyv, &rx00, &ry00, &rc00);
	  GETMAP(map, mpw, ixv, iyv+1, &rx01, &ry01, &rc01);
	  GETMAP(map, mpw, ixv+1, iyv, &rx10, &ry10, &rc10);
	  GETMAP(map, mpw, ixv+1, iyv+1, &rx11, &ry11, &rc11);
	  rx = rx00 * (rrx - 1.0) * (rry - 1.0)
	    - rx10 * rrx * (rry - 1.0) 
	    - rx01 * (rrx - 1.0) * rry
	    + rx11 * rrx * rry;
	  ry = ry00 * (rrx - 1.0) * (rry - 1.0)
	    - ry10 * rrx * (rry - 1.0) 
	    - ry01 * (rrx - 1.0) * rry
	    + ry11 * rrx * rry;
	  rx = lFactor * rx;
	  ry = lFactor * ry;
	  rrx = cpts[i].rx;
	  rry = cpts[i].ry;
	  distance = kFactor * (hypot(rrx-rx, rry-ry) - c.correspondenceThreshold);
	  Log("CPTS %d: ix=%f iy=%f rx=%f ry=%f rrx=%f rry=%f kF=%f dist=%f\n",
	      i, cpts[i].ix, cpts[i].iy, rx, ry, rrx, rry, kFactor, distance);
	  if (distance < 0.0)
	    cpts[i].energy = 0.0;
	  else
	    cpts[i].energy = distance;
	  correspondence += cpts[i].energy;
	}

      /* contribution from constraining map */
      constraining = 0.0;
      cth = c.constrainingThreshold / lFactor;
      scc = 0.0;
      if (constrainingMapFactor != 0)
	{
	  for (y = 0; y < mph; ++y)
	    for (x = 0; x < mpw; ++x)
	      {
		e = &MAP(map, mpw, x, y);
		if (e->c == 0.0)
		  continue;

		ec = &MAP(mapCons, mpw, x, y);
		if (ec->c == 0.0)
		  continue;
		distance = hypot(e->x - ec->x, e->y - ec->y) - cth;
		if (distance < 0.0)
		  distance = 0.0;
		scc += ec->c;
		constraining += ec->c * distance;
	      }
	  constraining = constraining / mpw / mph;
	}
      Log("init contrib from constraining map = %f %f %f\n", constraining, cth,
	  scc);

      energy = distortion * c.distortion - correlation + correspondence * c.correspondence +
	constraining * c.constraining;
      Log("Total energy of level %d map is %f (dist %f - correl %f + corresp %f + constr %f)\n",
	  level, energy, distortion * c.distortion, correlation,
	  correspondence * c.correspondence,
	  constraining * c.constraining);

      mSize = ((size_t) mph) * mpw;
      prop = (MapElement*) malloc(mSize * sizeof(MapElement));
      if (prop == NULL)
	{
	  sprintf(r.message, "Could not allocate prop array\n");
	  return;
	}
      memset(prop, 0, mSize * sizeof(MapElement));
      
      fflush(stdout);

      multiplier = 1 << level;
      goalCoverage = (c.quality * mpw) * mph * multiplier;
      actualCoverage = 0.0;

      moveCount = 0;
      acceptedMoveCount = 0;
      logMaxRadius = log(0.5);
      logMinRadius = log(0.02 / factor); /* no use going smaller than 2% of the pixel size */
      logRadiusRange = logMaxRadius - logMinRadius;
      udLimit = 0.5 * 0.5;
      diagLimit = 0.5;   // (sqrt(2.0) / 2.0) * (sqrt(2.0) / 2.0)

      memset(statLogRadius, 0, 21*sizeof(size_t));
      memset(statTheta, 0, 21*sizeof(size_t));
      memset(statDeltaE, 0, 21*sizeof(double));

      while (actualCoverage < goalCoverage)
	{
#if GRAPHICS
	  if (moveCount % (mpw * mph) == 0)
	    displayLevel = level;
#endif
#if 1
	  // go in sequential order for cache-friendliness
	  icx = moveCount % mpw; 
	  icy = (moveCount / mpw) % mph;
#endif
#if 0
	  icx = (int) floor(drand48() * mpw);
	  icy = (int) floor(drand48() * mph);
#endif
	  ++moveCount;

	  // make sure (icx,icy) is in the trimmed subregion
	  GETMAP(map, mpw, icx, icy, &cx, &cy, &cc);
	  if (cc == 0.0)
	    continue;
	  if (isnan(cx) || isnan(cy))
	    Error("ISNAN cx %f cy %f\n", cx, cy);

	  actualCoverage += 1.0;
	  rnd = drand48();
	  logRadius = rnd * logRadiusRange + logMinRadius;
	  radius = exp(logRadius);
	  theta = drand48() * 2.0 * M_PI;
	  cx += radius * cos(theta);
	  cy += radius * sin(theta);
#if DEBUG_MOVES
#if 0
	  Log("Testing move %d: cx = %f cy = %f radius = %f theta = %f\n",
 	      moveCount, cx, cy, radius, theta);
#endif
#endif

	  // check that the move will not distort the grid too much;
	  //  we check that the distance to each neighboring point
	  //  is not less than half the nominal distance
	  if (icx >= 1)
	    {
	      if (icy >= 1)
		{
		  e = &MAP(map, mpw, icx-1, icy-1);
		  if (e->c != 0.0)
		    {
		      dx = e->x - cx;
		      dy = e->y - cy;
		      d2 = dx*dx + dy*dy;
		      if (d2 < diagLimit)
			continue;
		    }
		}
	      if (icy < mph_minus_1)
		{
		  e = &MAP(map, mpw, icx-1, icy+1);
		  if (e->c != 0.0)
		    {
		      dx = e->x - cx;
		      dy = e->y - cy;
		      d2 = dx*dx + dy*dy;
		      if (d2 < diagLimit)
			continue;
		    }
		}
	      e = &MAP(map, mpw, icx-1, icy);
	      if (e->c != 0.0)
		{
		  dx = e->x - cx;
		  dy = e->y - cy;
		  d2 = dx*dx + dy*dy;
		  if (d2 < udLimit)
		    continue;
		}
	    }
	  if (icx < mpw_minus_1)
	    {
	      if (icy >= 1)
		{
		  e = &MAP(map, mpw, icx+1, icy-1);
		  if (e->c != 0.0)
		    {
		      dx = e->x - cx;
		      dy = e->y - cy;
		      d2 = dx*dx + dy*dy;
		      if (d2 < diagLimit)
			continue;
		    }
		}
	      if (icy < mph_minus_1)
		{
		  e = &MAP(map, mpw, icx+1, icy+1);
		  if (e->c != 0.0)
		    {
		      dx = e->x - cx;
		      dy = e->y - cy;
		      d2 = dx*dx + dy*dy;
		      if (d2 < diagLimit)
			continue;
		    }
		}
	      e = &MAP(map, mpw, icx+1, icy);
	      if (e->c != 0.0)
		{
		  dx = e->x - cx;
		  dy = e->y - cy;
		  d2 = dx*dx + dy*dy;
		  if (d2 < udLimit)
		    continue;
		}
	    }
	  if (icy >= 1)
	    {
	      e = &MAP(map, mpw, icx, icy-1);
	      if (e->c != 0.0)
		{
		  dx = e->x - cx;
		  dy = e->y - cy;
		  d2 = dx*dx + dy*dy;
		  if (d2 < udLimit)
		    continue;
		}
	    }
	  if (icy < mph_minus_1)
	    {
	      e = &MAP(map, mpw, icx, icy+1);
	      if (e->c != 0.0)
		{
		  dx = e->x - cx;
		  dy = e->y - cy;
		  d2 = dx*dx + dy*dy;
		  if (d2 < udLimit)
		    continue;
		}
	    }
	  
	  SETMAP(prop, mpw, icx, icy, cx, cy, 1.0);
	  
	  changeMinX = icx;
	  changeMaxX = icx;
	  changeMinY = icy;
	  changeMaxY = icy;

	  /* evaluate effect of that move on energy */

	  /* add in contribution from correlation */
	  dsir = 0.0;
	  dsr2 = 0.0;
	  dsr = 0.0;
	  dsi = 0.0;
	  dsi2 = 0.0;
	  cPoints = 0;
	  sx = (changeMinX - 1 + mox) * factor - imgox;
	  if (sx < 0)
	    sx = 0;
	  ex = (changeMaxX + 1 + mox) * factor - 1 - imgox;
	  if (ex >= iw)
	    ex = iw - 1;
	  sy = (changeMinY - 1 + moy) * factor - imgoy;
	  if (sy < 0)
	    sy = 0;
	  ey = (changeMaxY + 1 + moy) * factor - 1 - imgoy;
	  if (ey >= ih)
	    ey = ih - 1;

	  //	  printf("sx = %d ex = %d sy = %d ey = %d\n",
	  //		 sx, ex, sy, ey);
	  for (y = sy; y <= ey; ++y)
            for (x = sx; x <= ex; ++x)
	      {
#if MASKING
		if ((cimask[y*icmbpl + (x >> 3)] & (0x80 >> (x & 7))) == 0)
		  continue;
#endif
		xv = (x + 0.5 + imgox) / factor - mox;
		yv = (y + 0.5 + imgoy) / factor - moy;
		ixv = ((int) (xv + 2.0)) - 2;
		iyv = ((int) (yv + 2.0)) - 2;
		rrx = xv - ixv;
		rry = yv - iyv;

		if (ixv < 0 || ixv >= mpw_minus_1 ||
		    iyv < 0 || iyv >= mph_minus_1)
		  continue;

		if (MAP(prop, mpw, ixv, iyv).c == 0.0 &&
		    MAP(prop, mpw, ixv+1, iyv).c == 0.0 &&
		    MAP(prop, mpw, ixv, iyv+1).c == 0.0 &&
		    MAP(prop, mpw, ixv+1, iyv+1).c == 0.0)
		    continue;

		/* use bilinear interpolation to find value */
		iv = IMAGE(cimage, iw, x, y);
		if (iv < 0.0)
		  Error("Internal error: iv out of range: %f\n", iv);
		GETMAP(map, mpw, ixv, iyv, &rx00, &ry00, &rc00);
		GETMAP(map, mpw, ixv, iyv+1, &rx01, &ry01, &rc01);
		GETMAP(map, mpw, ixv+1, iyv, &rx10, &ry10, &rc10);
		GETMAP(map, mpw, ixv+1, iyv+1, &rx11, &ry11, &rc11);
		if (rc00 != 0.0 && rc01 != 0.0 && rc10 != 0.0 && rc11 != 0.0)
		  {
		    rx = rx00 * (rrx - 1.0) * (rry - 1.0)
		      - rx10 * rrx * (rry - 1.0) 
		      - rx01 * (rrx - 1.0) * rry
		      + rx11 * rrx * rry;
		    ry = ry00 * (rrx - 1.0) * (rry - 1.0)
		      - ry10 * rrx * (rry - 1.0) 
		      - ry01 * (rrx - 1.0) * rry
		      + ry11 * rrx * rry;
		    rx = factor * rx - 0.5;
		    ry = factor * ry - 0.5;
		    irx = ((int) (rx + 1.0)) - 1;
		    iry = ((int) (ry + 1.0)) - 1;
		    rrx = rx - irx;
		    rry = ry - iry;
		    if (irx >= 0 && irx < rw-1 && iry >= 0 && iry < rh-1)
		      {
#if MASKING
			if ((crmask[iry*rcmbpl + (irx >> 3)] & (0x80 >> (irx & 7))) != 0 &&
			    ((rrx <= 0.0) || (crmask[iry*rcmbpl + ((irx+1) >> 3)] & (0x80 >> ((irx+1) & 7))) != 0) &&
			    ((rry <= 0.0) || (crmask[(iry+1)*rcmbpl + (irx >> 3)] & (0x80 >> (irx & 7))) != 0) &&
			    ((rrx <= 0.0) || (rry <= 0.0) || (crmask[(iry+1)*rcmbpl +((irx+1) >> 3)] & (0x80 >> ((irx+1) & 7))) != 0))
			  {
#endif
			    r00 = IMAGE(cref, rw, irx, iry);
			    r01 = IMAGE(cref, rw, irx, iry + 1);
			    r10 = IMAGE(cref, rw, irx + 1, iry);
			    r11 = IMAGE(cref, rw, irx + 1, iry + 1);

			    rv = r00 * (rrx - 1.0) * (rry - 1.0)
			      - r10 * rrx * (rry - 1.0) 
			      - r01 * (rrx - 1.0) * rry
			      + r11 * rrx * rry;
			    dsi -= iv;
			    dsi2 -= iv * iv;
			    dsir -= iv * rv;
			    dsr2 -= rv * rv;
			    dsr -= rv;
			    cPoints -= 1;
#if MASKING
			  }
#endif
		      }
		  }

		/* use bilinear interpolation to find value */
		mp = (MAP(prop, mpw, ixv, iyv).c != 0.0) ? prop : map;
		GETMAP(mp, mpw, ixv, iyv, &rx00, &ry00, &rc00);
		mp = (MAP(prop, mpw, ixv, iyv+1).c != 0.0) ? prop : map;
		GETMAP(mp, mpw, ixv, iyv+1, &rx01, &ry01, &rc01);
		mp = (MAP(prop, mpw, ixv+1, iyv).c != 0.0) ? prop : map;
		GETMAP(mp, mpw, ixv+1, iyv, &rx10, &ry10, &rc10);
		mp = (MAP(prop, mpw, ixv+1, iyv+1).c != 0.0) ? prop : map;
		GETMAP(mp, mpw, ixv+1, iyv+1, &rx11, &ry11, &rc11);
		if (isnan(ry11))
		  Error("ISNAN ry11 %d %d %d %d %f %d\n",
			ixv, iyv, mpw, mph, MAP(prop, mpw, ixv+1, iyv+1).c, level);
		if (rc00 == 0.0 || rc01 == 0.0 || rc10 == 0.0 || rc11 == 0.0)
		  continue;

		rrx = xv - ixv;
		rry = yv - iyv;
		rx = rx00 * (rrx - 1.0) * (rry - 1.0)
		  - rx10 * rrx * (rry - 1.0) 
		  - rx01 * (rrx - 1.0) * rry
		  + rx11 * rrx * rry;
		ry = ry00 * (rrx - 1.0) * (rry - 1.0)
		  - ry10 * rrx * (rry - 1.0) 
		  - ry01 * (rrx - 1.0) * rry
		  + ry11 * rrx * rry;
		if (isnan(ry))
		  Error("Internal error: RY ISNAN %f %f %f %f %f %f\n%d %d %d %d\n",
			ry00, ry01, ry10, ry11, rrx, rry,
			ixv, iyv, mpw, mph);
		rx = factor * rx - 0.5;
		ry = factor * ry - 0.5;

		irx = ((int) (rx + 1.0)) - 1;
		iry = ((int) (ry + 1.0)) - 1;
		if (irx < 0 || irx >= rw-1 || iry < 0 || iry >= rh-1)
		  continue;
		rrx = rx - irx;
		rry = ry - iry;
		//		printf("irx = %d iry = %d rw = %d rh = %d factor = %d ixv = %d iyv = %d\n",
		//		       irx, iry, rw, rh, factor, ixv, iyv);
#if MASKING
		if ((crmask[iry*rcmbpl + (irx >> 3)] & (0x80 >> (irx & 7))) == 0 ||
		    rrx > 0.0 && (crmask[iry*rcmbpl + ((irx+1) >> 3)] & (0x80 >> ((irx+1) & 7))) == 0 ||
		    rry > 0.0 && (crmask[(iry+1)*rcmbpl + (irx >> 3)] & (0x80 >> (irx & 7))) == 0 ||
		    rrx > 0.0 && rry > 0.0 && (crmask[(iry+1)*rcmbpl +((irx+1) >> 3)] & (0x80 >> ((irx+1) & 7))) == 0)
		  continue;
#endif
		r00 = IMAGE(cref, rw, irx, iry);
		r01 = IMAGE(cref, rw, irx, iry + 1);
		r10 = IMAGE(cref, rw, irx + 1, iry);
		r11 = IMAGE(cref, rw, irx + 1, iry + 1);
		    
		rv = r00 * (rrx - 1.0) * (rry - 1.0)
		  - r10 * rrx * (rry - 1.0) 
		  - r01 * (rrx - 1.0) * rry
		  + r11 * rrx * rry;
		dsi += iv;
		dsi2 += iv * iv;
		dsir += iv * rv;
		dsr2 += rv * rv;
		dsr += rv;
		cPoints += 1;
	      }

	  if (nPoints < requiredPoints)
	    {
	      dsr -= (requiredPoints - nPoints) * 255.0;
	      dsr2 -= (requiredPoints - nPoints) * 255.0 * 255.0;
	    }
	  newPoints = nPoints + cPoints;
	  if (newPoints < requiredPoints)
	    {
	      dsr += (requiredPoints - newPoints) * 255.0;
	      dsr2 += (requiredPoints - newPoints) * 255.0 * 255.0;
	      newEffectivePoints = requiredPoints;
	    }
	  else
	    newEffectivePoints = newPoints;
	  newsi = si + dsi;
	  if (newsi < -0.1)
	    Error("Internal error: negative newsi\nsi = %f  dsi = %f  newsi = %f\n",
		  si, dsi, newsi);
	  newsi2 = si2 + dsi2;
	  newsir = sir + dsir;
	  newsr = sr + dsr;
	  newsr2 = sr2 + dsr2;
	  if (newEffectivePoints != 0)
	    {
	      newmi = newsi / newEffectivePoints;
	      newmr = newsr / newEffectivePoints;

	      newDenom = (newsi2 - 2.0 * newmi * newsi + newEffectivePoints * newmi * newmi) *
		(newsr2 - 2.0 * newmr * newsr + newEffectivePoints * newmr * newmr);
	      if (newDenom < 0.001)
		newCorr = -1000000.0;
	      else
		newCorr = (newsir - newmi * newsr - newmr * newsi + newEffectivePoints * newmi * newmr) / sqrt(newDenom);
	    }
	  else
	    newCorr = -1000000.0;
	  newCorrelation = newCorr;
	  if (newCorrelation > 1.1)
	    {
	      Log("mi %f mr %f si %f sr %f ih %d iw %d si2 %f sr2 %f denom %f sir %f r %f\n",
		  newmi, newmr, newsi, newsr, ih, iw, newsi2, newsr2, newDenom,
		  newsir, newCorr);
	      Log("nPoints = %d requiredPoints = %d effectivePoints = %d\n",
		  newPoints, requiredPoints, newEffectivePoints);
	      Error("Internal error: Level %d map has new correlation energy of %f\n",
		    level, newCorrelation);
	    }

	  /* add in contribution from distortion */
	  de = 0.0;
	  newDPoints = dPoints;
          for (y = changeMinY - 1; y <= changeMaxY; ++y)
	    {
              if (y < 0 || y >= mph_minus_1)
		continue;
              for (x = changeMinX - 1; x <= changeMaxX; ++x)
      	        {
	          if (x < 0 || x >= mpw_minus_1)
	            continue;

#if FOLDING
		  nirdisc[y*imbpl + (x >> 3)] &= ~(0x80 >> (x & 7));
		  nirdisc[y*imbpl + (x >> 3)] |= irdisc[y*imbpl + (x >> 3)] & (0x80 >> (x & 7));
		  if ((idisc[y*imbpl + (x >> 3)] & (0x80 >> (x & 7))) == 0)
		    continue;
#endif

		  upd0 = MAP(prop, mpw, x, y).c != 0.0;
		  upd1 = MAP(prop, mpw, x, y+1).c != 0.0;
		  upd2 = MAP(prop, mpw, x+1, y+1).c != 0.0;
		  upd3 = MAP(prop, mpw, x+1, y).c != 0.0;
		      
		  if (!upd0 && !upd1 && !upd2 && !upd3)
		    continue;
#if FOLDING
		  if ((idisc[y*imbpl + (x >> 3)] & (0x80 >> (x & 7))) == 0)
		    continue;
		  if ((irdisc[y*imbpl + (x >> 3)] & (0x80 >> (x & 7))) == 0)
		    continue;
#endif

		  GETMAP(map, mpw, x, y, &x0, &y0, &c0);
		  GETMAP(map, mpw, x, y+1, &x1, &y1, &c1);
		  GETMAP(map, mpw, x+1, y+1, &x2, &y2, &c2);
		  GETMAP(map, mpw, x+1, y, &x3, &y3, &c3);

		  if (c0 == 0.0 || c1 == 0.0 || c2 == 0.0 || c3 == 0.0)
		    continue;

		  nx0 = x0; ny0 = y0; nc0 = c0;
		  nx1 = x1; ny1 = y1; nc1 = c1;
		  nx2 = x2; ny2 = y2; nc2 = c2;
		  nx3 = x3; ny3 = y3; nc3 = c3;

		  if (upd0)
		    GETMAP(prop, mpw, x, y, &nx0, &ny0, &nc0);
		  if (upd1)
		    GETMAP(prop, mpw, x, y+1, &nx1, &ny1, &nc1);
		  if (upd2)
		    GETMAP(prop, mpw, x+1, y+1, &nx2, &ny2, &nc2);
		  if (upd3)
		    GETMAP(prop, mpw, x+1, y, &nx3, &ny3, &nc3);

#if FOLDING
		  rx = 0.25 * (nx0 + nx1 + nx2 + nx3);
		  ry = 0.25 * (ny0 + ny1 + ny2 + ny3);
		  ix = floor(rx);
		  iy = floor(ry);
		  if (ix >= 0 && ix < refw && iy >= 0 && iy < refh &&
		      (rdisc[iy*rmbpl + (ix >> 3)] & (0x80 >> (ix & 7))) == 0)
		    {
		      nirdisc[y*imbpl + (x >> 3)] &= ~(0x80 >> (x & 7));
		      continue;
		    }
#endif

		  area = - 0.5 *((x2 - x0) * (y3 - y1) -
				 (x3 - x1) * (y2 - y0));
		  l0 = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
		  l1 = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
		  l2 = sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2));
		  l3 = sqrt((x0 - x3) * (x0 - x3) + (y0 - y3) * (y0 - y3));
		  narea = - 0.5 *((nx2 - nx0) * (ny3 - ny1) -
				  (nx3 - nx1) * (ny2 - ny0));
		  nl0 = sqrt((nx1 - nx0) * (nx1 - nx0) +
			     (ny1 - ny0) * (ny1 - ny0));
		  nl1 = sqrt((nx2 - nx1) * (nx2 - nx1) +
			     (ny2 - ny1) * (ny2 - ny1));
		  nl2 = sqrt((nx3 - nx2) * (nx3 - nx2) +
			     (ny3 - ny2) * (ny3 - ny2));
		  nl3 = sqrt((nx0 - nx3) * (nx0 - nx3) +
			     (ny0 - ny3) * (ny0 - ny3));
		  ode = (area - 1.0) * (area - 1.0) +
		    (l0 - 1.0) * (l0 - 1.0) +
		    (l1 - 1.0) * (l1 - 1.0) +
		    (l2 - 1.0) * (l2 - 1.0) +
		    (l3 - 1.0) * (l3 - 1.0);
		  nde = (narea - 1.0) * (narea - 1.0) +
		    (nl0 - 1.0) * (nl0 - 1.0) +
		    (nl1 - 1.0) * (nl1 - 1.0) +
		    (nl2 - 1.0) * (nl2 - 1.0) +
		    (nl3 - 1.0) * (nl3 - 1.0);
		  de += nde - ode;
		  ++updcnt;
		}
	    }
	  // FIX to divide by adjusted dpoints
	  newDistortion = distortion + de / dPoints;

	  /* add in contribution from correspondence energy */
	  ce = 0.0;
	  for (i = 0; i < nCpts; ++i)
	    {
	      xv = cpts[i].ix / lFactor - mox;
	      yv = cpts[i].iy / lFactor - moy;
	      ixv = ((int) (xv + 2.0)) - 2;
	      iyv = ((int) (yv + 2.0)) - 2;
	      if (ixv < changeMinX-1 || ixv > changeMaxX ||
		  iyv < changeMinY-1 || iyv > changeMaxY)
		{
		  cpts[i].newEnergy = cpts[i].energy;
		  continue;
		}
	      rrx = xv - ixv;
	      rry = yv - iyv;
	      if (ixv < 0 || ixv >= mpw_minus_1 || iyv < 0 || iyv >= mph_minus_1)
		{
		  cpts[i].newEnergy = 1000000.0;
		  continue;
		}
	      mp = (MAP(prop, mpw, ixv, iyv).c != 0.0) ? prop : map;
	      GETMAP(mp, mpw, ixv, iyv, &rx00, &ry00, &rc00);
	      mp = (MAP(prop, mpw, ixv, iyv+1).c != 0.0) ? prop : map;
	      GETMAP(mp, mpw, ixv, iyv+1, &rx01, &ry01, &rc01);
	      mp = (MAP(prop, mpw, ixv+1, iyv).c != 0.0) ? prop : map;
	      GETMAP(mp, mpw, ixv+1, iyv, &rx10, &ry10, &rc10);
	      mp = (MAP(prop, mpw, ixv+1, iyv+1).c != 0.0) ? prop : map;
	      GETMAP(mp, mpw, ixv+1, iyv+1, &rx11, &ry11, &rc11);

	      rx = rx00 * (rrx - 1.0) * (rry - 1.0)
		- rx10 * rrx * (rry - 1.0) 
		- rx01 * (rrx - 1.0) * rry
		+ rx11 * rrx * rry;
	      ry = ry00 * (rrx - 1.0) * (rry - 1.0)
		- ry10 * rrx * (rry - 1.0) 
		- ry01 * (rrx - 1.0) * rry
		+ ry11 * rrx * rry;
	      rx = lFactor * rx;
	      ry = lFactor * ry;
	      rrx = cpts[i].rx;
	      rry = cpts[i].ry;
	      distance = kFactor * (hypot(rrx-rx, rry-ry) - c.correspondenceThreshold);
	      if (distance < 0.0)
		cpts[i].newEnergy = 0.0;
	      else
		cpts[i].newEnergy = distance;
	      ce += cpts[i].newEnergy - cpts[i].energy;
	    }
	  newCorrespondence = correspondence + ce;

	  /* add in contribution from constraining map */
	  de = 0.0;
	  if (constrainingMapFactor != 0)
	    {
	      cth = c.constrainingThreshold / lFactor;
	      for (y = changeMinY; y <= changeMaxY; ++y)
		{
		  if (y < 0 || y >= mph)
		    continue;
		  for (x = changeMinX; x <= changeMaxX; ++x)
		    {
		      if (x < -1 || x >= mpw)
			continue;

		      if (MAP(prop, mpw, x, y).c != 0.0)
			{
			  e = &MAP(map, mpw, x, y);
			  if (e->c == 0.0)
			    continue;
			  ec = &MAP(mapCons, mpw, x, y);
			  if (ec->c == 0.0)
			    continue;
			  ep = &MAP(prop, mpw, x, y);
			  distance = hypot(e->x - ec->x, e->y - ec->y) - cth;
			  if (distance < 0.0)
			    distance = 0.0;
			  oldEnergy = ec->c * distance;
			  distance = hypot(ep->x - ec->x, ep->y - ec->y) - cth;
			  if (distance < 0.0)
			    distance = 0.0;
			  de += ec->c * distance - oldEnergy;
			}
		    }
		}
	    }
	  newConstraining = constraining + de / mph / mpw;

	  newEnergy = newDistortion * c.distortion - newCorrelation +
	    newCorrespondence * c.correspondence + newConstraining * c.constraining;

	  /* accept or reject move */
	  accept = (newEnergy < energy);
	  for (y = changeMinY; y <= changeMaxY; ++y)
	    for (x = changeMinX; x <= changeMaxX; ++x)
	      if (MAP(prop, mpw, x, y).c != 0.0)
		{
		  if (x == 210 && y == 175)
		    Log("SCB 3: %d\n", accept);
		  MAP(prop, mpw, x, y).c = 0;
		  if (accept)
		    {
		      /* make the provisional state the new state */
		      GETMAP(prop, mpw, x, y, &rx00, &ry00, &rc00);
		      rc00 = 1.0;
		      SETMAP(map, mpw, x, y, rx00, ry00, rc00);
		    }
		}
	  if (accept)
	    {
#if DEBUG_MOVES
	      Log("Accepted move %d since deltaE = %f (dist %f - correl %f + corres %f + constr %f)\n",
		  moveCount, newEnergy - energy,
		  (newDistortion - distortion) * c.distortion,
		  newCorrelation - correlation,
		  (newCorrespondence - correspondence) * c.correspondence,
		  (newConstraining - constraining) * c.constraining);
	      Log("  cx = %f cy = %f radius = %f offset = %f theta = %f\n",
		  cx, cy, radius, offset, theta);
	      Log("  old energy = %f  estimated new energy = %f (dist %f - correl %f + corres %f + constr %f)\n",
		  energy,
		  newEnergy, newDistortion * c.distortion, newCorrelation,
		  newCorrespondence * c.correspondence,
		  newConstraining * c.constraining);
#endif
	      ++statLogRadius[(int) (20.0 * rnd)];
	      ++statTheta[(int) (20.0 * theta / (2.0 * M_PI))];
	      statDeltaE[(int) (20.0 * rnd)] += energy - newEnergy;

#if FOLDING
	      for (y = changeMinY - 1; y <= changeMaxY; ++y)
		{
		  FIX;
		  for (x = changeMinX - 1; x <= changeMaxX; ++x)
		    {
		      irdisc[y*imbpl + (x >> 3)] &= ~(0x80 >> (x & 7));
		      irdisc[y*imbpl + (x >> 3)] |= nirdisc[y*imbpl + (x >> 3)] & (0x80 >> (x & 7));
		    }
		}
#endif	      

	      for (i = 0; i < nCpts; ++i)
		cpts[i].energy = cpts[i].newEnergy;

	      nPoints = newPoints;
	      si = newsi;
	      si2 = newsi2;
	      sr = newsr;
	      sr2 = newsr2;
	      sir = newsir;
	      correlation = newCorrelation;
	      distortion = newDistortion;
	      correspondence = newCorrespondence;
	      energy = newEnergy;
	      ++acceptedMoveCount;
	      displayLevel = level;
	      //	      sleep(1);
	    }
	}

      Log("Level %d: after %d moves, and %d accepted moves...\n",
	  level, moveCount, acceptedMoveCount);
      Log("          energy is %f\n", energy);

      Log("Accepted shift move statistics:\n");
      Log("   range:  logRadius   theta  deltaE\n");
      for (i = 0; i < 20; ++i)
	Log(" %d-%d%%:  %10.5f  %10.5f  %12.8f  %10.6f\n",
	    5*i, 5*i+5,
	    (100.0 * statLogRadius[i]) / moveCount,
	    (100.0 * statTheta[i]) / moveCount,
	    statLogRadius[i] != 0 ?
	    statDeltaE[i] / statLogRadius[i] : 0.0,
	    statDeltaE[i]);
      Log("--- done with level %d ---\n", level);
      

      free(prop);
#if FOLDING
      free(cirdisc);
#endif

      if (level == c.outputLevel || c.writeAllMaps)
	{
	  warpedArray = (float *) malloc(ih * iw * sizeof(float));
	  validArray = (unsigned char *) malloc(iw * ih *
						sizeof(unsigned char));
	  ComputeWarpedImage(warpedArray, validArray,
			     iw, ih,
			     imgox, imgoy,
			     cref, crmask,
			     rw, rh,
			     map,
			     factor, mpw, mph, mox, moy);

	  correlationArray = (float *) malloc(ih * iw * sizeof(float));
	  ComputeCorrelation(correlationArray,
			     cimage, warpedArray, validArray,
			     iw, ih,
			     c.correlationHalfWidth);

	  // use the correlation as the confidence values
	  //   for the map entries
	  for (y = 0; y < mph; ++y)
	    for (x = 0; x < mpw; ++x)
	      {
		ixv = factor * (x + mox) - imgox;
		iyv = factor * (y + moy) - imgoy;
		MAP(map, mpw, x, y).c = correlationArray[iyv * iw + ixv];
	      }

	  WriteOutputMap(outputName, level, map, mpw, mph, mox, moy);
	  if (c.outputImagesName[0] != '\0')
	    WriteOutputImage(c.outputImagesName, c.outputLevel, warpedArray,
			     iw, ih, 1.0);
	  if (c.outputCorrelationName[0] != '\0')
	    WriteOutputImage(c.outputCorrelationName, c.outputLevel, correlationArray,
			     iw, ih, 256.0);

	  free(warpedArray);
	  free(validArray);
	  free(correlationArray);
	}
    }

  // write out the score for this mapping
  sprintf(fn, "%s.score", outputName);
  f = fopen(fn, "w");
  fprintf(f, "%f %f %f %f %f\n",
	  correlation
	  - c.distortion * distortion
	  - c.correspondence * correspondence
	  - c.constraining * constraining,
	  correlation, distortion, correspondence, constraining);
  fclose(f);

  memcpy(&(r.pair), &(t.pair), sizeof(Pair));
  r.updated = 1;
  r.correlation = correlation;
  r.distortion = distortion;
  r.correspondence = correspondence;
  r.constraining = constraining;
  r.message[0] = '\0';
}


void
WriteOutputMap (char *outputName, int level, MapElement *map,
		int mpw, int mph, int mox, int moy)
{
  char fn[PATH_MAX];
  FILE *f;
  float rx, ry, rc;
  int irx, iry;
  int x, y;
  int minX, maxX, minY, maxY;
  int nx, ny;
  MapElement *outMap;
  
  /* determine how much of the map has valid points in it */
  for (minX = 0; minX < mpw; ++minX)
    {
      for (y = 0; y < mph; ++y)
	if (MAP(map, mpw, minX, y).c != 0.0)
	  break;
      if (y < mph)
	break;
    }
  if (minX >= mpw)
    {
      Log("WriteOutputMap: no valid map points are present.\n");
      return;
    }
  for (maxX = mpw-1; maxX >= minX; --maxX)
    {
      for (y = 0; y < mph; ++y)
	if (MAP(map, mpw, maxX, y).c != 0.0)
	  break;
      if (y < mph)
	break;
    }
  for (minY = 0; minY < mph; ++minY)
    {
      for (x = 0; x < mpw; ++x)
	if (MAP(map, mpw, x, minY).c != 0.0)
	  break;
      if (x < mpw)
	break;
    }
  for (maxY = mph-1; maxY >= minY; --maxY)
    {
      for (x = 0; x < mpw; ++x)
	if (MAP(map, mpw, x, maxY).c != 0.0)
	  break;
      if (x < mpw)
	break;
    }
  nx = maxX - minX + 1;
  ny = maxY - minY + 1;

  /* write out resultant map */
  if (level > c.outputLevel)
    sprintf(fn, "%s.%0.2d.map", outputName, level);
  else
    sprintf(fn, "%s.map", outputName);

  f = fopen(fn, "w");
  if (f == NULL)
    Error("Could not open map file %s for writing\n", fn);
  fprintf(f, "M1\n");
  fprintf(f, "%d\n", level);
  fprintf(f, "%d %d\n", nx, ny);
  fprintf(f, "%d %d\n", minX + mox, minY + moy);
  fprintf(f, "%s %s\n", t.pair.imageName[0], t.pair.imageName[1]);

  outMap = (MapElement*) malloc(nx * ny * sizeof(MapElement));
  if (outMap == NULL)
    {
      sprintf(r.message, "Could not allocate outMap array\n");
      return;
    }
  for (y = 0; y < ny; ++y)
    for (x = 0; x < nx; ++x)
      {
	GETMAP(map, mpw, minX + x, minY + y, &rx, &ry, &rc);
	SETMAP(outMap, nx, x, y, rx, ry, rc);
    }
  fwrite(outMap, nx * ny * sizeof(MapElement), 1, f);
  fclose(f);
  free(outMap);
}


int
WriteOutputImage (char *outputName, int level, float *output,
		  int w, int h, float scale)
{
  int x, y;
  int bv;
  unsigned char *image;
  char fn[PATH_MAX];
  char extension[8];
  char msg[PATH_MAX+256];

  /* write out the output image */
  if (c.type == 't')
    strcpy(extension, "tif");
  else
    strcpy(extension, "pgm");

  if (level > c.outputLevel)
    sprintf(fn, "%s.%0.2d.%s", outputName, level, extension);
  else
    sprintf(fn, "%s.%s", c.outputCorrelationName, extension);
  image = (unsigned char *) malloc(w * h * sizeof(unsigned char));
  for (y = 0; y < h; ++y)
    for (x = 0; x < w; ++x)
      {
	if (output[y * w + x] <= 0.0)
	  bv = 0;
	else
	  {
	    bv = (int) floor(output[y * w + x] * scale);
	    if (bv > 255)
	      bv = 255;
	  }
	image[y * w + x] = bv;
      }
  if (!WriteImage(fn, image, w, h, UncompressedImage, msg))
    Error("Could not write output image %s:\n%s", fn, msg);
  free(image);
  return(1);
}

void
ComputeWarpedImage (float *warped, unsigned char *valid,
		    int w, int h,
		    int imgox, int imgoy,
		    float *image, unsigned char *mask,
		    int iw, int ih,
		    MapElement *map,
		    int mapFactor,
		    int mpw, int mph,
		    int mox, int moy)
{
  float xv, yv;
  int ixv, iyv;
  int irx, iry;
  float rrx, rry;
  float r00, r01, r10, r11;
  float rv;
  float rx00, rx01, rx10, rx11;
  float ry00, ry01, ry10, ry11;
  float rc00, rc01, rc10, rc11;
  float rx, ry, rc;
  int x, y;
  int bv;
  int mbpl;

  mbpl = (iw + 7) >> 3;
  memset(valid, 0, w * h * sizeof(unsigned char));

  for (y = 0; y < h; ++y)
    for (x = 0; x < w; ++x)
      {
	// use bilinear interpolation to find value                                   
	xv = (x + 0.5 + imgox) / mapFactor - mox;
	yv = (y + 0.5 + imgoy) / mapFactor - moy;
	ixv = ((int) (xv + 2.0)) - 2;
	iyv = ((int) (yv + 2.0)) - 2;
	rrx = xv - ixv;
	rry = yv - iyv;
	if (ixv < 0 || ixv >= mpw-1 ||
	    iyv < 0 || iyv >= mph-1)
	  {
	    warped[y * w + x] = 0.0;
	    continue;
	  }

	GETMAP(map, mpw, ixv, iyv, &rx00, &ry00, &rc00);
	GETMAP(map, mpw, ixv, iyv+1, &rx01, &ry01, &rc01);
	GETMAP(map, mpw, ixv+1, iyv, &rx10, &ry10, &rc10);
	GETMAP(map, mpw, ixv+1, iyv+1, &rx11, &ry11, &rc11);

	rx = rx00 * (rrx - 1.0) * (rry - 1.0)
	  - rx10 * rrx * (rry - 1.0)
	  - rx01 * (rrx - 1.0) * rry
	  + rx11 * rrx * rry;
	ry = ry00 * (rrx - 1.0) * (rry - 1.0)
	  - ry10 * rrx * (rry - 1.0)
	  - ry01 * (rrx - 1.0) * rry
	  + ry11 * rrx * rry;
	rc = rc00 * (rrx - 1.0) * (rry - 1.0)
	  - rc10 * rrx * (rry - 1.0)
	  - rc01 * (rrx - 1.0) * rry
	  + rc11 * rrx * rry;
	rx = mapFactor * rx;
	ry = mapFactor * ry;
	rc = rc;

	rx -= 0.5;
	ry -= 0.5;
	irx = (int) floor(rx);
	iry = (int) floor(ry);
	if (irx < 0 || irx >= iw - 1 ||
	    iry < 0 || iry >= ih - 1)
	  {
	    warped[y * w + x] = 0.0;
	    continue;
	  }
	rrx = rx - irx;
	rry = ry - iry;

#if MASKING
	if ((mask[iry*mbpl + (irx >> 8)] & (0x80 >> (irx & 7))) == 0 ||
	    rrx > 0.0 && (mask[iry*mbpl + ((irx+1) >> 3)] & (0x80 >> ((irx+1) & 7))) == 0 ||
	    rry > 0.0 && (mask[(iry+1)*mbpl + (irx >> 3)] & (0x80 >> (irx & 7))) == 0 ||
	    rrx > 0.0 && rry > 0.0 && (mask[(iry+1)*mbpl +((irx+1) >> 3)] & (0x80 >> ((irx+1) & 7))) == 0)
	  {
	    warped[y * w + x] = 0.0;
	    continue;
	  }
#endif

	r00 = image[iry * w + irx];
	r01 = image[(iry + 1) * w + irx];
	r10 = image[iry * w + (irx + 1)];
	r11 = image[(iry + 1) * w + irx + 1];
	rv = r00 * (rrx - 1.0) * (rry - 1.0)
	  - r10 * rrx * (rry - 1.0)
	  - r01 * (rrx - 1.0) * rry
	  + r11 * rrx * rry;
	bv = (int) ((0.5 + 0.5 * rc) * rv);
	if (bv < 0)
	  bv = 0;
	else if (bv > 255)
	  bv = 255;
	warped[y * w + x] = bv;
	valid[y * w + x] = 1;
      }
}

void
ComputeCorrelation (float *correlation,
		    float *a, float *b,
		    unsigned char *valid,
		    int w, int h,
		    int hw)
{
  int startN;
  double startSumA, startSumA2, startSumB, startSumB2, startSumAB;
  int n;
  double sumA, sumA2, sumB, sumB2, sumAB;
  double meanA, meanB;
  double denom;
  double av, bv;
  int x, y;
  int xc, yc;
  int i;
  int *lim;

  // then compute the correlation
  lim = (int*) malloc((hw+1) * sizeof(int));
  for (x = 0; x <= hw; ++x)
    lim[x] = (int) floor(sqrt((hw + 0.5) * (hw + 0.5) - x * x));

  // first compute the sums for (-1, -1)
  n = 0;
  sumA = 0.0;
  sumB = 0.0;
  sumA2 = 0.0;
  sumB2 = 0.0;
  sumAB = 0.0;
  for (y = 0; y < hw; ++y)
    for (x = 0; x < lim[y+1]; ++x)
      {
	if (!valid[y*w+x])
	  continue;
	av = a[y*w+x];
	bv = b[y*w+x];
	startSumA += av;
	startSumA2 += av * av;
	startSumB += bv;
	startSumB2 += bv * bv;
	startSumAB += av * bv;
	++startN;
      }
  // go through y
  for (yc = 0; yc < h; ++yc)
    {
      for (i = 1; i <= hw; ++i)
	{
	  // subtract old top
	  x = i - 1;
	  y = yc - lim[i] - 1;
	  if (y >= 0 && y < h && x < w && valid[y*w+x])
	    {
	      av = a[y*w+x];
	      bv = b[y*w+x];
	      startSumA -= av;
	      startSumA2 -= av * av;
	      startSumB -= bv;
	      startSumB2 -= bv * bv;
	      startSumAB -= av * bv;
	      --startN;
	    }

	  // add new bottom
	  y = yc + lim[i];
	  if (y < h && x < w && valid[y*w+x])
	    {
	      av = a[y*w+x];
	      bv = b[y*w+x];
	      startSumA += av;
	      startSumA2 += av * av;
	      startSumB += bv;
	      startSumB2 += bv * bv;
	      startSumAB += av * bv;
	      ++startN;
	    }
	}

      sumA = startSumA;
      sumA2 = startSumA2;
      sumB = startSumB;
      sumB2 = startSumB2;
      sumAB = startSumAB;
      n = startN;

      for (xc = 0; xc < w; ++xc)
	{
	  // subtract old left
	  x = xc - lim[0] - 1;
	  if (x >= 0 && x < w && valid[yc*w+x])
	    {
	      av = a[yc*w+x];
	      bv = b[yc*w+x];
	      sumA -= av;
	      sumA2 -= av * av;
	      sumB -= bv;
	      sumB2 -= bv * bv;
	      sumAB -= av * bv;
	      --n;
	    }

	  // add new right
	  x = xc + lim[0];
	  if (x >= 0 && x < w && valid[yc*w+x])
	    {
	      av = a[yc*w+x];
	      bv = b[yc*w+x];
	      sumA -= av;
	      sumA2 -= av * av;
	      sumB -= bv;
	      sumB2 -= bv * bv;
	      sumAB -= av * bv;
	      ++n;
	    }

	  for (i = 1; i <= hw; ++i)
	    {
	      // subtract old left
	      x = xc - lim[i] - 1;
	      if (x < 0)
		continue;
	      y = yc - i;
	      if (y >= 0 && valid[y*w+x])
		{
		  av = a[y*w+x];
		  bv = b[y*w+x];
		  sumA -= av;
		  sumA2 -= av * av;
		  sumB -= bv;
		  sumB2 -= bv * bv;
		  sumAB -= av * bv;
		  --n;
		}
	      y = yc + i;
	      if (y < h && valid[y*w+x])
		{
		  av = a[y*w+x];
		  bv = b[y*w+x];
		  sumA -= av;
		  sumA2 -= av * av;
		  sumB -= bv;
		  sumB2 -= bv * bv;
		  sumAB -= av * bv;
		  --n;
		}

	      // add new right
	      x = xc + lim[i];
	      if (x >= h)
		continue;
	      y = yc - i;
	      if (y >= 0 && valid[y*w+x])
		{
		  av = a[y*w+x];
		  bv = b[y*w+x];
		  sumA += av;
		  sumA2 += av * av;
		  sumB += bv;
		  sumB2 += bv * bv;
		  sumAB += av * bv;
		  ++n;
		}
	      y = yc + i;
	      if (y < h && valid[y*w+x])
		{
		  av = a[y*w+x];
		  bv = b[y*w+x];
		  sumA += av;
		  sumA2 += av * av;
		  sumB += bv;
		  sumB2 += bv * bv;
		  sumAB += av * bv;
		  ++n;
		}
	    }

	  // calculate the correlation for the disk surrounding (xc, yc)
	  if (n > 0)
	    {
	      meanA = sumA / n;
	      meanB = sumB / n;
	      denom = (sumA2 - 2.0 * meanA * sumA + n * meanA * meanA) *
		(sumB2 - 2.0 * meanB * sumB + n * meanB * meanB);
	      if (denom < 0.001)
		correlation[yc*w+xc] = -1000000.0;
	      else
		correlation[yc*w+xc] = (sumAB - meanA * sumB - meanB * sumA + n * meanA * meanB) / sqrt(denom);
	    }
	  else
	    correlation[yc*w+xc] = -1000000.0;
	}
    }
  free(lim);
}



void
PackContext ()
{
  int i;

  par_pkbyte((unsigned char) c.type);
  par_pkstr(c.imageBasename);
  par_pkstr(c.maskBasename);
  par_pkstr(c.discontinuityBasename);
  par_pkint(c.strictMasking);

  par_pkstr(c.cptsName);

  par_pkstr(c.inputMapName);
  par_pkstr(c.constrainingMapName);

  par_pkstr(c.outputMapBasename);
  par_pkstr(c.outputImagesName);
  par_pkstr(c.outputCorrelationName);

  par_pkstr(c.logBasename);

  par_pkint(c.startLevel);
  par_pkint(c.outputLevel);
  par_pkint(c.minResolution);
  par_pkint(c.depth);
  par_pkdouble(c.distortion);
  par_pkdouble(c.correspondence);
  par_pkdouble(c.correspondenceThreshold);
  par_pkdouble(c.constraining);
  par_pkdouble(c.constrainingThreshold);
  par_pkdouble(c.quality);
  par_pkdouble(c.minOverlap);

  par_pkint(c.writeAllMaps);
  par_pkint(c.update);
  par_pkint(c.partial);
}

void
UnpackContext ()
{
  c.type = (char) par_upkbyte();
  par_upkstr(c.imageBasename);
  par_upkstr(c.maskBasename);
  par_upkstr(c.discontinuityBasename);
  c.strictMasking = par_upkint();

  par_upkstr(c.cptsName);

  par_upkstr(c.inputMapName);
  par_upkstr(c.constrainingMapName);

  par_upkstr(c.outputMapBasename);
  par_upkstr(c.outputImagesName);
  par_upkstr(c.outputCorrelationName);

  par_upkstr(c.logBasename);

  c.startLevel = par_upkint();
  c.outputLevel = par_upkint();
  c.minResolution = par_upkint();
  c.depth = par_upkint();
  c.distortion = par_upkdouble();
  c.correspondence = par_upkdouble();
  c.correspondenceThreshold = par_upkdouble();
  c.constraining = par_upkdouble();
  c.constrainingThreshold = par_upkdouble();
  c.quality = par_upkdouble();
  c.minOverlap = par_upkdouble();

  c.writeAllMaps = par_upkint();
  c.update = par_upkint();
  c.partial = par_upkint();
}

void
PackPair (Pair *p)
{
  int imi;

  for (imi = 0; imi < 2; ++imi)
    {
      par_pkstr(p->imageName[imi]);
      par_pkint(p->imageMinX[imi]);
      par_pkint(p->imageMaxX[imi]);
      par_pkint(p->imageMinY[imi]);
      par_pkint(p->imageMaxY[imi]);
    }
  par_pkstr(p->pairName);
}

void
UnpackPair (Pair *p)
{
  char s[PATH_MAX];
  int imi;

  for (imi = 0; imi < 2; ++imi)
    {
      if (p->imageName[imi] != NULL)
	free(p->imageName[imi]);
      par_upkstr(s);
      p->imageName[imi] = (char*) malloc(strlen(s)+1);
      strcpy(p->imageName[imi], s);
      
      p->imageMinX[imi] = par_upkint();
      p->imageMaxX[imi] = par_upkint();
      p->imageMinY[imi] = par_upkint();
      p->imageMaxY[imi] = par_upkint();
    }

  if (p->pairName != NULL)
    free(p->pairName);
  par_upkstr(s);
  p->pairName = (char*) malloc(strlen(s)+1);
  strcpy(p->pairName, s);
}

void
PackTask ()
{
  PackPair(&(t.pair));
}

void
UnpackTask ()
{
  UnpackPair(&(t.pair));
}

void
PackResult ()
{
  PackPair(&(r.pair));
  par_pkint(r.updated);
  par_pkdouble(r.distortion);
  par_pkdouble(r.correlation);
  par_pkdouble(r.correspondence);
  par_pkdouble(r.constraining);
  par_pkstr(r.message);
}

void
UnpackResult ()
{
  UnpackPair(&(r.pair));
  r.updated = par_upkint();
  r.distortion = par_upkdouble();
  r.correlation = par_upkdouble();
  r.correspondence = par_upkdouble();
  r.constraining = par_upkdouble();
  par_upkstr(r.message);
}

size_t
CountBits (unsigned char *p, size_t n)
{
  size_t ii;
  unsigned char b;
  size_t sum = 0;

  for (ii = 0; ii < n; ++ii)
    {
      b = *p++;
      while (b != 0)
	{
	  if (b & 1)
	    ++sum;
	  b >>= 1;
	}
    }
  return(sum);
}

int
CreateDirectories (char *fn)
{
  int pos;
  char dn[PATH_MAX];
  char *slash;
  int len;
  struct stat statBuf;
  int hv;
  int i;

  for (pos = 0;
       (slash = strchr(&fn[pos], '/')) != NULL;
       pos = len + 1)
    {
      len = slash - fn;
      if (len == 0)
	continue;
      strncpy(dn, fn, len);
      dn[len] = '\0';

      /* check if in hash table */
      hv = 0;
      for (i = 0; i < len; ++i)
	hv = 239*hv + dn[i];
      hv &= DIR_HASH_SIZE-1;
      i = hv;
      while (dirHash[i] != NULL)
	{
	  if (strcmp(dirHash[i], dn) == 0)
	    break;
	  i = (i + 1) & (DIR_HASH_SIZE-1);
	  if (i == hv)
	    {
	      Log("Directory hash table is full!\n");
	      return(0);
	    }
	}
      if (dirHash[i] != NULL)
	continue;
      dirHash[i] = (char *) malloc(strlen(dn)+1);
      strcpy(dirHash[i], dn);

      if (stat(dn, &statBuf) == 0)
	{
	  if (S_ISDIR(statBuf.st_mode) ||
	      S_ISLNK(statBuf.st_mode))
	    continue;
	  Log("Output path component %s is not a directory\n", dn);
	  return(0);
	}
      if (errno != ENOENT)
	{
	  Log("Could not stat directory %s\n", dn);
	  return(0);
	}
      
      if (mkdir(dn, 0755) != 0)
	{
	  Log("Could not create directory %s\n", dn);
	  return(0);
	}
    }
  return(1);
}

#if FOLDING
void
FloodFill (int x, int y, int dmaw, int dmah,
	   unsigned char *disc, int w, int h,
	   int ox, int oy)
{
  int y1;
  int left, right;
  int sp;
  unsigned int *stack;
  int bpl, dmabpl;

  bpl = (w+7) >> 3;
  dmabpl = (dmaw+7) >> 3;
  sp = 0;
  stack = (unsigned int *) malloc(h * w * sizeof(int));
  stack[sp++] = y*w+x;
  while (sp > 0)
    {
      --sp;
      y = stack[sp] / w;
      x = stack[sp] - y * w;
      y1 = y;
      while (y1 >= 0 &&
             idisc[(y1 + oy)*bpl + ((x + ox) >> 3)] & (0x80 >> ((x + ox) & 7)) == 1 &&
             dma[y1*dmabpl + (x >> 3)] & (0x80 >> (x & 7)) != 0)
        --y1;
      ++y1;
      left = 0;
      right = 0;
      while (y1 < h &&
             background[y1*w + x] &&
             dma[y1*dmabpl + (x >> 3)] & (0x80 >> (x & 7)) != 0)
        {
          dma[y1*dmabpl + (x >> 3)] &= ~(0x80 >> (x & 7));
          if (!left &&
              x > 0 &&
	      background[y1*w + (x-1)] &&
	      (dma[y1*dmabpl + ((x-1) >> 3)] & (0x80 >> ((x-1) & 7))) != 0)
            {
              stack[sp++] = y1 * w + (x-1);
              left = 1;
            }
          else if (left &&
                   x > 0 &&
                   (!background[y1*w + (x-1)] ||
                    (dma[y1*dmabpl + ((x-1) >> 3)] & (0x80 >> ((x-1) & 7))) == 0))
            left = 0;

          if (!right &&
              x < w-1 &&
              background[y1*w + (x+1)] &&
              (dma[y1*dmabpl + ((x+1) >> 3)] & (0x80 >> ((x+1) & 7))) != 0)
            {
              stack[sp++] = y1 * w + (x+1);
              right = 1;
            }
          else if (right &&
                   x < w-1 &&
                   (!background[y1*w + (x+1)] ||
                    (dma[y1*dmabpl + ((x+1) >> 3)] & (0x80 >> ((x+1) & 7))) == 0))
            right = 0;
          ++y1;
        }
    }
  free(stack);
}
#endif

void Error (char *fmt, ...)
{
  va_list args;

  if (logFile != NULL)
    {
      va_start(args, fmt);
      fprintf(logFile, "%f: ERROR: ", GetCurrentTime());
      vfprintf(logFile, fmt, args);
      va_end(args);
      fflush(logFile);
    }
  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
  fflush(stderr);
  abort();
}

void Log (char *fmt, ...)
{
  va_list args;

  if (logFile == NULL)
    {
      int i;
      char logName[PATH_MAX];
      char logPrefix[PATH_MAX];
      if (c.logBasename[0] != '\0')
	strcpy(logPrefix, c.logBasename);
      else
	strcpy(logPrefix, "logs/");
      i = par_instance();
      if (i < 0)
	sprintf(logName, "%s00.log", logPrefix);
      else if (i < 100)
	sprintf(logName, "%s%0.2d.log", logPrefix, i);
      else
	sprintf(logName, "%s%d.log", logPrefix, i);
      logFile = fopen(logName, "w");
      if (logFile == NULL)
	Error("Could not open log file %s\n", logName);
    }

  va_start(args, fmt);
  fprintf(logFile, "%f: ", GetCurrentTime());
  vfprintf(logFile, fmt, args);
  va_end(args);
  fflush(logFile);
}

double
GetCurrentTime ()
{
  struct timeval current_time;

  gettimeofday(&current_time, NULL);
  return(((double) current_time.tv_sec) +
         ((double) current_time.tv_usec) / 1000000.0);
}
