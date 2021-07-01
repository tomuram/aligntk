/*
 *  warp3.c  -  registers two grayscale images, producing a map of
 *                the corresponding points in the two images
 *
 *  Copyright (c) 2006-2011 National Resource for Biomedical
 *                          Supercomputing,
 *                          Pittsburgh Supercomputing Center,
 *                          Carnegie Mellon University
 *
 *  Distribution of this code is prohibited without the prior written
 *  consent of the National Resource for Biomedical Supercomputing.
 *
 *  Acknowledgements:
 *     Development of this code was supported in part by
 *       NIH grant 5P41RR006009.
 *
 *  This program is distributed in the hope that it will
 *  be useful, but WITHOUT ANY WARRANTY; without even the
 *  implied warranty of MERCHANTABILITY or FITNESS FOR A
 *  PARTICULAR PURPOSE.  Neither Carnegie Mellon University
 *  nor any of the authors assume any liability for
 *  damages, incidental or otherwise, caused by the
 *  installation or use of this software.
 *
 *  CLINICAL APPLICATIONS ARE NOT RECOMMENDED, AND THIS
 *  SOFTWARE HAS NOT BEEN EVALUATED BY THE UNITED STATES
 *  FDA FOR ANY CLINICAL USE.
 *
 *  HISTORY
 *    2006  Written by Greg Hood (ghood@psc.edu)
 *    2009    Simplified move generator, and refined
 *             confidence value computation to take into
 *             account image and reference masks (ghood@psc.edu)
 */

#include <dirent.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <sched.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>
#include <unistd.h>

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

#define TRANSLATION_METHOD	0
#define RIGID_METHOD		1
#define AFFINE_METHOD		2
#define QUADRATIC_METHOD	3

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
  char outputWarpedBasename[PATH_MAX];
  char outputCorrelationBasename[PATH_MAX];
  char outputMaskBasename[PATH_MAX];
  char outputPairwiseMaskBasename[PATH_MAX];

  char logBasename[PATH_MAX];

  int startLevel;
  int outputLevel;
  int minResolution;
  int depth;
  int cptsMethod;
  double distortion;
  double correspondence;
  double correspondenceThreshold;
  double constraining;
  double constrainingThreshold;
  double constrainingConfidenceThreshold;
  double quality;
  double minOverlap;
  double trimMapSourceThreshold;
  double trimMapTargetThreshold;

  int correlationHalfWidth;
  int writeAllMaps;
  int update;
  int partial;
  int nWorkers;
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

unsigned char *outputMasks[MAX_LEVELS];/* masks of image 0 to specify
					  the desired extent of output map */
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
MapElement *inputMap = 0;
int mpwi, mphi;
int offxi, offyi;

int constrainingMapLevel;
int constrainingMapFactor;
MapElement *constrainingMap = 0;
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
int Init ();
void Compute (char *outputName, char *outputWarpedName, char *outputCorrelationName);
void ComputeWarpedImage (float *warped, unsigned char *valid,
			 int w, int h,
			 int imgox, int imgoy,
			 float *image, unsigned char *mask,
			 int iw, int ih,
			 int refox, int refoy,
			 MapElement *map,
			 int mapFactor,
			 int mpw, int mph,
			 int mox, int moy);
void ComputeCorrelation (float *correlation,
			 float *a, float *b,
			 unsigned char *valid,
			 int w, int h,
			 int hw);
void TrimOutputMap (MapElement *map, int mpw, int mph, int mox, int moy,
		    int factor,
		    unsigned int iw, unsigned int ih, int imgox, int imgoy,
		    unsigned int rw, unsigned int rh, int refox, int refoy,
		    unsigned char *cimask,
		    unsigned char *crmask);
void WriteOutputMap (char *outputName, int level, MapElement *map,
		     int mpw, int mph, int mox, int moy);
int WriteOutputImage (char *outputName, int level, float *output,
		      int w, int h, float scale);
int Extrapolate (double *prx, double *pry, double *prc,
		 int ix, int iy, float arrx, float arry,
		 MapElement* map, int mw, int mh, float threshold);
int Compare (const void *x, const void *y);
int SortBySlice (const void *x, const void *y);
int SortByEnergy (const void *x, const void *y);
int ParseRange (char *s, int *pos, int *minValue, int *maxValue);
int ParseValue (char *s, int *pos, int *value);
size_t CountBits (unsigned char *p, size_t n);
size_t CountIntersectionBits (unsigned char *p, unsigned char *q, size_t n);
int CreateDirectories (char *fn);
void CopyString (char **dst, char *src);
void SetMessage (char *fmt, ...);
void Error (char *fmt, ...);
void Log (char *fmt, ...);
char *GetTimestamp (char *timestamp, size_t size);
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
#if 1
  // TEMPORARY FOR DEBUGGING
  struct rlimit rlim;
  rlim.rlim_cur = 10000000000;
  rlim.rlim_max = 10000000000;
  setrlimit(RLIMIT_CORE, &rlim);
#endif
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
  char outputPairsFile[PATH_MAX];
  char outputSortedPairsFile[PATH_MAX];
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
  char line[LINE_LENGTH+1];
  FILE *opf;

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
  c.outputWarpedBasename[0] = '\0';
  c.outputCorrelationBasename[0] = '\0';
  c.logBasename[0] = '\0';
  c.startLevel = -1;
  c.outputLevel = 0;
  c.minResolution = 1;
  c.distortion = 1.0;
  c.correspondence = 1.0;
  c.correspondenceThreshold = 0.0;
  c.constraining = 1.0;
  c.constrainingThreshold = 0.0;
  c.constrainingConfidenceThreshold = 0.5;
  c.quality = 10.0;
  c.minOverlap = 20.0;
  c.correlationHalfWidth = 31;
  c.writeAllMaps = 0;
  c.depth = 0;
  c.cptsMethod = -1;
  c.update = 0;
  c.partial = 0;
  c.nWorkers = par_workers();
  c.trimMapSourceThreshold = 0.0;
  c.trimMapTargetThreshold = 0.0;
  r.pair.imageName[0]  = r.pair.imageName[1] = NULL;
  r.pair.pairName = NULL;
  r.message = NULL;
  forward = 1;
  pairsFile[0] = '\0';
  outputPairsFile[0] = '\0';
  outputSortedPairsFile[0] = '\0';
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
    else if (strcmp(argv[i], "-output_mask") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.outputMaskBasename, argv[i]);
      }
    else if (strcmp(argv[i], "-output_pairwise_mask") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.outputPairwiseMaskBasename, argv[i]);
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
    else if (strcmp(argv[i], "-output_warped") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.outputWarpedBasename, argv[i]);
      }
    else if (strcmp(argv[i], "-output_correlation") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.outputCorrelationBasename, argv[i]);
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
    else if (strcmp(argv[i], "-constraining_confidence_threshold") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%lf", &c.constrainingConfidenceThreshold) != 1)
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
    else if (strcmp(argv[i], "-min_overlap") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%lf", &c.minOverlap) != 1)
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
    else if (strcmp(argv[i], "-output_pairs") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(outputPairsFile, argv[i]);
      }
    else if (strcmp(argv[i], "-output_sorted_pairs") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(outputSortedPairsFile, argv[i]);
      }
    else if (strcmp(argv[i], "-translation") == 0)
      {
	if (c.cptsMethod >= 0)
	  {
	    error = 1;
	    break;
	  }
	c.cptsMethod = TRANSLATION_METHOD;
      }
    else if (strcmp(argv[i], "-rigid") == 0)
      {
	if (c.cptsMethod >= 0)
	  {
	    error = 1;
	    break;
	  }
	c.cptsMethod = RIGID_METHOD;
      }
    else if (strcmp(argv[i], "-affine") == 0)
      {
	if (c.cptsMethod >= 0)
	  {
	    error = 1;
	    break;
	  }
	c.cptsMethod = AFFINE_METHOD;
      }
    else if (strcmp(argv[i], "-quadratic") == 0)
      {
	if (c.cptsMethod >= 0)
	  {
	    error = 1;
	    break;
	  }
	c.cptsMethod = QUADRATIC_METHOD;
      }
    else if (strcmp(argv[i], "-trim_map_source_threshold") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%lf", &c.trimMapSourceThreshold) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-trim_map_target_threshold") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%lf", &c.trimMapTargetThreshold) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else error = 1;

  if (error)
    {
      fprintf(stderr, "Usage: warp3 -input file_prefix -output file_prefix\n");
      fprintf(stderr, "              [-resolution HORIZxVERT]\n");
      fprintf(stderr, "              [-distortion distortion_weight]\n");
      fprintf(stderr, "              [-correspondence correspondence_weight]\n");
      fprintf(stderr, "              [-correspondence_threshold pixels]\n");
      fprintf(stderr, "              [-translation]\n");
      fprintf(stderr, "              [-rigid]\n");
      fprintf(stderr, "              [-affine]\n");
      fprintf(stderr, "              [-quadratic]\n");
      fprintf(stderr, "              [-quality quality_factor]\n");
      fprintf(stderr, "              [-min_res minimum_resolution_in_pixels]\n");
      fprintf(stderr, "              [-trim_map_source_threshold]\n");
      fprintf(stderr, "              [-trim_map_target_threshold]\n");
      fprintf(stderr, "              [-depth delta_depth]\n");
      fprintf(stderr, "              [-all_maps]\n");
      fprintf(stderr, "              [-update]\n");
      fprintf(stderr, "              [-partial]\n");
      fprintf(stderr, "              [-pairs <pair_file>]\n");
      fprintf(stderr, "              [-input_map <input_map_prefix>]\n");
      fprintf(stderr, "              [-constraining_map <constraining_map_prefix>]\n");
      fprintf(stderr, "              [-constraining constraining_weight]\n");
      fprintf(stderr, "              [-constraining_threshold pixels]\n");
      fprintf(stderr, "              [-constraining_confidence_threshold score]\n");
      fprintf(stderr, "              [-strict_masking]\n");
      fprintf(stderr, "              [-min_overlap percent]\n");
      fprintf(stderr, "              [-output_pairs <output_pair_file>]\n");
      fprintf(stderr, "              [-output_sorted_pairs <output_pair_file>]\n");
      fprintf(stderr, "              [-logs <log_file_directory>]\n");
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
  if (c.cptsMethod < 0)
    c.cptsMethod = AFFINE_METHOD;

  f = fopen(pairsFile, "r");
  if (f == NULL)
    Error("Could not open pairs file %s\n", pairsFile);
  while (fgets(line, LINE_LENGTH, f) != NULL)
    {
      if (line[0] == '\0' || line[0] == '#')
	continue;
      if (sscanf(line, "%s %d %d %d %d %s %d %d %d %d %s",
		 imgn[0], &imgMinX[0], &imgMaxX[0], &imgMinY[0], &imgMaxY[0],
		 imgn[1], &imgMinX[1], &imgMaxX[1], &imgMinY[1], &imgMaxY[1],
		 pairn) != 11)
	{
	  if (sscanf(line, "%s %s %s", imgn[0], imgn[1], pairn) != 3)
	    Error("Invalid line in pairs file %s:\n%s\n", pairsFile, line);
	  imgMinX[0] = -1;
	  imgMaxX[0] = -1;
	  imgMinY[0] = -1;
	  imgMaxY[0] = -1;
	  imgMinX[1] = -1;
	  imgMaxX[1] = -1;
	  imgMinY[1] = -1;
	  imgMaxY[1] = -1;
	}
      if ((nPairs & 1023) == 0)
	pairs = (Pair *) realloc(pairs, (nPairs + 1024) * sizeof(Pair));
      for (imi = 0; imi < 2; ++imi)
	{
	  pairs[nPairs].imageName[imi] = NULL;
	  CopyString(&(pairs[nPairs].imageName[imi]), imgn[imi]);
	  pairs[nPairs].imageMinX[imi] = imgMinX[imi];
	  pairs[nPairs].imageMaxX[imi] = imgMaxX[imi];
	  pairs[nPairs].imageMinY[imi] = imgMinY[imi];
	  pairs[nPairs].imageMaxY[imi] = imgMaxY[imi];
	}
      pairs[nPairs].pairName = NULL;
      CopyString(&(pairs[nPairs].pairName), pairn);
      ++nPairs;
    }
  fclose(f);

  /* check that output directories are writeable */ 
  Log("MASTER checking output directories\n");
  sprintf(fn, "%sTEST.map", c.outputMapBasename);
  if (!CreateDirectories(fn))
    Error("Could not create directory for output maps: %s\n",
	  c.outputMapBasename);
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

  if (c.outputWarpedBasename[0] != '\0')
    {
      sprintf(fn, "%sTEST.pgm", c.outputWarpedBasename);
      if (!CreateDirectories(fn))
	Error("Could not create directory for output warped images: %s\n",
	      c.outputWarpedBasename);
      f = fopen(fn, "w");
      if (f == NULL)
	{
	  for (i = strlen(c.outputWarpedBasename)-1; i >= 0 && c.outputWarpedBasename[i] != '/'; --i) ;
	  if (i != 0)
	    strncpy(outputDirName, c.outputWarpedBasename, i+1);
	  else
	    outputDirName[0] = '.';
	  outputDirName[i+1] = '\0';
	  Error("Could not open test output image file %s --\n      does directory %s exist and is it writeable?\n", fn, outputDirName);
	}
      fclose(f);
      unlink(fn);
    }

  if (c.outputCorrelationBasename[0] != '\0')
    {
      sprintf(fn, "%sTEST.pgm", c.outputCorrelationBasename);
      if (!CreateDirectories(fn))
	Error("Could not create directory for output correlation: %s\n",
	      c.outputCorrelationBasename);
      f = fopen(fn, "w");
      if (f == NULL)
	{
	  for (i = strlen(c.outputCorrelationBasename)-1; i >= 0 && c.outputCorrelationBasename[i] != '/'; --i) ;
	  if (i != 0)
	    strncpy(outputDirName, c.outputCorrelationBasename, i+1);
	  else
	    outputDirName[0] = '.';
	  outputDirName[i+1] = '\0';
	  Error("Could not open test output image file %s --\n      does directory %s exist and is it writeable?\n", fn, outputDirName);
	}
      fclose(f);
      unlink(fn);
    }

  if (outputPairsFile[0] != '\0')
    {
      if (!CreateDirectories(outputPairsFile))
	Error("Could not create directory for output pairs file: %s\n",
	      outputPairsFile);
      f = fopen(outputPairsFile, "w");
      if (f == NULL)
	Error("Could not open output pairs file %s\n", outputPairsFile);
      fclose(f);
      unlink(outputPairsFile);
    }

  if (outputSortedPairsFile[0] != '\0')
    {
      if (!CreateDirectories(outputSortedPairsFile))
	Error("Could not create directory for output sorted pairs file: %s\n",
	      outputSortedPairsFile);
      f = fopen(outputSortedPairsFile, "w");
      if (f == NULL)
	Error("Could not open output sorted pairs file %s\n",
	      outputSortedPairsFile);
      fclose(f);
      unlink(outputSortedPairsFile);
    }

  if (summaryName[0] != '\0')
    {
      if (!CreateDirectories(summaryName))
	Error("Could not create directory for summary file: %s\n",
	      summaryName);
      f = fopen(summaryName, "w");
      if (f == NULL)
	Error("Could not open summary file %s\n",
	      outputSortedPairsFile);
      fclose(f);
      unlink(summaryName);
    }

  if (c.logBasename[0] != '\0')
    {
      sprintf(fn, "%sTEST.log", c.logBasename);
      if (!CreateDirectories(fn))
	Error("Could not create directory for test logfile: %s\n", fn);
      f = fopen(fn, "w");
      if (f == NULL)
	Error("Could not open test log file %s\n", fn);
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
      for (imi = 0; imi < 2; ++imi)
	{
	  CopyString(&(t.pair.imageName[imi]), pairs[pn].imageName[imi]);
	  t.pair.imageMinX[imi] = pairs[pn].imageMinX[imi];
	  t.pair.imageMaxX[imi] = pairs[pn].imageMaxX[imi];
	  t.pair.imageMinY[imi] = pairs[pn].imageMinY[imi];
	  t.pair.imageMaxY[imi] = pairs[pn].imageMaxY[imi];
	}
      CopyString(&(t.pair.pairName), pairs[pn].pairName);

      // make sure that output directories exist
      sprintf(fn, "%s%s.map", c.outputMapBasename, t.pair.pairName);
      if (!CreateDirectories(fn))
	continue;
      if (c.outputWarpedBasename[0] != '\0')
	{
	  sprintf(fn, "%s%s.pgm", c.outputWarpedBasename, t.pair.pairName);
	  if (!CreateDirectories(fn))
	    continue;
	}
      if (c.outputCorrelationBasename[0] != '\0')
	{
	  sprintf(fn, "%s%s.pgm", c.outputCorrelationBasename, t.pair.pairName);
	  if (!CreateDirectories(fn))
	    continue;
	}

      Log("Delegating pair %d\n", pn);
      par_delegate_task();
    }
  par_finish();

  if (outputPairsFile[0] != '\0')
    {
      qsort(results, nResults, sizeof(Result), SortBySlice);
      opf = fopen(outputPairsFile, "w");
      if (opf == NULL)
	Error("Could not open output pairs file %s for writing.\n",
	      outputPairsFile);
      for (i = 0; i < nResults; ++i)
	if (results[i].updated)
	  fprintf(opf, "%s %d %d %d %d %s %d %d %d %d %s\n",
		  results[i].pair.imageName[0],
		  results[i].pair.imageMinX[0],
		  results[i].pair.imageMaxX[0],
		  results[i].pair.imageMinY[0],
		  results[i].pair.imageMaxX[0],
		  results[i].pair.imageName[1],
		  results[i].pair.imageMinX[1],
		  results[i].pair.imageMaxX[1],
		  results[i].pair.imageMinY[1],
		  results[i].pair.imageMaxX[1],
		  results[i].pair.pairName);
      fclose(opf);
    }
  if (outputSortedPairsFile[0] != '\0')
    {
      qsort(results, nResults, sizeof(Result), SortByEnergy);
      opf = fopen(outputSortedPairsFile, "w");
      if (opf == NULL)
	Error("Could not open output sorted pairs file %s for writing.\n",
	      outputSortedPairsFile);
      for (i = 0; i < nResults; ++i)
	if (results[i].updated)
	  fprintf(opf, "%s %d %d %d %d %s %d %d %d %d %s\n",
		  results[i].pair.imageName[0],
		  results[i].pair.imageMinX[0],
		  results[i].pair.imageMaxX[0],
		  results[i].pair.imageMinY[0],
		  results[i].pair.imageMaxX[0],
		  results[i].pair.imageName[1],
		  results[i].pair.imageMinX[1],
		  results[i].pair.imageMaxX[1],
		  results[i].pair.imageMinY[1],
		  results[i].pair.imageMaxX[1],
		  results[i].pair.pairName);
      fclose(opf);
    }

  if (summaryName[0] != '\0')
    {
      qsort(results, nResults, sizeof(Result), SortBySlice);
      f = fopen(summaryName, "w");
      if (f == NULL)
	Error("Could not open summary output file %s\n", summaryName);

      fprintf(f, "Sorted by slice:\n");
      fprintf(f, "IMAGE    REFERENCE  CORRELATION   DISTORTION    CORRESPOND    CONSTRAIN     ENERGY\n");
      for (i = 0; i < nResults; ++i)
	fprintf(f, "%8s %8s %12.6f  %12.6f  %12.6f  %12.6f  %12.6f\n",
		results[i].pair.imageName[0],
		results[i].pair.imageName[1],
		results[i].correlation,
		results[i].distortion,
		results[i].correspondence,
		results[i].constraining,
		results[i].distortion * c.distortion - results[i].correlation +
		results[i].correspondence * c.correspondence + results[i].constraining * c.constraining);
      fprintf(f, "\n\nSorted by energy:\n");
      fprintf(f, "IMAGE    REFERENCE  CORRELATION   DISTORTION    CORRESPOND    CONSTRAIN     ENERGY\n");
      qsort(results, nResults, sizeof(Result), SortByEnergy);
      for (i = 0; i < nResults; ++i)
	fprintf(f, "%8s %8s  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f\n",
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
  int imi;

  if (r.message != NULL)
    Error("\nThe following error was encountered by one of the worker processes:%s\n", r.message);

  results = (Result*) realloc(results, (nResults + 1) * sizeof(Result));
  for (imi = 0; imi < 2; ++imi)
    {
      results[nResults].pair.imageName[imi] = NULL;
      CopyString(&(results[nResults].pair.imageName[imi]), r.pair.imageName[imi]);
      results[nResults].pair.imageMinX[imi] = r.pair.imageMinX[imi];
      results[nResults].pair.imageMaxX[imi] = r.pair.imageMaxX[imi];
      results[nResults].pair.imageMinY[imi] = r.pair.imageMinY[imi];
      results[nResults].pair.imageMaxY[imi] = r.pair.imageMaxY[imi];
    }
  results[nResults].pair.pairName = NULL;
  CopyString(&(results[nResults].pair.pairName), r.pair.pairName);
  results[nResults].updated = r.updated;
  results[nResults].distortion = r.distortion;
  results[nResults].correlation = r.correlation;
  results[nResults].correspondence = r.correspondence;
  results[nResults].constraining = r.constraining;
  results[nResults].message = NULL;
  CopyString(&(results[nResults].message), r.message);

  if ((nResults % 50) == 0 && nResults != 0)
    printf(" %d \n                   ", nResults);
  printf(".");
  ++nResults;
}

int
Extrapolate (double *prx, double *pry, double *prc,
	     int ix, int iy, float arrx, float arry,
	     MapElement* map, int mw, int mh, float threshold)
{
  float rx = 0.0;
  float ry = 0.0;
  float rc = 1.0;
  float totalWeight = 0.0;
  int maxDelta = ((int) ceil(threshold)) + 2;
  float rrx, rry;
  int dx, dy;
  int ixv, iyv;
  float rx00, rx01, rx10, rx11;
  float ry00, ry01, ry10, ry11;
  float rc00, rc01, rc10, rc11;
  float weight;
  float rcMin;

  for (dy = -maxDelta; dy <= maxDelta; ++dy)
    {
      iyv = iy + dy;
      if (iyv < 0 || iyv >= mh-1)
	continue;
      rry = arry - dy;

      for (dx = -maxDelta; dx <= maxDelta; ++dx)
	{
	  ixv = ix + dx;
	  if (ixv < 0 || ixv >= mw-1 ||
	      map[iyv*mw+ixv].c == 0.0 ||
	      map[(iyv+1)*mw+ixv].c == 0.0 ||
	      map[iyv*mw+ixv+1].c == 0.0 ||
	      map[(iyv+1)*mw+ixv+1].c == 0.0)
	    continue;

	  rrx = arrx - dx;

	  rx00 = map[iyv*mw+ixv].x;
	  ry00 = map[iyv*mw+ixv].y;
	  rc00 = map[iyv*mw+ixv].c;
	  rx01 = map[(iyv+1)*mw+ixv].x;
	  ry01 = map[(iyv+1)*mw+ixv].y;
	  rc01 = map[(iyv+1)*mw+ixv].c;
	  rx10 = map[iyv*mw+ixv+1].x;
	  ry10 = map[iyv*mw+ixv+1].y;
	  rc10 = map[iyv*mw+ixv+1].c;
	  rx11 = map[(iyv+1)*mw+ixv+1].x;
	  ry11 = map[(iyv+1)*mw+ixv+1].y;
	  rc11 = map[(iyv+1)*mw+ixv+1].c;

	  if (dx == 0 && dy == 0)
	    weight = M_SQRT2;
	  else
	    weight = 1.0 / hypot((float) dx, (float) dy);
	  rx += weight * (rx00 * (rrx - 1.0) * (rry - 1.0)
			  - rx10 * rrx * (rry - 1.0) 
			  - rx01 * (rrx - 1.0) * rry
			  + rx11 * rrx * rry);
	  ry += weight * (ry00 * (rrx - 1.0) * (rry - 1.0)
			  - ry10 * rrx * (rry - 1.0) 
			  - ry01 * (rrx - 1.0) * rry
			  + ry11 * rrx * rry);
	  rcMin = rc00;
	  if (rc01 < rcMin)
	    rcMin = rc01;
	  if (rc10 < rcMin)
	    rcMin = rc10;
	  if (rc11 < rcMin)
	    rcMin = rc11;
	  rc += weight * rcMin;
	  totalWeight += weight;
	}
    }
  if (totalWeight == 0.0)
    return(0);
  *prx = rx / totalWeight;
  *pry = ry / totalWeight;
  *prc = rc / totalWeight;
  return(1);
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
  r.message = NULL;
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
  size_t ombpl;
  float *warp;
  unsigned char *warpout;
  int level;
  int imw, imh;
  int omw, omh;
  int mpw, mph;
  int maskPresent;
  int outputMaskPresent;
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
  char outputName[PATH_MAX], outputWarpedName[PATH_MAX], outputCorrelationName[PATH_MAX];
  char outputMaskName[PATH_MAX];
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
      //      sprintf(imageName[imi], "%s%s.%s", c.imageBasename, t.pair.imageName[imi],
      //	      c.type == 't' ? "tif" : "pgm");
      sprintf(imageName[imi], "%s%s", c.imageBasename, t.pair.imageName[imi]);
      if (c.maskBasename[0] != '\0')
	sprintf(maskName[imi], "%s%s", c.maskBasename, t.pair.imageName[imi]);
      else
	maskName[imi][0] = '\0';
      if (c.discontinuityBasename[0] != '\0')
	sprintf(discontinuityName[imi], "%s%s", c.discontinuityBasename, t.pair.imageName[imi]);
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
  if (c.outputWarpedBasename[0] != '\0')
    sprintf(outputWarpedName, "%s%s", c.outputWarpedBasename, t.pair.pairName);
  else
    outputWarpedName[0] = '\0';
  if (c.outputCorrelationBasename[0] != '\0')
    sprintf(outputCorrelationName, "%s%s.pgm", c.outputCorrelationBasename, t.pair.pairName);
  else
    outputCorrelationName[0] = '\0';
  if (c.outputPairwiseMaskBasename[0] != '\0')
    sprintf(outputMaskName, "%s%s.pbm.gz", c.outputPairwiseMaskBasename,
	    t.pair.pairName);
  else if (c.outputMaskBasename[0] != '\0')
    sprintf(outputMaskName, "%s%s.pbm.gz", c.outputMaskBasename,
	    t.pair.imageName[0]);

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
  if (!computeMap && outputMaskName[0] != '\0')
    {
      if (stat(outputMaskName, &sb) == 0 &&
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

      for (imi = 0; imi < 2; ++imi)
	{
	  CopyString(&(r.pair.imageName[imi]), t.pair.imageName[imi]);
	  r.pair.imageMinX[imi] = t.pair.imageMinX[imi];
	  r.pair.imageMaxX[imi] = t.pair.imageMaxX[imi];
	  r.pair.imageMinY[imi] = t.pair.imageMinY[imi];
	  r.pair.imageMaxY[imi] = t.pair.imageMaxY[imi];
	}
      CopyString(&(r.pair.pairName), t.pair.pairName);
      r.updated = 0;
      CopyString(&(r.message), NULL);
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
	  for (imi = 0; imi < 2; ++imi)
	    {
	      CopyString(&(r.pair.imageName[imi]), t.pair.imageName[imi]);
	      r.pair.imageMinX[imi] = t.pair.imageMinX[imi];
	      r.pair.imageMaxX[imi] = t.pair.imageMaxX[imi];
	      r.pair.imageMinY[imi] = t.pair.imageMinY[imi];
	      r.pair.imageMaxY[imi] = t.pair.imageMaxY[imi];
	    }
	  CopyString(&(r.pair.pairName), t.pair.pairName);
	  r.updated = 0;
	  r.correlation = 0.0;
	  r.distortion = 0.0;
	  r.correspondence = 0.0;
	  r.constraining = 0.0;
	  CopyString(&(r.message), NULL);
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
		     errorMsg))
	{
	  SetMessage("Could not read image %s\n", imageName[imi]);
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
	  SetMessage("Could not allocate image arrays (%zd)\n",
		     imagePixels * sizeof(float));
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
			  errorMsg))
	    {
	      SetMessage("Could not read mask %s\n", maskName[imi]);
	      return;
	    }
	  Log("WORKER READ THE IMAGE MASK\n");

	  Log("imw[%d] = %d imh[%d] = %d imageWidth[%d][0] = %d imageHeight[%d][0] = %d\n",
	      imi, imw, imi, imh, imi, imageWidth[imi][0], imi, imageHeight[imi][0]);

	  if (imw != imageWidth[imi][0] ||
	      imh != imageHeight[imi][0])
	    {
	      SetMessage("Error: incorrect reference mask size: %s\n", maskName[imi]);
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
	      SetMessage("Could not allocate image arrays (%zd)\n",
			 imageHeight[imi][0]*imbpl);
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
			  errorMsg))
	    {
	      SetMessage("Could not read reference discontinuity mask %s\n", discontinuityName);
	      return;
	    }
	  Log("WORKER READ THE DISCONTINUITY MASK\n");
	  if (imw != imageWidth[imi][0] || imh != imageHeight[imi][0])
	    {
	      SetMessage("Error: incorrect reference discontinuity mask size: %s\n", discontinuityName);
	      return;
	    }
	}
#endif
    }

#if MASKING
  outputMasks[0] = NULL;
  if (outputMaskName[0] != '\0')
    {
      if (!ReadBitmap(outputMaskName, &outputMasks[0],
		      &omw, &omh,
		      t.pair.imageMinX[0], t.pair.imageMaxX[0],
		      t.pair.imageMinY[0], t.pair.imageMaxY[0],
		      errorMsg))
	{
	  Log("Could not read output mask %s\n", outputMaskName);
	  SetMessage("Could not read output mask %s\n", outputMaskName);
	  return;
	}
      Log("omw = %d omh = %d ones = %zd  %d %d %d %d\n", omw, omh,
	  CountBits(outputMasks[0], omh * ((omw + 7) >> 3)),
	  t.pair.imageMinX[0], t.pair.imageMaxX[0],
	  t.pair.imageMinY[0], t.pair.imageMaxY[0]);
      if (omw != imageWidth[0][0] ||
	  omh != imageHeight[0][0])
	{
	  Log("Output mask %s size does not match image size\n", outputMaskName);
	  SetMessage("Output mask %s size does not match image size\n",
		     outputMaskName);
	  return;
	}
      Log("WORKER READ THE OUTPUT MASK\n");
    }
#endif

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

  Compute(outputName, outputWarpedName, outputCorrelationName);

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
      if (outputMasks[level] != NULL)
	free(outputMasks[level]);
      if (maps[level] != NULL)
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
  size_t ombpl;
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

  memset(maps, 0, MAX_LEVELS * sizeof(MapElement*));
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
	      SetMessage("Could not allocate image arrays (%zd)\n",
			 imagePixels * sizeof(float));
	      return(0);
	    }

#if MASKING
	  masks[imi][level] = (unsigned char*) malloc(ih*imbpl);
	  if (masks[imi][level] == NULL)
	    {
	      SetMessage("Could not allocate mask arrays (%zd)\n",
			 ih * imbpl);
	      return(0);
	    }

	  iCount = (unsigned char*) malloc(imagePixels * sizeof(unsigned char));
	  if (iCount == NULL)
	    {
	      SetMessage("Could not allocate count arrays (%zd)\n",
			 imagePixels * sizeof(unsigned char));
	      return(0);
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
	      SetMessage("Could not allocate discontinuity arrays (%zd)\n",
			 imph * impbpl);
	      return(0);
	    }
	  memset(imdisc[imi][level], 0, imph * impbpl);

	  // allocate a bitmap array for one map square at this level
	  dmadim = 1 << level;
	  dmabpl = (dmadim + 7) >> 3;
	  dmasize = dmadim * dmabpl;
	  dma = (unsigned char *) malloc(dmasize);
	  if (dma == NULL)
	    {
	      SetMessage("Could not allocate arrays (%zd)\n", dmasize);
	      return(0);
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

#if MASKING      
      if (outputMasks[level-1] != NULL)
	{
	  iw = imageWidth[0][level];
	  ih = imageHeight[0][level];
	  ombpl = (iw + 7) >> 3;
	  outputMasks[level] = (unsigned char *) malloc(ih * ombpl);
	  m = outputMasks[level];
	  Log("Allocated outputMasks[%d] = %p of size %zd\n",
	      level, m, ih * ombpl);
	  memset(m, 0, ih * ombpl);
	  src_ih = imageHeight[0][level-1];
	  src_iw = imageWidth[0][level-1];
	  src_mask = outputMasks[level-1];
	  src_mbpl = (src_iw + 7) >> 3;
	  Log("At level %d, next lower level output_mask contains %d bits.\n",
	      level, CountBits(src_mask, src_ih * src_mbpl));
	  delta_x = 2*imageOffsetX[0][level] - imageOffsetX[0][level-1];
	  delta_y = 2*imageOffsetY[0][level] - imageOffsetY[0][level-1];
	  Log("delta_x = %d  delta_y = %d\n", delta_x, delta_y);
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
			  if (MASK(src_mask, src_mbpl, src_x, src_y))
			    SETMASKBIT(m, ombpl, x, y);
			}
		    }
		}
	    }
	  Log("At level %d, output_mask contains %d ones.\n",
	      level, CountBits(outputMasks[level], ih * ombpl));
	}
      else
	outputMasks[level] = NULL;
#endif
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

	  if (imi == 0 && outputMasks[level] != NULL)
	    maskCount = CountIntersectionBits(masks[imi][level],
					      outputMasks[level],
					      ih*imbpl);
	  else
	    maskCount = CountBits(masks[imi][level], ih*imbpl);
	  if (maskCount < minMaskCount[level])
	    minMaskCount[level] = maskCount;
	  Log("maskCount[%d][%d] = %ld (%dx%d = %ld) %f%%\n",
	      imi, level, maskCount,
	      iw, ih, imagePixels, 100.0 * ((float) maskCount) / imagePixels);
	}
      Log("minMaskCount[%d] = %d\n", level, minMaskCount[level]);
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
Compute (char *outputName, char *outputWarpedName, char *outputCorrelationName)
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
  size_t moveCount, goalMoveCount, acceptedMoveCount;
  double l0, l1, l2, l3;
  double nl0, nl1, nl2, nl3;
  size_t imbpl, rmbpl;
  size_t ombpl;
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
  int nValid;
  int nTrimmed;
  int linePos;
  char line[65536];
  float *warpedArray;
  unsigned char *validArray;
  float *correlationArray;
  int moveDebug;
  unsigned char *m;
  int count;
  int count1, count2, count3, count4;
  size_t k;
  double *nomArea;
  double *nomL0, *nomL1, *nomL2, *nomL3;
  double *nomThetaX, *nomThetaY;
  double thetaX, thetaY;
  double txd, tyd;
  double nThetaX, nThetaY;
  double ntxd, ntyd;
  int imi;
  int dnx, dny;
  int immbpl;
  unsigned char *inputMapMask;
  float *inputMapDist;
  double lb;

  for (imi = 0; imi < 2; ++imi)
    {
      CopyString(&(r.pair.imageName[imi]), t.pair.imageName[imi]);
      r.pair.imageMinX[imi] = t.pair.imageMinX[imi];
      r.pair.imageMaxX[imi] = t.pair.imageMaxX[imi];
      r.pair.imageMinY[imi] = t.pair.imageMinY[imi];
      r.pair.imageMaxY[imi] = t.pair.imageMaxY[imi];
    }
  CopyString(&(r.pair.pairName), t.pair.pairName);
  r.updated = 1;
  r.correlation = 0.0;
  r.distortion = 0.0;
  r.correspondence = 0.0;
  r.constraining = 0.0;
  CopyString(&(r.message), NULL);

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
	  SetMessage("Could not allocate map arrays (%zd)\n",
		     mSize * sizeof(MapElement));
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
	      SetMessage("Could not allocate constraining arrays (%zd)\n",
			 mSize * sizeof(MapElement));
	      return;
	    }
	}

      if (level == startLevel)
	if (inputMapFactor != 0 && nCpts < 2)
	  {
	    /* use the provided input map as the initial map */
	    /* create an array with coverage corresponding to the
	       startLevel, and resolution determined by the inputMapLevel;
	       fill in the entries with a 1 if the map has nonzero c;
	       fill in the other entries with 0
	     compute distance
	     let the threshold be a factor from sqrt(2)/2 to sqrt(2)
	    */
	    dnx = (imageWidth[0][0] + lFactor - 1) / lFactor;
	    dny = (imageHeight[0][0] + lFactor - 1) / lFactor;
	    dnx = dnx * lFactor / inputMapFactor + 1;
	    dny = dny * lFactor / inputMapFactor + 1;
	    immbpl = (dnx + 7) >> 3;
	    inputMapMask = (unsigned char *) malloc(dny * immbpl);
	    memset(inputMapMask, 0, dny * immbpl);
	    for (y = 0; y < mphi; ++y)
	      for (x = 0; x < mpwi; ++x)
		if (inputMap[y*mpwi+x].c > 0.0)
		  {
		    ixv = x + offxi;
		    iyv = y + offyi;
		    if (ixv >= 0 && ixv < dnx &&
			iyv >= 0 && iyv < dny)
		      inputMapMask[iyv*immbpl+(ixv >> 3)] |= 0x80 >> (ixv & 7);
		  }
	    inputMapDist = (float *) malloc(dny * dnx * sizeof(float));
	    computeDistance(CHESSBOARD_DISTANCE, dnx, dny,
			    inputMapMask, inputMapDist);
	    
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
		  if (ixv < -1 || ixv > dnx ||
		      iyv < -1 || iyv > dny)
		    Error("ixv (%d) or iyv (%d) out-of-range during initial position computation.\n",
			  ixv, iyv);
		  if (ixv < 0)
		    ixv = 0;
		  else if (ixv >= dnx)
		    ixv = dnx-1;
		  if (iyv < 0)
		    iyv = 0;
		  else if (iyv >= dny)
		    iyv = dny-1;
		  threshold = inputMapDist[iyv*dnx+ixv] + 2.0;
		      
		  ixv -= offxi;
		  iyv -= offyi;
		  if (ixv >= 0 && ixv < mpwi-1 &&
		      iyv >= 0 && iyv < mphi-1)
		    {
		      GETMAP(inputMap, mpwi, ixv, iyv, &rx00, &ry00, &rc00);
		      GETMAP(inputMap, mpwi, ixv, iyv+1, &rx01, &ry01, &rc01);
		      GETMAP(inputMap, mpwi, ixv+1, iyv, &rx10, &ry10, &rc10);
		      GETMAP(inputMap, mpwi, ixv+1, iyv+1, &rx11, &ry11, &rc11);
		      if (rc00 > 0.0 &&
			  rc01 > 0.0 &&
			  rc10 > 0.0 &&
			  rc11 > 0.0)
			{
			  rx = rx00 * (rrx - 1.0) * (rry - 1.0)
			    - rx10 * rrx * (rry - 1.0) 
			    - rx01 * (rrx - 1.0) * rry
			    + rx11 * rrx * rry;
			  ry = ry00 * (rrx - 1.0) * (rry - 1.0)
			    - ry10 * rrx * (rry - 1.0) 
			    - ry01 * (rrx - 1.0) * rry
			    + ry11 * rrx * rry;
			}
		      else
			{
			  if (!Extrapolate(&rx, &ry, &rc,
					   ixv, iyv, rrx, rry,
					   inputMap, mpwi, mphi, threshold))
			    Error("Could not extrapolate initial position\n");
			}
		    }
		  else
		    {
		      if (!Extrapolate(&rx, &ry, &rc,
				       ixv, iyv, rrx, rry,
				       inputMap, mpwi, mphi, threshold))
			Error("Could not extrapolate initial position\n");
		    }
		  SETMAP(map, mpw, x, y,
			 rx * inputMapFactor / lFactor,
			 ry * inputMapFactor / lFactor, 1.0);
		  Log("Setting map[%d,%d] to (%f %f)\n",
		      x, y,
		      rx * inputMapFactor / lFactor,
		      ry * inputMapFactor / lFactor);
		}
	    free(inputMapDist);
	    free(inputMapMask);
	    free(inputMap);
	    inputMap = NULL;
	  }
	else if (nCpts >= 4 && c.cptsMethod >= QUADRATIC_METHOD)
	  {
	    double a[6], b[6];
	    double *fx, *fy, *rx, *ry;
	    double xp, yp;

	    fx = (double *) malloc(nCpts * sizeof(double));
	    fy = (double *) malloc(nCpts * sizeof(double));
	    rx = (double *) malloc(nCpts * sizeof(double));
	    ry = (double *) malloc(nCpts * sizeof(double));

	    for (i = 0; i < nCpts; ++i)
	      {
		fx[i] = cpts[i].ix;
		fy[i] = cpts[i].iy;
		rx[i] = cpts[i].rx;
		ry[i] = cpts[i].ry;
	      }
	    ComputeMapping(a, b, fx, fy, rx, ry, nCpts, 6);
	    for (y = 0; y < mph; ++y)
	      {
		yp = lFactor * (y + moy);
		for (x = 0; x < mpw; ++x)
		  {
		    xp = lFactor * (x + mox);
		    SETMAP(map, mpw, x, y,
			   (a[4] * xp * xp + a[5] * yp * yp + a[3] * xp * yp +
			    a[1] * xp + a[2] * yp + a[0]) / lFactor,
			   (b[4] * xp * xp + b[5] * yp * yp + b[3] * xp * yp +
			    b[1] * xp + b[2] * yp + b[0]) / lFactor,
		       1.0);
		  }
	      }
	    Log("Setting initial map to a quadratic map:\n");
	    Log("  x' = %gx^2 %+gy^2 %+gxy %+gx %+gy %+g\n", a[4], a[5], a[3], a[1], a[2], a[0]);
	    Log("  y' = %gx^2 %+gy^2 %+gxy %+gx %+gy %+g\n", b[4], b[5], b[3], b[1], b[2], b[0]);

	    free(fx);
	    free(fy);
	    free(rx);
	    free(ry);
	  }
	else if (nCpts >= 2 && c.cptsMethod >= AFFINE_METHOD)
	  {
	    double a[6], b[6];
	    double *fx, *fy, *rx, *ry;
	    double xp, yp;

	    fx = (double *) malloc(nCpts * sizeof(double));
	    fy = (double *) malloc(nCpts * sizeof(double));
	    rx = (double *) malloc(nCpts * sizeof(double));
	    ry = (double *) malloc(nCpts * sizeof(double));

	    for (i = 0; i < nCpts; ++i)
	      {
		fx[i] = cpts[i].ix;
		fy[i] = cpts[i].iy;
		rx[i] = cpts[i].rx;
		ry[i] = cpts[i].ry;
	      }
	    ComputeMapping(a, b, fx, fy, rx, ry, nCpts, 3);
	    for (y = 0; y < mph; ++y)
	      {
		yp = lFactor * (y + moy);
		for (x = 0; x < mpw; ++x)
		  {
		    xp = lFactor * (x + mox);
		    SETMAP(map, mpw, x, y,
			   (a[1] * xp + a[2] * yp + a[0]) / lFactor,
			   (b[1] * xp + b[2] * yp + b[0]) / lFactor,
			   1.0);
		  }
	      }
	    Log("Setting initial map to an affine map:\n");
	    Log("  x' = %gx %+gy %+g\n", a[1], a[2], a[0]);
	    Log("  y' = %gx %+gy %+g\n", b[1], b[2], b[0]);

	    free(fx);
	    free(fy);
	    free(rx);
	    free(ry);
	  }
	else if (nCpts >= 2 && c.cptsMethod >= RIGID_METHOD)
	  {
	    Log("-rigid is currently unimplemented\n");
	    SetMessage("-rigid is currently unimplemented\n");
	    return;
	  }
	else if (nCpts >= 1 && c.cptsMethod >= TRANSLATION_METHOD)
	  {
	    double dx, dy;
	    
	    dx = dy = 0.0;
	    for (i = 0; i < nCpts; ++i)
	      {
		dx += cpts[i].ix - cpts[i].rx;
		dy += cpts[i].iy - cpts[i].ry;
	      }
	    dx = dx / nCpts;
	    dy = dy / nCpts;

	    /* use a simple translation map */
	    for (y = 0; y < mph; ++y)
	      for (x = 0; x < mpw; ++x)
		SETMAP(map, mpw, x, y,
		       ((float) (x + mox)) + dx / lFactor,
		       ((float) (y + moy)) + dy / lFactor,
		       1.0);
	    Log("Setting initial map to a translation of (%f, %f)\n",
		dx, dy);
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

      nValid = 0;
      for (y = 0; y < mph; ++y)
	for (x = 0; x < mpw; ++x)
	  if (MAP(map, mpw, x, y).c != 0.0)
	    ++nValid;
      Log("Starting with %d of %d valid map points.\n", nValid, mph * mpw);

#if MASKING
      /* trim off any map points that correspond to masked areas of
	 the outputMask or outputPairwiseMask */
      nTrimmed = 0;
      if (outputMasks[cLevel] != NULL)
	{
	  m = outputMasks[cLevel];
	  ombpl = (iw + 7) >> 3;
	  for (y = 0; y < mph; ++y)
	    {
	      iy = (y + moy) * factor - imgoy;
	      for (x = 0; x < mpw; ++x)
		{
		  ix = (x + mox) * factor - imgox;
		  valid = 0;
		  for (dy = -factor; dy < factor && !valid; ++dy) 
		    {
		      iyv = iy + dy;
		      if (iyv < 0 || iyv >= ih)
			continue;
		      for (dx = -factor; dx < factor && !valid; ++dx)
			{
			  ixv = ix + dx;
			  if (ixv < 0 || ixv >= iw)
			    continue;
			  if (MASK(m, ombpl, ixv, iyv))
			    valid = 1;
			}
		    }
		  if (!valid && MAP(map, mpw, x, y).c != 0.0)
		    {
		      ++nTrimmed;
		      MAP(map, mpw, x, y).c = 0.0;
		    }
		}
	    }
	}
      nValid -= nTrimmed;
      Log("Trimmed %d map points in stage 0; %d of %d remaining\n", nTrimmed, nValid, mph * mpw);
#endif

      /* trim off exterior of map */
      /* remove any map point that is farther than sqrt(2)*2^depth pixels
	 from an unmasked image pixel */
      nTrimmed = 0;
      dist = (float *) malloc(((size_t) ih) * iw * sizeof(float));
      if (dist == NULL)
	{
	  SetMessage("Could not allocate dist array (%zd)\n",
		     ((size_t) ih) * iw * sizeof(float));
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
	      if (IMAGE(dist, iw, ix, iy) > threshold && MAP(map, mpw, x, y).c != 0.0)
		{
		  ++nTrimmed;
#if 0
		  if (x == 210 && y == 175)
		    Log("SCB 0: %f %f ix=%d iy=%d mox=%d moy=%d iw=%d ih=%d mpw=%d mph=%d factor=%d imgox=%d imgoy=%d\n",
			IMAGE(dist, iw, ix, iy), threshold,
			ix, iy, mox, moy, iw, ih, mpw, mph,
			factor, imgox, imgoy);
#endif
		  MAP(map, mpw, x, y).c = 0.0;
		}
	    }
	}
      free(dist);
      nValid -= nTrimmed;
      Log("Trimmed %d map points in stage 1; %d of %d remaining\n", nTrimmed, nValid, mph * mpw);

      /* remove any map point that maps to a point farther than 2*2^depth
	 pixels from an unmasked reference pixel */
      nTrimmed = 0;
      border = 2*factor;
      // round border up to a multiple of 8
      if ((border & 7) != 0 || border == 0)
	border = (border + 8) & ~7;
      rpbw = rw + 2 * border;
      rpbh = rh + 2 * border;
      rpbbpl = (rpbw + 7) >> 3;
      rpbmask = (unsigned char *) malloc(rpbh * rpbbpl);
      if (rpbmask == NULL)
	{
	  SetMessage("Could not allocate rpbmask array (%zd)\n",
		     rpbh * rpbbpl);
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
	  SetMessage("Could not allocate dist array (rpbh=%d rpbw=%d)\n",
		     rpbh, rpbw);
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

#if 1
      Log("\nrpbw = %d  rpbh = %d  factor = %d  border = %d\n",
	  rpbw, rpbh, factor, border);
#endif
      for (y = 0; y < mph; ++y)
	for (x = 0; x < mpw; ++x)
	  {
	    if (MAP(map, mpw, x, y).c == 0.0)
	      continue;

	    GETMAP(map, mpw, x, y, &rx00, &ry00, &rc00);
	    xv = factor * rx00 - refox;
	    yv = factor * ry00 - refoy;
	    ixv = ((int) floor(xv)) + border;
	    iyv = ((int) floor(yv)) + border;

	    if (ixv < 0 || ixv >= rpbw ||
		iyv < 0 || iyv >= rpbh ||
		IMAGE(dist, rpbw, ixv, iyv) > threshold &&
		MAP(map, mpw, x, y).c != 0.0)
	      {
		++nTrimmed;
#if 1
		Log("TRIM %d %d %f %f %d %d %f %f\n",
		    x, y, xv, yv, ixv, iyv,
		    ixv >= 0 && ixv < rpbw && iyv >= 0 && iyv < rpbh ?
		    IMAGE(dist, rpbw, ixv, iyv) : -999.0,
		    threshold);
#endif
#if 0
		if (x == 210 && y == 175)
		  Log("SCB 1: %f %f\n",
		      IMAGE(dist, rpbw, ixv, iyv), threshold);
#endif
		MAP(map, mpw, x, y).c = 0.0;
	      }
	  }
      free(dist);
      free(rpbmask);
      nValid -= nTrimmed;
      Log("Trimmed %d map points in stage 2; %d of %d remaining\n", nTrimmed, nValid, mph * mpw);

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
#if 0
		if (x == 210 && y == 175)
		  Log("SCB 2: %d\n", nTrimmed);
#endif
		MAP(map, mpw, x, y).c = 0.0;
	      }
	  }
      nValid -= nTrimmed;
      Log("Trimmed %d map points in stage 3; %d of %d remaining\n", nTrimmed, nValid, mph * mpw);

      /* if insufficient map points are left, exit */
      if (nValid == 0)
	{
	  Log("No valid map points remaining!\n");
	  correlation = 0.0;
	  distortion = 0.0;
	  correspondence = 0.0;
	  constraining = 0.0;
	  goto writeScore;
	}

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

      /* compute the nominal areas and lengths for measuring distortion */
      nomArea = (double *) malloc(mph_minus_1 * mpw_minus_1 * sizeof(double));
      nomL0 = (double *) malloc(mph_minus_1 * mpw_minus_1 * sizeof(double));
      nomL1 = (double *) malloc(mph_minus_1 * mpw_minus_1 * sizeof(double));
      nomL2 = (double *) malloc(mph_minus_1 * mpw_minus_1 * sizeof(double));
      nomL3 = (double *) malloc(mph_minus_1 * mpw_minus_1 * sizeof(double));
      nomThetaX = (double *) malloc(mph_minus_1 * mpw_minus_1 * sizeof(double));
      nomThetaY = (double *) malloc(mph_minus_1 * mpw_minus_1 * sizeof(double));
      for (y = 0; y < mph_minus_1; ++y)
	for (x = 0; x < mpw_minus_1; ++x)
	  {
	    k = y * mpw_minus_1 + x;
	    GETMAP(map, mpw, x, y, &x0, &y0, &c0);
	    GETMAP(map, mpw, x, y+1, &x1, &y1, &c1);
	    GETMAP(map, mpw, x+1, y+1, &x2, &y2, &c2);
	    GETMAP(map, mpw, x+1, y, &x3, &y3, &c3);
	    if (c0 == 0.0 || c1 == 0.0 || c2 == 0.0 || c3 == 0.0)
	      {
		nomArea[k] = 1000000000.0;
		nomL0[k] = 1000000000.0;
		nomL1[k] = 1000000000.0;
		nomL2[k] = 1000000000.0;
		nomL3[k] = 1000000000.0;
		nomThetaX[k] = 0.0;
		nomThetaY[k] = 0.0;
		continue;
	      }
	    nomArea[k] = - 0.5 *((x2 - x0) * (y3 - y1) -
				 (x3 - x1) * (y2 - y0));
	    nomL0[k] = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
	    nomL1[k] = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
	    nomL2[k] = sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2));
	    nomL3[k] = sqrt((x0 - x3) * (x0 - x3) + (y0 - y3) * (y0 - y3));
	    nomThetaX[k] = atan2(y3 - y0,  x3 - x0);
	    nomThetaY[k] = atan2(y1 - y0,  x1 - x0);
	  }

      /* calculate energy of mapping */
      distortion = 0.0;
      correlation = 0.0;
      correspondence = 0.0;
      energy = 0.0;
      Log("energy set to %f\n", energy);

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
	    rx = factor * rx - 0.5 - refox;
	    ry = factor * ry - 0.5 - refoy;
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
	    thetaX = atan2(y3 - y0, x3 - x0);
	    thetaY = atan2(y1 - y0, x1 - x0);
	    
	    k = y * mpw_minus_1 + x;
#if 0
	    de += (area - nomArea[k]) * (area - nomArea[k]) +
	      (l0 - nomL0[k]) * (l0 - nomL0[k]) +
	      (l1 - nomL1[k]) * (l1 - nomL1[k]) +
	      (l2 - nomL2[k]) * (l2 - nomL2[k]) +
	      (l3 - nomL3[k]) * (l3 - nomL3[k]);
#else
	    txd = fabs(fmod(nomThetaX[k] + M_PI - thetaX, 2.0*M_PI) - M_PI);
	    tyd = fabs(fmod(nomThetaY[k] + M_PI - thetaY, 2.0*M_PI) - M_PI);
	    de += txd * txd + tyd * tyd +
	      (l0 - nomL0[k]) * (l0 - nomL0[k]) +
	      (l1 - nomL1[k]) * (l1 - nomL1[k]) +
	      (l2 - nomL2[k]) * (l2 - nomL2[k]) +
	      (l3 - nomL3[k]) * (l3 - nomL3[k]);
#endif
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
		if (ec->c < c.constrainingConfidenceThreshold)
		  continue;
		distance = hypot(e->x - ec->x, e->y - ec->y) - cth;
		if (distance < 0.0)
		  distance = 0.0;
		scc += ec->c;
		constraining += distance;
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
	  SetMessage("Could not allocate prop array (%zd)\n",
		     mSize * sizeof(MapElement));
	  return;
	}
      memset(prop, 0, mSize * sizeof(MapElement));
      
      fflush(stdout);

      multiplier = 1 << level;
      goalMoveCount = (size_t) ceil((c.quality * mpw) * mph * multiplier);
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

      while (moveCount < goalMoveCount)
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

	  moveDebug = 0;
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
		    rx = factor * rx - 0.5 - refox;
		    ry = factor * ry - 0.5 - refoy;
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
		rx = factor * rx - 0.5 - refox;
		ry = factor * ry - 0.5 - refoy;

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

		  k = y * mpw_minus_1 + x;
		  area = - 0.5 *((x2 - x0) * (y3 - y1) -
				 (x3 - x1) * (y2 - y0));
		  l0 = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
		  l1 = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
		  l2 = sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2));
		  l3 = sqrt((x0 - x3) * (x0 - x3) + (y0 - y3) * (y0 - y3));
		  thetaX = atan2(y3 - y0, x3 - x0);
		  thetaY = atan2(y1 - y0, x1 - x0);
		  txd = fabs(fmod(nomThetaX[k] + M_PI - thetaX, 2.0*M_PI) - M_PI);
		  tyd = fabs(fmod(nomThetaY[k] + M_PI - thetaY, 2.0*M_PI) - M_PI);
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
		  nThetaX = atan2(ny3 - ny0, nx3 - nx0);
		  nThetaY = atan2(ny1 - ny0, nx1 - nx0);
		  ntxd = fabs(fmod(nomThetaX[k] + M_PI - nThetaX, 2.0*M_PI) -
			      M_PI);
		  ntyd = fabs(fmod(nomThetaY[k] + M_PI - nThetaY, 2.0*M_PI) -
			      M_PI);
#if 0
		  ode = (area - nomArea[k]) * (area - nomArea[k]) +
		    (l0 - nomL0[k]) * (l0 - nomL0[k]) +
		    (l1 - nomL1[k]) * (l1 - nomL1[k]) +
		    (l2 - nomL2[k]) * (l2 - nomL2[k]) +
		    (l3 - nomL3[k]) * (l3 - nomL3[k]);
		  nde = (narea - nomArea[k]) * (narea - nomArea[k]) +
		    (nl0 - nomL0[k]) * (nl0 - nomL0[k]) +
		    (nl1 - nomL1[k]) * (nl1 - nomL1[k]) +
		    (nl2 - nomL2[k]) * (nl2 - nomL2[k]) +
		    (nl3 - nomL3[k]) * (nl3 - nomL3[k]);
#else
		  ode = txd * txd + tyd * tyd +
		    (l0 - nomL0[k]) * (l0 - nomL0[k]) +
		    (l1 - nomL1[k]) * (l1 - nomL1[k]) +
		    (l2 - nomL2[k]) * (l2 - nomL2[k]) +
		    (l3 - nomL3[k]) * (l3 - nomL3[k]);
		  nde = ntxd * ntxd + ntyd * ntyd +
		    (nl0 - nomL0[k]) * (nl0 - nomL0[k]) +
		    (nl1 - nomL1[k]) * (nl1 - nomL1[k]) +
		    (nl2 - nomL2[k]) * (nl2 - nomL2[k]) +
		    (nl3 - nomL3[k]) * (nl3 - nomL3[k]);
#endif
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
			  if (ec->c < c.constrainingConfidenceThreshold)
			    continue;
			  ep = &MAP(prop, mpw, x, y);
			  distance = hypot(e->x - ec->x, e->y - ec->y) - cth;
			  if (distance < 0.0)
			    distance = 0.0;
			  oldEnergy = distance;
			  distance = hypot(ep->x - ec->x, ep->y - ec->y) - cth;
			  if (distance < 0.0)
			    distance = 0.0;
			  de += distance - oldEnergy;
#if 0
			  if (x * lFactor == 32*64 && y * lFactor == 160*64)
			    {
			      moveDebug = 1;
			      Log("cons constr move: old = %f new = %f de = %f  cct = %f mpw = %d mph = %d cth = %f lFactor = %d\n",
				  oldEnergy, distance, de,
				  c.constrainingThreshold, mpw, mph,
				  cth, lFactor);
			    }
#endif
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
#if 0
		  if (x == 210 && y == 175)
		    Log("SCB 3: %d\n", accept);
#endif
		  MAP(prop, mpw, x, y).c = 0;
		  if (accept)
		    {
		      /* make the provisional state the new state */
		      GETMAP(prop, mpw, x, y, &rx00, &ry00, &rc00);
		      rc00 = 1.0;
		      SETMAP(map, mpw, x, y, rx00, ry00, rc00);
		    }
		}
#if 0
	  if (moveDebug)
	    {
	      Log("%s move %d since deltaE = %f\n",
		  accept ? "Accepted" : "Rejected",
		  moveCount, newEnergy - energy,
		      (newDistortion - distortion) * c.distortion,
		      newCorrelation - correlation,
		      (newCorrespondence - correspondence) * c.correspondence,
		      (newConstraining - constraining) * c.constraining);
	      Log("  old energy = %f (%f * dist %f - correl %f + %f * corres %f + %f * constr %f)\n",
		  energy,
		  c.distortion, distortion,
		  correlation,
		  c.correspondence, correspondence,
		  c.constraining, constraining);
	      Log("  new energy = %f (%f * dist %f - correl %f + %f * corres %f + %f * constr %f)\n",
		  newEnergy,
		  c.distortion, newDistortion,
		  newCorrelation,
		  c.correspondence, newCorrespondence,
		  c.constraining, newConstraining);
	      Log("  cx = %f cy = %f radius = %f offset = %f theta = %f\n",
		  cx, cy, radius, offset, theta);
	    }
#endif

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
	      constraining = newConstraining;
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
      free(nomArea);
      free(nomL0);
      free(nomL1);
      free(nomL2);
      free(nomL3);
      free(nomThetaX);
      free(nomThetaY);
#if FOLDING
      free(cirdisc);
#endif

      if (level == c.outputLevel || c.writeAllMaps)
	{
	  Log("Going to write output map at level %d\n", level);
	  count = 0;
	  for (y = 0; y < mph; ++y)
	    for (x = 0; x < mpw; ++x)
	      if (MAP(map, mpw, x, y).c >= 0.001)
		++count;
	  Log("map initially has %d valid points of %d possible\n", count, mph * mpw);
	  Log("mox = %d moy = %d imgox = %d imgoy = %d\n",
	      mox, moy, imgox, imgoy);
#if 1
	  warpedArray = (float *) malloc(ih * iw * sizeof(float));
	  validArray = (unsigned char *) malloc(iw * ih *
						sizeof(unsigned char));
	  ComputeWarpedImage(warpedArray, validArray,
			     iw, ih,
			     imgox, imgoy,
			     cref, crmask,
			     rw, rh,
			     refox, refoy,
			     map,
			     factor, mpw, mph, mox, moy);

	  correlationArray = (float *) malloc(ih * iw * sizeof(float));
	  ComputeCorrelation(correlationArray,
			     cimage, warpedArray, validArray,
			     iw, ih,
			     c.correlationHalfWidth);

	  // use the correlation as the confidence values
	  //   for the map entries
	  count1 = count2 = count3 = count4 = 0;
	  for (y = 0; y < mph; ++y)
	    for (x = 0; x < mpw; ++x)
	      {
		ixv = factor * (x + mox) - imgox;
		iyv = factor * (y + moy) - imgoy;
		if (ixv < 0 || ixv >= iw ||
		    iyv < 0 || iyv >= ih)
		  {
		    if (MAP(map, mpw, x, y).c >= 0.001)
		      {
			++count1;
			MAP(map, mpw, x, y).c = 0.001;
		      }
		    else
		      {
			++count2;
			MAP(map, mpw, x, y).c = 0.0;
		      }
		  }
		else
		  {
		    if (MAP(map, mpw, x, y).c >= 0.001)
		      {
			++count3;
			MAP(map, mpw, x, y).c = correlationArray[iyv * iw + ixv];
			if (MAP(map, mpw, x, y).c == 0.0)
			  MAP(map, mpw, x, y).c = 0.001;
		      }
		    else
		      {
			++count4;
			MAP(map, mpw, x, y).c = 0.0;
		      }
		  }
	      }
	  Log("map to be output stats: %d %d %d %d\n", count1, count2, count3, count4);
#endif

	  if (level == c.outputLevel)
	    TrimOutputMap(map, mpw, mph, mox, moy,
			  factor,
			  iw, ih, imgox, imgoy,
			  rw, rh, refox, refoy,
			  cimask,
			  crmask);

	  WriteOutputMap(outputName, level, map, mpw, mph, mox, moy);
	  if (outputWarpedName[0] != '\0')
	    WriteOutputImage(outputWarpedName, c.outputLevel, warpedArray,
			     iw, ih, 1.0);
	  if (outputCorrelationName[0] != '\0')
	    WriteOutputImage(outputCorrelationName, c.outputLevel, correlationArray,
			     iw, ih, 256.0);
	  
#if 1
	  free(warpedArray);
	  free(validArray);
	  free(correlationArray);
#endif
	}
    }

 writeScore:  
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

  r.correlation = correlation;
  r.distortion = distortion;
  r.correspondence = correspondence;
  r.constraining = constraining;
}


void
TrimOutputMap (MapElement *map, int mpw, int mph, int mox, int moy,
	       int factor,
	       unsigned int iw, unsigned int ih, int imgox, int imgoy,
	       unsigned int rw, unsigned int rh, int refox, int refoy,
	       unsigned char *cimask,
	       unsigned char *crmask)
{
  int x, y;
  int ixv, iyv;
  int ix, iy;
  float total;
  int nInside;
  unsigned char *validArray;
  size_t icmbpl, rcmbpl;
  float rx00, rx01, rx10, rx11;
  float ry00, ry01, ry10, ry11;
  float rc00, rc01, rc10, rc11;
  float x00, y00, x01, y01, x11, y11, x10, y10;
  float ix00, iy00, ix01, iy01, ix11, iy11, ix10, iy10;
  int minX, maxX, minY, maxY;
  float xv, yv;
  int inside0, inside1, inside2, inside3;

  // the validArray contains a 1 for every valid map quadrilateral;
  //  note the the mph-1 row and mpw-1 column are unused, but
  //  retained here for consistency
  validArray = (unsigned char *) malloc(mpw * mph * sizeof(unsigned char));
  memset(validArray, 0x00, mpw * mph * sizeof(unsigned char));
  for (y = 0; y < mph-1; ++y)
    for (x = 0; x < mpw-1; ++x)
      MAP(validArray, mpw, x, y) = 1;

  // first trim out map quadrilaterals that have previously been trimmed
  for (y = 0; y < mph; ++y)
    for (x = 0; x < mpw; ++x)
      if (MAP(map, mpw, x, y).c == 0.0)
	{
	  MAP(validArray, mpw, x, y) = 0;
	  if (x > 0)
	    MAP(validArray, mpw, x-1, y) = 0;
	  if (y > 0)
	    MAP(validArray, mpw, x, y-1) = 0;
	  if (x > 0 && y > 0)
	    MAP(validArray, mpw, x-1, y-1) = 0;
	}

  // then trim map quadrilaterals that don't contain sufficient valid pixels in
  //   the source image
  icmbpl = (iw + 7) >> 3;
  for (y = 0; y < mph-1; ++y)
    {
      iyv = factor * (y + moy) - imgoy;
      for (x = 0; x < mpw-1; ++x)
	{
	  if (MAP(validArray, mpw, x, y) == 0)
	    continue;
	  total = 0.0;
	  ixv = factor * (x + mox) - imgox;
	  for (iy = iyv; iy < iyv+factor; ++iy)
	    {
	      if (iy < 0 || iy >= ih)
		continue;
	      for (ix = ixv; ix < ixv+factor; ++ix)
		{
		  if (ix < 0 || ix >= iw)
		    continue;
		  if ((cimask[iy*icmbpl + (ix >> 3)] & (0x80 >> (ix & 7))) == 0)
		    continue;
		  total += 1.0;
		}
	    }
	  if (100.0 * total / (factor*factor) <= c.trimMapSourceThreshold)
	    MAP(validArray, mpw, x, y) = 0;
	}
    }

  rcmbpl = (rw + 7) >> 3;
  for (y = 0; y < mph-1; ++y)
    {
      iyv = factor * (y + moy) - imgoy;
      for (x = 0; x < mpw-1; ++x)
	{
	  if (MAP(validArray, mpw, x, y) == 0)
	    continue;

	  // get map element corners in clockwise order
	  GETMAP(map, mpw, x, y, &rx00, &ry00, &rc00);
	  GETMAP(map, mpw, x, y+1, &rx01, &ry01, &rc01);
	  GETMAP(map, mpw, x+1, y+1, &rx11, &ry11, &rc11);
	  GETMAP(map, mpw, x+1, y, &rx10, &ry10, &rc10);
	  x00 = factor * rx00;
	  y00 = factor * ry00;
	  x01 = factor * rx01;
	  y01 = factor * ry01;
	  x11 = factor * rx11;
	  y11 = factor * ry11;
	  x10 = factor * rx10;
	  y10 = factor * ry10;
	  ix00 = (int) floor(x00 - refox);
	  iy00 = (int) floor(y00 - refoy);
	  ix01 = (int) floor(x01 - refox);
	  iy01 = (int) floor(y01 - refoy);
	  ix11 = (int) floor(x11 - refox);
	  iy11 = (int) floor(y11 - refoy);
	  ix10 = (int) floor(x10 - refox);
	  iy10 = (int) floor(y10 - refoy);
	  ix00 = (int) floor(x00 - refox);
	  iy00 = (int) floor(y00 - refoy);
	  ix01 = (int) floor(x01 - refox);
	  iy01 = (int) floor(y01 - refoy);
	  ix11 = (int) floor(x11 - refox);
	  iy11 = (int) floor(y11 - refoy);
	  ix10 = (int) floor(x10 - refox);
	  iy10 = (int) floor(y10 - refoy);
	  
	  // get limits in reference space
	  minX = ix00;
	  maxX = ix00;
	  minY = iy00;
	  maxY = iy00;
	  if (ix01 < minX)
	    minX = ix01;
	  if (ix01 > maxX)
	    maxX = ix01;
	  if (iy01 < minY)
	    minY = iy01;
	  if (iy01 > maxY)
	    maxY = iy01;
	  if (ix11 < minX)
	    minX = ix11;
	  if (ix11 > maxX)
	    maxX = ix11;
	  if (iy11 < minY)
	    minY = iy11;
	  if (iy11 > maxY)
	    maxY = iy11;
	  if (ix10 < minX)
	    minX = ix10;
	  if (ix10 > maxX)
	    maxX = ix10;
	  if (iy10 < minY)
	    minY = iy10;
	  if (iy10 > maxY)
	    maxY = iy10;
	  
	  // go through all pixels
	  total = 0.0;
	  nInside = 0;
	  for (iy = minY; iy <= maxY; ++iy)
	    for (ix = minX; ix <= maxX; ++ix)
	      {
		// find pixel center
		xv = ix + refox + 0.5;
		yv = iy + refoy + 0.5;

		// check if it is inside map element quadrilateral
		inside0 = (x01 - x00)*(yv - y00) - (y01 - y00)*(xv - x00) < 0.0;
		inside1 = (x11 - x01)*(yv - y01) - (y11 - y01)*(xv - x01) < 0.0;
		inside2 = (x10 - x11)*(yv - y11) - (y10 - y11)*(xv - x11) < 0.0;
		inside3 = (x00 - x10)*(yv - y10) - (y00 - y10)*(xv - x10) < 0.0;
		if (!inside0 || !inside1 || !inside2 || !inside3)
		  continue;
		++nInside;

		// check if pixel is marked as valid
		if (ix < 0 || ix >= rw || iy < 0 || iy >= rh)
		  continue;
		if ((crmask[iy*rcmbpl + (ix >> 3)] & (0x80 >> (ix & 7))) == 0)
		    continue;
		  
		total += 1.0;
	      }
	}
      if (nInside == 0 || 100.0 * total / nInside <= c.trimMapTargetThreshold)
	MAP(validArray, mpw, x, y) = 0;
    }

  // mark map nodes as invalid if no quadrilateral containing them is valid
#if 0
  {
    char out[4096];

    Log("VALID ARRAY:\n");
    for (y = 0; y < mph; ++y)
      {
	memset(out, 0, 4096);
	for (x = 0; x < mpw; ++x)
	  sprintf(&out[2*x], "%d ", MAP(validArray, mpw, x, y));
	Log("    %s\n", out);
      }
  }
#endif

  for (y = 0; y < mph; ++y)
    for (x = 0; x < mpw; ++x)
      if (!(x < mpw-1 && y < mph-1 && MAP(validArray, mpw, x, y) ||
	    x > 0 && y < mph-1 && MAP(validArray, mpw, x-1, y) ||
	    x < mpw-1 && y > 0 && MAP(validArray, mpw, x, y-1) ||
	    x > 0 && y > 0 && MAP(validArray, mpw, x-1, y-1)))
	MAP(map, mpw, x, y).c = 0.0;

#if 0
  {
    char out[4096];
    int len;

    Log("CONF ARRAY:\n");
    for (y = 0; y < mph; ++y)
      {
	memset(out, 0, 4096);
	len = 0;
	for (x = 0; x < mpw; ++x)
	  len += sprintf(&out[len], "%f ", MAP(map, mpw, x, y).c);
	Log("    %s\n", out);
      }
  }
#endif

  free(validArray);
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
  char msg[PATH_MAX+256];
  int count;

  count = 0;
  for (y = 0; y < mph; ++y)
    for (x = 0; x < mpw; ++x)
      if (MAP(map, mpw, x, y).c != 0.0)
	++count;
  Log("WriteOutputMap: %d of %d\n", count, mph * mpw);

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

  outMap = (MapElement*) malloc(nx * ny * sizeof(MapElement));
  if (outMap == NULL)
    {
      SetMessage("Could not allocate outMap array (%zd)\n",
		 ((size_t) nx) * ny * sizeof(MapElement));
      return;
    }
  for (y = 0; y < ny; ++y)
    for (x = 0; x < nx; ++x)
      {
	GETMAP(map, mpw, minX + x, minY + y, &rx, &ry, &rc);
	SETMAP(outMap, nx, x, y, rx, ry, rc);
    }
  if (!WriteMap(fn, outMap, level,
		nx, ny,
		minX + mox, minY + moy,
		t.pair.imageName[0], t.pair.imageName[1],
		UncompressedMap,
		msg))
    Error("Could not write output map %s :\n%s\n",
	  fn, msg);

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
    sprintf(fn, "%s.%s", outputName, extension);
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
		    int w, int h,           /* of the warped reference */
		    int imgox, int imgoy,   /* of the warped reference */
		    float *image, unsigned char *mask,
		    int iw, int ih,         /* of the reference (image & mask) */
		    int refox, int refoy,   /* of the reference (image & mask) */
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
  unsigned char msg[PATH_MAX];
  size_t count1, count2, count3, count4, count5;

  mbpl = (iw + 7) >> 3;
  memset(valid, 0, w * h * sizeof(unsigned char));

  count1 = count2 = count3 = count4 = count5 = 0;
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
	    ++count1;
	    warped[y * w + x] = 0.0;
	    continue;
	  }

	GETMAP(map, mpw, ixv, iyv, &rx00, &ry00, &rc00);
	GETMAP(map, mpw, ixv, iyv+1, &rx01, &ry01, &rc01);
	GETMAP(map, mpw, ixv+1, iyv, &rx10, &ry10, &rc10);
	GETMAP(map, mpw, ixv+1, iyv+1, &rx11, &ry11, &rc11);

	if (rc00 == 0.0 || rc01 == 0.0 || rc10 == 0.0 || rc11 == 0.0)
	  {
	    ++count2;
	    warped[y * w + x] = 0.0;
	    continue;
	  }

	rx = rx00 * (rrx - 1.0) * (rry - 1.0)
	  - rx10 * rrx * (rry - 1.0)
	  - rx01 * (rrx - 1.0) * rry
	  + rx11 * rrx * rry;
	ry = ry00 * (rrx - 1.0) * (rry - 1.0)
	  - ry10 * rrx * (rry - 1.0)
	  - ry01 * (rrx - 1.0) * rry
	  + ry11 * rrx * rry;
	rx = mapFactor * rx - 0.5 - refox;
	ry = mapFactor * ry - 0.5 - refoy;
	irx = ((int) floor(rx + 1.0)) - 1;
	iry = ((int) floor(ry + 1.0)) - 1;
	if (irx < 0 || irx >= iw - 1 ||
	    iry < 0 || iry >= ih - 1)
	  {
	    ++count3;
	    warped[y * w + x] = 0.0;
	    continue;
	  }
	rrx = rx - irx;
	rry = ry - iry;

#if MASKING
	if ((mask[iry*mbpl + (irx >> 3)] & (0x80 >> (irx & 7))) == 0 ||
	    rrx > 0.0 && (mask[iry*mbpl + ((irx+1) >> 3)] & (0x80 >> ((irx+1) & 7))) == 0 ||
	    rry > 0.0 && (mask[(iry+1)*mbpl + (irx >> 3)] & (0x80 >> (irx & 7))) == 0 ||
	    rrx > 0.0 && rry > 0.0 && (mask[(iry+1)*mbpl +((irx+1) >> 3)] & (0x80 >> ((irx+1) & 7))) == 0)
	  {
	    ++count4;
	    warped[y * w + x] = 0.0;
	    continue;
	  }
#endif

	r00 = image[iry * iw + irx];
	r01 = image[(iry + 1) * iw + irx];
	r10 = image[iry * iw + (irx + 1)];
	r11 = image[(iry + 1) * iw + irx + 1];
	rv = r00 * (rrx - 1.0) * (rry - 1.0)
	  - r10 * rrx * (rry - 1.0)
	  - r01 * (rrx - 1.0) * rry
	  + r11 * rrx * rry;
	++count5;
	warped[y * w + x] = rv;
	valid[y * w + x] = 1;
      }
  Log("ComputeWarpedImage: %zd %zd %zd %zd %zd\n",
      count1, count2, count3, count4, count5);
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
  double checkSumA, checkSumA2, checkSumB, checkSumB2, checkSumAB;
  int checkN;
  int dx, dy, cx, cy;
  double checkCorr;
  int minimumSamples;
  size_t count1, count2, count3;

  // find the extents of one quadrant of a circle of the desired half-width;
  // also, let the minimum number of samples be about one-quarter of the
  //   the potential samples
  lim = (int*) malloc((hw+1) * sizeof(int));
  minimumSamples = 1;
  for (x = 0; x <= hw; ++x)
    {
      lim[x] = (int) floor(sqrt((hw + 0.5) * (hw + 0.5) - x * x));
      if (x > 0)
	minimumSamples += (lim[x] + 1);
    }
  Log("ComputeCorrelation: minimumSamples = %d\n", minimumSamples);

  // first compute the sums for (-1, -1)
  startSumA = 0.0;
  startSumB = 0.0;
  startSumA2 = 0.0;
  startSumB2 = 0.0;
  startSumAB = 0.0;
  startN = 0;
  for (y = 0; y < hw; ++y)
    for (x = 0; x < lim[y+1]; ++x)
      {
	if (x >= w || y >= h || !valid[y*w+x])
	  continue;
	av = a[y*w+x];
	bv = b[y*w+x];
	startSumA += av;
	startSumA2 += av * av;
	startSumB += bv;
	startSumB2 += bv * bv;
	startSumAB += av * bv;
	++startN;
	//	printf("INCR %0.4d %0.4d ADD\n", x, y);
      }
  count1 = count2 = count3 = 0;
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
	      //	      printf("INCR %0.4d %0.4d DROP\n", x, y);
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
	      //	      printf("INCR %0.4d %0.4d ADD\n", x, y);
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
#if 0
	  checkSumA = 0.0;
	  checkSumA2 = 0.0;
	  checkSumB = 0.0;
	  checkSumB2 = 0.0;
	  checkSumAB = 0.0;
	  checkN = 0;
	  for (dy = -31; dy <= 31; ++dy)
	    for (dx = -31; dx <= 31; ++dx)
	      {
		cx = xc + dx;
		cy = yc + dy;
		if (cx < 0 || cx >= w ||
		    cy < 0 || cy >= h)
		  continue;
		if ((dx * dx + dy * dy) > ((hw + 0.5) * (hw + 0.5)))
		  continue;
		if (!valid[cy*w+cx])
		  continue;
		av = a[cy*w+cx];
		bv = b[cy*w+cx];
		checkSumA += av;
		checkSumA2 += av * av;
		checkSumB += bv;
		checkSumB2 += bv * bv;
		checkSumAB += av * bv;
		//		printf("CHECK %0.4d %0.4d ADD\n", cx, cy);
		++checkN;
	      }

	  if (checkN >= minimumSamples)
	    {
	      meanA = checkSumA / checkN;
	      meanB = checkSumB / checkN;
	      denom = (checkSumA2 - 2.0 * meanA * checkSumA + checkN * meanA * meanA) *
		(checkSumB2 - 2.0 * meanB * checkSumB + checkN * meanB * meanB);
	      if (denom < 0.001)
		checkCorr = 0.0;
	      else
		checkCorr = (checkSumAB - meanA * checkSumB - meanB * checkSumA + checkN * meanA * meanB) / sqrt(denom);
	    }
	  else
	    checkCorr = 0.0;
#endif
#if 0
	  sumA = checkSumA;
	  sumA2 = checkSumA2;
	  sumB = checkSumB;
	  sumB2 = checkSumB2;
	  sumAB = checkSumAB;
	  n = checkN;
#endif

#if 1
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
	      //	      printf("INCR %0.4d %0.4d DROP\n", x, y);
	    }

	  // add new right
	  x = xc + lim[0];
	  if (x >= 0 && x < w && valid[yc*w+x])
	    {
	      av = a[yc*w+x];
	      bv = b[yc*w+x];
	      sumA += av;
	      sumA2 += av * av;
	      sumB += bv;
	      sumB2 += bv * bv;
	      sumAB += av * bv;
	      ++n;
	      //	      printf("INCR %0.4d %0.4d ADD\n", x, y);
	    }

	  for (i = 1; i <= hw; ++i)
	    {
	      // subtract old left
	      x = xc - lim[i] - 1;
	      if (x >= 0)
		{
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
		      //		      printf("INCR %0.4d %0.4d DROP\n", x, y);
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
		      //		      printf("INCR %0.4d %0.4d DROP\n", x, y);
		    }
		}

	      // add new right
	      x = xc + lim[i];
	      if (x < w)
		{
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
		      //		      printf("INCR %0.4d %0.4d ADD\n", x, y);
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
		      //		      printf("INCR %0.4d %0.4d ADD\n", x, y);
		    }
		}
	    }
#endif

	  // calculate the correlation for the disk surrounding (xc, yc)
	  if (n >= minimumSamples)
	    {
	      meanA = sumA / n;
	      meanB = sumB / n;
	      denom = (sumA2 - 2.0 * meanA * sumA + n * meanA * meanA) *
		(sumB2 - 2.0 * meanB * sumB + n * meanB * meanB);
	      if (denom < 0.001)
		{
		  correlation[yc*w+xc] = 0.0;
		  ++count1;
		}
	      else
		{
		  correlation[yc*w+xc] = (sumAB - meanA * sumB - meanB * sumA + n * meanA * meanB) / sqrt(denom);
		  ++count2;
		}
	    }
	  else
	    {
	      correlation[yc*w+xc] = 0.0;
	      ++count3;
	    }

#if 0
	  if (fabs(correlation[yc*w+xc] - checkCorr) / checkCorr > 0.001)
	    Error("mismatch at %d %d %f %f\n(%f %f %f %f %f %d)\n(%f %f %f %f %f %d)\n",
		  xc, yc, checkCorr, correlation[yc*w+xc],
		  checkSumA, checkSumA2, checkSumB, checkSumB2, checkSumAB, checkN,
		  sumA, sumA2, sumB, sumB2, sumAB, n);
#endif
	}
    }
  Log("ComputeCorrelation: %zd %zd %zd\n", count1, count2, count3);
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
  par_pkstr(c.outputWarpedBasename);
  par_pkstr(c.outputCorrelationBasename);
  par_pkstr(c.outputMaskBasename);
  par_pkstr(c.outputPairwiseMaskBasename);

  par_pkstr(c.logBasename);

  par_pkint(c.startLevel);

  par_pkint(c.outputLevel);
  par_pkint(c.minResolution);
  par_pkint(c.depth);
  par_pkint(c.cptsMethod);
  par_pkdouble(c.distortion);
  par_pkdouble(c.correspondence);
  par_pkdouble(c.correspondenceThreshold);
  par_pkdouble(c.constraining);
  par_pkdouble(c.constrainingThreshold);
  par_pkdouble(c.constrainingConfidenceThreshold);
  par_pkdouble(c.quality);
  par_pkdouble(c.minOverlap);
  par_pkdouble(c.trimMapSourceThreshold);
  par_pkdouble(c.trimMapTargetThreshold);
  par_pkint(c.correlationHalfWidth);
  par_pkint(c.writeAllMaps);
  par_pkint(c.update);
  par_pkint(c.partial);
  par_pkint(c.nWorkers);
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
  par_upkstr(c.outputWarpedBasename);
  par_upkstr(c.outputCorrelationBasename);
  par_upkstr(c.outputMaskBasename);
  par_upkstr(c.outputPairwiseMaskBasename);

  par_upkstr(c.logBasename);

  c.startLevel = par_upkint();
  c.outputLevel = par_upkint();
  c.minResolution = par_upkint();
  c.depth = par_upkint();
  c.cptsMethod = par_upkint();
  c.distortion = par_upkdouble();
  c.correspondence = par_upkdouble();
  c.correspondenceThreshold = par_upkdouble();
  c.constraining = par_upkdouble();
  c.constrainingThreshold = par_upkdouble();
  c.constrainingConfidenceThreshold = par_upkdouble();
  c.quality = par_upkdouble();
  c.minOverlap = par_upkdouble();
  c.trimMapSourceThreshold = par_upkdouble();
  c.trimMapTargetThreshold = par_upkdouble();

  c.correlationHalfWidth = par_upkint();
  c.writeAllMaps = par_upkint();
  c.update = par_upkint();
  c.partial = par_upkint();
  c.nWorkers = par_upkint();
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
      par_upkstr(s);
      CopyString(&(p->imageName[imi]), s);
      
      p->imageMinX[imi] = par_upkint();
      p->imageMaxX[imi] = par_upkint();
      p->imageMinY[imi] = par_upkint();
      p->imageMaxY[imi] = par_upkint();
    }

  par_upkstr(s);
  CopyString(&(p->pairName), s);
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
  if (r.message != NULL)
    par_pkstr(r.message);
  else
    par_pkstr("");
}

void
UnpackResult ()
{
  char s[PATH_MAX];

  UnpackPair(&(r.pair));
  r.updated = par_upkint();
  r.distortion = par_upkdouble();
  r.correlation = par_upkdouble();
  r.correspondence = par_upkdouble();
  r.constraining = par_upkdouble();
  par_upkstr(s);
  if (s[0] != '\0')
    CopyString(&(r.message), s);
  else
    CopyString(&(r.message), NULL);
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

size_t
CountIntersectionBits (unsigned char *p, unsigned char *q, size_t n)
{
  size_t ii;
  unsigned char b;
  size_t sum = 0;

  for (ii = 0; ii < n; ++ii)
    {
      b = (*p++) & (*q++);
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
      
      if (mkdir(dn, 0777) != 0)
	{
	  Log("Could not create directory %s\n", dn);
	  return(0);
	}
    }
  return(1);
}

void
CopyString (char **dst, char *src)
{
  if (*dst != NULL)
    {
      free(*dst);
      *dst = NULL;
    }
  if (src != NULL)
    {
      *dst = (char *) malloc(strlen(src) + 1);
      strcpy(*dst, src);
    }
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

void
SetMessage (char *fmt, ...)
{
  va_list args;
  char s[PATH_MAX];

  va_start(args, fmt);
  vsnprintf(s, PATH_MAX, fmt, args);
  va_end(args);
  CopyString(&(r.message), s);
}

void Error (char *fmt, ...)
{
  va_list args;
  char timestamp[32];

  if (logFile != NULL)
    {
      va_start(args, fmt);
      fprintf(logFile, "%s: ERROR: ", GetTimestamp(timestamp, 32));
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
  char timestamp[32];

  if (c.logBasename[0] == '\0')
    return;

  if (logFile == NULL)
    {
      int i;
      char logName[PATH_MAX];
      i = par_instance();
      if (i < 0)
	sprintf(logName, "%s0.log", c.logBasename);
      else if (c.nWorkers < 10)
	sprintf(logName, "%s%0.1d.log", c.logBasename, i);
      else if (c.nWorkers < 100)
	sprintf(logName, "%s%0.2d.log", c.logBasename, i);
      else if (c.nWorkers < 1000)
	sprintf(logName, "%s%0.3d.log", c.logBasename, i);
      else
	sprintf(logName, "%s%0.6d.log", c.logBasename, i);
      logFile = fopen(logName, "w");
      if (logFile == NULL)
	Error("Could not open log file %s\n", logName);
    }

  va_start(args, fmt);
  fprintf(logFile, "%s: ", GetTimestamp(timestamp, 32));
  vfprintf(logFile, fmt, args);
  va_end(args);
  fflush(logFile);
}

char *
GetTimestamp (char *timestamp, size_t size)
{
  struct timeval current_time;
  struct tm lt;
  size_t len;

  if (size < 24 ||
      gettimeofday(&current_time, NULL) != 0 ||
      strftime(timestamp, size, "%y%m%d %H:%M:%S",
	       localtime_r(&(current_time.tv_sec), &lt)) == 0)
    return(NULL);
  len = strlen(timestamp);
  if (size - len > 8)
    sprintf(&timestamp[len], ".%0.6d", (int) current_time.tv_usec);
  return(timestamp);
}
