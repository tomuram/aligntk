/*
 * align12.c  - align a set of maps using a relaxation algorithm
 *                on a system of springs
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
 *    2009  Development of multiresolution algorithm (ghood@psc.edu)
 */

#include <stdio.h>
#include <mpi.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <sched.h>

#include "imio.h"
#include "dt.h"

#define DEBUG	0
#define PDEBUG	0

#define LINE_LENGTH	256
#define MAX_SPEC        1024
#define CONSTRAINED	(1.0e+30)

typedef struct Step
{
  float level;      // actually an int, but we use a float to simplify
                    //  the MPI_Bcast
  float kIntra;
  float kInter;
  float kAbsolute;
  float threshold;
  float dampingFactor;
} Step;

typedef struct ModelPoint
{
  float x, y;
} ModelPoint;

typedef struct Point
{
  float x;    /* current x position of this Point */
  float y;    /* current y position of this Point */
  float fx;   /* force in x direction;
		 if fx > 0.5*CONSTRAINED, Point is not allowed to
		 move in the x direction */
  float fy;   /* force in y direction;
		 if fy > 0.5*CONSTRAINED, Point is not allowed to
		 move in the y direction */
} Point;

typedef struct Image
{
  int next;		/* index of next image in hash bucket */
  char *name;   	/* name of this image */
  int width, height;    /* width and height in pixels */
  float rotation, scale;/* initial rotation, scale, and translation */
  float tx, ty;         /*    to be applied to image to obtain starting */
                        /*    positions (expressed in pixels) */
  int fixed;		/* true if this image is fixed in position */
  int owner;		/* which process owns the image */
  int needed;		/* true if this image is needed in this process */ 
  unsigned char *sendTo;/* bitmask of the processes that this image should
			   be sent to */
  int nx, ny;           /* number of points in x and y at current level */
  Point **points;       /* points in the image at each level */
  char *modelName;      /* name of map to use for initial positions and
			   distances */
  struct IntraImageMap *map;
                        /* map to use for intra-image springs */
} Image;

typedef struct IntraImageMap
{
  char *name;           /* filename of map (or "" if uniform map) */
  struct IntraImageMap *next;  /* next in hash bucket */
  int nx, ny;           /* number of points in x and y at the ending level */
  ModelPoint *points;   /* points in the model at the starting level */
  int *nSprings;        /* number of springs at each level */
  struct IntraImageSpring **springs;
                        /* array of springs at each level */
} IntraImageMap;

typedef struct IntraImageSpring
{
  int index0;		/* index of first point in image */
  int index1;		/* index of second point in image */
  float nomD;           /* nominal distance between the points */
  float k;		/* spring constant for this entry */
} IntraImageSpring;

typedef struct InterImageMap
{
  char *name;           /* filename of map */
  int image0;           /* index of source image */
  int image1;           /* index of target (reference) image */
  float energyFactor;   /* either 0.0 or 1.0; determines whether this entry
			   contributes to the energy sum on this node */
  float k;		/* overall spring constant for this map */

  int *nStrips;          /* number of strips at each level */
  struct InterImageStrip **strips;     /* all strips at each level */
  int *nSprings;         /* number of springs at each level */
  struct InterImageSpring **springs;  /* all springs at each level */
} InterImageMap;

typedef struct InterImageStrip
{
  unsigned short nSprings; /* number of springs in this strip */
  unsigned short x0, y0;   /* source position of first spring */
  unsigned short x1, y1;   /* target position of first spring */
} InterImageStrip;

typedef struct InterImageSpring
{
  short dx, dy;         /* pixel offset for distance computations;
			   this is in units of 0.0001 * nominal spring spacing */
  unsigned char k;      /* spring constant (multiplied by 100) */
  unsigned char dxy1;   /* delta position in target:
			     delta_x = (dxy1 >> 4) - 8;
			     delta_y = (dxy1 & 0xf) - 8; */
  short irrx, irry;     /* ratio of goal offset to distances to adjacent points
			   in target map (multiplied by 65536); should be
			   in range -32768 to 32767 (corresponding to a ratio
			   of -0.5 to 0.5) */
} InterImageSpring;

typedef struct CommPhase
{
  int otherProcess;      /* the process to communicate with in this phase;
			    the lower-numbered process always sends first,
			    then receives */
  int nSendImages;	 /* number of images to send */
  int *sendImages;	 /* the images to send */
  int nSends;            /* number of MPI_Send's */
  int *nImagesToSend;    /* number of images to send for each MPI_Send */
  int nReceiveImages;    /* number of images to receive */
  int *receiveImages;    /* the images to receive */
  int nReceives;	 /* number of MPI_Recv's */
  int *floatsToReceive;  /* number of floats to receive for each MPI_Recv */
  int *nImagesToReceive; /* number of images to receive for each MPI_Recv */
} CommPhase;


int p;      /* rank of this process */
int np;     /* number of processes in this run */

char schedule[PATH_MAX];
int params[12];
float fParams[5];
char imageListName[PATH_MAX];
char mapListName[PATH_MAX];
char mapsName[PATH_MAX];
char initialMapsName[PATH_MAX];
char modelsName[PATH_MAX];
char outputName[PATH_MAX];
int outputIncremental = 0;
char outputGridName[PATH_MAX];
int outputGridWidth = 1024;
int outputGridHeight = 1024;
int outputGridOverlay = 0;
int outputGridSprings = 0;
int outputGridInterval = -1;
int outputGridFocusWidth = -1;
int outputGridFocusHeight = -1;
int outputGridFocusX = -1;
int outputGridFocusY = -1;
char outputGridFocusImage[PATH_MAX];
int outputGridSequence = 0;
char constraintName[PATH_MAX];
int nFixedImages = 0;
char **fixedImages = 0;
int showConstraints = 0;
int epochIterations = 512;
int foldRecovery = 0;

int nImages = 0;
Image *images = 0;
int *imageHashTable = 0;
IntraImageMap **mapHashTable = 0;
int myFirstImage, myLastImage;

int nMaps = 0;
InterImageMap *maps = 0;
int mapsSize = 0;
int mapsPos = 0;

int nPhases;
CommPhase *commPhases;
int bufferSize = 0;
float *buffer = 0;

FILE *logFile = 0;
float kAbsolute = 0.0;
float kIntra = 0.0;
float kInter = 0.0;
float fixedDamping = 0.0;
float threshold = 1.0;  /* ppm change per iteration */
int startLevel = -1;
int endLevel = -1;
int nLevels = -1;

int startFactor;   /* spacing at starting level */
int endFactor;	   /* spacing at end level */
int factor;        /* spacing at current level */

int nSteps = 0;
Step *steps = NULL;

int foldDetected = 0;
int processWithFold = -1;
char foldImage[PATH_MAX];
int foldWidth = 0;
int foldHeight = 0;
int foldX = 0;
int foldY = 0;

#if PDEBUG
int ioi = 0;
float poixv = 844.0;
float poiyv = 14038.0;
int poix, poiy;
#endif

#define DIR_HASH_SIZE	8192
char *dirHash[DIR_HASH_SIZE];

unsigned char colors[16][3] =
  {{100, 0, 0},
   {200, 0, 0},
   {255, 0, 0},
   {255, 50, 0},
   {255, 150, 0},

   {255, 255, 0},
   {160, 255, 0},
   {80, 255, 0},
   {0, 255, 0},
   {0, 255, 128},

   {0, 255, 255},
   {0, 128, 255},
   {0, 0, 255},
   {100, 0, 255},
   {128, 0, 200},

   {100, 0, 100}
  };

void Error (char *fmt, ...);
void Log (char *fmt, ...);
void Output (int level, int iter);
void OutputGrid (int level, int iter);
void hsv_to_rgb (unsigned char* r, unsigned char *g, unsigned char *b,
		 float h, float s, float v);
void DrawLine(unsigned char *img, float px0, float py0, float px1, float py1,
	      int imageWidth, int imageHeight,
	      unsigned char r, unsigned char g, unsigned char b);
void UpdateDeltas (int level, int iter);
void RefinePositions (int prevLevel, int level);
void PlanCommunications (int level);
void CommunicatePositions (int level);
unsigned int Hash (char *s);
unsigned int HashMap (char *s, int nx, int ny);
int CreateDirectories (char *fn);
int CheckForFolds (int level, int iter, int reportError);
void ComputeSpringConstants (int level, int iter);
int Extrapolate (float *prx, float *pry, int ix, int iy, float arrx, float arry,
		 MapElement* map, int mw, int mh, float threshold);

int
main (int argc, char **argv)
{
  char fn[PATH_MAX];
  FILE *f;
  int x, y, z;
  int i, j, k;
  int error;
  int iw, ih;
  int im;
  int n;
  int w, h;
  float rx, ry, rc;
  int irx, iry;
  float rrx, rry;
  float r00, r01, r10, r11;
  float xv, yv;
  Point *pt;
  int dx, dy;
  float *fx;
  float *fy;
  double totalEnergy;
  double prevTotalEnergy;
  double deltaEnergy;
  double energy;
  double epochInitialTotalEnergy;
  int iter;
  float deltaX, deltaY;
  float nomD;
  float d;
  float force;
  float kdx, kdy;
  float dampingFactor;
  MPI_Status status;
  int nDecrease;
  int nIncrease;
  float forceOverD;
  float dfx, dfy;
  unsigned char red, green, blue;
  float angle;
  float mag;
  int nPts;
  int nFloats;
  Point *p00, *p10, *p01, *p11;
  float kIntraThisSection;
  struct stat sb;
  int ixv, iyv;
  int ind;
  float rx00, rx01, rx10, rx11;
  float ry00, ry01, ry10, ry11;
  float rc00, rc01, rc10, rc11;
  time_t lastOutput;
  int cx, cy;
  float consX, consY;
  float maxF;
  float globalMaxF;
  char hostName[256];
  double basis;
  float maxStepX, maxStepY;
  int nx, ny, nz;
  int level;
  char termName[PATH_MAX];
  char triggerName[PATH_MAX];
  unsigned int hv;
  int found;
  int mLevel;
  int mw, mh;
  int mxMin, myMin;
  char imName0[PATH_MAX], imName1[PATH_MAX];
  int imageParamsSize;
  int imageParamsPos;
  float *imageParams;
  int imageNamesSize;
  int imageNamesPos;
  char *imageNames;
  int mapParamsSize;
  int mapParamsPos;
  float *mapParams;
  char imageName0[PATH_MAX], imageName1[PATH_MAX];
  char pairName[PATH_MAX];
  int modelNamesSize;
  int modelNamesPos;
  char *modelNames;
  char line[LINE_LENGTH];
  int nItems;
  char imageName[PATH_MAX];
  int imageNameLen;
  int width, height;
  double rotation, scale;
  double tx, ty;
  double weight;
  char modelName[PATH_MAX];
  int modelNameLen;
  int firstImage, lastImage;
  int mapNamesSize;
  int mapNamesPos;
  char *mapNames = 0;
  int imageNameLen0, imageNameLen1, pairNameLen;
  int totalLength;
  int image0, image1;
  IntraImageMap *iim;
  IntraImageSpring *iis;
  MapElement *modelMap;
  char msg[PATH_MAX+256];
  MapElement *map;
  int mFactor;
  int stripsSize;
  int springsSize;
  int springsPos;
  int prevX, prevY;
  int irrx, irry;
  int firstInStrip;
  InterImageMap *m;
  InterImageStrip *strip;
  InterImageSpring *s;
  int fixedImageNameSize;
  Point *pts0, *pts1;
  ModelPoint *ipt;
  float cost, sint;
  Point *pts;
  int phase;
  int op;
  int y0;
  int nx1, ny1;
  CommPhase *cp;
  int separation;
  int mp;
  int *phaseCount;
  int *globalPhaseCount;
  int outputRequestedIter;
  int terminationRequestedIter;
  int refinementRequestedIter;
  int controlFlag;
  ModelPoint *mpts;
  int ns;
  int *pnStrips;
  InterImageStrip **pStrips;
  int *pnSprings;
  InterImageSpring **pSprings;
  int nStrips;
  InterImageStrip *strips;
  InterImageSpring *springs;
  int mx, my;
  IntraImageSpring **piSprings;
  MapElement *initialMap;
  int nSprings;
  IntraImageSpring *iSprings;
  float msk;
  float sk;
  int sme;
  int prevLevel;
  int step;
  cpu_set_t cpumask;
  int dnx, dny;
  int mbpl;
  unsigned char *mask;
  float *dist;
  /* DECLS */

  /* initialize MPI */
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    Error("Could not do MPI_Init\n");
  if (MPI_Comm_size(MPI_COMM_WORLD, &np) != MPI_SUCCESS)
    Error("Could not do MPI_Comm_size\n");
  if (MPI_Comm_rank(MPI_COMM_WORLD, &p) != MPI_SUCCESS)
    Error("Could not do MPI_Comm_rank\n");
  //  if (MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN) != MPI_SUCCESS)
  //    Error("Could not set MPI_ERRORS_RETURN.\n");

  printf("I am node %d of %d\n", p, np);
  fflush(stdout);
  
  // FIX -- TEMPORARY HACK
  sleep(5);
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

  sprintf(fn, "logs/%0.2d.log", p);
  logFile = fopen(fn, "w");
  gethostname(hostName, 255);
  Log("Process %d is running on %s\n", p, hostName);
  memset(dirHash, 0, DIR_HASH_SIZE * sizeof(char*));
  if (p == 0)
    {
      error = 0;
      imageListName[0] = '\0';
      mapListName[0] = '\0';
      mapsName[0] = '\0';
      initialMapsName[0] = '\0';
      modelsName[0] = '\0';
      outputName[0] = '\0';
      outputGridName[0] = '\0';
      outputGridFocusImage[0] = '\0';
      constraintName[0] = '\0';
      schedule[0] = '\0';

      for (i = 0; i < argc; ++i)
	Log("ARGV[%d] = %s\n", i, argv[i]);
      for (i = 1; i < argc; ++i)
	if (strcmp(argv[i], "-image_list") == 0)
	  {
	    if (++i == argc)
	      {
		error = 1;
		break;
	      }
	    strcpy(imageListName, argv[i]);
	  }
	else if (strcmp(argv[i], "-map_list") == 0)
	  {
	    if (++i == argc)
	      {
		error = 1;
		break;
	      }
	    strcpy(mapListName, argv[i]);
	  }
	else if (strcmp(argv[i], "-maps") == 0)
	  {
	    if (++i == argc)
	      {
		error = 1;
		break;
	      }
	    strcpy(mapsName, argv[i]);
	  }
	else if (strcmp(argv[i], "-models") == 0)
	  {
	    if (++i == argc)
	      {
		error = 1;
		break;
	      }
	    strcpy(modelsName, argv[i]);
	  }
	else if (strcmp(argv[i], "-initial_maps") == 0)
	  {
	    if (++i == argc)
	      {
		error = 1;
		break;
	      }
	    strcpy(initialMapsName, argv[i]);
	  }
	else if (strcmp(argv[i], "-output") == 0)
	  {
	    if (++i == argc)
	      {
		error = 1;
		break;
	      }
	    strcpy(outputName, argv[i]);
	  }	    
	else if (strcmp(argv[i], "-incremental") == 0)
	  outputIncremental = 1;
	else if (strcmp(argv[i], "-output_grid") == 0)
	  {
	    if (++i == argc)
	      {
		error = 1;
		break;
	      }
	    strcpy(outputGridName, argv[i]);
	  }	    
	else if (strcmp(argv[i], "-grid_size") == 0)
	  {
	    if (++i == argc ||
		sscanf(argv[i], "%dx%d", &outputGridWidth, &outputGridHeight) != 2)
	      {
		error = 1;
		break;
	      }
	  }
	else if (strcmp(argv[i], "-overlay") == 0)
	  outputGridOverlay = 1;
	else if (strcmp(argv[i], "-show_springs") == 0)
	  outputGridSprings = 1;
	else if (strcmp(argv[i], "-interval") == 0)
	  {
	    if (++i == argc || sscanf(argv[i], "%d", &outputGridInterval) != 1)
	      {
		error = 1;
		break;
	      }
	  }
	else if (strcmp(argv[i], "-focus") == 0)
	  {
	    if (++i == argc ||
		sscanf(argv[i], "%dx%d+%d+%d,%s",
		       &outputGridFocusWidth, &outputGridFocusHeight,
		       &outputGridFocusX, &outputGridFocusY,
		       outputGridFocusImage) != 5)
	      {
		error = 1;
		break;
	      }
	  }
	else if (strcmp(argv[i], "-fixed") == 0)
	  {
	    if (++i == argc)
	      {
		error = 1;
		break;
	      }
	    ++nFixedImages;
	    fixedImages = realloc(fixedImages, nFixedImages * sizeof(char *));
	    fixedImages[nFixedImages-1] = (char *) malloc(strlen(argv[i]) + 1);
	    strcpy(fixedImages[nFixedImages-1], argv[i]);
	  }
	else if (strcmp(argv[i], "-schedule") == 0)
	  {
	    if (++i == argc ||
		sscanf(argv[i], "%s", schedule) != 1)
	      {
		error = 1;
		break;
	      }
	  }
	else if (strcmp(argv[i], "-constraints") == 0)
	  {
	    if (++i == argc)
	      {
		error = 1;
		break;
	      }
	    strcpy(constraintName, argv[i]);
	  }
	else if (strcmp(argv[i], "-fold_recovery") == 0)
	  foldRecovery = 1;
	else error = 1;

      if (error)
	{
	  fprintf(stderr, "Usage: align12\n");
	  fprintf(stderr, "              -image_list images_file\n");
	  fprintf(stderr, "              -map_list maps_file\n");
	  fprintf(stderr, "              -output output_prefix\n");
	  fprintf(stderr, "              -schedule schedule_file\n");
	  fprintf(stderr, "              [-output_grid output_grid_prefix]\n");
	  fprintf(stderr, "              [-grid_size widthxheight]\n");
	  fprintf(stderr, "              [-overlay]\n");
	  fprintf(stderr, "              [-interval iterations]\n");
	  fprintf(stderr, "              [-focus x,y,image_name]\n");
	  fprintf(stderr, "              [-fixed image_name]\n");
	  fprintf(stderr, "              [-constraints constraints_prefix]\n");
	  fprintf(stderr, "              [-fold_recovery]\n");
	  exit(1);
	}
      
      /* images_file contains one line for each image to be aligned:
            image_name rotation scale tx ty [map_file]
	    image_name rotation scale tx ty [map_file]
            ...
	
         maps_file contains one line for each map between images:
            image_name0 image_name1 map_file weight
            image_name0 image_name1 map_file weight
            ...

	 constraint files are based on the image name with a .con
	 suffix; they contains a list of points:
	    x y constrained_location_x constrained_location_y
	    ...

	 a schedule file contains 1 line per step in the format:
            level kIntra kInter kAbsolute [threshold [dampingFactor]]

      */

      /* check that at least minimal parameters were supplied */
      if (imageListName[0] == '\0')
	Error("-image_list must be specified\n");
      if (mapListName[0] == '\0')
	Error("-map_list must be specified\n");
      if (outputName[0] == '\0' && outputGridName[0] == '\0')
	Error("At least one of -output or -outputgrid must be specified.\n");
      if (schedule[0] == '\0')
	Error("-schedule must be specified\n");
      
      /* read in the schedule */
      f = fopen(schedule, "r");
      if (f == NULL)
	Error("Could not open schedule file %s\n", schedule);
      while (fgets(line, LINE_LENGTH, f) != NULL)
	{
	  kAbsolute = 0.0;
	  threshold = 1.0;
	  dampingFactor = 0.0;
	  nItems = sscanf(line, "%d%f%f%f%f%f",
			  &level, &kIntra, &kInter,
			  &kAbsolute, &threshold, &dampingFactor);
	  if (nItems < 3)
	    Error("Malformed line in %s:\n%s\n", schedule, line);
	  if (nSteps > 0 && level > steps[nSteps-1].level)
	    Error("Levels must monotonically decrease in schedule.\n");
	  steps = (Step *) realloc(steps, (nSteps+1) * sizeof(Step));
	  steps[nSteps].level = level;
	  steps[nSteps].kIntra = kIntra;
	  steps[nSteps].kInter = kInter;
	  steps[nSteps].kAbsolute = kAbsolute;
	  steps[nSteps].threshold = threshold;
	  steps[nSteps].dampingFactor = dampingFactor;
	  ++nSteps;
	}
      fclose(f);
      if (nSteps == 0)
	Error("At least one step must be listed in schedule file %s\n",
	      schedule);
    }

  /* broadcast the info */
  if (MPI_Bcast(mapsName, PATH_MAX, MPI_CHAR, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(initialMapsName, PATH_MAX, MPI_CHAR, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(modelsName, PATH_MAX, MPI_CHAR, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(outputName, PATH_MAX, MPI_CHAR, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&outputIncremental, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(outputGridName, PATH_MAX, MPI_CHAR, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&outputGridWidth, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&outputGridHeight, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&outputGridOverlay, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&outputGridSprings, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&outputGridInterval, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&outputGridFocusWidth, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&outputGridFocusHeight, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&outputGridFocusX, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&outputGridFocusY, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(outputGridFocusImage, PATH_MAX, MPI_CHAR, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(constraintName, PATH_MAX, MPI_CHAR, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&nFixedImages, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&nSteps, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("Broadcast of parameters failed.\n");

  if (p != 0)
    {
      fixedImages = (char **) malloc(nFixedImages * sizeof(char *));
      steps = (Step *) malloc(nSteps * sizeof(Step));
    }
  for (i = 0; i < nFixedImages; ++i)
    {
      if (p == 0)
	fixedImageNameSize = strlen(fixedImages[i]) + 1;
      else
	fixedImageNameSize = 0;
      if (MPI_Bcast(&fixedImageNameSize, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
	Error("Broadcast of fixedImageNameSize failed\n");
      if (p != 0)
	fixedImages[i] = (char *) malloc(fixedImageNameSize);
      if (MPI_Bcast(fixedImages[i], fixedImageNameSize, MPI_CHAR, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
	Error("Broadcast of fixedImages failed\n");
    }      
  if (MPI_Bcast(steps, nSteps * sizeof(Step) / sizeof(float), MPI_FLOAT, 0,
		MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("Broadcast of steps failed\n");
  
  if (p == 0)
    {
      /* read the images file */
      f = fopen(imageListName, "r");
      if (f == NULL)
	Error("Could not open file %s for reading\n", imageListName);
      nImages = 0;
      imageParamsSize = 0;
      imageParamsPos = 0;
      imageParams = 0;
      imageNamesSize = 0;
      imageNamesPos = 0;
      imageNames = 0;
      modelNamesSize = 0;
      modelNamesPos = 0;
      modelNames = 0;
      while (fgets(line, LINE_LENGTH, f) != NULL)
	{
	  nItems = sscanf(line, "%s%d%d%lf%lf%lf%lf%s",
			  imageName, &width, &height,
			  &rotation, &scale, &tx, &ty,
			  modelName);
	  if (nItems != 7 && nItems != 8)
	    Error("Malformed line in %s:\n%s\n", imageListName, line);
	  while (imageParamsPos + 6 >= imageParamsSize)
	    {
	      imageParamsSize = (imageParamsSize > 0) ? imageParamsSize * 2 : 1024;
	      imageParams = (float *) realloc(imageParams, imageParamsSize * sizeof(float));
	    }
	  imageParams[imageParamsPos] = width;
	  imageParams[imageParamsPos+1] = height;
	  imageParams[imageParamsPos+2] = rotation;
	  imageParams[imageParamsPos+3] = scale;
	  imageParams[imageParamsPos+4] = tx;
	  imageParams[imageParamsPos+5] = ty;
	  imageParamsPos += 6;

	  imageNameLen = strlen(imageName);
	  while (imageNamesPos + imageNameLen + 1 >= imageNamesSize)
	    {
	      imageNamesSize = (imageNamesSize > 0) ? imageNamesSize * 2 : 1024;
	      imageNames = (char *) realloc(imageNames, imageNamesSize);
	    }
	  strcpy(&imageNames[imageNamesPos], imageName);
	  imageNamesPos += imageNameLen + 1;

	  if (nItems == 7)
	    modelName[0] = '\0';
	  modelNameLen = strlen(modelName);
	  while (modelNamesPos + modelNameLen + 1 >= modelNamesSize)
	    {
	      modelNamesSize = (modelNamesSize > 0) ? modelNamesSize * 2 : 1024;
	      modelNames = (char *) realloc(modelNames, modelNamesSize);
	    }
	  strcpy(&modelNames[modelNamesPos], modelName);
	  modelNamesPos += modelNameLen + 1;

	  ++nImages;
	}
      fclose(f);

      imageParamsSize = imageParamsPos;
      imageParams = (float *) realloc(imageParams, imageParamsSize * sizeof(float));
      imageNamesSize = imageNamesPos;
      imageNames = (char *) realloc(imageNames, imageNamesSize);
      modelNamesSize = modelNamesPos;
      modelNames = (char *) realloc(modelNames, modelNamesSize);
    }

  /* broadcast the image info */
  if (MPI_Bcast(&nImages, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&imageParamsSize, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&imageNamesSize, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&modelNamesSize, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("Broadcast of nImages, imageParamsSize, imageNamesSize, or modelNamesSize failed.\n");

  if (p != 0)
    {
      imageParams = (float *) malloc(imageParamsSize * sizeof(float));
      imageNames = (char *) malloc(imageNamesSize * sizeof(char));
      modelNames = (char *) malloc(modelNamesSize * sizeof(char));
    }
  if (MPI_Bcast(imageParams, imageParamsSize, MPI_FLOAT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(imageNames, imageNamesSize, MPI_CHAR, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(modelNames, modelNamesSize, MPI_CHAR, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("Broadcast of imageParams, imageNames, and modelNames failed.\n");

  images = (Image *) malloc(nImages * sizeof(Image));
  imageHashTable = (int *) malloc(nImages * sizeof(int));
  for (i = 0; i < nImages; ++i)
    imageHashTable[i] = -1;
  imageNamesPos = 0;
  modelNamesPos = 0;
  startLevel = steps[0].level;
  endLevel = steps[nSteps-1].level;
  nLevels = startLevel - endLevel + 1;
  startFactor = 1 << startLevel;
  endFactor = 1 << endLevel;
  for (i = 0; i < nImages; ++i)
    {
      images[i].name = &imageNames[imageNamesPos];
      imageNamesPos += strlen(&imageNames[imageNamesPos]) + 1;
      images[i].width = imageParams[6*i];
      images[i].height = imageParams[6*i+1];
      images[i].rotation = imageParams[6*i+2];
      images[i].scale = imageParams[6*i+3];
      images[i].tx = imageParams[6*i+4];
      images[i].ty = imageParams[6*i+5];
      images[i].fixed = 0;
      images[i].owner = -1;
      images[i].needed = 0;
      images[i].sendTo = NULL;
      images[i].nx = (images[i].width + startFactor - 1) / startFactor + 1;
      images[i].ny = (images[i].height + startFactor - 1) / startFactor + 1;
      images[i].points = 0;
      images[i].modelName = &modelNames[modelNamesPos];
      modelNamesPos += strlen(&modelNames[modelNamesPos]) + 1;
      images[i].map = 0;
      hv = Hash(images[i].name) % nImages ;
      images[i].next = imageHashTable[hv];
      imageHashTable[hv] = i;
    }
  free(imageParams);

  for (i = 0; i < np; ++i)
    {
      firstImage = (i * nImages) / np;
      lastImage = ((i+1) * nImages / np) - 1;
      for (j = firstImage; j <= lastImage; ++j)
	{
	  images[j].owner = i;
	  images[j].needed = 0;
	  if (i == p)
	    {
	      images[j].sendTo = (unsigned char *)
		malloc((np + 7) >> 3);
	      memset(images[j].sendTo, 0,
		     (np + 7) >> 3);
	    }
	}
    }
  myFirstImage = (p * nImages) / np;
  myLastImage = ((p+1) * nImages / np) - 1;
  for (i = 0; i < nFixedImages; ++i)
    {
      hv = Hash(fixedImages[i]) % nImages ;
      for (j = imageHashTable[hv]; j >= 0; j = images[j].next)
	if (strcmp(images[j].name, fixedImages[i]) == 0)
	  {
	    images[j].fixed = 1;
	    break;
	  }
      if (j < 0)
	Error("Could not find fixed image %s in set of images.\n",
	      fixedImages[i]);
    }

  Log("On node %d first = %d last = %d (nz = %d)\n",
      p, myFirstImage, myLastImage, nImages);

  /* set trigger and termination files */
  sprintf(triggerName, "%strigger", outputName);
  sprintf(termName, "%sterm", outputName);

  if (p == 0)
    {
      /* read the maps file */
      f = fopen(mapListName, "r");
      if (f == NULL)
	Error("Could not open file %s for reading\n", mapListName);

      nMaps = 0;
      mapNamesSize = 0;
      mapNamesPos = 0;
      mapNames = 0;
      mapParamsSize = 0;
      mapParamsPos = 0;
      mapParams = 0;
      while (fgets(line, LINE_LENGTH, f) != NULL)
	{
	  nItems = sscanf(line, "%s%s%s%lf",
			  imageName0, imageName1, pairName, &weight);
	  if (nItems != 3 && nItems != 4)
	    Error("Malformed line in %s:\n%s\n", mapListName, line);
	  if (nItems == 3)
	    weight = 1.0;

	  imageNameLen0 = strlen(imageName0);
	  imageNameLen1 = strlen(imageName1);
	  pairNameLen = strlen(pairName);
	  totalLength = imageNameLen0 + imageNameLen1 + pairNameLen;
	  while (mapNamesPos + totalLength + 3 >= mapNamesSize)
	    {
	      mapNamesSize = (mapNamesSize > 0) ? mapNamesSize * 2 : 1024;
	      mapNames = (char *) realloc(mapNames, mapNamesSize);
	    }
	  strcpy(&mapNames[mapNamesPos], imageName0);
	  mapNamesPos += imageNameLen0 + 1;
	  strcpy(&mapNames[mapNamesPos], imageName1);
	  mapNamesPos += imageNameLen1 + 1;
	  strcpy(&mapNames[mapNamesPos], pairName);
	  mapNamesPos += pairNameLen + 1;

	  while (mapParamsPos + 1 >= mapParamsSize)
	    {
	      mapParamsSize = (mapParamsSize > 0) ? mapParamsSize * 2 : 1024;
	      mapParams = (float *) realloc(mapParams, mapParamsSize * sizeof(float));
	    }
	  mapParams[mapParamsPos++] = weight;

	  ++nMaps;
	}
      fclose(f);
      mapParamsSize = mapParamsPos;
      mapParams = (float *) realloc(mapParams, mapParamsSize * sizeof(float));
      mapNamesSize = mapNamesPos;
      mapNames = (char *) realloc(mapNames, mapNamesSize);
    }

  /* broadcast the map info */
  if (MPI_Bcast(&nMaps, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&mapParamsSize, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&mapNamesSize, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("Broadcast of nMaps, mapParamsSize, or mapNamesSize failed.\n");
  if (p != 0)
    {
      mapParams = (float *) malloc(mapParamsSize * sizeof(float));
      mapNames = (char *) malloc(mapNamesSize * sizeof(char));
    }
  if (MPI_Bcast(mapParams, mapParamsSize, MPI_FLOAT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(mapNames, mapNamesSize, MPI_CHAR, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("Broadcast of mapParams and mapNames failed.\n");

  /* initialize the commPhases array */
  commPhases = (CommPhase*) malloc(2 * np * sizeof(CommPhase));
  for (i = 0; i < 2 * np; ++i)
    {
      cp = &(commPhases[i]);
      separation = i >> 1;
      if (separation == 0)
	cp->otherProcess = -1;
      else if (((p / separation) ^ i) & 1)
	cp->otherProcess = p - separation;
      else
	cp->otherProcess = p + separation;
      if (cp->otherProcess < 0 || cp->otherProcess >= np)
	cp->otherProcess = -1;
      cp->nSendImages = 0;
      cp->sendImages = NULL;
      cp->nSends = 0;
      cp->nImagesToSend = NULL;
      cp->nReceiveImages = 0;
      cp->receiveImages = NULL;
      cp->nReceives = 0;
      cp->floatsToReceive = NULL;
      cp->nImagesToReceive = NULL;
    }

  if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("MPI_Barrier after cpi failed.\n");

  /* go through the mapNames array and pull out the relevant maps */
  mapNamesPos = 0;
  mapParamsPos = 0;
  mapsSize = 0;
  mapsPos = 0;
  for (i = 0; i < nMaps; ++i)
    {
      if (mapNamesPos >= mapNamesSize ||
	  mapParamsPos >= mapParamsSize)
	Error("Internal error: mapParams or MapNames too short.\n");
      strcpy(imageName0, &mapNames[mapNamesPos]);
      imageNameLen0 = strlen(imageName0);
      mapNamesPos += imageNameLen0 + 1;
      strcpy(imageName1, &mapNames[mapNamesPos]);
      imageNameLen1 = strlen(imageName1);
      mapNamesPos += imageNameLen1 + 1;
      strcpy(pairName, &mapNames[mapNamesPos]); 
      pairNameLen = strlen(pairName);
      mapNamesPos += pairNameLen + 1;
      weight = mapParams[mapParamsPos++];

      /* check if this node holds the source of
	 the map */
      found = 0;
      hv = Hash(imageName0) % nImages;
      for (j = imageHashTable[hv]; j >= 0; j = images[j].next)
	if (strcmp(images[j].name, imageName0) == 0)
	  {
	    image0 = j;
	    if (j >= myFirstImage && j <= myLastImage)
	      found = 1;
	    break;
	  }

      /* check if this node holds the destination
	 of the map */
      hv = Hash(imageName1) % nImages;
      for (j = imageHashTable[hv]; j >= 0; j = images[j].next)
	if (strcmp(images[j].name, imageName1) == 0)
	  {
	    image1 = j;
	    if (j >= myFirstImage && j <= myLastImage)
	      found = 1;
	    break;
	  }

      if (found)
	{
	  images[image0].needed = 1;
	  images[image1].needed = 1;

	  /* create map table entry */
	  if (mapsPos >= mapsSize)
	    {
	      mapsSize = (mapsSize > 0) ? mapsSize * 2 : 1024;
	      maps = (InterImageMap *) realloc(maps, mapsSize *
					       sizeof(InterImageMap));
	    }
	  m = &(maps[mapsPos]);
	  m->name = (char *) malloc(pairNameLen+1);
	  strcpy(m->name, pairName);
	  m->image0 = image0;
	  m->image1 = image1;
	  m->energyFactor = (image0 >= myFirstImage &&
			     image0 <= myLastImage) ? 1.0 : 0.0;
	  m->k = 1.0;
	  m->nStrips = (int *) malloc(nLevels * sizeof(int));
	  memset(m->nStrips, 0, nLevels * sizeof(int));
	  m->strips = (InterImageStrip**)
	    malloc(nLevels * sizeof(InterImageStrip*));
	  memset(m->strips, 0, nLevels * sizeof(InterImageStrip*));
	  m->nSprings = (int *) malloc(nLevels * sizeof(int));
	  memset(m->nSprings, 0, nLevels * sizeof(int));
	  m->springs = (InterImageSpring**)
	    malloc(nLevels * sizeof(InterImageSpring*));
	  memset(m->springs, 0, nLevels * sizeof(InterImageSpring*));

	  if (images[image0].owner != images[image1].owner)
	    {
	      if (images[image0].owner == p)
		{
		  op = images[image1].owner;
		  images[image0].sendTo[op >> 3] |= 0x80 >> (op & 7);
		}
	      else
		{
		  op = images[image0].owner;
		  images[image1].sendTo[op >> 3] |= 0x80 >> (op & 7);
		}
	    }
	  ++mapsPos;
	}
    }
  nMaps = mapsPos;
  maps = (InterImageMap *) realloc(maps, nMaps * sizeof(InterImageMap));

  /* allocate the memory for holding the image points at each level */
  for (i = 0; i < nImages; ++i)
    {
      if (images[i].owner != p && !images[i].needed)
	continue;
      images[i].points = (Point *) malloc(nLevels * sizeof(Point *));
      for (level = startLevel; level >= endLevel; --level)
	{
	  factor = 1 << level;
	  nx = (images[i].width + factor - 1) / factor + 1;
	  ny = (images[i].height + factor - 1) / factor + 1;
	  images[i].points[startLevel - level] =
	    (Point *) malloc(nx * ny * sizeof(Point));
	  memset(images[i].points[startLevel - level], 0,
		 nx * ny * sizeof(Point));
	}
    }

  /* construct the list of images to be sent and received
     at each communication phase */
  for (i = 0; i < nImages; ++i)
    if (images[i].owner == p)
      {
	for (op = 0; op < np; ++op)
	  if (images[i].sendTo[op >> 3] & (0x80 >> (op & 7)))
	    {
	      separation = abs(op - p);
	      mp = p;
	      if (op < mp)
		mp = op;
	      phase = 2 * separation + ((mp / separation) & 1);
	      cp = &(commPhases[phase]);
	      ++(cp->nSendImages);
	      cp->sendImages = (int *) realloc(cp->sendImages,
					       cp->nSendImages * sizeof(int));
	      cp->sendImages[cp->nSendImages-1] = i;
	    }
	nFloats = ((images[i].width + endFactor - 1) / endFactor + 1) *
	  ((images[i].height + endFactor - 1) / endFactor + 1) * 2;
	if (nFloats > bufferSize)
	  bufferSize = nFloats;
      }
    else if (images[i].needed)
      {
	op = images[i].owner;
	separation = abs(op - p);
	mp = p;
	if (op < mp)
	  mp = op;
	phase = 2 * separation + ((mp / separation) & 1);
	cp = &(commPhases[phase]);
	++(cp->nReceiveImages);
	cp->receiveImages = (int *) realloc(cp->receiveImages,
					    cp->nReceiveImages * sizeof(int));
	cp->receiveImages[cp->nReceiveImages-1] = i;
	nFloats = ((images[i].width + endFactor - 1) / endFactor + 1) *
	  ((images[i].height + endFactor - 1) / endFactor + 1) * 2;
	if (nFloats > bufferSize)
	  bufferSize = nFloats;
      }
  if (bufferSize < 4096*1024)
    bufferSize = 4096*1024;
  buffer = (float *) malloc(bufferSize * sizeof(float));
  
  /* eliminate unnecessary phases */
  phaseCount = (int *) malloc(2 * np * sizeof(int));
  memset(phaseCount, 0, 2 * np * sizeof(int));
  globalPhaseCount = (int *) malloc(2 * np * sizeof(int));
  for (phase = 0; phase < 2 * np; ++phase)
    {
      phaseCount[phase] += commPhases[phase].nSendImages;
      phaseCount[phase] += commPhases[phase].nReceiveImages;
    }
  if (MPI_Allreduce(phaseCount, globalPhaseCount, 2 * np,
		    MPI_INT, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("MPI_Allreduce failed.\n");
  j = 0;
  for (i = 0; i < 2 * np; ++i)
    if (globalPhaseCount[i] > 0 && i != j)
      {
	memcpy(&commPhases[j], &commPhases[i], sizeof(CommPhase));
	++j;
      }
  nPhases = j;
  commPhases = (CommPhase *) realloc(commPhases, nPhases * sizeof(CommPhase));
  free(phaseCount);
  free(globalPhaseCount);

  /* initialize the models that are needed */
  mapHashTable = (IntraImageMap **) malloc(nImages * sizeof(IntraImageMap*));
  for (i = 0; i < nImages; ++i)
    mapHashTable[i] = NULL;
  modelMap = NULL;
  for (i = 0; i < nImages; ++i)
    {
      mx = (images[i].width + endFactor - 1) / endFactor + 1;
      my = (images[i].height + endFactor - 1) / endFactor + 1;
      hv = HashMap(images[i].modelName, mx, my) % nImages;
      for (iim = mapHashTable[hv]; iim != NULL; iim = iim->next)
	if (strcmp(iim->name, images[i].modelName) == 0 &&
	    iim->nx == mx && iim->ny == my)
	  break;
      if (iim == NULL)
	{
	  /* create the initial map */
	  iim = (IntraImageMap*) malloc(sizeof(IntraImageMap));
	  iim->name = (char *) malloc(strlen(images[i].modelName) + 1);
	  strcpy(iim->name, images[i].modelName);
	  iim->next = mapHashTable[hv];
	  mapHashTable[hv] = iim;
	  if (iim->name[0] != '\0')
	    {
	      /* read the map from a file */
	      sprintf(fn, "%s%s.map", modelsName, iim->name);
	      if (!ReadMap(fn, &modelMap, &mLevel,
			   &mw, &mh, &mxMin, &myMin,
			   imName0, imName1,
			   msg))
		Error("Could not read map %s:\n  error: %s\n",
		      fn, msg);
	      if (mxMin != 0 || myMin != 0)
		Error("Model map %s is only a partial map.\n",
		      fn);
	    }
	  iim->nSprings = (int *) malloc(nLevels * sizeof(int));
	  memset(iim->nSprings, 0, nLevels * sizeof(int));
	  iim->springs = (IntraImageSpring**)
	    malloc(nLevels * sizeof(IntraImageSpring*));
	  memset(iim->springs, 0, nLevels * sizeof(IntraImageSpring*));
	  for (level = startLevel; level >= endLevel; --level)
	    {
	      factor = 1 << level;
	      nx = (images[i].width + factor - 1) / factor + 1;
	      ny = (images[i].height + factor - 1) / factor + 1;
	      mpts = (ModelPoint*) malloc(nx * ny * sizeof(ModelPoint));
	      pnSprings = &(iim->nSprings[startLevel - level]);
	      piSprings = &(iim->springs[startLevel - level]);
	      if (iim->name[0] == '\0')
		/* use the identity map */
		for (y = 0; y < ny; ++y)
		  for (x = 0; x < nx; ++x)
		    {
		      mpts[y*nx+x].x = x;
		      mpts[y*nx+x].y = y;
		    }
	      else
		{
		  scale = ((double) factor) / (1 << mLevel);
#if PDEBUG
		  poix = (int) floor(poixv / factor + 0.5);
		  poiy = (int) floor(poiyv / factor + 0.5);
#endif
		  for (y = 0; y < ny; ++y)
		    for (x = 0; x < nx; ++x)
		      {
			xv = scale * x;
			yv = scale * y;
			ixv = (int) floor(xv);
			iyv = (int) floor(yv);
			rrx = xv - ixv;
			rry = yv - iyv;
			while (ixv >= mw-1)
			  {
			    --ixv;
			    rrx += 1.0;
			  }
			while (iyv >= mh-1)
			  {
			    --iyv;
			    rry += 1.0;
			  }
			rx00 = modelMap[iyv*mw+ixv].x;
			ry00 = modelMap[iyv*mw+ixv].y;
			rx01 = modelMap[(iyv+1)*mw+ixv].x;
			ry01 = modelMap[(iyv+1)*mw+ixv].y;
			rx10 = modelMap[iyv*mw+ixv+1].x;
			ry10 = modelMap[iyv*mw+ixv+1].y;
			rx11 = modelMap[(iyv+1)*mw+ixv+1].x;
			ry11 = modelMap[(iyv+1)*mw+ixv+1].y;
			rx = rx00 * (rrx - 1.0) * (rry - 1.0)
			  - rx10 * rrx * (rry - 1.0) 
			  - rx01 * (rrx - 1.0) * rry
			  + rx11 * rrx * rry;
			ry = ry00 * (rrx - 1.0) * (rry - 1.0)
			  - ry10 * rrx * (rry - 1.0) 
			  - ry01 * (rrx - 1.0) * rry
			  + ry11 * rrx * rry;
			mpts[y*nx+x].x = rx / scale;
			mpts[y*nx+x].y = ry / scale;
#if DEBUG
#if PDEBUG
			if (i == ioi && x == poix && y == poiy)
#endif
			  Log("iim pts (%lu) %d %d = %f %f\n",
			      mpts, x, y,
			      mpts[y*nx+x].x,
			      mpts[y*nx+x].y);
#endif
		      }
		}
	      // creates the springs at this level
	      *pnSprings = (ny - 1) * (nx - 1) * 4 +
		(ny - 1) + (nx - 1);
	      *piSprings = (IntraImageSpring*) malloc(*pnSprings *
						      sizeof(IntraImageSpring));
	      iis = *piSprings;
	      for (y = 0; y < ny; ++y)
		for (x = 0; x < nx; ++x)
		  {
		    if (x < nx-1 && y > 0)
		      {
			iis->index0 = y*nx+x;
			iis->index1 = (y-1)*nx+x+1;
			++iis;
		      }
		    if (x < nx-1)
		      {
			iis->index0 = y*nx+x;
			iis->index1 = y*nx+x+1;
			++iis;
		      }
		    if (x < nx-1 && y < ny-1)
		      {
			iis->index0 = y*nx+x;
			iis->index1 = (y+1)*nx+x+1;
			++iis;
		      }
		    if (y < ny-1)
		      {
			iis->index0 = y*nx+x;
			iis->index1 = (y+1)*nx+x;
			++iis;
		      }
		  }

	      iis = *piSprings;
	      for (j = 0; j < *pnSprings; ++j)
		{
		  iis->nomD = hypot(mpts[iis->index0].x -
				    mpts[iis->index1].x,
				    mpts[iis->index0].y -
				    mpts[iis->index1].y);
		  iis->k = 1.0;
		  ++iis;
		}
	      if (level == startLevel)
		iim->points = mpts;
	      else
		free(mpts);
	    }
	}
      images[i].map = iim;
    }
  if (modelMap != NULL)
    {
      free(modelMap);
      modelMap = NULL;
    }

  /* initialize the image points */
  initialMap = NULL;
  for (i = 0; i < nImages; ++i)
    {
      if (images[i].owner != p && !images[i].needed)
	continue;
      nx = images[i].nx;
      ny = images[i].ny;
      images[i].points[0] = (Point*) malloc(nx * ny * sizeof(Point));
      pt = images[i].points[0];
      if (images[i].owner != p)
	{
	  memset(pt, 0, nx * ny * sizeof(Point));
	  continue;
	}
      if (initialMapsName[0] != '\0')
	{
	  sprintf(fn, "%s%s.map", initialMapsName, images[i].name);
	  if (!ReadMap(fn, &initialMap, &mLevel,
		       &mw, &mh, &mxMin, &myMin,
		       imName0, imName1,
		       msg))
	    Error("Could not read map %s:\n  error: %s\n",
	      fn, msg);
	  if (mxMin != 0 || myMin != 0)
	    Error("Initial map %s is only a partial map.\n",
		  fn);
	  scale = ((double) startFactor) / (1 << mLevel);
	  for (y = 0; y < ny; ++y)
	    for (x = 0; x < nx; ++x)
	      {
		xv = scale * x;
		yv = scale * y;
		ixv = (int) floor(xv);
		iyv = (int) floor(yv);
		rrx = xv - ixv;
		rry = yv - iyv;
		while (ixv >= mw-1)
		  {
		    --ixv;
		    rrx += 1.0;
		  }
		while (iyv >= mh-1)
		  {
		    --iyv;
		    rry += 1.0;
		  }
		rx00 = initialMap[iyv*mw+ixv].x;
		ry00 = initialMap[iyv*mw+ixv].y;
		rx01 = initialMap[(iyv+1)*mw+ixv].x;
		ry01 = initialMap[(iyv+1)*mw+ixv].y;
		rx10 = initialMap[iyv*mw+ixv+1].x;
		ry10 = initialMap[iyv*mw+ixv+1].y;
		rx11 = initialMap[(iyv+1)*mw+ixv+1].x;
		ry11 = initialMap[(iyv+1)*mw+ixv+1].y;
		rx = rx00 * (rrx - 1.0) * (rry - 1.0)
		  - rx10 * rrx * (rry - 1.0) 
		  - rx01 * (rrx - 1.0) * rry
		  + rx11 * rrx * rry;
		ry = ry00 * (rrx - 1.0) * (rry - 1.0)
		  - ry10 * rrx * (rry - 1.0) 
		  - ry01 * (rrx - 1.0) * rry
		  + ry11 * rrx * rry;
		pt->x = rx / scale;
		pt->y = ry / scale;
		pt->fx = 0.0;
		pt->fy = 0.0;
		++pt;
	      }
	}
      else
	{
	  iim = images[i].map;
	  mpts = iim->points;
	  cost = cos(images[i].rotation) * images[i].scale;
	  sint = sin(images[i].rotation) * images[i].scale;
#if PDEBUG
	  poix = (int) floor(poixv / startFactor + 0.5);
	  poiy = (int) floor(poiyv / startFactor + 0.5);
#endif
	  for (y = 0; y < ny; ++y)
	    for (x = 0; x < nx; ++x)
	      {
		pt->x =  cost * mpts[y*nx+x].x + sint * mpts[y*nx+x].y
		  + images[i].tx / startFactor;
		pt->y =  -sint * mpts[y*nx+x].x + cost * mpts[y*nx+x].y
		  + images[i].ty / startFactor;
		pt->fx = 0.0;
		pt->fy = 0.0;
#if DEBUG
#if PDEBUG
		if (i == ioi && x == poix && y == poiy)
#endif
		  Log("Initializing position of %d(%d,%d) to (%f,%f)\n",
		      i, x, y, pt->x, pt->y);
#endif
		++pt;
	      }
	}
    }
  if (initialMap != NULL)
    {
      free(initialMap);
      initialMap = NULL;
    }

  /* initialize the maps that are needed */
  map = NULL;
  for (i = 0; i < nMaps; ++i)
    {
      sprintf(fn, "%s%s.map", mapsName, maps[i].name);
      if (!ReadMap(fn, &map, &mLevel,
		   &mw, &mh, &mxMin, &myMin,
		   imName0, imName1,
		   msg))
	Error("Could not read map %s:\n  error: %s\n",
	      fn, msg);
      mFactor = 1 << mLevel;
      m = &maps[i];

      /* create an array with coverage corresponding to the
	 startLevel, and resolution determined by the mLevel
	 fill in the entries with a 1 if the map has nonzero c
	 fill in the other entries with 0
	 compute distance
	 let the threshold at each level be one-half of the
	   distance from grid point to grid point in terms
	   of the map spacing
	 at each level determine if the distance is less than
	   or equal to the threshold; if so, make a spring
      */
      factor = 1 << startLevel;
      dnx = (images[m->image0].width + factor - 1) / factor + 1;
      dny = (images[m->image0].height + factor - 1) / factor + 1;
      dnx = dnx * factor / mFactor;
      dny = dny * factor / mFactor;
      mbpl = (dnx + 7) >> 3;
      mask = (unsigned char *) malloc(dny * mbpl);
      memset(mask, 0, dny * mbpl);
      for (y = 0; y < mh; ++y)
	for (x = 0; x < mw; ++x)
	  if (map[y*mw+x].c > 0.0)
	    {
	      ixv = x + mxMin;
	      iyv = y + myMin;
	      if (ixv >= 0 && ixv < dnx &&
		  iyv >= 0 && iyv < dny)
		mask[iyv*mbpl+(ixv >> 3)] |= 0x80 >> (ixv & 7);
	    }
      dist = (float *) malloc(dny * dnx * sizeof(float));
      computeDistance(EUCLIDEAN_DISTANCE, dnx, dny, mask, dist);

      for (level = startLevel; level >= endLevel; --level)
	{
	  factor = 1 << level;
	  threshold = ((float) factor) / M_SQRT2 / (1 << mLevel);

	  pnStrips = &(m->nStrips[startLevel - level]);
	  pStrips = &(m->strips[startLevel - level]);
	  pnSprings = &(m->nSprings[startLevel - level]);
	  pSprings = &(m->springs[startLevel - level]);

	  *pnStrips = 0;
	  *pStrips = NULL;
	  *pnSprings = 0;
	  *pSprings = NULL;

	  nx = (images[m->image0].width + factor - 1) / factor + 1;
	  ny = (images[m->image0].height + factor - 1) / factor + 1;
	  nx1 = (images[m->image1].width + factor - 1) / factor + 1;
	  ny1 = (images[m->image1].height + factor - 1) / factor + 1;
#if PDEBUG
	  poix = (int) floor(poixv / factor + 0.5);
	  poiy = (int) floor(poiyv / factor + 0.5);
	  Log("at level %d poix = %d poiy = %d  nx=%d ny=%d nx1=%d ny1=%d\n",
	      level, poix, poiy, nx, ny, nx1, ny1);
#endif
	  stripsSize = 0;
	  springsSize = 0;
	  for (y = 0; y < ny; ++y)
	    {
	      firstInStrip = 1;
	      yv = ((float) (y * factor)) / mFactor;
	      for (x = 0; x < nx; ++x)
		{
		  xv = ((float) (x * factor)) / mFactor;
		  ixv = (int) floor(xv);
		  iyv = (int) floor(yv);
#if DEBUG
#if PDEBUG
		  if (m->image0 == ioi && x == poix && y == poiy)
#endif
		    Log("cons spring at level %d x = %d y = %d xv = %f yv = %f ixv = %d iyv = %d dnx = %d dny = %d dist = %f thresh = %f\n",
			level, x, y, xv, yv, ixv, iyv,
			dnx, dny,
			ixv >= 0 && ixv < dnx && iyv >= 0 && iyv < dny ?
			dist[iyv*dnx+ixv] : -999.0,
			threshold);
#endif
		  if (ixv < 0 || ixv >= dnx ||
		      iyv < 0 || iyv >= dny ||
		      dist[iyv*dnx+ixv] > threshold)
		    {
		      firstInStrip = 1;
		      continue;
		    }
		  rrx = xv - ixv;
		  rry = yv - iyv;
		  ixv -= mxMin;
		  iyv -= myMin;
		  if (ixv >= 0 && ixv < mw-1 &&
		      iyv >= 0 && iyv < mh-1 &&
		      map[iyv*mw+ixv].c > 0.0 &&
		      map[(iyv+1)*mw+ixv].c > 0.0 &&
		      map[iyv*mw+ixv+1].c > 0.0 &&
		      map[(iyv+1)*mw+ixv+1].c > 0.0)
		    {
		      rx00 = map[iyv*mw+ixv].x;
		      ry00 = map[iyv*mw+ixv].y;
		      rx01 = map[(iyv+1)*mw+ixv].x;
		      ry01 = map[(iyv+1)*mw+ixv].y;
		      rx10 = map[iyv*mw+ixv+1].x;
		      ry10 = map[iyv*mw+ixv+1].y;
		      rx11 = map[(iyv+1)*mw+ixv+1].x;
		      ry11 = map[(iyv+1)*mw+ixv+1].y;
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
		      if (!Extrapolate(&rx, &ry, ixv, iyv, rrx, rry,
				       map, mw, mh, threshold))
			{
			  firstInStrip = 1;
			  continue;
			}
		    }
		  rx = rx * mFactor / factor;
		  ry = ry * mFactor / factor;
		  irx = floor(rx + 0.5 + 0.5 / 65536.0);
		  iry = floor(ry + 0.5 + 0.5 / 65536.0);
#if DEBUG
#if PDEBUG
		  if (m->image0 == ioi && x == poix && y == poiy ||
		      m->image1 == ioi && irx == poix && iry == poiy)
#endif
		    Log("spring at lvl = %d from im %d x = %d y = %d xv = %f yv = %f ixv = %d iyv = %d rrx = %f rry = %f to im %d rx = %f ry = %f irx = %d iry = %d\n",
			level, m->image0, x, y, xv, yv, ixv, iyv, rrx, rry,
			m->image1, rx, ry, irx, iry);
#endif

		  if (rx < -0.5 || rx >= nx1 - 0.5 - 0.5 / 65536.0 ||
		      ry < -0.5 || ry >= ny1 - 0.5 - 0.5 / 65536.0)
		    {
		      firstInStrip = 1;
		      continue;
		    }
		  rrx = rx - irx;
		  rry = ry - iry;
		  
		  if (firstInStrip)
		    {
		      firstInStrip = 0;
		      if (++*pnStrips > stripsSize)
			{
			  stripsSize = (stripsSize == 0) ? ny : 2 * stripsSize;
			  *pStrips = (InterImageStrip*)
			    realloc(*pStrips, stripsSize*sizeof(InterImageStrip));
			}
		      strip = &((*pStrips)[*pnStrips - 1]);
		      strip->nSprings = 1;
		      strip->x0 = x;
		      strip->y0 = y;
		      strip->x1 = irx;
		      strip->y1 = iry;
		      prevX = irx;
		      prevY = iry;
		    }
		  else
		    ++(strip->nSprings);
		  if (++*pnSprings > springsSize)
		    {
		      springsSize = (springsSize == 0) ? nx*ny : 2 * springsSize;
		      *pSprings = (InterImageSpring*)
			realloc(*pSprings, springsSize*sizeof(InterImageSpring));
		    }		  
		  s = &((*pSprings)[*pnSprings-1]);
		  s->dx = 0;
		  s->dy = 0;
		  s->k = 255;
		  dx = irx - prevX;
		  dy = iry - prevY;
		  if (dx < -8 || dx > 7 ||
		      dy < -8 || dy > 7)
		    Error("Map is too divergent for internal representation: dx = %d dy = %d\n",
			  dx, dy);
		  s->dxy1 = ((dx + 8) << 4) | (dy + 8);
		  irrx = floor(65536.0 * rrx + 0.5);
		  irry = floor(65536.0 * rry + 0.5);
		  if (irrx < -32768 || irrx > 32767)
		    Error("Internal error: irrx is out-of-range: %d\n", irrx);
		  if (irry < -32768 || irry > 32767)
		    Error("Internal error: irry is out-of-range: %d\n", irry);
		  s->irrx = irrx;
		  s->irry = irry;
		  prevX = irx;
		  prevY = iry;
		}
	    }
	  if (*pnStrips != 0)
	    *pStrips = (InterImageStrip*)
	      realloc(*pStrips, *pnStrips * sizeof(InterImageStrip));
	  if (*pnSprings != 0)
	    *pSprings = (InterImageSpring*)
	      realloc(*pSprings, *pnSprings * sizeof(InterImageSpring));
	  else
	    Error("No springs were generated from map %s\n", m->name);
	}
      free(dist);
      free(mask);
    }
  free(map);
  map = NULL;

#if 0
  /* check what constraints are present */
  if (constraintName[0] != '\0')
    {
      for (i = myFirstImage; i <= myLastImage; ++i)
	{
	  sprintf(fn, "%s%s.con", constraintName, images[i].name);
	  f = fopen(fn, "r");
	  if (f != NULL)
	    {
	      pts = images[i].points[0];
	      j = 0;
	      while (fscanf(f, "%d %d %f %f", &cx, &cy, &consX, &consY) == 4)
		{
		  if (cx < 0 || cx >= nx || cy < 0 || cy >= ny)
		    continue;
		  pt[cy*nx + cx].x = consX;
		  pt[cy*nx + cx].fx = (j + 1) * CONSTRAINED;
		  pt[cy*nx + cx].y = consY;
		  pt[cy*nx + cx].fy = (j + 1) * CONSTRAINED;
		  ++j;
		}
	      fclose(f);
	    }
	}
    }
#endif

  prevLevel = -1;
  for (step = 0; step < nSteps; ++step)
    {
    restartStep:
      level = steps[step].level;
      kIntra = steps[step].kIntra;
      kInter = steps[step].kInter;
      kAbsolute = steps[step].kAbsolute;
      threshold = steps[step].threshold;
      if (steps[step].dampingFactor > 0.0)
	dampingFactor = steps[step].dampingFactor;
      else
	dampingFactor = 0.001;

      if (level != prevLevel)
	{
	  if (prevLevel > 0)
	    RefinePositions(prevLevel, level);
	  PlanCommunications(level);
	  prevLevel = level;
	}

      factor = 1 << level;
      maxStepX = 0.1 * factor;
      maxStepY = 0.1 * factor;
      outputRequestedIter = -1;
      terminationRequestedIter = -1;
      refinementRequestedIter = -1;

      if (p == 0)
	Log("Starting iterations for step %d (level %d)\n", step, level);
      epochInitialTotalEnergy = 1.0e+30;
      prevTotalEnergy = 1.0e+30;
      lastOutput = 0;
      nDecrease = 0;
      nIncrease = 0;
      foldDetected = 0;
      for (iter = 0; ; ++iter)
	{
#if DEBUG
	  Log("Starting iteration %d\n", iter);
#endif
	  CommunicatePositions(level);

	  /* output the grids if requested */
	  if (outputGridName[0] != '\0' &&
	      outputGridInterval > 0 &&
	      iter % outputGridInterval == 0)
	    OutputGrid(level, iter);

	  /* check for folds periodically */
	  if (iter % epochIterations == 0)
	    {
	      foldDetected = CheckForFolds(level, iter,
					   !outputIncremental && !foldRecovery);
	      if (foldDetected)
		processWithFold = p;
	      else
		processWithFold = np;
	      if (MPI_Allreduce(MPI_IN_PLACE, &processWithFold, 1, MPI_INT,
				MPI_MIN, MPI_COMM_WORLD) != MPI_SUCCESS)
		Error("MPI_Allreduce: could not check for folds.\n");
	      if (processWithFold < np)
		{
		  foldDetected = 1;

		  if (foldRecovery)
		    {
		      RecoverFromFold();
		      goto

		      // determine which InterImageMap is most distorted
		      //   in the region of the fold, and set the spring
		      //   constants to 0 in the vicinity of the fold
		      
		    }

		  break;
		}
	    }
		
	  /* update all deltas periodically */
	  if (iter % epochIterations == 0)
	    UpdateDeltas(level, iter);

	  /* update all forces, also computing energy */
	  energy = 0.0;
	  basis = energy;
	  for (i = myFirstImage; i <= myLastImage; ++i)
	    {
	      nx = images[i].nx;
	      ny = images[i].ny;
#if PDEBUG
	      poix = (int) floor(poixv / factor + 0.5);
	      poiy = (int) floor(poiyv / factor + 0.5);
#endif
	      if (kAbsolute > 0.0)
		{
		  /* compute absolute location forces */
		  pt = images[i].points[startLevel - level];
		  for (y = 0; y < ny; ++y)
		    for (x = 0; x < nx; ++x)
		      {
			deltaX = ((float) x) - pt->x;
			deltaY = ((float) y) - pt->y;
			if (pt->fx < 0.5 * CONSTRAINED)
			  pt->fx = kAbsolute * deltaX;
			if (pt->fy < 0.5 * CONSTRAINED)
			  pt->fy = kAbsolute * deltaY;
			++pt;
			energy += kAbsolute * (deltaX * deltaX + deltaY * deltaY);
		      }
		}
	      else
		{
		  /* zero all forces */
		  pt = images[i].points[startLevel - level];
		  for (y = 0; y < ny; ++y)
		    for (x = 0; x < nx; ++x)
		      {
			if (pt->fx < 0.5 * CONSTRAINED)
			  pt->fx = 0.0;
			if (pt->fy < 0.5 * CONSTRAINED)
			  pt->fy = 0.0;
			++pt;
		      }
		}

	      /* add in intra-section forces */
	      iim = images[i].map;
	      pts = images[i].points[startLevel - level];
	      nSprings = iim->nSprings[startLevel - level];
	      iSprings = iim->springs[startLevel - level];
	      for (k = 0; k < nSprings; ++k)
		{
		  iis = &(iSprings[k]);
		  deltaX = pts[iis->index1].x - pts[iis->index0].x;
		  deltaY = pts[iis->index1].y - pts[iis->index0].y;
		  d = sqrt(deltaX * deltaX + deltaY * deltaY);
		  sk = kIntra * iis->k;
		  force = sk * (d - iis->nomD);
		  energy += force * (d - iis->nomD);
		  if (d != 0.0)
		    {
		      forceOverD = force / d;
		      dfx = forceOverD * deltaX;
		      dfy = forceOverD * deltaY;
		      pts[iis->index0].fx += dfx;
		      pts[iis->index0].fy += dfy;
		      pts[iis->index1].fx -= dfx;
		      pts[iis->index1].fy -= dfy;
		    }
#if DEBUG
#if PDEBUG
		  if (i == ioi && iis->index0 % images[i].nx == poix &&
		      iis->index0 / images[i].nx == poiy)
#endif
		    Log("intraforce: %d(%d,%d) - %d(%d,%d): (%f %f) to (%f %f) dist %f nom %f force (%f %f)\n",
		      i, iis->index0 % images[i].nx, iis->index0 / images[i].nx,
		      i, iis->index1 % images[i].nx, iis->index1 / images[i].nx,
		      pts[iis->index0].x, pts[iis->index0].y,
		      pts[iis->index1].x, pts[iis->index1].y,
		      d, iis->nomD,
		      d != 0.0 ? dfx : 1000000000.0,
		      d != 0.0 ? dfy : 1000000000.0);
#endif
		}
	    }
	  if (p == 0 && iter % 100 == 0)
	    Log("intra-energy = %f  (kIntra = %f)\n", energy - basis, kIntra);
	  //	  pts0 = images[0].points[startLevel-level];
	  //	  nx = images[0].nx;
	  //	  if (iter == 0)
	  //	    Log("PTEST %f %f\n", pts[2*nx].x, pts[2*nx].y);

	  /* add in inter-image forces (one spring method) */
	  basis = energy;
	  for (i = 0; i < nMaps; ++i)
	    {
	      m = &maps[i];
	      msk = kInter * m->k;
	      sme = m->energyFactor != 0.0;
	      nStrips = m->nStrips[startLevel - level];
	      strips = m->strips[startLevel - level];
	      springs = m->springs[startLevel - level];
	      pts0 = images[m->image0].points[startLevel - level];
	      pts1 = images[m->image1].points[startLevel - level];
	      nx = images[m->image0].nx;
	      nx1 = images[m->image1].nx;
	      springsPos = 0;
	      for (j = 0; j < nStrips; ++j)
		{
		  strip = &(strips[j]);
		  ns = strip->nSprings;
		  x = strip->x0;
		  y = strip->y0;
		  irx = strip->x1;
		  iry = strip->y1;
		  for (k = 0; k < ns; ++k, ++x)
		    {
		      s = &(springs[springsPos++]);
		      irx += ((s->dxy1) >> 4) - 8;
		      iry += ((s->dxy1) & 0xf) - 8;
		      xv = pts1[iry*nx1+irx].x + 0.0001 * s->dx;
		      yv = pts1[iry*nx1+irx].y + 0.0001 * s->dy;
		      deltaX = xv - pts0[y*nx + x].x;
		      deltaY = yv - pts0[y*nx + x].y;
		      sk = (s->k / 255.0) * msk;
		      kdx = sk * deltaX;
		      kdy = sk * deltaY;
#if DEBUG
#if PDEBUG
		      if (m->image0 == ioi && x == poix && y == poiy ||
			  m->image1 == ioi && irx == poix && iry == poiy)
#endif
		      Log("interforce %d(%d,%d) - %d(%d,%d): (%f %f) to (%f %f) (%d %d) (%f %f) (%f %f) (%f %f)%\n",
			  m->image0, x, y,
			  m->image1, irx, iry,
			  pts0[y*nx+x].x, pts0[y*nx+x].y,
			  xv, yv,
			  irx, iry,
			  pts1[iry*nx1+irx].x, pts1[iry*nx1+irx].y,
			  0.0001 * s->dx, 0.0001 * s->dy,
			  kdx, kdy);
#endif
		      pts0[y*nx + x].fx += kdx;
		      pts0[y*nx + x].fy += kdy;
		      if (sme)
			energy += sk * (deltaX * deltaX + deltaY * deltaY);
		      pts1[iry*nx1 + irx].fx -= kdx;
		      pts1[iry*nx1 + irx].fy -= kdy;
		    }
		}
	    }
	  if (p == 0 && iter % 100 == 0)
	    Log("inter-energy = %f  (kInter = %f)\n", energy - basis, kInter);

	  maxF = 0.0;
	  for (i = myFirstImage; i <= myLastImage; ++i)
	    {
	      if (images[i].fixed)
		continue;
	      nPts = images[i].nx * images[i].ny;
	      pt = images[i].points[startLevel-level];
	      for (k = 0; k < nPts; ++k, ++pt)
		{
		  if (pt->fx < 0.5 * CONSTRAINED)
		    if (pt->fy < 0.5 * CONSTRAINED)
		      force = dampingFactor * hypot(pt->fx, pt->fy);
		    else
		      force = dampingFactor * fabsf(pt->fx);
		  else 
		    if (pt->fy < 0.5 * CONSTRAINED)
		      force = dampingFactor * fabsf(pt->fy);
		    else
		      continue;
		  if (force > maxF)
		    maxF = force;
		}
	    }

	  /* find global maximum */
	  if (MPI_Allreduce(&maxF, &globalMaxF, 1, MPI_FLOAT, MPI_MAX,
			    MPI_COMM_WORLD) != MPI_SUCCESS)
	    Error("Could not find global maximum force\n");

	  /* update all positions */
	  if (globalMaxF > 0.5)
	    scale = dampingFactor * 0.5 / globalMaxF;
	  else
	    scale = dampingFactor;
	  //	  printf("df = %f scale = %f msx = %f msy = %f\n",
	  //		 dampingFactor, scale, maxStepX, maxStepY);

	  for (i = myFirstImage; i <= myLastImage; ++i)
	    {
	      if (images[i].fixed)
		continue;
	      nPts = images[i].nx * images[i].ny;
	      pt = images[i].points[startLevel-level];
	      for (k = 0; k < nPts; ++k)
		{
		  if (pt->fx < 0.5 * CONSTRAINED)
		    {
		      deltaX = scale * pt->fx;
		      if (deltaX > maxStepX)
			deltaX = maxStepX;
		      pt->x += deltaX;

		    }
		  if (pt->fy < 0.5 * CONSTRAINED)
		    {
		      deltaY = scale * pt->fy;
		      if (deltaY > maxStepY)
			deltaY = maxStepY;
		      pt->y += deltaY;
		    }
#if DEBUG
#if PDEBUG
		  if (i == ioi && k % images[i].nx == poix &&
		      k / images[i].nx == poiy)
#endif
		  Log("moving point %d(%d,%d) at (%f %f) by (%f %f) to (%f %f)\n",
		      i, k % images[i].nx, k / images[i].nx,
		      pt->x - deltaX, pt->y - deltaY,
		      deltaX, deltaY,
		      pt->x, pt->y);
#endif
		  ++pt;
		}
	    }

	  if (MPI_Allreduce(&energy, &totalEnergy, 1, MPI_DOUBLE,
			    MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS)
	    Error("Could not sum up energies.\n");
	  /*	  if (p == 0 && (iter % 1000) == 999) */
	  if (p == 0 && iter % 100 == 0)
	    Log("After %d iterations, total energy is %f  (df = %f mgf = %f)\n",
		iter, totalEnergy, dampingFactor, globalMaxF);
	  deltaEnergy = totalEnergy - prevTotalEnergy;
	  prevTotalEnergy = totalEnergy;
	  if (iter % epochIterations == 0)
	    epochInitialTotalEnergy = totalEnergy;

	  // check periodically for termination conditions
	  if (iter % epochIterations == epochIterations-1)
	    {
	      if (p == 0)
		{
		  controlFlag = 0;
		  if (triggerName[0] != '\0' &&
		      stat(triggerName, &sb) == 0 &&
		      sb.st_mtime > lastOutput)
		    {
		      Log("Update of file %s forced write of output images.\n",
			  triggerName);
		      controlFlag |= 4;
		      lastOutput = sb.st_mtime;
		    }
		  if (termName[0] != '\0' &&
		      (f = fopen(termName, "r")) != NULL)
		    {
		      fclose(f);
		      Log("Presence of file %s forcing termination.\n", termName);
		      controlFlag |= 2;
		    }
		  if (epochInitialTotalEnergy - totalEnergy <
		      totalEnergy * threshold * 0.000001 * epochIterations ||
		      totalEnergy < 0.000001)
		    controlFlag |= 1;
		}
	      if (MPI_Bcast(&controlFlag, 1, MPI_INT, 0, MPI_COMM_WORLD) !=
		  MPI_SUCCESS)
		Error("Could not broadcast control flag.\n");
	      if (controlFlag & 4)
		outputRequestedIter = iter;
	      if (controlFlag & 2)
		terminationRequestedIter = iter;
	      if (controlFlag & 1)
		refinementRequestedIter = iter;
	    }

	  if (deltaEnergy > 0.0)
	    {
	      nDecrease = 0;
	      ++nIncrease;
	      if (nIncrease > 1000 ||
		  nIncrease > 10 && deltaEnergy > totalEnergy ||
		  deltaEnergy > 1000.0 * totalEnergy)
		Error("Alignment is diverging instead of converging.\n");
	      
	      if ((iter < 100 || nIncrease < 2) && fixedDamping == 0.0)
		{
		  dampingFactor *= 0.5;
		  if (dampingFactor < 0.000000001)
		    Error("dampingFactor became too small\n");
		}
	    }
	  else
	    {
	      ++nDecrease;
	      nIncrease = 0;
	      if (fixedDamping == 0.0 && dampingFactor < 0.5)
		  dampingFactor *= 1.01;
	    }

	  /* 32 is chosen as a preferred time to output or terminate because,
	     given that we increase the damping by 1% on each iteration, and decrease
	     it by a factor of 2 when instability sets in, the number of iterations
	     between instabilities is about 70, and we want to output/terminate
	     when the state is not near an instability */
	  if (outputRequestedIter >= 0 &&
	      (nDecrease == 32 || iter > outputRequestedIter + 128))
	    {
	      Output(level, iter+1);
	      outputRequestedIter = -1;
	    }
	  if (terminationRequestedIter >= 0 &&
	      (nDecrease == 32 || iter > terminationRequestedIter + 128))
	    break;
	  if (refinementRequestedIter >= 0 &&
	      (nDecrease == 32 || iter > refinementRequestedIter + 128))
	    {
	      Log("Refinement (requested at iter %d, nDecrease = %d)\n",
		  refinementRequestedIter, nDecrease);
	      break;
	    }
	}

      if (foldDetected)
	{
	  Log("Detection of fold causing program to exit from step %d\n",
	      step);
	  break;
	}

      if (p == 0)
	Log("Finished alignment at step %d (level %d).\n", step, level);
      
      if (terminationRequestedIter >= 0)
	break;

      if (outputIncremental && outputName[0] != '\0')
	{
	  Log("OUTPUTTING SECTIONS\n");
	  Output(level, -1);
	}
      if (outputIncremental && outputGridName[0] != '\0')
	{
	  Log("OUTPUTTING GRIDS\n");
	  OutputGrid(level, -1);
	}
    }

  if (!outputIncremental && outputName[0] != '\0')
    {
      Log("OUTPUTTING SECTIONS\n");
      Output(level, -1);
    }
  if (!outputIncremental && outputGridName[0] != '\0')
    {
      Log("OUTPUTTING GRIDS\n");
      OutputGrid(level, -1);
    }

  if (foldDetected && outputGridName[0] != '\0')
    {
      Log("OUTPUTTING FOLDED GRID\n");

      if (MPI_Bcast(foldImage, PATH_MAX, MPI_CHAR,
		    processWithFold, MPI_COMM_WORLD) != MPI_SUCCESS ||
	  MPI_Bcast(&foldWidth, 1, MPI_INT,
		    processWithFold, MPI_COMM_WORLD) != MPI_SUCCESS ||
	  MPI_Bcast(&foldHeight, 1, MPI_INT,
		    processWithFold, MPI_COMM_WORLD) != MPI_SUCCESS ||
	  MPI_Bcast(&foldX, 1, MPI_INT,
		    processWithFold, MPI_COMM_WORLD) != MPI_SUCCESS ||
	  MPI_Bcast(&foldY, 1, MPI_INT,
		    processWithFold, MPI_COMM_WORLD) != MPI_SUCCESS)
	Error("Could not broadcast fold position.\n");
      strcpy(outputGridFocusImage, foldImage);
      outputGridFocusWidth = foldWidth;
      outputGridFocusHeight = foldHeight;
      outputGridFocusX = foldX;
      outputGridFocusY = foldY;
      outputGridSequence = 0;
      OutputGrid(level, -2);
    }

  Log("FINALIZING\n");
  MPI_Finalize();
  fclose(logFile);
  return(0);
}

void
Output (int level, int iter)
{
  int x, y;
  int i;
  char fn[PATH_MAX];
  char dirName[PATH_MAX];
  Point *pt;
  int nx, ny;
  MapElement *map;
  char msg[PATH_MAX+256];

  if (iter < 0)
    sprintf(dirName, "%s", outputName);
  else
    sprintf(dirName, "%sl%0.2di%0.6d",
	    outputName, level, iter);
  for (i = myFirstImage; i <= myLastImage; ++i)
    {
      nx = images[i].nx;
      ny = images[i].ny;
      Log("Going to output section %s  nx = %d ny = %d\n",
	  images[i].name, nx, ny);

      if (outputName[0] != '\0')
	{
	  sprintf(fn, "%s/%s.map", dirName, images[i].name);
	  if (!CreateDirectories(fn))
	    Error("Could not create directories for %s\n", fn);
	  map = malloc(nx * ny * sizeof(MapElement));
	  pt = images[i].points[startLevel-level];
	  for (y = 0; y < ny; ++y)
	    for (x = 0; x < nx; ++x)
	      {
		map[y * nx + x].x = pt->x;
		map[y * nx + x].y = pt->y;
		map[y * nx + x].c = 1.0;
		++pt;
	      }
	  if (!WriteMap(fn, map, level, nx, ny, 0, 0, images[i].name, "align12",
			UncompressedMap, msg))
	    Error("Could not write map file %s: %s\n", fn, msg);
	  free(map);
	}
    }
}

void
OutputGrid (int level, int iter)
{
  char gridDirName[PATH_MAX];  
  int x, y;
  int i, j, k;
  char fn[PATH_MAX];
  FILE *f;
  Point *pt;
  float rx, ry;
  int irx, iry;
  float rrx, rry;
  float rx00, rx01, rx10, rx11;
  float ry00, ry01, ry10, ry11;
  int nx, ny;
  int ix, iy;
  int dx, dy;
  int xp, yp;
  unsigned char *grid;
  int factor;
  MapElement *map;
  char msg[PATH_MAX+256];
  float nomD, dist;
  float hue;
  unsigned char r, g, b;
  int n;
  int hv;
  float xv, yv;
  int ixv, iyv;
  float bMinX, bMaxX, bMinY, bMaxY;
  int code;
  static float minX, maxX, minY, maxY;
  static float scale;
  static float offsetX, offsetY;

  Log("OutputGrid called at level = %d iter = %d\n", level, iter);
  factor = 1 << level;
  if (outputGridSequence == 0)
    {
      /* determine the scaling and offset for the output grids */
      minX = 1.0e30;
      maxX = -1.0e30;
      minY = 1.0e30;
      maxY = -1.0e30;
      if (outputGridFocusImage[0] != '\0')
	{
	  hv = Hash(outputGridFocusImage) % nImages;
	  for (i = imageHashTable[hv]; i >= 0; i = images[i].next)
	    if (strcmp(images[i].name, outputGridFocusImage) == 0)
	      {
		if (images[i].owner != p)
		  break;
		nx = images[i].nx;
		ny = images[i].ny;
		pt = images[i].points[startLevel-level];
		xv = ((double) outputGridFocusX) / factor;
		yv = ((double) outputGridFocusY) / factor;
		ixv = (int) floor(xv);
		iyv = (int) floor(yv);
		rrx = xv - ixv;
		rry = yv - iyv;
		while (ixv < 0)
		  {
		    ++ixv;
		    rrx -= 1.0;
		  }
		while (ixv >= nx-1)
		  {
		    --ixv;
		    rrx += 1.0;
		  }
		while (iyv < 0)
		  {
		    ++iyv;
		    rry -= 1.0;
		  }
		while (iyv >= ny-1)
		  {
		    --iyv;
		    rry += 1.0;
		  }
		rx00 = pt[iyv*nx+ixv].x;
		ry00 = pt[iyv*nx+ixv].y;
		rx01 = pt[(iyv+1)*nx+ixv].x;
		ry01 = pt[(iyv+1)*nx+ixv].y;
		rx10 = pt[iyv*nx+ixv+1].x;
		ry10 = pt[iyv*nx+ixv+1].y;
		rx11 = pt[(iyv+1)*nx+ixv+1].x;
		ry11 = pt[(iyv+1)*nx+ixv+1].y;
		rx = rx00 * (rrx - 1.0) * (rry - 1.0)
		  - rx10 * rrx * (rry - 1.0) 
		  - rx01 * (rrx - 1.0) * rry
		  + rx11 * rrx * rry;
		ry = ry00 * (rrx - 1.0) * (rry - 1.0)
		  - ry10 * rrx * (rry - 1.0) 
		  - ry01 * (rrx - 1.0) * rry
		  + ry11 * rrx * rry;
		minX = rx - outputGridFocusWidth / 2.0 / factor;
		maxX = rx + outputGridFocusWidth / 2.0 / factor;
		minY = ry - outputGridFocusHeight / 2.0 / factor;
		maxY = ry + outputGridFocusHeight / 2.0 / factor; 
		break;
	      }
	  if (i < 0)
	    Error("Could not find focus image %s\n", outputGridFocusImage);
	}
      else
	for (i = myFirstImage; i <= myLastImage; ++i)
	  {
	    nx = images[i].nx;
	    ny = images[i].ny;
	    pt = images[i].points[startLevel-level];
	    for (y = 0; y < ny; ++y)
	      for (x = 0; x < nx; ++x)
		{
		  if (pt->x < minX)
		    minX = pt->x;
		  if (pt->x > maxX)
		    maxX = pt->x;
		  if (pt->y < minY)
		    minY = pt->y;
		  if (pt->y > maxY)
		    maxY = pt->y;
		  ++pt;
		}
	  }
      if (MPI_Allreduce(MPI_IN_PLACE, &minX, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD) != MPI_SUCCESS ||
	  MPI_Allreduce(MPI_IN_PLACE, &maxX, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD) != MPI_SUCCESS ||	  
	  MPI_Allreduce(MPI_IN_PLACE, &minY, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD) != MPI_SUCCESS ||
	  MPI_Allreduce(MPI_IN_PLACE, &maxY, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD) != MPI_SUCCESS)
	Error("MPI_Allreduce: could not compute grid limits.\n");
      
      if (minX == maxX || minY == maxY)
	Error("maps have 0 extent in x or y\n");
      minX *= factor;
      maxX *= factor;
      minY *= factor;
      maxY *= factor;
      //      printf("After allreduce: %f %f %f %f\n", minX, maxX, minY, maxY);
      scale = 0.95 * outputGridWidth / (maxX - minX);
      if (0.95 * outputGridHeight / (maxY - minY) < scale)
	scale = 0.95 * outputGridHeight / (maxY - minY);
      offsetX = 0.5 * outputGridWidth - scale * 0.5 * (maxX + minX);
      offsetY = 0.5 * outputGridHeight - scale * 0.5 * (maxY + minY);
      printf("ogw = %d ogh = %d scale = %f offsetX = %f offsetY = %f\n",
	     outputGridWidth, outputGridHeight, scale, offsetX, offsetY);
    }
  if (outputGridOverlay || iter < 0)
    {
      if (iter == -2)
	sprintf(gridDirName, "%sfold", outputGridName);
      else
	sprintf(gridDirName, "%s", outputGridName);
    }
  else
    sprintf(gridDirName, "%sl%0.2di%0.6d",
	    outputGridName, level, iter);

  n = outputGridWidth * outputGridHeight * 3;
  grid = (unsigned char *) malloc(n);
  if (outputGridOverlay)
    memset(grid, 255, n);

  bMinX = (minX - 0.25 * (maxX - minX)) / factor;
  bMaxX = (maxX + 0.25 * (maxX - minX)) / factor;
  bMinY = (minY - 0.25 * (maxY - minY)) / factor;
  bMaxY = (maxY + 0.25 * (maxY - minY)) / factor;
  
  for (i = myFirstImage; i <= myLastImage; ++i)
    {
      if (!outputGridOverlay)
	memset(grid, 255, n);
      nx = images[i].nx;
      ny = images[i].ny;
      pt = images[i].points[startLevel-level];
      if (pt[0].x < bMinX && pt[nx-1].x < bMinX &&
	  pt[(ny-1)*nx].x < bMinX && pt[(ny-1)*nx+nx-1].x < bMinX)
	continue;
      if (pt[0].x > bMaxX && pt[nx-1].x > bMaxX &&
	  pt[(ny-1)*nx].x > bMaxX && pt[(ny-1)*nx+nx-1].x > bMaxX)
	continue;
      if (pt[0].y < bMinY && pt[nx-1].y < bMinY &&
	  pt[(ny-1)*nx].y < bMinY && pt[(ny-1)*nx+nx-1].y < bMinY)
	continue;
      if (pt[0].y > bMaxY && pt[nx-1].y > bMaxY &&
	  pt[(ny-1)*nx].y > bMaxY && pt[(ny-1)*nx+nx-1].y > bMaxY)
	continue;
      for (y = 0; y < ny; ++y)
	for (x = 0; x < nx; ++x)
	  for (dy = -1; dy <= 0; ++dy)
	    {
	      yp = y + dy;
	      if (yp < 0)
		continue;
	      for (dx = -1; dx <= 1; ++dx)
		{
		  xp = x + dx;
		  if (xp < 0)
		    continue;
		  if (xp >= nx || dy == 0 && dx >= 0)
		    break;
			
#if 0
		  nomD = (dx * dy != 0) ? M_SQRT2 : 1.0;
		  dist = hypot(pt[y*nx+x].x - pt[yp*nx+xp].x,
			       pt[y*nx+x].y - pt[yp*nx+xp].y);
		  hue = 4.0 * (dist / nomD) - 2.0;
		  if (hue < 0.0)
		    hue = 0.0;
		  if (hue > 5.0)
		    hue = 5.0;
		  hsv_to_rgb(&r, &g, &b, hue, 1.0, 1.0);
#else
		  r = colors[i & 15][0];
		  g = colors[i & 15][1];
		  b = colors[i & 15][2];
#endif

		  DrawLine(grid,
			   scale * factor * pt[y*nx+x].x + offsetX,
			   scale * factor * pt[y*nx+x].y + offsetY,
			   scale * factor * pt[yp*nx+xp].x + offsetX,
			   scale * factor * pt[yp*nx+xp].y + offsetY,
			   outputGridWidth, outputGridHeight,
			   r, g, b);
		}
	    }

      for (y = 0; y < ny; ++y)
	for (x = 0; x < nx; ++x)
	  {
	    if (pt[y*nx+x].fx > 0.5 * CONSTRAINED ||
		pt[y*nx+x].fy < 0.5 * CONSTRAINED)
	      {
		ix = (int) (scale * factor * pt[y*nx+x].x + offsetX);
		iy = (int) (scale * factor * pt[y*nx+x].y + offsetY);
		for (dy = -1; dy <= 1; ++dy)
		  for (dx = -1; dx <= 1; ++dx)
		    {
		      xp = ix + dx;
		      yp = iy + dy;
		      if (xp >= 0 && xp < outputGridWidth &&
			  yp >= 0 && yp < outputGridHeight)
			{
			  grid[3*(yp * outputGridWidth + xp)] = 0;
			  grid[3*(yp * outputGridWidth + xp) + 1] = 0;
			  grid[3*(yp * outputGridWidth + xp) + 2] = 0;
			}
		    }
	      }
	  }

      if (!outputGridOverlay)
	{
	  sprintf(fn, "%s/%s.pnm", gridDirName, images[i].name);
	  if (!CreateDirectories(fn))
	    Error("Could not create directories for %s\n", fn);
	  f = fopen(fn, "w");
	  if (f == NULL)
	    Error("Could not open output grid file %s\n", fn);
	  fprintf(f, "P6\n%d %d\n%d\n", outputGridWidth, outputGridHeight, 255);
	  if (fwrite(grid, sizeof(unsigned char), n, f) != n)
	    Error("Could not write to grid file %s\n", fn);
	  fclose(f);
	}
    }

  if (outputGridOverlay)
    {
      // for some reason MPI_Reduce give an error with MPI_IN_PLACE
      if ((code = MPI_Allreduce(MPI_IN_PLACE, grid,  n, MPI_UNSIGNED_CHAR,
				MPI_MIN, MPI_COMM_WORLD)) != MPI_SUCCESS)
	Error("MPI_reduce of grid image failed: %d\n", code);
	  
      if (p == 0)
	{
	  sprintf(fn, "%s%0.6d.pnm", gridDirName, outputGridSequence);
	  if (!CreateDirectories(fn))
	    Error("Could not create directories for %s\n", fn);
	  f = fopen(fn, "w");
	  if (f == NULL)
	    Error("Could not open output grid file %s\n", fn);
	  fprintf(f, "P6\n%d %d\n%d\n", outputGridWidth, outputGridHeight, 255);
	  if (fwrite(grid, sizeof(unsigned char), n, f) != n)
	    Error("Could not write to grid file %s\n", fn);
	  fclose(f);
	}
    }
  free(grid);
  ++outputGridSequence;
}

void Error (char *fmt, ...)
{
  va_list args;

  if (logFile != NULL)
    {
      va_start(args, fmt);
      fprintf(logFile, "%f: ERROR: ", MPI_Wtime());
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

  va_start(args, fmt);
  fprintf(logFile, "%f: ", MPI_Wtime());
  vfprintf(logFile, fmt, args);
  va_end(args);
  fflush(logFile);
}

void hsv_to_rgb (unsigned char* r, unsigned char *g, unsigned char *b,
		 float h, float s, float v)
{ 
  /* H is given on [0, 6] or UNDEFINED. S and V are given on [0, 1]. 
     RGB are each returned on [0, 255]. */
#define RETURN_RGB(rd,gr,bl)	{ *r = floor(255.999999 * rd); *g = floor(255.999999 * gr); *b = floor(255.999999 * bl); return; }
  
  float m, n, f;
  int i; 
  
  i = floor(h);
  f = h - i; 
  if ( !(i&1) )
    f = 1 - f; /* if i is even */
  m = v * (1 - s); 
  n = v * (1 - s * f); 
  
  switch (i)
    { 
    case 6: 
    case 0: RETURN_RGB(v, n, m); 
    case 1: RETURN_RGB(n, v, m); 
    case 2: RETURN_RGB(m, v, n);
    case 3: RETURN_RGB(m, n, v); 
    case 4: RETURN_RGB(n, m, v); 
    case 5: RETURN_RGB(v, m, n); 
    } 
} 

void
DrawLine(unsigned char *img, float px0, float py0, float px1, float py1,
	 int imageWidth, int imageHeight,
	 unsigned char r, unsigned char g, unsigned char b)
{
  int x0, y0, x1, y1;
  int steep;
  int tmp;
  int error;
  int yStep;
  int x, y;
  int dx, dy;
  float dist;
  float hue;
  unsigned char *p;

  x0 = (int) floor(px0);
  y0 = (int) floor(py0);
  x1 = (int) floor(px1);
  y1 = (int) floor(py1);

  // use the Bresenham Algorithm for lines
  //    (derived from the Wikipedia entry on
  //     "Bresenham's line algorithm")
  steep = abs(y1 - y0) > abs(x1 - x0);
  if (steep)
    {
      tmp = x0;
      x0 = y0;
      y0 = tmp;
      tmp = x1;
      x1 = y1;
      y1 = tmp;
    }
  if (x0 > x1)
    {
      tmp = x0;
      x0 = x1;
      x1 = tmp;
      tmp = y0;
      y0 = y1;
      y1 = tmp;
    }
  dx = x1 - x0;
  dy = abs(y1 - y0);
  error = -dx / 2;
  y = y0;
  if (y0 < y1)
    yStep = 1;
  else
    yStep = -1;
  for (x = x0; x <= x1; ++x)
    {
      if (steep)
	{
	  if (y >= 0 && y < imageWidth &&
	      x >= 0 && x < imageHeight)
	    {
	      p = &img[3*(x * imageWidth + y)];
	      *p++ = r;
	      *p++ = g;
	      *p++ = b;
	    }
	}
      else
	{
	  if (x >= 0 && x < imageWidth &&
	      y >= 0 && y < imageHeight)
	    {
	      p = &img[3*(y * imageWidth + x)];
	      *p++ = r;
	      *p++ = g;
	      *p++ = b;
	    }
	}
      error += dy;
      if (error > 0)
	{
	  y += yStep;
	  error -= dx;
	}
    }
}

void
UpdateDeltas (int level, int iter)
{
  InterImageMap *m;
  int i, j, k;
  int x, y;
  float rrx, rry;
  Point *pts1;
  float rx00, rx01, rx10, rx11;
  float ry00, ry01, ry10, ry11;
  float rx, ry;
  int nx, ny;
  int irx, iry;
  int ix, iy;
  int nx1, ny1;
  InterImageSpring *s;
  InterImageStrip *strip;
  float deltaX, deltaY;
  int springsPos;
  int ns;
  int nStrips;
  InterImageStrip *strips;
  InterImageSpring *springs;

  for (i = 0; i < nMaps; ++i)
    {
      m = &maps[i];
      nx1 = images[m->image1].nx;
      ny1 = images[m->image1].ny;
      pts1 = images[m->image1].points[startLevel-level];

      springsPos = 0;
      nStrips = m->nStrips[startLevel - level];
      strips = m->strips[startLevel - level];
      springs = m->springs[startLevel - level];
      for (j = 0; j < nStrips; ++j)
	{
	  strip = &(strips[j]);
	  ns = strip->nSprings;
	  x = strip->x0;
	  y = strip->y0;
	  ix = strip->x1;
	  iy = strip->y1;
	  for (k = 0; k < ns; ++k)
	    {
	      s = &(springs[springsPos++]);
	      ix += ((s->dxy1) >> 4) - 8;
	      iy += ((s->dxy1) & 0xf) - 8;
	      irx = ix;
	      iry = iy;
	      rrx = s->irrx / 65536.0;
	      rry = s->irry / 65536.0;
	      if (rrx >= 0.0)
		{
		  if (irx >= nx1-1)
		    {
		      --irx;
		      rrx += 1.0;
		    }
		}
	      else
		{
		  if (irx > 0)
		    {
		      --irx;
		      rrx += 1.0;
		    }
		}
	      if (rry >= 0.0)
		{
		  if (iry >= ny1-1)
		    {
		      --iry;
		      rry += 1.0;
		    }
		}
	      else
		{
		  if (iry > 0)
		    {
		      --iry;
		      rry += 1.0;
		    }
		}
	      rx00 = pts1[iry*nx1 + irx].x;
	      rx01 = pts1[(iry+1)*nx1+irx].x;
	      rx10 = pts1[iry*nx1 + irx+1].x;
	      rx11 = pts1[(iry+1)*nx1+irx+1].x;
	      ry00 = pts1[iry*nx1 + irx].y;
	      ry01 = pts1[(iry+1)*nx1+irx].y;
	      ry10 = pts1[iry*nx1 + irx+1].y;
	      ry11 = pts1[(iry+1)*nx1+irx+1].y;
	      rx = rx00 * (rrx - 1.0) * (rry - 1.0)
		- rx10 * rrx * (rry - 1.0) 
		- rx01 * (rrx - 1.0) * rry
		+ rx11 * rrx * rry;
	      ry = ry00 * (rrx - 1.0) * (rry - 1.0)
		- ry10 * rrx * (rry - 1.0) 
		- ry01 * (rrx - 1.0) * rry
		+ ry11 * rrx * rry;
	      deltaX = 10000.0 * (rx - pts1[iy*nx1+ix].x) + 0.5;
	      deltaY = 10000.0 * (ry - pts1[iy*nx1+ix].y) + 0.5;
	      if (deltaX < -32768.0 || deltaX >= 32768.0 ||
		  deltaY < -32768.0 || deltaY >= 32768.0)
		Error("Internal error: pixel offsets too large: dx = %f dy = %f\n",
		      deltaX, deltaY);
	      s->dx = (int) floor(deltaX);
	      s->dy = (int) floor(deltaY);
	    }
	}
    }
}

void
RefinePositions (int prevLevel, int level)
{
  int i;
  int factor;
  int factor1;
  int nx, ny;
  int nx1, ny1;
  int x, y;
  float xv, yv;
  int ixv, iyv;
  float rrx, rry;
  float rx00, rx01, rx10, rx11;
  float ry00, ry01, ry10, ry11;
  float rx, ry;
  Point *pts, *pts1;
  float dFactor;

  factor = 1 << prevLevel;
  factor1 = 1 << level;
  dFactor = 1 << (prevLevel - level);
  for (i = 0; i < nImages; ++i)
    {
      if (images[i].owner != p && !images[i].needed)
	continue;
      nx = (images[i].width + factor - 1) / factor + 1;
      ny = (images[i].height + factor - 1) / factor + 1;
      nx1 = (images[i].width + factor1 - 1) / factor1 + 1;
      ny1 = (images[i].height + factor1 - 1) / factor1 + 1;
      images[i].nx = nx1;
      images[i].ny = ny1;
      pts = images[i].points[startLevel-prevLevel];
      pts1 = images[i].points[startLevel-level];
      if (images[i].owner != p)
	{
	  memset(pts1, 0, nx1 * ny1 * sizeof(Point));
	  continue;
	}
      for (y = 0; y < ny1; ++y)
	for (x = 0; x < nx1; ++x)
	  {
	    xv = x / dFactor;
	    yv = y / dFactor;
	    ixv = (int) floor(xv);
	    iyv = (int) floor(yv);
	    rrx = xv - ixv;
	    rry = yv - iyv;
	    while (ixv >= nx-1)
	      {
		--ixv;
		rrx += 1.0;
	      }
	    while (iyv >= ny-1)
	      {
		--iyv;
		rry += 1.0;
	      }
	    rx00 = pts[iyv * nx + ixv].x;
	    ry00 = pts[iyv * nx + ixv].y;
	    rx01 = pts[(iyv+1) * nx + ixv].x;
	    ry01 = pts[(iyv+1) * nx + ixv].y;
	    rx10 = pts[iyv * nx + ixv + 1].x;
	    ry10 = pts[iyv * nx + ixv + 1].y;
	    rx11 = pts[(iyv+1) * nx + ixv + 1].x;
	    ry11 = pts[(iyv+1) * nx + ixv + 1].y;
	    rx = rx00 * (rrx - 1.0) * (rry - 1.0)
	      - rx10 * rrx * (rry - 1.0) 
	      - rx01 * (rrx - 1.0) * rry
	      + rx11 * rrx * rry;
	    ry = ry00 * (rrx - 1.0) * (rry - 1.0)
	      - ry10 * rrx * (rry - 1.0) 
	      - ry01 * (rrx - 1.0) * rry
	      + ry11 * rrx * rry;

	    pts1[y * nx1 + x].x = dFactor * rx;
	    pts1[y * nx1 + x].y = dFactor * ry;
	    pts1[y * nx1 + x].fx = 0.0;
	    pts1[y * nx1 + x].fy = 0.0;
	  }
    }
}

void
PlanCommunications (int level)
{
  int phase;
  CommPhase *cp;
  int bufferPos;
  int i, j, k;
  int nFloats;

  Log("Planning communication for level %d\n", level);
  for (phase = 0; phase < nPhases; ++phase)
    {
      // schedule the communication of positions between processors
      cp = &(commPhases[phase]);

      // first, the sends
      cp->nSends = 0;
      bufferPos = 0;
      k = 0;
      for (i = 0; i < cp->nSendImages; ++i)
	{
	  j = cp->sendImages[i];
	  nFloats = images[j].nx * images[j].ny * 2;
	  if (bufferPos + nFloats > bufferSize)
	    {
	      ++(cp->nSends);
	      cp->nImagesToSend = (int *) realloc(cp->nImagesToSend,
						  cp->nSends * sizeof(int));
	      cp->nImagesToSend[cp->nSends-1] = k;
	      Log("In phase %d will send %d images to %d in a message of %d floats\n",
		  phase, k, cp->otherProcess, bufferPos);
	      k = 0;
	      bufferPos = 0;
	    }
	  bufferPos += nFloats;
	  ++k;
	}
      if (k > 0)
	{
	  ++(cp->nSends);
	  cp->nImagesToSend = (int *) realloc(cp->nImagesToSend,
					      cp->nSends * sizeof(int));
	  cp->nImagesToSend[cp->nSends-1] = k;
	  Log("In phase %d will send %d images to %d in a final message of %d floats\n",
	      phase, k, cp->otherProcess, bufferPos);
	}

      // next, the receives
      cp->nReceives = 0;
      bufferPos = 0;
      k = 0;
      for (i = 0; i < cp->nReceiveImages; ++i)
	{
	  j = cp->receiveImages[i];
	  nFloats = images[j].nx * images[j].ny * 2;
	  if (bufferPos + nFloats > bufferSize)
	    {
	      ++(cp->nReceives);
	      cp->floatsToReceive = (int *) realloc(cp->floatsToReceive,
						    cp->nReceives * sizeof(int));
	      cp->nImagesToReceive = (int *) realloc(cp->nImagesToReceive,
						  cp->nReceives * sizeof(int));
	      cp->floatsToReceive[cp->nReceives-1] = bufferPos;
	      cp->nImagesToReceive[cp->nReceives-1] = k;
	      Log("In phase %d will receive %d images from %d in a message of %d floats\n",
		  phase, k, cp->otherProcess, bufferPos);
	      k = 0;
	      bufferPos = 0;
	    }
	  bufferPos += nFloats;
	  ++k;
	}
      if (k > 0)
	{
	  ++(cp->nReceives);
	  cp->floatsToReceive = (int *) realloc(cp->floatsToReceive,
						cp->nReceives * sizeof(int));
	  cp->nImagesToReceive = (int *) realloc(cp->nImagesToReceive,
						 cp->nReceives * sizeof(int));
	  cp->floatsToReceive[cp->nReceives-1] = bufferPos;
	  cp->nImagesToReceive[cp->nReceives-1] = k;
	  Log("In phase %d will receive %d images from %d in a final message of %d floats\n",
	      phase, k, cp->otherProcess, bufferPos);
	}
    }
}


void
CommunicatePositions (int level)
{
  /* communicate image positions to other nodes */
  //      (min(a,b) / sep) & 1 == 0  : first phase
  //      (min(a,b) / sep) & 1 == 1  : second phase
  
  //	  0-1,2-3,4-5,6-7,8-9,10-11,12-13      0,2,4,6,8,10,12
  //        1-2,3-4,5-6,7-8,9-10,11-12,13-14   1,3,5,7,9,11,13
  //      sep 1, 

  //      0-2,1-3,4-6,5-7,8-10,9-11,12-14      0,1,4,5,8,9,12
  //      2-4,3-5,6-8,7-9,10-12,11-13          2,3,6,7,10,11
  //      sep 2

  //      0-3,1-4,2-5,6-9,7-10,8-11            0,1,2,6,7,8
  //      3-6,4-7,5-8,9-12,10-13,11-14         3,4,5,9,10,11
  
  //      0-4,1-5,2-6,3-7,8-12,9-13,10-14
  //      4-8,5-9,6-10,7-11

  //      0-5,1-6,2-7,3-8,4-9
  //      5-10,6-11,7-12,8-13,9-14

  //      0-6,1-7,2-8,3-9,4-10,5-11
  //      6-12,7-13,8-14

  //      0-7,1-8,2-9,3-10,4-11,5-12,6-13
  //      7-14

  //      0-8,1-9,2-10,3-11,4-12,5-13,6-14
  //      none

  //      0-9,1-10,2-11,3-12,4-13,5-14
  //      none

  //      0-10,1-11,2-12,3-13,4-14
  //      none

  //      0-11,1-12,2-13,3-14
  //      none

  //      0-12,1-13,2-14
  //      none

  //      0-13,1-14
  //      none

  //      0-14
  //      none
  int phase;
  int op;
  int subPhase;
  int i, j, k;
  int jj;
  int bufferPos;
  int nPts;
  int its;
  Point *pts;
  MPI_Status status;

  for (phase = 0; phase < nPhases; ++phase)
    {
      if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS)
	Error("MPI_Barrier() failed.\n");
      op = commPhases[phase].otherProcess;
      if (op < 0)
	continue;
      for (subPhase = 0; subPhase < 2; ++subPhase)
	{
	  if (subPhase ^ (p < op))
	    {
	      /* send */
	      for (i = 0; i < commPhases[phase].nSends; ++i)
		{
		  jj = 0;
		  bufferPos = 0;
		  for (j = 0; j < commPhases[phase].nImagesToSend[i]; ++j, ++jj)
		    {
		      its = commPhases[phase].sendImages[jj];
		      nPts = images[its].nx * images[its].ny;
		      pts = images[its].points;
		      for (k = 0; k < nPts; ++k)
			{
			  buffer[bufferPos++] = pts[k].x;
			  buffer[bufferPos++] = pts[k].y;
			}
		    }
		  if (MPI_Send(buffer, bufferPos, MPI_FLOAT,
			       op, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
		    Error("Could not send to process %d\n", op);
		}
	    }
	  else
	    {
	      /* receive */
	      for (i = 0; i < commPhases[phase].nReceives; ++i)
		{
		  if (MPI_Recv(buffer, commPhases[phase].floatsToReceive[i], MPI_FLOAT,
			       op, 0, MPI_COMM_WORLD, &status) != MPI_SUCCESS)
		    Error("Could not receive from process %d\n", op);
		  jj = 0;
		  bufferPos = 0;
		  for (j = 0; j < commPhases[phase].nImagesToReceive[i]; ++j, ++jj)
		    {
		      its = commPhases[phase].receiveImages[jj];
		      nPts = images[its].nx * images[its].ny;
		      pts = images[its].points;
		      for (k = 0; k < nPts; ++k)
			{
			  pts[k].x =  buffer[bufferPos++];
			  pts[k].y =  buffer[bufferPos++];
			}
		    }
		}
	    }
	}
    }
}

unsigned int
Hash (char *s)
{
  unsigned int v = 0;
  char *p = s;
  while (*p != '\0')
    v = 37 * v + *p++;
  return(v);
}

unsigned int
HashMap (char *s, int nx, int ny)
{
  unsigned int v = 0;
  char *p = s;
  while (*p != '\0')
    v = 37 * v + *p++;
  v = 37 * v + nx;
  v = 37 * v + ny;
  return(v);
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
      
      if (mkdir(dn, 0755) != 0 && errno != EEXIST)
	{
	  Log("Could not create directory %s\n", dn);
	  return(0);
	}
    }
  return(1);
}

int
CheckForFolds (int level, int iter, int reportError)
{
  int i;
  int nx, ny;
  int x, y;
  float xv, yv;
  Point *pts;
  float ax, ay;
  float bx, by;

  for (i = myFirstImage; i <= myLastImage; ++i)
    {
      nx = images[i].nx;
      ny = images[i].ny;
      pts = images[i].points;
      for (y = 0; y < ny; ++y)
	for (x = 0; x < nx; ++x)
	  {
	    xv = pts[y * nx + x].x;
	    yv = pts[y * nx + x].y;
	    if (x > 0)
	      {
		ax = pts[y * nx + x - 1].x - xv;
		ay = pts[y * nx + x - 1].y - yv;
		if (y > 0)
		  {
		    bx = pts[(y-1) * nx + x].x - xv;
		    by = pts[(y-1) * nx + x].y - yv;
		    if (ax * by - ay * bx <= 0.0)
		      goto foldError;
		  }
		if (y < ny - 1)
		  {
		    bx = pts[(y+1) * nx + x].x - xv;
		    by = pts[(y+1) * nx + x].y - yv;
		    if (ax * by - ay * bx >= 0.0)
		      goto foldError;
		  }
	      }
	    if (x < nx - 1)
	      {
		ax = pts[y * nx + x + 1].x - xv;
		ay = pts[y * nx + x + 1].y - yv;
		if (y > 0)
		  {
		    bx = pts[(y-1) * nx + x].x - xv;
		    by = pts[(y-1) * nx + x].y - yv;
		    if (ax * by - ay * bx >= 0.0)
		      goto foldError;
		  }
		if (y < ny - 1)
		  {
		    bx = pts[(y+1) * nx + x].x - xv;
		    by = pts[(y+1) * nx + x].y - yv;
		    if (ax * by - ay * bx <= 0.0)
		      goto foldError;
		  }
	      }
	  }
    }
  return(0);

 foldError:
  strcpy(foldImage, images[i].name);
  foldWidth = (1 << level) * 64;
  foldHeight = (1 << level) * 64;
  foldX = (1 << level) * ((double) x);
  foldY = (1 << level) * ((double) y);
  if (reportError)
    Error("Fold detected in map of image %s at (%f %f) (%d(%d %d)) at iteration %d\n",
	images[i].name,
	(1 << level) * ((double) x), (1 << level) * ((double) y),
	level, x, y, iter);
  else
    Log("Fold detected in map of image %s at (%f %f) (%d(%d %d)) at iteration %d\n",
	images[i].name,
	(1 << level) * ((double) x), (1 << level) * ((double) y),
	level, x, y, iter);
  return(1);
}

int
Extrapolate (float *prx, float *pry, int ix, int iy, float arrx, float arry,
	     MapElement* map, int mw, int mh, float threshold)
{
  float rx = 0.0;
  float ry = 0.0;
  float totalWeight = 0.0;
  int maxDelta = ((int) ceil(threshold)) + 2;
  float rrx, rry;
  int dx, dy;
  int ixv, iyv;
  float rx00, rx01, rx10, rx11;
  float ry00, ry01, ry10, ry11;
  float weight;

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
	  rx01 = map[(iyv+1)*mw+ixv].x;
	  ry01 = map[(iyv+1)*mw+ixv].y;
	  rx10 = map[iyv*mw+ixv+1].x;
	  ry10 = map[iyv*mw+ixv+1].y;
	  rx11 = map[(iyv+1)*mw+ixv+1].x;
	  ry11 = map[(iyv+1)*mw+ixv+1].y;

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
	  totalWeight += weight;
	}
    }
  if (totalWeight == 0.0)
    return(0);
  *prx = rx / totalWeight;
  *pry = ry / totalWeight;
  return(1);
}
