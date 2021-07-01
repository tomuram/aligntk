/*
 * align12.c  - align a set of maps using a relaxation algorithm
 *                on a system of springs
 *
 *   Copyright (c) 2006-2010 Pittsburgh Supercomputing Center
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
#include <limits.h>
#include <sys/resource.h>

#include "imio.h"
#include "dt.h"

#define DEBUG	0
#define PDEBUG	0

#define LINE_LENGTH	256
#define MAX_SPEC        1024
#define CONSTRAINED	(1.0e+30)
#define UNSPECIFIED	(1.0e+30)

typedef struct Step
{
  float level;      // actually an int, but we use a float to simplify
                    //  the MPI_Bcast
  float kIntra;
  float kInter;
  float kAbsolute;
  int minIter;
  float threshold;
  float dampingFactor;
  int clampX;
  int clampY;
} Step;

typedef struct Point
{
  float x, y;
} Point;

typedef struct Node
{
  float x;    /* current x position of this Node;
	         if x > 0.5*UNSPECIFIED, then Node is not valid
	         for the image */
  float y;    /* current y position of this Node */
  float fx;   /* force in x direction;
		 if fx > 0.5*CONSTRAINED, Node is not allowed to
		 move in the x direction */
  float fy;   /* force in y direction;
		 if fy > 0.5*CONSTRAINED, Node is not allowed to
		 move in the y direction */
} Node;

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
  int nx, ny;           /* number of nodes in x and y at current level */
  Node *nodes;          /* nodes in the image at current level */
  Node *initialNodes;   /* the initial node positions for the current step;
			   used only if recovery from a fold is necessary */
  Point **absolutePositions;
                        /* absolute positions at each level for each node */
  char *modelName;      /* name of map to use for initial positions and
			   distances */
  struct IntraImageMap *map;
                        /* map to use for intra-image springs */
  unsigned char *mask;  /* mask at the endLevel of blocks containing valid pixels */
  float kFactor;        /* factor to be applied to spring constants for
			   intra-image springs of this image */
  int marked;		/* true if we have marked this image as near a fold */
} Image;

typedef struct IntraImageMap
{
  char *name;           /* filename of map (or "" if uniform map) */
  struct IntraImageMap *next;  /* next in hash bucket */
  int nx, ny;           /* number of nodes in x and y at the ending level */
  Point *points;   /* points in the model at the starting level */
  int *nSprings;        /* number of springs at each level */
  struct IntraImageSpring **springs;
                        /* array of springs at each level */
} IntraImageMap;

typedef struct IntraImageSpring
{
  int index0;		/* index of first point in image */
  int index1;		/* index of second point in image */
  float nomD;           /* nominal distance between the nodes */
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
			   this is in units of 0.00005 * nominal spring spacing */
  unsigned char k;      /* spring constant (multiplied by 100) */
  unsigned char dxy1;   /* delta position in target:
			     delta_x = (dxy1 >> 4) - 8;
			     delta_y = (dxy1 & 0xf) - 8; */
  short irrx, irry;     /* goal offset from nearest grid point in
			   target map; this is in units of 0.00005 *
			   nominal spring spacing */
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

typedef struct SpringForce
{
  int delta;             /* ordering of map producing this spring */
  float theta;
  float force;
} SpringForce;


int p;      /* rank of this process */
int np;     /* number of processes in this run */

char schedule[PATH_MAX];
int params[12];
float fParams[5];
char imageListName[PATH_MAX];
char mapListName[PATH_MAX];
char mapsName[PATH_MAX];
char initialMapsName[PATH_MAX];
char absoluteMapsName[PATH_MAX];
char modelsName[PATH_MAX];
char masksName[PATH_MAX];
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
int outputGridFocusDepth = -1;
int outputGridFocusX = -1;
int outputGridFocusY = -1;
char outputGridFocusImage[PATH_MAX];
int outputGridSequence = 0;
int outputFoldMaps = 0;
char constraintName[PATH_MAX];
int nFixedImages = 0;
char **fixedImages = 0;
int showConstraints = 0;
int epochIterations = 512;

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

int startFactor;        /* spacing at starting level */
int endFactor;	        /* spacing at end level */
int factor;             /* spacing at current level */

int nSteps = 0;
Step *steps = NULL;

int foldRecovery = 0;
int foldRecoveryCount = 0;
int foldDetected = 0;
int minImageWithFold = -1;
int foldImage = -1;
int foldWidth = 0;
int foldHeight = 0;
int foldX = 0;
int foldY = 0;
float foldAbsoluteX = 0.0;
float foldAbsoluteY = 0.0;
float foldRadius = 32.0;

int *exfoldImage = 0;
int *exfoldWidth = 0;
int *exfoldHeight = 0;
int *exfoldX = 0;
int *exfoldY = 0;

int minimizeArea = 0;

#if PDEBUG
int ioi = 28;
float poixv = 10432.0;
float poiyv = 4032.0;

//int ioi = 29;
//float poixv = 11328.0;
//float poiyv = 4160.0;

//int ioi = 9;
//float poixv = 312.0;
//float poiyv = 15192.0;

// int ioi = 10;
// float poixv = 1272.0;
// float poiyv = 848.0;

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
void OutputFoldMaps (int level, int iter);
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
int CheckForFolds (int level, int iter);
void ComputeSpringConstants (int level, int iter);
int Extrapolate (float *prx, float *pry, float *prc,
		 int ix, int iy, float arrx, float arry,
		 MapElement* map, int mw, int mh, float threshold);
void RecoverFromFold (int level);
double PerpDist (double x1, double y1, double x2, double y2, double s);
void PointProjectionOnLine (Point p0, Point p1, Point q, Point *proj);
double Angle (double x0, double y0,
	      double x1, double y1,
	      double x2, double y2);
int ComparePoints (const void *p0, const void *p1);
double Cross (double x0, double y0,
	      double x1, double y1,
	      double x2, double y2);
void ConvexHull (int n, Point *pts,
		 int *hn, Point *hpts);
double MinimumAreaRectangle (int n, Point *pts, Point *corners);
int CompareSpringForces (const void *p0, const void *p1);

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
  Node *node;
  int dx, dy;
  float *fx;
  float *fy;
  double totalEnergy;
  double prevTotalEnergy;
  double deltaEnergy;
  double energy;
  double epochInitialTotalEnergy;
  double epochFinalTotalEnergy;
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
  int nNodes;
  int nFloats;
  Node *p00, *p10, *p01, *p11;
  float kIntraThisImage;
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
  double kFactor;
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
  Node *nodes0, *nodes1;
  Point *ipt;
  double cost, sint;
  Node *nodes;
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
  Point *mpts;
  float *mptsc;
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
  int imageWithFold;
  int processWithFold;
  int minIter;
  int clampX, clampY;
  Point *absPos;
  int rnx, rny;
  int valid;
  int ix, iy;
  int eFactor;
  int rFactor;
  int count;
  int count1, count2;
  unsigned char *imask;
  int imbpl;
  int maskWidth, maskHeight;
  int ixMin, ixMax, iyMin, iyMax;
  /* DECLS */

#if 0
  /* temp for testing */
  int npts = 32;
  Point pts[32];
  Point hpts[64];
  Point corners[4];
  int hn;
  int upperLeftCorner;
  double upperLeftDistance;
  double distance;
  double theta;
  double a[6];
  Point v0, v1, v2;
  double theta2;
  srand48(34237L);
  for (i = 0; i < npts; ++i)
    {
      pts[i].x = drand48();
      pts[i].y = drand48();
      printf("%d: %f %f\n", i, pts[i].x, pts[i].y);
    }
  pts[0].x = 1.9999;
  pts[0].y = -1.5;
  ConvexHull(npts, pts, &hn, hpts);
  printf("hn = %d\n", hn);
  for (i = 0; i < hn; ++i)
    printf("%i: %f %f\n", i, hpts[i].x, hpts[i].y);
  printf("area = %f\n", MinimumAreaRectangle(hn, hpts, corners));
  for (i = 0; i < 4; ++i)
    printf("Corner %d: %f %f\n", i, corners[i].x, corners[i].y);

  upperLeftCorner = -1;
  upperLeftDistance = 1.0e30;
  for (i = 0; i < 4; ++i)
    {
      distance = hypot(corners[i].x, corners[i].y);
      if (distance < upperLeftDistance)
	{
	  upperLeftDistance = distance;
	  upperLeftCorner = i;
	}
    }
  /* we want to translate the upper left to 0, then
     rotate so that the line from the upper left to the
     upper right ((upperLeftCorner + 3) % 4) is along the
     x-axis */
  theta = -atan2(corners[(upperLeftCorner + 3) % 4].y -
		corners[upperLeftCorner].y,
		corners[(upperLeftCorner + 3) % 4].x -
		corners[upperLeftCorner].x);
  printf("theta = %f\n", theta);
  theta2 = -atan2(corners[(upperLeftCorner + 1) % 4].y -
		corners[upperLeftCorner].y,
		corners[(upperLeftCorner + 1) % 4].x -
		corners[upperLeftCorner].x);
  printf("theta2 = %f\n", theta2);
  v0.x = corners[(upperLeftCorner + 3) % 4].x -
    corners[upperLeftCorner].x;
  v0.y = corners[(upperLeftCorner + 3) % 4].y -
    corners[upperLeftCorner].y;
  printf("v0 = %f %f\n", v0.x, v0.y);
  v1.x = corners[(upperLeftCorner + 1) % 4].x -
    corners[upperLeftCorner].x;
  v1.y = corners[(upperLeftCorner + 1) % 4].y -
    corners[upperLeftCorner].y;
  printf("v1 = %f %f\n", v1.x, v1.y);
  v2.x = corners[(upperLeftCorner + 2) % 4].x -
    corners[upperLeftCorner].x;
  v2.y = corners[(upperLeftCorner + 2) % 4].y -
    corners[upperLeftCorner].y;
  printf("v2 = %f %f\n", v2.x, v2.y);

  cost = cos(theta);
  sint = sin(theta);
  a[0] = cost;
  a[1] = -sint;
  a[2] = a[0] * (-corners[upperLeftCorner].x) +
    a[1] * (-corners[upperLeftCorner].y);
  a[3] = sint;
  a[4] = cost;
  a[5] = a[3] * (-corners[upperLeftCorner].x) +
    a[4] * (-corners[upperLeftCorner].y);
  
  printf("Xform[0-2]: %f %f %f\n", a[0], a[1], a[2]);
  printf("Xform[3-5]: %f %f %f\n", a[3], a[4], a[5]);

  for (i = 0; i < 4; ++i)
    printf("xcorner %d: %f %f\n", i,
	   a[0] * corners[i].x + a[1] * corners[i].y + a[2],
	   a[3] * corners[i].x + a[4] * corners[i].y + a[5]);
  exit(0);
#endif

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
  
#if 1
  // TEMPORARY FOR DEBUGGING
  struct rlimit rlim;
  rlim.rlim_cur = 10000000000;
  rlim.rlim_max = 10000000000;
  setrlimit(RLIMIT_CORE, &rlim);
#endif

#if 0
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
#endif

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
      absoluteMapsName[0] = '\0';
      modelsName[0] = '\0';
      masksName[0] = '\0';
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
	else if (strcmp(argv[i], "-absolute_maps") == 0)
	  {
	    if (++i == argc)
	      {
		error = 1;
		break;
	      }
	    strcpy(absoluteMapsName, argv[i]);
	  }
	else if (strcmp(argv[i], "-masks") == 0)
	  {
	    if (++i == argc)
	      {
		error = 1;
		break;
	      }
	    strcpy(masksName, argv[i]);
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
	  {
	    if (++i == argc ||
		sscanf(argv[i], "%d", &foldRecovery) != 1)
	      {
		error = 1;
		break;
	      }
	  }
	else if (strcmp(argv[i], "-fold_radius") == 0)
	  {
	    if (++i == argc ||
		sscanf(argv[i], "%f", &foldRadius) != 1)
	      {
		error = 1;
		break;
	      }
	  }
	else if (strcmp(argv[i], "-minimize_area") == 0)
	  minimizeArea = 1;
	else if (strcmp(argv[i], "-output_fold_maps") == 0)
	  outputFoldMaps = 1;
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
	  fprintf(stderr, "              [-fold_recovery count]\n");
	  fprintf(stderr, "              [-output_fold_maps]\n");
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
           level kIntra kInter [kAbsolute [min_iter [threshold [dampingFactor [clampX clampY]]]]]

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
	  minIter = 2000;
	  threshold = 1.0;
	  dampingFactor = 0.0;
	  clampX = 0;
	  clampY = 0;
	  nItems = sscanf(line, "%d%f%f%f%d%f%f%d%d",
			  &level, &kIntra, &kInter, &kAbsolute,
			  &minIter, &threshold, &dampingFactor,
			  &clampX, &clampY);
	  if (nItems < 3)
	    Error("Malformed line in %s:\n%s\n", schedule, line);
	  if (nSteps > 0 && level > steps[nSteps-1].level)
	    Error("Levels must monotonically decrease in schedule.\n");
	  steps = (Step *) realloc(steps, (nSteps+1) * sizeof(Step));
	  steps[nSteps].level = level;
	  steps[nSteps].kIntra = kIntra;
	  steps[nSteps].kInter = kInter;
	  steps[nSteps].kAbsolute = kAbsolute;
	  steps[nSteps].minIter = minIter;
	  steps[nSteps].threshold = threshold;
	  steps[nSteps].dampingFactor = dampingFactor;
	  steps[nSteps].clampX = clampX;
	  steps[nSteps].clampY = clampY;
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
      MPI_Bcast(absoluteMapsName, PATH_MAX, MPI_CHAR, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(modelsName, PATH_MAX, MPI_CHAR, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(masksName, PATH_MAX, MPI_CHAR, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
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
      MPI_Bcast(&nSteps, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&foldRecovery, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&foldRadius, 1, MPI_FLOAT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&minimizeArea, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&outputFoldMaps, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
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

  if (foldRecovery > 0)
    {
      exfoldImage = (int *) malloc(foldRecovery * sizeof(int));
      exfoldWidth = (int *) malloc(foldRecovery * sizeof(int));
      exfoldHeight = (int *) malloc(foldRecovery * sizeof(int));
      exfoldX = (int *) malloc(foldRecovery * sizeof(int));
      exfoldY = (int *) malloc(foldRecovery * sizeof(int));
    }
  
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
	  modelName[0] = '\0';
	  kFactor = 1.0;
	  nItems = sscanf(line, "%s%d%d%lf%lf%lf%lf%s%lf",
			  imageName, &width, &height,
			  &rotation, &scale, &tx, &ty,
			  modelName, &kFactor);
	  if (nItems != 7 && nItems != 8 && nItems != 9)
	    Error("Malformed line in %s:\n%s\n", imageListName, line);
	  while (imageParamsPos + 7 >= imageParamsSize)
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
	  imageParams[imageParamsPos+6] = kFactor;
	  imageParamsPos += 7;

	  imageNameLen = strlen(imageName);
	  while (imageNamesPos + imageNameLen + 1 >= imageNamesSize)
	    {
	      imageNamesSize = (imageNamesSize > 0) ? imageNamesSize * 2 : 1024;
	      imageNames = (char *) realloc(imageNames, imageNamesSize);
	    }
	  strcpy(&imageNames[imageNamesPos], imageName);
	  imageNamesPos += imageNameLen + 1;

	  if (strcmp(modelName, "-") == 0)
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
      images[i].width = imageParams[7*i];
      images[i].height = imageParams[7*i+1];
      images[i].rotation = imageParams[7*i+2];
      images[i].scale = imageParams[7*i+3];
      images[i].tx = imageParams[7*i+4];
      images[i].ty = imageParams[7*i+5];
      images[i].kFactor = imageParams[7*i+6];
      images[i].fixed = 0;
      images[i].owner = -1;
      images[i].needed = 0;
      images[i].sendTo = NULL;
      images[i].nx = (images[i].width + startFactor - 1) / startFactor + 1;
      images[i].ny = (images[i].height + startFactor - 1) / startFactor + 1;
      images[i].nodes = 0;
      images[i].initialNodes = 0;
      images[i].absolutePositions = 0;
      images[i].mask = 0;
      images[i].marked = 0;
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

  Log("Going to initialize the commPhases array\n");
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
  Log("Going to pull out relevant maps\n");
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

  /* construct the list of images to be sent and received
     at each communication phase */
  Log("Constructing list of images to be sent and received\n");
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
  Log("Eliminating unnecessary phases\n");
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
  Log("Initializing models\n");
  mapHashTable = (IntraImageMap **) malloc(nImages * sizeof(IntraImageMap*));
  for (i = 0; i < nImages; ++i)
    mapHashTable[i] = NULL;
  modelMap = NULL;
  for (i = 0; i < nImages; ++i)
    {
      if (images[i].owner != p)
	continue;
      Log("Working on model for image %s (name %s)\n", images[i].name,
	  images[i].modelName);
      mx = (images[i].width + endFactor - 1) / endFactor + 1;
      my = (images[i].height + endFactor - 1) / endFactor + 1;
      hv = HashMap(images[i].modelName, mx, my) % nImages;
      Log("mx = %d my = %d nImages = %d hv = %d\n",
	  mx, my, nImages, hv);
      for (iim = mapHashTable[hv]; iim != NULL; iim = iim->next)
	if (strcmp(iim->name, images[i].modelName) == 0 &&
	    iim->nx == mx && iim->ny == my)
	  break;
      if (iim == NULL)
	{
	  /* create the initial map */
	  Log("Creating iim for image %s\n", images[i].name);
	  iim = (IntraImageMap*) malloc(sizeof(IntraImageMap));
	  iim->name = (char *) malloc(strlen(images[i].modelName) + 1);
	  strcpy(iim->name, images[i].modelName);
	  iim->nx = mx;
	  iim->ny = my;
	  iim->next = mapHashTable[hv];
	  mapHashTable[hv] = iim;
	  if (iim->name[0] != '\0')
	    {
	      /* read the map from a file */
	      sprintf(fn, "%s%s.map", modelsName, iim->name);
	      Log("Going to read map %s\n", fn);
	      if (!ReadMap(fn, &modelMap, &mLevel,
			   &mw, &mh, &mxMin, &myMin,
			   imName0, imName1,
			   msg))
		Error("Could not read map %s:\n  error: %s\n",
		      fn, msg);
	      Log("Finished reading map %s\n", fn);
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
	      mpts = (Point*) malloc(nx * ny * sizeof(Point));
	      mptsc = (float*) malloc(nx * ny * sizeof(float));
	      pnSprings = &(iim->nSprings[startLevel - level]);
	      piSprings = &(iim->springs[startLevel - level]);
	      if (iim->name[0] == '\0')
		/* use the identity map */
		for (y = 0; y < ny; ++y)
		  for (x = 0; x < nx; ++x)
		    {
		      mpts[y*nx+x].x = x;
		      mpts[y*nx+x].y = y;
		      mptsc[y*nx+x] = 1.0;
		    }
	      else
		{
		  scale = ((double) factor) / (1 << mLevel);
		  threshold = scale;
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
			ixv -= mxMin;
			iyv -= myMin;
			if (ixv >= 0 && ixv < mw-1 &&
			    iyv >= 0 && iyv < mh-1 &&
			    modelMap[iyv*mw+ixv].c > 0.0 &&
			    modelMap[(iyv+1)*mw+ixv].c > 0.0 &&
			    modelMap[iyv*mw+ixv+1].c > 0.0 &&
			    modelMap[(iyv+1)*mw+ixv+1].c > 0.0)
			  {
			    rx00 = modelMap[iyv*mw+ixv].x;
			    ry00 = modelMap[iyv*mw+ixv].y;
			    rc00 = modelMap[iyv*mw+ixv].c;
			    rx01 = modelMap[(iyv+1)*mw+ixv].x;
			    ry01 = modelMap[(iyv+1)*mw+ixv].y;
			    rc01 = modelMap[(iyv+1)*mw+ixv].c;
			    rx10 = modelMap[iyv*mw+ixv+1].x;
			    ry10 = modelMap[iyv*mw+ixv+1].y;
			    rc10 = modelMap[iyv*mw+ixv+1].c;
			    rx11 = modelMap[(iyv+1)*mw+ixv+1].x;
			    ry11 = modelMap[(iyv+1)*mw+ixv+1].y;
			    rc11 = modelMap[(iyv+1)*mw+ixv+1].c;
			    rx = rx00 * (rrx - 1.0) * (rry - 1.0)
			      - rx10 * rrx * (rry - 1.0) 
			      - rx01 * (rrx - 1.0) * rry
			      + rx11 * rrx * rry;
			    ry = ry00 * (rrx - 1.0) * (rry - 1.0)
			      - ry10 * rrx * (rry - 1.0) 
			      - ry01 * (rrx - 1.0) * rry
			      + ry11 * rrx * rry;
			    rc = rc00;
			    if (rc01 < rc)
			      rc = rc01;
			    if (rc10 < rc)
			      rc = rc10;
			    if (rc11 < rc)
			      rc = rc11;
			    if (rc < 0.0)
			      Error("rc is negative! %f\n", rc);
			  }
			else
			  {
			    if (!Extrapolate(&rx, &ry, &rc,
					     ixv, iyv, rrx, rry,
					     modelMap, mw, mh, threshold))
			      {
				rx = 0.0;
				ry = 0.0;
				rc = 0.0;
			      }
			  }
			mpts[y*nx+x].x = rx / scale;
			mpts[y*nx+x].y = ry / scale;
			mptsc[y*nx+x] = rc;
#if DEBUG
#if PDEBUG
			if (i == ioi && x == poix && y == poiy)
#endif
			  Log("iim pts (%lu) %d %d = %f %f  ixv = %d iyv = %d rrx = %f rry = %f\n",
			      mpts, x, y,
			      mpts[y*nx+x].x,
			      mpts[y*nx+x].y,
			      ixv, iyv, rrx, rry);
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
		  iis->k = mptsc[iis->index0] < mptsc[iis->index1] ?
		    mptsc[iis->index0] : mptsc[iis->index1];
		  ++iis;
		}
	      if (level == startLevel)
		iim->points = mpts;
	      else
		free(mpts);
	      free(mptsc);
	    }
	  if (modelMap != NULL)
	    {
	      free(modelMap);
	      modelMap = NULL;
	    }
	}
      images[i].map = iim;
    }

  /* initialize the image masks */
  Log("Initializing image masks\n");
  eFactor = 1 << endLevel;
  rFactor = (1 << startLevel) / eFactor;
  Log("Going to initialize the image masks (eFactor = %d rFactor = %d)\n",
      eFactor, rFactor);
  for (i = 0; i < nImages; ++i)
    {
      if (images[i].owner != p && !images[i].needed)
	continue;
      if (masksName[0] != '\0')
	{
	  /* read the mask from a file */
	  sprintf(fn, "%s%s.pbm", masksName, images[i].name);
	  Log("Reading mask %s\n", fn);
	  mask = NULL;
	  if (ReadBitmap(fn, &mask, &maskWidth, &maskHeight,
			 -1, -1, -1, -1,
			 msg))
	    {
	      Log("Read mask %s  -- width=%d height=%d\n",
		  fn, maskWidth, maskHeight);
	      if (maskWidth != images[i].width ||
		  maskHeight != images[i].height)
		Error("Mask %s is not same size as image %s\n",
		      fn, images[i].name);
	    }
	  else
	    {
	      Log("WARNING: could not read mask %s -- assuming unmasked image.\n",
		  fn);
	      maskWidth = images[i].width;
	      maskHeight = images[i].height;
	      mbpl = (maskWidth + 7) >> 3;
	      mask = (unsigned char *) malloc(maskHeight * mbpl);
	      memset(mask, 0xff, maskHeight * mbpl);
	    }
	  mbpl = (maskWidth + 7) >> 3;
	  rnx = (maskWidth + eFactor - 1) / eFactor;
	  rny = (maskHeight + eFactor - 1) / eFactor;
	  imbpl = (rnx + 7) >> 3;
	  imask = (unsigned char *) malloc(rny * imbpl *
					   sizeof(unsigned char));
	  memset(imask, 0, rny * imbpl * sizeof(unsigned char));
	  count = 0;
	  for (y = 0; y < rny; ++y)
	    for (x = 0; x < rnx; ++x)
	      {
		valid = 0;
		for (iy = 0; !valid && iy < eFactor; ++iy)
		  {
		    iyv = y * eFactor + iy;
		    if (iyv >= maskHeight)
		      break;
		    for (ix = 0; ix < eFactor; ++ix)
		      {
			ixv = x * eFactor + ix;
			if (ixv >= maskWidth)
			  break;
			if (mask[iyv * mbpl + (ixv >> 3)] & (0x80 >> (ixv & 7)))
			  {
			    valid = 1;
			    break;
			  }
		      }
		  }
		if (valid)
		  {
		    imask[y*imbpl + (x >> 3)] |= 0x80 >> (x & 7);
		    ++count;
		  }
	      }
	  Log("Initialized mask for %s with %d of %d (%d x %d) bits set.\n",
	      images[i].name, count, rnx*rny, rnx, rny);
	  images[i].mask = imask;
	  free(mask);
	}
      else
	images[i].mask = NULL;
    }

  Log("Before node init barrier\n");
  if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("MPI_Barrier after cpi failed.\n");
  Log("After node init barrier\n");

  /* initialize the image nodes */
  initialMap = NULL;
  for (i = 0; i < nImages; ++i)
    {
      if (images[i].owner != p && !images[i].needed)
	continue;
      nx = images[i].nx;
      ny = images[i].ny;
      images[i].nodes = (Node*) malloc(nx * ny * sizeof(Node));
      node = images[i].nodes;
      if (images[i].owner != p)
	{
	  memset(node, 0, nx * ny * sizeof(Node));
	  continue;
	}

      // reset the nodes
      for (y = 0; y < ny; ++y)
	for (x = 0; x < nx; ++x)
	  {
	    node[y * nx + x].x = UNSPECIFIED;
	    node[y * nx + x].y = UNSPECIFIED;
	    node[y * nx + x].fx = 0.0;
	    node[y * nx + x].fy = 0.0;
	  }

      // mark the valid map nodes
      Log("Marking the valid map nodes for image %s\n", images[i].name);
      rnx = (images[i].width + eFactor - 1) / eFactor;
      rny = (images[i].height + eFactor - 1) / eFactor;
      imbpl = (rnx + 7) >> 3;
      imask = images[i].mask;
      for (y = 0; y < ny; ++y)
	for (x = 0; x < nx; ++x)
	  {
	    valid = (imask == NULL);
	    for (iy = 0; !valid && iy < rFactor; ++iy)
	      {
		iyv = y * rFactor + iy;
		if (iyv >= rny)
		  break;
		for (ix = 0; ix < rFactor; ++ix)
		  {
		    ixv = x * rFactor + ix;
		    if (ixv >= rnx)
		      break;
		    if (imask[iyv * imbpl + (ixv >> 3)] & (0x80 >> (ixv & 7)))
		      {
			valid = 1;
			break;
		      }
		  }
	      }
	    if (valid)
	      {
		node[y * nx + x].x = 0.0;
		node[y * nx + x].y = 0.0;
		if (x < nx - 1)
		  {
		    node[y * nx + x + 1].x = 0.0;
		    node[y * nx + x + 1].y = 0.0;
		  }
		if (y < ny - 1)
		  {
		    node[(y + 1) * nx + x].x = 0.0;
		    node[(y + 1) * nx + x].y = 0.0;
		  }
		if (x < nx - 1 && y < ny - 1)
		  {
		    node[(y + 1) * nx + x + 1].x = 0.0;
		    node[(y + 1) * nx + x + 1].y = 0.0;
		  }
	      }
	  }
      Log("Finished marking the valid map nodes for image %s\n", images[i].name);

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
	    for (x = 0; x < nx; ++x, ++node)
	      {
		if (node->x > 0.5 * UNSPECIFIED)
		  continue;
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
		node->x = rx / scale;
		node->y = ry / scale;
		node->fx = 0.0;
		node->fy = 0.0;
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
	    for (x = 0; x < nx; ++x, ++node)
	      {
		if (node->x > 0.5 * UNSPECIFIED)
		  continue;
		node->x =  cost * mpts[y*nx+x].x + sint * mpts[y*nx+x].y
		  + images[i].tx / startFactor;
		node->y =  -sint * mpts[y*nx+x].x + cost * mpts[y*nx+x].y
		  + images[i].ty / startFactor;
		node->fx = 0.0;
		node->fy = 0.0;
#if DEBUG
#if PDEBUG
		if (i == ioi && x == poix && y == poiy)
#endif
		  Log("Initializing position of %d(%d,%d) to (%f,%f)\n",
		      i, x, y, node->x, node->y);
#endif
	      }
	}
    }
  if (initialMap != NULL)
    {
      free(initialMap);
      initialMap = NULL;
    }

  Log("Initializing maps\n");
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
	 let the threshold be a factor from sqrt(2)/2 to sqrt(2)
           of the distance from grid point to grid point in terms
	   of the map spacing; we let this be sqrt(2) at the
           startLevel to always have at least 4 springs; and reduce
           the value to sqrt(2)/2 at the endLevel)
	 at each level determine if the distance is less than
	   or equal to the threshold; if so, make a spring
      */
      factor = 1 << startLevel;
      dnx = (images[m->image0].width + factor - 1) / factor;
      dny = (images[m->image0].height + factor - 1) / factor;
      dnx = dnx * factor / mFactor + 1;
      dny = dny * factor / mFactor + 1;
      mbpl = (dnx + 7) >> 3;
      mask = (unsigned char *) malloc(dny * mbpl);
      memset(mask, 0, dny * mbpl);
      count1  = 0;
      count2 = 0;
      for (y = 0; y < mh; ++y)
	for (x = 0; x < mw; ++x)
	  {
	    if (strcmp(m->name, "z3803_z3807") == 0)
	      Log("map[%d*mw+%d].x = %f .y = %f .c = %f\n",
		  y, x, map[y*mw+x].x, map[y*mw+x].y, map[y*mw+x].c);
	    if (map[y*mw+x].c > 0.0)
	    {
	      ixv = x + mxMin;
 	      iyv = y + myMin;
	      ++count2;
	      if (ixv >= 0 && ixv < dnx &&
		  iyv >= 0 && iyv < dny)
		{
		  mask[iyv*mbpl+(ixv >> 3)] |= 0x80 >> (ixv & 7);
		  ++count1;
		}
	    }
	  }
      Log("For map %s count1 = %d count2 = %d mw=%d mh=%d mw*mh = %d dnx*dny = %d  dnx=%d dny=%d mxMin=%d myMin=%d\n",
	  m->name, count1, count2, mw, mh, mw * mh, dnx * dny,
	  dnx, dny, mxMin, myMin);
      dist = (float *) malloc(dny * dnx * sizeof(float));
      computeDistance(CHESSBOARD_DISTANCE, dnx, dny, mask, dist);
#if 0
      if (strcmp(m->name, "z3803_z3807") == 0)
	{
	  for (iyv = 0; iyv < dny; ++iyv)
	    for (ixv = 0; ixv < dnx; ++ixv)
	      Log("dist[%d,%d] = %f\n",
		  ixv, iyv, dist[iyv*dnx+ixv]);
	}
#endif

      imask = images[m->image0].mask;
      rnx = (images[m->image0].width + eFactor - 1) / eFactor;
      rny = (images[m->image0].height + eFactor - 1) / eFactor;
      imbpl = (rnx + 7) >> 3;
      for (level = startLevel; level >= endLevel; --level)
	{
	  factor = 1 << level;
	  threshold = ((float) factor) / mFactor;

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
	  count1 = 0;
	  count2 = 0;
	  for (y = 0; y < ny; ++y)
	    {
	      firstInStrip = 1;
	      iyMin = (y - 1) * factor / eFactor;
	      if (iyMin < 0)
		iyMin = 0;
	      iyMax = (y + 1) * factor / eFactor - 1;
	      if (iyMax >= rny)
		iyMax = rny - 1;
	      yv = ((float) (y * factor)) / mFactor;
	      for (x = 0; x < nx; ++x)
		{
		  // check if the mask supports having a spring
		  //    at this node
		  
		  ixMin = (x - 1) * factor / eFactor;
		  if (ixMin < 0)
		    ixMin = 0;
		  ixMax = (x + 1) * factor / eFactor - 1;
		  if (ixMax >= rnx)
		    ixMax = rnx - 1;
		  valid = (imask == NULL);
		  for (iyv = iyMin; iyv <= iyMax && !valid; ++iyv)
		    for (ixv = ixMin; ixv <= ixMax && !valid; ++ixv)
		      if (imask[iyv*imbpl + (ixv >> 3)] & (0x80 >> (ixv & 7)))
			valid = 1;
#if DEBUG
#if PDEBUG
		  if (m->image0 == ioi && x == poix && y == poiy)
#endif
		    Log("cons point at level %d x = %d y = %d ixMin = %d ixMax = %d iyMin = %d iyMax = %d  valid = %d  imask = %p rnx = %d rny = %d\n",
			level, x, y,
			ixMin, ixMax, iyMin, iyMax,
			valid, imask, rnx, rny);
#endif		    
		  if (!valid)
		    {
		      firstInStrip = 1;
		      continue;
		    }
		  ++count1;

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
#if 0
		      if (strcmp(m->name, "z3803_z3807") == 0)
			Log("dist[%d*dnx+%d]=%f > %f\n",
			    iyv, ixv, dist[iyv*dnx+ixv], threshold);
#endif
		      firstInStrip = 1;
		      continue;
		    }
		  ++count2;
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
#if DEBUG
#if PDEBUG
		      if (m->image0 == ioi && x == poix && y == poiy)
#endif
			Log("RXY: ixv = %d iyv = %d rx00=%f ry00=%f rx01=%f ry01=%f rx10=%f ry10=%f rx11=%f ry11=%f rrx=%f rry=%f rx=%f ry=%f\n",
			    ixv, iyv, rx00, ry00, rx01, ry01, rx10, ry10, rx11, ry11, rrx, rry, rx, ry);
#endif
		    }
		  else
		    {
#if DEBUG
#if PDEBUG
		      if (m->image0 == ioi && x == poix && y == poiy)
#endif
		      
			Log("Attempting to extrapolate\n");
#endif
		      if (!Extrapolate(&rx, &ry, &rc,
				       ixv, iyv, rrx, rry,
				       map, mw, mh, threshold))
			{
			  firstInStrip = 1;
			  continue;
			}
#if DEBUG
#if PDEBUG
		      if (m->image0 == ioi && x == poix && y == poiy)
#endif
			Log("Extrapolated: ixv=%d iyv=%d rrx=%f rry=%f mw=%d mh=%d rx=%f ry=%f\n",
			    ixv, iyv, rrx, rry, mw, mh, rx, ry);
#endif
		    }
		  rx = rx * mFactor / factor;
		  ry = ry * mFactor / factor;
		  irx = floor(rx + 0.5);
		  if (irx < 0)
		    irx = 0;
		  if (irx >= nx1)
		    irx = nx1 - 1;
		  iry = floor(ry + 0.5);
		  if (iry < 0)
		    iry = 0;
		  if (iry >= ny1)
		    iry = ny1 - 1; 
#if DEBUG
#if PDEBUG
		  if (m->image0 == ioi && x == poix && y == poiy ||
		      m->image1 == ioi && irx == poix && iry == poiy)
#endif
		    Log("spring at lvl = %d from im %d x = %d y = %d xv = %f yv = %f ixv = %d iyv = %d rrx = %f rry = %f to im %d rx = %f ry = %f irx = %d iry = %d\n",
			level, m->image0, x, y, xv, yv, ixv, iyv, rrx, rry,
			m->image1, rx, ry, irx, iry);
#endif

		  rrx = rx - irx;
		  rry = ry - iry;
		  if (rrx < -1.6383 || rrx > 1.6383 ||
		      rry < -1.6383 || rry > 1.6383)
		    {
		      firstInStrip = 1;
		      continue;
		    }
		  
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
		  irrx = floor(20000.0 * rrx + 0.5);
		  irry = floor(20000.0 * rry + 0.5);
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
	  Log("For map %s at level %d, count1 = %d  count2 = %d  total = %d\n",
	      m->name, level, count1, count2, nx * ny);
	  if (*pnStrips != 0)
	    *pStrips = (InterImageStrip*)
	      realloc(*pStrips, *pnStrips * sizeof(InterImageStrip));
	  if (*pnSprings != 0)
	    *pSprings = (InterImageSpring*)
	      realloc(*pSprings, *pnSprings * sizeof(InterImageSpring));
#if 0
	  else
	    Error("No springs were generated from map %s at level %d\n",
		  m->name, level);
#endif
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
	      nodes = images[i].nodes;
	      j = 0;
	      while (fscanf(f, "%d %d %f %f", &cx, &cy, &consX, &consY) == 4)
		{
		  if (cx < 0 || cx >= nx || cy < 0 || cy >= ny)
		    continue;
		  node[cy*nx + cx].x = consX;
		  node[cy*nx + cx].fx = (j + 1) * CONSTRAINED;
		  node[cy*nx + cx].y = consY;
		  node[cy*nx + cx].fy = (j + 1) * CONSTRAINED;
		  ++j;
		}
	      fclose(f);
	    }
	}
    }
#endif

  /* compute absolute positions if needed */
  for (step = 0; step < nSteps; ++step)
    {
      level = steps[step].level;
      if (level > startLevel || level < endLevel)
	continue;
      kAbsolute = steps[step].kAbsolute;
      if (kAbsolute == 0.0)
	continue;
      factor = 1 << level;
      for (i = myFirstImage; i <= myLastImage; ++i)
	{
	  if (images[i].absolutePositions == NULL)
	    {
	      images[i].absolutePositions = (Point **) malloc(nLevels *
							      sizeof(Point*));
	      memset(images[i].absolutePositions, 0, nLevels * sizeof(Point*));
	    }
	  if (images[i].absolutePositions[startLevel - level] == NULL)
	    {
	      nx = (images[i].width + factor - 1) / factor + 1;
	      ny = (images[i].height + factor - 1) / factor + 1;
	      //	      Log("Allocating absPos[%d][%d] with size %d %d\n", i, level, nx, ny);
	      images[i].absolutePositions[startLevel - level] = (Point *)
		malloc(nx * ny * sizeof(Point));
	    }
	}
    }
  for (i = myFirstImage; i <= myLastImage; ++i)
    {
      if (images[i].absolutePositions == NULL)
	continue;

      /* read in absolute positions map */
      map = NULL;
      sprintf(fn, "%s%s.map", absoluteMapsName, images[i].name);
      if (!ReadMap(fn, &map, &mLevel,
		   &mw, &mh, &mxMin, &myMin,
		   imName0, imName1,
		   msg))
	{
	  Log("WARNING: Could not open %s so will not use absolute map.\n",
	      fn);
	  free(images[i].absolutePositions);
	  images[i].absolutePositions = NULL;
	  continue;
	}
      mFactor = 1 << mLevel;

      /* create an array with coverage corresponding to the
	 startLevel, and resolution determined by the mLevel
	 fill in the entries with a 1 if the map has nonzero c
	 fill in the other entries with 0
	 compute distance
	 let the threshold be a factor from sqrt(2)/2 to sqrt(2)
           of the distance from grid point to grid point in terms
	   of the map spacing; we let this be sqrt(2) at the
           startLevel to always have at least 4 springs; and reduce
           the value to sqrt(2)/2 at the endLevel)
	 at each level determine if the distance is less than
	   or equal to the threshold; if so, make a spring
      */
      factor = 1 << startLevel;
      dnx = (images[i].width + factor - 1) / factor;
      dny = (images[i].height + factor - 1) / factor;
      dnx = dnx * factor / mFactor + 1;
      dny = dny * factor / mFactor + 1;
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
      computeDistance(CHESSBOARD_DISTANCE, dnx, dny, mask, dist);

      for (level = startLevel; level >= endLevel; --level)
	{
	  factor = 1 << level;

	  if (startLevel == endLevel)
	    threshold = ((float) factor) / (1 << mLevel);
	  else
	    threshold = ((float) factor)
	      * exp(M_LN2 * (level - endLevel) / (startLevel - endLevel))
	      / 2.0
	      / (1 << mLevel);
	  nx = (images[i].width + factor - 1) / factor + 1;
	  ny = (images[i].height + factor - 1) / factor + 1;
	  //	  Log("Initializing absPos[%d][%d] with size %d %d\n", i, level, nx, ny);
	  absPos = images[i].absolutePositions[startLevel - level];
	  for (y = 0; y < ny; ++y)
	    {
	      yv = ((float) (y * factor)) / mFactor;
	      for (x = 0; x < nx; ++x)
		{
		  absPos->x = UNSPECIFIED;
		  absPos->y = UNSPECIFIED;
		  xv = ((float) (x * factor)) / mFactor;
		  ixv = (int) floor(xv);
		  iyv = (int) floor(yv);
#if DEBUG
#if PDEBUG
		  if (i == ioi && x == poix && y == poiy)
#endif
		    Log("cons abs spring at level %d x = %d y = %d xv = %f yv = %f ixv = %d iyv = %d dnx = %d dny = %d dist = %f thresh = %f\n",
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
		      ++absPos;
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
#if DEBUG
		      Log("RXY: ixv = %d iyv = %d rx00=%f ry00=%f rx01=%f ry01=%f rx10=%f ry10=%f rx11=%f ry11=%f rrx=%f rry=%f rx=%f ry=%f\n",
			  ixv, iyv, rx00, ry00, rx01, ry01, rx10, ry10, rx11, ry11, rrx, rry, rx, ry);
#endif
		    }
		  else
		    {
#if DEBUG
		      Log("Attempting to extrapolate\n");
#endif
		      if (!Extrapolate(&rx, &ry, &rc,
				       ixv, iyv, rrx, rry,
				       map, mw, mh, threshold))
			{
			  ++absPos;
			  continue;
			}
#if DEBUG
		      Log("Extrapolated: ixv=%d iyv=%d rrx=%f rry=%f mw=%d mh=%d rx=%f ry=%f\n",
			  ixv, iyv, rrx, rry, mw, mh, rx, ry);
#endif
		    }
		  absPos->x = rx * mFactor / factor;
		  absPos->y = ry * mFactor / factor;
		  ++absPos;
		}
	    }
	}
      free(dist);
      free(mask);
      free(map);
      map = NULL;
    }

  prevLevel = -1;
  for (step = 0; step < nSteps; ++step)
    {
      level = steps[step].level;
      if (p == 0)
	Log("Starting iterations for step %d (level %d)\n", step, level);

      if (level != prevLevel)
	{
	  if (prevLevel >= 0)
	    RefinePositions(prevLevel, level);
	  PlanCommunications(level);
	  prevLevel = level;
	}

      if (foldRecovery > 0)
	{
	  // make a backup copy of the point positions in case
	  //   we have to restart this step
	  for (i = myFirstImage; i <= myLastImage; ++i)
	    {
	      nx = images[i].nx;
	      ny = images[i].ny;
	      images[i].initialNodes = (Node *) realloc(images[i].initialNodes,
							  nx * ny * sizeof(Node));
	      memcpy(images[i].initialNodes, images[i].nodes,
		     nx * ny * sizeof(Node));
	    }
	}

      kIntra = steps[step].kIntra;
      kInter = steps[step].kInter;
      kAbsolute = steps[step].kAbsolute;
      threshold = steps[step].threshold;
      factor = 1 << level;

    restartStep:
      if (steps[step].dampingFactor > 0.0)
	{
	  fixedDamping = steps[step].dampingFactor;
	  dampingFactor = fixedDamping;
	}
      else
	{
	  fixedDamping = 0.0;
	  dampingFactor = 0.001;
	}

      outputRequestedIter = -1;
      terminationRequestedIter = -1;
      refinementRequestedIter = -1;

      epochInitialTotalEnergy = 1.0e+30;
      epochFinalTotalEnergy = 1.0e+30;
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
	      imageWithFold = CheckForFolds(level, iter);
	      if (MPI_Allreduce(&imageWithFold, &minImageWithFold, 1, MPI_INT,
				MPI_MIN, MPI_COMM_WORLD) != MPI_SUCCESS)
		Error("MPI_Allreduce: could not check for folds.\n");
	      if (minImageWithFold < nImages)
		{
		  foldDetected = 1;

		  if (foldRecovery > 0)
		    {
		      if (foldRecoveryCount >= foldRecovery)
			Error("Maximum number of fold recoveries exceeded.\n");
		      RecoverFromFold(level);
		      ++foldRecoveryCount;
		      goto restartStep;
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
	      if (images[i].absolutePositions != NULL)
		{
		  /* compute absolute location forces */
		  node = images[i].nodes;
		  absPos = images[i].absolutePositions[startLevel - level];
		  for (y = 0; y < ny; ++y)
		    for (x = 0; x < nx; ++x, ++node, ++absPos)
		      {
			if (node->x > 0.5 * UNSPECIFIED)
			  continue;
			if (absPos->x < 0.5 * UNSPECIFIED)
			  deltaX = absPos->x - node->x;
			else
			  deltaX = 0.0;
			if (absPos->y < 0.5 * UNSPECIFIED)
			  deltaY = absPos->y - node->y;
			else
			  deltaY = 0;
			if (node->fx < 0.5 * CONSTRAINED)
			  node->fx = kAbsolute * deltaX;
			if (node->fy < 0.5 * CONSTRAINED)
			  node->fy = kAbsolute * deltaY;
			energy += kAbsolute * (deltaX * deltaX + deltaY * deltaY);
			if (isinf(energy) || isnan(energy))
			  abort();
			//			if (x == 0 && y == 0)
			//			  Log("ENERGY COMP %f %f %f %f %f %f %f %f\n",
			//			      kAbsolute, deltaX, deltaY, absPos->x, absPos->y,
			//			      node->x, node->y, energy);
		      }
		}
	      else
		{
		  /* zero all forces */
		  node = images[i].nodes;
		  for (y = 0; y < ny; ++y)
		    for (x = 0; x < nx; ++x, ++node)
		      {
			if (node->fx < 0.5 * CONSTRAINED)
			  node->fx = 0.0;
			if (node->fy < 0.5 * CONSTRAINED)
			  node->fy = 0.0;
		      }
		}

	      /* add in intra-section forces */
	      kIntraThisImage = kIntra * images[i].kFactor;
	      iim = images[i].map;
	      nodes = images[i].nodes;
	      nSprings = iim->nSprings[startLevel - level];
	      iSprings = iim->springs[startLevel - level];
	      for (k = 0; k < nSprings; ++k)
		{
		  iis = &(iSprings[k]);
		  if (nodes[iis->index0].x > 0.5 * UNSPECIFIED ||
		      nodes[iis->index1].x > 0.5 * UNSPECIFIED)
		    continue;
		  deltaX = nodes[iis->index1].x - nodes[iis->index0].x;
		  deltaY = nodes[iis->index1].y - nodes[iis->index0].y;
		  d = sqrt(deltaX * deltaX + deltaY * deltaY);
		  sk = kIntraThisImage * iis->k;
		  force = sk * (d - iis->nomD);
		  energy += force * (d - iis->nomD);
		  if (isinf(energy) || isnan(energy))
		    abort();
		  if (d != 0.0)
		    {
		      forceOverD = force / d;
		      dfx = forceOverD * deltaX;
		      dfy = forceOverD * deltaY;
		      nodes[iis->index0].fx += dfx;
		      nodes[iis->index0].fy += dfy;
		      nodes[iis->index1].fx -= dfx;
		      nodes[iis->index1].fy -= dfy;
		    }
#if DEBUG
#if PDEBUG
		  if (i == ioi && iis->index0 % images[i].nx == poix &&
		      iis->index0 / images[i].nx == poiy)
#endif
		    Log("intraforce: %d(%d,%d) - %d(%d,%d): (%f %f) to (%f %f) dist %f nom %f force (%f %f)\n",
		      i, iis->index0 % images[i].nx, iis->index0 / images[i].nx,
		      i, iis->index1 % images[i].nx, iis->index1 / images[i].nx,
		      nodes[iis->index0].x, nodes[iis->index0].y,
		      nodes[iis->index1].x, nodes[iis->index1].y,
		      d, iis->nomD,
		      d != 0.0 ? dfx : 1000000000.0,
		      d != 0.0 ? dfy : 1000000000.0);
#endif
		}
	    }
	  if (p == 0 && iter % 100 == 0)
	    {
	      Log("intra-energy = %f  (kIntra = %f)\n", energy - basis, kIntra);
	      Log("energy = %f basis = %f\n", energy, basis);
	    }
	  //	  nodes0 = images[0].nodes;
	  //	  nx = images[0].nx;
	  //	  if (iter == 0)
	  //	    Log("PTEST %f %f\n", nodes[2*nx].x, nodes[2*nx].y);

	  /* add in inter-image forces (one spring method) */
	  basis = energy;
	  for (i = 0; i < nMaps; ++i)
	    {
	      m = &maps[i];
	      msk = kInter * m->k;
	      if (msk == 0.0)
		continue;
	      sme = m->energyFactor != 0.0;
	      nStrips = m->nStrips[startLevel - level];
	      strips = m->strips[startLevel - level];
	      springs = m->springs[startLevel - level];
	      nodes0 = images[m->image0].nodes;
	      if (nodes0 == NULL)
		abort();
	      nodes1 = images[m->image1].nodes;
	      nx = images[m->image0].nx;
	      nx1 = images[m->image1].nx;
	      springsPos = 0;
	      for (j = 0; j < nStrips; ++j)
		{
		  if (nodes0 == NULL)
		    abort();
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
		      if (nodes0[y*nx+x].x > 0.5 * UNSPECIFIED ||
			  nodes1[iry*nx1+irx].x > 0.5 * UNSPECIFIED)
			continue;
		      xv = nodes1[iry*nx1+irx].x + 0.00005 * s->dx;
		      yv = nodes1[iry*nx1+irx].y + 0.00005 * s->dy;
		      deltaX = xv - nodes0[y*nx + x].x;
		      deltaY = yv - nodes0[y*nx + x].y;
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
			  nodes0[y*nx+x].x, nodes0[y*nx+x].y,
			  xv, yv,
			  irx, iry,
			  nodes1[iry*nx1+irx].x, nodes1[iry*nx1+irx].y,
			  0.00005 * s->dx, 0.00005 * s->dy,
			  kdx, kdy);
#endif
		      if (nodes0 == NULL)
			abort();
		      nodes0[y*nx + x].fx += kdx;
		      nodes0[y*nx + x].fy += kdy;
		      if (sme)
			energy += sk * (deltaX * deltaX + deltaY * deltaY);
		      if (isinf(energy) || isnan(energy))
			abort();
		      nodes1[iry*nx1 + irx].fx -= kdx;
		      nodes1[iry*nx1 + irx].fy -= kdy;
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
	      nNodes = images[i].nx * images[i].ny;
	      node = images[i].nodes;
	      for (k = 0; k < nNodes; ++k, ++node)
		{
		  if (node->x > 0.5 * UNSPECIFIED)
		    continue;
		  if (node->fx < 0.5 * CONSTRAINED)
		    if (node->fy < 0.5 * CONSTRAINED)
		      force = dampingFactor * hypot(node->fx, node->fy);
		    else
		      force = dampingFactor * fabsf(node->fx);
		  else 
		    if (node->fy < 0.5 * CONSTRAINED)
		      force = dampingFactor * fabsf(node->fy);
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
	  if (steps[step].clampX)
	    maxStepX = 0.0;
	  else
	    maxStepX = 0.1 * factor;
	  if (steps[step].clampY)
	    maxStepY = 0.0;
	  else
	    maxStepY = 0.1 * factor;
	  for (i = myFirstImage; i <= myLastImage; ++i)
	    {
	      if (images[i].fixed)
		continue;
	      nNodes = images[i].nx * images[i].ny;
	      node = images[i].nodes;
	      for (k = 0; k < nNodes; ++k, ++node)
		{
		  if (node->x > 0.5 * UNSPECIFIED)
		    continue;
		  if (node->fx < 0.5 * CONSTRAINED)
		    {
		      deltaX = scale * node->fx;
		      if (deltaX > maxStepX)
			deltaX = maxStepX;
		      else if (deltaX < -maxStepX)
			deltaX = -maxStepX;
		      node->x += deltaX;

		    }
		  if (node->fy < 0.5 * CONSTRAINED)
		    {
		      deltaY = scale * node->fy;
		      if (deltaY > maxStepY)
			deltaY = maxStepY;
		      else if (deltaY < -maxStepY)
			deltaY = -maxStepY;
		      node->y += deltaY;
		    }
#if DEBUG
#if PDEBUG
		  if (i == ioi && k % images[i].nx == poix &&
		      k / images[i].nx == poiy)
#endif
		  Log("moving point %d(%d,%d) at (%f %f) by (%f %f) to (%f %f)\n",
		      i, k % images[i].nx, k / images[i].nx,
		      node->x - deltaX, node->y - deltaY,
		      deltaX, deltaY,
		      node->x, node->y);
#endif
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
		  if (iter >= steps[step].minIter &&
		      (epochInitialTotalEnergy - totalEnergy <
		       totalEnergy * threshold * 0.000001 * epochIterations ||
		       epochFinalTotalEnergy - totalEnergy <
		       totalEnergy * threshold * 0.000001 * epochIterations ||
		       totalEnergy < 0.000001))
		    controlFlag |= 1;
		}
	      epochFinalTotalEnergy = totalEnergy;
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
		  if (dampingFactor < 0.000001)
		    {
		      Log("Damping factor = %f... ending iterations for this step.\n",
			  dampingFactor);
		      goto finalFoldCheck;
		    }
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
	    Log("Refinement (requested at iter %d, nDecrease = %d)\n",
		refinementRequestedIter, nDecrease);
	  else
	    continue;

	finalFoldCheck:
	  imageWithFold = CheckForFolds(level, iter);
	  if (MPI_Allreduce(&imageWithFold, &minImageWithFold, 1, MPI_INT,
			    MPI_MIN, MPI_COMM_WORLD) != MPI_SUCCESS)
	    Error("MPI_Allreduce: could not check for folds.\n");
	  if (minImageWithFold < nImages)
	    {
	      foldDetected = 1;
	      
	      if (foldRecovery > 0)
		{
		  if (foldRecoveryCount >= foldRecovery)
		    Error("Maximum number of fold recoveries exceeded.\n");
		  RecoverFromFold(level);
		  ++foldRecoveryCount;
		  goto restartStep;
		}
	    }
	  break;
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

  if (foldDetected)
    {
      processWithFold = images[minImageWithFold].owner;
      if (MPI_Bcast(&foldImage, 1, MPI_INT,
		    processWithFold, MPI_COMM_WORLD) != MPI_SUCCESS ||
	  MPI_Bcast(&foldWidth, 1, MPI_INT,
		    processWithFold, MPI_COMM_WORLD) != MPI_SUCCESS ||
	  MPI_Bcast(&foldHeight, 1, MPI_INT,
		    processWithFold, MPI_COMM_WORLD) != MPI_SUCCESS ||
	  MPI_Bcast(&foldX, 1, MPI_INT,
		    processWithFold, MPI_COMM_WORLD) != MPI_SUCCESS ||
	  MPI_Bcast(&foldY, 1, MPI_INT,
		    processWithFold, MPI_COMM_WORLD) != MPI_SUCCESS ||
	  MPI_Bcast(&foldAbsoluteX, 1, MPI_FLOAT,
		    processWithFold, MPI_COMM_WORLD) != MPI_SUCCESS ||
	  MPI_Bcast(&foldAbsoluteY, 1, MPI_FLOAT,
		    processWithFold, MPI_COMM_WORLD) != MPI_SUCCESS)
	Error("Could not broadcast fold position.\n");

      if (outputGridName[0] != '\0')
	{
	  Log("OUTPUTTING FOLDED GRID\n");
	  strcpy(outputGridFocusImage, images[foldImage].name);
	  outputGridFocusWidth = foldWidth;
	  outputGridFocusHeight = foldHeight;
	  //	  outputGridFocusDepth = 2;
	  outputGridFocusX = foldX;
	  outputGridFocusY = foldY;
	  outputGridSequence = 0;
	  OutputGrid(level, -2);
	}

      if (outputFoldMaps)
	{
	  Log("OUTPUTTING FOLDED MAPS\n");
	  OutputFoldMaps(level, iter);
	}
    }

  if (foldRecoveryCount > 0)
    {
      for (i = 0; i < foldRecoveryCount; ++i)
	{
	  strcpy(outputGridFocusImage, images[exfoldImage[i]].name);
	  outputGridFocusWidth = exfoldWidth[i];
	  outputGridFocusHeight = exfoldHeight[i];
	  outputGridFocusX = exfoldX[i];
	  outputGridFocusY = exfoldY[i];
	  outputGridSequence = 0;
	  Log("ogfw = %d ogfh = %d\n", outputGridFocusWidth, outputGridFocusHeight);
	  OutputGrid(endLevel, -3 - i);
	}
      if (p == 0)
	{
	  sprintf(fn, "%sfolds.txt", outputGridName);
	  f = fopen(fn, "w");
	  fprintf(f, "%d\n", foldRecoveryCount);
	  fclose(f);
	}
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
  Node *nodes;
  Node *node;
  int nx, ny;
  MapElement *map;
  char msg[PATH_MAX+256];
  double a[6];
  int nPts;
  Point *pts;
  int nHullPts;
  Point *hullPts;
  int nCumulativeHullPts;
  Point *cumulativeHullPts;
  int nNewCumulativeHullPts;
  Point *newCumulativeHullPts;
  int nGlobalPts;
  Point *globalPts;
  int nGlobalHullPts;
  Point *globalHullPts;
  int *nCumulativeHullPtsFromProcess;
  int *nFloats;
  int *displacements;
  double distance;
  double theta;
  double cost, sint;
  int upperLeftCorner;
  double upperLeftDistance;
  Point corners[4];

  if (minimizeArea)
    {
      /* determine the transformation to apply to reduce the output area */
      cumulativeHullPts = NULL;
      nCumulativeHullPts = 0;
      for (i = myFirstImage; i <= myLastImage; ++i)
	{
	  nx = images[i].nx;
	  ny = images[i].ny;

	  /* compute the convex hull for just this image */
	  nPts = 0;
	  pts = (Point*) malloc(nx * ny * sizeof(Point));
	  node = images[i].nodes;
	  for (y = 0; y < ny; ++y)
	    for (x = 0; x < nx; ++x, ++node)
	      if (node->x < 0.5 * UNSPECIFIED)
		{
		  pts[nPts].x = node->x;
		  pts[nPts].y = node->y;
		  ++nPts;
		}
	  hullPts = (Point *) malloc(2 * nPts * sizeof(Point));
	  ConvexHull(nPts, pts, &nHullPts, hullPts);
	  free(pts);

	  /* compute the cumulative convex hull for this process */
	  if (cumulativeHullPts != NULL)
	    {
	      cumulativeHullPts = (Point *)
		realloc(cumulativeHullPts,
			(nCumulativeHullPts + nHullPts) * sizeof(Point)); 
	      memcpy(&cumulativeHullPts[nCumulativeHullPts],
		     hullPts,
		     nHullPts * sizeof(Point));
	      nCumulativeHullPts += nHullPts;
	      free(hullPts);
	      hullPts = NULL;

	      newCumulativeHullPts = (Point *)
		malloc(2 * nCumulativeHullPts * sizeof(Point));
	      ConvexHull(nCumulativeHullPts, cumulativeHullPts,
			 &nNewCumulativeHullPts, newCumulativeHullPts);
	      free(cumulativeHullPts);
	      cumulativeHullPts = newCumulativeHullPts;
	      nCumulativeHullPts = nNewCumulativeHullPts;
	    }
	  else
	    {
	      cumulativeHullPts = hullPts;
	      nCumulativeHullPts = nHullPts;
	      hullPts = NULL;
	    }
	}

      /* compute the cumulative convex hull over all processes */
      if (p == 0)
	nCumulativeHullPtsFromProcess = (int *) malloc(np * sizeof(int));
      else
	nCumulativeHullPtsFromProcess = NULL;
      if (MPI_Gather(&nCumulativeHullPts, 1, MPI_INT,
		     nCumulativeHullPtsFromProcess, 1, MPI_INT,
		     0, MPI_COMM_WORLD) != MPI_SUCCESS)
	Error("Could not gather nCumulativeHullPtsFromProcess\n");
      if (p == 0)
	{
	  nFloats = (int *) malloc(np * sizeof(int));
	  displacements = (int *) malloc(np * sizeof(int));
	  nGlobalPts = 0;
	  for (i = 0; i < np; ++i)
	    {
	      nGlobalPts += nCumulativeHullPtsFromProcess[i]; 
	      nFloats[i] = 2 * nCumulativeHullPtsFromProcess[i];
	      if (i > 0)
		displacements[i] = displacements[i-1] +
		  2 * nCumulativeHullPtsFromProcess[i-1];
	      else
		displacements[i] = 0;
	    }
	  free(nCumulativeHullPtsFromProcess);
	  globalPts = (Point *) malloc(nGlobalPts * sizeof(Point));
	}
      else
	{
	  globalPts = NULL;
	  nFloats = NULL;
	  displacements = NULL;
	}
      if (MPI_Gatherv(cumulativeHullPts, 2*nCumulativeHullPts, MPI_FLOAT,
		      globalPts, nFloats, displacements, MPI_FLOAT,
		      0, MPI_COMM_WORLD) != MPI_SUCCESS)
	Error("Could not gather cumulativeHullPts\n");
      if (p == 0)
	{
	  free(nFloats);
	  free(displacements);
	  globalHullPts = (Point *) malloc(2 * nGlobalPts * sizeof(Point));
	  ConvexHull(nGlobalPts, globalPts,
		     &nGlobalHullPts, globalHullPts);
	  free(globalPts);

	  /* compute the transformation to minimize the area */
	  MinimumAreaRectangle(nGlobalHullPts, globalHullPts, corners);
	  free(globalHullPts);

	  /* use the first image to pick the upper left */
	  nodes = images[0].nodes;
	  node = NULL;
	  for (y = 0; node == NULL && y < nx+ny; ++y)
	    for (x = 0; node == NULL && x < nx && x <= y; ++x)
	      if (y-x < ny && nodes[(y-x)*nx+x].x < 0.5 * UNSPECIFIED)
		node = &nodes[(y-x)*nx+x];
	  if (node == NULL)
	    Error("Could not find upper left corner in first image.\n");

	  upperLeftCorner = -1;
	  upperLeftDistance = 1.0e30;
	  for (i = 0; i < 4; ++i)
	    {
	      distance = hypot(node->x - corners[i].x,
			       node->y - corners[i].y);
	      if (distance < upperLeftDistance)
		{
		  upperLeftDistance = distance;
		  upperLeftCorner = i;
		}
	    }
	  /* we want to translate the upper left to 0, then
	     rotate so that the line from the upper left to the
	     upper right ((upperLeftCorner + 3) % 4) is along the
	     x-axis */
	  theta = -atan2(corners[(upperLeftCorner + 3) % 4].y -
			 corners[upperLeftCorner].y,
			 corners[(upperLeftCorner + 3) % 4].x -
			 corners[upperLeftCorner].x);
	  cost = cos(theta);
	  sint = sin(theta);
	  a[0] = cost;
	  a[1] = -sint;
	  a[2] = a[0] * (-corners[upperLeftCorner].x) +
	    a[1] * (-corners[upperLeftCorner].y);
	  a[3] = sint;
	  a[4] = cost;
	  a[5] = a[3] * (-corners[upperLeftCorner].x) +
	    a[4] * (-corners[upperLeftCorner].y);
	}
      if (MPI_Bcast(a, 6, MPI_DOUBLE,
		    0, MPI_COMM_WORLD) != MPI_SUCCESS)
	Error("Could not broadcast output transformation\n");
    }
  else
    {
      a[0] = 1.0;
      a[1] = 0.0;
      a[2] = 0.0;
      a[3] = 0.0;
      a[4] = 1.0;
      a[5] = 0.0;
    }

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
	  node = images[i].nodes;
	  for (y = 0; y < ny; ++y)
	    for (x = 0; x < nx; ++x)
	      {
		if (node->x < 0.5 * UNSPECIFIED)
		  {
		    map[y * nx + x].x = a[0] * node->x + a[1] * node->y + a[2];
		    map[y * nx + x].y = a[3] * node->x + a[4] * node->y + a[5];
		    map[y * nx + x].c = 1.0;
		  }
		else
		  {
		    map[y * nx + x].x = 0.0;
		    map[y * nx + x].y = 0.0;
		    map[y * nx + x].c = 0.0;
		  }
		++node;
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
  Node *node;
  float rx, ry;
  int irx, iry;
  float rrx, rry;
  float rx00, rx01, rx10, rx11;
  float ry00, ry01, ry10, ry11;
  int nx, ny;
  int ix, iy;
  int dx, dy;
  int xp, yp;
  unsigned char *grid, *combinedGrid;
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
  float globalMinX, globalMaxX, globalMinY, globalMaxY;
  int fIndex;
  int ci;
  static float minX, maxX, minY, maxY;
  static float scale;
  static float offsetX, offsetY;

  Log("OutputGrid called at level = %d iter = %d\n", level, iter);
  factor = 1 << level;
  fIndex = -1;
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
		fIndex = i;
		if (images[i].owner != p)
		  break;
		nx = images[i].nx;
		ny = images[i].ny;
		node = images[i].nodes;
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
		rx00 = node[iyv*nx+ixv].x;
		ry00 = node[iyv*nx+ixv].y;
		rx01 = node[(iyv+1)*nx+ixv].x;
		ry01 = node[(iyv+1)*nx+ixv].y;
		rx10 = node[iyv*nx+ixv+1].x;
		ry10 = node[iyv*nx+ixv+1].y;
		rx11 = node[(iyv+1)*nx+ixv+1].x;
		ry11 = node[(iyv+1)*nx+ixv+1].y;
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
		Log("rx = %f ry = %f ogfw = %d ofgh = %d factor = %d\n",
		    rx, ry, outputGridFocusWidth, outputGridFocusHeight, factor);
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
	    node = images[i].nodes;
	    for (y = 0; y < ny; ++y)
	      for (x = 0; x < nx; ++x, ++node)
		{
		  if (node->x > 0.5 * UNSPECIFIED)
		    continue;
		  if (node->x < minX)
		    minX = node->x;
		  if (node->x > maxX)
		    maxX = node->x;
		  if (node->y < minY)
		    minY = node->y;
		  if (node->y > maxY)
		    maxY = node->y;
		}
	  }
      if (MPI_Allreduce(&minX, &globalMinX, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD) != MPI_SUCCESS ||
	  MPI_Allreduce(&maxX, &globalMaxX, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD) != MPI_SUCCESS ||
	  MPI_Allreduce(&minY, &globalMinY, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD) != MPI_SUCCESS ||
	  MPI_Allreduce(&maxY, &globalMaxY, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD) != MPI_SUCCESS)
	Error("MPI_Allreduce: could not compute grid limits.\n");
      
      if (globalMinX == globalMaxX || globalMinY == globalMaxY)
	Error("maps have 0 extent in x or y\n");
      minX = globalMinX * factor;
      maxX = globalMaxX * factor;
      minY = globalMinY * factor;
      maxY = globalMaxY * factor;
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
      if (iter < -2)
	sprintf(gridDirName, "%sexfold%0.3d", outputGridName, -(iter + 3));
      else if (iter == -2)
	sprintf(gridDirName, "%sfold", outputGridName);
      else
	sprintf(gridDirName, "%s", outputGridName);
    }
  else
    sprintf(gridDirName, "%sl%0.2di%0.6d",
	    outputGridName, level, iter);

  n = outputGridWidth * outputGridHeight * 3;
  grid = (unsigned char *) malloc(n);
  combinedGrid = (unsigned char *) malloc(n);
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
      if (outputGridFocusDepth >= 0 && fIndex >= 0 &&
	  abs(i - fIndex) > outputGridFocusDepth)
	continue;
      if (outputGridFocusDepth >= 0)
	ci = (i - fIndex) & 15;
      else
	ci = i & 15;
      nx = images[i].nx;
      ny = images[i].ny;
      node = images[i].nodes;
#if 0
      if (node[0].x < bMinX && node[nx-1].x < bMinX &&
	  node[(ny-1)*nx].x < bMinX && node[(ny-1)*nx+nx-1].x < bMinX)
	continue;
      if (node[0].x > bMaxX && node[nx-1].x > bMaxX &&
	  node[(ny-1)*nx].x > bMaxX && node[(ny-1)*nx+nx-1].x > bMaxX)
	continue;
      if (node[0].y < bMinY && node[nx-1].y < bMinY &&
	  node[(ny-1)*nx].y < bMinY && node[(ny-1)*nx+nx-1].y < bMinY)
	continue;
      if (node[0].y > bMaxY && node[nx-1].y > bMaxY &&
	  node[(ny-1)*nx].y > bMaxY && node[(ny-1)*nx+nx-1].y > bMaxY)
	continue;
#endif
      for (y = 0; y < ny; ++y)
	for (x = 0; x < nx; ++x)
	  {
	    if (node[y*nx+x].x > 0.5 * UNSPECIFIED)
	      continue;
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
		    
		    if (node[yp*nx+xp].x > 0.5 * UNSPECIFIED ||
			node[y*nx+x].x < bMinX && node[yp*nx+xp].x < bMinX ||
			node[y*nx+x].x > bMaxX && node[yp*nx+xp].x > bMaxX ||
			node[y*nx+x].y < bMinY && node[yp*nx+xp].y < bMinY ||
			node[y*nx+x].y > bMaxY && node[yp*nx+xp].y > bMaxY)
		      continue;

#if 0
		    nomD = (dx * dy != 0) ? M_SQRT2 : 1.0;
		    dist = hypot(node[y*nx+x].x - node[yp*nx+xp].x,
				 node[y*nx+x].y - node[yp*nx+xp].y);
		    hue = 4.0 * (dist / nomD) - 2.0;
		    if (hue < 0.0)
		      hue = 0.0;
		    if (hue > 5.0)
		      hue = 5.0;
		    hsv_to_rgb(&r, &g, &b, hue, 1.0, 1.0);
#else
		    r = colors[ci][0];
		    g = colors[ci][1];
		    b = colors[ci][2];
#endif

		    DrawLine(grid,
			     scale * factor * node[y*nx+x].x + offsetX,
			     scale * factor * node[y*nx+x].y + offsetY,
			     scale * factor * node[yp*nx+xp].x + offsetX,
			     scale * factor * node[yp*nx+xp].y + offsetY,
			     outputGridWidth, outputGridHeight,
			     r, g, b);
		  }
	      }
	  }

      for (y = 0; y < ny; ++y)
	for (x = 0; x < nx; ++x)
	  {
	    if (node[y*nx+x].x > 0.5 * UNSPECIFIED)
	      continue;
	    if (node[y*nx+x].fx > 0.5 * CONSTRAINED ||
		node[y*nx+x].fy < 0.5 * CONSTRAINED)
	      {
		ix = (int) (scale * factor * node[y*nx+x].x + offsetX);
		iy = (int) (scale * factor * node[y*nx+x].y + offsetY);
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
      if ((code = MPI_Allreduce(grid, combinedGrid, n, MPI_UNSIGNED_CHAR,
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
	  if (fwrite(combinedGrid, sizeof(unsigned char), n, f) != n)
	    Error("Could not write to grid file %s\n", fn);
	  fclose(f);
	}
    }
  free(grid);
  free(combinedGrid);
  ++outputGridSequence;
}

void OutputFoldMaps (int level, int iter)
{
  int i, j, k, q, r;
  int x, y;
  int nx, ny;
  Node *nodes;
  int minX, maxX, minY, maxY;
  float xv, yv;
  float ax, ay, bx, by;
  float xMin, xMax, yMin, yMax;
  unsigned char *valid;
  int snx, sny;
  MapElement *submap;
  float kdx, kdy;
  InterImageMap *m;
  int mark;
  int ix, iy;
  int delta;
  char fn[PATH_MAX];
  char msg[PATH_MAX+256];
  int placed;
  int occupied[7];
  SpringForce *springForces;
  float msk;
  float sk;
  int irx, iry;
  int *code;
  int nStrips;
  InterImageStrip *strips;
  InterImageSpring *springs;
  Node *nodes0, *nodes1;
  int nx0, ny0, nx1, ny1;
  int dir;
  int springsPos;
  InterImageStrip *strip;
  InterImageSpring *s;
  int ns;
  float deltaX, deltaY;
  int js;

  xMin = foldAbsoluteX - 8.0;
  xMax = foldAbsoluteX + 8.0;
  yMin = foldAbsoluteY - 8.0;
  yMax = foldAbsoluteY + 8.0;
  for (i = myFirstImage; i <= myLastImage; ++i)
    {
      /* only output images within a given distance of fold */
      if (abs(i - foldImage) > 2)
	continue;
      nx = images[i].nx;
      ny = images[i].ny;
      nodes = images[i].nodes;
      valid = (unsigned char *) malloc(nx * ny * sizeof(unsigned char));
      memset(valid, 0, nx * ny * sizeof(unsigned char));
      minX = 1000000000;
      maxX = -1000000000;
      minY = 1000000000;
      maxY = -1000000000;
      for (y = 0; y < ny-1; ++y)
	for (x = 0; x < nx-1; ++x)
	  {
	    if (nodes[y*nx+x].x > 0.5 * UNSPECIFIED ||
		nodes[y*nx+x+1].x > 0.5 * UNSPECIFIED ||
		nodes[(y+1)*nx+x].x > 0.5 * UNSPECIFIED ||
		nodes[(y+1)*nx+x+1].x > 0.5 * UNSPECIFIED)
	      continue;

	    if (nodes[y*nx+x].x < xMin || nodes[y*nx+x].x > xMax ||
		nodes[y*nx+x].y < yMin || nodes[y*nx+x].y > yMax ||
		nodes[y*nx+x+1].x < xMin || nodes[y*nx+x+1].x > xMax ||
		nodes[y*nx+x+1].y < yMin || nodes[y*nx+x+1].y > yMax ||
		nodes[(y+1)*nx+x].x < xMin || nodes[(y+1)*nx+x].x > xMax ||
		nodes[(y+1)*nx+x].y < yMin || nodes[(y+1)*nx+x].y > yMax ||
		nodes[(y+1)*nx+x+1].x < xMin || nodes[(y+1)*nx+x+1].x > xMax ||
		nodes[(y+1)*nx+x+1].y < yMin || nodes[(y+1)*nx+x+1].y > yMax)
	      continue;

	    /* eliminate the folded elements from the map */
	    xv = nodes[y * nx + x].x;
	    yv = nodes[y * nx + x].y;
	    if (x > 0)
	      {
		ax = nodes[y * nx + x - 1].x - xv;
		ay = nodes[y * nx + x - 1].y - yv;
		if (y > 0)
		  {
		    bx = nodes[(y-1) * nx + x].x - xv;
		    by = nodes[(y-1) * nx + x].y - yv;
		    if (ax < 0.5 * UNSPECIFIED &&
			bx < 0.5 * UNSPECIFIED &&
			ax * by - ay * bx <= 0.0)
		      continue;
		  }
		if (y < ny - 1)
		  {
		    bx = nodes[(y+1) * nx + x].x - xv;
		    by = nodes[(y+1) * nx + x].y - yv;
		    if (ax < 0.5 * UNSPECIFIED &&
			bx < 0.5 * UNSPECIFIED &&
			ax * by - ay * bx >= 0.0)
		      continue;
		  }
	      }
	    if (x < nx - 1)
	      {
		ax = nodes[y * nx + x + 1].x - xv;
		ay = nodes[y * nx + x + 1].y - yv;
		if (y > 0)
		  {
		    bx = nodes[(y-1) * nx + x].x - xv;
		    by = nodes[(y-1) * nx + x].y - yv;
		    if (ax < 0.5 * UNSPECIFIED &&
			bx < 0.5 * UNSPECIFIED &&
			ax * by - ay * bx >= 0.0)
		      continue;
		  }
		if (y < ny - 1)
		  {
		    bx = nodes[(y+1) * nx + x].x - xv;
		    by = nodes[(y+1) * nx + x].y - yv;
		    if (ax < 0.5 * UNSPECIFIED &&
			bx < 0.5 * UNSPECIFIED &&
			ax * by - ay * bx <= 0.0)
		      continue;
		  }
	      }

	    valid[y*nx+x] = 1;
	    if (x < minX)
	      minX = x;
	    if (x > maxX)
	      maxX = x;
	    if (y < minY)
	      minY = y;
	    if (y > maxY)
	      maxY = y;
	  }
      if (maxX < minX | maxY < minY)
	{
	  free(valid);
	  continue;
	}

      snx = maxX - minX + 2;
      sny = maxY - minY + 2;
      Log("Allocating submap of size %d x %d\n", snx, sny);
      submap = malloc(snx * sny * sizeof(MapElement));
      memset(submap, 0, snx * sny * sizeof(MapElement));
      springForces = (SpringForce*) malloc(snx * sny * 7 * sizeof(SpringForce));
      memset(springForces, 0, snx * sny * 7 * sizeof(SpringForce));
      for (j = 0; j < nMaps; ++j)
	{
	  m = &maps[j];
	  if (m->image0 == i)
	    q = m->image1;
	  else if (m->image1 == i)
	    q = m->image0;
	  else
	    continue;
	  if (q < i)
	    delta = -1;
	  else
	    delta = 1;
	  for (k = 0; k < nMaps; ++k)
	    {
	      if (maps[k].image0 == i)
		r = maps[k].image1;
	      else if (maps[k].image1 == i)
		r = maps[k].image0;
	      else continue;
	      if (q < i && r < i && q < r)
		--delta;
	      if (q > i && r > i && q > r)
		++delta;
	    }
	  if (delta < -3)
	    delta = 1;
	  else if (delta < 0)
	    delta = 4 + delta;
	  else if (delta > 4)
	    delta = 7;
	  else
	    delta = delta + 3;

#if 1
	  msk = kInter * m->k;
	  nStrips = m->nStrips[startLevel - level];
	  strips = m->strips[startLevel - level];
	  springs = m->springs[startLevel - level];
	  nodes0 = images[m->image0].nodes;
	  if (nodes0 == NULL)
	    abort();
	  nodes1 = images[m->image1].nodes;
	  nx0 = images[m->image0].nx;
	  nx1 = images[m->image1].nx;
	  springsPos = 0;
	  for (js = 0; js < nStrips; ++js)
	    {
	      if (nodes0 == NULL)
		abort();
	      strip = &(strips[js]);
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
		  if (nodes0[y*nx0+x].x > 0.5 * UNSPECIFIED ||
		      nodes1[iry*nx1+irx].x > 0.5 * UNSPECIFIED)
		    continue;
		  xv = nodes1[iry*nx1+irx].x + 0.00005 * s->dx;
		  yv = nodes1[iry*nx1+irx].y + 0.00005 * s->dy;
		  deltaX = xv - nodes0[y*nx0 + x].x;
		  deltaY = yv - nodes0[y*nx0 + x].y;
		  sk = (s->k / 255.0) * msk;
		  kdx = sk * deltaX;
		  kdy = sk * deltaY;
		  mark = 0;
		  if  (m->image0 == i &&
		       x >= minX && x <= maxX &&
		       y >= minY && y <= maxY)
		    {
		      mark = 1;
		      ix = x - minX;
		      iy = y - minY;
		    }
		  else if (m->image1 == i &&
			   irx >= minX && irx <= maxX &&
			   iry >= minY && iry <= maxY)
		    {
		      mark = 1;
		      ix = irx - minX;
		      iy = iry - minY;
		      kdx = -kdx;
		      kdy = -kdy;
		    }
		  if (mark)
		    {
		      if (ix < 0 || ix >= snx)
			abort();
		      if (iy < 0 || iy >= sny)
			abort();
		      for (q = 0; q < 7; ++q)
			if (springForces[(iy * snx + ix) * 7 + q].delta == 0)
			  {
			    springForces[(iy * snx + ix) * 7 + q].delta = delta;
			    //			    if (delta != 3)
			    //			      abort();
			    springForces[(iy * snx + ix) * 7 + q].theta = atan2f(kdy, kdx);
			    springForces[(iy * snx + ix) * 7 + q].force = hypot(kdx, kdy);
			    break;
			}
		    }
		}
	    }
#endif
	}
      code = (int *) malloc(snx * sny * sizeof(int));
      for (y = 0; y < sny; ++y)
	for (x = 0; x < snx; ++x)
	  {
	    code[y*snx + x] = 1;
	    qsort(&springForces[(y * snx + x) * 7], 7, sizeof(SpringForce), CompareSpringForces);
	    memset(occupied, 0, 7*sizeof(int));
	    for (k = 0; k < 7; ++k)
	      {
		delta = springForces[(y * snx + x) * 7 + k].delta;
		if (delta == 0)
		  break;
		dir = (int) floor((springForces[(y * snx + x) * 7 + k].theta * 180.0 / M_PI + 330.0) / 60.0);
		dir = dir % 6;
		placed = 0;
		for (q = 0; !placed && q < 4; ++q)
		  {
		    if (!occupied[(dir+q) % 6])
		      {
			code[y * snx + x] |= delta << (((dir+q) % 6) * 3 + 1);
			occupied[(dir+q) % 6] = 1;
			placed = 1;
		      }
		    else if (!occupied[(dir-q+6) % 6])
		      {
			code[y * snx + x] |= delta << (((dir-q+6) % 6) * 3 + 1);
			occupied[(dir-q+6) % 6] = 1;
			placed = 1;
		      }
		  }
		if (!placed && !occupied[6])
		  {
		    code[y * snx + x] |= delta << (6 * 3 + 1);
		    occupied[6] = 1;
		  }
	      }
	  }

      for (y = 0; y < sny-1; ++y)
	for (x = 0; x < snx-1; ++x)
	  {
	    if (!valid[(y+minY)*nx + x + minX])
	      continue;
	    submap[y*snx + x].x = nodes[(y+minY)*nx + x + minX].x;
	    submap[y*snx + x].y = nodes[(y+minY)*nx + x + minX].y;
	    submap[y*snx + x].c = (float) code[y*snx + x];
	    submap[y*snx + x + 1].x = nodes[(y+minY)*nx + x + 1 + minX].x;
	    submap[y*snx + x + 1].y = nodes[(y+minY)*nx + x + 1 + minX].y;
	    submap[y*snx + x + 1].c = (float) code[y*snx + x + 1];
	    submap[(y+1)*snx + x].x = nodes[(y+1+minY)*nx + x + minX].x;
	    submap[(y+1)*snx + x].y = nodes[(y+1+minY)*nx + x + minX].y;
	    submap[(y+1)*snx + x].c = (float) code[(y+1)*snx + x];
	    submap[(y+1)*snx + x + 1].x = nodes[(y+1+minY)*nx + x + 1 + minX].x;
	    submap[(y+1)*snx + x + 1].y = nodes[(y+1+minY)*nx + x + 1 + minX].y;
	    submap[(y+1)*snx + x + 1].c = (float) code[(y+1)*snx + x + 1];
	  }
      sprintf(fn, "%sfold.%s.map", outputGridName, images[i].name);
      if (!WriteMap(fn, submap, level, snx, sny, minX, minY,
		    images[i].name, "align12_fold",
		    UncompressedMap, msg))
	Error("Could not write map file %s: %s\n", fn, msg);
      free(submap);
      free(code);
      free(valid);
      free(springForces);
    }
}

int
CompareSpringForces (const void *p0, const void *p1)
{
  if (((const SpringForce *) p0)->force > ((const SpringForce *) p1)->force)
    return(-1);
  else if (((const SpringForce *) p0)->force > ((const SpringForce *) p1)->force)
    return(1);
  else
    return(0);
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
  Node *nodes1;
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
      nodes1 = images[m->image1].nodes;

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
	      rrx = s->irrx / 20000.0;
	      rry = s->irry / 20000.0;
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
	      rx00 = nodes1[iry*nx1 + irx].x;
	      rx01 = nodes1[(iry+1)*nx1+irx].x;
	      rx10 = nodes1[iry*nx1 + irx+1].x;
	      rx11 = nodes1[(iry+1)*nx1+irx+1].x;
	      ry00 = nodes1[iry*nx1 + irx].y;
	      ry01 = nodes1[(iry+1)*nx1+irx].y;
	      ry10 = nodes1[iry*nx1 + irx+1].y;
	      ry11 = nodes1[(iry+1)*nx1+irx+1].y;
	      if (rx00 > 0.5 * UNSPECIFIED ||
		  rx01 > 0.5 * UNSPECIFIED ||
		  rx10 > 0.5 * UNSPECIFIED ||
		  rx11 > 0.5 * UNSPECIFIED)
		{
		  s->dx = 0;
		  s->dy = 0;
		  s->k = 0;
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
	      deltaX = 20000.0 * (rx - nodes1[iy*nx1+ix].x) + 0.5;
	      deltaY = 20000.0 * (ry - nodes1[iy*nx1+ix].y) + 0.5;
	      if (deltaX >= -32768.0 && deltaX < 32768.0 &&
		  deltaY >= -32768.0 && deltaY < 32768.0)
		{
		  s->dx = (int) floor(deltaX);
		  s->dy = (int) floor(deltaY);
		}
	      else
		{
		  s->dx = 0;
		  s->dy = 0;
		  s->k = 0;
		}
	    }
	}
    }
}

void
RefinePositions (int prevLevel, int level)
{
  int i;
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
  Node *nodes, *nodes1;
  float dFactor;
  int eFactor;
  int rFactor;
  unsigned char *imask;
  int rnx, rny;
  int imbpl;
  int valid;
  int ix, iy;
  int count;

  factor1 = 1 << level;
  dFactor = 1 << (prevLevel - level);
  eFactor = 1 << endLevel;
  rFactor = factor1 / eFactor;
  for (i = 0; i < nImages; ++i)
    {
      if (images[i].owner != p && !images[i].needed)
	continue;
      nx = images[i].nx;
      ny = images[i].ny;
      nx1 = (images[i].width + factor1 - 1) / factor1 + 1;
      ny1 = (images[i].height + factor1 - 1) / factor1 + 1;
      images[i].nx = nx1;
      images[i].ny = ny1;
      nodes = images[i].nodes;
      nodes1 = (Node*) malloc(nx1 * ny1 * sizeof(Node));
      images[i].nodes = nodes1;
      if (images[i].owner != p)
	{
	  memset(nodes1, 0, nx1 * ny1 * sizeof(Node));
	  free(nodes);
	  continue;
	}
      for (y = 0; y < ny1; ++y)
	for (x = 0; x < nx1; ++x)
	  {
	    nodes1[y * nx1 + x].x = UNSPECIFIED;
	    nodes1[y * nx1 + x].y = UNSPECIFIED;
	    nodes1[y * nx1 + x].fx = 0.0;
	    nodes1[y * nx1 + x].fy = 0.0;
	  }
      // mark the valid map nodes
      rnx = (images[i].width + eFactor - 1) / eFactor;
      rny = (images[i].height + eFactor - 1) / eFactor;
      imbpl = (rnx + 7) >> 3;
      imask = images[i].mask;
      for (y = 0; y < ny1; ++y)
	for (x = 0; x < nx1; ++x)
	  {
	    valid = (imask == NULL);
	    for (iy = 0; !valid && iy < rFactor; ++iy)
	      {
		iyv = y * rFactor + iy;
		if (iyv >= rny)
		  break;
		for (ix = 0; ix < rFactor; ++ix)
		  {
		    ixv = x * rFactor + ix;
		    if (ixv >= rnx)
		      break;
		    if (imask[iyv * imbpl + (ixv >> 3)] & (0x80 >> (ixv & 7)))
		      {
			valid = 1;
			break;
		      }
		  }
	      }
	    if (valid)
	      {
		nodes1[y * nx1 + x].x = 0.0;
		nodes1[y * nx1 + x].y = 0.0;
		if (x < nx1 - 1)
		  {
		    nodes1[y * nx1 + x + 1].x = 0.0;
		    nodes1[y * nx1 + x + 1].y = 0.0;
		  }
		if (y < ny1 - 1)
		  {
		    nodes1[(y + 1) * nx1 + x].x = 0.0;
		    nodes1[(y + 1) * nx1 + x].y = 0.0;
		  }
		if (x < nx1 - 1 && y < ny1 - 1)
		  {
		    nodes1[(y + 1) * nx1 + x + 1].x = 0.0;
		    nodes1[(y + 1) * nx1 + x + 1].y = 0.0;
		  }
	      }
	  }

      count = 0;
      for (y = 0; y < ny1; ++y)
	for (x = 0; x < nx1; ++x)
	  {
	    if (nodes1[y * nx1 + x].x > 0.5 * UNSPECIFIED)
	      continue;
	    ++count;
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
	    rx00 = nodes[iyv * nx + ixv].x;
	    ry00 = nodes[iyv * nx + ixv].y;
	    rx01 = nodes[(iyv+1) * nx + ixv].x;
	    ry01 = nodes[(iyv+1) * nx + ixv].y;
	    rx10 = nodes[iyv * nx + ixv + 1].x;
	    ry10 = nodes[iyv * nx + ixv + 1].y;
	    rx11 = nodes[(iyv+1) * nx + ixv + 1].x;
	    ry11 = nodes[(iyv+1) * nx + ixv + 1].y;
	    if (rx00 > 0.5 * UNSPECIFIED ||
		rx01 > 0.5 * UNSPECIFIED && rry != 0.0 ||
		rx10 > 0.5 * UNSPECIFIED && rrx != 0.0 ||
		rx11 > 0.5 * UNSPECIFIED && rrx != 0.0 && rry != 0.0)
	      Error("RefinePositions: could not interpolate from invalid nodes (image = %s x = %d y = %d nx1 = %d ny1 = %d ixv = %d, iyv = %d nx = %d ny = %d xv = %f yv = %f rrx = %f rry = %f dFactor = %f)\n",
		    images[i].name, x, y, nx1, ny1,
		    ixv, iyv, nx, ny,
		    xv, yv, rrx, rry, dFactor);
	    rx = rx00 * (rrx - 1.0) * (rry - 1.0)
	      - rx10 * rrx * (rry - 1.0) 
	      - rx01 * (rrx - 1.0) * rry
	      + rx11 * rrx * rry;
	    ry = ry00 * (rrx - 1.0) * (rry - 1.0)
	      - ry10 * rrx * (rry - 1.0) 
	      - ry01 * (rrx - 1.0) * rry
	      + ry11 * rrx * rry;

	    nodes1[y * nx1 + x].x = dFactor * rx;
	    nodes1[y * nx1 + x].y = dFactor * ry;
	    nodes1[y * nx1 + x].fx = 0.0;
	    nodes1[y * nx1 + x].fy = 0.0;
	  }
      Log("Refined map of image %s (%d of %d nodes valid)\n",
	  images[i].name, count, nx1*ny1);
      free(nodes);
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
  int nNodes;
  int its;
  Node *nodes;
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
		      nNodes = images[its].nx * images[its].ny;
		      nodes = images[its].nodes;
		      for (k = 0; k < nNodes; ++k)
			{
			  buffer[bufferPos++] = nodes[k].x;
			  buffer[bufferPos++] = nodes[k].y;
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
		      nNodes = images[its].nx * images[its].ny;
		      nodes = images[its].nodes;
		      for (k = 0; k < nNodes; ++k)
			{
			  nodes[k].x =  buffer[bufferPos++];
			  nodes[k].y =  buffer[bufferPos++];
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
CheckForFolds (int level, int iter)
{
  int i;
  int nx, ny;
  int x, y;
  float xv, yv;
  Node *nodes;
  float ax, ay;
  float bx, by;
  int ix, iy;

  for (i = myFirstImage; i <= myLastImage; ++i)
    {
      nx = images[i].nx;
      ny = images[i].ny;
      nodes = images[i].nodes;
      for (y = 0; y < ny; ++y)
	for (x = 0; x < nx; ++x)
	  {
	    xv = nodes[y * nx + x].x;
	    if (xv > 0.5 * UNSPECIFIED)
	      continue;
	    yv = nodes[y * nx + x].y;
	    if (x > 0)
	      {
		ax = nodes[y * nx + x - 1].x - xv;
		ay = nodes[y * nx + x - 1].y - yv;
		if (y > 0)
		  {
		    bx = nodes[(y-1) * nx + x].x - xv;
		    by = nodes[(y-1) * nx + x].y - yv;
		    if (ax < 0.5 * UNSPECIFIED &&
			bx < 0.5 * UNSPECIFIED &&
			ax * by - ay * bx <= 0.0)
		      goto foldError;
		  }
		if (y < ny - 1)
		  {
		    bx = nodes[(y+1) * nx + x].x - xv;
		    by = nodes[(y+1) * nx + x].y - yv;
		    if (ax < 0.5 * UNSPECIFIED &&
			bx < 0.5 * UNSPECIFIED &&
			ax * by - ay * bx >= 0.0)
		      goto foldError;
		  }
	      }
	    if (x < nx - 1)
	      {
		ax = nodes[y * nx + x + 1].x - xv;
		ay = nodes[y * nx + x + 1].y - yv;
		if (y > 0)
		  {
		    bx = nodes[(y-1) * nx + x].x - xv;
		    by = nodes[(y-1) * nx + x].y - yv;
		    if (ax < 0.5 * UNSPECIFIED &&
			bx < 0.5 * UNSPECIFIED &&
			ax * by - ay * bx >= 0.0)
		      goto foldError;
		  }
		if (y < ny - 1)
		  {
		    bx = nodes[(y+1) * nx + x].x - xv;
		    by = nodes[(y+1) * nx + x].y - yv;
		    if (ax < 0.5 * UNSPECIFIED &&
			bx < 0.5 * UNSPECIFIED &&
			ax * by - ay * bx <= 0.0)
		      goto foldError;
		  }
	      }
	  }
    }
  return(nImages);

 foldError:
  foldImage = i;
  foldWidth = (1 << level) * 32;
  foldHeight = (1 << level) * 32;
  foldX = (1 << level) * x;
  foldY = (1 << level) * y;
  foldAbsoluteX = nodes[y * nx + x].x;
  foldAbsoluteY = nodes[y * nx + x].y;

#if 0
  Log("DUMPING FOLDED MAP for %s:\n", images[i].name);
  for (iy = 0; iy < ny; ++iy)
    for (ix = 0; ix < nx; ++ix)
      Log("x = %d y = %d:  mapx = %f mapy = %f\n",
	  ix, iy, nodes[iy*nx+ix].x, nodes[iy*nx+ix].y);
#endif

  Log("Fold detected in map of image %s at (%f %f) (%d(%d %d)) at iteration %d\n",
      images[i].name,
      (1 << level) * ((double) x), (1 << level) * ((double) y),
      level, x, y, iter);
  return(i);
}

int
Extrapolate (float *prx, float *pry, float *prc,
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

	  //	  Log("extrap rx00=%f rx01=%f rx10=%f rx11=%f\n",
	  //	      rx00, rx01, rx10, rx11);
	  //	  Log("extrap ry00=%f ry01=%f ry10=%f ry11=%f\n",
	  //	      ry00, ry01, ry10, ry11);
	  //	  Log("extrap term: dx=%d dy=%d rrx=%f rry=%f weight=%f drx=%f dry=%f\n",
	  //	      dx, dy, rrx, rry, weight,
	  //	      rx00 * (rrx - 1.0) * (rry - 1.0)
	  //	      - rx10 * rrx * (rry - 1.0) 
	  //	      - rx01 * (rrx - 1.0) * rry
	  //	      + rx11 * rrx * rry,
	  //	      ry00 * (rrx - 1.0) * (rry - 1.0)
	  //	      - ry10 * rrx * (rry - 1.0) 
	  //	      - ry01 * (rrx - 1.0) * rry
	  //	      + ry11 * rrx * rry);

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

void
RecoverFromFold (int level)
{
  float myWorstK;
  int myWorstMap;
  float worstK;
  MapElement *map;
  InterImageMap *m;
  int nStrips;
  InterImageStrip *strips;
  InterImageStrip *strip;
  InterImageSpring *springs;
  InterImageSpring *s;
  Node *nodes0;
  int i, j, k;
  int ns;
  char *worstMapName;
  int worstMapNameSize;
  int imageNameSize;
  float area;
  float l0, l1, l2, l3;
  float xv, yv;
  int nx, ny;
  int nx1;
  int springsPos;
  int dPoints;
  int x, y;
  int irx, iry;
  float x0, x1, x2, x3;
  float y0, y1, y2, y3;
  int processWithWorstMap;
  int minProcessWithWorstMap;
  int lvl;
  int dLevel;
  unsigned char *cancel;
  int processWithFold;
  unsigned char *myVector, *markedVector;
  int nExpansions;
  int expansion;

  // broadcast the fold location
  processWithFold = images[minImageWithFold].owner;
  if (MPI_Bcast(&foldImage, 1, MPI_INT,
		processWithFold, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&foldWidth, 1, MPI_INT,
		processWithFold, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&foldHeight, 1, MPI_INT,
		processWithFold, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&foldX, 1, MPI_INT,
		processWithFold, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&foldY, 1, MPI_INT,
		processWithFold, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&foldAbsoluteX, 1, MPI_FLOAT,
		processWithFold, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&foldAbsoluteX, 1, MPI_FLOAT,
		processWithFold, MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("Could not broadcast fold location\n");

  // decide which nearby InterImageMap has the lowest weight

  // mark the folded image
  for (i = 0; i < nImages; ++i)
    images[i].marked = (i == foldImage);
  myVector = (unsigned char *) malloc(nImages);
  markedVector = (unsigned char *) malloc(nImages);

  nExpansions = 2;
  for (expansion = 0; expansion < nExpansions; ++expansion)
    {
      memset(myVector, 0, nImages);
      // find all maps involving a marked image
      myWorstK = 1.0e+30;
      myWorstMap = -1;
      for (i = 0; i < nMaps; ++i)
	{
	  m = &maps[i];
	  if (images[m->image0].marked ||
	      images[m->image1].marked)
	    {
	      if (m->k < myWorstK)
		{
		  myWorstMap = i;
		  myWorstK = m->k;
		}
	      myVector[m->image0] = 1;
	      myVector[m->image1] = 1;
	    }
	}

      // also mark all the images connected by a map
      if (MPI_Allreduce(myVector, markedVector, nImages, MPI_BYTE,
			MPI_BOR, MPI_COMM_WORLD) != MPI_SUCCESS)
	Error("MPI_Allreduce: could not 'or' marked vector.\n");
      for (i = 0; i < nImages; ++i)
	if (markedVector[i])
	  images[i].marked = 1;
    }
  free(myVector);
  free(markedVector);

  // find the worst map
  if (MPI_Allreduce(&myWorstK, &worstK, 1, MPI_FLOAT,
		    MPI_MIN, MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("MPI_Allreduce: could not find min k.\n");
  if (myWorstK == worstK)
    processWithWorstMap = p;
  else
    processWithWorstMap = np;
  if (MPI_Allreduce(&processWithWorstMap, &minProcessWithWorstMap,
		    1, MPI_INT, MPI_MIN, MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("MPI_Allreduce: could not find process with worst map.\n");
  if (p == minProcessWithWorstMap)
    worstMapNameSize = strlen(maps[myWorstMap].name) + 1;
  if (MPI_Bcast(&worstMapNameSize, 1, MPI_INT,
		minProcessWithWorstMap, MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("Could not broadcast worst map name size\n");
  
  worstMapName = (char *) malloc(worstMapNameSize);
  if (p == minProcessWithWorstMap)
    {
      strcpy(worstMapName, maps[myWorstMap].name);
      exfoldImage[foldRecoveryCount] = maps[myWorstMap].image0;
      exfoldWidth[foldRecoveryCount] = foldWidth;
      exfoldHeight[foldRecoveryCount] = foldHeight;
      exfoldX[foldRecoveryCount] = foldX;
      exfoldY[foldRecoveryCount] = foldY;
    }
  if (MPI_Bcast(worstMapName, worstMapNameSize, MPI_CHAR,
		minProcessWithWorstMap, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&exfoldImage[foldRecoveryCount], 1, MPI_INT,
		minProcessWithWorstMap, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&exfoldWidth[foldRecoveryCount], 1, MPI_INT,
		minProcessWithWorstMap, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&exfoldHeight[foldRecoveryCount], 1, MPI_INT,
		minProcessWithWorstMap, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&exfoldX[foldRecoveryCount], 1, MPI_INT,
		minProcessWithWorstMap, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&exfoldY[foldRecoveryCount], 1, MPI_INT,
		minProcessWithWorstMap, MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("Could not broadcast fold information\n");

  Log("exfoldWidth[%d] = %d\n", foldRecoveryCount, exfoldWidth[foldRecoveryCount]);

  if (p == 0)
    Log("Discarding map %s with weight %f to try to eliminate fold\n",
	worstMapName, worstK);

  for (i = 0; i < nMaps; ++i)
    if (strcmp(maps[i].name, worstMapName) == 0)
      maps[i].k = 0.0;
  free(worstMapName);

  // restore point positions
  for (i = myFirstImage; i <= myLastImage; ++i)
    {
      nx = images[i].nx;
      ny = images[i].ny;
      memcpy(images[i].nodes, images[i].initialNodes,
	     nx * ny * sizeof(Node));
    }
}

void
ConvexHull (int n, Point *pts,
	    int *hn, Point *hpts)
{
  int i;
  int k = 0;
  int t;

  // Sort points lexicographically
  qsort(pts, n, sizeof(Point), ComparePoints);
 
  // Build lower hull
  for (i = 0; i < n; ++i)
    {
      while (k >= 2 && Cross(hpts[k-2].x, hpts[k-2].y,
			     hpts[k-1].x, hpts[k-1].y,
			     pts[i].x, pts[i].y) <= 0.0)
	--k;
      hpts[k++] = pts[i];
    }
 
  // Build upper hull
  t = k + 1;
  for (i = n - 2; i >= 0; --i)
    {
      while (k >= t && Cross(hpts[k-2].x, hpts[k-2].y,
			     hpts[k-1].x, hpts[k-1].y,
			     pts[i].x, pts[i].y) <= 0.0)
	--k;
      hpts[k++] = pts[i];
    }

  *hn = k-1;
}

int
ComparePoints (const void *p0, const void *p1)
{
  if (((const Point *) p0)->x < ((const Point *) p1)->x)
    return(-1);
  else if (((const Point *) p0)->x > ((const Point *) p1)->x)
    return(1);
  else if (((const Point *) p0)->y < ((const Point *) p1)->y)
    return(-1);
  else if (((const Point *) p0)->y > ((const Point *) p1)->y)
    return(1);
  else
    return(0);
}


double
Cross (double x0, double y0,
       double x1, double y1,
       double x2, double y2)
{
  return ((x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0));
}



double
MinimumAreaRectangle (int n, Point *p, Point *corners)
{
  int i, j, k, m;
  double alpha, beta, gamma, delta;
  double slope;
  double d1, d2;
  double A, Amin;
  int start;
  Point p2;

  /* convert pts to have the point with the smallest y coordinate first */
  start = 0;
  for (i = 1; i < n; ++i)
    if (p[i].y < p[start].y ||
	p[i].y == p[start].y && p[i].x < p[start].x)
      start = i;
  printf("start = %d\n", start);
  
  alpha = 0.0;
  beta = Angle(p[start].x, p[start].y,
	       p[(start + 1) % n].x, p[(start + 1) % n].y,
	       p[(start + 2) % n].x, p[(start + 2) % n].y);
  gamma = beta;
  delta = beta;
  j = k = m = 1;
  Amin = 1.0e+30;
  /* find area of encasing rectangle anchored on each edge */
  for (i = 0; i < n; ++i)
    {
      /* find angle of rotation of next edge of polygon */
      if (i > 0)
	alpha += Angle(p[(start + i + n - 1) % n].x, p[(start + i + n - 1) % n].y,
		       p[(start + i) % n].x, p[(start + i) % n].y,
		       p[(start + i + 1) % n].x, p[(start + i + 1) % n].y);

      /* find a vertex on first perpendicular line of support */
      while (beta < alpha + M_PI / 2.0)
	{
	  ++j;
	  beta += Angle(p[(start + j + n - 1) % n].x, p[(start + j + n - 1) % n].y,
			p[(start + j) % n].x, p[(start + j) % n].y,
			p[(start + j + 1) % n].x, p[(start + j + 1) % n].y);
	}

      /* find a vertex on parallel line of support */
      while (gamma < alpha + M_PI)
	{
	  ++k;
	  gamma += Angle(p[(start + k + n - 1) % n].x, p[(start + k + n - 1) % n].y,
			 p[(start + k) % n].x, p[(start + k) % n].y,
			 p[(start + k + 1) % n].x, p[(start + k + 1) % n].y);
	}
      
      /* find a vertex on second perpendicular line of support */
      while (delta < alpha + 3.0 * M_PI / 2.0)
	{
	  ++m;
	  delta += Angle(p[(start + m + n - 1) % n].x, p[(start + m + n - 1) % n].y,
			 p[(start + m) % n].x, p[(start + m) % n].y,
			 p[(start + m + 1) % n].x, p[(start + m + 1) % n].y);
	}

      /* find distances between parallel and perpendicular lines of support */
      if (p[(start + i + 1) % n].x == p[(start + i) % n].x)
	{
	  d1 = fabs(p[(start + k) % n].x - p[(start + i) % n].x);
	  d2 = fabs(p[(start + m) % n].y - p[(start + j) % n].y);
	}
      else if (p[(start + i + 1) % n].y == p[(start + i) % n].y)
	{
	  d1 = fabs(p[(start + k) % n].y - p[(start + i) % n].y);
	  d2 = fabs(p[(start + m) % n].x - p[(start + j) % n].x);
	}
      else
	{
	  slope = (p[(start + i + 1) % n].y- p[(start + i) %n].y) /
	    (p[(start + i + 1) % n].x - p[(start + i) % n].x);
	  d1 = PerpDist(p[(start + i) % n].x, p[(start + i) % n].y,
			p[(start + k) % n].x, p[(start + k) % n].y,
			slope);
	  d2 = PerpDist(p[(start + j) % n].x, p[(start + j) % n].y,
			p[(start + m) % n].x, p[(start + m) % n].y,
			-1.0 / slope);
	}

      /* compute area of encasing rectangle anchored on current edge */
      A = d1 * d2;
      printf("i = %d  d1 = %f d2 = %f A = %f\n", i, d1, d2, A);
      if (i == 0 || A < Amin)
	{
	  Amin = A;
	  printf("PTi = %f %f\n", p[(start + i) % n].x, p[(start + i) % n].y);
	  printf("PTi+1 = %f %f\n", p[(start + i + 1) % n].x, p[(start + i + 1) % n].y);
	  printf("PTj = %f %f\n", p[(start + j) % n].x, p[(start + j) % n].y);
	  printf("PTk = %f %f\n", p[(start + k) % n].x, p[(start + k) % n].y);
	  printf("PTm = %f %f\n", p[(start + m) % n].x, p[(start + m) % n].y);
	  p2.x = p[(start + k) % n].x + (p[(start + i + 1) % n].x -
					 p[(start + i) % n].x);
	  p2.y = p[(start + k) % n].y + (p[(start + i + 1) % n].y -
					 p[(start + i) % n].y);
	  printf("PT2 = %f %f\n", p2.x, p2.y);

	  PointProjectionOnLine(p[(start + i) % n],
				p[(start + i + 1) % n],
				p[(start + j) % n],
				&corners[0]);
	  PointProjectionOnLine(p[(start + i) % n],
				p[(start + i + 1) % n],
				p[(start + m) % n],
				&corners[1]);
	  PointProjectionOnLine(p[(start + k) % n],
				p2,
				p[(start + m) % n],
				&corners[2]);
	  PointProjectionOnLine(p[(start + k) % n],
				p2,
				p[(start + j) % n],
				&corners[3]);
	}
    }
  return(Amin);
}

void
PointProjectionOnLine (Point p0, Point p1, Point q, Point *proj)
{
  double dx, dy;
  double u;
  
  dx = p1.x - p0.x;
  dy = p1.y - p0.y;
  u = ((q.x - p0.x) * dx + (q.y - p0.y) * dy) / (dx * dx + dy * dy);
  proj->x = p0.x + u * dx;
  proj->y = p0.y + u * dy;
}

double
PerpDist (double x1, double y1, double x2, double y2, double s)
{
  return(fabs(s * (x2 - x1) - y2 + y1) / sqrt(s*s + 1.0));
}

double
Angle (double x0, double y0,
       double x1, double y1,
       double x2, double y2)
{
  double a = atan2(y2 - y1, x2 - x1) - atan2(y1 - y0, x1 - x0);
  if (a > M_PI)
    a -= 2.0 * M_PI;
  else if (a < -M_PI)
    a += 2.0 * M_PI;
  return(a);
}
