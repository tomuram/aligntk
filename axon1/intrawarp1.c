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
#include <mpi.h>

#include "par.h"

extern int ReadImage (char *filename, unsigned char **pixels,
		      int *width, int *height,
		      int minX, int maxX, int minY, int maxY,
		      char *error);
extern int ReadBitmap (char *filename, unsigned char **bitmap,
		       int *width, int *height,
		       int minX, int maxX, int minY, int maxY,
		       char *error);

#define MAX_LEVELS	32       /* max image size (w or h) is 2^MAX_LEVELS */
#define MAP(map,x,y)	map[((y)+1)*((size_t) mpwp) + (x) + 1]
#define MAP1(map,x,y)	map[((y)+1)*((size_t) mpw1p) + (x) + 1]
#define MAPINP(map,x,y)	map[((y)+1)*((size_t) mpwi) + (x) + 1]
#define IIMAGE(i,x,y)   i[(y)*((size_t) iw) + (x)]
#define RIMAGE(i,x,y)   i[(y)*((size_t) rw) + (x)]
#define REFIMAGE(i,x,y) i[(y)*((size_t) refw) + (x)]
#define LINE_LENGTH	255

typedef struct Context {
  /* NOTE: any new fields added to this struct should
     be also added to PackContext and UnpackContext */ 
  char type;
  char imageBasename[PATH_MAX];
  char referenceBasename[PATH_MAX];

  char cptsName[PATH_MAX];
  char inputMapName[PATH_MAX];

  char outputMapBasename[PATH_MAX];
  char outputImagesName[PATH_MAX];
  char outputCorrelationName[PATH_MAX];

  int outputLevel;
  int minResolution;
  int depth;
  double distortion;
  double correspondence;
  double correspondenceThreshold;
  double quality;

  int correlationHalfWidth;
  int correlationHalfHeight;
  int writeAllMaps;
  double afterTime;
} Context;

typedef struct Pair {
  char *imageName;
  int imageMinX, imageMaxX, imageMinY, imageMaxY;
  char *refName;
  int refMinX, refMaxX, refMinY, refMaxY;
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
  char *message;
} Result;

typedef struct CPoint {
  int ix, iy;
  int rx, ry;
  double energy;
  double newEnergy;
} CPoint;

/* GLOBAL VARIABLES FOR MASTER */
int nResults = 0;
Result* results = 0;

/* GLOBAL VARIABLES FOR MASTER & WORKER */
Context c;
Task t;
Result r;
FILE *logFile = NULL;

/* GLOBAL VARIABLES FOR MASTER & WORKER */
Task t;

/* GLOBAL VARIABLES FOR WORKER */
unsigned char *ref_in;
unsigned char *image_in;

float* mapping;
int nLevels;
float* refs[MAX_LEVELS];       /* reference image to warp to (at various
			         resolutions) */
float* images[MAX_LEVELS];     /* image to be warped (at various resolutions) */
float* warped[MAX_LEVELS];     /* warped reference image (at various resol.) */

int imageWidth[MAX_LEVELS], imageHeight[MAX_LEVELS];
int referenceWidth[MAX_LEVELS], referenceHeight[MAX_LEVELS];

float* xmaps[MAX_LEVELS];
float* ymaps[MAX_LEVELS];
int inputMapFactor;
float *inputMapX, *inputMapY, *inputMapC;
int mpwi, mphi;
int mpwip, mphip;
int mpwc, mphc;
int mpwcp, mphcp;


#if 0
float* isig[MAX_LEVELS];
float* rsig[MAX_LEVELS];
#endif

int nCpts = 0;
CPoint *cpts = 0;

#define DIR_HASH_SIZE	8192
char *dirHash[DIR_HASH_SIZE];

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
void WriteMap (char *outputName, int level);
void Compute(char *outputName);
int Compare (const void *x, const void *y);
int SortBySlice (const void *x, const void *y);
int SortByEnergy (const void *x, const void *y);
int ParseRange (char *s, int *pos, int *minValue, int *maxValue);
int ParseValue (char *s, int *pos, int *value);
size_t CountBits (unsigned char *p, size_t n);
int CreateDirectories (char *fn);
void Error (char *fmt, ...);
void Log (char *fmt, ...);

int
main (int argc, char **argv, char **envp)
{
  par_process(argc, argv, envp,
              (void (*)()) MasterTask, MasterResult,
              WorkerContext, WorkerTask, NULL,
              PackContext, UnpackContext,
              PackTask, UnpackTask,
              PackResult, UnpackResult);
}

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
  char outputDirName[PATH_MAX];
  char tc;
  unsigned int w, h;
  int m;
  int xRes, yRes;
  int minS, maxS;
  int forward;
  char afterFileName[PATH_MAX];
  struct stat statBuf;
  int pn;
  char imgn[PATH_MAX];
  char refn[PATH_MAX];
  char pairn[PATH_MAX];
  int imgMinX, imgMaxX, imgMinY, imgMaxY;
  int refMinX, refMaxX, refMinY, refMaxY;
  unsigned int iw, ih;
  int imageSlice, refSlice;

  error = 0;
  c.type = '\0';
  c.imageBasename[0] = '\0';
  c.referenceBasename[0] = '\0';
  c.inputMapName[0] = '\0';
  c.cptsName[0] = '\0';
  c.outputMapBasename[0] = '\0';
  c.outputImagesName[0] = '\0';
  c.outputCorrelationName[0] = '\0';
  c.outputLevel = 0;
  c.minResolution = 1;
  c.distortion = 1.0;
  c.correspondence = 1.0;
  c.correspondenceThreshold = 0.0;
  c.quality = 10.0;
  c.correlationHalfWidth = -1;
  c.correlationHalfHeight = -1;
  c.writeAllMaps = 0;
  c.depth = 0;
  r.pair.imageName = NULL;
  r.pair.refName = NULL;
  r.pair.pairName = NULL;
  r.message = (char *) malloc(PATH_MAX + 1024);
  forward = 1;
  pairsFile[0] = '\0';
  nPairs = 0;
  pairs = 0;
  minSlice = 0;
  maxSlice = 1000000000;
  afterFileName[0] = '\0';
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
    else if (strcmp(argv[i], "-output") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.outputMapBasename, argv[i]);
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
    else if (strcmp(argv[i], "-pairs") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(pairsFile, argv[i]);
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
    else if (strcmp(argv[i], "-correlation_half_height") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%d", &c.correlationHalfHeight) != 1)
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
    else if (strcmp(argv[i], "-after") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(afterFileName, argv[i]);
      }
    else if (strcmp(argv[i], "-reference") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.referenceBasename, argv[i]);
      }
    else if (strcmp(argv[i], "-tif") == 0)
      c.type = 't';
    else error = 1;

  if (error)
    {
      fprintf(stderr, "Usage: intrawarp -input file_prefix -output file_prefix\n");
      fprintf(stderr, "              [-skip skip_file] [-slice slice_range]\n");
      fprintf(stderr, "              [-resolution HORIZxVERT]\n");
      fprintf(stderr, "              [-distortion distortion_weight]\n");
      fprintf(stderr, "              [-correspondence correspondence_weight]\n");
      fprintf(stderr, "              [-quality quality_factor]\n");
      fprintf(stderr, "              [-min_res minimum_resolution_in_pixels]\n");
      fprintf(stderr, "              [-depth delta_depth]\n");
      fprintf(stderr, "              [-all_maps]\n");
      fprintf(stderr, "              [-after file]\n");
      fprintf(stderr, "              [-pairs <pair_file>]\n");
      fprintf(stderr, "              [-reference <reference_prefix>]\n");
      fprintf(stderr, "              [-input_map <input_map_prefix>]\n");
      fprintf(stderr, "   where ranges are expressed as: integer\n");
      fprintf(stderr, "                              or: integer-integer\n");
      fprintf(stderr, "                              or: integer-\n");
      fprintf(stderr, "                              or: -integer\n");
      exit(1);
    }

  Log("MASTER starting on node %d\n", par_instance());

  /* check that at least minimal parameters were supplied */
  if (c.imageBasename[0] == '\0' || c.outputMapBasename[0] == '\0' ||
      pairsFile[0] == '\0')
    Error("-input, -output, and -pairs parameters must be specified.\n");

  if (afterFileName[0] != '\0')
    {
      if (stat(afterFileName, &statBuf) != 0)
	Error("Cannot access the -after file: %s\n", afterFileName);
      c.afterTime = (double) (statBuf.st_mtime);
    }
  else
    c.afterTime = 0.0;

  f = fopen(pairsFile, "r");
  if (f == NULL)
    Error("Could not open pairs file %s\n", pairsFile);
  while (fscanf(f, "%s %d %d %d %d %s %d %d %d %d %s",
		imgn, &imgMinX, &imgMaxX, &imgMinY, &imgMaxY,
		refn, &refMinX, &refMaxX, &refMinY, &refMaxY,
		pairn) == 11)
    {
      if ((nPairs & 1023) == 0)
	pairs = (Pair *) realloc(pairs, (nPairs + 1024) * sizeof(Pair));
      pairs[nPairs].imageName = (char *) malloc(strlen(imgn)+1);
      strcpy(pairs[nPairs].imageName, imgn);
      pairs[nPairs].imageMinX = imgMinX;
      pairs[nPairs].imageMaxX = imgMaxX;
      pairs[nPairs].imageMinY = imgMinY;
      pairs[nPairs].imageMaxY = imgMaxY;
      pairs[nPairs].refName = (char *) malloc(strlen(refn)+1);
      strcpy(pairs[nPairs].refName, refn);
      pairs[nPairs].refMinX = refMinX;
      pairs[nPairs].refMaxX = refMaxX;
      pairs[nPairs].refMinY = refMinY;
      pairs[nPairs].refMaxY = refMaxY;
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
  t.pair.imageMinX = t.pair.imageMaxX = t.pair.imageMinY = t.pair.imageMaxY = -1;
  t.pair.refMinX = t.pair.refMaxX = t.pair.refMinY = t.pair.refMaxY = -1;
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

  sprintf(fn, "%squality.out", c.outputMapBasename);
  f = fopen(fn, "w");
  qsort(results, nResults, sizeof(Result), SortBySlice);
  fprintf(f, "Sorted by slice:\n");
  fprintf(f, "         IMAGE     REFERENCE   CORRELATION    DISTORTION    CORRESPOND        ENERGY\n");
  for (i = 0; i < nResults; ++i)
    fprintf(f, "%s  %s  %12.6f  %12.6f  %12.6f  %12.6f\n",
	    results[i].pair.imageName,
	    results[i].pair.refName,
	    results[i].correlation,
	    results[i].distortion,
	    results[i].correspondence,
	    results[i].distortion * c.distortion - results[i].correlation +
	    results[i].correspondence * c.correspondence);
  fprintf(f, "\n\nSorted by energy:\n");
  fprintf(f, "         IMAGE     REFERENCE   CORRELATION    DISTORTION    CORRESPOND        ENERGY\n");
  qsort(results, nResults, sizeof(Result), SortByEnergy);
  for (i = 0; i < nResults; ++i)
    fprintf(f, "%s  %s  %12.6f  %12.6f  %12.6f  %12.6f\n",
	    results[i].pair.imageName,
	    results[i].pair.refName,
	    results[i].correlation,
	    results[i].distortion,
	    results[i].correspondence,
	    results[i].distortion * c.distortion - results[i].correlation +
	    results[i].correspondence * c.correspondence);
  fclose(f);

  printf(" %d\nAll slices completed.\n", nResults);
  exit(0);
}

void
MasterResult ()
{
  if (r.message[0] != '\0')
    Error("\nThe following error was encountered by one of the worker processes:\n");

  if (!r.updated)
    return;

  results = (Result*) realloc(results, (nResults + 1) * sizeof(Result));
  results[nResults] = r;
  results[nResults].message = NULL;
  r.pair.imageName = NULL;
  r.pair.refName = NULL;
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
  return(strcmp(rx->pair.imageName, ry->pair.imageName));
}

int
SortByEnergy (const void *x, const void *y)
{
  Result *rx, *ry;
  double ex, ey;
  rx = (Result *) x;
  ry = (Result *) y;
  ex = rx->distortion * c.distortion - rx->correlation + rx->correspondence * c.correspondence;
  ey = ry->distortion * c.distortion - ry->correlation + rx->correspondence * c.correspondence;
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
  r.pair.imageName = NULL;
  r.pair.refName = NULL;
  r.pair.pairName = NULL;
  r.message = (char *) malloc(PATH_MAX + 1024);
  r.message[0] = '\0';
  t.pair.imageName = NULL;
  t.pair.refName = NULL;
  t.pair.pairName = NULL;
}

void
WorkerTask ()
{
  FILE *f;
  unsigned int w, h;
  unsigned int iw, ih;
  unsigned int rw, rh;
  char image_filename[256];
  char reference_filename[256];
  int i;
  int n;
  float *xmap;
  float *ymap;
  int mpwp;
  char fn[PATH_MAX];
  char tc;
  int m;
  int bp;  /* bytes per pixel */
  int mx, my;
  int x, y;
  int irx, iry;
  double rx, ry;
  double rrx, rry;
  float *warp;
  unsigned char *warpout;
  int level;
  int imw, imh;
  int rmw, rmh;
  int mpw, mph;
  int imx, imy;
  int refSlice;
  int refx, refy;
  char line[LINE_LENGTH+1];
  struct stat sb;
  int update;
  size_t imagePixels;
  size_t referencePixels;
  size_t ii;
  int mapHeader[3];
  float *mapX, *mapY, *mapC;
  char imageName[PATH_MAX];
  char refName[PATH_MAX];
  char cptsName[PATH_MAX], inputMapName[PATH_MAX];
  char outputName[PATH_MAX], outputImagesName[PATH_MAX], outputCorrelationName[PATH_MAX];
  int cw, ch;
  int factor;
  double si, si2, sr, sr2, sir;
  int nPoints;
  int dx, dy;
  int px, py;
  unsigned char *correlation;
  double xv, yv;
  int ixv, iyv;
  double rx00, ry00, rx01, ry01, rx10, ry10, rx11, ry11;
  double r00, r01, r10, r11;
  double rv, iv;
  double mi, mr;
  double denom;
  float *cimage, *cref;

  Log("WORKER starting on node %d\n", par_instance());

  /* construct filenames */
  sprintf(imageName, "%s%s.%s", c.imageBasename, t.pair.imageName,
	  c.type == 't' ? "tif" : "pgm");
  if (c.referenceBasename[0] != '\0')
    sprintf(refName, "%s%s.%s", c.referenceBasename, t.pair.refName,
	    c.type == 't' ? "tif" : "pgm");
  else
    sprintf(refName, "%s%s.%s", c.imageBasename, t.pair.refName,
	    c.type == 't' ? "tif" : "pgm");
  if (c.cptsName[0] != '\0')
    sprintf(cptsName, "%s%s.cpts", c.cptsName, t.pair.pairName);
  else
    cptsName[0] = '\0';
  if (c.inputMapName[0] != '\0')
    sprintf(inputMapName, "%s%s.map", c.inputMapName, t.pair.pairName);
  else
    inputMapName[0] = '\0';
  sprintf(outputName, "%s%s", c.outputMapBasename, t.pair.pairName);
  if (c.outputImagesName[0] != '\0')
    sprintf(outputImagesName, "%s%s.pgm", c.outputImagesName, t.pair.pairName);
  else
    outputImagesName[0] = '\0';
  if (c.outputCorrelationName[0] != '\0')
    sprintf(outputCorrelationName, "%s%s.pgm", c.outputCorrelationName, t.pair.pairName);
  else
    outputCorrelationName[0] = '\0';

  /* check if we can skip this task because of the -after option */
  update = 0;
  if (c.afterTime == 0.0)
    update = 1;
  if (!update)
    {
      if (stat(imageName, &sb) == 0 &&
	  (double) (sb.st_mtime) > c.afterTime)
	update = 1;
    }
  if (!update)
    {
      if (stat(refName, &sb) == 0 &&
	  (double) (sb.st_mtime) > c.afterTime)
	update = 1;
    }
  if (!update && cptsName[0] != '\0')
    {
      if (stat(cptsName, &sb) == 0 &&
	  (double) (sb.st_mtime) > c.afterTime)
	update = 1;
    }
  if (!update && inputMapName[0] != '\0')
    {
      if (stat(inputMapName, &sb) == 0 &&
	  (double) (sb.st_mtime) > c.afterTime)
	update = 1;
    }
  if (!update)
    {
      memcpy(&(r.pair), &(t.pair), sizeof(Pair));
      r.updated = 0;
      r.correlation = 0.0;
      r.distortion = 0.0;
      r.correspondence = 0.0;
      r.message[0] = '\0';
      return;
    }

  Log("STARTING WORKER TASK\n");

  inputMapX = NULL;
  inputMapY = NULL;
  inputMapC = NULL;

  Log("WORKER reading ref %s\n", refName);

  /* get the reference image */
  if (!ReadImage(refName, &ref_in,
		 &referenceWidth[0], &referenceHeight[0],
		 t.pair.refMinX, t.pair.refMaxX,
		 t.pair.refMinY, t.pair.refMaxY,
		 r.message))
    {
      sprintf(r.message, "Could not read reference image %s\n", refName);
      return;
    }
  Log("WORKER READ REFERENCE: %s\n", refName);

  /* read in the image to warp */
  Log("WORKER READING IMAGE: %s\n", imageName);
  if (!ReadImage(imageName, &image_in,
		 &imageWidth[0], &imageHeight[0],
		 t.pair.imageMinX, t.pair.imageMaxX,
		 t.pair.imageMinY, t.pair.imageMaxY,
		 r.message))
    {
      sprintf(r.message, "Could not read image %s\n", imageName);
      return;
    }
  Log("WORKER READ IMAGE: %s\n", imageName);

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
	      sscanf(line, "%d%d%d%d", &imx, &imy, &refx, &refy);
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


  imagePixels = ((size_t) imageWidth[0]) * imageHeight[0];
  referencePixels = ((size_t) referenceWidth[0]) * referenceHeight[0];
  images[0] = (float *) malloc(imagePixels * sizeof(float));
  Log("images[0] = %llx  imp = %d\n", (long long) images[0], (int) imagePixels);
  for (ii = 0; ii < imagePixels; ++ii)
    images[0][ii] = image_in[ii];
  free(image_in);
  refs[0] = (float *) malloc(referencePixels * sizeof(float));
  for (ii = 0; ii < referencePixels; ++ii)
    refs[0][ii] = ref_in[ii];
  free(ref_in);

  /* get the input map, if present */
  if (inputMapName[0] != '\0')
    {
      f = fopen(inputMapName, "r");
      if (f != NULL)
	{
	  if (fread(mapHeader, sizeof(int), 3, f) != 3)
	    Error("Could not read map header from file %s", fn);
	  mpwi = mapHeader[0];
	  mphi = mapHeader[1];
	  mpwip = mpwi + 2;
	  mphip = mphi + 2;
	  inputMapX = (float *) malloc(mpwip * mphip * sizeof(float));
	  inputMapY = (float *) malloc(mpwip * mphip * sizeof(float));
	  inputMapC = (float *) malloc(mpwip * mphip * sizeof(float));
	  if (fread(inputMapX, sizeof(float), mpwip * mphip, f) != mpwip * mphip ||
	      fread(inputMapY, sizeof(float), mpwip * mphip, f) != mpwip * mphip ||
	      fread(inputMapC, sizeof(float), mpwip * mphip, f) != mpwip * mphip)
	    Error("Could not read map from file %s\n", fn);
	  fclose(f);
	  if (mpwi > imageWidth[0])
	    Error("Input map resolution cannot be greater than image resolution\n");
	  inputMapFactor = 1 << ((int) floor(log(((float) imageWidth[0]) / mpwi) / log(2.0) + 0.5));
	  Log("Read map %s of size %d x %d;  inputMapFactor = %d\n", fn, mpwi, mphi, inputMapFactor);
	}
      else
	{
	  Log("No input map found for slice %s\n", t.pair.imageName);
	  inputMapFactor = 0;
	}
    }
  else
    inputMapFactor = 0;

  if (!Init())
    {
      Log("Init was unsuccessful.\n");
      return;
    }
  Compute(outputName);

  if (outputImagesName[0] != '\0')
    {
      /* write out the warped image */
      warp = warped[c.outputLevel];
      iw = imageWidth[c.outputLevel];
      ih = imageHeight[c.outputLevel];
      imagePixels = ((size_t) iw) * ih;
      warpout = (unsigned char*) malloc(imagePixels);
      f = fopen(outputImagesName, "w");
      if (f == NULL)
	{
	  sprintf(r.message, "Could not open image file %s for writing.\n", fn);
	  return;
	}
      fprintf(f, "P5\n%d %d\n255\n", iw, ih);
      for (y = 0; y < ih; ++y)
	for (x = 0; x < iw; ++x)
	  warpout[y*((size_t) iw) + x] = (unsigned char) warp[y * ((size_t) iw) + x];
      if (fwrite(warpout, 1, imagePixels, f) != imagePixels)
	{
	  fclose(f);
	  sprintf(r.message, "Could not write to image file %s\n", fn);
	  return;
	}
      fclose(f);
      free(warpout);
    }

  if (outputCorrelationName[0] != '\0')
    {
      /* write out the correlation image */
      iw = imageWidth[0];
      ih = imageHeight[0];
      imagePixels = ((size_t) iw) * ih;
      rw = referenceWidth[0];
      rh = referenceHeight[0];
      cimage = images[0];
      cref = refs[0];
      correlation = (unsigned char*) malloc(iw*ih);
      f = fopen(outputCorrelationName, "w");
      if (f == NULL)
	{
	  sprintf(r.message, "Could not open image file %s for writing.\n", fn);
	  return;
	}
      fprintf(f, "P5\n%d %d\n255\n", iw, ih);
      if (c.correlationHalfWidth >= 0)
	if (c.correlationHalfHeight >= 0)
	  {
	    cw = c.correlationHalfWidth;
	    ch = c.correlationHalfHeight;
	  }
	else
	  cw = ch = c.correlationHalfWidth;
      else
	if (c.correlationHalfHeight >= 0)
	  cw = ch = c.correlationHalfHeight;
	else
	  cw = ch = 7;
      factor = 1 << c.outputLevel;
      for (y = 0; y < ih; ++y)
	{
	  si = 0.0;
	  si2 = 0.0;
	  sr = 0.0;
	  sr2 = 0.0;
	  sir = 0.0;
	  nPoints = 0;
	  for (x = -cw; x < iw; ++x)
	    {
	      for (i = -1; i <= 1; i += 2)
		{
		  px = x + i * cw;
		  if (px < 0 || px >= iw)
		    continue;
		  for (dy = -ch; dy <= ch; ++dy)
		    {
		      py = y + dy;
		      if (py < 0 || py >= ih)
			continue;
		    
		      /* use bilinear interpolation to find value */
		      xv = (2.0 * px - factor + 1) / (2.0 * factor);
		      yv = (2.0 * py - factor + 1) / (2.0 * factor);
		      ixv = ((int) (xv + 2.0)) - 2;
		      iyv = ((int) (yv + 2.0)) - 2;
		      rrx = xv - ixv;
		      rry = yv - iyv;
		      if (ixv < -1 || ixv > mpw || ixv == mpw && rrx > 0.0 ||
			  iyv < -1 || iyv > mph || iyv == mph && rry > 0.0)
			continue;
		      rx00 = MAP(xmap, ixv, iyv);
		      ry00 = MAP(ymap, ixv, iyv);
		      if (iyv < mph)
			{
			  rx01 = MAP(xmap, ixv, iyv + 1);
			  ry01 = MAP(ymap, ixv, iyv + 1);
			}
		      else
			{
			  rx01 = 0.0;
			  ry01 = 0.0;
			}
		      if (ixv < mpw)
			{
			  rx10 = MAP(xmap, ixv + 1, iyv);
			  ry10 = MAP(ymap, ixv + 1, iyv);
			}
		      else
			{
			  rx10 = 0.0;
			  ry10 = 0.0;
			}
		      if (ixv < mpw && iyv < mph)
			{
			  rx11 = MAP(xmap, ixv + 1, iyv + 1);
			  ry11 = MAP(ymap, ixv + 1, iyv + 1);
			}
		      else
			{
			  rx11 = 0.0;
			  ry11 = 0.0;
			}

		      rx = rx00 * (rrx - 1.0) * (rry - 1.0)
			- rx10 * rrx * (rry - 1.0) 
			- rx01 * (rrx - 1.0) * rry
			+ rx11 * rrx * rry;
		      ry = ry00 * (rrx - 1.0) * (rry - 1.0)
			- ry10 * rrx * (rry - 1.0) 
			- ry01 * (rrx - 1.0) * rry
			+ ry11 * rrx * rry;
		      rx = factor * rx + (factor - 1.0) / 2.0;
		      ry = factor * ry + (factor - 1.0) / 2.0;
		      irx = ((int) (rx + 1.0)) - 1;
		      iry = ((int) (ry + 1.0)) - 1;
		      rrx = rx - irx;
		      rry = ry - iry;
		      if (irx < 0 || irx >= rw || irx == rw-1 && rrx > 0.0 ||
			  iry < 0 || iry >= rh || iry == rh-1 && rry > 0.0)
			continue;

		      r00 = RIMAGE(cref, irx, iry);
		      if (iry < rh-1)
			r01 = RIMAGE(cref, irx, iry + 1);
		      else
			r01 = 0.0;
		      if (irx < rw-1)
			{
			  r10 = RIMAGE(cref, irx + 1, iry);
			  if (iry < rh-1)
			    r11 = RIMAGE(cref, irx + 1, iry + 1);
			  else
			    r11 = 0.0;
			}
		      else
			{
			  r10 = 0.0;
			  r11 = 0.0;
			}
		      rv = r00 * (rrx - 1.0) * (rry - 1.0)
			- r10 * rrx * (rry - 1.0) 
			- r01 * (rrx - 1.0) * rry
			+ r11 * rrx * rry;
		      iv = IIMAGE(cimage, x, y);
		      nPoints += i;
		      si += i * iv;
		      si2 += i * iv * iv;
		      sr += i * rv;
		      sr2 += i * rv * rv;
		      sir += i * iv * rv;
		    }
		}
	      if (x < 0)
		continue;
	      if (nPoints <= 0)
		{
		  correlation[y*w+x] = 0;
		  continue;
		}
	      mi = si / nPoints;
	      mr = sr / nPoints;
	      denom = (si2 - 2.0 * mi * si + nPoints * mi * mi) *
		(sr2 - 2.0 * mr * sr + nPoints * mr * mr);
	      if (denom < 0.001)
		{
		  correlation[y*w+x] = 0;
		  continue;
		}
	      correlation[y*iw+x] = (int) floor(((sir - mi * sr - mr * si + nPoints * mi * mr) /
						sqrt(denom)) * 255.999);
	    }
	  }
      if (fwrite(correlation, 1, imagePixels, f) != imagePixels)
	{
	  fclose(f);
	  sprintf(r.message, "Could not write to image file %s\n", fn);
	  return;
	}
      fclose(f);
      free(correlation);
    }

  for (level = 0; level < nLevels; ++level)
    {
      free(images[level]);
      free(refs[level]);
      free(xmaps[level]);
      free(ymaps[level]);
    }
}

int
Init ()
{
  int level;
  unsigned int iw, ih;
  unsigned int rw, rh;
  float *image, *ref;
  size_t imagePixels, referencePixels;
  size_t ii;
  int i;
  int x, y;
  int b;
  int dx, dy;
  unsigned char *iCount, *rCount;
  size_t impw, imph;

  /* make the full resolution warped image */
  imagePixels = ((size_t) imageWidth[0]) * imageHeight[0];
  if (c.outputImagesName[0] != '\0')
    warped[0] = (float *) malloc(imagePixels * sizeof(float));
  Log("iw = %d ih = %d pixels = %d\n", imageWidth[0], imageHeight[0], imagePixels);
  Log("rw = %d rh = %d\n", referenceWidth[0], referenceHeight[0]);

#if 0
  /* make the full resolution image significance array by taking the standard
     deviation of the pixel values over a circular region centered on each pixel */
  isig = (float *) malloc(imagePixels * sizeof(float));
#endif    

  iw = imageWidth[0];
  ih = imageHeight[0];
  rw = referenceWidth[0];
  rh = referenceHeight[0];

  /* create hierarchical representations of images */
  for (level = 1;
       imageWidth[level-1] > 1 || imageHeight[level-1] > 1 ||
	 referenceWidth[level-1] > 1 || referenceHeight[level-1] > 1;
       ++level)
    {
      Log("creating level %d\n", level);
      iw = imageWidth[level] = (imageWidth[level-1] + 1) >> 1;
      ih = imageHeight[level] = (imageHeight[level-1] + 1) >> 1;
      rw = referenceWidth[level] = (referenceWidth[level-1] + 1) >> 1;
      rh = referenceHeight[level] = (referenceHeight[level-1] + 1) >> 1;
      imagePixels = ((size_t) imageWidth[level]) * imageHeight[level];
      referencePixels = ((size_t) referenceWidth[level]) * referenceHeight[level];
      impw = iw + 2;
      imph = ih + 2;

      Log("imagePixels = %d referencePixels = %d\n", imagePixels,
	  referencePixels);
      images[level] = (float*) malloc(imagePixels * sizeof(float));
      refs[level] = (float*) malloc(referencePixels * sizeof(float));
      if (images[level] == NULL || refs[level] == NULL)
	{
	  sprintf(r.message, "Could not allocate image arrays\n");
	  return(0);
	}
      if (c.outputImagesName[0] != '\0')
	{
	  warped[level] = (float*) malloc(imagePixels * sizeof(float));
	  if (warped[level] == NULL)
	    {
	      sprintf(r.message, "Could not allocate warped array\n");
	      return(0);
	    }
	}

      image = images[level];
      ref = refs[level];
      for (ii = 0; ii < imagePixels; ++ii)
	image[ii] = 0.0;
      for (ii = 0; ii < referencePixels; ++ii)
	ref[ii] = 0.0;

      ii = 0;
      for (y = 0; y < imageHeight[level-1]; ++y)
	for (x = 0; x < imageWidth[level-1]; ++x)
	  {
	    image[(y >> 1) * ((size_t) iw) + (x >> 1)] += images[level-1][ii];
	    ++ii;
	  }

      ii = 0;
      for (y = 0; y < referenceHeight[level-1]; ++y)
	for (x = 0; x < referenceWidth[level-1]; ++x)
	  {
	    ref[(y >> 1) * ((size_t) rw) + (x >> 1)] += refs[level-1][ii];
	    ++ii;
	  }

      for (ii = 0; ii < imagePixels; ++ii)
	image[ii] *= 1.0 / 4.0;
      for (ii = 0; ii < referencePixels; ++ii)
	ref[ii] *= 1.0 / 4.0;
    }
  nLevels = level;
  Log("nLevels = %d\n", nLevels);

  if (c.outputLevel < 0 || c.outputLevel >= nLevels)
    Error("outputLevel (%d) is invalid for image (nLevels = %d)\n",
	  c.outputLevel, nLevels);

  Log("Constructed the multi-resolution image set.\n");
  return(1);
}
  
void
Compute (char *outputName)
{
  int level;
  float *xprop, *yprop;
  unsigned char* changed;
  int i;
  size_t ii;
  int x, y;
  double x00, x01, x10, x11;
  double y00, y01, y10, y11;
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
  double cx, cy;
  int xdir, ydir;
  int ix, iy;
  double dxp, dyp;
  double mrd;
  double nx, ny;
  double ode, nde;
  double newDistortion, newCorrelation, newCorrespondence;;
  double newE;
  double aCorrelation, aDistortion, aCorrespondence;
  int upd0, upd1, upd2, upd3;
  double iv, rv;
  double distortion, correlation, correspondence;
  double e;
  double corr;
  int move;
  int accept;
  double prob;
  double ran;
  float *image;
  float *ref;
  float *warp;
  int mpwp, mphp;
  int mpw1, mph1;
  int mpw1p, mph1p;
  float *xmap, *ymap;
  float *xmap1, *ymap1;
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
  size_t n;
  size_t updcnt;
  double x0, x1, x2, x3;
  double y0, y1, y2, y3;
  double area;
  double nx0, nx1, nx2, nx3;
  double ny0, ny1, ny2, ny3;
  double narea;
  size_t rawMoveCount, moveCount, acceptedMoveCount;
  double l0, l1, l2, l3;
  double nl0, nl1, nl2, nl3;
  double goalCoverage;
  double actualCoverage;
  int cLevel;
  int nb;
  float *cimage;
  float *cref;
  double rx00, rx01, rx10, rx11, ry00, ry01, ry10, ry11;
  size_t nPoints;
  size_t requiredPoints;
  int factor;
  int sx, sy, ex, ey;
  int mpw, mph;
  int ixv, iyv;
  double xv, yv;
  double ea;
  double check_sir, check_si2, check_sr2, check_si, check_sr;
  double check_sum;
  size_t check_points;
  size_t check_effectivePoints;
  double check_mi, check_mr, check_denom, check_r;
  double minOverlap;
  double sqrtMaxMoveRadius, sqrtMaxMoveOffset;
  double kFactor;
  double lFactor;
  double distance;
  double ce;
  double ces;
  size_t mSize;
  double cth, sth;
  double thetaRange;
  size_t statSqrtRadius[21];
  size_t statSqrtOffset[21];
  size_t statTheta[21];
  size_t statX[21];
  size_t statY[21];
  size_t nShiftMoves;
  int ix0, iy0;
  double rc;
  double rc00, rc01, rc10, rc11;
  unsigned int refw, refh;
  int multiplier;
  int startLevel;
  char fn[PATH_MAX];
  FILE *f;
  double logMinRadius, logMaxRadius;
  double minDist2;
  double d2;

  minOverlap = 0.25;
  kFactor = 1.0 / sqrt((double) (((size_t) referenceWidth[0]) * referenceWidth[0] +
	                         ((size_t) referenceHeight[0]) * referenceHeight[0]));
  Log("Warping started with %d total resolution levels\n", nLevels);

  if (inputMapFactor != 0)
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
      image = images[level];
      ref = refs[level];
      mpw = imageWidth[level];
      mph = imageHeight[level];
      mpwp = mpw + 2;
      mphp = mph + 2;
      refw = referenceWidth[level];
      refh = referenceHeight[level];

      cLevel = level - c.depth;
      if (cLevel < 0)
	cLevel = 0;
      while (cLevel > 0 &&
	     (imageWidth[cLevel] < c.minResolution ||
	      imageHeight[cLevel] < c.minResolution))
	--cLevel;
      nb = level - cLevel;
      factor = (1 << nb);
      ih = imageHeight[cLevel];
      iw = imageWidth[cLevel];
      rh = referenceHeight[cLevel];
      rw = referenceWidth[cLevel];
      cimage = images[cLevel];
      cref = refs[cLevel];
      lFactor = (1 << level);

      mSize = ((size_t) mpwp) * mphp;
      xmaps[level] = (float*) malloc(mSize * sizeof(float));
      ymaps[level] = (float*) malloc(mSize * sizeof(float));
      xmap = xmaps[level];
      ymap = ymaps[level];

      if (level == startLevel)
	if (inputMapFactor != 0)
	  {
	    /* use the provided input map as the initial map */
	    for (y = -1; y <= mph; ++y)
	      for (x = -1; x <= mpw; ++x)
		{
		  MAP(xmap, x, y) = MAPINP(inputMapX, x, y);
		  MAP(ymap, x, y) = MAPINP(inputMapY, x, y);
		}
	  }
	else
	  {
	    /* use the identity map */
	    for (i = 0; i < 9; ++i)
	      {
		xmap[i] = (i % 3) - 1;
		ymap[i] = (i / 3) - 1;
	      }
	    Log("Setting initial map to the identity map (mpw=%d mph=%d)\n",
		mpw, mph);
	  }
      else
	{
	  /* create initial map from next higher level map */
	  mpw1 = imageWidth[level+1];
	  mph1 = imageHeight[level+1];
	  mpw1p = mpw1 + 2;
	  mph1p = mph1 + 2;
	  xmap1 = xmaps[level+1];
	  ymap1 = ymaps[level+1];

	  /* interpolate to construct the current level map */
	  for (y = -1; y <= mph; ++y)
	    for (x = -1; x <= mpw; ++x)
	      {
		x00 = MAP1(xmap1, ((x + 1) >> 1) - 1, ((y + 1) >> 1) - 1);
		x10 = MAP1(xmap1, ((x + 1) >> 1), ((y + 1) >> 1) - 1);
		x01 = MAP1(xmap1, ((x + 1) >> 1) - 1, ((y + 1) >> 1));
		x11 = MAP1(xmap1, ((x + 1) >> 1), ((y + 1) >> 1));
		xd = (3.0 - 2.0 * (x & 1)) / 4.0;
		yd = (3.0 - 2.0 * (y & 1)) / 4.0;

		xi = x00 * (xd - 1.0) * (yd - 1.0)
		  - x10 * xd * (yd - 1.0)
		  - x01 * (xd - 1.0) * yd
		  + x11 * xd * yd;
		MAP(xmap, x, y) = 2.0 * xi + 1.0 / 2.0;

		y00 = MAP1(ymap1, ((x + 1) >> 1) - 1, ((y + 1) >> 1) - 1);
		y10 = MAP1(ymap1, ((x + 1) >> 1), ((y + 1) >> 1) - 1);
		y01 = MAP1(ymap1, ((x + 1) >> 1) - 1, ((y + 1) >> 1));
		y11 = MAP1(ymap1, ((x + 1) >> 1), ((y + 1) >> 1));

		yi = y00 * (xd - 1.0) * (yd - 1.0)
		  - y10 * xd * (yd - 1.0)
		  - y01 * (xd - 1.0) * yd
		  + y11 * xd * yd;
		MAP(ymap, x, y) = 2.0 * yi + 1.0 / 2.0;
	      }
	  Log("Interpolated level %d map (mpw=%d mph=%d + boundary) to obtain level %d map (mpw=%d mph=%d)\n",
	      level+1, mpw1, mph1, level, mpw, mph);
	}

      /* calculate energy of mapping */
      distortion = 0.0;
      correlation = 0.0;
      correspondence = 0.0;
      e = 0.0;

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
	    /* use bilinear interpolation to find value */
	    xv = (2.0 * x - factor + 1) / (2.0 * factor);
	    yv = (2.0 * y - factor + 1) / (2.0 * factor);
	    ixv = ((int) (xv + 2.0)) - 2;
	    iyv = ((int) (yv + 2.0)) - 2;
	    rrx = xv - ixv;
	    rry = yv - iyv;
	    if (ixv < -1 || ixv > mpw || ixv == mpw && rrx > 0.0 ||
		iyv < -1 || iyv > mph || iyv == mph && rry > 0.0)
	      continue;
	    rx00 = MAP(xmap, ixv, iyv);
	    ry00 = MAP(ymap, ixv, iyv);
	    if (iyv < mph)
	      {
		rx01 = MAP(xmap, ixv, iyv + 1);
		ry01 = MAP(ymap, ixv, iyv + 1);
	      }
	    else
	      {
		rx01 = 0.0;
		ry01 = 0.0;
	      }
	    if (ixv < mpw)
	      {
		rx10 = MAP(xmap, ixv + 1, iyv);
		ry10 = MAP(ymap, ixv + 1, iyv);
	      }
	    else
	      {
		rx10 = 0.0;
		ry10 = 0.0;
	      }
	    if (ixv < mpw && iyv < mph)
	      {
		rx11 = MAP(xmap, ixv + 1, iyv + 1);
		ry11 = MAP(ymap, ixv + 1, iyv + 1);
	      }
	    else
 {
		rx11 = 0.0;
		ry11 = 0.0;
	      }

	    rx = rx00 * (rrx - 1.0) * (rry - 1.0)
	      - rx10 * rrx * (rry - 1.0) 
	      - rx01 * (rrx - 1.0) * rry
	      + rx11 * rrx * rry;
	    ry = ry00 * (rrx - 1.0) * (rry - 1.0)
	      - ry10 * rrx * (rry - 1.0) 
	      - ry01 * (rrx - 1.0) * rry
	      + ry11 * rrx * rry;
	    rx = factor * rx + (factor - 1.0) / 2.0;
	    ry = factor * ry + (factor - 1.0) / 2.0;
	    irx = ((int) (rx + 1.0)) - 1;
	    iry = ((int) (ry + 1.0)) - 1;
	    rrx = rx - irx;
	    rry = ry - iry;
	    if (irx < 0 || irx >= rw || irx == rw-1 && rrx > 0.0 ||
		iry < 0 || iry >= rh || iry == rh-1 && rry > 0.0)
	      continue;

	    ++nPoints;
	    r00 = RIMAGE(cref, irx, iry);
	    if (iry < rh-1)
	      r01 = RIMAGE(cref, irx, iry + 1);
	    else
	      r01 = 0.0;
	    if (irx < rw-1)
	      {
		r10 = RIMAGE(cref, irx + 1, iry);
		if (iry < rh-1)
		  r11 = RIMAGE(cref, irx + 1, iry + 1);
		else
		  r11 = 0.0;
	      }
	    else
	      {
		r10 = 0.0;
		r11 = 0.0;
	      }
	    rv = r00 * (rrx - 1.0) * (rry - 1.0)
	      - r10 * rrx * (rry - 1.0) 
	      - r01 * (rrx - 1.0) * rry
	      + r11 * rrx * rry;
	    if (isnan(rv))
	      Error("ISNAN internal error\nrx = %f irx = %d ry = %f iry = %d rrx = %f rry = %f\n",
		    rx, irx, ry, iry, rrx, rry);
	    iv = IIMAGE(cimage, x, y);
	    si += iv;
	    si2 += iv * iv;
	    sr += rv;
	    sr2 += rv * rv;
	    sir += iv * rv;
	  }
      requiredPoints = ((size_t) iw) * ih * minOverlap;
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
      dPoints = 0;
      de = 0.0;
      for (y = -1; y < mph; ++y)
	for (x = -1; x < mpw; ++x)
	  {

	    x0 = MAP(xmap, x, y);
	    y0 = MAP(ymap, x, y);
	    x1 = MAP(xmap, x, y+1);
	    y1 = MAP(ymap, x, y+1);
	    x2 = MAP(xmap, x+1, y+1);
	    y2 = MAP(ymap, x+1, y+1);
	    x3 = MAP(xmap, x+1, y);
	    y3 = MAP(ymap, x+1, y);

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
	  xv = (2.0 * cpts[i].ix - lFactor + 1) / (2.0 * lFactor);
	  yv = (2.0 * cpts[i].iy - lFactor + 1) / (2.0 * lFactor);
	  ixv = ((int) (xv + 2.0)) - 2;
	  iyv = ((int) (yv + 2.0)) - 2;
	  rrx = xv - ixv;
	  rry = yv - iyv;
	  if (ixv < -1 || ixv > mpw || ixv == mpw && rrx > 0.0 ||
	      iyv < -1 || iyv > mph || iyv == mph && rry > 0.0)
	    {
	      cpts[i].energy = 1000000.0;
	      correspondence += cpts[i].energy;
	      continue;
	    }
	  rx00 = MAP(xmap, ixv, iyv);
	  ry00 = MAP(ymap, ixv, iyv);
	  if (iyv < mph)
	    {
	      rx01 = MAP(xmap, ixv, iyv + 1);
	      ry01 = MAP(ymap, ixv, iyv + 1);
	    }
	  else
	    {
	      rx01 = 0.0;
	      ry01 = 0.0;
	    }
	  if (ixv < mpw)
	    {
	      rx10 = MAP(xmap, ixv + 1, iyv);
	      ry10 = MAP(ymap, ixv + 1, iyv);
	    }
	  else
	    {
	      rx10 = 0.0;
	      ry10 = 0.0;
	    }
	  if (ixv < mpw && iyv < mph)
	    {
	      rx11 = MAP(xmap, ixv + 1, iyv + 1);
	      ry11 = MAP(ymap, ixv + 1, iyv + 1);
	    }
	  else
	    {
	      rx11 = 0.0;
	      ry11 = 0.0;
	    }

	  rx = rx00 * (rrx - 1.0) * (rry - 1.0)
	    - rx10 * rrx * (rry - 1.0) 
	    - rx01 * (rrx - 1.0) * rry
	    + rx11 * rrx * rry;
	  ry = ry00 * (rrx - 1.0) * (rry - 1.0)
	    - ry10 * rrx * (rry - 1.0) 
	    - ry01 * (rrx - 1.0) * rry
	    + ry11 * rrx * rry;
	  rx = lFactor * rx + (lFactor - 1.0) / 2.0;
	  ry = lFactor * ry + (lFactor - 1.0) / 2.0;
	  rrx = cpts[i].rx;
	  rry = cpts[i].ry;
	  distance = kFactor * (sqrt((rrx-rx) * (rrx-rx) + (rry-ry) * (rry-ry)) -
				c.correspondenceThreshold);
	  if (distance < 0.0)
	    cpts[i].energy = 0.0;
	  else
	    cpts[i].energy = distance;
	  correspondence += cpts[i].energy;
#if 0
	  printf("cci-A %d ix %d iy %d rx %d ry %d\n",
		 i, cpts[i].ix, cpts[i].iy, cpts[i].rx, cpts[i].ry);
	  printf("cci-B %d xv %f yv %f rx00 %f ry00 %f rx %f ry %f rrx %f rry %f dist %f en %f corr %f\n",
		 i, xv, yv, rx00, ry00, rx, ry, rrx, rry, distance, cpts[i].energy,
		 correspondence);
#endif
	}

      e = distortion * c.distortion - correlation + correspondence * c.correspondence;
      Log("Total energy of level %d map is %f (dist %f - correl %f + corresp %f)\n",
	  level, e, distortion * c.distortion, correlation,
	  correspondence * c.correspondence);

      mSize = ((size_t) mphp) * mpwp;
      xprop = (float*) malloc(mSize * sizeof(float));
      yprop = (float*) malloc(mSize * sizeof(float));
      changed = (unsigned char*) malloc(mSize * sizeof(unsigned char));
      memset(changed, 0, mSize * sizeof(unsigned char));
      
      fflush(stdout);

      multiplier = 1 << level;
      if (multiplier > 4096)
	multiplier = 4096;
      goalCoverage = (c.quality * mpw) * mph * multiplier;
      Log("for level %d goalCoverage is %f\n", level, goalCoverage);
      actualCoverage = 0.0;

      rawMoveCount = 0;
      moveCount = 0;
      acceptedMoveCount = 0;
      nShiftMoves = 0;

      memset(statSqrtRadius, 0, 21*sizeof(size_t));
      memset(statSqrtOffset, 0, 21*sizeof(size_t));
      memset(statTheta, 0, 21*sizeof(size_t));
      memset(statX, 0, 21*sizeof(size_t));
      memset(statY, 0, 21*sizeof(size_t));

      logMinRadius = log(0.01 / factor); /* no use going smaller than 1% of the pixel size */
      logMaxRadius = log(0.5); /* never go more than half the distance to the next map point */
      while (actualCoverage < goalCoverage)
	{
	  ++rawMoveCount;

	  /* choose point in mapping to distort */
	  /*    center(x,y) radius(r) offset(x,y) */
#if 0
	  icx = ((int) floor(drand48() * (mpw+2))) - 1;
	  icy = ((int) floor(drand48() * (mph+2))) - 1;
#else
	  icx = (rawMoveCount % (mpw+2)) - 1; 
	  icy = ((rawMoveCount / (mpw+2)) % (mph+2)) - 1;
#endif
	  cx = MAP(xmap, icx, icy);
	  cy = MAP(ymap, icx, icy);

	  /* choose size of move to make */
	  minDist2 = 1.0e+30;
	  if (icx >= 0)
	    {
	      if (icy >= 0)
		{
		  dx = MAP(xmap, icx-1, icy-1) - cx;
		  dy = MAP(ymap, icx-1, icy-1) - cy;
		  d2 = dx*dx + dy*dy;
		  if (d2 < minDist2)
		    minDist2 = d2;
		}
	      if (icy < mph)
		{
		  dx = MAP(xmap, icx-1, icy+1) - cx;
		  dy = MAP(ymap, icx-1, icy+1) - cy;
		  d2 = dx*dx + dy*dy;
		  if (d2 < minDist2)
		    minDist2 = d2;
		}
	      dx = MAP(xmap, icx-1, icy) - cx;
	      dy = MAP(ymap, icx-1, icy) - cy;
	      d2 = dx*dx + dy*dy;
	      if (d2 < minDist2)
		minDist2 = d2;
	    }
	  if (icx < mpw)
	    {
	      if (icy >= 0)
		{
		  dx = MAP(xmap, icx+1, icy-1) - cx;
		  dy = MAP(ymap, icx+1, icy-1) - cy;
		  d2 = dx*dx + dy*dy;
		  if (d2 < minDist2)
		    minDist2 = d2;
		}
	      if (icy < mph)
		{
		  dx = MAP(xmap, icx+1, icy+1) - cx;
		  dy = MAP(ymap, icx+1, icy+1) - cy;
		  d2 = dx*dx + dy*dy;
		  if (d2 < minDist2)
		    minDist2 = d2;
		}
	      dx = MAP(xmap, icx+1, icy) - cx;
	      dy = MAP(ymap, icx+1, icy) - cy;
	      d2 = dx*dx + dy*dy;
	      if (d2 < minDist2)
		minDist2 = d2;
	    }
	  if (icy >= 0)
	    {
	      dx = MAP(xmap, icx, icy-1) - cx;
	      dy = MAP(ymap, icx, icy-1) - cy;
	      d2 = dx*dx + dy*dy;
	      if (d2 < minDist2)
		minDist2 = d2;
	    }
	  if (icy < mph)
	    {
	      dx = MAP(xmap, icx, icy+1) - cx;
	      dy = MAP(ymap, icx, icy+1) - cy;
	      d2 = dx*dx + dy*dy;
	      if (d2 < minDist2)
		minDist2 = d2;
	    }
	  logMaxRadius = 0.5 * log(minDist2) - M_LN2;
	  if (logMaxRadius < logMinRadius)
	    continue;
	  radius = exp(drand48() * (logMaxRadius - logMinRadius) + logMinRadius);
	  theta = drand48() * 2.0 * M_PI;
	  
	  actualCoverage += 1.0;
	  ++moveCount;

	  MAP(xprop, icx, icy) = cx + radius * cos(theta);
	  MAP(yprop, icx, icy) = cy + radius * sin(theta);
	  MAP(changed, icx, icy) = 1;
	  
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
	  dPoints = 0;
	  sx = (changeMinX - 1) * factor + (factor - 1) / 2;
	  if (sx < 0)
	    sx = 0;
	  ex = (changeMaxX + 1) * factor + factor / 2;
	  if (ex >= iw)
	    ex = iw - 1;
	  sy = (changeMinY - 1) * factor + (factor - 1) / 2;
	  if (sy < 0)
	    sy = 0;
	  ey = (changeMaxY + 1) * factor + factor / 2;
	  if (ey >= ih)
	    ey = ih - 1;

	  /*	  printf("sx = %d %d %d %d  factor = %d\n", sx, ex, sy, ey, factor); */
	  for (y = sy; y <= ey; ++y)
            for (x = sx; x <= ex; ++x)
	      {
		xv = (2.0 * x - factor + 1) / (2.0 * factor);
		yv = (2.0 * y - factor + 1) / (2.0 * factor);
		ixv = ((int) (xv + 2.0)) - 2;
		iyv = ((int) (yv + 2.0)) - 2;
		rrx = xv - ixv;
		rry = yv - iyv;
		if (ixv < -1 || ixv > mpw || ixv == mpw && rrx > 0.0 ||
		    iyv < -1 || iyv > mph || iyv == mph && rry > 0.0)
		  continue;

		if (!MAP(changed, ixv, iyv) &&
		    !MAP(changed, ixv+1, iyv) &&
		    !MAP(changed, ixv, iyv+1) &&
		    !MAP(changed, ixv+1, iyv+1))
		    continue;

		/* use bilinear interpolation to find value */
		iv = IIMAGE(cimage, x, y);
		if (iv < 0.0)
		  Error("Internal error: iv out of range: %f\n", iv);
		rx00 = MAP(xmap, ixv, iyv);
		ry00 = MAP(ymap, ixv, iyv);
		if (iyv < mph)
		  {
		    rx01 = MAP(xmap, ixv, iyv + 1);
		    ry01 = MAP(ymap, ixv, iyv + 1);
		  }
		else
		  {
		    rx01 = 0.0;
		    ry01 = 0.0;
		  }
		if (ixv < mpw)
		  {
		    rx10 = MAP(xmap, ixv + 1, iyv);
		    ry10 = MAP(ymap, ixv + 1, iyv);
		  }
		else
		  {
		    rx10 = 0.0;
		    ry10 = 0.0;
		  }
		if (ixv < mpw && iyv < mph)
		  {
		    rx11 = MAP(xmap, ixv + 1, iyv + 1);
		    ry11 = MAP(ymap, ixv + 1, iyv + 1);
		  }
		else
		  {
		    rx11 = 0.0;
		    ry11 = 0.0;
		  }
		rx = rx00 * (rrx - 1.0) * (rry - 1.0)
		  - rx10 * rrx * (rry - 1.0) 
		  - rx01 * (rrx - 1.0) * rry
		  + rx11 * rrx * rry;
		ry = ry00 * (rrx - 1.0) * (rry - 1.0)
		  - ry10 * rrx * (rry - 1.0) 
		  - ry01 * (rrx - 1.0) * rry
		  + ry11 * rrx * rry;
		rx = factor * rx + (factor - 1.0) / 2.0;
		ry = factor * ry + (factor - 1.0) / 2.0;
		irx = ((int) (rx + 1.0)) - 1;
		iry = ((int) (ry + 1.0)) - 1;
		rrx = rx - irx;
		rry = ry - iry;
		if (irx >= 0 && (irx < rw-1 || irx == rw-1 && rrx == 0.0) &&
		    iry >= 0 && (iry < rh-1 || iry == rh-1 && rry == 0.0))
		  {
			r00 = RIMAGE(cref, irx, iry);
			if (iry < rh-1)
			  r01 = RIMAGE(cref, irx, iry + 1);
			else
			  r01 = 0.0;
			if (irx < rw-1)
			  {
			    r10 = RIMAGE(cref, irx + 1, iry);
			    if (iry < rh-1)
			      r11 = RIMAGE(cref, irx + 1, iry + 1);
			    else
			      r11 = 0.0;
			  }
			else
			  {
			    r10 = 0.0;
			    r11 = 0.0;
			  }
			rv = r00 * (rrx - 1.0) * (rry - 1.0)
			  - r10 * rrx * (rry - 1.0) 
			  - r01 * (rrx - 1.0) * rry
			  + r11 * rrx * rry;
			dsi -= iv;
			dsi2 -= iv * iv;
			dsir -= iv * rv;
			dsr2 -= rv * rv;
			dsr -= rv;
			dPoints -= 1;
		  }
		
		/* use bilinear interpolation to find value */
		if (MAP(changed, ixv, iyv))
		  {
		    rx00 = MAP(xprop, ixv, iyv);
		    ry00 = MAP(yprop, ixv, iyv);
		  }
		else
		  {
		    rx00 = MAP(xmap, ixv, iyv);
		    ry00 = MAP(ymap, ixv, iyv);
		  }
		if (iyv < mph)
		  if (MAP(changed, ixv, iyv+1))
		    {
		      rx01 = MAP(xprop, ixv, iyv + 1);
		      ry01 = MAP(yprop, ixv, iyv + 1);
		    }
		  else
		    {
		      rx01 = MAP(xmap, ixv, iyv + 1);
		      ry01 = MAP(ymap, ixv, iyv + 1);
		    }
		else
		  {
		    rx01 = 0.0;
		    ry01 = 0.0;
		  }
		if (ixv < mpw)
		  if (MAP(changed, ixv+1, iyv))
		    {
		      rx10 = MAP(xprop, ixv + 1, iyv);
		      ry10 = MAP(yprop, ixv + 1, iyv);
		    }
		  else
		    {
		      rx10 = MAP(xmap, ixv + 1, iyv);
		      ry10 = MAP(ymap, ixv + 1, iyv);
		    }
		else
		  {
		    rx10 = 0.0;
		    ry10 = 0.0;
		  }
		if (ixv < mpw && iyv < mph)
		  if (MAP(changed, ixv+1, iyv+1))
		    {
		      rx11 = MAP(xprop, ixv + 1, iyv + 1);
		      ry11 = MAP(yprop, ixv + 1, iyv + 1);
		    }
		  else
		    {
		      rx11 = MAP(xmap, ixv + 1, iyv + 1);
		      ry11 = MAP(ymap, ixv + 1, iyv + 1);
		    }
		else
		  {
		    rx11 = 0.0;
		    ry11 = 0.0;
		  }
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
		  Error("Internal error: RY ISNAN %f %f %f %f %f %f\n",
			ry00, ry01, ry10, ry11, rrx, rry);
		rx = factor * rx + (factor - 1.0) / 2.0;
		ry = factor * ry + (factor - 1.0) / 2.0;

		irx = ((int) (rx + 1.0)) - 1;
		iry = ((int) (ry + 1.0)) - 1;
		rrx = rx - irx;
		rry = ry - iry;
		if (irx < 0 || irx >= rw || irx == rw-1 && rrx > 0.0 ||
		    iry < 0 || iry >= rh || iry == rh-1 && rry > 0.0)
		  continue;
		r00 = RIMAGE(cref, irx, iry);
		if (iry < rh-1)
		  r01 = RIMAGE(cref, irx, iry + 1);
		else
		  r01 = 0.0;
		if (irx < rw-1)
		  {
		    r10 = RIMAGE(cref, irx + 1, iry);
		    if (iry < rh-1)
		      r11 = RIMAGE(cref, irx + 1, iry + 1);
		    else
		      r11 = 0.0;
		  }
		else
		  {
		    r10 = 0.0;
		    r11 = 0.0;
		  }
		rv = r00 * (rrx - 1.0) * (rry - 1.0)
		  - r10 * rrx * (rry - 1.0) 
		  - r01 * (rrx - 1.0) * rry
		  + r11 * rrx * rry;
		dsi += iv;
		dsi2 += iv * iv;
		dsir += iv * rv;
		dsr2 += rv * rv;
		dsr += rv;
		dPoints += 1;
	      }

	  if (nPoints < requiredPoints)
	    {
	      dsr -= (requiredPoints - nPoints) * 255.0;
	      dsr2 -= (requiredPoints - nPoints) * 255.0 * 255.0;
	    }
	  newPoints = nPoints + dPoints;
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
              if (y < -1 || y >= mph)
		continue;
              for (x = changeMinX - 1; x <= changeMaxX; ++x)
      	        {
	          if (x < -1 || x >= mpw)
	            continue;

		  upd0 = MAP(changed, x, y);
		  upd1 = MAP(changed, x, y+1);
		  upd2 = MAP(changed, x+1, y+1);
		  upd3 = MAP(changed, x+1, y);
		      
		  if (!upd0 && !upd1 && !upd2 && !upd3)
		    continue;

		  x0 = MAP(xmap, x, y);
		  y0 = MAP(ymap, x, y);
		  x1 = MAP(xmap, x, y+1);
		  y1 = MAP(ymap, x, y+1);
		  x2 = MAP(xmap, x+1, y+1);
		  y2 = MAP(ymap, x+1, y+1);
		  x3 = MAP(xmap, x+1, y);
		  y3 = MAP(ymap, x+1, y);
		  if (upd0)
		    {
		      nx0 = MAP(xprop, x, y);
		      ny0 = MAP(yprop, x, y);
		    }
		  else
		    {
		      nx0 = x0;
		      ny0 = y0;
		    }
		  if (upd1)
		    {
		      nx1 = MAP(xprop, x, y+1);
		      ny1 = MAP(yprop, x, y+1);
		    }
		  else
		    {
		      nx1 = x1;
		      ny1 = y1;
		    }
		  if (upd2)
		    {
		      nx2 = MAP(xprop, x+1, y+1);
		      ny2 = MAP(yprop, x+1, y+1);
		    }
		  else
		    {
		      nx2 = x2;
		      ny2 = y2;
		    }
		  if (upd3)
		    {
		      nx3 = MAP(xprop, x+1, y);
		      ny3 = MAP(yprop, x+1, y);
		    }
		  else
		    {
		      nx3 = x3;
		      ny3 = y3;
		    }

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
	  newDistortion = distortion + de / (mph+1) / (mpw+1);

	  /* add in contribution from correspondence energy */
	  ce = 0.0;
	  for (i = 0; i < nCpts; ++i)
	    {
	      xv = (2.0 * cpts[i].ix - lFactor + 1) / (2.0 * lFactor);
	      yv = (2.0 * cpts[i].iy - lFactor + 1) / (2.0 * lFactor);
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
	      if (ixv < -1 || ixv > mpw || ixv == mpw && rrx > 0.0 ||
		  iyv < -1 || iyv > mph || iyv == mph && rry > 0.0)
		{
		  cpts[i].newEnergy = 1000000.0;
		  continue;
		}
	      if (MAP(changed, ixv, iyv))
		{
		  rx00 = MAP(xprop, ixv, iyv);
		  ry00 = MAP(yprop, ixv, iyv);
		}
	      else
		{
		  rx00 = MAP(xmap, ixv, iyv);
		  ry00 = MAP(ymap, ixv, iyv);
		}
	      if (iyv < mph)
		if (MAP(changed, ixv, iyv+1))
		  {
		    rx01 = MAP(xprop, ixv, iyv + 1);
		    ry01 = MAP(yprop, ixv, iyv + 1);
		  }
		else
		  {
		    rx01 = MAP(xmap, ixv, iyv + 1);
		    ry01 = MAP(ymap, ixv, iyv + 1);
		  }
	      else
		{
		  rx01 = 0.0;
		  ry01 = 0.0;
		}
	      if (ixv < mpw)
		if (MAP(changed, ixv+1, iyv))
		  {
		    rx10 = MAP(xprop, ixv + 1, iyv);
		    ry10 = MAP(yprop, ixv + 1, iyv);
		  }
		else
		  {
		    rx10 = MAP(xmap, ixv + 1, iyv);
		    ry10 = MAP(ymap, ixv + 1, iyv);
		  }
	      else
		{
		  rx10 = 0.0;
		  ry10 = 0.0;
		}
	      if (ixv < mpw && iyv < mph)
		if (MAP(changed, ixv+1, iyv+1))
		  {
		    rx11 = MAP(xprop, ixv + 1, iyv + 1);
		    ry11 = MAP(yprop, ixv + 1, iyv + 1);
		  }
		else
		  {
		    rx11 = MAP(xmap, ixv + 1, iyv + 1);
		    ry11 = MAP(ymap, ixv + 1, iyv + 1);
		  }
	      else
		{
		  rx11 = 0.0;
		  ry11 = 0.0;
		}

	      rx = rx00 * (rrx - 1.0) * (rry - 1.0)
		- rx10 * rrx * (rry - 1.0) 
		- rx01 * (rrx - 1.0) * rry
		+ rx11 * rrx * rry;
	      ry = ry00 * (rrx - 1.0) * (rry - 1.0)
		- ry10 * rrx * (rry - 1.0) 
		- ry01 * (rrx - 1.0) * rry
		+ ry11 * rrx * rry;
	      rx = lFactor * rx + (lFactor - 1.0) / 2.0;
	      ry = lFactor * ry + (lFactor - 1.0) / 2.0;
	      rrx = cpts[i].rx;
	      rry = cpts[i].ry;
	      distance = kFactor * (sqrt((rrx-rx) * (rrx-rx) + (rry-ry) * (rry-ry)) -
				    c.correspondenceThreshold);
	      if (distance < 0.0)
		cpts[i].newEnergy = 0.0;
	      else
		cpts[i].newEnergy = distance;
	      ce += cpts[i].newEnergy - cpts[i].energy;
	    }
	  newCorrespondence = correspondence + ce;

	  newE = newDistortion * c.distortion - newCorrelation +
	    newCorrespondence * c.correspondence;

	  /* accept or reject move */
	  accept = (newE < e);
	  /*	  printf("MAPCH = %d  %f %f\n", MAP(changed, 1, 1),
		  MAP(xmap, 1, 1), MAP(xprop, 1, 1)); */
	  for (y = changeMinY; y <= changeMaxY; ++y)
	    for (x = changeMinX; x <= changeMaxX; ++x)
	      if (MAP(changed, x, y))
		{
		  MAP(changed, x, y) = 0;
		  if (accept)
		    {
		      /* make the provisional state the new state */
		      /*		      printf("Moving point (%d,%d) from (%f, %f) to (%f, %f)\n",
					      x, y,
					      MAP(xmap, x, y), MAP(ymap, x, y),
					      MAP(xprop, x, y), MAP(yprop, x, y)); */
		      MAP(xmap, x, y) = MAP(xprop, x, y);
		      MAP(ymap, x, y) = MAP(yprop, x, y);
		    }
		}
	  if (accept)
	    {
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
	      e = newE;
	      ++acceptedMoveCount;
	    }
	}

      Log("Level %d: after %d raw moves, %d moves, and %d accepted moves...\n",
	  level, rawMoveCount, moveCount, acceptedMoveCount);
      Log("          energy is %f\n", e);
#if 0
      Log("Accepted shift move statistics:\n");
      Log("   range:  sqrt(radius) sqrt(offset)    theta       x        y\n");
      for (i = 0; i < 20; ++i)
	    Log(" %d-%d%%:  %10.5f  %10.5f  %10.5f  %10.5f  %10.5f\n",
		5*i, 5*i+5,
		(100.0 * statSqrtRadius[i]) / nShiftMoves,
		(100.0 * statSqrtOffset[i]) / nShiftMoves,
		(100.0 * statTheta[i]) / nShiftMoves,
		(100.0 * statX[i]) / nShiftMoves,
		(100.0 * statY[i]) / nShiftMoves);
#endif
      Log("--- done with level %d ---\n", level);

#if 0
      if (level > 7)
	{
	  printf("Level %d map:\n", level);
	  for (y = 0; y < mph; ++y)
	    for (x = 0; x < mpw; ++x)
	      {
		printf("%d %d: %f %f\n",
		       y, x, MAP(xmap, x, y), MAP(ymap, x, y));
	      }
	}
#endif

      free(xprop);
      free(yprop);
      free(changed);

      /* compute the warped image */
      if (c.outputImagesName[0] != '\0')
        {
          warp = warped[level];
          for (y = 0; y < mph; ++y)
            for (x = 0; x < mpw; ++x)
              {
	        /* use bilinear interpolation to find value */
	        rx = MAP(xmap, x, y);
	        ry = MAP(ymap, x, y);
	        irx = (int) floor(rx);
	        iry = (int) floor(ry);
	        rrx = rx - irx;
	        rry = ry - iry;

		if (irx >= 0 && irx < refw)
		  {
		    if (iry >= 0 && iry < refh)
		      r00 = REFIMAGE(ref, irx, iry);
		    else
		      r00 = 0.0;
		    if (iry >= -1 && iry < refh-1)
		      r01 = REFIMAGE(ref, irx, iry + 1);
		    else
		      r01 = 0.0;
		  }
		else
		  {
		    r00 = 0.0;
		    r01 = 0.0;
		  }
		if (irx >= -1 && irx < refw-1)
		  {
		    if (iry >= 0 && iry < refh)
		      r10 = REFIMAGE(ref, irx + 1, iry);
		    else
		      r10 = 0.0;
		    if (iry >= -1 && iry < refh-1)
		      r11 = REFIMAGE(ref, irx + 1, iry + 1);
		    else
		      r11 = 0.0;
		  }
		else
		  {
		    r10 = 0.0;
		    r11 = 0.0;
		  }
	        rv = r00 * (rrx - 1.0) * (rry - 1.0)
	          - r10 * rrx * (rry - 1.0) 
	          - r01 * (rrx - 1.0) * rry
	          + r11 * rrx * rry;
	        warp[y * mpw + x] = rv;
	      }
          /* printf("Computed warped image for level %d\n", level); */
	}

      if (level == c.outputLevel || c.writeAllMaps)
	WriteMap(outputName, level);
    }

  // write out the score for this mapping
  sprintf(fn, "%s.score", outputName);
  f = fopen(fn, "w");
  fprintf(f, "%f %f %f %f\n",
	  correlation
	  - c.distortion * distortion
	  - c.correspondence * correspondence,
	  correlation, distortion, correspondence);
  fclose(f);

  memcpy(&(r.pair), &(t.pair), sizeof(Pair));
  r.updated = 1;
  r.correlation = correlation;
  r.distortion = distortion;
  r.correspondence = correspondence;
  r.message[0] = '\0';
}

void
WriteMap (char *outputName, int level)
{
  char fn[PATH_MAX];
  int mpw;
  int mph;
  int refSlice;
  float *xmap;
  float *ymap;
  float *cmap;
  size_t cmapCount1;
  size_t cmapCount2;
  size_t cmapCount3;
  int mpwp, mphp;
  int x, y;
  FILE *f;
  float rx, ry;
  int irx, iry;
  size_t mSize;
  unsigned int rw, rh;
  
  /* write out resultant map */
  if (level > c.outputLevel)
    sprintf(fn, "%s.%0.2d.map", outputName, level);
  else
    sprintf(fn, "%s.map", outputName);

  f = fopen(fn, "w");
  if (f == NULL)
    Error("Could not open map file %s for writing\n", fn);
  mpw = imageWidth[level];
  mph = imageHeight[level];
  mpwp = mpw + 2;
  mphp = mph + 2;
  mSize = ((size_t) mphp) * mpwp;
  rw = referenceWidth[level];
  rh = referenceHeight[level];
  xmap = xmaps[level];
  ymap = ymaps[level];
  cmap = (float *) malloc(mSize*sizeof(float));
  cmapCount1 = 0;
  cmapCount2 = 0;
  cmapCount3 = 0;
  for (y = -1; y <= mph; ++y)
    {    
      for (x = -1; x <= mpw; ++x)
	if (y < 0 || y >= mph ||
	    x < 0 || x >= mpw)
	  MAP(cmap, x, y) = 0.0;
	else
	  {
	    rx = MAP(xmap, x, y);
	    ry = MAP(ymap, x, y);
	    irx = ((int) (rx + 1.5)) - 1;
	    iry = ((int) (ry + 1.5)) - 1;
	    if (irx >= 0 && irx < rw && iry >= 0 && iry < rh)
	      {
		++cmapCount1;
	        ++cmapCount2;
		++cmapCount3;
		MAP(cmap, x, y) = 1.0;
	      }
	    else
	      MAP(cmap, x, y) = 0.0;
	  }
    }
  /* printf("CMAPCOUNT for %s is %ld %ld %ld\n", fn, cmapCount1,
     cmapCount2, cmapCount3); */

  fprintf(f, "M1\n");
  fprintf(f, "%d %d\n", mpw, mph);
  fprintf(f, "%d %d\n", t.pair.imageMinX, t.pair.imageMinY);
  fprintf(f, "%s %s\n", t.pair.imageName, t.pair.refName);
  fwrite(xmap, mSize*sizeof(float), 1, f);
  fwrite(ymap, mSize*sizeof(float), 1, f);
  fwrite(cmap, mSize*sizeof(float), 1, f);
  fclose(f);
  free(cmap);
}

void
PackContext ()
{
  int i;

  par_pkbyte((unsigned char) c.type);
  par_pkstr(c.imageBasename);

  par_pkstr(c.referenceBasename);

  par_pkstr(c.cptsName);

  par_pkstr(c.inputMapName);

  par_pkstr(c.outputMapBasename);
  par_pkstr(c.outputImagesName);
  par_pkstr(c.outputCorrelationName);

  par_pkint(c.outputLevel);
  par_pkint(c.minResolution);
  par_pkint(c.depth);
  par_pkdouble(c.distortion);
  par_pkdouble(c.correspondence);
  par_pkdouble(c.correspondenceThreshold);
  par_pkdouble(c.quality);

  par_pkint(c.writeAllMaps);
  par_pkdouble(c.afterTime);
}

void
UnpackContext ()
{
  c.type = (char) par_upkbyte();
  par_upkstr(c.imageBasename);

  par_upkstr(c.referenceBasename);

  par_upkstr(c.cptsName);

  par_upkstr(c.inputMapName);

  par_upkstr(c.outputMapBasename);
  par_upkstr(c.outputImagesName);
  par_upkstr(c.outputCorrelationName);

  c.outputLevel = par_upkint();
  c.minResolution = par_upkint();
  c.depth = par_upkint();
  c.distortion = par_upkdouble();
  c.correspondence = par_upkdouble();
  c.correspondenceThreshold = par_upkdouble();
  c.quality = par_upkdouble();

  c.writeAllMaps = par_upkint();
  c.afterTime = par_upkdouble();
}

void
PackPair (Pair *p)
{
  par_pkstr(p->imageName);
  par_pkint(p->imageMinX);
  par_pkint(p->imageMaxX);
  par_pkint(p->imageMinY);
  par_pkint(p->imageMaxY);
  par_pkstr(p->refName);
  par_pkint(p->refMinX);
  par_pkint(p->refMaxX);
  par_pkint(p->refMinY);
  par_pkint(p->refMaxY);
  par_pkstr(p->pairName);
}

void
UnpackPair (Pair *p)
{
  char s[PATH_MAX];

  if (p->imageName != NULL)
    free(p->imageName);
  par_upkstr(s);
  p->imageName = (char*) malloc(strlen(s)+1);
  strcpy(p->imageName, s);

  p->imageMinX = par_upkint();
  p->imageMaxX = par_upkint();
  p->imageMinY = par_upkint();
  p->imageMaxY = par_upkint();

  if (p->refName != NULL)
    free(p->refName);
  par_upkstr(s);
  p->refName = (char*) malloc(strlen(s)+1);
  strcpy(p->refName, s);

  p->refMinX = par_upkint();
  p->refMaxX = par_upkint();
  p->refMinY = par_upkint();
  p->refMaxY = par_upkint();

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

  if (logFile == NULL)
    {
      int i;
      char logName[128];
      i = par_instance();
      if (i < 0)
	sprintf(logName, "logs/00.log");
      else if (i < 100)
	sprintf(logName, "logs/%0.2d.log", i);
      else
	sprintf(logName, "logs/%d.log", i);
      logFile = fopen(logName, "w");
      if (logFile == NULL)
	Error("Could not open log file %s\n", logName);
    }

  va_start(args, fmt);
  fprintf(logFile, "%f: ", MPI_Wtime());
  vfprintf(logFile, fmt, args);
  va_end(args);
  fflush(logFile);
}

