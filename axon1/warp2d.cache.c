#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <unistd.h>
#include <limits.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <mpi.h>

#include "par.h"

#define CHECKING	0
#define TRACKING	0
#define DEBUG_MOVES	0
#define MASKING		1

#define MAX_LEVELS	32       /* max image size (w or h) is 2^MAX_LEVELS */
#define MAP(map,x,y)	map[((y)+1)*((size_t) mpwp) + (x) + 1]
#define MAP1(map,x,y)	map[((y)+1)*((size_t) mpw1p) + (x) + 1]
#define MAPINP(map,x,y)	map[((y)+1)*((size_t) mpwi) + (x) + 1]
#define MAPCON(map,x,y)	map[((y)+1)*((size_t) mpwcp) + (x) + 1]
#define IIMAGE(i,x,y)   i[(y)*((size_t) iw) + (x)]
#define RIMAGE(i,x,y)   i[(y)*((size_t) rw) + (x)]
#define REFIMAGE(i,x,y) i[(y)*((size_t) refw) + (x)]
#define LINE_LENGTH	255

typedef struct Context {
  /* NOTE: any new fields added to this struct should
     be also added to PackContext and UnpackContext */ 
  char type;
  char inputName[PATH_MAX];
  int inputDigits;

  char maskName[PATH_MAX];
  int maskDigits;

  char referenceName[PATH_MAX];
  char referenceMaskName[PATH_MAX];

  char cptsName[PATH_MAX];

  char inputMapName[PATH_MAX];
  char constrainingMapName[PATH_MAX];

  int pairwiseNaming;
  char outputName[PATH_MAX];
  char outputImagesName[PATH_MAX];
  char outputCorrelationName[PATH_MAX];

  int outputLevel;
  int minResolution;
  int depth;
  double distortion;
  double correspondence;
  double correspondenceThreshold;
  double quality;

  int writeAllMaps;
  double afterTime;
} Context;

typedef struct Task {
  /* NOTE: any new fields added to this struct should
     be also added to PackTask and UnpackTask */ 
  int refSlice;
  int imageSlice;
} Task;

typedef struct Result {
  int refSlice;
  int imageSlice;
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
unsigned char *ref_mask_in;
unsigned char *image_in;
unsigned char *image_mask_in;

float* mapping;
int nLevels;
float* refs[MAX_LEVELS];       /* reference image to warp to (at various
			         resolutions) */
float* images[MAX_LEVELS];     /* image to be warped (at various resolutions) */
float* warped[MAX_LEVELS];     /* warped reference image (at various resol.) */

int imageWidth[MAX_LEVELS], imageHeight[MAX_LEVELS];
int referenceWidth[MAX_LEVELS], referenceHeight[MAX_LEVELS];

int imageMaskWidth, imageMaskHeight;
int imageMaskFactor;
int referenceMaskWidth, referenceMaskHeight;
int referenceMaskFactor;
float* xmaps[MAX_LEVELS];
float* ymaps[MAX_LEVELS];
unsigned char *imask[MAX_LEVELS];
unsigned char *rmask[MAX_LEVELS];
unsigned char *eimask[MAX_LEVELS];
unsigned char *ermask[MAX_LEVELS];
size_t imaskCount[MAX_LEVELS];
size_t rmaskCount[MAX_LEVELS];
size_t minMaskCount[MAX_LEVELS];
size_t eimaskCount[MAX_LEVELS];
size_t ermaskCount[MAX_LEVELS];
int inputMapFactor;
float *inputMapX, *inputMapY, *inputMapC;
int mpwi, mphi;
int mpwip, mphip;
int constrainingMapFactor;
float *constrainingMapX, *constrainingMapY, *constrainingMapC;
int mpwc, mphc;
int mpwcp, mphcp;


#if 0
float* isig[MAX_LEVELS];
float* rsig[MAX_LEVELS];
#endif

int nCpts = 0;
CPoint *cpts = 0;

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
void Init();
void WriteMap (int level);
void Compute();
int Compare (const void *x, const void *y);
int SortBySlice (const void *x, const void *y);
int SortByEnergy (const void *x, const void *y);
int ParseRange (char *s, int *pos, int *minValue, int *maxValue);
int ParseValue (char *s, int *pos, int *value);
int ReadHeader (FILE *f, char *tc, unsigned int *w, unsigned int *h, int *m);
size_t CountBits (unsigned char *p, size_t n);
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
  int nSkipSlices;
  int maxSkipSlices;
  int *skipSlices;
  int len;
  char skipFile[PATH_MAX];
  char pairsFile[PATH_MAX];
  int nPairs;
  int *pairs;
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
  int skipIndex;
  int xRes, yRes;
  int minS, maxS;
  int forward;
  char afterFileName[PATH_MAX];
  struct stat statBuf;
  int pn;
  int imgn, refn;
  unsigned int iw, ih;

  error = 0;
  c.type = '\0';
  c.inputName[0] = '\0';
  c.maskName[0] = '\0';
  c.referenceName[0] = '\0';
  c.referenceMaskName[0] = '\0';
  c.inputMapName[0] = '\0';
  c.cptsName[0] = '\0';
  c.constrainingMapName[0] = '\0';
  c.pairwiseNaming = 0;
  c.outputName[0] = '\0';
  c.outputImagesName[0] = '\0';
  c.outputCorrelationName[0] = '\0';
  c.outputLevel = 0;
  c.minResolution = 1;
  c.distortion = 1.0;
  c.correspondence = 1.0;
  c.correspondenceThreshold = 0.0;
  c.quality = 10.0;
  c.writeAllMaps = 0;
  c.depth = 0;
  r.message = (char *) malloc(PATH_MAX + 1024);
  forward = 1;
  skipFile[0] = '\0';
  pairsFile[0] = '\0';
  nPairs = 0;
  pairs = 0;
  minSlice = 0;
  maxSlice = 1000000000;
  afterFileName[0] = '\0';
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
	strcpy(c.inputName, argv[i]);
      }
    else if (strcmp(argv[i], "-mask") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.maskName, argv[i]);
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
	strcpy(c.outputName, argv[i]);
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
    else if (strcmp(argv[i], "-skip") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(skipFile, argv[i]);
      }
    else if (strcmp(argv[i], "-slice") == 0)
      {
	pos = 0;
	if (++i == argc ||
	    !ParseRange(argv[i], &pos, &minS, &maxS) ||
	    argv[i][pos] != '\0')
	  {
	    error = 1;
	    break;
	  }
	if (minS <= maxS)
	  {
	    minSlice = minS;
	    maxSlice = maxS;
	  }
	else
	  {
	    minSlice = maxS;
	    maxSlice = minS;
	  }
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
	strcpy(c.referenceName, argv[i]);
      }
    else if (strcmp(argv[i], "-reference_mask") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.referenceMaskName, argv[i]);
      }
    else error = 1;

  if (error)
    {
      fprintf(stderr, "Usage: warp -input file_prefix -output file_prefix\n");
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
      fprintf(stderr, "              [-reference_mask <reference_mask_prefix>]\n");
      fprintf(stderr, "              [-input_map <input_map_prefix>]\n");
      fprintf(stderr, "   where ranges are expressed as: integer\n");
      fprintf(stderr, "                              or: integer-integer\n");
      fprintf(stderr, "                              or: integer-\n");
      fprintf(stderr, "                              or: -integer\n");
      exit(1);
    }

  /* check that at least minimal parameters were supplied */
  if (c.inputName[0] == '\0' || c.outputName[0] == '\0')
    Error("Both -input and -output parameters must be specified.\n");

  /* get a sorted list of the slices to skip */
  nSkipSlices = 0;
  skipSlices = (int *) malloc(sizeof(int));
  maxSkipSlices = 1;
  if (skipFile[0] != '\0')
    {
      f = fopen(skipFile, "r");
      if (f == NULL)
	Error("Could not open skip file %s\n", skipFile);
      while (fscanf(f, "%d", &n) == 1)
	{
	  if (++nSkipSlices > maxSkipSlices)
	    {
	      maxSkipSlices *= 2;
	      skipSlices = realloc(skipSlices, maxSkipSlices * sizeof(int));
	    }
	  skipSlices[nSkipSlices-1] = n;
	}
      if (!feof(f))
	Error("Read error while reading skip file %s\n", skipFile);
      
      fclose(f);
      qsort(skipSlices, nSkipSlices, sizeof(int), Compare);
    }    

  /* check what images are present */
  if (strlen(c.inputName) > 4 &&
      strcmp(&c.inputName[strlen(c.inputName) - 4], ".pgm") == 0)
    {
      f = fopen(c.inputName, "r");
      if (f == NULL)
	Error("Could not open image file %s", c.inputName);
	
      if (!ReadHeader(f, &tc, &iw, &ih, &m))
	{
	  fclose(f);
	  Error("Invalid image file header: %s\n", c.inputName);
	}
      if (c.type == '\0')
	c.type = tc;
      else if (tc != c.type)
	{
	  fclose(f);
	  Error("Inconsistent mix of pgm and ppm image files: %s\n",
		c.inputName);
	}
      if (m != 255)
	{
	  fclose(f);
	  Error("Image file %s does not have byte-size image components.\n",
		c.inputName);
	}
      fclose(f);
    }
  else
    {
      for (i = strlen(c.inputName)-1; i >= 0 && c.inputName[i] != '/'; --i) ;
      if (i >= 0)
	{
	  strncpy(inputDirName, c.inputName, i+1);
	  inputDirName[i+1] = '\0';
	  strcpy(inputPrefix, &c.inputName[i+1]);
	}
      else
	{
	  strcpy(inputDirName, "./");
	  strcpy(inputPrefix, c.inputName);
	}
      
      /* read the directory to look for files of the form nameNNNN.ppm  (or .pgm) */
      Log("Scanning input directory %s\n", inputDirName);
      dir = opendir(inputDirName);
      if (dir == NULL)
	Error("Could not open directory %s for input files\n");
      inputPrefixLen = strlen(inputPrefix);
      minN = 1000000000;
      maxN = -1;
      c.inputDigits = -1;
      while ((de = readdir(dir)) != NULL)
	{
	  if (strncmp(de->d_name, inputPrefix, inputPrefixLen) == 0 &&
	      (len = strlen(de->d_name)) > 4 &&
	      strcmp(&(de->d_name[len-4]), ".pgm") == 0)
	    {
	      if (sscanf(&(de->d_name[inputPrefixLen]), "%d%n",
			 &n, &digitLen) < 1)
		continue;
	      if (c.inputDigits < 0)
		c.inputDigits = digitLen;
	      else if (digitLen != c.inputDigits)
		Error("Inconsistent number of digits in input file names: %s\n", de->d_name);
	      if (n < minN)
		{
		  minN = n;
		  
		  /* read the header to determine image size */
		  sprintf(fn, "%s%s", inputDirName, de->d_name);
		  f = fopen(fn, "r");
		  if (f == NULL)
		    Error("Could not open image file %s", de->d_name);
		  if (!ReadHeader(f, &tc, &iw, &ih, &m))
		    {
		      fclose(f);
		      Error("Invalid image file header: %s\n", de->d_name);
		    }
		  if (c.type == '\0')
		    c.type = tc;
		  else if (tc != c.type)
		    {
		      fclose(f);
		      Error("Inconsistent mix of pgm and ppm image files: %s\n", de->d_name);
		    }
		  if (m != 255)
		    {
		      fclose(f);
		      Error("Image file %s does not have byte-size image components.\n", de->d_name);
		    }
		  fclose(f);
		}
	      if (n > maxN)
		maxN = n;
	    }
	}
      closedir(dir);
      if (maxN < 0)
	Error("No input image files found.\n");
      Log("minN = %d maxN = %d\n", minN, maxN);
    }

  if (c.maskName[0] != '\0')
    {
      /* check what masks are present */
      for (i = strlen(c.maskName)-1; i >= 0 && c.maskName[i] != '/'; --i) ;
      if (i >= 0)
	{
	  strncpy(maskDirName, c.maskName, i+1);
	  maskDirName[i+1] = '\0';
	  strcpy(maskPrefix, &c.maskName[i+1]);
	}
      else
	{
	  strcpy(maskDirName, "./");
	  strcpy(maskPrefix, c.maskName);
	}
      
      /* read the directory to look for files of the form nameNNNN.pbm */
      dir = opendir(maskDirName);
      if (dir == NULL)
	Error("Could not open directory %s for mask files\n");
      maskPrefixLen = strlen(maskPrefix);
      maskMinN = 1000000000;
      maskMaxN = -1;
      c.maskDigits = -1;
      while ((de = readdir(dir)) != NULL)
	{
	  if (strncmp(de->d_name, maskPrefix, maskPrefixLen) == 0 &&
	      (len = strlen(de->d_name)) > 4 &&
	      strcmp(&(de->d_name[len-4]), ".pbm") == 0)
	    {
	      if (sscanf(&(de->d_name[maskPrefixLen]), "%d%n",
			 &n, &digitLen) < 1)
		continue;
	      if (c.maskDigits < 0)
		c.maskDigits = digitLen;
	      else if (digitLen != c.maskDigits)
		Error("Inconsistent number of digits in mask file names: %s\n", de->d_name);
	      if (n < maskMinN)
		{
		  maskMinN = n;
		  
		  /* read the header to determine image size */
		  sprintf(fn, "%s%s", maskDirName, de->d_name);
		  f = fopen(fn, "r");
		  if (f == NULL)
		    Error("Could not open mask file %s", fn);
		  if (!ReadHeader(f, &tc, &iw, &ih, &m))
		    {
		      fclose(f);
		      Error("Invalid image file header: %s\n", fn);
		    }
		  if (tc != '4')
		    {
		      fclose(f);
		      Error("Bad .pbm header in mask file: %s\n", fn);
		    }
		  fclose(f);
		}
	      if (n > maskMaxN)
		maskMaxN = n;
	    }
	}
      closedir(dir);
      if (maskMaxN < 0)
	Error("No mask image files found.\n");
      if (maskMaxN < minN || maskMinN > maxN)
	Error("No corresponding mask image files found.\n");
    }
  else
    {
      c.maskDigits = -1;
    }

  if (afterFileName[0] != '\0')
    {
      if (stat(afterFileName, &statBuf) != 0)
	Error("Cannot access the -after file: %s\n", afterFileName);
      c.afterTime = (double) (statBuf.st_mtime);
    }
  else
    c.afterTime = 0.0;

  if (pairsFile[0] != '\0')
    {
      f = fopen(pairsFile, "r");
      if (f == NULL)
	Error("Could not open pairs file %s\n", pairsFile);
      while (fscanf(f, "%d %d", &imgn, &refn) == 2)
	{
	  if ((nPairs & 1023) == 0)
	    pairs = (int *) realloc(pairs, (nPairs + 1024) * 2 * sizeof(int));
	  pairs[2*nPairs] = imgn;
	  pairs[2*nPairs+1] = refn;
	  ++nPairs;
	}
      fclose(f);
      c.pairwiseNaming = 1;
    }

  /* check that output directories are writeable */ 
  sprintf(fn, "%sTEST.map", c.outputName);
  f = fopen(fn, "w");
  if (f == NULL)
    {
      for (i = strlen(c.outputName)-1; i >= 0 && c.outputName[i] != '/'; --i) ;
      if (i != 0)
	strncpy(outputDirName, c.outputName, i+1);
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

  par_set_context();

  /* for all slices */
  printf("Processing slices: ");
  fflush(stdout);
  if (pairs != 0)
    {
      for (pn = 0; pn < nPairs; ++pn)
	{
	  t.imageSlice = pairs[2*pn];
	  t.refSlice = pairs[2*pn+1];
	  if (t.imageSlice < minSlice || t.imageSlice > maxSlice ||
	      t.refSlice < minSlice || t.refSlice > maxSlice)
	    continue;
	  par_delegate_task();
	}
    }
  else if (strlen(c.inputName) > 4 &&
	   strcmp(&c.inputName[strlen(c.inputName) - 4], ".pgm") == 0)
    {
      t.imageSlice = 0;
      skipIndex = 0;
      for (t.refSlice = minN; t.refSlice <= maxN; ++t.refSlice)
	{
	  if (t.refSlice < minSlice || t.refSlice > maxSlice)
	    continue;
	  while (skipIndex < nSkipSlices && skipSlices[skipIndex] < t.refSlice)
	    ++skipIndex;
	  if (skipIndex < nSkipSlices && skipSlices[skipIndex] == t.refSlice)
	    continue;
	  par_delegate_task();
	}
    }
  else if (strlen(c.referenceName) > 4 &&
	   strcmp(&c.referenceName[strlen(c.referenceName) - 4], ".pgm") == 0)
    {
      t.refSlice = 0;
      skipIndex = 0;
      for (t.imageSlice = minN; t.imageSlice <= maxN; ++t.imageSlice)
	{
	  if (t.imageSlice < minSlice || t.imageSlice > maxSlice)
	    continue;
	  while (skipIndex < nSkipSlices && skipSlices[skipIndex] < t.imageSlice)
	    ++skipIndex;
	  if (skipIndex < nSkipSlices && skipSlices[skipIndex] == t.imageSlice)
	    continue;
	  par_delegate_task();
	}
    }
  else if (c.referenceName[0] != '\0')
    {
      skipIndex = 0;
      for (t.imageSlice = minN; t.imageSlice <= maxN; ++t.imageSlice)
	{
	  if (t.imageSlice < minSlice || t.imageSlice > maxSlice)
	    continue;
	  while (skipIndex < nSkipSlices && skipSlices[skipIndex] < t.imageSlice)
	    ++skipIndex;
	  if (skipIndex < nSkipSlices && skipSlices[skipIndex] == t.imageSlice)
	    continue;
	  t.refSlice = t.imageSlice;
	  par_delegate_task();
	}
    }
  else if (forward)
    {
      skipIndex = nSkipSlices - 1;
      t.refSlice = -1;
      for (t.imageSlice = maxN; t.imageSlice >= minN; --t.imageSlice)
	{
	  if (t.imageSlice < minSlice || t.imageSlice > maxSlice)
	    continue;
	  while (skipIndex >= 0 && skipSlices[skipIndex] > t.imageSlice)
	    --skipIndex;
	  if (skipIndex >= 0 && skipSlices[skipIndex] == t.imageSlice)
	    continue;
	  if (t.refSlice >= 0)
	    par_delegate_task();
	  t.refSlice = t.imageSlice;
	}
    }
  else
    {
      skipIndex = 0;
      t.refSlice = -1;
      for (t.imageSlice = minN; t.imageSlice <= maxN; ++t.imageSlice)
	{
	  if (t.imageSlice < minSlice || t.imageSlice > maxSlice)
	    continue;
	  while (skipIndex < nSkipSlices && skipSlices[skipIndex] < t.imageSlice)
	    ++skipIndex;
	  if (skipIndex < nSkipSlices && skipSlices[skipIndex] == t.imageSlice)
	    continue;
	  if (t.refSlice >= 0)
	    par_delegate_task();
	  t.refSlice = t.imageSlice;
	}
    }
  par_finish();

  sprintf(fn, "%squality.out", c.outputName);
  f = fopen(fn, "w");
  qsort(results, nResults, sizeof(Result), SortBySlice);
  fprintf(f, "Sorted by slice:\n");
  fprintf(f, "         IMAGE     REFERENCE   CORRELATION    DISTORTION    CORRESPOND        ENERGY\n");
  for (i = 0; i < nResults; ++i)
    fprintf(f, "%012.*d  %012.*d  %12.6f  %12.6f  %12.6f  %12.6f\n",
	    c.inputDigits, results[i].imageSlice,
	    c.inputDigits, results[i].refSlice,
	    results[i].correlation,
	    results[i].distortion,
	    results[i].correspondence,
	    results[i].distortion * c.distortion - results[i].correlation +
	    results[i].correspondence * c.correspondence);
  fprintf(f, "\n\nSorted by energy:\n");
  fprintf(f, "         IMAGE     REFERENCE   CORRELATION    DISTORTION    CORRESPOND        ENERGY\n");
  qsort(results, nResults, sizeof(Result), SortByEnergy);
  for (i = 0; i < nResults; ++i)
    fprintf(f, "%012.*d  %012.*d  %12.6f  %12.6f  %12.6f  %12.6f\n",
	    c.inputDigits, results[i].imageSlice,
	    c.inputDigits, results[i].refSlice,
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
  results[nResults].message = 0;

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
  if (rx->imageSlice < ry->imageSlice)
    return(-1);
  else if (rx->imageSlice == ry->imageSlice)
    return(0);
  else
    return(1);
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
  r.message = (char *) malloc(PATH_MAX + 1024);
  r.message[0] = '\0';
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
  size_t imbpl, rmbpl;
  size_t imbpl1, rmbpl1;
  float *warp;
  unsigned char *warpout;
  int level;
  unsigned int imw, imh;
  unsigned int rmw, rmh;
  int mpw, mph;
  int maskPresent;
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

  /* check if we can skip this task because of the -after option */
  update = 0;
  if (c.afterTime == 0.0)
    update = 1;
  if (!update)
    {
      if (strlen(c.referenceName) > 4 &&
	  strcmp(&c.referenceName[strlen(c.referenceName) - 4], ".pgm") == 0)
	strcpy(fn, c.referenceName);
      else
	sprintf(fn, "%s%0.*d.%s",
		(c.referenceName[0] != '\0') ? c.referenceName : c.inputName,
		c.inputDigits, t.refSlice,
		c.type == '5' ? "pgm" : "ppm");
      if (stat(fn, &sb) != 0 || (double) (sb.st_mtime) > c.afterTime)
	update = 1;
    }
  if (!update && (c.referenceName[0] == '\0' && c.maskName[0] != '\0' ||
		  c.referenceName[0] != '\0' && c.referenceMaskName[0] != '\0'))
    {
      if (c.referenceName[0] != '\0')
	sprintf(fn, "%s%0.*d.pbm", c.referenceMaskName, c.maskDigits, t.refSlice);
      else
	sprintf(fn, "%s%0.*d.pbm", c.maskName, c.maskDigits, t.refSlice);
      if (stat(fn, &sb) == 0 && (double) (sb.st_mtime) > c.afterTime)
	update = 1;
    }
  if (!update)
    {
      if (strlen(c.inputName) > 4 &&
	  strcmp(&c.inputName[strlen(c.inputName) - 4], ".pgm") == 0)
	strcpy(fn, c.inputName);
      else
	sprintf(fn, "%s%0.*d.%s", c.inputName, c.inputDigits, t.imageSlice,
		c.type == '5' ? "pgm" : "ppm");
      if (stat(fn, &sb) != 0 || (double) (sb.st_mtime) > c.afterTime)
	update = 1;
    }	  
  if (!update && c.maskName[0] != '\0')
    {
      sprintf(fn, "%s%0.*d.pbm", c.maskName, c.maskDigits, t.imageSlice);
      if (stat(fn, &sb) == 0 && (double) (sb.st_mtime) > c.afterTime)
	update = 1;
    }
  if (!update && c.cptsName[0] != '\0')
    {
      sprintf(fn, "%s%0.*d.pts", c.cptsName, c.inputDigits, t.imageSlice);
      if (stat(fn, &sb) == 0 && (double) (sb.st_mtime) > c.afterTime)
	{
	  f = fopen(fn, "r");
	  if (f != NULL)
	    {
	      while (fgets(line, LINE_LENGTH, f) != NULL)
		{
		  sscanf(line, "%d%d%d%d%d", &imx, &imy, &refSlice, &refx, &refy);
		  if (refSlice == t.refSlice)
		    {
		      update = 1;
		      break;
		    }
		}
	      fclose(f);
	    }
	}
    }
  if (!update)
    {
      r.refSlice = t.refSlice;
      r.imageSlice = t.imageSlice;
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

  constrainingMapX = NULL;
  constrainingMapY = NULL;
  constrainingMapC = NULL;

  /* get the reference image */
  if (strlen(c.referenceName) > 4 &&
      strcmp(&c.referenceName[strlen(c.referenceName) - 4], ".pgm") == 0)
    strcpy(fn, c.referenceName);
  else
    sprintf(fn, "%s%0.*d.%s",
	    (c.referenceName[0] != '\0') ? c.referenceName : c.inputName,
	    c.inputDigits, t.refSlice,
	    c.type == '5' ? "pgm" : "ppm");
  f = fopen(fn, "r");
  if (f == NULL)
    {
      sprintf(r.message, "Could not open image file %s\n", fn);
      return;
    }
  Log("WORKER OPENED %s\n", fn);
  if (!ReadHeader(f, &tc, &w, &h, &m))
    {
      sprintf(r.message, "Image file not binary pgm or ppm: %s\n", fn);
      fclose(f);
      return;
    }
  referenceWidth[0] = w;
  referenceHeight[0] = h;
  if (tc != c.type)
    {
      fclose(f);
      sprintf(r.message, "Inconsistent image file type: %s\n", fn);
      return;
    }
  if (c.type == '5')
    bp = 1;
  else
    bp = 3;
  referencePixels = ((size_t) referenceWidth[0]) * referenceHeight[0];
  ref_in = (unsigned char *) malloc(referencePixels*bp);
  if (fread(ref_in, 1, referencePixels*bp, f) != referencePixels*bp)
    {
      fclose(f);
      sprintf(r.message, "Image file %s apparently truncated\n", fn);
      return;
    }
  fclose(f);
  Log("WORKER READ THE REFERENCE\n");

#if MASKING
  maskPresent = 0;
  if (c.referenceName[0] == '\0' && c.maskName[0] != '\0' ||
      c.referenceName[0] != '\0' && c.referenceMaskName[0] != '\0')
    {
      /* get the mask for the reference image */
      if (c.referenceName[0] != '\0')
	{
	  if (strlen(c.referenceMaskName) > 4 &&
	      strcmp(&c.referenceMaskName[strlen(c.referenceMaskName) - 4], ".pbm") == 0)
	    strcpy(fn, c.referenceMaskName);
	  else
	    sprintf(fn, "%s%0.*d.pbm", c.referenceMaskName, c.maskDigits, t.refSlice);
	}
      else
	{
	  if (strlen(c.maskName) > 4 &&
	      strcmp(&c.maskName[strlen(c.maskName) - 4], ".pbm") == 0)
	    strcpy(fn, c.maskName);
	  else
	    sprintf(fn, "%s%0.*d.pbm", c.maskName, c.maskDigits, t.refSlice);
	}
      f = fopen(fn, "r");
      if (f != NULL)
	{
	  if (!ReadHeader(f, &tc, &rmw, &rmh, &m) || tc != '4')
	    {
	      sprintf(r.message, "Mask file not binary pbm: %s\n", fn);
	      fclose(f);
	      return;
	    }
	  Log("WORKER READ THE MASK HEADER\n");
	  referenceMaskWidth = rmw;
	  referenceMaskHeight = rmh;
	  referenceMaskFactor = referenceWidth[0] / rmw;
	  if (referenceMaskFactor * rmw != referenceWidth[0] ||
	      referenceMaskFactor * rmh != referenceHeight[0])
	    {
	      sprintf(r.message, "Size of reference is not a multiple of reference mask: %s\n", fn);
	      fclose(f);
	      return;
	    }
	  rmbpl = (rmw + 7) >> 3;
	  ref_mask_in = (unsigned char *) malloc(rmh * rmbpl);
	  if (fread(ref_mask_in, 1, rmh * rmbpl, f) != rmh * rmbpl)
	    {
	      fclose(f);
	      sprintf(r.message, "Mask file %s apparently truncated\n", fn);
	      return;
	    }
	  fclose(f);
	  Log("WORKER READ THE REFERENCE MASK\n");
	  
	  /* apply the mask */
	  if (referenceMaskFactor == 1)
	    {
	      for (my = 0; my < rmh; ++my)
		for (mx = 0; mx < rmw; ++mx)
		  if (!(ref_mask_in[my*rmbpl+(mx>>3)] & (0x80 >> (mx & 7))))
		    if (bp == 1)
		      ref_in[my*((size_t) rmw) + mx] = 0;
		    else
		      memset(&ref_in[bp*(my*((size_t) rmw)+mx)], 0, bp);
	    }
	  else
	    for (my = 0; my < rmh; ++my)
	      for (mx = 0; mx < rmw; ++mx)
		if (!(ref_mask_in[my*rmbpl+(mx>>3)] & (0x80 >> (mx & 7))))
		  for (y = 0; y < referenceMaskFactor; ++y)
		    memset(&ref_in[(my*referenceMaskFactor*((size_t) rmw)+mx)*
				   referenceMaskFactor*bp],
			   0, referenceMaskFactor*bp);
	  Log("WORKER APPLIED THE REFERENCE MASK\n");
	  maskPresent = 1;
	}
    }
  if (!maskPresent)
    {
      rmbpl = (referenceWidth[0] + 7) >> 3;
      ref_mask_in = (unsigned char *) malloc(referenceHeight[0]*rmbpl);
      memset(ref_mask_in, 0xff, referenceHeight[0]*rmbpl);
      referenceMaskFactor = 1;
      referenceMaskWidth = referenceWidth[0];
      referenceMaskHeight = referenceHeight[0];
    }
#endif

  /* read in the image to warp */
  if (strlen(c.inputName) > 4 &&
      strcmp(&c.inputName[strlen(c.inputName) - 4], ".pgm") == 0)
    strcpy(fn, c.inputName);
  else
    sprintf(fn, "%s%0.*d.%s", c.inputName, c.inputDigits, t.imageSlice,
	    c.type == '5' ? "pgm" : "ppm");
  f = fopen(fn, "r");
  if (f == NULL)
    {
      sprintf(r.message, "Could not open image file %s\n", fn);
      return;
    }
  if (!ReadHeader(f, &tc, &iw, &ih, &m))
    {
      sprintf(r.message, "Image file not binary pgm or ppm: %s\n", fn);
      fclose(f);
      return;
    }
  imageWidth[0] = iw;
  imageHeight[0] = ih;
  imagePixels = ((size_t) iw) * ih;
  image_in = (unsigned char *) malloc(imagePixels*bp);
  if (fread(image_in, 1, imagePixels*bp, f) != imagePixels*bp)
    {
      fclose(f);
      sprintf(r.message, "Image file %s apparently truncated\n", fn);
      return;
    }
  fclose(f);
  Log("WORKER READ THE IMAGE\n");

#if MASKING
  maskPresent = 0;
  if (c.maskName[0] != '\0')
    {
      /* get the mask for the image */
      if (strlen(c.maskName) > 4 &&
	   strcmp(&c.maskName[strlen(c.maskName) - 4], ".pbm") == 0)
	strcpy(fn, c.maskName);
      else
	sprintf(fn, "%s%0.*d.pbm", c.maskName, c.maskDigits, t.imageSlice);
      f = fopen(fn, "r");
      if (f != NULL)
	{
	  if (!ReadHeader(f, &tc, &imw, &imh, &m) || tc != '4')
	    {
	      sprintf(r.message, "Mask file not binary pbm: %s\n", fn);
	      fclose(f);
	      return;
	    }
	  imageMaskWidth = imw;
	  imageMaskHeight = imh;
	  imageMaskFactor = imageWidth[0] / imw;
	  if (imageMaskFactor * imw != imageWidth[0] ||
	      imageMaskFactor * imh != imageHeight[0])
	    {
	      sprintf(r.message, "Size of image is not a multiple of image mask: %s\n", fn);
	      fclose(f);
	      return;
	    }
	  imbpl = (imw + 7) >> 3;
	  image_mask_in = (unsigned char *) malloc(imh*imbpl);
	  if (fread(image_mask_in, 1, imh*imbpl, f) != imh*imbpl)
	    {
	      fclose(f);
	      sprintf(r.message, "Mask file %s apparently truncated\n", fn);
	      return;
	    }
	  fclose(f);
	  Log("WORKER READ THE IMAGE MASK\n");

	  /* apply the mask */
	  if (imageMaskFactor == 1)
	    {
	      for (my = 0; my < imh; ++my)
		for (mx = 0; mx < imw; ++mx)
		  if (!(image_mask_in[my*imbpl+(mx>>3)] & (0x80 >> (mx & 7))))
		    if (bp == 1)
		      image_in[my*((size_t) imw)+mx] = 0;
		    else
		      memset(&image_in[bp*(my*((size_t) imw)+mx)], 0, bp);
	    }
	  else
	    for (my = 0; my < imh; ++my)
	      for (mx = 0; mx < imw; ++mx)
		if (!(image_mask_in[my*imbpl+(mx>>3)] & (0x80 >> (mx & 7))))
		  for (y = 0; y < imageMaskFactor; ++y)
		    memset(&image_in[(my*imageMaskFactor*((size_t) imw)+mx)*
				     imageMaskFactor*bp],
			   0, imageMaskFactor*bp);
	  Log("WORKER APPLIED THE IMAGE MASK\n");
	  maskPresent = 1;
	}
    }
  if (!maskPresent)
    {
      imbpl = (imageWidth[0] + 7) >> 3;
      image_mask_in = (unsigned char *) malloc(imageHeight[0]*imbpl);
      memset(image_mask_in, 0xff, imageHeight[0]*imbpl);
      imageMaskFactor = 1;
      imageMaskWidth = imageWidth[0];
      imageMaskHeight = imageHeight[0];
    }
#endif

  /* read in the correspondence points */
  if (cpts != NULL)
    {
      nCpts = 0;
      free(cpts);
      cpts = NULL;
    }
  if (c.cptsName[0] != '\0')
    {
      sprintf(fn, "%s%0.*d.pts", c.cptsName, c.inputDigits, t.imageSlice);
      f = fopen(fn, "r");
      if (f != NULL)
	{
	  while (fgets(line, LINE_LENGTH, f) != NULL)
	    {
	      sscanf(line, "%d%d%d%d%d", &imx, &imy, &refSlice, &refx, &refy);
	      if (refSlice != t.refSlice)
		continue;
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


  images[0] = (float *) malloc(imagePixels * sizeof(float));
  for (ii = 0; ii < imagePixels; ++ii)
    images[0][ii] = image_in[ii];
  free(image_in);
  refs[0] = (float *) malloc(referencePixels * sizeof(float));
  for (ii = 0; ii < referencePixels; ++ii)
    refs[0][ii] = ref_in[ii];
  free(ref_in);

  /* get the input map, if present */
  if (c.inputMapName[0] != '\0')
    {
      sprintf(fn, "%s%0.*d.map", c.inputMapName, c.inputDigits, t.imageSlice);
      f = fopen(fn, "r");
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
	  if (mpwi > w)
	    Error("Input map resolution cannot be greater than image resolution\n");
	  inputMapFactor = 1 << ((int) floor(log(((float) w) / mpwi) / log(2.0) + 0.5));
	  Log("Read map %s of size %d x %d;  inputMapFactor = %d\n", fn, mpwi, mphi, inputMapFactor);
	}
      else
	{
	  Log("No input map found for slice %0.*d\n", c.inputDigits, t.imageSlice);
	  inputMapFactor = 0;
	}
    }
  else
    inputMapFactor = 0;

  /* get the constraining map, if present */
  if (c.constrainingMapName[0] != '\0')
    {
      sprintf(fn, "%s%0.*d.map", c.constrainingMapName, c.inputDigits, t.imageSlice);
      f = fopen(fn, "r");
      if (f != NULL)
	{
	  if (fread(mapHeader, sizeof(int), 3, f) != 3)
	    Error("Could not read map header from file %s", fn);
	  mpwc = mapHeader[0];
	  mphc = mapHeader[1];
	  mpwcp = mpwc + 2;
	  mphcp = mphc + 2;
	  constrainingMapX = (float *) malloc(mpwcp * mphcp * sizeof(float));
	  constrainingMapY = (float *) malloc(mpwcp * mphcp * sizeof(float));
	  constrainingMapC = (float *) malloc(mpwcp * mphcp * sizeof(float));
	  if (fread(constrainingMapX, sizeof(float), mpwcp * mphcp, f) != mpwcp * mphcp ||
	      fread(constrainingMapY, sizeof(float), mpwcp * mphcp, f) != mpwcp * mphcp ||
	      fread(constrainingMapC, sizeof(float), mpwcp * mphcp, f) != mpwcp * mphcp)
	    Error("Could not read map from file %s\n", fn);
	  fclose(f);
	  if (mpwc > w)
	    Error("Constraining map resolution cannot be greater than image resolution\n");
	  constrainingMapFactor = 1 << ((int) floor(log(((float) w) / mpwc) / log(2.0) + 0.5));
	  Log("Read map %s of size %d x %d;  constrainingMapFactor = %d\n", fn, mpwc, mphc, constrainingMapFactor);
	}
      else
	{
	  Log("No constraining map found for slice %0.*d\n", c.inputDigits, t.imageSlice);
	  constrainingMapFactor = 0;
	}
    }
  else
    constrainingMapFactor = 0;

  Init();
  Compute();

  if (c.outputImagesName[0] != '\0')
    {
      /* write out the warped image */
      warp = warped[c.outputLevel];
      iw = imageWidth[c.outputLevel];
      ih = imageHeight[c.outputLevel];
      imagePixels = ((size_t) iw) * ih;
      warpout = (unsigned char*) malloc(imagePixels);
      if (c.pairwiseNaming)
	sprintf(fn, "%s%0.*d-%0.*d.pgm", c.outputImagesName,
		c.inputDigits, t.imageSlice,
		c.inputDigits, t.refSlice);
      else
	sprintf(fn, "%s%0.*d.pgm", c.outputImagesName,
		c.inputDigits, t.imageSlice);
      f = fopen(fn, "w");
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

  if (c.outputCorrelationName[0] != '\0')
    {
      /* write out the correlation image */
#if 0
      /* FIX */
      warp = warped[0];
      w = width[0];
      h = height[0];
      correlation = (unsigned char*) malloc(w*h);
      if (c.pairwiseNaming)
	sprintf(fn, "%s%0.*d-%0.*d.pgm", c.outputCorrelationName,
		c.inputDigits, t.imageSlice,
		c.inputDigits, t.refSlice);
      else
	sprintf(fn, "%s%0.*d.pgm", c.outputCorrelationName,
		c.inputDigits, t.imageSlice);
      f = fopen(fn, "w");
      if (f == NULL)
	{
	  sprintf(r.message, "Could not open image file %s for writing.\n", fn);
	  return;
	}
      fprintf(f, "P5\n%d %d\n255\n", w, h);
      for (y = 0; y < h; ++y)
	for (x = 0; x < w; ++x)
	  correlation[y*w+x] = FIX;
      if (fwrite(correlation, 1, nPixels*sizeof(float), f) != nPixels*sizeof(float))
	{
	  fclose(f);
	  sprintf(r.message, "Could not write to image file %s\n", fn);
	  return;
	}
      fclose(f);
      free(correlation);
#endif
    }

  if (constrainingMapFactor != 0)
    {
      free(constrainingMapX);
      free(constrainingMapY);
      free(constrainingMapC);
    }
#if MASKING
  free(ref_mask_in);
  free(image_mask_in);
#endif
  for (level = 0; level < nLevels; ++level)
    {
      free(images[level]);
      free(refs[level]);
#if MASKING
      free(imask[level]);
      free(rmask[level]);
      free(eimask[level]);
      free(ermask[level]);
#endif
      free(xmaps[level]);
      free(ymaps[level]);
    }
}

void
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
  size_t imbpl, rmbpl;
  size_t imbpl1, rmbpl1;
  unsigned char *new_mask;
  int b;
  int dx, dy;
  unsigned char *iCount, *rCount;

  /* make the full resolution warped image */
  imagePixels = ((size_t) imageWidth[0]) * imageHeight[0];
  if (c.outputImagesName[0] != '\0')
    warped[0] = (float *) malloc(imagePixels * sizeof(float));

#if 0
  /* make the full resolution image significance array by taking the standard
     deviation of the pixel values over a circular region centered on each pixel */
  isig = (float *) malloc(imagePixels * sizeof(float));
#endif    

#if MASKING
  /* make the full resolution mask */
  imbpl1 = (imageMaskWidth + 7) >> 3;
  rmbpl1 = (referenceMaskWidth + 7) >> 3;
  imbpl = (imageWidth[0] + 7) >> 3;
  rmbpl = (referenceWidth[0] + 7) >> 3;
  imask[0] = (unsigned char *) malloc(imageHeight[0] * imbpl);
  rmask[0] = (unsigned char *) malloc(referenceHeight[0] * rmbpl);
  memset(imask[0], 0, imageHeight[0] * imbpl);
  memset(rmask[0], 0, referenceHeight[0] * rmbpl);
  iw = imageWidth[0];
  ih = imageHeight[0];
  /*  printf("iw = %d ih = %d\n", iw, ih); */
  for (y = 0; y < ih; ++y)
    for (x = 0; x < iw; ++x)
      if (image_mask_in[(y / imageMaskFactor) * imbpl1 + ((x / imageMaskFactor) >> 3)] &
	  (0x80 >> ((x / imageMaskFactor) & 7)))
	imask[0][y * imbpl + (x >> 3)] |= 0x80 >> (x & 7);
  rw = referenceWidth[0];
  rh = referenceHeight[0];
  /*  printf("rw = %d rh = %d\n", rw, rh); */
  for (y = 0; y < rh; ++y)
    for (x = 0; x < rw; ++x)
      if (ref_mask_in[(y / referenceMaskFactor) * rmbpl1 + ((x / referenceMaskFactor) >> 3)] &
	  (0x80 >> ((x / referenceMaskFactor) & 7)))
	rmask[0][y * rmbpl + (x >> 3)] |= 0x80 >> (x & 7);
  Log("imageMaskFactor = %d %ld  referenceMaskFactor = %d %ld\n",
      imageMaskFactor, (long) imbpl,
      referenceMaskFactor, (long) rmbpl);
  Log("icount = %ld\n", CountBits(image_mask_in, ih * imbpl / imageMaskFactor));
  Log("ocount = %ld\n", CountBits(imask[0], ih * imbpl));
#endif

  /* create hierarchical representations of images */
  for (level = 1;
       imageWidth[level-1] > 1 || imageHeight[level-1] > 1 ||
	 referenceWidth[level-1] > 1 || referenceHeight[level-1] > 1;
       ++level)
    {
      iw = imageWidth[level] = (imageWidth[level-1] + 1) >> 1;
      ih = imageHeight[level] = (imageHeight[level-1] + 1) >> 1;
      rw = referenceWidth[level] = (referenceWidth[level-1] + 1) >> 1;
      rh = referenceHeight[level] = (referenceHeight[level-1] + 1) >> 1;
      imagePixels = ((size_t) imageWidth[level]) * imageHeight[level];
      referencePixels = ((size_t) referenceWidth[level]) * referenceHeight[level];
      imbpl = (iw + 7) >> 3;
      imbpl1 = (imageWidth[level-1] + 7) >> 3;
      rmbpl = (rw + 7) >> 3;
      rmbpl1 = (referenceWidth[level-1] + 7) >> 3;
      images[level] = (float*) malloc(imagePixels * sizeof(float));
      refs[level] = (float*) malloc(referencePixels * sizeof(float));
      if (c.outputImagesName[0] != '\0')
        warped[level] = (float*) malloc(imagePixels * sizeof(float));
#if MASKING
      imask[level] = (unsigned char*) malloc(ih*imbpl);
      rmask[level] = (unsigned char*) malloc(rh*rmbpl);
      iCount = (unsigned char*) malloc(imagePixels * sizeof(unsigned char));
      rCount = (unsigned char*) malloc(referencePixels * sizeof(unsigned char));
#endif
	
      image = images[level];
      ref = refs[level];
      for (ii = 0; ii < imagePixels; ++ii)
	image[ii] = 0.0;
      for (ii = 0; ii < referencePixels; ++ii)
	ref[ii] = 0.0;
#if MASKING
      memset(iCount, 0, imagePixels*sizeof(unsigned char));
      memset(rCount, 0, referencePixels*sizeof(unsigned char));
#endif

      ii = 0;
      for (y = 0; y < imageHeight[level-1]; ++y)
	for (x = 0; x < imageWidth[level-1]; ++x)
	  {
	    image[(y >> 1) * ((size_t) iw) + (x >> 1)] += images[level-1][ii];
	    ++ii;
#if MASKING
	    if ((imask[level-1][y*imbpl1 + (x >> 3)] & (0x80 >> (x & 7))) != 0)
	      ++iCount[(y >> 1) * ((size_t) iw) + (x >> 1)];
#endif
	  }

      ii = 0;
      for (y = 0; y < referenceHeight[level-1]; ++y)
	for (x = 0; x < referenceWidth[level-1]; ++x)
	  {
	    ref[(y >> 1) * ((size_t) rw) + (x >> 1)] += refs[level-1][ii];
	    ++ii;
#if MASKING
	    if ((rmask[level-1][y*rmbpl1 + (x >> 3)] & (0x80 >> (x & 7))) != 0)
	      ++rCount[(y >> 1) * ((size_t) rw) + (x >> 1)];
#endif
	  }

      for (ii = 0; ii < imagePixels; ++ii)
	image[ii] *= 1.0 / 4.0;
      for (ii = 0; ii < referencePixels; ++ii)
	ref[ii] *= 1.0 / 4.0;

#if MASKING
      memset(imask[level], 0, ih*imbpl);
      memset(rmask[level], 0, rh*rmbpl);
      for (y = 0; y < ih; ++y)
	for (x = 0; x < iw; ++x)
	  if (iCount[y * ((size_t) iw) + x] == 4)
	    imask[level][y * imbpl + (x >> 3)] |= 0x80 >> (x & 7);
      for (y = 0; y < rh; ++y)
	for (x = 0; x < rw; ++x)
	  if (rCount[y * ((size_t) rw) + x] == 4)
	    rmask[level][y * rmbpl + (x >> 3)] |= 0x80 >> (x & 7);
      free(iCount);
      free(rCount);
#endif
    }
  nLevels = level;

  if (c.outputLevel < 0 || c.outputLevel >= nLevels)
    Error("outputLevel (%d) is invalid for image (nLevels = %d)\n",
	  c.outputLevel, nLevels);

#if MASKING
  for (level = 0; level < nLevels; ++level)
    {
      iw = imageWidth[level];
      ih = imageHeight[level];
      rw = referenceWidth[level];
      rh = referenceHeight[level];
      imagePixels = ((size_t) iw) * ih;
      referencePixels = ((size_t) rw) * rh;
      imbpl = (iw + 7) >> 3;
      rmbpl = (rw + 7) >> 3;
      eimask[level] = (unsigned char *) malloc(ih * imbpl);
      ermask[level] = (unsigned char *) malloc(rh * rmbpl);

      memcpy(eimask[level], imask[level], ih * imbpl);
      memcpy(ermask[level], rmask[level], rh * rmbpl);

#if 0
     //  this code is disabled since we are not using the extended masks
     //      in Compute
      new_mask = (unsigned char*) malloc(h * mbpl);
      memset(new_mask, 0, h * mbpl);
      for (i = 0; i < 17; ++i)
	{
	  for (y = 0; y < h; ++y)
	    for (x = 0; x < w; ++x)
	      {
		b = 0;
		for (dy = -1; dy <= 1; ++dy)
		  if (y + dy >= 0 && y + dy < h)
		    {
		      for (dx = -1; dx <= 1; ++dx)
			if (x + dx >= 0 && x + dx < w &&
			    (eimask[level][(y + dy)*mbpl + ((x + dx) >> 3)] &
			     (0x80 >> ((x + dx) & 7))) != 0)
			  {
			    b = 1;
			    break;
			  }
		      if (b)
			break;
		    }
		if (b)
		  new_mask[y*mbpl + (x >> 3)] |= 0x80 >> (x & 7);
	      }
	  memcpy(eimask[level], new_mask, h*mbpl);

	  memset(new_mask, 0, h * mbpl);
	  for (y = 0; y < h; ++y)
	    for (x = 0; x < w; ++x)
	      {
		b = 0;
		for (dy = -1; dy <= 1; ++dy)
		  if (y + dy >= 0 && y + dy < h)
		    {
		      for (dx = -1; dx <= 1; ++dx)
			if (x + dx >= 0 && x + dx < w &&
			    (ermask[level][(y + dy)*mbpl + ((x + dx) >> 3)] &
			     (0x80 >> ((x + dx) & 7))) != 0)
			  {
			    b = 1;
			    break;
			  }
		      if (b)
			break;
		    }
		if (b)
		  new_mask[y*mbpl + (x >> 3)] |= 0x80 >> (x & 7);
	      }
	  memcpy(ermask[level], new_mask, h*mbpl);
	}
      free(new_mask);
#endif

      imaskCount[level] = CountBits(imask[level], ih*imbpl);
      rmaskCount[level] = CountBits(rmask[level], rh*rmbpl);
      minMaskCount[level] = imaskCount[level];
      if (rmaskCount[level] < minMaskCount[level])
	minMaskCount[level] = rmaskCount[level];
      eimaskCount[level] = CountBits(eimask[level], ih*imbpl);
      ermaskCount[level] = CountBits(ermask[level], rh*rmbpl);
      Log("imaskCount[%d] = %ld (%dx%d = %ld) %f%%\n",
	  level, imaskCount[level],
	  iw, ih, imagePixels, 100.0 * ((float) imaskCount[level]) / imagePixels);
      Log("rmaskCount[%d] = %ld (%dx%d = %ld) %f%%\n",
	  level, rmaskCount[level],
	  rw, rh, referencePixels, 100.0 * ((float) rmaskCount[level]) / referencePixels);
      Log("eimaskCount[%d] = %ld (%dx%d = %ld) %f%%\n",
	  level, eimaskCount[level],
	  iw, ih, imagePixels, 100.0 * ((float) eimaskCount[level]) / imagePixels);
      Log("ermaskCount[%d] = %ld (%dx%d = %ld) %f%%\n",
	  level, ermaskCount[level],
	  rw, rh, referencePixels, 100.0 * ((float) ermaskCount[level]) / referencePixels);
    }
#endif

  Log("Constructed the multi-resolution image set.\n");
}
  
void
Compute ()
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
  int nMoves;
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
  int ixv, iyv;
  double xv, yv;
#if CHECKING
  float *bugOldRV;
  float *bugNewRV;
  unsigned char *bugChecked;
#endif
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
  int rotationalMove;
  double thetaRange;
  size_t statRotationalTheta[21];
  size_t statRotationalShiftX[21];
  size_t statRotationalShiftY[21];
  size_t nRotationalMoves;
  size_t statSqrtRadius[21];
  size_t statSqrtOffset[21];
  size_t statTheta[21];
  size_t statX[21];
  size_t statY[21];
  size_t nShiftMoves;
  float *xMapForced, *yMapForced, *cMapForced;
  int ix0, iy0;
  double rc;
  double rc00, rc01, rc10, rc11;
  unsigned int refw, refh;
  int multiplier;
  int startLevel;

  minOverlap = 0.5;
  kFactor = 1.0 / sqrt((double) (((size_t) referenceWidth[0]) * referenceWidth[0] +
	                         ((size_t) referenceHeight[0]) * referenceHeight[0]));
  xMapForced = NULL;
  yMapForced = NULL;
  cMapForced = NULL;
#if CHECKING
  bugOldRV = NULL;
  bugNewRV = NULL;
  bugChecked = NULL;
#endif
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
      imbpl = (mpw + 7) >> 3;
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
#if MASKING
      cimask = imask[cLevel];
      crmask = rmask[cLevel];
#endif
      icmbpl = (iw + 7) >> 3;
      rcmbpl = (rw + 7) >> 3;
      cimage = images[cLevel];
      cref = refs[cLevel];
      lFactor = (1 << level);

      mSize = ((size_t) mpwp) * mphp;
      xmaps[level] = (float*) malloc(mSize * sizeof(float));
      ymaps[level] = (float*) malloc(mSize * sizeof(float));
      xmap = xmaps[level];
      ymap = ymaps[level];
#if MASKING
      mask = imask[level];
#endif
      if (constrainingMapFactor != 0)
	{
	  if (xMapForced != NULL)
	    free(xMapForced);
	  if (yMapForced != NULL)
	    free(yMapForced);
	  if (cMapForced != NULL)
	    free(cMapForced);
	  xMapForced = (float *) malloc(mSize * sizeof(float));
	  yMapForced = (float *) malloc(mSize * sizeof(float));
	  cMapForced = (float *) malloc(mSize * sizeof(float));
	}

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
#if 0
	  printf("MAP AT LEVEL %d (imf = %d)\n", level, constrainingMapFactor);
	  for (y = -1; y <= mph; ++y)
	    for (x = -1; x <= mpw; ++x)
	      printf("xmap(%d, %d) = %f  ymap(%d, %d) = %f\n",
		     x, y, MAP(xmap, x, y), x, y, MAP(ymap, x, y));
	  printf("\n");
#endif	  

	}

      /* force the map if a constraining map was given */
      if (constrainingMapFactor != 0)
	for (y = -1; y <= mph; ++y)
	  for (x = -1; x <= mpw; ++x)
	    {
	      rx = factor * x + (factor - 1.0) / 2.0;
	      ry = factor * y + (factor - 1.0) / 2.0;
	      x0 = (2.0 * rx - constrainingMapFactor + 1) / (2.0 * constrainingMapFactor);
	      y0 = (2.0 * ry - constrainingMapFactor + 1) / (2.0 * constrainingMapFactor);
	      ix0 = (int) floor(x0);
	      iy0 = (int) floor(y0);
	      rrx = x0 - ix0;
	      rry = y0 - iy0;
	      if (ix0 < -1 || ix0 > mpwc || ix0 == mpwc && rrx > 0.0 ||
		  iy0 < -1 || iy0 > mphc || iy0 == mphc && rry > 0.0)
		{
		  rc = 0.0;
		  rx = 0.0;
		  ry = 0.0;
		}
	      else
		{
		  rx00 = MAPCON(constrainingMapX, ix0, iy0);
		  ry00 = MAPCON(constrainingMapY, ix0, iy0);
		  rc00 = MAPCON(constrainingMapC, ix0, iy0);
		  if (iy0 < mphc)
		    {
		      rx01 = MAPCON(constrainingMapX, ix0, iy0 + 1);
		      ry01 = MAPCON(constrainingMapY, ix0, iy0 + 1);
		      rc01 = MAPCON(constrainingMapC, ix0, iy0 + 1);
		    }
		  else
		    {
		      rx01 = 0.0;
		      ry01 = 0.0;
		      rc01 = 0.0;
		    }
		  if (ix0 < mpwc)
		    {
		      rx10 = MAPCON(constrainingMapX, ix0 + 1, iy0);
		      ry10 = MAPCON(constrainingMapY, ix0 + 1, iy0);
		      rc10 = MAPCON(constrainingMapC, ix0 + 1, iy0);
		    }
		  else
		    {
		      rx10 = 0.0;
		      ry10 = 0.0;
		      rc10 = 0.0;
		    }
		  if (ix0 < mpwc && iy0 < mphc)
		    {
		      rx11 = MAPCON(constrainingMapX, ix0 + 1, iy0 + 1);
		      ry11 = MAPCON(constrainingMapY, ix0 + 1, iy0 + 1);
		      rc11 = MAPCON(constrainingMapC, ix0 + 1, iy0 + 1);
		    }
		  else
		    {
		      rx11 = 0.0;
		      ry11 = 0.0;
		      rc11 = 0.0;
		    }
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
		}

	      x0 = constrainingMapFactor * rx + (constrainingMapFactor - 1.0) / 2.0;
	      y0 = constrainingMapFactor * ry + (constrainingMapFactor - 1.0) / 2.0;
	      rx = (2.0 * x0 - factor + 1) / (2.0 * factor);
	      ry = (2.0 * y0 - factor + 1) / (2.0 * factor);

	      MAP(xMapForced, x, y) = rx;
	      MAP(yMapForced, x, y) = ry;
	      MAP(cMapForced, x, y) = rc;
	      MAP(xmap, x, y) = rc * rx + (1.0 - rc) * MAP(xmap, x, y);
	      MAP(ymap, x, y) = rc * ry + (1.0 - rc) * MAP(ymap, x, y);
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
#if CHECKING
      if (bugOldRV != NULL)
	free(bugOldRV);
      if (bugNewRV != NULL)
	free(bugNewRV);
      if (bugChecked != NULL)
	free(bugChecked);
      bugOldRV = (float *) malloc(((size_t) iw) * ih * sizeof(float));
      bugNewRV = (float *) malloc(((size_t) iw) * ih * sizeof(float));
      bugChecked = (unsigned char *) malloc(((size_t) iw) * ih);
#endif
      for (y = 0; y < ih; ++y)
	for (x = 0; x < iw; ++x)
	  {
#if CHECKING
	    bugNewRV[y*((size_t) iw) + x] = 0.0;
#endif
#if MASKING
	    if ((cimask[y*icmbpl + (x >> 3)] & (0x80 >> (x & 7))) == 0)
	      continue;
#endif
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
#if MASKING
	    if ((crmask[iry*rcmbpl + (irx >> 3)] & (0x80 >> (irx & 7))) == 0 ||
		rrx > 0.0 && (crmask[iry*rcmbpl + ((irx+1) >> 3)] & (0x80 >> ((irx+1) & 7))) == 0 ||
		rry > 0.0 && (crmask[(iry+1)*rcmbpl + (irx >> 3)] & (0x80 >> (irx & 7))) == 0 ||
		rrx > 0.0 && rry > 0.0 && (crmask[(iry+1)*rcmbpl +((irx+1) >> 3)] & (0x80 >> ((irx+1) & 7))) == 0)
	      continue;
#endif
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
#if CHECKING	    
	    bugNewRV[y*iw + x] = rv;
#endif
	    iv = IIMAGE(cimage, x, y);
	    si += iv;
	    si2 += iv * iv;
	    sr += rv;
	    sr2 += rv * rv;
	    sir += iv * rv;
	  }
#if MASKING
      requiredPoints = (size_t) ceil(minMaskCount[cLevel] * minOverlap);
#else
      requiredPoints = ((size_t) iw) * ih * minOverlap;
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
      de = 0.0;
      for (y = -1; y < mph; ++y)
	for (x = -1; x < mpw; ++x)
	  {
#if 0
	    if ((mask[y*imbpl + (x >> 3)] & (0x80 >> (x & 7))) == 0)
	      continue;
#endif
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
	  }
      distortion = de / (mph+1) / (mpw+1);

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

      sqrtMaxMoveRadius = sqrt(0.25 * (nLevels - level - 1) * (nLevels - level - 1) + 2.0);
      sqrtMaxMoveOffset = 0.3 * sqrt(0.01 * (nLevels - level - 1) * (nLevels - level - 1) + 0.25);
      multiplier = 1 << level;
      if (multiplier > 4096)
	multiplier = 4096;
      goalCoverage = (c.quality * mpw) * mph * multiplier;
      actualCoverage = 0.0;

      rawMoveCount = 0;
      moveCount = 0;
      acceptedMoveCount = 0;
      nRotationalMoves = 0;
      nShiftMoves = 0;

      memset(statRotationalTheta, 0, 21*sizeof(size_t));
      memset(statRotationalShiftX, 0, 21*sizeof(size_t));
      memset(statRotationalShiftY, 0, 21*sizeof(size_t));
      memset(statSqrtRadius, 0, 21*sizeof(size_t));
      memset(statSqrtOffset, 0, 21*sizeof(size_t));
      memset(statTheta, 0, 21*sizeof(size_t));
      memset(statX, 0, 21*sizeof(size_t));
      memset(statY, 0, 21*sizeof(size_t));

      while (actualCoverage < goalCoverage)
	{
	  ++rawMoveCount;

	  if (level >= nLevels - 2 && drand48() < 0.5)
	    {
	      /* try a rotational move */
	      rotationalMove = 1;
	      if (level == nLevels - 1)
		thetaRange = 2.0 * M_PI;  /* arbitrary rotation */
	      else
		thetaRange = 0.25 * M_PI; /* small rotation */
	      theta = (drand48() - 0.5) * thetaRange;
	      dxp = drand48() - 0.5;
	      dyp = drand48() - 0.5;

	      ++moveCount;
	      cth = cos(theta);
	      sth = sin(theta);
	      if (constrainingMapFactor != 0)
		{
		  changeMinX = mpw + 1;
		  changeMaxX = -2;
		  changeMinY = mph + 1;
		  changeMaxY = -2;
		  for (y = -1; y <= mph; ++y)
		    for (x = -1; x <= mpw; ++x)
		      {
			rc = MAP(cMapForced, x, y);
			if (rc >= 1.0)
			  continue;
			x1 = cth * MAP(xmap, x, y) +
			  sth * MAP(ymap, x, y) + dxp;
			y1 = -sth * MAP(xmap, x, y) +
			  cth * MAP(ymap, x, y) + dyp;
			MAP(xprop, x, y) = rc * MAP(xMapForced, x, y) + (1.0 - rc) * x1;
			MAP(yprop, x, y) = rc * MAP(yMapForced, x, y) + (1.0 - rc) * y1;
			MAP(changed, x, y) = 1;
			if (x < changeMinX)
			  changeMinX = x;
			if (x > changeMaxX)
			  changeMaxX = x;
			if (y < changeMinY)
			  changeMinY = y;
			if (y > changeMaxY)
			  changeMaxY = y;
		      }
		}
	      else
		{
		  for (y = -1; y <= mph; ++y)
		    for (x = -1; x <= mpw; ++x)
		      {
			MAP(xprop, x, y) = cth * MAP(xmap, x, y) +
			  sth * MAP(ymap, x, y) + dxp;
			MAP(yprop, x, y) = -sth * MAP(xmap, x, y) +
			  cth * MAP(ymap, x, y) + dyp;
			MAP(changed, x, y) = 1;
		      }
		  changeMinX = -1;
		  changeMaxX = mpw;
		  changeMinY = -1;
		  changeMaxY = mph;
		}
	    }
	  else
	    {
	      /* choose size of move to make */
	      rotationalMove = 0;
	      sqrtRadius = drand48() * sqrtMaxMoveRadius;
	      radius = sqrtRadius * sqrtRadius;
	      sqrtOffset = drand48() * sqrtMaxMoveOffset;
	      offset = sqrtOffset * sqrtOffset;
	      if (offset > 0.5 * radius)
		continue;
	      theta = drand48() * 2.0 * M_PI;

	      /* choose region in mapping to distort */
	      /*    center(x,y) radius(r) offset(x,y) */
	      icx = (int) floor(drand48() * mpw);
	      icy = (int) floor(drand48() * mph);
	      cx = MAP(xmap, icx, icy);
	      cy = MAP(ymap, icx, icy);

	      /*	  printf("Considering map shift %d at (%f, %f) of radius %f, offset %f,\n    and theta %f ... ", rawMoveCount, cx, cy, radius, offset, theta); */
	      if (M_PI * radius * radius < ((double) mpw) * mph)
		actualCoverage += M_PI * radius * radius;
	      else
		actualCoverage += ((double) mpw) * mph;
	      /*	  printf("  actualCoverage = %f  goalCoverage = %f\n",
			  actualCoverage, goalCoverage); */

#if 0
	      if ((eimask[level][icy*imbpl + (icx >> 3)] & (0x80 >> (icx & 7))) == 0)
		continue;
#endif
	      ++moveCount;

	      /* shift the map */
	      /*	  memcpy(xprop, xmap, hp * wp * sizeof(float));
			  memcpy(yprop, ymap, hp * wp * sizeof(float)); */

	      changeMinX = mpw + 1;
	      changeMaxX = -2;
	      changeMinY = mph + 1;
	      changeMaxY = -2;
	      cth = cos(theta);
	      sth = sin(theta);
	      n = 0;
	      for (xdir = -1; xdir <= 1; xdir += 2)
		for (ix = ((xdir-1) >> 1); ; ix += xdir)
		  {
		    x = icx + ix;
		    if (x < -1)
		      if (xdir == -1)
			break;
		      else
			continue;
		    if (x > mpw)
		      if (xdir == 1)
			break;
		      else
			continue;
		    for (ydir = -1; ydir <= 1; ydir += 2)
		      for (iy = ((ydir-1) >> 1); ; iy += ydir)
			{
			  y = icy + iy;
			  if (y < -1)
			    if (ydir == -1)
			      break;
			    else
			      continue;
			  if (y > mph)
			    if (ydir == 1)
			      break;
			    else
			      continue;

			  if (constrainingMapFactor != 0)
			    {
			      rc = MAP(cMapForced, x, y);
			      if (rc >= 1.0)
				continue;
			    }
			  else
			    rc = 0.0;

			  dx = MAP(xmap, x, y) - cx;
			  dy = MAP(ymap, x, y) - cy;
			  if (dx * dx + dy * dy >= 4.0 * radius * radius)
			    break;
			  if (dx * dx + dy * dy >= radius * radius)
			    continue;
			  MAP(changed, x, y) = 1;
			  ++n;
			  if (x < changeMinX)
			    changeMinX = x;
			  if (x > changeMaxX)
			    changeMaxX = x;
			  if (y < changeMinY)
			    changeMinY = y;
			  if (y > changeMaxY)
			    changeMaxY = y;
			  dxp = cth * dx + sth * dy;
			  dyp = -sth * dx + cth * dy;
			  mrd = (radius - sqrt(dxp * dxp + dyp * dyp)) / radius;
			  dxp += mrd * offset;
			  nx = cth * dxp - sth * dyp;
			  ny = sth * dxp + cth * dyp;
			  x1 = cx + nx;
			  y1 = cy + ny;
			  if (rc > 0.0)
			    {
			      MAP(xprop, x, y) = rc * MAP(xMapForced, x, y) + (1.0 - rc) * x1;
			      MAP(yprop, x, y) = rc * MAP(yMapForced, x, y) + (1.0 - rc) * y1;
			    }
			  else
			    {			  
			      MAP(xprop, x, y) = x1;
			      MAP(yprop, x, y) = y1;
			    }
			}
		  }
	      /*	  printf("(moving %d pts) ", n); */
	      /*	  printf("changeMin = %d %d %d %d\n",
			  changeMinX, changeMaxX, changeMinY, changeMaxY); */
	    }

	  /* evaluate effect of that move on energy */

	  /* add in contribution from correlation */
#if CHECKING
	  memcpy(bugOldRV, bugNewRV,
		 ((size_t) imageWidth[cLevel]) * imageHeight[cLevel] * sizeof(float));
	  memset(bugChecked, 0, ((size_t) iw) * ih);
#endif
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
	  /* TEMPORARY */
	  /*	  sx = 0;
	  ex = w - 1;
	  sy = 0;
	  ey = h - 1; */

	  /*	  printf("sx = %d %d %d %d  factor = %d\n", sx, ex, sy, ey, factor); */
	  for (y = sy; y <= ey; ++y)
            for (x = sx; x <= ex; ++x)
	      {
#if MASKING
		if ((cimask[y*icmbpl + (x >> 3)] & (0x80 >> (x & 7))) == 0)
		  continue;
#endif
		xv = (2.0 * x - factor + 1) / (2.0 * factor);
		yv = (2.0 * y - factor + 1) / (2.0 * factor);
		ixv = ((int) (xv + 2.0)) - 2;
		iyv = ((int) (yv + 2.0)) - 2;
		rrx = xv - ixv;
		rry = yv - iyv;
#if TRACKING
		if (x == 2 && y == 2)
		  printf("xv = %f yv = %f\n", xv, yv);
#endif
		if (ixv < -1 || ixv > mpw || ixv == mpw && rrx > 0.0 ||
		    iyv < -1 || iyv > mph || iyv == mph && rry > 0.0)
		  continue;
#if 1
		if (!MAP(changed, ixv, iyv) &&
		    !MAP(changed, ixv+1, iyv) &&
		    !MAP(changed, ixv, iyv+1) &&
		    !MAP(changed, ixv+1, iyv+1))
		    continue;
#endif
#if CHECKING
		bugChecked[y * iw + x] = 1;
#endif
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
#if TRACKING
		if (x == 2 && y == 2)
		  printf("1: rx = %f ry = %f\n", rx, ry);
#endif
		if (irx >= 0 && (irx < rw-1 || irx == rw-1 && rrx == 0.0) &&
		    iry >= 0 && (iry < rh-1 || iry == rh-1 && rry == 0.0))
		  {
#if MASKING
		    if ((crmask[iry*rcmbpl + (irx >> 3)] & (0x80 >> (irx & 7))) != 0 &&
			((rrx <= 0.0) || (crmask[iry*rcmbpl + ((irx+1) >> 3)] & (0x80 >> ((irx+1) & 7))) != 0) &&
			((rry <= 0.0) || (crmask[(iry+1)*rcmbpl + (irx >> 3)] & (0x80 >> (irx & 7))) != 0) &&
			((rrx <= 0.0) || (rry <= 0.0) || (crmask[(iry+1)*rcmbpl +((irx+1) >> 3)] & (0x80 >> ((irx+1) & 7))) != 0))
		      {
#endif
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
#if TRACKING
			if (x == 2 && y == 2)
			  printf("Subtracting rv = %f\n", rv);
#endif
#if CHECKING
			bugOldRV[y*iw+x] -= rv;
#endif
#if TRACKING
			if (x == 2 && y == 2)
			  printf("yielding %f\n", bugOldRV[y*iw+x]);
#endif
#if MASKING
		      }
#endif
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

#if TRACKING
		if (x == 2 && y == 2)
		  printf("2: rx = %f ry = %f\n", rx, ry);
#endif
		irx = ((int) (rx + 1.0)) - 1;
		iry = ((int) (ry + 1.0)) - 1;
		rrx = rx - irx;
		rry = ry - iry;
		if (irx < 0 || irx >= rw || irx == rw-1 && rrx > 0.0 ||
		    iry < 0 || iry >= rh || iry == rh-1 && rry > 0.0)
		  continue;
#if MASKING
		if ((crmask[iry*rcmbpl + (irx >> 3)] & (0x80 >> (irx & 7))) == 0 ||
		    rrx > 0.0 && (crmask[iry*rcmbpl + ((irx+1) >> 3)] & (0x80 >> ((irx+1) & 7))) == 0 ||
		    rry > 0.0 && (crmask[(iry+1)*rcmbpl + (irx >> 3)] & (0x80 >> (irx & 7))) == 0 ||
		    rrx > 0.0 && rry > 0.0 && (crmask[(iry+1)*rcmbpl +((irx+1) >> 3)] & (0x80 >> ((irx+1) & 7))) == 0)
		  continue;
#endif
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
#if TRACKING
		if (x == 2 && y == 2)
		  printf("Adding rv = %f\n", rv);
#endif
#if CHECKING
		bugOldRV[y*iw+x] += rv;
#endif
#if TRACKING
		if (x == 2 && y == 2)
		  printf("yielding %f\n", bugOldRV[y*iw+x]);
#endif
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
	  updcnt = 0;
          for (y = changeMinY - 1; y <= changeMaxY; ++y)
	    {
              if (y < -1 || y > mph)
		continue;
              for (x = changeMinX - 1; x <= changeMaxX; ++x)
      	        {
	          if (x < -1 || x > mpw)
	            continue;
#if 0
		  if ((mask[y*imbpl + (x >> 3)] & (0x80 >> (x & 7))) == 0)
		    continue;
#endif
		  if (x < mpw && y < mph)
		    {
		      upd0 = MAP(changed, x, y);
		      upd1 = MAP(changed, x, y+1);
		      upd2 = MAP(changed, x+1, y+1);
		      upd3 = MAP(changed, x+1, y);
		      
		      if (upd0 || upd1 || upd2 || upd3)
			{
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
#if DEBUG_MOVES
	      Log("Accepted move %d since deltaE = %f (dist %f - correl %f + corres %f)\n",
		  moveCount, newE - e,
		  (newDistortion - distortion) * c.distortion,
		  newCorrelation - correlation,
		  (newCorrespondence - correspondence) * c.correspondence);
	      Log("  cx = %f cy = %f radius = %f offset = %f theta = %f\n",
		  cx, cy, radius, offset, theta);
	      Log("  estimated new energy = %f (dist %f - correl %f + corres %f)\n",
		  newE, newDistortion * c.distortion, newCorrelation,
		  newCorrespondence * c.correspondence);
	      /*	      printf("  move area: %d %d %d %d\n", changeMinX, changeMinY,
			      changeMaxX, changeMaxY); */
#endif
	      if (rotationalMove)
		{
		  ++statRotationalTheta[(int) (20.0 * (theta / thetaRange + 0.5))];
		  ++statRotationalShiftX[(int) (20.0 * (dxp + 0.5))];
		  ++statRotationalShiftY[(int) (20.0 * (dyp + 0.5))];
		  ++nRotationalMoves;
		}
	      else
		{
		  ++statSqrtRadius[(int) (20.0 * (sqrtRadius / sqrtMaxMoveRadius))];
		  ++statSqrtOffset[(int) (20.0 * (sqrtOffset / sqrtMaxMoveOffset))];
		  ++statTheta[(int) (20.0 * theta / (2.0 * M_PI))];
		  ++statX[(int) ((20.0 * icx) / mpw)];
		  ++statY[(int) ((20.0 * icy) / mph)];
		  ++nShiftMoves;
		}

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

#if 0
	      /* recalculate correspondence energy */
	      ces = 0.0;
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
		      ces += 1000000.0;
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
		    ce = 0.0;
		  else
		    ce = distance;
		  ces += ce;
#if 0
		  printf("cci-Ap %d ix %d iy %d rx %d ry %d\n",
			 i, cpts[i].ix, cpts[i].iy, cpts[i].rx, cpts[i].ry);
		  printf("cci-Bp %d xv %f yv %f rx00 %f ry00 %f rx %f ry %f rrx %f rry %f dist %f en %f corr %f\n",
			 i, xv, yv, rx00, ry00, rx, ry, rrx, rry, distance, ce,
			 ces);
#endif
		}
	      if (fabs(ces - correspondence) > 0.000001)
		{
		  Log("CORRESPONDENCE FAILURE!!!  (ces = %f corres = %f\n",
		      ces, correspondence);
		}
#endif

#if CHECKING

      /* calculate energy of mapping */
      ea = 0.0;

      /* contribution from correlation */
      check_sir = 0.0;
      check_si2 = 0.0;
      check_sr2 = 0.0;
      check_si = 0.0;
      check_sr = 0.0;
      check_sum = 0.0;
      check_points = 0;
      for (y = 0; y < ih; ++y)
	for (x = 0; x < iw; ++x)
	  {
	    bugNewRV[y * iw + x] = 0.0;
#if MASKING
	    if ((cimask[y*icmbpl + (x >> 3)] & (0x80 >> (x & 7))) == 0)
	      continue;
#endif
	    /* use bilinear interpolation to find value */
	    xv = (2.0 * x - factor + 1) / (2.0 * factor);
	    yv = (2.0 * y - factor + 1) / (2.0 * factor);
	    ixv = ((int) (xv + 2.0)) - 2;
	    iyv = ((int) (yv + 2.0)) - 2;
	    rrx = xv - ixv;
	    rry = yv - iyv;
#if TRACKING
	    if (x == 2 && y == 2)
	      printf("xv = %f yv = %f\n", xv, yv);
#endif
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
#if TRACKING
	    if (x == 2 && y == 2)
	      printf("rx = %f ry = %f\n", rx, ry);
#endif
	    if (irx < 0 || irx >= rw || irx == rw-1 && rrx > 0.0 ||
		iry < 0 || iry >= rh || iry == rh-1 && rry > 0.0)
	      continue;

#if MASKING
	    if ((crmask[iry*rcmbpl + (irx >> 3)] & (0x80 >> (irx & 7))) == 0 ||
		rrx > 0.0 && (crmask[iry*rcmbpl + ((irx+1) >> 3)] & (0x80 >> ((irx+1) & 7))) == 0 ||
		rry > 0.0 && (crmask[(iry+1)*rcmbpl + (irx >> 3)] & (0x80 >> (irx & 7))) == 0 ||
		rrx > 0.0 && rry > 0.0 && (crmask[(iry+1)*rcmbpl +((irx+1) >> 3)] & (0x80 >> ((irx+1) & 7))) == 0)
	      continue;
#endif
	    ++check_points;
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
	      Error("ISNAN internal error:\n  rx = %f irx = %d ry = %f iry = %d rrx = %f rry = %f\n",
		    rx, irx, ry, iry, rrx, rry);
	    iv = IIMAGE(cimage, x, y);
	    bugNewRV[y * iw + x] = rv;
	    check_sr += rv;
	    check_si += iv;
	    check_si2 += iv * iv;
	    check_sr2 += rv * rv;
	    check_sir += iv * rv;
	  }
      if (check_points < requiredPoints)
	{
	  check_sr += (requiredPoints - check_points) * 255.0;
	  check_sr2 += (requiredPoints - check_points) * 255.0 * 255.0;
	  check_effectivePoints = requiredPoints;
	}
      else
	check_effectivePoints = check_points;
      if (check_effectivePoints != 0)
	{
	  check_mi = check_si / check_effectivePoints;
	  check_mr = check_sr / check_effectivePoints;
	  check_denom = (check_si2 - 2.0 * check_mi * check_si +
			 check_effectivePoints * check_mi * check_mi) *
	    (check_sr2 - 2.0 * check_mr * check_sr + check_effectivePoints *
	     check_mr * check_mr);
	  if (check_denom < 0.001)
	    check_r = -1000000.0;
	  else
	    check_r = (check_sir - check_mi * check_sr - check_mr * check_si +
		       check_effectivePoints * check_mi * check_mr) /
	      sqrt(check_denom);
	}
      else
	check_r = -1000000.0;
      aCorrelation = check_r;

      for (y = 0; y < ih; ++y)
	for (x = 0; x < iw; ++x)
	  {
	    if (fabs(bugNewRV[y * iw + x] - bugOldRV[y * iw + x]) > 0.1)
	      Error("Internal error: found disparate location (%d %d) old: %f new: %f checked: %d\n",
		    x, y, bugOldRV[y * iw + x], bugNewRV[y * iw + x],
		    bugChecked[y * iw + x]);
	  }

      /* contribution from distortion */
      de = 0.0;
      for (y = -1; y < mph; ++y)
	for (x = -1; x < mpw; ++x)
	  {
#if 0
	    if ((mask[y*imbpl + (x >> 3)] & (0x80 >> (x & 7))) == 0)
	      continue;
#endif
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
	  }
      aDistortion = de / (mph+1) / (mpw+1);
      ea = aDistortion * c.distortion - aCorrelation + aCorrespondence * c.correspondence;
#if DEBUG_MOVES
      Log("  actual new energy is %f (dist %f - correl %f + corres %f)\n",
	  ea, aDistortion * c.distortion, aCorrelation,
	  aCorrespondence * c.correspondence);
      Log("si = %f sr = %f\n", check_si, check_sr);
#endif

      if (fabs(check_si - si) > 0.000001 ||
	  fabs(check_sr - sr) > 0.000001)
	Error("CHECK_SI/SR ERROR\n");
      if (fabs(ea - e) > 0.000001)
	Error("Energy computation error\n");
      else
	{
#if DEBUG_MOVES
	  Log("Energy computation OK  %f %f %f\n", ea, e, fabs(ea-e));
#endif
	}
#endif
	    }
	}

      Log("Level %d: after %d raw moves, %d moves, and %d accepted moves...\n",
	  level, rawMoveCount, moveCount, acceptedMoveCount);
      Log("          energy is %f\n", e);
      if (level >= nLevels - 2)
	if (nRotationalMoves != 0)
	  {
            Log("Accepted rotational move statistics:\n");
	    Log("   range:  theta     shiftX     shiftY\n");
	    for (i = 0; i < 20; ++i)
	      Log(" %3d-%3d%%:  %10.5f  %10.5f  %10.5f\n",
		  5*i, 5*i+5,
		  (100.0 * statRotationalTheta[i]) / nRotationalMoves,
		  (100.0 * statRotationalShiftX[i]) / nRotationalMoves,
		  (100.0 * statRotationalShiftY[i]) / nRotationalMoves);
          }
	else
	  Log("No accepted rotational moves\n");
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
	WriteMap(level);
    }
  r.refSlice = t.refSlice;
  r.imageSlice = t.imageSlice;
  r.updated = 1;
  r.correlation = correlation;
  r.distortion = distortion;
  r.correspondence = correspondence;
  r.message[0] = '\0';
}

void
WriteMap (int level)
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
  size_t imbpl;
  int x, y;
  FILE *f;
  float rx, ry;
  int irx, iry;
  size_t mSize;
  unsigned int rw, rh;
  
  /* write out resultant map */
  if (level > c.outputLevel)
    if (c.pairwiseNaming)
      sprintf(fn, "%s%0.*d-%0.*d.%0.2d.map", c.outputName, c.inputDigits,
	      t.imageSlice, c.inputDigits, t.refSlice, level);
    else
      sprintf(fn, "%s%0.*d.%0.2d.map", c.outputName, c.inputDigits,
	      t.imageSlice, level);
  else
    if (c.pairwiseNaming)
      sprintf(fn, "%s%0.*d-%0.*d.map", c.outputName,
	      c.inputDigits, t.imageSlice, c.inputDigits, t.refSlice);
    else
      sprintf(fn, "%s%0.*d.map", c.outputName, c.inputDigits, t.imageSlice);
  f = fopen(fn, "w");
  if (f == NULL)
    Error("Could not open map file %s for writing\n", fn);
  mpw = imageWidth[level];
  mph = imageHeight[level];
  mpwp = mpw + 2;
  mphp = mph + 2;
  mSize = ((size_t) mphp) * mpwp;
  imbpl = (mpw + 7) >> 3;
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
#if MASKING
      for (x = -1; x <= mpw; ++x)
	if (y < 0 || y >= mph ||
	    x < 0 || x >= mpw ||
	    imask[level][y*imbpl+(x>>3)] & (0x80 >> (x & 7)) == 0)
	  MAP(cmap, x, y) = 0.0;
	else
	  {
	    rx = MAP(xmap, x, y);
	    ry = MAP(ymap, x, y);
	    irx = ((int) (rx + 1.5)) - 1;
	    iry = ((int) (ry + 1.5)) - 1;
	    if (irx >= 0 && irx < mpw && iry >= 0 && iry < mph)
	      {
		++cmapCount1;
	      }

	    if (irx >= 0 && irx < mpw && iry >= 0 && iry < mph &&
		(rmask[level][iry*imbpl + (irx>>3)] & (0x80 >> (irx & 7))) != 0)
	      ++cmapCount2;

	    if (irx >= 0 && irx < mpw &&
		iry >= 0 && iry < mph &&
		(rmask[level][iry*imbpl + (irx>>3)] & (0x80 >> (irx & 7))) != 0)
	      {
		MAP(cmap, x, y) = 1.0;
		++cmapCount3;
	      }
	    else
	      MAP(cmap, x, y) = 0.0;
	  }
#else
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
#endif
    }
  /* printf("CMAPCOUNT for %s is %ld %ld %ld\n", fn, cmapCount1,
     cmapCount2, cmapCount3); */

  fwrite(&mpw, sizeof(int), 1, f);
  fwrite(&mph, sizeof(int), 1, f);
  fwrite(&t.refSlice, sizeof(int), 1, f);
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
  par_pkstr(c.inputName);
  par_pkint(c.inputDigits);

  par_pkstr(c.maskName);
  par_pkint(c.maskDigits);

  par_pkstr(c.referenceName);
  par_pkstr(c.referenceMaskName);

  par_pkstr(c.cptsName);

  par_pkstr(c.constrainingMapName);

  par_pkint(c.pairwiseNaming);
  par_pkstr(c.outputName);
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
  par_upkstr(c.inputName);
  c.inputDigits = par_upkint();

  par_upkstr(c.maskName);
  c.maskDigits = par_upkint();

  par_upkstr(c.referenceName);
  par_upkstr(c.referenceMaskName);

  par_upkstr(c.cptsName);

  par_upkstr(c.constrainingMapName);

  c.pairwiseNaming = par_upkint();
  par_upkstr(c.outputName);
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
PackTask ()
{
  par_pkint(t.refSlice);
  par_pkint(t.imageSlice);
}

void
UnpackTask ()
{
  t.refSlice = par_upkint();
  t.imageSlice = par_upkint();
}

void
PackResult ()
{
  par_pkint(r.refSlice);
  par_pkint(r.imageSlice);
  par_pkint(r.updated);
  par_pkdouble(r.distortion);
  par_pkdouble(r.correlation);
  par_pkdouble(r.correspondence);
  par_pkstr(r.message);
}

void
UnpackResult ()
{
  r.refSlice = par_upkint();
  r.imageSlice = par_upkint();
  r.updated = par_upkint();
  r.distortion = par_upkdouble();
  r.correlation = par_upkdouble();
  r.correspondence = par_upkdouble();
  par_upkstr(r.message);
}

int
ReadHeader (FILE *f, char *tc, unsigned int *w, unsigned int *h, int *m)
{
  int c;
  int v;

  c = fgetc(f);
  while (c == '#')
    {
      while ((c = fgetc(f)) != EOF && c != '\n') ;
    }
  if (c != 'P')
    return(0);
  c = fgetc(f);
  if (c != '4' && c != '5' && c != '6')
    return(0);
  *tc = c;
  c = fgetc(f);

  while (c == ' ' || c == '\t' || c == '\n' || c == '#')
    if (c == '#')
      {
	while ((c = fgetc(f)) != EOF && c != '\n') ;
      }
    else
      c = fgetc(f);
  if (c >= '0' && c <= '9')
    {
      v = 0;
      while (c != EOF && c >= '0' && c <= '9')
	{
	  v = 10 * v + (c - '0');
	  c = fgetc(f);
	}
      *w = v;
    }
  else
    return(0);

  while (c == ' ' || c == '\t' || c == '\n' || c == '#')
    if (c == '#')
      {
	while ((c = fgetc(f)) != EOF && c != '\n') ;
      }
    else
      c = fgetc(f);
  if (c >= '0' && c <= '9')
    {
      v = 0;
      while (c != EOF && c >= '0' && c <= '9')
	{
	  v = 10 * v + (c - '0');
	  c = fgetc(f);
	}
      *h = v;
    }
  else
    return(0);

  if (*tc == '4')
    return(1);

  while (c == ' ' || c == '\t' || c == '\n' || c == '#')
    if (c == '#')
      while ((c = fgetc(f)) != EOF && c != '\n') ;
    else
      c = fgetc(f);
  if (c >= '0' && c <= '9')
    {
      v = 0;
      while (c != EOF && c >= '0' && c <= '9')
	{
	  v = 10 * v + (c - '0');
	  c = fgetc(f);
	}
      *m = v;
    }
  else
    return(0);

  return(1);
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

