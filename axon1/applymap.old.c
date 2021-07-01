#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <dirent.h>
#include <stdarg.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "csa/csa.h"
#undef __GNUC__
#include "csa/nan.h"
#define __GNUC__
#include "ptinpoly.h"

#define CHECKING	0
#define TRACKING	0
#define MASKING		1

#define N_LEVELS	16       /* max image size (w or h) is 2^N_LEVELS */
#define MAP(map,x,y)	map[((y)+1)*mpwp + (x) + 1]
#define NEW_MAP(map,x,y)	map[((y)+1)*new_mpwp + (x) + 1]
#define IMAGE(i,x,y)    (x >= 0 && x < w && y >= 0 && y < h ? i[(y)*w + (x)] : background)
#define LINE_LENGTH	255

int forward;
int noCrop;

char inputName[PATH_MAX];
int inputDigits;
int inputWidth;
int inputHeight;
int inputBits;

char outputName[PATH_MAX];

int scale;
int mapFactor;
char mapName[PATH_MAX];

float* map;                  /* map to use */
float* image;                /* image to be warped */
float* output;               /* warped output image */

int nProcessed = 0;

#if 0
/* the following was included to gain access to some internal fields
   of a csa struct for debugging */
struct csa {
    double xmin;
    double xmax;
    double ymin;
    double ymax;

    int npoints;
    point** points;
    int npointsallocated;

    int nstd;
    double** std;
    int nstdallocated;

    /*
     * squarization 
     */
    int ni;
    int nj;
    double h;
    struct square*** squares;          /* square* [j][i] */

    int npt;                    /* Number of Primary Triangles */
    struct triangle** pt;              /* Primary Triangles -- triangle* [npt] */
    int nincreased;             /* Number of sub-datasets thinned */
    int nthinned;               /* Number of sub-datasets increased */
    int norder[4];              /* Number of fittings of given interpolation
                                 * order */

    /*
     * algorithm parameters 
     */
    int npmin;                  /* minimal number of points locally involved
                                 * in spline calculation (normally = 3) */
    int npmax;                  /* maximal number of points locally involved
                                 * in spline calculation (required > 10,
                                 * recommended 20 < npmax < 60) */
    double k;                   /* relative tolerance multiple in fitting
                                 * spline coefficients: the higher this
                                 * value, the higher degree of the locally
                                 * fitted spline (recommended 80 < k < 200) */
    int nppc;                   /* average number of points per cell */
};
#endif

int GridTest( register pGridSet p_gs, double point[2] );


/* FORWARD DECLARATIONS */
int Compare (const void *x, const void *y);
int ParseRange (char *s, int *pos, int *minValue, int *maxValue);
int ParseValue (char *s, int *pos, int *value);
int ReadHeader (FILE *f, char *tc, int *w, int *h, int *m);
void Error (char *fmt, ...);

int
main (int argc, char **argv, char **envp)
{
  int i;
  size_t ii;
  int n;
  int nSkipSlices;
  int maxSkipSlices;
  int *skipSlices;
  int len;
  char skipFile[PATH_MAX];
  int minSlice, maxSlice;
  int error;
  int pos;
  DIR *dir;
  struct dirent *de;
  char fn[PATH_MAX];
  char inputDirName[PATH_MAX];
  char inputPrefix[PATH_MAX];
  int inputPrefixLen;
  int minN, maxN;
  int digitLen;
  char tc;
  int w, h;
  int m;
  int skipIndex;
  int minS, maxS;
  float *xmap;
  float *ymap;
  float *cmap;
  size_t mpwp, mphp;
  int x, y;
  int irx, iry;
  double rx, ry;
  double rrx, rry;
  float *warp;
  unsigned char *warpout;
  int mpw, mph;
  float *image;
  double rx00, rx01, rx10, rx11, ry00, ry01, ry10, ry11;
  int imageSlice;
  unsigned char *input;
  int mapHeader[3];
  double r00, r01, r10, r11;
  double rv;
  FILE *f;
  size_t nPixels;
  float background;
  unsigned short v;
  double pt[2];
  int maxMapRes;
  int gridResolution;
  int nVertices;
  double (*polygon)[2];
  GridSet gridSet;
  csa* csx;
  csa* csy;
  point csa_pt;
  point *xPoints;
  point *yPoints;
  size_t oPixels;
  float minX, minY, maxX, maxY;
  int oMinX, oMinY, oMaxX, oMaxY;
  size_t oWidth, oHeight;
  int oxmin, oymin, oxmax, oymax;
  double xv, yv;
  int ixv, iyv;
  int new_mpw, new_mph;
  size_t new_mpwp, new_mphp;
  float *new_xmap;
  float *new_ymap;
  float *new_cmap;
  int lmyx, lmyy;

  error = 0;
  forward = 1;
  noCrop = 0;
  inputName[0] = '\0';
  outputName[0] = '\0';
  scale = 1;
  mapFactor = 1;
  mapName[0] = '\0';
  skipFile[0] = '\0';
  minSlice = 0;
  maxSlice = 1000000000;
  background = 0.0;
  for (i = 1; i < argc; ++i)
    if (strcmp(argv[i], "-forward") == 0)
      forward = 1;
    else if (strcmp(argv[i], "-reverse") == 0)
      forward = 0;
    else if (strcmp(argv[i], "-nocrop") == 0)
      noCrop = 1;
    else if (strcmp(argv[i], "-input") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(inputName, argv[i]);
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
    else if (strcmp(argv[i], "-map") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(mapName, argv[i]);
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
    else if (strcmp(argv[i], "-factor") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d", &mapFactor) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-scale") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d", &scale) != 1)
	  {
	    error = 1;
	    break;
	  }
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
    else error = 1;

  if (error)
    {
      fprintf(stderr, "Usage: applymap -input file_prefix -output file_prefix\n");
      fprintf(stderr, "              -map map_file_prefix\n");
      fprintf(stderr, "              [-skip skip_file] [-slice slice_range]\n");
      fprintf(stderr, "   where ranges are expressed as: integer\n");
      fprintf(stderr, "                              or: integer-integer\n");
      fprintf(stderr, "                              or: integer-\n");
      fprintf(stderr, "                              or: -integer\n");
      exit(1);
    }

  /* check that at least minimal parameters were supplied */
  if (inputName[0] == '\0' || outputName[0] == '\0' || mapName[0] == '\0')
    {
      fprintf(stderr, "-input and -output and -map parameters must be specified.\n");
      exit(1);
    }

  /* get a sorted list of the slices to skip */
  nSkipSlices = 0;
  skipSlices = (int *) malloc(sizeof(int));
  maxSkipSlices = 1;
  if (skipFile[0] != '\0')
    {
      f = fopen(skipFile, "r");
      if (f == NULL)
	{
	  fprintf(stderr, "Could not open skip file %s\n", skipFile);
	  exit(1);
	}
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
	{
	  fprintf(stderr, "Read error while reading skip file %s\n", skipFile);
	  exit(1);
	}
      
      fclose(f);
      qsort(skipSlices, nSkipSlices, sizeof(int), Compare);
    }    

  /* check what images are present */
  if (strlen(inputName) > 4 &&
      strcmp(&inputName[strlen(inputName) - 4], ".pgm") == 0)
    {
      f = fopen(inputName, "r");
      if (f == NULL)
	Error("Could not open image file %s", inputName);
      if (!ReadHeader(f, &tc, &inputWidth, &inputHeight, &m))
	Error("Invalid image file header: %s\n", inputName);
      if (tc != '5')
	Error("Expected binary pgm file: %s\n", inputName);
      if (m != 255 && m != 65535)
	Error("Expected 8-bit or 16-bit pgm file: %s\n", inputName);
      if (m == 255)
	inputBits = 8;
      else
	inputBits = 16;
      fclose(f);
    }
  else
    {
      for (i = strlen(inputName)-1; i >= 0 && inputName[i] != '/'; --i) ;
      if (i >= 0)
	{
	  strncpy(inputDirName, inputName, i+1);
	  inputDirName[i+1] = '\0';
	  strcpy(inputPrefix, &inputName[i+1]);
	}
      else
	{
	  strcpy(inputDirName, "./");
	  strcpy(inputPrefix, inputName);
	}
      
      /* read the directory to look for files of the form nameNNNN.ppm  (or .pgm) */
      printf("Scanning input directory %s\n", inputDirName);
      dir = opendir(inputDirName);
      if (dir == NULL)
	{
	  fprintf(stderr, "Could not open directory %s for input files\n",
		  inputDirName);
	  exit(1);
	}
      inputPrefixLen = strlen(inputPrefix);
      minN = 1000000000;
      maxN = -1;
      inputDigits = -1;
      while ((de = readdir(dir)) != NULL)
	{
	  if (strncmp(de->d_name, inputPrefix, inputPrefixLen) == 0 &&
	      (len = strlen(de->d_name)) > 4 &&
	      strcmp(&(de->d_name[len-4]), ".pgm") == 0)
	    {
	      if (sscanf(&(de->d_name[inputPrefixLen]), "%d%n",
			 &n, &digitLen) < 1)
		continue;
	      if (inputDigits < 0)
		inputDigits = digitLen;
	      else if (digitLen != inputDigits)
		Error("Inconsistent number of digits in input file names: %s\n", de->d_name);
	      if (n < minN)	
	{
		  minN = n;
		  
		  /* read the header to determine image size */
		  sprintf(fn, "%s%s", inputDirName, de->d_name);
		  f = fopen(fn, "r");
		  if (f == NULL)
		    Error("Could not open image file %s", de->d_name);
		  if (!ReadHeader(f, &tc, &inputWidth, &inputHeight, &m))
		    Error("Invalid image file header: %s\n", de->d_name);
		  if (tc != '5' || (m != 255 && m != 65535))
		    Error("Expected 8-bit or 16-bit pgm file: %s\n", de->d_name);
		  if (m == 255)
		    inputBits = 8;
		  else
		    inputBits = 16;
		  fclose(f);
		}
	      if (n > maxN)
		maxN = n;
	    }
	}
      closedir(dir);
      if (maxN < 0)
	{
	  fprintf(stderr, "No input image files found.\n");
	  exit(1);
	}
      printf("minN = %d maxN = %d bits = %d\n", minN, maxN, inputBits);
    }

  if (forward && noCrop)
    {
      printf("Previewing maps: ");
      fflush(stdout);

      oMinX = 1000000000;
      oMaxX = -1000000000;
      oMinY = 1000000000;
      oMaxY = -1000000000;
      skipIndex = 0;
      for (imageSlice = minN; imageSlice <= maxN; ++imageSlice)
	{
	  if (imageSlice < minSlice || imageSlice > maxSlice)
	    continue;
	  while (skipIndex < nSkipSlices && skipSlices[skipIndex] < imageSlice)
	    ++skipIndex;
	  if (skipIndex < nSkipSlices && skipSlices[skipIndex] == imageSlice)
	    continue;

	  if (strlen(inputName) > 4 &&
	      strcmp(&inputName[strlen(inputName) - 4], ".pgm") == 0)
	    strcpy(fn, inputName);
	  else
	    sprintf(fn, "%s%0.*d.pgm", inputName, inputDigits, imageSlice);
	  f = fopen(fn, "r");
	  if (f == NULL)
	    Error("Could not open image file %s\n", fn);
	  if (!ReadHeader(f, &tc, &w, &h, &m))
	    Error("Could not read header of %s\n", fn);
	  fclose(f);

	  if (strlen(mapName) > 4 &&
	      strcmp(&mapName[strlen(mapName) - 4], ".map") == 0)
	    strcpy(fn, mapName);
	  else
	    sprintf(fn, "%s%0.*d.map", mapName, inputDigits, imageSlice);
	  f = fopen(fn, "r");
	  if (f == NULL)
	    Error("Could not open map file %s\n", fn);
	  if (fread(mapHeader, sizeof(int), 3, f) != 3)
	    Error("Could not read map header from file %s\n", fn);
	  mpw = mapHeader[0];
	  mph = mapHeader[1];
	  mpwp = mpw + 2;
	  mphp = mph + 2;
	  xmap = (float *) malloc(mpwp*mphp*sizeof(float));
	  ymap = (float *) malloc(mpwp*mphp*sizeof(float));
	  cmap = (float *) malloc(mpwp*mphp*sizeof(float));
	  if (fread(xmap, sizeof(float), mpwp*mphp, f) != mpwp*mphp ||
	      fread(ymap, sizeof(float), mpwp*mphp, f) != mpwp*mphp ||
	      fread(cmap, sizeof(float), mpwp*mphp, f) != mpwp*mphp)
	    Error("Map file %s apparently truncated\n", fn);
	  fclose(f);

#if 0
	  /* crop map */
	  new_mpw = mpw - 512;
	  new_mph = mph - 512;
	  new_mpwp = new_mpw + 2;
	  new_mphp = new_mph + 2;
	  new_xmap = (float *) malloc(new_mpwp*new_mphp*sizeof(float));
	  new_ymap = (float *) malloc(new_mpwp*new_mphp*sizeof(float));
	  new_cmap = (float *) malloc(new_mpwp*new_mphp*sizeof(float));
	  for (y = -1; y <= new_mph; ++y)
	    for (x = -1; x < new_mpw; ++x)
	      {
		NEW_MAP(new_xmap, x, y) = MAP(xmap, x, y);
		NEW_MAP(new_ymap, x, y) = MAP(ymap, x, y);
		NEW_MAP(new_cmap, x, y) = MAP(cmap, x, y);
	      }
	  free(xmap);
	  free(ymap);
	  free(cmap);
	  mpw = new_mpw;
	  mph = new_mph;
	  mpwp = new_mpwp;
	  mphp = new_mphp;
	  xmap = new_xmap;
	  ymap = new_ymap;
	  cmap = new_cmap;
#endif

	  minX = 1000000000.0;
	  maxX = -1000000000.0;
	  minY = 1000000000.0;
	  maxY = -1000000000.0;
	  lmyx = -1000;
	  lmyy = -1000;
	  for (x = -1; x <= mpw; ++x)
	    {
	      MAP(xmap, x, -1) = (2.0 * MAP(xmap, x, -1) - scale + 1) /
		(2.0 * scale);
	      MAP(ymap, x, -1) = (2.0 * MAP(ymap, x, -1) - scale + 1) /
		(2.0 * scale);
	      if (MAP(xmap, x, -1) < minX)
		minX = MAP(xmap, x, -1);
	      if (MAP(xmap, x, -1) > maxX)
		maxX = MAP(xmap, x, -1);
	      if (MAP(ymap, x, -1) < minY)
		minY = MAP(ymap, x, -1);
	      if (MAP(ymap, x, -1) > maxY)
		{
		  maxY = MAP(ymap, x, -1);
		  lmyx = x;
		  lmyy = -1;
		}

	      MAP(xmap, x, mph) = (2.0 * MAP(xmap, x, mph) - scale + 1) /
		(2.0 * scale);
	      MAP(ymap, x, mph) = (2.0 * MAP(ymap, x, mph) - scale + 1) /
		(2.0 * scale);
	      if (MAP(xmap, x, mph) < minX)
		minX = MAP(xmap, x, mph);
	      if (MAP(xmap, x, mph) > maxX)
		maxX = MAP(xmap, x, mph);
	      if (MAP(ymap, x, mph) < minY)
		minY = MAP(ymap, x, mph);
	      if (MAP(ymap, x, mph) > maxY)
		{
		  maxY = MAP(ymap, x, mph);
		  lmyx = x;
		  lmyy = mph;
		}
	    }

	  for (y = 0; y < mph; ++y)
	    {
	      MAP(xmap, -1, y) = (2.0 * MAP(xmap, -1, y) - scale + 1) /
		(2.0 * scale);
	      MAP(ymap, -1, y) = (2.0 * MAP(ymap, -1, y) - scale + 1) /
		(2.0 * scale);
	      if (MAP(xmap, -1, y) < minX)
		minX = MAP(xmap, -1, y);
	      if (MAP(xmap, -1, y) > maxX)
		maxX = MAP(xmap, -1, y);
	      if (MAP(ymap, -1, y) < minY)
		minY = MAP(ymap, -1, y);
	      if (MAP(ymap, -1, y) > maxY)
		{
		  maxY = MAP(ymap, -1, y);
		  lmyx = -1;
		  lmyy = y;
		}

	      MAP(xmap, mpw, y) = (2.0 * MAP(xmap, mpw, y) - scale + 1) /
		(2.0 * scale);
	      MAP(ymap, mpw, y) = (2.0 * MAP(ymap, mpw, y) - scale + 1) /
		(2.0 * scale);
	      if (MAP(xmap, mpw, y) < minX)
		minX = MAP(xmap, mpw, y);
	      if (MAP(xmap, mpw, y) > maxX)
		maxX = MAP(xmap, mpw, y);
	      if (MAP(ymap, mpw, y) < minY)
		minY = MAP(ymap, mpw, y);
	      if (MAP(ymap, mpw, y) > maxY)
		{
		  maxY = MAP(ymap, mpw, y);
		  lmyx = mpw;
		  lmyy = y;
		}
	    }

	  free(xmap);
	  free(ymap);
	  free(cmap);

	  oxmin = (int) floor(mapFactor * minX + (mapFactor - 1.0) / 2.0);
	  oymin = (int) floor(mapFactor * minY + (mapFactor - 1.0) / 2.0);
	  oxmax = (int) ceil(mapFactor * maxX + (mapFactor - 1.0) / 2.0);
	  oymax = (int) ceil(mapFactor * maxY + (mapFactor - 1.0) / 2.0);
	  printf("mapFactor = %d oxmin = %d oymin = %d oxmax = %d oymax = %d\n",
		 mapFactor, oxmin, oymin, oxmax, oymax);
	  printf("lmyx = %d lmyy = %d\n", lmyx, lmyy);
	  if (oxmin < oMinX)
	    oMinX = oxmin;
	  if (oymin < oMinY)
	    oMinY = oymin;
	  if (oxmax > oMaxX)
	    oMaxX = oxmax;
	  if (oymax > oMaxY)
	    oMaxY = oymax;

	  if ((nProcessed % 50) == 0 && nProcessed != 0)
	    printf(" %d\n                ", nProcessed);
	  printf(".");
	  fflush(stdout);
	  ++nProcessed;
	}
      printf("\nAll slices previewed.\n");

      oWidth = oMaxX - oMinX + 1;
      oHeight = oMaxY - oMinY + 1;
      printf("output width = %lu (%d to %d) output height = %lu (%d to %d)\n\n",
	     oWidth, oMinX, oMaxX,
	     oHeight, oMinY, oMaxY);
    }
  else
    {
      oMinX = 0;
      oMinY = 0;
      oMaxX = inputWidth-1;
      oMaxY = inputHeight-1;
      oWidth = inputWidth;
      oHeight = inputHeight;
    }

  /* for all slices */
  nProcessed = 0;
  printf("Processing slices: ");
  fflush(stdout);

  skipIndex = 0;
  for (imageSlice = minN; imageSlice <= maxN; ++imageSlice)
    {
      if (imageSlice < minSlice || imageSlice > maxSlice)
	continue;
      while (skipIndex < nSkipSlices && skipSlices[skipIndex] < imageSlice)
	++skipIndex;
      if (skipIndex < nSkipSlices && skipSlices[skipIndex] == imageSlice)
	continue;

      /* read in the image to warp */
      if (strlen(inputName) > 4 &&
	  strcmp(&inputName[strlen(inputName) - 4], ".pgm") == 0)
	strcpy(fn, inputName);
      else
	sprintf(fn, "%s%0.*d.pgm", inputName, inputDigits, imageSlice);
      f = fopen(fn, "r");
      if (f == NULL)
	Error("Could not open image file %s\n", fn);
      if (!ReadHeader(f, &tc, &w, &h, &m))
	Error("Could not read header of %s\n", fn);
      nPixels = ((size_t) w) * h;
      input = (unsigned char *) malloc(nPixels*inputBits/8);
      if (fread(input, 1, nPixels*inputBits/8, f) != nPixels*inputBits/8)
	Error("Image file %s apparently truncated\n", fn);
      fclose(f);
      image = (float *) malloc(w * h * sizeof(float));
      if (inputBits == 8)
	for (ii = 0; ii < nPixels; ++ii)
	  image[ii] = input[ii];
      else
	for (ii = 0; ii < nPixels; ++ii)
	  image[ii] = (input[2*ii] << 8) | input[2*ii+1];
      free(input);

      /* get the map for the image */
      if (strlen(mapName) > 4 &&
	  strcmp(&mapName[strlen(mapName) - 4], ".map") == 0)
	strcpy(fn, mapName);
      else
	sprintf(fn, "%s%0.*d.map", mapName, inputDigits, imageSlice);
      f = fopen(fn, "r");
      if (f == NULL)
	Error("Could not open map file %s\n", fn);
      if (fread(mapHeader, sizeof(int), 3, f) != 3)
	Error("Could not read map header from file %s\n", fn);
      mpw = mapHeader[0];
      mph = mapHeader[1];
      if (mapFactor == 0)
	{
	  mapFactor = w / mpw;
	  if (mpw * mapFactor != w || mph * mapFactor != h)
	    Error("Image resolution not a multiple of map resolution.\n");
	}
      mpwp = mpw + 2;
      mphp = mph + 2;
      xmap = (float *) malloc(mpwp*mphp*sizeof(float));
      ymap = (float *) malloc(mpwp*mphp*sizeof(float));
      cmap = (float *) malloc(mpwp*mphp*sizeof(float));
      if (fread(xmap, sizeof(float), mpwp*mphp, f) != mpwp*mphp ||
	  fread(ymap, sizeof(float), mpwp*mphp, f) != mpwp*mphp ||
	  fread(cmap, sizeof(float), mpwp*mphp, f) != mpwp*mphp)
	Error("Map file %s apparently truncated\n", fn);
      fclose(f);

#if 0
	  /* crop map */
	  new_mpw = mpw - 512;
	  new_mph = mph - 512;
	  new_mpwp = new_mpw + 2;
	  new_mphp = new_mph + 2;
	  new_xmap = (float *) malloc(new_mpwp*new_mphp*sizeof(float));
	  new_ymap = (float *) malloc(new_mpwp*new_mphp*sizeof(float));
	  new_cmap = (float *) malloc(new_mpwp*new_mphp*sizeof(float));
	  for (y = -1; y <= new_mph; ++y)
	    for (x = -1; x < new_mpw; ++x)
	      {
		NEW_MAP(new_xmap, x, y) = MAP(xmap, x, y);
		NEW_MAP(new_ymap, x, y) = MAP(ymap, x, y);
		NEW_MAP(new_cmap, x, y) = MAP(cmap, x, y);
	      }
	  free(xmap);
	  free(ymap);
	  free(cmap);
	  mpw = new_mpw;
	  mph = new_mph;
	  mpwp = new_mpwp;
	  mphp = new_mphp;
	  xmap = new_xmap;
	  ymap = new_ymap;
	  cmap = new_cmap;
#endif

      for (y = -1; y <= mph; ++y)
	for (x = -1; x <= mpw; ++x)
	  {
	    MAP(xmap, x, y) = (2.0 * MAP(xmap, x, y) - scale + 1) /
	      (2.0 * scale);
	    MAP(ymap, x, y) = (2.0 * MAP(ymap, x, y) - scale + 1) /
	      (2.0 * scale);
	  }
      output = (float *) malloc(oWidth * oHeight * sizeof(float));

      if (forward)
	{
	  maxMapRes = mpw;
	  if (mph > maxMapRes)
	    maxMapRes = mph;
	  gridResolution = maxMapRes / 16;
	  if (gridResolution < 4)
	    gridResolution = 4;
	  nVertices = mpw + mph + mpw + mph + 4;
	  polygon = (double (*)[2]) malloc(nVertices * sizeof(double[2]));
	  i = 0;
	  for (x = 0; x < mpw; ++x)
	    {
	      polygon[i][0] = 0.5 * (MAP(xmap, x, 0) + MAP(xmap, x, -1));
	      polygon[i][1] = 0.5 * (MAP(ymap, x, 0) + MAP(ymap, x, -1));
	      ++i;
	    }
	  polygon[i][0] = 0.5 * (MAP(xmap, mpw-1, 0) + MAP(xmap, mpw, -1));
	  polygon[i][1] = 0.5 * (MAP(ymap, mpw-1, 0) + MAP(ymap, mpw, -1));
	  ++i;
	  for (y = 0; y < mph; ++y)
	    {
	      polygon[i][0] = 0.5 * (MAP(xmap, mpw-1, y) + MAP(xmap, mpw, y));
	      polygon[i][1] = 0.5 * (MAP(ymap, mpw-1, y) + MAP(ymap, mpw, y));
	      ++i;
	    }
	  polygon[i][0] = 0.5 * (MAP(xmap, mpw-1, mph-1) + MAP(xmap, mpw, mph));
	  polygon[i][1] = 0.5 * (MAP(ymap, mpw-1, mph-1) + MAP(ymap, mpw, mph));
	  ++i;
	  for (x = mpw-1; x >= 0; --x)
	    {
	      polygon[i][0] = 0.5 * (MAP(xmap, x, mph-1) + MAP(xmap, x, mph));
	      polygon[i][1] = 0.5 * (MAP(ymap, x, mph-1) + MAP(ymap, x, mph));
	      ++i;
	    }
	  polygon[i][0] = 0.5 * (MAP(xmap, 0, mph-1) + MAP(xmap, -1, mph));
	  polygon[i][1] = 0.5 * (MAP(ymap, 0, mph-1) + MAP(ymap, -1, mph));
	  ++i;
	  for (y = mph-1; y >= 0; --y)
	    {
	      polygon[i][0] = 0.5 * (MAP(xmap, 0, y) + MAP(xmap, -1, y));
	      polygon[i][1] = 0.5 * (MAP(ymap, 0, y) + MAP(ymap, -1, y));
	      ++i;
	    }
	  polygon[i][0] = 0.5 * (MAP(xmap, 0, 0) + MAP(xmap, -1, -1));
	  polygon[i][1] = 0.5 * (MAP(ymap, 0, 0) + MAP(ymap, -1, -1));
	  ++i;
	  /*	  printf("Before GridSetup\n");
		  fflush(stdout); */
	  GridSetup(polygon, nVertices, gridResolution, &gridSet);
	  /*	  printf("After GridSetup\n");
		  fflush(stdout); */

	  /*  use Bivariate Cubic Spline approximations to compute the
	      original coordinates of the reconstructed dataset */
	  xPoints = (point *) malloc(mpwp * mphp * sizeof(point));
	  yPoints = (point *) malloc(mpwp * mphp * sizeof(point));
	  i = 0;
	  for (y = -1; y <= mph; ++y)
	    for (x = -1; x <= mpw; ++x)
	      {
		xPoints[i].x = MAP(xmap, x, y);
		xPoints[i].y = MAP(ymap, x, y);
		xPoints[i].z = x;
		yPoints[i].x = MAP(xmap, x, y);
		yPoints[i].y = MAP(ymap, x, y);
		yPoints[i].z = y;
		++i;
	      }

	  /*	  printf("Before csa_create #1\n");
		  fflush(stdout); */
	  csx = csa_create();
	  csa_setnpmin(csx, 20);
	  csa_setnpmax(csx, 80);
	  csa_setk(csx, 80);
	  csa_setnppc(csx, 30);
	  csa_addpoints(csx, mpwp * mphp, xPoints);
	  csa_calculatespline(csx);
	  /*	  printf("After csa_create #1\n");
		  fflush(stdout); */


	  /*	  printf("Before csa_create #2\n");
		  fflush(stdout); */
	  csy = csa_create();
	  csa_setnpmin(csy, 20);
	  csa_setnpmax(csy, 80);
	  csa_setk(csy, 80);
	  csa_setnppc(csy, 30);
	  csa_addpoints(csy, mpwp * mphp, yPoints);
	  csa_calculatespline(csy);
	  /*	  printf("After csa_create #2\n");
		  fflush(stdout); */

	  for (y = oMinY; y <= oMaxY; ++y)
	    {
	    for (x = oMinX; x <= oMaxX; ++x)
	      {
#if 0
		pt[0] = (2.0 * x - mapFactor + 1) / (2.0 * mapFactor);
		pt[1] = (2.0 * y - mapFactor + 1) / (2.0 * mapFactor);
#endif
		pt[0] = (2.0 * x - scale + 1) / (2.0 * scale);
		pt[1] = (2.0 * y - scale + 1) / (2.0 * scale);
		
		if (!GridTest(&gridSet, pt))
		  {
		    output[(y - oMinY) * oWidth + (x - oMinX)] = 0;
		    continue;
		  }

		csa_pt.x = pt[0];
		csa_pt.y = pt[1];
		csa_approximatepoint(csx, &csa_pt);
		rx = (csa_pt.z != NaN) ?
		  (mapFactor * csa_pt.z + (mapFactor - 1.0) / 2.0) : -1.0;
		csa_approximatepoint(csy, &csa_pt);
		ry = (csa_pt.z != NaN) ?
		  (mapFactor * csa_pt.z + (mapFactor - 1.0) / 2.0) : -1.0;

		irx = (int) floor(rx);
		iry = (int) floor(ry);
		rrx = rx - irx;
		rry = ry - iry;
		if (irx < 0 || irx >= w || irx == w-1 && rrx > 0.0 ||
		    iry < 0 || iry >= h || iry == h-1 && rry > 0.0)
		  {
		    output[(y - oMinY) * oWidth + (x - oMinX)] = 0;
		    continue;
		  }
		r00 = IMAGE(image, irx, iry);
		r01 = IMAGE(image, irx, iry + 1);
		r10 = IMAGE(image, irx + 1, iry);
		r11 = IMAGE(image, irx + 1, iry + 1);
		rv = r00 * (rrx - 1.0) * (rry - 1.0)
		  - r10 * rrx * (rry - 1.0) 
		  - r01 * (rrx - 1.0) * rry
		  + r11 * rrx * rry;
		output[(y - oMinY) * oWidth + (x - oMinX)] = rv;
	      }
	    /*	    printf(".");
		    fflush(stdout); */
	    }
	  GridCleanup(&gridSet);
	  csa_destroy(csx);
	  csa_destroy(csy);
	  free(xPoints);
	  free(yPoints);
	  free(polygon);
	}
      else
	{
	  /* compute the warped image */
	  for (y = oMinY; y <= oMaxY; ++y)
	    for (x = oMinX; x <= oMaxX; ++x)
	      {
		/* use bilinear interpolation to find value */
		xv = (2.0 * x - mapFactor + 1) / (2.0 * mapFactor);
		yv = (2.0 * y - mapFactor + 1) / (2.0 * mapFactor);
		ixv = (int) floor(xv);
		iyv = (int) floor(yv);
		rrx = xv - ixv;
		rry = yv - iyv;
		if (ixv < -1 || ixv > mpw || ixv == mpw && rrx > 0.0 ||
		    iyv < -1 || iyv > mph || iyv == mph && rry > 0.0)
		  {
		    output[(y - oMinY) * oWidth + (x - oMinX)] = 0;
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
		rx = mapFactor * rx + (mapFactor - 1.0) / 2.0;
		ry = mapFactor * ry + (mapFactor - 1.0) / 2.0;
		irx = (int) floor(rx);
		iry = (int) floor(ry);
		rrx = rx - irx;
		rry = ry - iry;
		if (irx < 0 || irx >= w || irx == w-1 && rrx > 0.0 ||
		    iry < 0 || iry >= h || iry == h-1 && rry > 0.0)
		  {
		    output[(y - oMinY) * oWidth + (x - oMinX)] = 0;
		    continue;
		  }

		r00 = IMAGE(image, irx, iry);
		r01 = IMAGE(image, irx, iry + 1);
		r10 = IMAGE(image, irx + 1, iry);
		r11 = IMAGE(image, irx + 1, iry + 1);
		rv = r00 * (rrx - 1.0) * (rry - 1.0)
		  - r10 * rrx * (rry - 1.0) 
		  - r01 * (rrx - 1.0) * rry
		  + r11 * rrx * rry;
		output[(y - oMinY) * oWidth + (x - oMinX)] = rv;

	      }
	}

      /* write out the warped image */
      oPixels = oWidth * oHeight;
      warpout = (unsigned char*) malloc(oPixels * inputBits / 8);
      if (strlen(outputName) > 4 &&
	  strcmp(&outputName[strlen(outputName) - 4], ".pgm") == 0)
	strcpy(fn, outputName);
      else
	sprintf(fn, "%s%0.*d.pgm", outputName, inputDigits, imageSlice);
      f = fopen(fn, "w");
      if (f == NULL)
	Error("Could not open image file %s for writing.\n", fn);
      fprintf(f, "P5\n%lu %lu\n%d\n", oWidth, oHeight, (1 << inputBits) - 1);
      if (inputBits == 8)
	for (ii = 0; ii < oPixels; ++ii)
	  if (output[ii] < 0.0)
	    warpout[ii] = 0;
	  else if (output[ii] >= 255.0)
	    warpout[ii] = 255;
	  else
	    warpout[ii] = (int) floor(output[ii] + 0.5);
      else
	for (ii = 0; ii < oPixels; ++ii)
	  {
	    if (output[ii] < 0.0)
	      v = 0;
	    else if (output[ii] >= 65535.0)
	      v = 65535;
	    else
	      v = (int) floor(output[ii]+0.5);
	    warpout[2*ii] = v & 0xff;
	    warpout[2*ii+1] = (v >> 8) & 0xff;
	  }
      
      if (fwrite(warpout, inputBits / 8, oPixels, f) != oPixels)
	Error("Could not write to image file %s\n", fn);
      fclose(f);

      free(image);
      free(xmap);
      free(ymap);
      free(cmap);
      free(output);
      free(warpout);

      if ((nProcessed % 50) == 0 && nProcessed != 0)
	printf(" %d\n                ", nProcessed);
      printf(".");
      fflush(stdout);
      ++nProcessed;
    }

  printf("\nAll slices completed.\n");
  exit(0);
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

int
ReadHeader (FILE *f, char *tc, int *w, int *h, int *m)
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

void Error (char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
  fflush(stderr);
  abort();
}
