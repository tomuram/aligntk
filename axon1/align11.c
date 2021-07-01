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

#define DEBUG	0

#define CONSTRAINED	(1.0e+30)
#define MAPIN(map,x,y)	map[((y)+1)*((size_t) mpwip) + (x) + 1]

typedef struct Point
{
  float x;
  float y;
  float fx;   /* if fx > 0.5*CONSTRAINED, Point is not allowed to
		 move in the x direction */
  float fy;   /* if fy > 0.5*CONSTRAINED, Point is not allowed to
		 move in the y direction */
} Point;

typedef struct IntraTableEntry
{
  int index1;
  int index2;
  float nomD;
  float k;
} IntraTableEntry;

typedef struct InterTableEntry
{
  Point *pt1;
  Point *pt2;
  float rrx, rry;    /* ratio of goal offset to distances to adjacent points
			in target map; should be in range -0.5 to 0.5 */
  float dx, dy;      /* pixel offset for distance computations */
  float energyFactor;
  int s1;
  int s2;
} InterTableEntry;

int p;      /* rank of this process */
int np;     /* number of processes in this run */

int params[12];
float fParams[5];
int *imageWidth, *imageHeight;
char inputName[PATH_MAX];
int inputDigits;
char alphaName[PATH_MAX];
char inputMapName[PATH_MAX];
char outputMapName[PATH_MAX];
char outputGridName[PATH_MAX];
int outputDigits;
char constraintName[PATH_MAX];
int constraintDigits;
char sliceMagnificationsName[PATH_MAX];
int nSkipSlices;
int *skipSlices;
int nMaps;
int *mapNameLen = 0;
char **mapName = 0;
int *mapDigits = 0;
int *mapForward = 0;
int *mapMin = 0;
int *mapMax = 0;
int fixedSlice = -1;
int minIter = 0;
char termName[PATH_MAX];
char triggerName[PATH_MAX];
float *sliceMagnifications;
int showConstraints = 0;

Point **pts = 0;
int *nIntraTableEntries;
IntraTableEntry **intraTables = 0;
int nInterTableEntries = 0;
InterTableEntry *interTable = 0;
FILE *logFile = 0;
float kAbsolute = 0.0;
float kIntra = 1.0;
float kInter = 0.03;
float tighten = 1.0;
float fixedDamping = 0.0;
int spacing = 4;
int mapFactor = 0;
int outputVolumes = 0;
int *actualSlice = 0;
int myFirstSlice, myLastSlice;

unsigned char *field;

unsigned char colors[11][3] =
  {{153, 0, 0},     // dark red
   {255, 255, 0},   // yellow
   {0, 255, 0},     // green
   {0, 255, 255},   // turquoise
   {255, 153, 0},   // orange

   {0, 51, 255},    // bright blue
   {204, 0, 255},   // bright purple
   {255, 204, 102}, // gold
   {0, 204, 255},   // light blue
   {255, 102, 153}, // pink

   {255, 0, 0}      // red
  };

int Compare (const void *x, const void *y);
int ComparePoints (const void *x, const void *y);
int ParseRange (char *s, int *pos, int *minValue, int *maxValue);
int ParseValue (char *s, int *pos, int *value);
int ReadHeader (FILE *f, char *tc, int *w, int *h, int *m);
void Error (char *fmt, ...);
void Log (char *fmt, ...);
void Output (int iter);
void hsv_to_rgb (unsigned char* r, unsigned char *g, unsigned char *b,
		 float h, float s, float v);
void DrawLine(unsigned char *img, float px0, float py0, float px1, float py1,
	      float nomD, int imageWidth, int imageHeight);
void UpdateDeltas (int iter);

int
main (int argc, char **argv)
{
  char fn[PATH_MAX];
  FILE *f;
  int x, y, z;
  int pixel;
  DIR *dir;
  struct dirent *de;
  int minZ, maxZ;
  int digitLen;
  int myMinZ, myMaxZ;
  int i, j, k;
  int error;
  int maxSkipSlices;
  char skipFile[PATH_MAX];
  int pos;
  int minS, maxS;
  int minSlice, maxSlice;
  char inputDirName[PATH_MAX];
  char inputPrefix[PATH_MAX];
  int inputPrefixLen;
  char mapDirName[PATH_MAX];
  char mapPrefix[PATH_MAX];;
  int mapPrefixLen;
  char constraintDirName[PATH_MAX];
  char constraintPrefix[PATH_MAX];
  int constraintPrefixLen;
  char tc;
  int iw, ih;
  int im;
  int n;
  int m;
  int w, h;
  float ***mapX;
  float ***mapY;
  float ***mapC;
  int **mnx;
  int **mny;
  int iz;
  int skipIndex;
  int nextZ;
  float rx, ry, rc;
  int irx, iry;
  float rrx, rry;
  float r00, r01, r10, r11;
  float xv, yv;
  float rv;
  Point *pt;
  int dx, dy;
  int wp, hp;
  double tw;
  double sum;
  float *fx;
  float *fy;
  int *slice;
  int sz, ez;
  float *mx, *my, *mc;
  int mnx2, mny2;
  double totalEnergy;
  double prevTotalEnergy;
  double deltaEnergy;
  double energy;
  double energyCorrection;
  double totalEnergyCorrection;
  int iter;
  int sepB, sepF, sep;
  float deltaX, deltaY;
  int xp, yp;
  float nomD;
  float d;
  float force;
  float s;
  float kdx, kdy;
  float dampingFactor;
  MPI_Status status;
  int mapHeader[3];
  int nTimes;
  int nDecrease;
  int nIncrease;
  float forceOverD;
  float dfx, dfy;
  unsigned char red, green, blue;
  float angle;
  float mag;
  int nPts;
  int *nFloats;
  IntraTableEntry *it;
  InterTableEntry *ite;
  Point *p00, *p10, *p01, *p11;
  float kIntraThisSlice;
  int terminate;
  int writeOutput;
  struct stat sb;
  int ixv, iyv;
  int ind;
  float rx00, rx01, rx10, rx11;
  float ry00, ry01, ry10, ry11;
  float rc00, rc01, rc10, rc11;
  time_t lastOutput;
  int cx, cy;
  float consX, consY;
  IntraTableEntry *intraTable;
  float maxF;
  float globalMaxF;
  float scale;
  char hostName[256];
  double basis;
  float maxStepX, maxStepY;
  float pk;
  int mpwi, mphi;
  int mpwip, mphip;
  int inputMapFactor;
  float x0, y0;
  int ix0, iy0;
  float min_rx, max_rx, min_ry, max_ry;
  int mpf;
  int totalPts;
  int nSlices;
  int nx, ny, nz;
  int nx2, ny2;
  int nxj, nyj;
  
  /* initialize MPI */
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    Error("Could not do MPI_Init\n");
  if (MPI_Comm_size(MPI_COMM_WORLD, &np) != MPI_SUCCESS)
    Error("Could not do MPI_Comm_size\n");
  if (MPI_Comm_rank(MPI_COMM_WORLD, &p) != MPI_SUCCESS)
    Error("Could not do MPI_Comm_rank\n");
  printf("I am node %d of %d\n", p, np);
  fflush(stdout);
  
  sprintf(fn, "logs/%0.2d.log", p);
  logFile = fopen(fn, "w");
  gethostname(hostName, 255);
  Log("Process %d is running on %s\n", p, hostName);
  if (p == 0)
    {
      error = 0;
      alphaName[0] = '\0';
      inputMapName[0] = '\0';
      outputMapName[0] = '\0';
      outputGridName[0] = '\0';
      triggerName[0] = '\0';
      termName[0] = '\0';
      constraintName[0] = '\0';
      sliceMagnificationsName[0] = '\0';
      nMaps = 0;
      mapName = NULL;
      skipFile[0] = '\0';
      minSlice = 0;
      maxSlice = 1000000000;
      fixedSlice = -1;
      minIter = 0;

      for (i = 0; i < argc; ++i)
	Log("ARGV[%d] = %s\n", i, argv[i]);
      for (i = 1; i < argc; ++i)
	if (strcmp(argv[i], "-input") == 0)
	  {
	    if (++i == argc)
	      {
		error = 1;
		break;
	      }
	    strcpy(inputName, argv[i]);
	  }
	else if (strcmp(argv[i], "-map") == 0 ||
	    strcmp(argv[i], "-partialmap") == 0)
	  {
	    if (++i == argc)
	      {
		error = 1;
		break;
	      }
	    ++nMaps;
	    mapNameLen = (int *) realloc(mapNameLen, nMaps*sizeof(int));
	    mapName = (char **) realloc(mapName, nMaps*sizeof(char*));
	    mapDigits = (int *) realloc(mapDigits, nMaps*sizeof(int));
	    mapForward = (int *) realloc(mapForward, nMaps*sizeof(int));
	    mapMin = (int *) realloc(mapMin, nMaps*sizeof(int));
	    mapMax = (int *) realloc(mapMax, nMaps*sizeof(int));
	    mapNameLen[nMaps-1] = strlen(argv[i]);
	    mapName[nMaps-1] = (char *) malloc(mapNameLen[nMaps-1]+1);
	    strcpy(mapName[nMaps-1], argv[i]);
	  }
	else if (strcmp(argv[i], "-inputmap") == 0)
	  {
	    if (++i == argc)
	      {
		error = 1;
		break;
	      }
	    strcpy(inputMapName, argv[i]);
	  }
	else if (strcmp(argv[i], "-outputmap") == 0)
	  {
	    if (++i == argc)
	      {
		error = 1;
		break;
	      }
	    strcpy(outputMapName, argv[i]);
	  }	    
	else if (strcmp(argv[i], "-outputgrid") == 0)
	  {
	    if (++i == argc)
	      {
		error = 1;
		break;
	      }
	    strcpy(outputGridName, argv[i]);
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
	else if (strcmp(argv[i], "-fixed") == 0)
	  {
	    if (++i == argc ||
		sscanf(argv[i], "%d", &fixedSlice) != 1)
	      {
		error = 1;
		break;
	      }
	  }
	else if (strcmp(argv[i], "-iter") == 0)
	  {
	    if (++i == argc ||
		sscanf(argv[i], "%d", &minIter) != 1)
	      {
		error = 1;
		break;
	      }
	  }
	else if (strcmp(argv[i], "-term") == 0)
	  {
	    if (++i == argc)
	      {
		error = 1;
		break;
	      }
	    strcpy(termName, argv[i]);
	  }
	else if (strcmp(argv[i], "-trigger") == 0)
	  {
	    if (++i == argc)
	      {
		error = 1;
		break;
	      }
	    strcpy(triggerName, argv[i]);
	  }
	else if (strcmp(argv[i], "-damping") == 0)
	  {
	    if (++i == argc ||
		sscanf(argv[i], "%f", &fixedDamping) != 1)
	      {
		error = 1;
		break;
	      }
	  }
	else if (strcmp(argv[i], "-intra") == 0)
	  {
	    if (++i == argc ||
		sscanf(argv[i], "%f", &kIntra) != 1)
	      {
		error = 1;
		break;
	      }
	  }
	else if (strcmp(argv[i], "-inter") == 0)
	  {
	    if (++i == argc ||
		sscanf(argv[i], "%f", &kInter) != 1)
	      {
		error = 1;
		break;
	      }
	  }
	else if (strcmp(argv[i], "-absolute") == 0)
	  {
	    if (++i == argc ||
		sscanf(argv[i], "%f", &kAbsolute) != 1)
	      {
		error = 1;
		break;
	      }
	  }
	else if (strcmp(argv[i], "-tighten") == 0)
	  {
	    if (++i == argc ||
		sscanf(argv[i], "%f", &tighten) != 1)
	      {
		error = 1;
		break;
	      }
	  }
	else if (strcmp(argv[i], "-spacing") == 0)
	  {
	    if (++i == argc ||
		sscanf(argv[i], "%d", &spacing) != 1)
	      {
		error = 1;
		break;
	      }
	  }
	else if (strcmp(argv[i], "-volume") == 0)
	  outputVolumes = 1;
	else if (strcmp(argv[i], "-constraints") == 0)
	  {
	    if (++i == argc)
	      {
		error = 1;
		break;
	      }
	    strcpy(constraintName, argv[i]);
	  }
	else if (strcmp(argv[i], "-slice_magnifications") == 0)
	  {
	    if (++i == argc)
	      {
		error = 1;
		break;
	      }
	    strcpy(sliceMagnificationsName, argv[i]);
	  }
	else if (strcmp(argv[i], "-show_constraints") == 0)
	  showConstraints = 1;
	else error = 1;

      if (error)
	{
	  fprintf(stderr, "Usage: align11\n");
	  fprintf(stderr, "              [-skip skip_file] [-slice slice_range]\n");
	  fprintf(stderr, "              [-map map_file_prefix] [-fixed slice_number]\n");
	  fprintf(stderr, "              [-inputmap map_prefix]\n");
	  fprintf(stderr, "              [-outputmap map_prefix]\n");
	  fprintf(stderr, "              [-constraints constraints_prefix]\n");
	  fprintf(stderr, "              [-slice_magnification file]\n");
	  fprintf(stderr, "              [-intra weight] [-inter weight] [-absolute weight]\n");
	  fprintf(stderr, "              [-tighten tightening_factor]\n");
	  fprintf(stderr, "              [-damping damping_factor]\n");
	  fprintf(stderr, "   where ranges are expressed as: integer\n");
	  fprintf(stderr, "                              or: integer-integer\n");
	  fprintf(stderr, "                              or: integer-\n");
	  fprintf(stderr, "                              or: -integer\n");
	  exit(1);
	}
      
      /* check that at least minimal parameters were supplied */
      if (inputName[0] == '\0')
	Error("-input must be specified\n");
      if (outputMapName[0] == '\0' && outputGridName[0] == '\0')
	Error("At least one of -outputmap or -outputgrid must be specified.\n");
      if (nMaps == 0)
	Error("At least one -map argument must be specified.\n");

      /* set default trigger and termination files */
      if (triggerName[0] == '\0')
	sprintf(triggerName, "%strigger", outputMapName);
      if (termName[0] == '\0')
	sprintf(termName, "%sterm", outputMapName);

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
		  skipSlices = (int*) realloc(skipSlices,
					      maxSkipSlices * sizeof(int));
		}
	      skipSlices[nSkipSlices-1] = n;
	    }
	  if (!feof(f))
	    Error("Read error while reading skip file %s\n", skipFile);
	  fclose(f);
	  qsort(skipSlices, nSkipSlices, sizeof(int), Compare);
	}    

      /* check what images are present */
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

      /* read the directory to look for files of the form nameNNNN.pgm */
      dir = opendir(inputDirName);
      if (dir == NULL)
	Error("Could not open directory %s for image files\n");
      inputPrefixLen = strlen(inputPrefix);
      minZ = 1000000000;
      maxZ = -1;
      inputDigits = -1;
      while ((de = readdir(dir)) != NULL)
	{
	  if (strncmp(de->d_name, inputPrefix, inputPrefixLen) != 0 ||
	      sscanf(&(de->d_name[inputPrefixLen]), "%d%n",
		     &n, &digitLen) < 1 ||
	      strcmp(&(de->d_name[inputPrefixLen+digitLen]), ".pgm") != 0)
	    continue;
	  if (inputDigits < 0)
	    inputDigits = digitLen;
	  else if (digitLen != inputDigits)
	    Error("Inconsistent number of digits in image file names: %s\n",
		  de->d_name);
	  if (n < minSlice || n > maxSlice)
	    continue;
	  if (n < minZ)
	    minZ = n;
	  if (n > maxZ)
	    maxZ = n;
	}
      closedir(dir);
      nSlices = maxZ - minZ + 1;
      if (nSlices <= 0)
	Error("No input image files found.\n");
      if (nSlices == 1)
	Error("Not enough slices selected to do alignment!\n");

      /* read the image sizes */
      slice = (int *) malloc(nSlices * sizeof(int));
      for (i = 0; i < nSlices; ++i)
	slice[i] = -1;
      actualSlice = (int *) malloc(nSlices * sizeof(int));
      imageWidth = (int *) malloc(nSlices * sizeof(int));
      imageHeight = (int *) malloc(nSlices * sizeof(int));
      skipIndex = 0;
      nz = 0;
      for (z = minZ; z <= maxZ; ++z)
	{
	  /* see if we should skip this slice */
	  while (skipIndex < nSkipSlices && skipSlices[skipIndex] < z)
	    ++skipIndex;
	  if (skipIndex < nSkipSlices && skipSlices[skipIndex] == z)
	    continue;
	  
	  /* read the header to determine image size */
	  sprintf(fn, "%s%0.*d.pgm", inputName, inputDigits, z);
	  f = fopen(fn, "r");
	  if (f == NULL)
	    Error("Could not open image file %s", fn);
	  if (!ReadHeader(f, &tc, &w, &h, &m))
	    Error("Could not read header of image file %s\n", fn);
	  if (m != 255)
	    Error("Image file %s does not have byte-size image components.", fn);
	  fclose(f);

	  slice[z - minZ] = nz;
	  actualSlice[nz] = z;
	  imageWidth[nz] = w;
	  imageHeight[nz] = h;
	  ++nz;
	}
      
      /* check that the map files are present */
      mapFactor = 0;
      for (m = 0; m < nMaps; ++m)
	{
	  for (i = strlen(mapName[m])-1;
	       i >= 0 && mapName[m][i] != '/';
	       --i) ;
	  if (i >= 0)
	    {
	      strncpy(mapDirName, mapName[m], i+1);
	      mapDirName[i+1] = '\0';
	      strcpy(mapPrefix, &mapName[m][i+1]);
	    }
	  else
	    {
	      strcpy(mapDirName, ".");
	      strcpy(mapPrefix, mapName[m]);
	    }

	  /* read the directory to look for files of the form nameNNNN.map */
	  dir = opendir(mapDirName);
	  if (dir == NULL)
	    Error("Could not open directory %s for map files\n", mapDirName);
	  mapPrefixLen = strlen(mapPrefix);
	  mapMin[m] = 1000000000;
	  mapMax[m] = -1;
	  mapDigits[m] = -1;
	  mapForward[m] = -1;
	  while ((de = readdir(dir)) != NULL)
	    {
	      if (strncmp(de->d_name, mapPrefix, mapPrefixLen) != 0 ||
		  sscanf(&(de->d_name[mapPrefixLen]), "%d%n",
			 &n, &digitLen) < 1 ||
		  strcmp(&(de->d_name[mapPrefixLen+digitLen]), ".map") != 0)
		continue;
	      if (mapDigits[m] < 0)
		mapDigits[m] = digitLen;
	      else if (digitLen != mapDigits[m])
		Error("Inconsistent number of digits in map file names: %s\n", de->d_name);
	      if (n < minZ || n > maxZ || slice[n - minZ] < 0)
		continue;
	      if (n < mapMin[m] || n > mapMax[m])
		{
		  if (n < mapMin[m])
		    mapMin[m] = n;
		  if (n > mapMax[m])
		    mapMax[m] = n;
		  
		  /* read the header to determine the map reference */
		  sprintf(fn, "%s%s", mapDirName, de->d_name);
		  f = fopen(fn, "r");
		  if (f == NULL)
		    Error("Could not open map file %s", fn);
		  if (fread(mapHeader, sizeof(int), 3, f) != 3)
		    Error("Could not read header of map file %s\n", fn);
		  fclose(f);
		  mapForward[m] = (mapHeader[2] > n);
		  mpf = 1 << ((int) floor(log(((float) imageWidth[slice[n - minZ]]) / mapHeader[0]) / log(2.0) + 0.5));
		  if (mapFactor == 0)
		    mapFactor = mpf;
		  else if (mpf != mapFactor)
		    Error("Maps must all be at the same scale.\n");
		  if (mapHeader[2] >= minSlice && mapHeader[2] <= maxSlice)
		    {
		      if (mapHeader[2] < mapMin[m])
			mapMin[m] = mapHeader[2];
		      if (mapHeader[2] > mapMax[m])
			mapMax[m] = mapHeader[2];
		    }
		}
	    }
	}

      /* check what constraints are present */
      for (i = strlen(constraintName)-1; i >= 0 &&
	     constraintName[i] != '/'; --i) ;
      if (i >= 0)
	{
	  strncpy(constraintDirName, constraintName, i+1);
	  constraintDirName[i+1] = '\0';
	  strcpy(constraintPrefix, &constraintName[i+1]);
	}
      else
	{
	  strcpy(constraintDirName, "./");
	  strcpy(constraintPrefix, constraintName);
	}
      /* read the directory to look for files of the form nameNNNN.con */
      dir = opendir(constraintDirName);
      if (dir == NULL)
	Error("Could not open directory %s for constraint files\n");
      constraintPrefixLen = strlen(constraintPrefix);
      constraintDigits = -1;
      while ((de = readdir(dir)) != NULL)
	{
	  if (strncmp(de->d_name, constraintPrefix,
		      constraintPrefixLen) != 0 ||
	      sscanf(&(de->d_name[constraintPrefixLen]), "%d%n",
		     &n, &digitLen) < 1 ||
	      strcmp(&(de->d_name[constraintPrefixLen+digitLen]),
		     ".con") != 0)
	    continue;
	  if (constraintDigits < 0)
	    constraintDigits = digitLen;
	  else if (digitLen != constraintDigits)
	    Error("Inconsistent number of digits in constraint file names: %s\n",
		  de->d_name);
	}
      closedir(dir);

      sliceMagnifications = (float *) malloc(nz * sizeof(float));
      for (i = 0; i < nz; ++i)
	sliceMagnifications[i] = 1.0;
      if (sliceMagnificationsName[0] != '\0')
	{
	  f = fopen(sliceMagnificationsName, "r");
	  if (f == NULL)
	    Error("Could not open slice magnifications file %s for reading\n",
		  sliceMagnificationsName);
	  while (fscanf(f, "%d %f", &z, &mag) == 2)
	    if (z >= minZ && z <= maxZ && slice[z - minZ] >= 0)
	      sliceMagnifications[slice[z - minZ]] = mag;
	  fclose(f);
	}


      params[0] = minZ;
      params[1] = maxZ;
      params[2] = nz;
      params[3] = nMaps;
      params[4] = fixedSlice;
      params[5] = spacing;
      params[6] = mapFactor;
      params[7] = outputVolumes;
      params[8] = constraintDigits;
      params[9] = showConstraints;
      params[10] = minIter;
      params[11] = inputDigits;

      fParams[0] = kAbsolute;
      fParams[1] = kIntra;
      fParams[2] = kInter;
      fParams[3] = tighten;
      fParams[4] = fixedDamping;
    }

  /* broadcast the info */
  if (MPI_Bcast(inputMapName, PATH_MAX, MPI_CHAR, 0, MPI_COMM_WORLD) !=
      MPI_SUCCESS)
    Error("Broadcast of input map name failed.");
  if (MPI_Bcast(outputMapName, PATH_MAX, MPI_CHAR, 0, MPI_COMM_WORLD) !=
      MPI_SUCCESS)
    Error("Broadcast of output map name failed.");
  if (MPI_Bcast(outputGridName, PATH_MAX, MPI_CHAR, 0, MPI_COMM_WORLD) !=
      MPI_SUCCESS)
    Error("Broadcast of output map name failed.");
  if (MPI_Bcast(constraintName, PATH_MAX, MPI_CHAR, 0, MPI_COMM_WORLD) !=
      MPI_SUCCESS)
    Error("Broadcast of constraint name failed.");
  if (MPI_Bcast(params, 12, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("Broadcast of image parameters failed.");
  if (MPI_Bcast(fParams, 5, MPI_FLOAT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("Broadcast of image float parameters failed.");
  minZ = params[0];
  maxZ = params[1];
  nz = params[2];
  nMaps = params[3];
  fixedSlice = params[4];
  spacing = params[5];
  mapFactor = params[6];
  outputVolumes = params[7];
  constraintDigits = params[8];
  showConstraints = params[9];
  minIter = params[10];
  inputDigits = params[11];
  kAbsolute = fParams[0];
  kIntra = fParams[1];
  kInter = fParams[2];
  tighten = fParams[3];
  fixedDamping = fParams[4];
  Log("minZ = %d maxZ = %d  nM = %d at %d\n",
      minZ, maxZ, nMaps, p);

  nSlices = maxZ - minZ + 1;
  if (p != 0)
    {
      slice = (int *) malloc(nSlices * sizeof(int));
      actualSlice = (int *) malloc(nz * sizeof(int));
      imageWidth = (int *) malloc(nz * sizeof(int));
      imageHeight = (int *) malloc(nz * sizeof(int));
      sliceMagnifications = (float *) malloc(nz * sizeof(float));
    }
  if (MPI_Bcast(slice, nSlices, MPI_INT,
		0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(actualSlice, nz, MPI_INT,
		0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(imageWidth, nz, MPI_INT,
		0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(imageHeight, nz, MPI_INT,
		0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(sliceMagnifications, nz, MPI_FLOAT,
		0, MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("Broadcast of slice info failed.");

  myFirstSlice = (p * nz) / np;
  myLastSlice = ((p+1) * nz / np) - 1;
  Log("On node %d first = %d last = %d (nz = %d)\n",
      p, myFirstSlice, myLastSlice, nz);

  pts = (Point **) malloc(nz * sizeof(Point *));
  memset(pts, 0, nz * sizeof(Point *));
  nFloats = malloc(nz * sizeof(int));
  memset(nFloats, 0, nz * sizeof(int));
  sz = (myFirstSlice > 0) ? (myFirstSlice-1) : 0;
  ez = (myLastSlice < nz-1) ? (myLastSlice+1) : (nz-1);
  totalPts = 0;
  for (i = sz; i <= ez; ++i)
    {
      nx = imageWidth[i] / spacing;
      ny = imageHeight[i] / spacing;
      nx2 = nx + 2;
      ny2 = ny + 2;
      nPts = nx * ny;
      totalPts += nPts;
      nFloats[i] = (sizeof(Point) / sizeof(float)) * nPts;
      pts[i] = (Point *) malloc(nPts * sizeof(Point));
      pt = pts[i];
      if (inputMapName[0] != '\0')
	{
	  sprintf(fn, "%s%0.*d.map", inputMapName, inputDigits, actualSlice[i]);
	  f = fopen(fn, "r");
	  if (f == NULL)
	    Error("Could not open input map file %s\n", fn);
	  if (fread(mapHeader, sizeof(int), 3, f) != 3)
	    Error("Could not read map header from file %s\n", fn);
	  mpwi = mapHeader[0];
	  mphi = mapHeader[1];
	  inputMapFactor = 1 << ((int) floor(log(((float) imageWidth[i]) / mpwi) / log(2.0) + 0.5));
	  mpwip = mpwi + 2;
	  mphip = mphi + 2;
	  mx = (float *) malloc(mpwip * mphip * sizeof(float));
	  my = (float *) malloc(mpwip * mphip * sizeof(float));
	  if (fread(mx, sizeof(float), mpwip*mphip, f) != mpwip*mphip ||
	      fread(my, sizeof(float), mpwip*mphip, f) != mpwip*mphip)
	    Error("Could not read map from file %s  (mpw = %d mph = %d)\n",
		  fn, mpwi, mphi);

	  min_rx = 1000000.0;
	  min_ry = 1000000.0;
	  max_rx = -1000000.0;
	  max_ry = -1000000.0;
	  for (y = 0; y < ny; ++y)
	    for (x = 0; x < nx; ++x)
	      {
		rx = spacing * x + (spacing - 1.0) / 2.0;
		ry = spacing * y + (spacing - 1.0) / 2.0;
		x0 = (2.0 * rx - inputMapFactor + 1) / (2.0 * inputMapFactor);
		y0 = (2.0 * ry - inputMapFactor + 1) / (2.0 * inputMapFactor);
		ix0 = (int) floor(x0);
		iy0 = (int) floor(y0);
		rrx = x0 - ix0;
		rry = y0 - iy0;
		if (ix0 >= mpwi)
		  {
		    rrx += ix0 - mpwi + 1;
		    ix0 = mpwi-1;
		  }
		if (iy0 >= mphi)
		  {
		    rry += iy0 - mphi + 1;
		    iy0 = mphi-1;
		  }
		rx00 = MAPIN(mx, ix0, iy0);
		ry00 = MAPIN(my, ix0, iy0);
		rx01 = MAPIN(mx, ix0, iy0 + 1);
		ry01 = MAPIN(my, ix0, iy0 + 1);
		rx10 = MAPIN(mx, ix0 + 1, iy0);
		ry10 = MAPIN(my, ix0 + 1, iy0);
		rx11 = MAPIN(mx, ix0 + 1, iy0 + 1);
		ry11 = MAPIN(my, ix0 + 1, iy0 + 1);
		rx = rx00 * (rrx - 1.0) * (rry - 1.0)
		  - rx10 * rrx * (rry - 1.0) 
		  - rx01 * (rrx - 1.0) * rry
		  + rx11 * rrx * rry;
		ry = ry00 * (rrx - 1.0) * (rry - 1.0)
		  - ry10 * rrx * (rry - 1.0) 
		  - ry01 * (rrx - 1.0) * rry
		  + ry11 * rrx * rry;
		  
		if (rx < min_rx)
		  min_rx = rx;
		if (rx > max_rx)
		  max_rx = rx;
		if (ry < min_ry)
		  min_ry = ry;
		if (ry > max_ry)
		  max_ry = ry;
		pt->x = rx;
		pt->y = ry;
		pt->fx = 0.0;
		pt->fy = 0.0;
		++pt;
	      }
	  Log("For slice %d (%d) mpw = %d mph = %d mf = %d (%f < rx < %f) (%f < ry < %f)\n",
	      i, actualSlice[i], mpwi, mphi, inputMapFactor, min_rx, max_rx, min_ry, max_ry);
	  free(mx);
	  free(my);
	}
      else
	for (y = 0; y < ny; ++y)
	  for (x = 0; x < nx; ++x)
	    {
	      pt->x = spacing * x + (spacing - 1.0) / 2.0;
	      pt->y = spacing * y + (spacing - 1.0) / 2.0;
	      pt->fx = 0.0;
	      pt->fy = 0.0;
	      ++pt;
	    }
      if (constraintName[0] != '\0')
	{
	  sprintf(fn, "%s%0.*d.con", constraintName, constraintDigits,
		  actualSlice[i]);
	  f = fopen(fn, "r");
	  if (f != NULL)
	    {
	      pt = pts[i];
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

  if (p != 0)
    {
      mapNameLen = (int*) malloc(nMaps*sizeof(int));
      mapName = (char**) malloc(nMaps*sizeof(char*));
      mapDigits = (int*) malloc(nMaps*sizeof(int));
      mapForward = (int*) malloc(nMaps*sizeof(int));
      mapMin = (int*) malloc(nMaps*sizeof(int));
      mapMax = (int*) malloc(nMaps*sizeof(int));
    }
  if (MPI_Bcast(mapNameLen, nMaps, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("Broadcast of mapNameLen failed.");
  for (m = 0; m < nMaps; ++m)
    {
      if (p != 0)
	mapName[m] = (char *) malloc(mapNameLen[m]+1);
      if (MPI_Bcast(mapName[m], mapNameLen[m]+1, MPI_BYTE, 0, MPI_COMM_WORLD) !=
	  MPI_SUCCESS)
	Error("Broadcast of mapName failed.");
    }
  if (MPI_Bcast(mapDigits, nMaps, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("Broadcast of mapDigits failed.");
  if (MPI_Bcast(mapForward, nMaps, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("Broadcast of mapForward failed.");
  if (MPI_Bcast(mapMin, nMaps, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("Broadcast of mapMin failed.");
  if (MPI_Bcast(mapMax, nMaps, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("Broadcast of mapMax failed.");

  mapX = (float ***) malloc(nMaps * sizeof(float**));
  mapY = (float ***) malloc(nMaps * sizeof(float**));
  mapC = (float ***) malloc(nMaps * sizeof(float**));
  mnx = (int **) malloc(nMaps * sizeof(int*));
  mny = (int **) malloc(nMaps * sizeof(int*));
  memset(mapX, 0, nMaps*sizeof(float**));
  memset(mapY, 0, nMaps*sizeof(float**));
  memset(mapC, 0, nMaps*sizeof(float**));
  memset(mnx, 0, nMaps*sizeof(int*));
  memset(mny, 0, nMaps*sizeof(int*));

  Log("GOING THROUGH MAPS\n");
  /* go through the maps one by one */
  for (m = 0; m < nMaps; ++m)
    {
      Log("LOADING MAP %d\n", m);

      Log("ALLOCATING\n");
      mapX[m] = (float **) malloc(nz * sizeof(float*));
      mapY[m] = (float **) malloc(nz * sizeof(float*));
      mapC[m] = (float **) malloc(nz * sizeof(float*));
      mnx[m] = (int *) malloc(nz * sizeof(int));
      mny[m] = (int *) malloc(nz * sizeof(int));
      memset(mapX[m], 0, nz * sizeof(float*));
      memset(mapY[m], 0, nz * sizeof(float*));
      memset(mapC[m], 0, nz * sizeof(float*));
      memset(mnx[m], 0, nz * sizeof(int));
      memset(mny[m], 0, nz * sizeof(int));

      Log("READING\n");
      /* read in the map data */
      for (iz = mapMin[m]; iz <= mapMax[m]; ++iz)
	{
	  if (mapForward[m])
	    z = iz;
	  else
	    z = mapMax[m] - iz + mapMin[m];
	  i = slice[z - minZ];
	  if (i < sz || i > ez)
	    continue;

	  /* check if the map is present */
	  sprintf(fn, "%s%0.*d.map", mapName[m], mapDigits[m], z);
	  Log("Attempting to read map %s into slot %d\n", fn, i);
	  f = fopen(fn, "r");
	  if (f == NULL)
	    {
	      Log("Could not open map file %s -- ignoring\n", fn);
	      continue;
	    }

	  /* check if the map is between locally-resident slices */
	  if (fread(mapHeader, sizeof(int), 3, f) != 3)
	    Error("Could not read map header from file %s\n", fn);
	  if (mapHeader[2] < actualSlice[sz] ||
	      mapHeader[2] > actualSlice[ez])
	    {
	      Log("Map reference is not in selected slices -- ignoring map\n");
	      fclose(f);
	      continue;
	    }
	  if (mapForward[m])
	    {
	      if (i >= nz-1 || mapHeader[2] != actualSlice[i+1])
		Error("Map files do not correspond to requested slices. (%d %d %d %d %d\n",
		      mapForward[m], i, nz, mapHeader[2], actualSlice[i+1]);
	    }
	  else
	    {
	      if (i == 0 || mapHeader[2] != actualSlice[i-1])
		Error("Map files do not correspond to requested slices. (%d %d %d %d %d\n",
		      mapForward[m], i, nz, mapHeader[2], actualSlice[i-1]);
	    }

	  /* allocate storage */
	  mnx[m][i] = mapHeader[0];
	  mny[m][i] = mapHeader[1];
	  mnx2 = mnx[m][i] + 2;
	  mny2 = mny[m][i] + 2;
	  mx = (float *) malloc(mnx2 * mny2 * sizeof(float));
	  my = (float *) malloc(mnx2 * mny2 * sizeof(float));
	  mc = (float *) malloc(mnx2 * mny2 * sizeof(float));
	  mapX[m][i] = mx;
	  mapY[m][i] = my;
	  mapC[m][i] = mc;

	  if (fread(mx, sizeof(float), mnx2*mny2, f) != mnx2*mny2 ||
	      fread(my, sizeof(float), mnx2*mny2, f) != mnx2*mny2 ||
	      fread(mc, sizeof(float), mnx2*mny2, f) != mnx2*mny2)
	    Error("Could not read map from file %s\n", fn);
	  fclose(f);
	  Log("Successfully read map %s\n", fn);
	}
    }
  
  if (p == 0)
    Log("Building tables\n");
  intraTables = (IntraTableEntry **) malloc(nz * sizeof(IntraTableEntry*));
  nIntraTableEntries = (int *) malloc(nz * sizeof(int));
  for (i = sz; i <= ez; ++i)
    {
      nx = imageWidth[i] / spacing;
      ny = imageHeight[i] / spacing;
      nPts = nx * ny;
      mag = sliceMagnifications[i];
      Log("Slice: %d mag = %f\n", i, mag);
      for (j = sz; j < i; ++j)
	if (sliceMagnifications[j] == mag &&
	    nx == imageWidth[j] / spacing &&
	    ny == imageHeight[j] / spacing)
	  {
	    Log("Using table %d for table %d\n", j, i);
	    nIntraTableEntries[i] = nIntraTableEntries[j];
	    intraTables[i] = intraTables[j];
	    break;
	  }
      if (j >= i)
	{
	  Log("Constructing table %d\n", i);
	  intraTable = (IntraTableEntry*)
	    malloc(nPts * 4 * sizeof(IntraTableEntry));
	  k = 0;
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
		      it = &intraTable[k++];
		      it->index1 = y * nx + x;
		      it->index2 = yp * nx + xp;
		      it->nomD = mag * spacing * ((dx * dy != 0) ? M_SQRT2 : 1.0);
		      it->k = kIntra;
		    }
		}
	  nIntraTableEntries[i] = k;
	  Log("Constructed  (k = %d)\n", k);
	  intraTables[i] = intraTable;
	}
    }
  maxStepX = 0.1 * spacing;
  maxStepY = 0.1 * spacing;

  interTable = (InterTableEntry*) malloc(nMaps * totalPts *
					 sizeof(InterTableEntry));
  Log("sz = %d ez = %d myFirstSlice = %d myLastSlice = %d\n",
      sz, ez, myFirstSlice, myLastSlice);
  k = 0;
  for (i = sz; i <= ez; ++i)
    {
      nx = imageWidth[i] / spacing;
      ny = imageHeight[i] / spacing;
      sepB = (i != 0) ? (actualSlice[i] - actualSlice[i-1]) : 1;
      sepF = (i != nz-1) ? (actualSlice[i+1] - actualSlice[i]) : 1;
      for (m = 0; m < nMaps; ++m)
	{
	  mx = mapX[m][i];
	  my = mapY[m][i];
	  mc = mapC[m][i];
	  if (mx == 0 || my == 0 || mc == 0)
	    continue;
	  mnx2 = mnx[m][i] + 2;
	  mny2 = mnx[m][i] + 2;
	  if (mapForward[m])
	    {
	      if (i > myLastSlice)
		continue;
	      j = i + 1;
	      if (j >= nz)
		continue;
	      sep = sepF;
	      Log("ITC F i = %d j = %d sep = %d\n", i, j, sep);
	    }
	  else
	    {
	      if (i < myFirstSlice)
		continue;
	      j = i - 1;
	      if (j < 0)
		continue;
	      sep = sepB;
	      Log("ITC B i = %d j = %d sep = %d\n", i, j, sep);
	    }
	  nxj = imageWidth[j] / spacing;
	  nyj = imageHeight[j] / spacing;

	  for (y = 0; y < ny; ++y)
	    for (x = 0; x < nx; ++x)
	      {
		xv = spacing * x + (spacing - 1.0) / 2.0;
		yv = spacing * y + (spacing - 1.0) / 2.0;
		xv = (2.0 * xv - mapFactor + 1) / (2.0 * mapFactor);
		yv = (2.0 * yv - mapFactor + 1) / (2.0 * mapFactor);
		ixv = (int) floor(xv);
		iyv = (int) floor(yv);
		if (ixv < -1 || iyv < -1)
		  continue;
		rrx = xv - ixv;
		rry = yv - iyv;
		ind = (iyv + 1) * mnx2 + ixv + 1;
		rx00 = mx[ind];
		ry00 = my[ind];
		rc00 = mc[ind];
		if (iyv < mny[m][i])
		  {
		    ind = (iyv + 2) * mnx2 + ixv + 1;
		    rx01 = mx[ind];
		    ry01 = my[ind];
		    rc01 = mc[ind];
		  }
		else
		  continue;
		if (ixv < mnx[m][i])
		  {
		    ind = (iyv + 1) * mnx2 + ixv + 2;
		    rx10 = mx[ind];
		    ry10 = my[ind];
		    rc10 = mc[ind];
		  }
		else
		  continue;
		if (ixv < mnx[m][i] && iyv < mny[m][i])
		  {
		    ind = (iyv + 2) * mnx2 + ixv + 2;
		    rx11 = mx[ind];
		    ry11 = my[ind];
		    rc11 = mc[ind];
		  }
		else
		  continue;
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
		s = rc / sep;
		if (s == 0.0)
		  continue;
		rx = mapFactor * rx + (mapFactor - 1.0) / 2.0;
		ry = mapFactor * ry + (mapFactor - 1.0) / 2.0;
		rx = (2.0 * rx - spacing + 1) / (2.0 * spacing);
		ry = (2.0 * ry - spacing + 1) / (2.0 * spacing);
		irx = (int) floor(rx + 0.5);
		iry = (int) floor(ry + 0.5);
		if (irx < 0 || irx >= nxj ||
		    iry < 0 || iry >= nyj)
		  continue;
	      	ite = &interTable[k++];
		ite->pt1 = &(pts[i][y * nx + x]);
		ite->pt2 = &(pts[j][iry * nxj + irx]);
		ite->rrx = rx - irx;
		ite->rry = ry - iry;
		ite->dx = 0.0;
		ite->dy = 0.0;
		ite->energyFactor = (i >= myFirstSlice && i <= myLastSlice) ? 1 : 0;
		ite->s1 = i;
		ite->s2 = j;
	      }
	  free(mx);
	  free(my);
	  free(mc);
	  mapX[m][i] = 0;
	  mapY[m][i] = 0;
	  mapC[m][i] = 0;
	}
    }
  nInterTableEntries = k;
  if (p == 0)
    Log("nInterTableEntries = %d\n", nInterTableEntries);

  if (p == 0)
    Log("Starting iterations\n");
  totalEnergy = 0.0;
  deltaEnergy = 1.0e+30;
  lastOutput = 0;
  nTimes = 0;
  nDecrease = 0;
  nIncrease = 0;
  if (fixedDamping != 0.0)
    dampingFactor = fixedDamping;
  else
    dampingFactor = 0.5;
  for (iter = 0; nTimes < 10; ++iter)
    {
      prevTotalEnergy = totalEnergy;
      energy = 0.0;
      energyCorrection = 0.0;

      /* update all deltas every 500 iterations */
      if (iter % 500 == 0)
	UpdateDeltas(iter);

      /* update all forces */
      for (i = myFirstSlice; i <= myLastSlice; ++i)
	{
	  nx = imageWidth[i] / spacing;
	  ny = imageHeight[i] / spacing;
	  
	  if (tighten == 1.0)
	    kIntraThisSlice = kIntra;
	  else
	    kIntraThisSlice = kIntra * ((i - nz / 2.0) * (i - nz / 2.0) *
					4.0 * (tighten - 1.0) / nz / nz + 1.0);

	  if (kAbsolute > 0.0)
	    {
	      /* compute absolute location forces */
	      pt = pts[i];
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
	      pt = pts[i];
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

	  /* add in intra-slice forces */
	  pt = pts[i];
	  basis = energy;
	  for (k = 0; k < nIntraTableEntries[i]; ++k)
	    {
	      it = &(intraTables[i][k]);
	      deltaX = (pt + it->index2)->x - (pt + it->index1)->x;
	      deltaY = (pt + it->index2)->y - (pt + it->index1)->y;
	      d = sqrt(deltaX * deltaX + deltaY * deltaY);
	      force = it->k * (d - it->nomD);
	      energy += force * (d - it->nomD);
#if DEBUG
	      if (iter == 1000000)
	      Log("intra %d: d %f nomD %f delta %f\n", k, d, it->nomD, d - it->nomD);
#endif
	      if (d != 0.0)
		{
		  forceOverD = force / d;
		  dfx = forceOverD * deltaX;
		  dfy = forceOverD * deltaY;
#if DEBUG
	      if (iter == 1000000)
		  Log("applying intraforce %f %f to %d-%d\n", dfx, dfy, i, it->index1);
	      if (iter == 1000000)
		  Log("applying intraforce %f %f to %d-%d\n", -dfx, -dfy, i, it->index2);
#endif
		  (pt + it->index1)->fx += dfx;
		  (pt + it->index1)->fy += dfy;
		  (pt + it->index2)->fx -= dfx;
		  (pt + it->index2)->fy -= dfy;
		}
	    }
#if DEBUG
	      if (iter == 1000000)
	  Log("Intraslice energy for slice %d = %f\n", i, energy - basis);
#endif
	}

      /* add in inter-slice forces (one spring method) */
      basis = energy;
      for (k = 0; k < nInterTableEntries; ++k)
	{
	  ite = &interTable[k];
	  xv = ite->pt2->x + ite->dx;
	  yv = ite->pt2->y + ite->dy;
	  deltaX = xv - ite->pt1->x;
	  deltaY = yv - ite->pt1->y;
#if DEBUG
	      if (iter == 1000000)
	  Log("inter %d: xv %f (%f + %f) x %f dx = %f  yv %f (%f + %f) y %f dy = %f\n",
	      k, xv, ite->pt2->x, ite->dx, ite->pt1->x, deltaX,
	      yv, ite->pt2->y, ite->dy, ite->pt1->y, deltaY);
#endif
	  kdx = kInter * deltaX;
	  kdy = kInter * deltaY;
#if DEBUG
	      if (iter == 1000000)
	  Log("applying interforce %f %f to %d-%d\n", kdx, kdy, ite->s1, ite->pt1 - pts[ite->s1]);
#endif
	  ite->pt1->fx += kdx;
	  ite->pt1->fy += kdy;
	  energy += ite->energyFactor * kInter * (deltaX * deltaX +
						  deltaY * deltaY);
#if DEBUG
	      if (iter == 1000000)
	  Log("applying interforce %f %f to %d-%d\n", -kdx, -kdy, ite->s2, ite->pt2 - pts[ite->s2]);
#endif
	  ite->pt2->fx -= kdx;
	  ite->pt2->fy -= kdy;
	}
#if DEBUG
	      if (iter == 1000000)
      Log("Interslice energy = %f\n", energy - basis);
#endif

      /* find maximum force */
      maxF = 0.0;
      for (i = myFirstSlice; i <= myLastSlice; ++i)
	{
	  if (fixedSlice >= 0 && actualSlice[i] == fixedSlice)
	    continue;
	  nx = imageWidth[i] / spacing;
	  ny = imageHeight[i] / spacing;
	  nPts = nx * ny;
	  pt = pts[i];
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
#if 1
      if (globalMaxF > 0.5)
	scale = dampingFactor * 0.5 / globalMaxF;
      else
	scale = dampingFactor;
#endif
#if 0
	scale = dampingFactor;
#endif
      for (i = myFirstSlice; i <= myLastSlice; ++i)
	{
	  //	  if (p == 0)
	  //	    Log("update of slice %d (actual %d)  fixedSlice = %d pt->fx = %f pt->fy = %f\n",
	  //		i, actualSlice[i], fixedSlice, pts[i]->fx, pts[i]->fy);
	      
	  if (fixedSlice >= 0 && actualSlice[i] == fixedSlice)
	    continue;
	  nx = imageWidth[i] / spacing;
	  ny = imageHeight[i] / spacing;
	  nPts = nx * ny;
	  pt = pts[i];
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
	      if (iter == 1000000)
	      Log("Net force on %d-%d is %f %f pos-adj %f %f\n",
		  i, k,
		  pt->fx, pt->fy, deltaX, deltaY);
#endif
	      ++pt;
	    }
	}

      /* transfer the boundary slices to adjacent nodes */
      /* 0<-->1, 2<-->3, 4<-->5, ... */
      if ((p & 1) == 0)
	{
	  /* send, then receive */
	  if (p != np-1)
	    {
	      if (MPI_Send(pts[myLastSlice], nFloats[myLastSlice], MPI_FLOAT,
			   p+1, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
		Error("Could not send positions forward (even-to-odd)\n");
	      if (MPI_Recv(pts[myLastSlice+1], nFloats[myLastSlice+1], MPI_FLOAT,
			   p+1, 0, MPI_COMM_WORLD, &status) != MPI_SUCCESS)
		Error("Could not receive positions backward (even-to-odd)\n");
	    }
	  
	}
      else
	{
	  /* receive, then send */
	  if (MPI_Recv(pts[myFirstSlice-1], nFloats[myFirstSlice-1], MPI_FLOAT,
		       p-1, 0, MPI_COMM_WORLD, &status) != MPI_SUCCESS)
	    Error("Could not receive positions forward (even-to-odd)\n");
	  if (MPI_Send(pts[myFirstSlice], nFloats[myFirstSlice], MPI_FLOAT,
		       p-1, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
	    Error("Could not send positions backward (even-to-odd)\n");
	}
      /* 1<-->2, 3<-->4, 5<-->6, ... */
      if ((p & 1) == 0)
	{
	  /* receive, then send */
	  if (p != 0)
	    {
	      if (MPI_Recv(pts[myFirstSlice-1], nFloats[myFirstSlice-1], MPI_FLOAT,
			   p-1, 0, MPI_COMM_WORLD, &status) != MPI_SUCCESS)
		Error("Could not receive positions backward (odd-to-even)\n");
	      if (MPI_Send(pts[myFirstSlice], nFloats[myFirstSlice], MPI_FLOAT,
			   p-1, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
		Error("Could not send positions forward (odd-to-even)\n");
	    }
	}
      else
	{
	  /* send, then receive */
	  if (p != np-1)
	    {
	      if (MPI_Send(pts[myLastSlice], nFloats[myLastSlice], MPI_FLOAT,
			   p+1, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
		Error("Could not send positions backward (odd-to-even)\n");
	      if (MPI_Recv(pts[myLastSlice+1], nFloats[myLastSlice+1], MPI_FLOAT,
			   p+1, 0, MPI_COMM_WORLD, &status) != MPI_SUCCESS)
		Error("Could not receive positions forward (odd-to-even)\n");
	    }
	}

#if DEBUG
	      if (iter == 1000000)
      Log("Local energy is %f\n", energy);
#endif
      if (MPI_Reduce(&energy, &totalEnergy, 1, MPI_DOUBLE, MPI_SUM,
		     0, MPI_COMM_WORLD) != MPI_SUCCESS)
	Error("Could not sum up energies.\n");
      if (MPI_Bcast(&totalEnergy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD) !=
	  MPI_SUCCESS)
	Error("Could not broadcast total energy.\n");
      if (MPI_Reduce(&energyCorrection, &totalEnergyCorrection, 1, MPI_DOUBLE, MPI_SUM,
		     0, MPI_COMM_WORLD) != MPI_SUCCESS)
	Error("Could not sum up energy corrections.\n");
      if (MPI_Bcast(&totalEnergyCorrection, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD) !=
	  MPI_SUCCESS)
	Error("Could not broadcast total energy.\n");
      if (p == 0 && (iter % 1000) == 999)
	/* if (p == 0) */
	Log("After %d iterations, total energy is %f  (df = %f mgf = %f tec = %f)\n",
	    iter+1, totalEnergy, dampingFactor, globalMaxF, totalEnergyCorrection);
      deltaEnergy = totalEnergy - totalEnergyCorrection - prevTotalEnergy;
      /* 32 is chosen as a time to test whether to output or terminate because,
         given that we increase the damping by 1% on each iteration, and decrease
         it by a factor of 2 when instability sets in, the number of iterations
         between instabilities is about 70, and we want to output/terminate
         when the state is not near an instability */
      if (nDecrease == 32)
	{
	  if (p == 0)
	    {
	      writeOutput = 0;
	      if (triggerName[0] != '\0' &&
		  stat(triggerName, &sb) == 0 &&
		  sb.st_mtime > lastOutput)
		{
		  Log("Update of file %s forced write of output images.\n",
		      triggerName);
		  writeOutput = 1;
		  lastOutput = sb.st_mtime;
		}
	    }
	  if (MPI_Bcast(&writeOutput, 1, MPI_INT, 0, MPI_COMM_WORLD) !=
	      MPI_SUCCESS)
	    Error("Could not broadcast writeOutput flag.\n");
	  if (writeOutput)
	    Output(iter+1);
	  if (p == 0)
	    {
	      terminate = 0;
	      if (termName[0] != '\0' &&
		  (f = fopen(termName, "r")) != NULL)
		{
		  fclose(f);
		  Log("Presence of file %s forced termination.\n", termName);
		  terminate = 1;
		}
	    }
	  if (MPI_Bcast(&terminate, 1, MPI_INT, 0, MPI_COMM_WORLD) !=
	      MPI_SUCCESS)
	    Error("Could not broadcast termination flag.\n");
	  if (terminate)
	    break;
	}
#if 0
      if (iter > minIter && fabs(deltaEnergy) < 0.000001 * totalEnergy)
	++nTimes;
      else
	nTimes = 0;
#else
      if (iter > minIter)
	++nTimes;
      else
	nTimes = 0;
#endif
#if 1
      if (deltaEnergy > 0.0)
	{
	  nDecrease = 0;
	  ++nIncrease;
	  if (p == 0)
	    {
#if 0
	      Log("Increase in energy: %d deltaEnergy %f totalEnergy %f\n",
		     nIncrease, deltaEnergy, totalEnergy);
#endif
	    }
	  if (nIncrease > 100 ||
	      nIncrease > 10 && deltaEnergy > totalEnergy ||
	      deltaEnergy > 1000.0 * totalEnergy)
	    Error("Alignment is diverging instead of converging.\n");
	      
	  if ((iter < 100 || nIncrease < 2) && fixedDamping == 0.0)
	    {
	      dampingFactor *= 0.5;
#if 0
	      if (p == 0)
		Log("Reducing dampingFactor to %f\n", dampingFactor);
#endif
	      if (dampingFactor < 0.000000001)
		Error("dampingFactor became too small\n");
	    }
	}
      else
	{
	  ++nDecrease;
	  nIncrease = 0;
#if 1
	  if (fixedDamping == 0.0 && dampingFactor < 0.5)
#endif
#if 0
	  if (fixedDamping == 0.0)
#endif
	    dampingFactor *= 1.01;
	}
#endif
    }
  if (p == 0)
    Log("Finished alignment.\n");

  Output(-1);

  Log("FINALIZING\n");
  MPI_Finalize();
  Log("DONE!\n");
  fclose(logFile);
  return(0);
}

void
Output (int iter)
{
  int x, y;
  int i, j, k;
  char fn[PATH_MAX];
  FILE *f;
  Point *pt;
  float rx, ry;
  int irx, iry;
  float rrx, rry;
  float r00, r01, r10, r11;
  float rv;
  int background;
  int nx, ny;
  int nx2, ny2;
  float *mx, *my, *mc;
  char tc;
  int iw, ih;
  int im;
  int mapHeader[3];
  int ix, iy;
  int dx, dy;
  int xp, yp;
  unsigned char *grid;

  Log("OUTPUTTING SLICES\n");
  /* for each slice in the full resolution dataset */
  outputDigits = inputDigits;
  if (p == 0)
    {
      /* create any necessary directories */
      if (outputMapName[0] != '\0')
	{
	  if (iter < 0)
	    sprintf(fn, "%sfinal", outputMapName);
	  else
	    sprintf(fn, "%s%0.6d", outputMapName, iter);
	  if (mkdir(fn, 0755) != 0)
	    Error("Could not create directory %s\n", fn);
	}
      if (outputGridName[0] != '\0')
	{
	  if (iter < 0)
	    sprintf(fn, "%sfinal", outputGridName);
	  else
	    sprintf(fn, "%s%0.6d", outputGridName, iter);
	  if (mkdir(fn, 0755) != 0)
	    Error("Could not create directory %s\n", fn);
	}
    }
  if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("MPI_Barrier in Output() failed.\n");
  for (i = myFirstSlice; i <= myLastSlice; ++i)
    {
      nx = imageWidth[i] / spacing;
      ny = imageHeight[i] / spacing;
      nx2 = nx + 2;
      ny2 = ny + 2;
      Log("Going to output slice %d (%d) nx = %d ny = %d pts[i] = %llx\n",
	  i, actualSlice[i], nx, ny, pts[i]);

      if (outputMapName[0] != '\0')
	{
	  if (iter < 0)
	    sprintf(fn, "%sfinal/%0.*d.map", outputMapName, outputDigits,
		    actualSlice[i]);
	  else
	    sprintf(fn, "%s%0.6d/%0.*d.map", outputMapName, iter, outputDigits,
		    actualSlice[i]);
	  f = fopen(fn, "w");
	  if (f == NULL)
	    Error("Could not open output map file %s\n", fn);
	  mx = (float *) malloc(nx2 * ny2 * sizeof(float));
	  my = (float *) malloc(nx2 * ny2 * sizeof(float));
	  mc = (float *) malloc(nx2 * ny2 * sizeof(float));
	  memset(mx, 0, nx2 * ny2 * sizeof(float));
	  memset(my, 0, nx2 * ny2 * sizeof(float));
	  memset(mc, 0, nx2 * ny2 * sizeof(float));
	  pt = pts[i];
	  for (y = 0; y < ny; ++y)
	    for (x = 0; x < nx; ++x)
	      {
		mx[(y + 1) * nx2 + (x + 1)] = pt->x;
		my[(y + 1) * nx2 + (x + 1)] = pt->y;
		mc[(y + 1) * nx2 + (x + 1)] = 1.0;
		++pt;
	      }

	  /* extrapolate periphery */
	  for (x = 0; x < nx; ++x)
	    {
	      mx[x + 1] = 2.0 * mx[nx2 + x + 1]  - mx[2*nx2 + x + 1];
	      my[x + 1] = 2.0 * my[nx2 + x + 1]  - my[2*nx2 + x + 1];
	      mc[x + 1] = 1.0;

	      mx[(ny + 1) * nx2 + x + 1] =
		2.0 * mx[ny * nx2 + x + 1] - mx[(ny - 1) * nx2 + x + 1];
	      my[(ny + 1) * nx2 + x + 1] =
		2.0 * my[ny * nx2 + x + 1] - my[(ny - 1) * nx2 + x + 1];
	      mc[(ny + 1) * nx2 + x + 1] = 1.0;
	    }
	  for (y = 0; y < ny; ++y)
	    {
	      mx[(y + 1) * nx2] =
		2.0 * mx[(y + 1) * nx2 + 1]  - mx[(y + 1) * nx2 + 2];
	      my[(y + 1) * nx2] =
		2.0 * my[(y + 1) * nx2 + 1]  - my[(y + 1) * nx2 + 2];
	      mc[(y + 1) * nx2] = 1.0;

	      mx[(y + 1) * nx2 + nx + 1] =
		2.0 * mx[(y + 1) * nx2 + nx] -
		mx[(y + 1) * nx2 + nx - 1];
	      my[(y + 1) * nx2 + nx + 1] =
		2.0 * my[(y + 1) * nx2 + nx] -
		my[(y + 1) * nx2 + nx - 1];
	      mc[(y + 1) * nx2 + nx + 1] = 1.0;
	    }
	  mx[0] = 0.5 * (2.0 * mx[1] - mx[2] +
			 2.0 * mx[nx2] - mx[2*nx2]);
	  my[0] = 0.5 * (2.0 * my[1] - my[2] +
			 2.0 * my[nx2] - my[2*nx2]);
	  mc[0] = 1.0;
	  mx[nx + 1] = 0.5 * (2.0 * mx[nx] - mx[nx-1] +
			      2.0 * mx[nx2 + nx + 1] - mx[2*nx2 + nx + 1]);
	  my[nx + 1] = 0.5 * (2.0 * my[nx] - my[nx-1] +
			      2.0 * my[nx2 + nx + 1] - my[2*nx2 + nx + 1]);
	  mc[nx + 1] = 1.0;
	  mx[(ny + 1) * nx2] =
	    0.5 * (2.0 * mx[ny * nx2] - mx[(ny - 1) * nx2] +
		   2.0 * mx[(ny + 1) * nx2 + 1] - mx[(ny + 1) * nx2 + 2]);
	  my[(ny + 1) * nx2] =
	    0.5 * (2.0 * my[ny * nx2] - my[(ny - 1) * nx2] +
		   2.0 * my[(ny + 1) * nx2 + 1] - my[(ny + 1) * nx2 + 2]);
	  mc[(ny + 1) * nx2] = 1.0;
	  mx[(ny + 1) * nx2 + nx + 1] =
	    0.5 * (2.0 * mx[ny * nx2 + nx + 1] - mx[(ny - 1) * nx2 + nx + 1] +
		   2.0 * mx[(ny + 1) * nx2 + nx] - mx[(ny + 1) * nx2 + nx - 1]);
	  my[(ny + 1) * nx2 + nx + 1] =
	    0.5 * (2.0 * my[ny * nx2 + nx + 1] - my[(ny - 1) * nx2 + nx + 1] +
		   2.0 * my[(ny + 1) * nx2 + nx] - my[(ny + 1) * nx2 + nx - 1]);
	  mc[(ny + 1) * nx2 + nx + 1] = 1.0;

	  mapHeader[0] = nx;
	  mapHeader[1] = ny;
	  mapHeader[2] = actualSlice[i];
	  if (fwrite(mapHeader, sizeof(int), 3, f) != 3 ||
	      fwrite(mx, sizeof(float), nx2*ny2, f) != nx2*ny2 ||
	      fwrite(my, sizeof(float), nx2*ny2, f) != nx2*ny2 ||
	      fwrite(mc, sizeof(float), nx2*ny2, f) != nx2*ny2)
	    Error("Could not write map file %s\n", fn);
	  fclose(f);
	  free(mx);
	  free(my);
	  free(mc);
	}

      if (outputGridName[0] != '\0')
	{
	  grid = (unsigned char *) malloc(imageWidth[i] * imageHeight[i] * 3);
	  memset(grid, 0, imageWidth[i] * imageHeight[i] * 3);
	  pt = pts[i];
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
			
		      DrawLine(grid,
			       pt[y*nx+x].x, pt[y*nx+x].y,
			       pt[yp*nx+xp].x, pt[yp*nx+xp].y,
			       sliceMagnifications[i] * spacing *
			       ((dx * dy != 0) ? M_SQRT2 : 1.0),
			       imageWidth[i], imageHeight[i]);
		    }
		}

	  for (y = 0; y < ny; ++y)
	    for (x = 0; x < nx; ++x)
	      {
		if (pt[y*nx+x].fx > 0.5 * CONSTRAINED ||
		    pt[y*nx+x].fy < 0.5 * CONSTRAINED)
		  {
		    ix = (int) (pt[y*nx+x].x + 0.5);
		    iy = (int) (pt[y*nx+x].y + 0.5);
		    for (dy = -1; dy <= 1; ++dy)
		      for (dx = -1; dx <= 1; ++dx)
			{
			  xp = ix + dx;
			  yp = iy + dy;
			  if (xp >= 0 && xp < imageWidth[i] &&
			      yp >= 0 && yp < imageHeight[i])
			    {
			      grid[3*(yp * imageWidth[i] + xp)] = 255;
			      grid[3*(yp * imageWidth[i] + xp) + 1] = 255;
			      grid[3*(yp * imageWidth[i] + xp) + 2] = 255;
			    }
			}
		  }
	      }

	  if (iter < 0)
	    sprintf(fn, "%sfinal/%0.*d.pnm", outputGridName, outputDigits,
		    actualSlice[i]);
	  else
	    sprintf(fn, "%s%0.6d/%0.*d.pnm", outputGridName, iter, outputDigits,
		    actualSlice[i]);
	  f = fopen(fn, "w");
	  if (f == NULL)
	    Error("Could not open output grid file %s\n", fn);
	  fprintf(f, "P6\n%d %d\n%d\n", imageWidth[i], imageHeight[i], 255);
	  if (fwrite(grid, sizeof(unsigned char),
		     3*imageWidth[i]*imageHeight[i], f) !=
	      3*imageWidth[i]*imageHeight[i])
	    Error("Could not write to grid file %s\n", fn);
	  fclose(f);
	  free(grid);
	}
    }

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
  fprintf(logFile, "%f: ERROR: ", MPI_Wtime());
  vfprintf(logFile, fmt, args);
  va_end(args);
  fflush(logFile);
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
	 float nomD, int imageWidth, int imageHeight)
{
  int x0, y0, x1, y1;
  int steep;
  int tmp;
  int error;
  int yStep;
  int x, y;
  int dx, dy;
  float dist;
  unsigned char r, g, b;
  float hue;
  unsigned char *p;

  x0 = (int) (px0 + 0.5);
  y0 = (int) (py0 + 0.5);
  x1 = (int) (px1 + 0.5);
  y1 = (int) (py1 + 0.5);

  dist = hypot(x1 - x0, y1 - y0);
  hue = 4.0 * (dist / nomD) - 2.0;
  if (hue < 0.0)
    hue = 0.0;
  if (hue > 5.0)
    hue = 5.0;
  hsv_to_rgb(&r, &g, &b, hue, 1.0, 1.0);

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
UpdateDeltas (int iter)
{
  int i, j, k;
  int x, y;
  InterTableEntry *ite;
  float rrx, rry;
  Point *pt;
  float rx00, rx01, rx10, rx11;
  float ry00, ry01, ry10, ry11;
  float rx, ry;
  int nx, ny;

  for (i = 0; i < nInterTableEntries; ++i)
    {
      ite = &interTable[i];
      pt = pts[ite->s2];
      k = ite->pt2 - pt;
      nx = imageWidth[ite->s2] / spacing;
      ny = imageHeight[ite->s2] / spacing;
      y = k / nx;
      x = k - y * nx;
      rrx = ite->rrx;
      rry = ite->rry;
      if (rrx >= 0.0)
	{
	  if (x >= nx-1)
	    {
	      --x;
	      rrx += 1.0;
	    }
	}
      else
	{
	  if (x > 0)
	    {
	      --x;
	      rrx += 1.0;
	    }
	}
      if (rry >= 0.0)
	{
	  if (y >= ny-1)
	    {
	      --y;
	      rry += 1.0;
	    }
	}
      else
	{
	  if (y > 0)
	    {
	      --y;
	      rry += 1.0;
	    }
	}
      rx00 = pt[y*nx + x].x;
      rx01 = pt[(y+1)*nx+x].x;
      rx10 = pt[y*nx + x+1].x;
      rx11 = pt[(y+1)*nx+x+1].x;
      ry00 = pt[y*nx + x].y;
      ry01 = pt[(y+1)*nx+x].y;
      ry10 = pt[y*nx + x+1].y;
      ry11 = pt[(y+1)*nx+x+1].y;
      rx = rx00 * (rrx - 1.0) * (rry - 1.0)
	- rx10 * rrx * (rry - 1.0) 
	- rx01 * (rrx - 1.0) * rry
	+ rx11 * rrx * rry;
      ry = ry00 * (rrx - 1.0) * (rry - 1.0)
	- ry10 * rrx * (rry - 1.0) 
	- ry01 * (rrx - 1.0) * rry
	+ ry11 * rrx * rry;
      ite->dx = rx - ite->pt2->x;
      ite->dy = ry - ite->pt2->y;
#if DEBUG
      if (iter >= 999000)
      Log("new deltas (s = %d, x = %d y = %d): %d %d %f %f\n",
	  ite->s2, k/nx, k - (k/nx)*nx, x, y, ite->dx, ite->dy);
#endif
    }
}
