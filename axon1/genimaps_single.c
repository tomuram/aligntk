#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <errno.h>
#include <stdarg.h>

#include "par.h"
#include "imio.h"
#include "invert.h"

#define LINE_LENGTH	255

typedef struct Context {
  /* NOTE: any new fields added to this struct should
     be also added to PackContext and UnpackContext */ 
  char inputName[PATH_MAX];
  char masksName[PATH_MAX];
  char outputName[PATH_MAX];
} Context;

typedef struct Task {
  /* NOTE: any new fields added to this struct should
     be also added to PackTask and UnpackTask */ 
  char imageName[PATH_MAX];
} Task;

typedef struct Result {
  char message[PATH_MAX + 1024];
} Result;

typedef struct Node
{
  float x;			/* value at node (may be either black or white value */
  float fx;			/* force applied to value */
} Node;

typedef struct Spring
{
  unsigned int node0;		/* first node */
  unsigned int node1;		/* second node */
  float k;  			/* spring constant */
  float offset;			/* how much node1 should be greater than node0 */
} Spring;

/* GLOBAL VARIABLES FOR MASTER */
int nProcessed = 0;

/* GLOBAL VARIABLES FOR MASTER & WORKER */
Context c;
Task t;
Result r;

int level = 6;
float intraframeK = 1.0;
float levelK = 1.0;
float threshold = 1.0;

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
void Error (char *fmt, ...);
int CreateDirectories (char *fn);

int
main (int argc, char **argv, char **envp)
{
  par_process(argc, argv, envp,
              (void (*)()) MasterTask, MasterResult,
              WorkerContext, WorkerTask, NULL,
              PackContext, UnpackContext,
              PackTask, UnpackTask,
              PackResult, UnpackResult);
  return(0);
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
  int error;
  FILE *f;
  char inputDirName[PATH_MAX];
  char inputPrefix[PATH_MAX];
  int inputPrefixLen;
  char fn[PATH_MAX];
  int w, h;
  char imageListName[PATH_MAX];
  char line[LINE_LENGTH+1];
  int nItems;

  error = 0;
  c.inputName[0] = '\0';
  c.masksName[0] = '\0';
  c.outputName[0] = '\0';
  imageListName[0] = '\0';
  for (i = 1; i < argc; ++i)
    if (strcmp(argv[i], "-input") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-input error\n");
	    break;
	  }
	strcpy(c.inputName, argv[i]);
      }
    else if (strcmp(argv[i], "-masks") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-input error\n");
	    break;
	  }
	strcpy(c.masksName, argv[i]);
      }
    else if (strcmp(argv[i], "-output") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-output error\n");
	    break;
	  }
	strcpy(c.outputName, argv[i]);
      }
    else if (strcmp(argv[i], "-image_list") == 0)
      {
	if (++i == argc)
	  {
	    fprintf(stderr, "-image_list error\n");
	    error = 1;
	    break;
	  }
	strcpy(imageListName, argv[i]);
      }
    else
      {
	fprintf(stderr, "Unknown option: %s\n", argv[i]);
	error = 1;
      }
      

  if (error)
    {
      fprintf(stderr, "Usage: reduce -input file_prefix -output file_prefix\n");
      fprintf(stderr, "              -image_list image_list_file\n");
      fprintf(stderr, "              [-factor reduction_factor]\n");
      exit(1);
    }

  /* check that at least minimal parameters were supplied */
  if (c.inputName[0] == '\0' || c.outputName[0] == '\0' ||
      imageListName[0] == '\0')
    {
      fprintf(stderr, "-input, -output, and -image_list parameters must be specified.\n");
      exit(1);
    }

  par_set_context();

  /* for all slices */
  printf("Processing slices: ");
  fflush(stdout);

  /* read the images file */
  f = fopen(imageListName, "r");
  if (f == NULL)
    Error("Could not open file %s for reading\n", imageListName);
  while (fgets(line, LINE_LENGTH, f) != NULL)
    {
      len = strlen(line);
      if (len > 0 && line[len-1] == '\n')
	line[len-1] = '\0';
      if (line[0] == '\0' || line[0] == '#')
	continue;
      nItems = sscanf(line, "%s", t.imageName);
      if (nItems != 1)
	Error("Malformed line in %s:\n%s\n", imageListName, line);
      par_delegate_task();
    }
  fclose(f);
  par_finish();

  printf(" %d\nAll slices completed.\n", nProcessed);
  exit(0);
}

void
MasterResult ()
{
  if (r.message[0] != '\0')
    {
      fprintf(stderr, "\nThe following error was encountered by one of the worker processes:\n%s\n", r.message);
      exit(1);
    }
  if ((nProcessed % 50) == 0 && nProcessed != 0)
    printf(" %d \n                   ", nProcessed);
  printf(".");
  fflush(stdout);
  ++nProcessed;
}


/* WORKER PROCEDURES */

void
WorkerContext ()
{
  /* nothing necessary here */
}

void
WorkerTask ()
{
  int imageWidth, imageHeight;
  unsigned char *image;
  unsigned char *mask;
  char sdirName[PATH_MAX];
  char fn[PATH_MAX];
  FILE *f;
  int nRows, nCols;
  int row, col, quad;
  int subrow, subcol;
  float blackLevel, whiteLevel;
  int nxNodes, nyNodes;
  Node *nodes;
  int nSprings;
  int springsSize;
  Spring *springs;
  int iw, ih;
  int mw, mh;
  int mbpl;
  int ix, iy;
  int dx, dy;
  int x, y;
  int i;
  int frameBlack, frameGray, frameWhite;
  int nValidPixels;
  int requiredPixels;
  int cumulativePixels;
  double sum;
  int nPixels;
  Spring *s;
  double prevEnergy;
  double energy;
  int iter;
  double delta;
  double force;
  double maxF;
  double scale;
  int epochIterations;
  double epochInitialEnergy;
  short *histograms;
  char msg[PATH_MAX + 1024];
  short *h;
  int fh[256], nh[256];
  int mode;
  int terminationRequestedIter;
  double dampingFactor;
  double maxStep;
  double deltaEnergy;
  int nDecrease, nIncrease;
  char imgName[PATH_MAX], refName[PATH_MAX];
  int spacing;
  int histogram[256];
  long long sectionHistogram[256];
  int sectionBlack, sectionWhite;
  int range;
  int nLevelSpringsBlack, nLevelSpringsWhite;
  int nIntraSpringsBlack, nIntraSpringsWhite;
  MapElement *map;

  /* get the image */
  sprintf(fn, "%s%s", c.inputName, t.imageName);
  if (!ReadImage(fn, &image, &imageWidth, &imageHeight, -1, -1, -1, -1, msg))
    Error("Could not read image file %s:\n%s\n", fn, msg);
  mbpl = (imageWidth + 7) / 8;

  /* get the mask */
  if (c.masksName[0] != '\0')
    {
      sprintf(fn, "%s%s.pbm", c.masksName, t.imageName);
      if (!ReadBitmap(fn, &mask, &mw, &mh, -1, -1, -1, -1, msg))
	{
	  fprintf(stderr, "Could not read mask from file %s:\n%s\n",
		  fn, msg);
	  exit(1);
	}		      
      if (mw != imageWidth || mh != imageHeight)
	{
	  fprintf(stderr, "Unexpected mask size in file %s\n", fn);
	  exit(1);
	}
    }
  else
    {
      mask = (unsigned char *) malloc(imageHeight * mbpl);
      memset(mask, 0xff, imageHeight * mbpl);
    }

  spacing = 1 << level;
  nxNodes = (imageWidth + spacing-1) / spacing + 1;
  nyNodes = (imageHeight + spacing-1) / spacing + 1;
  nodes = (Node*) malloc(2 * nyNodes * nxNodes * sizeof(Node));
  nLevelSpringsBlack = 0;
  nLevelSpringsWhite = 0;
  nIntraSpringsBlack = 0;
  nIntraSpringsWhite = 0;
  nSprings = 0;
  springsSize = 0;
  springs = NULL;
  histograms = (short *) malloc((nyNodes - 1) * (nxNodes - 1) * 256 * sizeof(short));

  memset(sectionHistogram, 0, 256*sizeof(long long));
  for (y = 0; y < imageHeight; ++y)
    for (x = 0; x < imageWidth; ++x)
      if (mask[y*mbpl + (x >> 3)] & (0x80 >> (x & 7)))
	++sectionHistogram[image[y*imageWidth + x]];

  /* choose section-wide black and white levels */
  mode = 0;
  for (i = 1; i < 256; ++i)
    if (sectionHistogram[i] > sectionHistogram[mode])
      mode = i;
  /* the black value will be the highest value with fewer than 1% of
     number of pixels of the mode */
  for (i = mode; i > 0; --i)
    if (sectionHistogram[i] < 0.01 * sectionHistogram[mode])
      break;
  sectionBlack = i;
  for (i = mode; i < 255; ++i)
    if (sectionHistogram[i] < 0.01 * sectionHistogram[mode])
      break;
  sectionWhite = i;

  /* increase coverage by 20% in each direction to ensure we cover
     all potentially valid pixels */
  range = sectionWhite - sectionBlack;
  sectionBlack = sectionBlack - range / 5;
  if (sectionBlack < 0)
    sectionBlack = 0;
  sectionWhite = sectionWhite + range / 5;
  if (sectionWhite > 255)
    sectionWhite = 255;
  printf("section black = %d  section white = %d\n",
	 sectionBlack, sectionWhite);

  /* make a histogram for each spacing x spacing pixel region of the image */
  memset(histograms, 0, (nyNodes - 1) * (nxNodes - 1) * 256 * sizeof(short));
  for (iy = 0; iy < nyNodes-1; ++iy)
    for (ix = 0; ix < nxNodes-1; ++ix)
      {
	h = &histograms[(iy * (nxNodes-1) + ix) * 256];
	for (dy = 0; dy < spacing; ++dy)
	  {
	    y = iy *spacing + dy;
	    if (y >= imageHeight)
	      break;
	    for (dx = 0; dx < spacing; ++dx)
	      {
		x = ix * spacing + dx;
		if (x >= imageWidth)
		  break;
		if (mask[y*mbpl+(x >> 3)] & (0x80 >> (x & 7)))
		  ++h[image[y*iw+x]];
	      }
	  }
      }


  /* now make a histogram for the entire frame */
  memset(fh, 0, 256*sizeof(int));
  for (iy = 0; iy < nyNodes-1; ++iy)
    for (ix = 0; ix < nxNodes-1; ++ix)
      {
	h = &histograms[(iy * (nxNodes-1) + ix) * 256];
	for (i = 0; i < 256; ++i)
	  fh[i] += h[i];
      }
  
  /* choose image-wide black and white levels */
  /* first find the mode */
  mode = 0;
  for (i = 1; i < 256; ++i)
    if (fh[i] > fh[mode])
      mode = i;
  /* the black value will be the highest value < mode that
     is used < 1% of the mode in the frame */
  for (i = mode; i > 0; --i)
    if (fh[i] < 0.01 * fh[mode])
      break;
  frameBlack = i;
  if (frameBlack < sectionBlack)
    frameBlack = sectionBlack;
  /* the white value will be the lowest value > mode that
     is used < 1% of the mode in the frame */
  for (i = mode; i < 255; ++i)
    if (fh[i] < 0.01 * fh[mode])
      break;
  frameWhite = i;
  if (frameWhite > sectionWhite)
    frameWhite = sectionWhite;
  frameGray = (frameBlack + frameWhite) / 2;
  
  /* go through all node positions */
  for (iy = 0; iy < nyNodes; ++iy)
    for (ix = 0; ix < nxNodes; ++ix)
      {
	/* set the initial values to the frame-wide values */
	nodes[iy * nxNodes + ix].x = frameBlack;
	nodes[nyNodes*nxNodes + (iy * nxNodes + ix)].x = frameWhite;

	/* make a spring to adjacent nodes */
	for (dy = -1; dy <= 0; ++dy)
	  for (dx = -1; dx <= 1; ++dx)
	    {
	      if (dx == 0 && dy == 0)
		break;
	      x = ix + dx;
	      y = iy + dy;
	      if (x < 0 || x >= nxNodes ||
		  y < 0)
		continue;
	      
	      if (nSprings+2 > springsSize)
		{
		  springsSize = (springsSize > 0) ? (springsSize + 1024*1024) : 1024*1024;
		  springs = (Spring*) realloc(springs, springsSize * sizeof(Spring));
		}
	      /* spring for black level */
	      springs[nSprings].node0 = y * nxNodes + x;
	      springs[nSprings].node1 = iy * nxNodes + ix;
	      springs[nSprings].k = intraframeK;
	      springs[nSprings].offset = 0.0;
	      ++nIntraSpringsBlack;
	      ++nSprings;
	      /* spring for white level */
	      springs[nSprings].node0 = nyNodes*nxNodes + (y * nxNodes + x);
	      springs[nSprings].node1 = nyNodes*nxNodes + (iy * nxNodes + ix);
	      springs[nSprings].k = intraframeK;
	      springs[nSprings].offset = 0.0;
	      ++nIntraSpringsWhite;
	      ++nSprings;
	    }

	/* compute the histogram of intensity values in a 256x256 square region
	   surrounding this node */
	memset(nh, 0, 256*sizeof(int));
	for (dy = -2; dy < 2; ++dy)
	  {
	    y = iy + dy;
	    if (y < 0 || y >= nyNodes-1)
	      continue;
	    for (dx = -2; dx < 2; ++dx)
	      {
		x = ix + dx;
		if (x < 0 || x >= nxNodes-1)
		  continue;
		h = &histograms[(y * (nxNodes-1) + x) * 256];
		for (i = 0; i < 256; ++i)
		  nh[i] += h[i];
	      }				
	  }

	/* choose the black value as the average of the lowest 1% of
	   valid pixels that are darker than gray; 2048 pixels
	   are required to consider this a usable estimate */
	nValidPixels = 0;
	for (i = frameBlack; i <= frameGray; ++i)
	  nValidPixels += nh[i];
	if (nValidPixels >= 2048)
	  {
	    requiredPixels = nValidPixels / 100;
	    cumulativePixels = 0;
	    sum = 0.0;
	    for (i = frameBlack; i <= frameGray; ++i)
	      {
		cumulativePixels += nh[i];
		if (cumulativePixels >= requiredPixels)
		  {
		    nPixels = nh[i] - (cumulativePixels - requiredPixels);
		    sum += nPixels * i;
		    break;
		  }
		nPixels = nh[i];
		sum += nPixels * i;
	      }
	    blackLevel = sum / requiredPixels;
	    if (nSprings+1 > springsSize)
	      {
		springsSize = (springsSize > 0) ? (springsSize + 1024*1024) : 1024*1024;
		springs = (Spring*) realloc(springs, springsSize * sizeof(Spring));
	      }
	    /* spring for black level */
	    springs[nSprings].node0 = iy * nxNodes + ix;
	    springs[nSprings].node1 = 0xffffffff;
	    springs[nSprings].k = levelK;
	    springs[nSprings].offset = blackLevel;
	    ++nLevelSpringsBlack;
	    ++nSprings;
	    
	    /* also set initial value */
	    nodes[iy * nxNodes + ix].x = blackLevel;
	  }

			
	/* choose the white value as the average of the top 1% of
	   valid pixels that are lighter than gray; 2048 pixels
	   are required to consider this a usable estimate */
	nValidPixels = 0;
	for (i = frameGray+1; i <= frameWhite; ++i)
	  nValidPixels += nh[i];
	if (nValidPixels >= 2048)
	  {
	    requiredPixels = nValidPixels / 100;
	    cumulativePixels = 0;
	    sum = 0.0;
	    for (i = frameWhite; i > frameGray; --i)
	      {
		cumulativePixels += nh[i];
		if (cumulativePixels >= requiredPixels)
		  {
		    nPixels = nh[i] - (cumulativePixels - requiredPixels);
		    sum += nPixels * i;
		    break;
		  }
		nPixels = nh[i];
		sum += nPixels * i;
	      }
	    whiteLevel = sum / requiredPixels;
	    if (nSprings+1 > springsSize)
	      {
		springsSize = (springsSize > 0) ? (springsSize + 1024*1024) : 1024*1024;
		springs = (Spring*) realloc(springs, springsSize * sizeof(Spring));
	      }
	    /* spring for white level */
	    springs[nSprings].node0 = nyNodes*nxNodes + (iy * nxNodes + ix);
	    springs[nSprings].node1 = 0xffffffff;
	    springs[nSprings].k = levelK;
	    springs[nSprings].offset = whiteLevel;
	    ++nLevelSpringsWhite;
	    ++nSprings;
	    /* also set initial value */
	    nodes[nyNodes*nxNodes + (iy * nxNodes + ix)].x = whiteLevel;
	  }
      }

  printf("nLevelSpringsBlack = %d nLevelSpringsWhite = %d\n",
	 nLevelSpringsBlack, nLevelSpringsWhite);
  printf("nIntraSpringsBlack = %d nIntraSpringsWhite = %d\n",
	 nIntraSpringsBlack, nIntraSpringsWhite);

  /* relax the spring system */
  dampingFactor = 0.1;
  epochIterations = 512;
  epochInitialEnergy = 1.0e+30;
  prevEnergy = 1.0e+30;
  terminationRequestedIter = -1;
  nIncrease = 0;
  nDecrease = 0;
  for (iter = 0; ; ++iter)
    {
      /* initialize all forces to zero */
      energy = 0.0;
      for (iy = 0; iy < nyNodes; ++iy)
	for (ix = 0; ix < nxNodes; ++ix)
	  {
	    nodes[iy * nxNodes + ix].fx = 0.0;
	    nodes[nyNodes*nxNodes + (iy * nxNodes + ix)].fx = 0.0;
	  }

      /* update all forces, also computing energy */
      for (i = 0; i < nSprings; ++i)
	{
	  s = &springs[i];
	  if (s->node1 == 0xffffffff)
	    {
	      delta = s->offset - nodes[s->node0].x;
	      force = s->k * delta;
	      nodes[s->node0].fx += force;
	      energy += force * delta;
	    }
	  else
	    {
	      delta = nodes[s->node1].x - nodes[s->node0].x - s->offset;
	      force = s->k * delta;
	      nodes[s->node0].fx += force;
	      nodes[s->node1].fx -= force;
	      energy += force * delta;
	    }
	}

      maxF = 0.0;
      for (iy = 0; iy < nyNodes; ++iy)
	for (ix = 0; ix < nxNodes; ++ix)
	  {
	    if (nodes[iy * nxNodes + ix].fx > maxF)
	      maxF = nodes[iy * nxNodes + ix].fx;
	    if (nodes[nyNodes*nxNodes + (iy * nxNodes + ix)].fx > maxF)
	      maxF = nodes[nyNodes*nxNodes + (iy * nxNodes + ix)].fx;
	  }

      /* update all positions */
      if (maxF > 0.5)
	scale = dampingFactor * 0.5 / maxF;
      else
	scale = dampingFactor;
      maxStep = 0.1;
      for (iy = 0; iy < nyNodes; ++iy)
	for (ix = 0; ix < nxNodes; ++ix)
	  {
	    delta = scale * nodes[iy * nxNodes + ix].fx;
	    if (delta > maxStep)
	      delta = maxStep;
	    nodes[iy * nxNodes + ix].x += delta;
	    delta = scale * nodes[nyNodes*nxNodes + (iy * nxNodes + ix)].fx;
	    if (delta > maxStep)
	      delta = maxStep;
	    nodes[nyNodes*nxNodes + (iy * nxNodes + ix)].x += delta;
	  }

      if (iter % 100 == 0)
	printf("After %d iterations, total energy is %f  (df = %f mgf = %f)\n",
	       iter, energy, dampingFactor, maxF);
      deltaEnergy = energy - prevEnergy;
      prevEnergy = energy;
      if (iter % epochIterations == 0)
	epochInitialEnergy = energy;

      // check periodically for termination conditions
      if (iter % epochIterations == epochIterations-1)
	{
	  if (epochInitialEnergy - energy <
	      energy * threshold * 0.000001 * epochIterations ||
	      energy < 0.000001)
	    terminationRequestedIter = iter;
	}

      if (deltaEnergy > 0.0)
	{
	  nDecrease = 0;
	  ++nIncrease;
	  if (nIncrease > 1000 ||
	      nIncrease > 10 && deltaEnergy > energy ||
	      deltaEnergy > 1000.0 * energy)
	    {
	      fprintf(stderr, "System is diverging instead of converging.\n");
	      exit(1);
	    }
	      
	  if (iter < 100 || nIncrease < 2)
	    {
	      dampingFactor *= 0.5;
	      if (dampingFactor < 0.000000001)
		{
		  fprintf(stderr, "dampingFactor became too small\n");
		  exit(1);
		}
	    }
	}
      else
	{
	  ++nDecrease;
	  nIncrease = 0;
	  if (dampingFactor < 0.5)
	    dampingFactor *= 1.01;
	}

      /* 32 is chosen as a preferred time to terminate because,
	 given that we increase the damping by 1% on each iteration, and decrease
	 it by a factor of 2 when instability sets in, the number of iterations
	 between instabilities is about 70, and we want to terminate
	 when the state is not near an instability */
      if (terminationRequestedIter >= 0 &&
	  (nDecrease == 32 || iter > terminationRequestedIter + 128))
	{
	  printf("Termination (requested at iter %d, nDecrease = %d)\n",
		 terminationRequestedIter, nDecrease);
	  break;
	}
    }
  printf("Finished intensity adjustment at iteration %d.\n", iter);

  /* output the intensity map */
  map = (MapElement*) malloc(nxNodes * nyNodes * sizeof(MapElement));
  for (iy = 0; iy < nyNodes; ++iy)
    for (ix = 0; ix < nxNodes; ++ix)
      {
	map[iy*nxNodes+ix].x = nodes[iy*nxNodes + ix].x;
	map[iy*nxNodes+ix].y = nodes[nyNodes*nxNodes + iy*nxNodes + ix].x;
	map[iy*nxNodes+ix].c = 1.0;
      }

  /* write out the map */
  sprintf(fn, "%s%s.map", c.outputName, t.imageName);
  if (!CreateDirectories(fn))
    {
      fprintf(stderr, "Could not create directories for %s\n", fn);
      exit(1);
    }
  sprintf(imgName, "%s%s", c.inputName, t.imageName);
  if (!WriteMap(fn, map, level, nxNodes, nyNodes, 0, 0, imgName, imgName,
		UncompressedMap, msg))
    {
      fprintf(stderr, "Could not write intensity map %s:\n%s\n",
	      fn, msg);
      exit(1);
    }
  free(map);
  free(nodes);
  free(springs);
  free(histograms);
  free(image);
  free(mask);

  r.message[0] = '\0';
}

void
PackContext ()
{
  par_pkstr(c.inputName);
  par_pkstr(c.masksName);
  par_pkstr(c.outputName);
}

void
UnpackContext ()
{
  par_upkstr(c.inputName);
  par_upkstr(c.masksName);
  par_upkstr(c.outputName);
}

void
PackTask ()
{
  par_pkstr(t.imageName);
}

void
UnpackTask ()
{
  par_upkstr(t.imageName);
}

void
PackResult ()
{
  par_pkstr(r.message);
}

void
UnpackResult ()
{
  par_upkstr(r.message);
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
	      fprintf(stderr, "Directory hash table is full!\n");
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
	  fprintf(stderr, "Output path component %s is not a directory\n", dn);
	  return(0);
	}
      if (errno != ENOENT)
	{
	  fprintf(stderr, "Could not stat directory %s\n", dn);
	  return(0);
	}
      
      if (mkdir(dn, 0755) != 0)
	{
	  fprintf(stderr, "Could not create directory %s\n", dn);
	  return(0);
	}
    }
  return(1);
}
