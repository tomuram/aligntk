/*
 * genimaps.c  - generate a set of intensity maps using a relaxation
 *                algorithm on a system of springs
 *
 *  Copyright (c) 2009-2013 National Resource for Biomedical
 *                          Supercomputing,
 *                          Pittsburgh Supercomputing Center,
 *                          Carnegie Mellon University
 *
 *  This file is part of AlignTK.
 *
 *  AlignTK is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  AlignTK is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with AlignTK.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Acknowledgements:
 *     Development of this code was supported in part by
 *       NIH NCRR grant 5P41RR006009 and
 *       NIH NIGMS grant P41GM103712
 *
 *  HISTORY
 *    2009  Development of intrasection_genimaps.c (ghood@psc.edu)
 *    2011  Generalized and embedded in a framework similar
 *              to that of align12.c (ghood@psc.edu)
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
#include "invert.h"
#include "dt.h"

#define DEBUG	0
#define PDEBUG	0

#define LINE_LENGTH	256

typedef struct Node
{
  float x;			/* value at node (may be either black or white value */
  float fx;			/* force applied to value */
} Node;

typedef struct Spring
{
  Node *node0;		/* first node */
  Node *node1;		/* second node */
  float k;  		/* spring constant */
  float offset;		/* how much node1 should be greater than node0 */
} Spring;

typedef struct Image
{
  int next;		/* index of next image in hash bucket */
  char *name;   	/* name of this image */
  int width, height;    /* width and height in pixels */
  int owner;		/* which process owns the image */
  int needed;		/* true if this image is needed in this process */ 
  unsigned char *sendTo;/* bitmask of the processes that this image should
			   be sent to */
  int nx, ny;           /* number of nodes in x and y */
  int black;		/* black level for this image */
  int gray;		/* gray level for this image */
  int white; 		/* white level for this image */
  Node *nodes;          /* nodes in the image */
  unsigned char *image; /* image pixels */
  unsigned char *mask;  /* mask of valid pixels */
  unsigned int *histogram;   /* histogram of pixel values for the entire image */
  unsigned short *nodeHistograms;
} Image;

typedef struct InterImageMap
{
  char *name;           /* filename of map */
  int image0;           /* index of source image */
  int image1;           /* index of target (reference) image */
  float energyFactor;   /* either 0.0 or 1.0; determines whether this entry
			   contributes to the energy sum on this node */
  float k;		/* overall spring constant for this map */
} InterImageMap;

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

char imagesName[PATH_MAX];
char imageListName[PATH_MAX];
char mapListName[PATH_MAX];
char mapsName[PATH_MAX];
char masksName[PATH_MAX];
char outputName[PATH_MAX];
int epochIterations = 512;

int nImages = 0;
Image *images = 0;
int *imageHashTable = 0;
int myFirstImage, myLastImage;

int nMaps = 0;
InterImageMap *maps = 0;
int mapsSize = 0;
int mapsPos = 0;

int nPhases;
CommPhase *commPhases;
int bufferSize = 0;
float *buffer = 0;

int springsSize = 0;
int nSprings = 0;
Spring *springs = NULL;

FILE *logFile = 0;

float levelK = 1.0;
float intraimageK = 1.0;
float interimageK = 10.0;
float threshold = 1.0;  /* ppm change per iteration */

int level = 6;
int spacing;		/* 2^level = pixels between intensity map points in x and y */

#define DIR_HASH_SIZE	8192
char *dirHash[DIR_HASH_SIZE];

void Error (char *fmt, ...);
void Log (char *fmt, ...);
void PlanCommunications ();
void CommunicateIntensities ();
unsigned int Hash (char *s);
unsigned int HashMap (char *s, int nx, int ny);
int CreateDirectories (char *fn);


int
main (int argc, char **argv)
{
  char fn[PATH_MAX];
  FILE *f;
  int x, y, z;
  int i, j, k;
  int error;
  int im;
  int n;
  float rx, ry, rc;
  int irx, iry;
  float rrx, rry;
  float xv, yv;
  int dx, dy;
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
  float maxF;
  float globalMaxF;
  char hostName[256];
  double basis;
  int nx, ny, nz;
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
  unsigned char *image;
  unsigned char *image0, *image1;
  int i0, i1;
  char msg[PATH_MAX+256];
  MapElement *map;
  InverseMap *iMap;
  int mFactor;
  int prevX, prevY;
  InterImageMap *m;
  Node *nodes0, *nodes1;
  Node *nodes;
  int phase;
  int op;
  CommPhase *cp;
  int separation;
  int mp;
  int *phaseCount;
  int *globalPhaseCount;
  int terminationRequestedIter;
  int ns;
  int mx, my;
  float msk;
  float sk;
  size_t mbpl, mbpl0, mbpl1;
  unsigned char *mask, *mask0, *mask1;
  int minIter;
  int rnx, rny;
  int valid;
  int ix, iy;
  int count;
  int maskWidth, maskHeight;
  int ixMin, ixMax, iyMin, iyMax;

  unsigned short *nh;
  unsigned int *ih;
  unsigned int rh[256];
  double dgh[256];
  double globalHistogram[256];

  long long globalTotal;
  long long imageTotal;
  long long total;
  int globalBlack, globalWhite;
  int imageWidth, imageHeight;
  unsigned short *nodeHistograms;
  int mode;
  int range;
  int imageBlack, imageGray, imageWhite;
  int nLevelSpringsBlack, nLevelSpringsWhite;
  int nIntraSpringsBlack, nIntraSpringsWhite;
  int nInterSpringsBlack, nInterSpringsWhite;
  int nValidPixels;
  int requiredPixels;
  int cumulativePixels;

  int nx0, ny0, nx1, ny1;
  int iw0, ih0, iw1, ih1;
  float iv, rv;
  float r00, r01, r10, r11;
  double sumBlack0, sumBlack1, sumWhite0, sumWhite1;
  double weightBlack0, weightBlack1, weightWhite0, weightWhite1;
  int blackCount0, whiteCount0, blackCount1, whiteCount1;
  int image0Black, image0Gray, image0White;
  int image1Black, image1Gray, image1White;
  float diff0to1;
  int ix1, iy1;
  float maxStep;
  char imgName[PATH_MAX];
  float blackLevel, whiteLevel;
  double delta;
  int nPixels;
  double sum;
  float x0, y0;
  Spring *s;

  int ni, nix, niy, nlv;
  float lvl0, lvl1, lvl;

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
  
  sprintf(fn, "logs/%.2d.log", p);
  logFile = fopen(fn, "w");
  gethostname(hostName, 255);
  Log("Process %d is running on %s\n", p, hostName);
  memset(dirHash, 0, DIR_HASH_SIZE * sizeof(char*));
  if (p == 0)
    {
      error = 0;
      imagesName[0] = '\0';
      imageListName[0] = '\0';
      mapListName[0] = '\0';
      mapsName[0] = '\0';
      masksName[0] = '\0';
      outputName[0] = '\0';

      for (i = 0; i < argc; ++i)
	Log("ARGV[%d] = %s\n", i, argv[i]);
      for (i = 1; i < argc; ++i)
	if (strcmp(argv[i], "-images") == 0)
	  {
	    if (++i == argc)
	      {
		error = 1;
		break;
	      }
	    strcpy(imagesName, argv[i]);
	  }
	else if (strcmp(argv[i], "-image_list") == 0)
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
	else if (strcmp(argv[i], "-level") == 0)
	  {
	    if (++i == argc ||
		sscanf(argv[i], "%d", &level) != 1)
	      {
		error = 1;
		break;
	      }
	  }
	else
	  {
	    fprintf(stderr, "Unrecognized option: %s\n", argv[i]);
	    error = 1;
	  }

      if (error)
	{
	  fprintf(stderr, "Usage: genimaps\n");
	  fprintf(stderr, "              -image_list images_file\n");
	  fprintf(stderr, "              -images images_prefix\n");
	  fprintf(stderr, "              -output output_prefix\n");
	  fprintf(stderr, "             [-map_list maps_file]\n");
	  exit(1);
	}
      
	
      /*    maps_file contains one line for each map between images:
            image_name0 image_name1 map_file weight
            image_name0 image_name1 map_file weight
            ...
      */

      /* check that at least minimal parameters were supplied */
      if (imageListName[0] == '\0')
	Error("-image_list must be specified\n");
      if (imagesName[0] == '\0')
	Error("-images must be specified.\n");
      if (outputName[0] == '\0')
	Error("-output must be specified.\n");
    }

  /* broadcast the info */
  if (MPI_Bcast(imagesName, PATH_MAX, MPI_CHAR, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(mapsName, PATH_MAX, MPI_CHAR, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(masksName, PATH_MAX, MPI_CHAR, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(outputName, PATH_MAX, MPI_CHAR, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&level, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("Broadcast of parameters failed.\n");
  spacing = 1 << level;

  if (p == 0)
    {
      /* read the images file */
      Log("imageListName = %s\n", imageListName);
      f = fopen(imageListName, "r");
      if (f == NULL)
	Error("Could not open file %s for reading\n", imageListName);
      nImages = 0;
      imageNamesSize = 0;
      imageNamesPos = 0;
      imageNames = 0;
      while (fgets(line, LINE_LENGTH, f) != NULL)
	{
	  if (line[0] == '\0' || line[0] == '#')
	    continue;
	  width = -1;
	  height = -1;
	  nItems = sscanf(line, "%s", imageName);
	  if (nItems != 1)
	    Error("Malformed line in %s:\n%s\n", imageListName, line);

	  imageNameLen = strlen(imageName);
	  while (imageNamesPos + imageNameLen + 1 >= imageNamesSize)
	    {
	      imageNamesSize = (imageNamesSize > 0) ? imageNamesSize * 2 : 1024;
	      imageNames = (char *) realloc(imageNames, imageNamesSize);
	    }
	  strcpy(&imageNames[imageNamesPos], imageName);
	  imageNamesPos += imageNameLen + 1;

	  ++nImages;
	}
      fclose(f);

      imageNamesSize = imageNamesPos;
      imageNames = (char *) realloc(imageNames, imageNamesSize);
    }

  /* broadcast the image info */
  if (MPI_Bcast(&nImages, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&imageNamesSize, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("Broadcast of nImages or imageNamesSize failed.\n");

  if (p != 0)
    imageNames = (char *) malloc(imageNamesSize * sizeof(char));
  if (MPI_Bcast(imageNames, imageNamesSize, MPI_CHAR, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("Broadcast of imageNames failed.\n");

  images = (Image *) malloc(nImages * sizeof(Image));
  imageHashTable = (int *) malloc(nImages * sizeof(int));
  for (i = 0; i < nImages; ++i)
    imageHashTable[i] = -1;

  imageNamesPos = 0;
  for (i = 0; i < nImages; ++i)
    {
      images[i].name = &imageNames[imageNamesPos];
      imageNamesPos += strlen(&imageNames[imageNamesPos]) + 1;
      images[i].width = 0;
      images[i].height = 0;
      images[i].owner = -1;
      images[i].needed = 0;
      images[i].sendTo = NULL;
      images[i].nx = 0;
      images[i].ny = 0;
      images[i].black = 0;
      images[i].gray = 0;
      images[i].white = 0;
      images[i].nodes = 0;
      images[i].image = NULL;
      images[i].mask = NULL;
      images[i].histogram = NULL;
      images[i].nodeHistograms = NULL;
      hv = Hash(images[i].name) % nImages ;
      images[i].next = imageHashTable[hv];
      imageHashTable[hv] = i;
    }

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

  Log("On node %d first = %d last = %d (nz = %d)\n",
      p, myFirstImage, myLastImage, nImages);

  if (p == 0)
    {
      nMaps = 0;
      mapNamesSize = 0;
      mapNamesPos = 0;
      mapNames = 0;
      mapParamsSize = 0;
      mapParamsPos = 0;
      mapParams = 0;
      if (mapListName[0] != '\0')
	{
	  /* read the maps file */
	  f = fopen(mapListName, "r");
	  if (f == NULL)
	    Error("Could not open file %s for reading\n", mapListName);
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
      i0 = -1;
      hv = Hash(imageName0) % nImages;
      for (j = imageHashTable[hv]; j >= 0; j = images[j].next)
	if (strcmp(images[j].name, imageName0) == 0)
	  {
	    i0 = j;
	    if (j >= myFirstImage && j <= myLastImage)
	      found = 1;
	    break;
	  }
      if (i0 < 0)
	Error("Could not find source image for map %s\n", pairName);

      /* check if this node holds the destination
	 of the map */
      i1 = -1;
      hv = Hash(imageName1) % nImages;
      for (j = imageHashTable[hv]; j >= 0; j = images[j].next)
	if (strcmp(images[j].name, imageName1) == 0)
	  {
	    i1 = j;
	    if (j >= myFirstImage && j <= myLastImage)
	      found = 1;
	    break;
	  }
      if (i1 < 0)
	Error("Could not find destination image for map %s\n", pairName);

      if (found)
	{
	  images[i0].needed = 1;
	  images[i1].needed = 1;

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
	  m->image0 = i0;
	  m->image1 = i1;
	  m->energyFactor = (i0 >= myFirstImage &&
			     i0 <= myLastImage) ? 1.0 : 0.0;
	  m->k = weight;

	  if (images[i0].owner != images[i1].owner)
	    {
	      if (images[i0].owner == p)
		{
		  op = images[i1].owner;
		  images[i0].sendTo[op >> 3] |= 0x80 >> (op & 7);
		}
	      else
		{
		  op = images[i0].owner;
		  images[i1].sendTo[op >> 3] |= 0x80 >> (op & 7);
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
	nFloats = 2 * images[i].nx * images[i].ny;
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
	nFloats = 2 * images[i].nx * images[i].ny;
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

  /* read the images and masks; also build the histograms */
  Log("Reading images and masks\n");
  Log("spacing = %d\n", spacing);
  memset(dgh, 0, 256*sizeof(double));
  for (i = 0; i < nImages; ++i)
    {
      if (images[i].owner != p && !images[i].needed)
	continue;
      sprintf(fn, "%s%s", imagesName, images[i].name);
      Log("Reading image %s\n", fn);
      if (!ReadImage(fn, &image, &imageWidth, &imageHeight,
		     -1, -1, -1, -1, msg))
	Error("Could not read image file %s:\n%s\n", fn, msg);
      images[i].image = image;
      images[i].width = imageWidth;
      images[i].height = imageHeight;
      mask = NULL;
      if (masksName[0] != '\0')
	{
	  /* read the mask from a file */
	  sprintf(fn, "%s%s.pbm", masksName, images[i].name);
	  Log("Reading mask %s\n", fn);
	  if (ReadBitmap(fn, &mask, &maskWidth, &maskHeight,
			 -1, -1, -1, -1,
			 msg))
	    {
	      Log("Read mask %s  -- width=%d height=%d\n",
		  fn, maskWidth, maskHeight);
	      if (maskWidth != imageWidth ||
		  maskHeight != imageHeight)
		Error("Mask %s is not same size as image %s  %d %d %d %d\n",
		      fn, images[i].name,
		      maskWidth, maskHeight, imageWidth, imageHeight);
	      mbpl = (maskWidth + 7) >> 3;
	    }
	  else
	    Log("WARNING: could not read mask %s -- assuming unmasked image.\n",
		fn);
	}
      if (mask == NULL)
	{
	  maskWidth = images[i].width;
	  maskHeight = images[i].height;
	  mbpl = (maskWidth + 7) >> 3;
	  mask = (unsigned char *) malloc(maskHeight * mbpl);
	  memset(mask, 0xff, maskHeight * mbpl);
	}
      images[i].mask = mask;
      images[i].nx = nx = (imageWidth + spacing-1) / spacing + 1;
      images[i].ny = ny = (imageHeight + spacing-1) / spacing + 1;

      /* make a histogram for each node -- a (spacing x spacing) pixel region of the image */
      images[i].nodeHistograms = nodeHistograms =
	(unsigned short *) malloc((ny - 1) * (nx - 1) * 256 * sizeof(unsigned short));
      memset(nodeHistograms, 0, (ny - 1) * (nx - 1) * 256 * sizeof(unsigned short));
      for (iy = 0; iy < ny-1; ++iy)
	for (ix = 0; ix < nx-1; ++ix)
	  {
	    nh = &nodeHistograms[(iy * (nx-1) + ix) * 256];
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
		      ++nh[image[y*imageWidth+x]];
		  }
	      }
	  }

      /* combine these into a histogram for the entire image */
      images[i].histogram = ih =
	(unsigned int *) malloc(256 * sizeof(unsigned int));
      memset(ih, 0, 256*sizeof(unsigned int));
      for (iy = 0; iy < ny-1; ++iy)
	for (ix = 0; ix < nx-1; ++ix)
	  {
	    nh = &nodeHistograms[(iy * (nx-1) + ix) * 256];
	    for (j = 0; j < 256; ++j)
	      ih[j] += nh[j];
	  }
      
      /* combine with this process's portion of global histogram */
      if (images[i].owner == p)
	for (j = 0; j < 256; ++j)
	  dgh[j] += ih[j];
    }
  if (MPI_Allreduce(dgh, globalHistogram, 256, MPI_DOUBLE,
		    MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("Could not sum up global histogram.\n");

  /* choose global black and white levels */
  globalTotal = 0;
  for (j = 0; j < 256; ++j)
    globalTotal += globalHistogram[j];
  globalBlack = -1;
  globalWhite = -1;
  total = 0;
  for (j = 0; j < 256; ++j)
    {
      total += globalHistogram[j];
      if (globalBlack < 0 && total >= globalTotal / 100)
	globalBlack = j;
      if (globalWhite < 0 && total >= 99 * globalTotal / 100)
	globalWhite = j;
    }

  /* increase coverage by 20% in each direction to ensure we cover
     all potentially valid pixels */
  range = globalWhite - globalBlack;
  globalBlack = globalBlack - range / 5;
  if (globalBlack < 0)
    globalBlack = 0;
  globalWhite = globalWhite + range / 5;
  if (globalWhite > 255)
    globalWhite = 255;
  Log("global black = %d  global white = %d\n",
      globalBlack, globalWhite);

  Log("Before node init barrier\n");
  if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("MPI_Barrier after cpi failed.\n");
  Log("After node init barrier\n");

  /* initialize the image nodes */
  for (i = 0; i < nImages; ++i)
    {
      if (images[i].owner != p && !images[i].needed)
	continue;

      /* choose image-wide black and white levels */
      /* first find the mode */
      ih = images[i].histogram;

      imageTotal = 0;
      for (j = 0; j < 256; ++j)
	imageTotal += ih[j];
      imageBlack = -1;
      imageWhite = -1;
      total = 0;
      for (j = 0; j < 256; ++j)
	{
	  total += ih[j];
	  if (imageBlack < 0 && total > imageTotal / 100)
	    imageBlack = j;
	  if (imageWhite < 0 && total > 99 * imageTotal / 100)
	    imageWhite = j;
	}
      if (imageBlack < globalBlack)
	imageBlack = globalBlack;
      if (imageWhite > globalWhite)
	imageWhite = globalWhite;
      images[i].black = imageBlack;
      images[i].white = imageWhite;
      imageGray = (imageBlack + imageWhite) / 2;
      images[i].gray = imageGray;
      Log("For image %s chose black=%d gray=%d white=%d\n",
	  images[i].name, images[i].black, images[i].gray, images[i].white);

      nx = images[i].nx;
      ny = images[i].ny;
      nodeHistograms = images[i].nodeHistograms;
      images[i].nodes = nodes = (Node*) malloc(2 * ny * nx * sizeof(Node));
      if (images[i].owner != p)
	{
	  memset(nodes, 0, 2 * nx * ny * sizeof(Node));
	  continue;
	}

      Log("Constructing intra-image springs for %s\n", images[i].name);
      nLevelSpringsBlack = 0;
      nLevelSpringsWhite = 0;
      nIntraSpringsBlack = 0;
      nIntraSpringsWhite = 0;
      for (iy = 0; iy < ny; ++iy)
	for (ix = 0; ix < nx; ++ix)
	  {
	    /* set the initial node values to the image-wide values */
	    nodes[iy * nx + ix].x = imageBlack;
	    nodes[ny*nx + (iy * nx + ix)].x = imageWhite;

	    if (nSprings+2 > springsSize)
	      {
		springsSize = (springsSize > 0) ? (springsSize + 1024*1024) : 1024*1024;
		springs = (Spring*) realloc(springs, springsSize * sizeof(Spring));
	      }
	    /* spring for black level */
	    springs[nSprings].node0 = &nodes[iy * nx + ix];
	    springs[nSprings].node1 = NULL;
	    springs[nSprings].k = levelK;
	    springs[nSprings].offset = imageBlack;
	    ++nLevelSpringsBlack;
	    ++nSprings;
	    
	    /* spring for white level */
	    springs[nSprings].node0 = &nodes[ny*nx + (iy * nx + ix)];
	    springs[nSprings].node1 = NULL;
	    springs[nSprings].k = levelK;
	    springs[nSprings].offset = imageWhite;
	    ++nLevelSpringsWhite;
	    ++nSprings;

	    /* make a spring to adjacent nodes */
	    for (dy = -1; dy <= 0; ++dy)
	      for (dx = -1; dx <= 1; ++dx)
		{
		  if (dx == 0 && dy == 0)
		    break;
		  x = ix + dx;
		  y = iy + dy;
		  if (x < 0 || x >= nx || y < 0)
		    continue;
		
		  if (nSprings+2 > springsSize)
		    {
		      springsSize = (springsSize > 0) ? (springsSize + 1024*1024) : 1024*1024;
		      springs = (Spring*) realloc(springs, springsSize * sizeof(Spring));
		    }
		  /* spring for black level */
		  springs[nSprings].node0 = &nodes[y * nx + x];
		  springs[nSprings].node1 = &nodes[iy * nx + ix];
		  springs[nSprings].k = intraimageK;
		  springs[nSprings].offset = 0.0;
		  ++nIntraSpringsBlack;
		  ++nSprings;
		
		  /* spring for white level */
		  springs[nSprings].node0 = &nodes[ny*nx + (y * nx + x)];
		  springs[nSprings].node1 = &nodes[ny*nx + (iy * nx + ix)];
		  springs[nSprings].k = intraimageK;
		  springs[nSprings].offset = 0.0;
		  ++nIntraSpringsWhite;
		  ++nSprings;
		}
	  }
      Log("nLevelSpringsBlack = %d nLevelSpringsWhite = %d\n",
	  nLevelSpringsBlack, nLevelSpringsWhite);
      Log("nIntraSpringsBlack = %d nIntraSpringsWhite = %d\n",
	  nIntraSpringsBlack, nIntraSpringsWhite);
    }

  Log("Initializing maps\n");
  /* initialize the maps that are needed */
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
      iMap = InvertMap(map, mw, mh);
      nInterSpringsBlack = 0;
      nInterSpringsWhite = 0;
      Log("Constructing inter-image springs for map %s\n", maps[i].name);
      
      //		  ++steps[1];
      i0 = maps[i].image0;
      nx0 = images[i0].nx;
      ny0 = images[i0].ny;
      iw0 = images[i0].width;
      ih0 = images[i0].height;
      image0 = images[i0].image;
      mask0 = images[i0].mask;
      mbpl0 = (iw0 + 7) >> 3;
      nodes0 = images[i0].nodes;
      image0Black = images[i0].black;
      image0Gray = images[i0].gray;
      image0White = images[i0].white;

      i1 = maps[i].image1;
      nx1 = images[i1].nx;
      ny1 = images[i1].ny;
      iw1 = images[i1].width;
      ih1 = images[i1].height;
      image1 = images[i1].image;
      mask1 = images[i1].mask;
      mbpl1 = (iw1 + 7) >> 3;
      nodes1 = images[i1].nodes;
      image1Black = images[i1].black;
      image1Gray = images[i1].gray;
      image1White = images[i1].white;

      /* go through all nodes */
      for (iy = 0; iy < ny0; ++iy)
	for (ix = 0; ix < nx0; ++ix)
	  {
	    //			++steps[2];
	    /* find the corresponding point in the other image */
	    x0 = ix * spacing;
	    y0 = iy * spacing;
	    xv = x0 / mFactor;
	    yv = y0 / mFactor;
	    ixv = (int) floor(xv);
	    iyv = (int) floor(yv);
	    rrx = xv - ixv;
	    rry = yv - iyv;
	    ixv -= mxMin;
	    iyv -= myMin;
	    if (ixv < 0 || ixv >= mw-1 ||
		iyv < 0 || iyv >= mh-1)
	      continue;
	    //			++steps[3];
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
	    rx = rx * mFactor / spacing;
	    ry = ry * mFactor / spacing;

	    /* pick the closest grid point */
	    ix1 = floor(rx + 0.5);
	    iy1 = floor(ry + 0.5);
	    if (ix1 < 0 || ix1 >= nx1 ||
		iy1 < 0 || iy1 >= ny1)
	      continue;
	    //			++steps[4];

	    /* go through all pixels in a (spacing x spacing) square centered on
	       (ix,iy) */
	    blackCount0 = 0;
	    whiteCount0 = 0;
	    sumBlack0 = 0.0;
	    sumWhite0 = 0.0;
	    weightBlack0 = 0.0;
	    weightWhite0 = 0.0;
	    for (dy = -spacing/2; dy < spacing/2; ++dy)
	      for (dx = -spacing/2; dx < spacing/2; ++dx)
		{
		  //			      ++steps[5];
		  x = ix * spacing + dx;
		  y = iy * spacing + dy;
		  if (x < 0 || x >= iw0 ||
		      y < 0 || y >= ih0)
		    continue;
		  //			      ++steps[6];
		  /* if pixel is not a valid pixel, ignore */
		  if ((mask0[y*mbpl0+(x>>3)] & (0x80 >> (x & 7))) == 0)
		    continue;
		  //			      ++steps[7];
		  if (image0[y*iw0 + x] < imageBlack ||
		      image0[y*iw0 + x] > imageWhite)
		    continue;
		  //			      ++steps[8];

		  /* find corresponding pixel value in the other image,
		     and compute difference */
		  xv = ((float) x) / mFactor;
		  yv = ((float) y) / mFactor;
		  ixv = (int) floor(xv);
		  iyv = (int) floor(yv);
		  rrx = xv - ixv;
		  rry = yv - iyv;
		  ixv -= mxMin;
		  iyv -= myMin;
		  if (ixv < 0 || ixv >= mw-1 ||
		      iyv < 0 || iyv >= mh-1)
		    continue;
		  //			      ++steps[9];
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
 		  rx = rx * mFactor;
		  ry = ry * mFactor;

		  /* interpolate to find pixel value */
		  irx = (int) floor(rx);
		  iry = (int) floor(ry);
		  if (irx < 0 || irx >= iw1-1 ||
		      iry < 0 || iry >= ih1-1)
		    continue;
		  //			      ++steps[10];
		  if ((mask1[iry*mbpl1 + (irx >> 3)] & (0x80 >> (irx & 7))) == 0 ||
		      (mask1[(iry+1)*mbpl1 + (irx >> 3)] & (0x80 >> (irx & 7))) == 0 ||
		      (mask1[iry*mbpl1 + ((irx+1) >> 3)] & (0x80 >> ((irx+1) & 7))) == 0 ||
		      (mask1[(iry+1)*mbpl1 + ((irx+1) >> 3)] & (0x80 >> ((irx+1) & 7))) == 0)
		    continue;
		  //			      ++steps[11];
		  r00 = image1[iry * iw1 + irx];
		  r01 = image1[(iry+1) * iw1 + irx];
		  r10 = image1[iry * iw1 + irx+1];
		  r11 = image1[(iry+1) * iw1 + irx+1];
		  if (r00 < image1Black || r00 > image1White ||
		      r01 < image1Black || r01 > image1White ||
		      r10 < image1Black || r10 > image1White ||
		      r11 < image1Black || r11 > image1White)
		    continue;
		  //			      ++steps[12];
		  rv = r00 * (rrx - 1.0) * (rry - 1.0)
		    - r10 * rrx * (rry - 1.0) 
		    - r01 * (rrx - 1.0) * rry
		    + r11 * rrx * rry;

		  iv = image0[y*iw0 + x];

		  lvl0 = (iv - image0Black) / (image0White - image0Black);
		  if (lvl0 < 0.0)
		    lvl0 = 0.0;
		  if (lvl0 > 1.0)
		    lvl0 = 1.0;
		  lvl1 = (rv - image1Black) / (image1White - image1Black);
		  if (lvl1 < 0.0)
		    lvl1 = 0.0;
		  if (lvl1 > 1.0)
		    lvl1 = 1.0;
		  lvl = 0.5 * (lvl0 + lvl1);

		  weight = 1.0 - lvl;
		  weightBlack0 += weight;
		  sumBlack0 += weight * (rv - iv);
		  ++blackCount0;

		  weight = lvl;
		  weightWhite0 += weight;
		  sumWhite0 += weight * (rv - iv);
		  ++whiteCount0;
		}

	    /* now go through all pixels in a (spacing x spacing) square centered
	       on the corresponding node */
	    //			++steps[13];
	    blackCount1 = 0;
	    whiteCount1 = 0;
	    sumBlack1 = 0.0;
	    sumWhite1 = 0.0;
	    weightBlack1 = 0.0;
	    weightWhite1 = 0.0;
	    for (dy = -spacing/2; dy < spacing/2; ++dy)
	      for (dx = -spacing/2; dx < spacing/2; ++dx)
		{
		  //			      ++steps[14];
		  x = ix1 * spacing + dx;
		  y = iy1 * spacing + dy;
		  if (x < 0 || x >= iw1 ||
		      y < 0 || y >= ih1)
		    continue;
		  //			      ++steps[15];
		  /* if pixel is not a valid pixel, ignore */
		  if ((mask1[y*mbpl1+(x>>3)] & (0x80 >> (x & 7))) == 0)
		    continue;
		  //			      ++steps[16];

		  /* Find corresponding pixel value in the first image,
		     and compute difference */
		  xv = ((float) x) / mFactor;
		  yv = ((float) y) / mFactor;
		  if (!Invert(iMap, &rx, &ry, xv, yv))
		    continue;
		  //			      ++steps[17];
		  rx = (rx + mxMin) * mFactor;
		  ry = (ry + myMin) * mFactor;

		  /* interpolate to find pixel value */
		  irx = (int) floor(rx);
		  iry = (int) floor(ry);
		  if (irx < 0 || irx >= iw0-1 ||
		      iry < 0 || iry >= ih0-1)
		    continue;
		  //			      ++steps[18];
		  if ((mask0[iry*mbpl0 + (irx >> 3)] & (0x80 >> (irx & 7))) == 0 ||
		      (mask0[(iry+1)*mbpl0 + (irx >> 3)] & (0x80 >> (irx & 7))) == 0 ||
		      (mask0[iry*mbpl0 + ((irx+1) >> 3)] & (0x80 >> ((irx+1) & 7))) == 0 ||
		      (mask0[(iry+1)*mbpl0 + ((irx+1) >> 3)] & (0x80 >> ((irx+1) & 7))) == 0)
		    continue;
		  //			      ++steps[19];
		  r00 = image0[iry * iw0 + irx];
		  r01 = image0[(iry+1) * iw0 + irx];
		  r10 = image0[iry * iw0 + irx+1];
		  r11 = image0[(iry+1) * iw0 + irx+1];
		  rv = r00 * (rrx - 1.0) * (rry - 1.0)
		    - r10 * rrx * (rry - 1.0) 
		    - r01 * (rrx - 1.0) * rry
		    + r11 * rrx * rry;
		  
		  iv = image1[y*iw1 + x];

		  lvl0 = (rv - image0Black) / (image0White - image0Black);
		  if (lvl0 < 0.0)
		    lvl0 = 0.0;
		  if (lvl0 > 1.0)
		    lvl0 = 1.0;
		  lvl1 = (iv - image1Black) / (image1White - image1Black);
		  if (lvl1 < 0.0)
		    lvl1 = 0.0;
		  if (lvl1 > 1.0)
		    lvl1 = 1.0;
		  lvl = 0.5 * (lvl0 + lvl1);

		  weight = 1.0 - lvl;
		  weightBlack1 += weight;
		  sumBlack1 += weight * (iv - rv);
		  ++blackCount1;

		  weight = lvl;
		  weightWhite1 += weight;
		  sumWhite1 += weight * (iv - rv);
		  ++whiteCount1;
		}		    

	    if (weightBlack0 > 0.0 || weightBlack1 > 0.0)
	      {
		diff0to1 = (sumBlack0 + sumBlack1) /
		  (weightBlack0 + weightBlack1);

		/* construct a spring representing the average difference of black values */
		if (nSprings+1 > springsSize)
		  {
		    springsSize = (springsSize > 0) ? (springsSize + 1024*1024) : 1024*1024;
		    springs = (Spring*) realloc(springs, springsSize * sizeof(Spring));
		  }
#if 0
		if (i0 == 2 && ix == 41 && iy == 34)
		  Log("Creating blk SPRING %d to %d %d %d   %f\n", nSprings, i1, ix1, iy1, diff0to1);
		if (i1 == 2 && ix1 == 41 && iy1 == 34)
		  Log("Creating blk SPRING %d from %d %d %d  %f\n", nSprings, i0, ix, iy, diff0to1);
#endif
		springs[nSprings].node0 = &nodes0[iy * nx0 + ix];
		springs[nSprings].node1 = &nodes1[iy1 * nx1 + ix1];
		springs[nSprings].k = interimageK * (weightBlack0 + weightBlack1) / (blackCount0 + blackCount1);
		springs[nSprings].offset = diff0to1;
		++nInterSpringsBlack;
		++nSprings;
		//			    ++steps[21];
	      }

	    if (weightWhite0 > 0.0 || weightWhite1 > 0.0)
	      {
		diff0to1 = (sumWhite0 + sumWhite1) /
		  (weightWhite0 + weightWhite1);

		/* construct a spring representing the average difference of black values */
		if (nSprings+1 > springsSize)
		  {
		    springsSize = (springsSize > 0) ? (springsSize + 1024*1024) : 1024*1024;
		    springs = (Spring*) realloc(springs, springsSize * sizeof(Spring));
		  }
#if 0
		if (i0 == 2 && ix == 41 && iy == 34)
		  Log("Creating wht SPRING %d to %d %d %d  %f\n", nSprings, i1, ix1, iy1, diff0to1);
		if (i1 == 2 && ix1 == 41 && iy1 == 34)
		  Log("Creating wht SPRING %d from %d %d %d  %f\n", nSprings, i0, ix, iy, diff0to1);
#endif
		springs[nSprings].node0 = &nodes0[ny0*nx0 + (iy * nx0 + ix)];
		springs[nSprings].node1 = &nodes1[ny1*nx1 + (iy1 * nx1 + ix1)];
		springs[nSprings].k = interimageK * (weightWhite0 + weightWhite1) / (whiteCount0 + whiteCount1);
		springs[nSprings].offset = diff0to1;
		++nInterSpringsWhite;
		++nSprings;
		//			    ++steps[21];
	      }
	  }
      Log("nInterSpringsBlack = %d nInterSpringsWhite = %d\n",
	  nInterSpringsBlack, nInterSpringsWhite);

      FreeInverseMap(iMap);
      free(map);
    }

  /* set up communications with other nodes */
  PlanCommunications();

  /* relax the spring system */
  Log("Starting relaxation iterations.\n");
  dampingFactor = 0.1;
  epochIterations = 512;
  terminationRequestedIter = -1;
  nIncrease = 0;
  nDecrease = 0;
  epochInitialTotalEnergy = 1.0e+30;
  epochFinalTotalEnergy = 1.0e+30;
  prevTotalEnergy = 1.0e+30;
  nDecrease = 0;
  nIncrease = 0;
  for (iter = 0; ; ++iter)
    {
#if DEBUG
      Log("Starting iteration %d\n", iter);
#endif
      CommunicateIntensities();

      /* update all forces, also computing energy */
      energy = 0.0;
      basis = energy;
      for (i = myFirstImage; i <= myLastImage; ++i)
	{
	  nx = images[i].nx;
	  ny = images[i].ny;
	  nodes = images[i].nodes;

	  /* zero all forces */
	  for (y = 0; y < ny; ++y)
	    for (x = 0; x < nx; ++x)
	      {
		if (nodes[y * nx + x].x < -512.0 ||
		    nodes[y * nx + x].x > 512.0 ||
		    nodes[ny*nx + (y * nx + x)].x < -512.0 ||
		    nodes[ny*nx + (y * nx + x)].x > 512.0)
		  {
		    Log("out-of-bounds value: %s %d %d %f %f\n",
			images[i].name, x, y,
			nodes[y * nx + x].x,
			nodes[ny*nx + (y * nx + x)].x);
		    exit(1);
		  }
		nodes[y * nx + x].fx = 0.0;
		nodes[ny*nx + (y * nx + x)].fx = 0.0;
	      }

	}

      //      Log("FX 0: %f\n", images[2].nodes[98*74+34*98+41].fx);

      /* update all forces, also computing energy */
      for (i = 0; i < nSprings; ++i)
	{
	  s = &springs[i];
	  if (s->node1 == NULL)
	    {
	      if (isnan(s->offset) || isnan(s->node0->x) || isnan(s->k))
		{
		  Log("invalid value: %f %f %f %d\n",
		      s->offset, s->node0->x, s->k, i);
		  exit(1);
		}
	      delta = s->offset - s->node0->x;
	      force = s->k * delta;
	      s->node0->fx += force;
	      energy += force * delta;
	    }
	  else
	    {
	      if (isnan(s->offset) || isnan(s->node0->x) || isnan(s->node1->x) || isnan(s->k))
		{
		  Log("invalid value: %f %f %f %f %d\n",
		      s->offset, s->node0->x, s->node1->x, s->k, i);
		  exit(1);
		}
	      delta = s->node1->x - s->node0->x - s->offset;
	      force = s->k * delta;
	      s->node0->fx += force;
	      s->node1->fx -= force;
	      energy += force * delta;
#if 0
	      if ((i & ~1) == 234452)
		Log("iter %d spr %d: %f %f %f %f  %f %f  %f %f\n",
		    iter, i,
		    s->offset, s->node0->x, s->node1->x, s->k,
		    delta, force,
		    s->node0->fx,
		    s->node1->fx);
#endif
	    }
	}

      //      Log("FX 1: %f\n", images[2].nodes[98*74+34*98+41].fx);

      maxF = 0.0;
      for (i = myFirstImage; i <= myLastImage; ++i)
	{
	  nx = images[i].nx;
	  ny = images[i].ny;
	  nodes = images[i].nodes;
	  for (iy = 0; iy < ny; ++iy)
	    for (ix = 0; ix < nx; ++ix)
	      {
#if 0
		if (i == 2 && ix == 41 && iy == 34)
		  {
		    Log("iter %d PT x = %f fx = %f nx = %d ny = %d\n",
			iter,
			nodes[iy*nx+ix].x,
			nodes[iy*nx+ix].fx,
			nx, ny);
		  }
#endif
		if (nodes[iy * nx + ix].fx > maxF)
		  {
		    maxF = nodes[iy * nx + ix].fx;
		    ni = i;
		    nix = ix;
		    niy = iy;
		    nlv = 0;
		  }
		if (nodes[ny*nx + (iy * nx + ix)].fx > maxF)
		  {
		    maxF = nodes[ny*nx + (iy * nx + ix)].fx;
		    ni = i;
		    nix = ix;
		    niy = iy;
		    nlv = 1;
		  }
	      }
	}
      if (iter % 100 == 0)
	Log("Iter %d: Node maxF occurred at %s %d %d %d: %f\n",
	    iter,
	    images[ni].name,
	    nix, niy, nlv, maxF);
      //      Log("FX 2: %f\n", images[2].nodes[98*74+34*98+41].fx);

      /* find global maximum */
      if (MPI_Allreduce(&maxF, &globalMaxF, 1, MPI_FLOAT, MPI_MAX,
			MPI_COMM_WORLD) != MPI_SUCCESS)
	Error("Could not find global maximum force\n");

      /* update all positions */
      if (globalMaxF > 0.5)
	scale = dampingFactor * 0.5 / globalMaxF;
      else
	scale = dampingFactor;
      maxStep = 0.1;
      for (i = myFirstImage; i <= myLastImage; ++i)
	{
	  nx = images[i].nx;
	  ny = images[i].ny;
	  nodes = images[i].nodes;
	  for (iy = 0; iy < ny; ++iy)
	    for (ix = 0; ix < nx; ++ix)
	      {
		delta = scale * nodes[iy * nx + ix].fx;
		if (delta > maxStep)
		  delta = maxStep;
		nodes[iy * nx + ix].x += delta;
#if 0
		if (i == 2 && ix == 41 && iy == 34)
		  Log("iter %d ADJ blk %f %f %f %f %f\n",
		      iter, scale, nodes[iy*nx+ix].fx, delta,
		      maxStep, nodes[iy*nx+ix].x);
#endif
		delta = scale * nodes[ny*nx + (iy * nx + ix)].fx;
		if (delta > maxStep)
		  delta = maxStep;
		nodes[ny*nx + (iy * nx + ix)].x += delta;
#if 0
		if (i == 2 && ix == 41 && iy == 34)
		  Log("iter %d ADJ wht %f %f %f %f %f\n",
		      iter, scale, nodes[ny*nx + iy*nx+ix].fx, delta,
		      maxStep, nodes[ny*nx + iy*nx+ix].x);
#endif
	      }
	}

      //      Log("FX 3: %f\n", images[2].nodes[98*74+34*98+41].fx);

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
	      if (epochInitialTotalEnergy - totalEnergy <
		  totalEnergy * threshold * 0.000001 * epochIterations ||
		  epochFinalTotalEnergy - totalEnergy <
		  totalEnergy * threshold * 0.000001 * epochIterations ||
		  totalEnergy < 0.000001)
		terminationRequestedIter = iter;
	    }
	  epochFinalTotalEnergy = totalEnergy;
	  if (MPI_Bcast(&terminationRequestedIter, 1, MPI_INT, 0, MPI_COMM_WORLD) !=
	      MPI_SUCCESS)
	    Error("Could not broadcast control flag.\n");
	}

      if (deltaEnergy > 0.0)
	{
	  nDecrease = 0;
	  ++nIncrease;
	  if (nIncrease > 1000 ||
	      nIncrease > 10 && deltaEnergy > totalEnergy ||
	      deltaEnergy > 1000.0 * totalEnergy)
	    Error("Alignment is diverging instead of converging.\n");
	      
	  if (iter < 100 || nIncrease < 2)
	    {
	      dampingFactor *= 0.5;
	      if (dampingFactor < 0.000001)
		{
		  Log("Damping factor = %f... ending iterations.\n",
		      dampingFactor);
		  break;
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

      /* 32 is chosen as a preferred time to output or terminate because,
	 given that we increase the damping by 1% on each iteration, and decrease
	 it by a factor of 2 when instability sets in, the number of iterations
	 between instabilities is about 70, and we want to output/terminate
	 when the state is not near an instability */
      if (terminationRequestedIter >= 0 &&
	  (nDecrease == 32 || iter > terminationRequestedIter + 128))
	break;
    }
  if (p == 0)
    Log("Finished intensity adjustment at iteration %d.\n", iter);
      
  /* output the intensity maps */
  for (i = myFirstImage; i <= myLastImage; ++i)
    {
      nx = images[i].nx;
      ny = images[i].ny;
      nodes = images[i].nodes;
      map = (MapElement*) malloc(nx * ny * sizeof(MapElement));
      for (iy = 0; iy < ny; ++iy)
	for (ix = 0; ix < nx; ++ix)
	  {
	    map[iy*nx+ix].x = nodes[iy*nx + ix].x;
	    map[iy*nx+ix].y = nodes[ny*nx + iy*nx + ix].x;
	    if (map[iy*nx+ix].x >= map[iy*nx+ix].y)
	      map[iy*nx+ix].y = map[iy*nx+ix].x + 1.0;
	    map[iy*nx+ix].c = 1.0;
	  }

      /* write out the map */
      sprintf(fn, "%s%s.map", outputName, images[i].name);
      if (!CreateDirectories(fn))
	{
	  fprintf(stderr, "Could not create directories for %s\n", fn);
	  exit(1);
	}
      sprintf(imgName, "%s%s", imagesName, images[i].name);
      if (!WriteMap(fn, map, level, nx, ny, 0, 0, imgName, imgName,
		    UncompressedMap, msg))
	{
	  fprintf(stderr, "Could not write intensity map %s:\n%s\n",
		  fn, msg);
	  exit(1);
	}
    }

  Log("FINALIZING\n");
  MPI_Finalize();
  fclose(logFile);
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

void
PlanCommunications ()
{
  int phase;
  CommPhase *cp;
  int bufferPos;
  int i, j, k;
  int nFloats;

  Log("Planning communication\n");
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
	  nFloats = 2 * images[j].nx * images[j].ny;
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
	  nFloats = 2 * images[j].nx * images[j].ny;
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
CommunicateIntensities ()
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
		      nNodes = 2 * images[its].nx * images[its].ny;
		      nodes = images[its].nodes;
		      for (k = 0; k < nNodes; ++k)
			buffer[bufferPos++] = nodes[k].x;
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
		      nNodes = 2 * images[its].nx * images[its].ny;
		      nodes = images[its].nodes;
		      for (k = 0; k < nNodes; ++k)
			nodes[k].x =  buffer[bufferPos++];
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
      
      if (mkdir(dn, 0777) != 0 && errno != EEXIST)
	{
	  Log("Could not create directory %s\n", dn);
	  return(0);
	}
    }
  return(1);
}
