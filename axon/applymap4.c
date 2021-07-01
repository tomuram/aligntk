/*
 *  applymap4.c  -  applies maps to images
 *
 *  Copyright (c) 2009-2011 National Resource for Biomedical
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
 *    2009  Written by Greg Hood (ghood@psc.edu)
 *    2011  Fixed out-of-range array indices when rendering multi-gigabyte
 *             images (ghood@psc.edu)
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <dirent.h>
#include <stdarg.h>
#include <errno.h>
#include <sched.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "imio.h"
#include "invert.h"
#include "dt.h"

#define LINE_LENGTH		255
#define MAX_LABEL_LENGTH	255
#define FONT_FILE		"/usr/users/6/ghood/warp/LTYPEB.pgm"

/* TYPES */
typedef struct Image
{
  int next;		/* index of next image in hash bucket */
  char *name;   	/* name of this image */
  int width, height;    /* width and height in pixels */
  unsigned char *image; /* the image bytes (NULL if not loaded) */
  unsigned char *mask;  /* the image mask, if present */
  unsigned char *dist;  /* the distance array -- each element holds the distance
			   from the corresponding image pixel to the closest
			   masked pixel (in units of 1/64 pixels) */
  float minX, maxX;     /* bounds of the area this image covers in the */
  float minY, maxY;	/*   final image */
  int needed;		/* true if this image is needed for the current set of
			   tiles */
  int mapBytes;         /* size of maps in bytes */
  int mLevel;		/* map level */
  int mw, mh;		/* map width and height */
  int mxMin, myMin;	/* map offset in x and y */
  MapElement *map;      /* map of this image into the final image */
  InverseMap *invMap;     /* inverse map that translates points in the final image
			   into points in this image */
  /* intensity map is optional */
  int imapLevel;	/* intensity map level */
  int imapw, imaph;	/* intensity map width and height */
  MapElement *imap;	/* intensity map */
  MapElement *targetMap;/* map of where pixels from this image ended up
			   in target image */
  time_t mtime;         /* the modification time for this image or its map */
} Image;

/* GLOBAL VARIABLES */
int resume = 0;
char imageListName[PATH_MAX];
char imageName[PATH_MAX];
char imagesName[PATH_MAX];
char extension[PATH_MAX];
char masksName[PATH_MAX];
char mapsName[PATH_MAX];
char imapsName[PATH_MAX];
char outputName[PATH_MAX];
char sourceMapName[PATH_MAX];
char targetMapsName[PATH_MAX];
int overlay = 0;
int blend = 1;
int margin = -1;
int tileWidth = -1;
int tileHeight = -1;
int memoryLimit = 1024;   /* 1 GB */
int tree = 0;
float blackValue = 0.0;
float whiteValue = 255.0;
float mapScale = 1.0;
float imapScale = 1.0;
float maskScale = 1.0;
int compress = 0;
int reductionFactor = 1;
int sourceMapLevel = 6;
int targetMapsLevel = 6;
int update = 0;
float rotation = 0.0;
float rotationX = 0.0;
float rotationY = 0.0;

int nImages = 0;
Image *images = 0;
int *imageHashTable = 0;
unsigned char *canvas = 0;
unsigned short *weight = 0;
size_t canvasSize = 0;
size_t canvasWidth, canvasHeight;
size_t oldCanvasWidth = 0;
int canvasMinX, canvasMinY;
size_t weightWidth, weightHeight;
int weightMinX, weightMinY;
size_t tw, th;
int oMinX, oMinY, oMaxX, oMaxY;
size_t oWidth, oHeight;
float range;
unsigned char *out;
size_t outWidth, outHeight;
int outMinY;
int nProcessed = 0;
size_t imageMem = 0;
int sourceMapSize;
MapElement *sourceMap = 0;
size_t sourceMapWidth, sourceMapHeight;
int sourceMapFactor;
int sourceMapMask;
int targetMapsFactor;
char label[MAX_LABEL_LENGTH+1];
char labelName[MAX_LABEL_LENGTH+1];
unsigned char *font = 0;
int fontWidth, fontHeight;
float cosRot, sinRot;

#define DIR_HASH_SIZE	8192
char *dirHash[DIR_HASH_SIZE];

/* FORWARD DECLARATIONS */
void PaintImage (int i, int minX, int maxX, int minY, int maxY);
void WriteTiles (int col, int startRow, int endRow, char *iName);
void Error (char *fmt, ...);
unsigned int Hash (char *s);
int CreateDirectories (char *fn);
void PrintUsage ();

int
main (int argc, char **argv, char **envp)
{
  int i, j;
  int n;
  int error;
  char fn[PATH_MAX];
  MapElement *map;
  int x, y;
  float minX, minY, maxX, maxY;
  int iMinX, iMinY, iMaxX, iMaxY;
  float xv, yv;
  float spacing;
  float ax, bx, ay, by;
  int mLevel;
  int mw, mh;
  int mxMin, myMin;
  char imName0[PATH_MAX];
  char imName1[PATH_MAX];
  char msg[PATH_MAX+256];
  FILE *f;
  int imagesSize;
  char line[LINE_LENGTH+1];
  int dir;
  int nItems;
  int width, height;
  int rows, cols;
  int ixv, iyv;
  float rx, ry;
  float rrx, rry;
  float rx00, rx01, rx10, rx11, ry00, ry01, ry10, ry11;
  float rv;
  cpu_set_t cpumask;
  int nOutputImages;
  int startImage, endImage;
  int imapLevel;
  int imapw, imaph;
  int imapXMin, imapYMin;
  char imapName0[PATH_MAX], imapName1[PATH_MAX];
  MapElement *imap;
  int regionWidth, regionHeight, regionOffsetX, regionOffsetY;
  int labelWidth, labelHeight, labelOffsetX, labelOffsetY;
  size_t *increase, *decrease;
  size_t memoryRequired;
  size_t maxMemoryRequired;
  int startX, endX;
  int startY, endY;
  int tx, ty;
  int tCol;
  int ix, iy;
  int dx, dy;
  int cx;
  int nx, ny;
  int oi;
  int hs;
  int hi;
  int startRow, endRow;
  int iv;
  int sum;
  int offset;
  int a;
  int valid;
  struct stat sb;
  int updateOutput;
  int labelMinX, labelMaxX, labelMinY, labelMaxY;
  int charMinX, charMaxX;
  int ci;
  float sx, ex, sy, ey;
  int isx, iex, isy, iey;
  float v;
  float wx, wy, ws;
  int len;
  int ilv;
  float rxp, ryp;

  error = 0;
  imageListName[0] = '\0';
  imageName[0] = '\0';
  imagesName[0] = '\0';
  masksName[0] = '\0';
  mapsName[0] = '\0';
  imapsName[0] = '\0';
  outputName[0] = '\0';
  sourceMapName[0] = '\0';
  targetMapsName[0] = '\0';
  labelName[0] = '\0';
  strcpy(extension, "tif");
  regionWidth = regionHeight = regionOffsetX = regionOffsetY = -1;
  labelWidth = labelHeight = labelOffsetX = labelOffsetY = -1;
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
    else if (strcmp(argv[i], "-image") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(imageName, argv[i]);
      }
    else if (strcmp(argv[i], "-images") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(imagesName, argv[i]);
      }
    else if (strcmp(argv[i], "-pgm") == 0)
      strcpy(extension, "pgm");
    else if (strcmp(argv[i], "-masks") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(masksName, argv[i]);
      }
    else if (strcmp(argv[i], "-mask_scale") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%f", &maskScale) != 1)
	  {
	    error = 1;
	    break;
	  }
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
    else if (strcmp(argv[i], "-map_scale") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%f", &mapScale) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-imaps") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(imapsName, argv[i]);
      }
    else if (strcmp(argv[i], "-imap_scale") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%f", &imapScale) != 1)
	  {
	    error = 1;
	    break;
	  }
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
    else if (strcmp(argv[i], "-overlay") == 0)
      overlay = 1;
    else if (strcmp(argv[i], "-blend") == 0)
      blend = 1;
    else if (strcmp(argv[i], "-mosaic") == 0)
      blend = 0;
    else if (strcmp(argv[i], "-margin") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d", &margin) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-tile") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%dx%d",
				  &tileWidth, &tileHeight) != 2)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-memory") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d", &memoryLimit) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-tree") == 0)
      tree = 1;
    else if (strcmp(argv[i], "-black") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%f", &blackValue) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-white") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%f", &whiteValue) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-compress") == 0)
      compress = 1;
    else if (strcmp(argv[i], "-resume") == 0)
      resume = 1;
    else if (strcmp(argv[i], "-region") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%dx%d%d%d",
				  &regionWidth, &regionHeight,
				  &regionOffsetX, &regionOffsetY) != 4)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-reduction") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d", &reductionFactor) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-source_map") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(sourceMapName, argv[i]);
      }
    else if (strcmp(argv[i], "-source_map_level") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d", &sourceMapLevel) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-target_maps") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(targetMapsName, argv[i]);
      }
    else if (strcmp(argv[i], "-target_maps_level") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d", &targetMapsLevel) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-update") == 0)
      update = 1;
    else if (strcmp(argv[i], "-label") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(labelName, argv[i]);
      }
    else if (strcmp(argv[i], "-label_location") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%dx%d%d%d",
				  &labelWidth, &labelHeight,
				  &labelOffsetX, &labelOffsetY) != 4)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-rotation") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%f", &rotation) != 1)
	  {
	    error = 1;
	    break;
	  }
	cosRot = cos(rotation * M_PI / 180.0);
	sinRot = sin(rotation * M_PI / 180.0);
      }
    else if (strcmp(argv[i], "-rotation_center") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%f,%f",
				  &rotationX, &rotationY) != 2)
	  {
	    error = 1;
	    break;
	  }
      }
    else error = 1;

  if (error)
    {
      fprintf(stderr, "Usage: applymap4 -image_list list_file -images image_prefix\n");
      fprintf(stderr, "              -maps map_prefix -output file_prefix\n");
      fprintf(stderr, "              [-map_scale scaling_factor]\n");
      fprintf(stderr, "              [-masks mask_prefix]\n");
      fprintf(stderr, "              [-masks_scale scaling_factor]\n");
      fprintf(stderr, "              [-imaps imaps_prefix]\n");
      fprintf(stderr, "              [-imap_scale scaling_factor]\n");
      fprintf(stderr, "              [-blend]\n");
      fprintf(stderr, "              [-mosaic]\n");
      fprintf(stderr, "              [-margin margin_thickness_in_pixels]\n");
      fprintf(stderr, "              [-tile tilewidthxtileheight]\n");
      fprintf(stderr, "              [-memory memory_limit_in_MB]\n");
      fprintf(stderr, "              [-tree]\n");
      fprintf(stderr, "              [-black black_value]\n");
      fprintf(stderr, "              [-white white_value]\n");
      fprintf(stderr, "              [-compress]\n");
      fprintf(stderr, "              [-region WxH+X+Y]\n");
      fprintf(stderr, "              [-reduction reduction_factor]\n");
      fprintf(stderr, "              [-resume]\n");
      fprintf(stderr, "              [-rotation CCW_rotation_in_degrees]\n");
      fprintf(stderr, "              [-rotation_center x,y]\n");
      fprintf(stderr, "              [-source_map map_name]\n");
      fprintf(stderr, "              [-source_map_level level]\n");
      fprintf(stderr, "              [-target_maps map_name_prefix]\n");
      fprintf(stderr, "              [-target_maps_level level]\n");
      fprintf(stderr, "              [-update]\n");
      fprintf(stderr, "              [-label WxH+Y+Y]\n");
      exit(1);
    }

  /* check that at least minimal parameters were supplied */
  if ((imageListName[0] == '\0' && imageName[0] == '\0') ||
      imagesName[0] == '\0' ||
      mapsName[0] == '\0' ||
      outputName[0] == '\0')
    {
      fprintf(stderr, "-image_list (or -image), -images, -maps, and -output parameters must be specified.\n");
      exit(1);
    }
  range = whiteValue - blackValue;
  if (range == 0.0)
    Error("White value cannot be same as black value\n");

  /* load the font (if necessary) */
  if (labelWidth > 0)
    {
      if (!ReadImage(FONT_FILE, &font, &fontWidth, &fontHeight, -1, -1, -1, -1, msg))
	Error("Could not load font: %s\n", msg);
      fontHeight = fontHeight / 95;
    }

  if (imageName[0] != '\0')
    {
      images = (Image *) malloc(sizeof(Image));
      images[0].name = (char *) malloc(strlen(imageName) + 1);
      strcpy(images[0].name, imageName);
      images[0].width = -1;
      images[0].height = -1;
      images[0].image = NULL;
      images[0].mask = NULL;
      images[0].dist = NULL;
      images[0].map = NULL;
      images[0].invMap = NULL;
      images[0].imap = NULL;
      images[0].targetMap = NULL;
      nImages = 1;
      imagesSize = 1;
    }
  else
    {
      /* read the image list */
      nImages = 0;
      imagesSize = 0;
      f = fopen(imageListName, "r");
      if (f == NULL)
	Error("Could not open file %s for reading\n", imageListName);
      while (fgets(line, LINE_LENGTH, f) != NULL)
	{
	  if (line[0] == '\0' || line[0] == '#')
	    continue;
	  width = -1;
	  height = -1;
	  nItems = sscanf(line, "%s%d%d", imageName, &width, &height);
	  if (nItems != 1 && nItems != 3)
	    Error("Malformed line in %s:\n%s\n", imageListName, line);

	  if (nImages >= imagesSize)
	    {
	      imagesSize = (imagesSize > 0) ? imagesSize * 2 : 64;
	      images = (Image*) realloc(images, imagesSize * sizeof(Image));
	    }
	  images[nImages].name = (char*) malloc(strlen(imageName) + 1);
	  strcpy(images[nImages].name, imageName);
	  images[nImages].width = width;
	  images[nImages].height = height;
	  images[nImages].image = NULL;
	  images[nImages].mask = NULL;
	  images[nImages].dist = NULL;
	  images[nImages].map = NULL;
	  images[nImages].invMap = NULL;
	  images[nImages].imap = NULL;
	  images[nImages].targetMap = NULL;
	  ++nImages;
	}
      fclose(f);
      images = (Image *) realloc(images, nImages * sizeof(Image));
      printf("nImages = %d\n", nImages);
    }
  for (i = 0; i < nImages; ++i)
    {
      /* check that the image exists and get its modification time */
      /* FIX:  TEMPORARY HACK */
      sprintf(fn, "%s%s.%s", imagesName, images[i].name, extension);
      if (stat(fn, &sb) != 0)
	{
	  sprintf(fn, "%s%s.tif", imagesName, images[i].name);
	  if (stat(fn, &sb) != 0)
	    Error("Could not stat file %s\n", fn);
	}
      if (S_ISDIR(sb.st_mode))
	Error("Image %s is a directory.\n", fn);
      images[i].mtime = sb.st_mtime;

      /* read the image size */
      if (images[i].width < 0 || images[i].height < 0)
	{
	  if (!ReadImageSize(fn, &width, &height, msg))
	    Error("Could not determine image size of %s:\n%s\n", fn, msg);
	  images[i].width = width;
	  images[i].height = height;
	}
    }

  printf("Previewing maps: ");
  fflush(stdout);
  oMinX = 1000000000;
  oMaxX = -1000000000;
  oMinY = 1000000000;
  oMaxY = -1000000000;
  map = NULL;
  imap = NULL;
  sourceMapFactor = 1 << sourceMapLevel;
  sourceMapMask = sourceMapFactor - 1;
  targetMapsFactor = 1 << targetMapsLevel;
  for (i = 0; i < nImages; ++i)
    {
      sprintf(fn, "%s%s.map", mapsName, images[i].name);
      if (stat(fn, &sb) != 0)
	Error("Could not stat map %s\n", fn);
      if (sb.st_mtime > images[i].mtime)
	images[i].mtime = sb.st_mtime;
      if (!ReadMap(fn, &map, &mLevel,
		   &mw, &mh, &mxMin, &myMin,
		   imName0, imName1,
		   msg))
	Error("Could not read map %s:\n  error: %s\n",
	      fn, msg);

      spacing = (1 << mLevel) * mapScale;
      minX = 1000000000;
      maxX = -1000000000;
      minY = 1000000000;
      maxY = -1000000000;
      for (y = 0; y < mh-1; ++y)
	for (x = 0; x < mw-1; ++x)
	  {
	    if (map[y*mw+x].c == 0.0 ||
		map[y*mw+x+1].c == 0.0 ||
		map[(y+1)*mw+x].c == 0.0 ||
		map[(y+1)*mw+x+1].c == 0.0)
	      continue;
	    for (dy = 0; dy < 2; ++dy)
	      for (dx = 0; dx < 2; ++dx)
		{
		  rx = map[(y+dy)*mw+x+dx].x * spacing;
		  ry = map[(y+dy)*mw+x+dx].y * spacing;
		  if (rotation != 0.0)
		    {
		      rxp = rx - rotationX;
		      ryp = ry - rotationY;
		      rx = cosRot * rxp + sinRot * ryp + rotationX;
		      ry = -sinRot * rxp + cosRot * ryp + rotationY;
		    }
		  if (rx < minX)
		    minX = rx;
		  if (rx > maxX)
		    maxX = rx;
		  if (ry < minY)
		    minY = ry;
		  if (ry > maxY)
		    maxY = ry;
		}
	  }
      free(map);
      images[i].mapBytes = mw * mh * (sizeof(MapElement) + sizeof(InverseMapElement) +
				      sizeof(unsigned char)) +
	5 * 5 * sizeof(long long) +
	sizeof(InverseMap) +
	256;

      printf("IMAGE BOUNDS[%d] = %f %f %f %f\n", i, minX, maxX, minY, maxY);
      images[i].minX = minX;
      images[i].maxX = maxX;
      images[i].minY = minY;
      images[i].maxY = maxY;
      
      if (minX < oMinX)
	oMinX = (int) floor(minX);
      if (maxX > oMaxX)
	oMaxX = (int) ceil(maxX);
      if (minY < oMinY)
	oMinY = (int) floor(minY);
      if (maxY > oMaxY)
	oMaxY = (int) ceil(maxY);

      if (imapsName[0] != '\0')
	{
	  sprintf(fn, "%s%s.map", imapsName, images[i].name);
	  if (!ReadMap(fn, &imap, &imapLevel,
		       &imapw, &imaph, &imapXMin, &imapYMin,
		       imapName0, imapName1,
		       msg))
	    Error("Could not read map %s:\n  error: %s\n",
		  fn, msg);
	  images[i].mapBytes += imapw * imaph * sizeof(MapElement);
	  free(imap);
	}

      if (targetMapsName[0] != '\0')
	images[i].mapBytes +=
	  ((images[i].width + targetMapsFactor - 1) / targetMapsFactor) *
	  ((images[i].height + targetMapsFactor - 1) / targetMapsFactor) *
	  sizeof(MapElement);

      if ((nProcessed % 50) == 0 && nProcessed != 0)
	printf(" %d\n    ", nProcessed);
      printf(".");
      fflush(stdout);
      ++nProcessed;
    }

  printf("\nAll images previewed.\n");
  oWidth = oMaxX - oMinX + 1;
  oHeight = oMaxY - oMinY + 1;
  printf("output width = %zu (%d to %d) output height = %zu (%d to %d)\n\n",
	 oWidth, oMinX, oMaxX,
	 oHeight, oMinY, oMaxY);

  if (regionWidth > 0)
    {
      oWidth = regionWidth;
      oHeight = regionHeight;
      oMinX = regionOffsetX;
      oMinY = regionOffsetY;
    }

  oWidth = ((oWidth + reductionFactor - 1) / reductionFactor) *
    reductionFactor;
  oHeight = ((oHeight + reductionFactor - 1) / reductionFactor) *
    reductionFactor;
  oMaxX = oMinX + ((int) oWidth) - 1;
  oMaxY = oMinY + ((int) oHeight) - 1;

  if (overlay)
    nOutputImages = 1;
  else
    nOutputImages = nImages;

  if (tileWidth > 0)
    {
      cols = (oWidth / reductionFactor + tileWidth - 1) / tileWidth;
      tw = tileWidth;
    }
  else
    {
      cols = 1;
      tw = oWidth / reductionFactor;
    }
  if (tileHeight > 0)
    {
      rows = (oHeight / reductionFactor + tileHeight - 1) / tileHeight;
      th = tileHeight;
    }
  else
    {
      rows = 1;
      th = oHeight / reductionFactor;
    }
  if (tileWidth > 0 || tileHeight > 0)
    printf("Tiling output into %d rows and %d columns\n",
	   rows, cols);

  increase = (size_t*) malloc(oWidth * sizeof(size_t));
  decrease = (size_t*) malloc(oWidth * sizeof(size_t));
  if (increase == 0 || decrease == 0)
    Error("malloc of increase/decrease failed; errno = %d\n", errno);

  for (oi = 0; oi < nOutputImages; ++oi)
    {
      if (overlay)
	{
	  startImage = 0;
	  endImage = nImages-1;
	  if (labelWidth > 0)
	    strncpy(label, (labelName[0] != '\0') ? labelName : outputName,
		    MAX_LABEL_LENGTH);
	}
      else
	{
	  startImage = oi;
	  endImage = oi;
	  if (labelWidth > 0)
	    snprintf(label, MAX_LABEL_LENGTH, "%s%s",
		     labelName, images[oi].name);
	}
      if (labelWidth > 0)
	{
	  label[MAX_LABEL_LENGTH] = '\0';
	  labelMinX = labelOffsetX;
	  labelMaxX = labelMinX + strlen(label) * labelWidth - 1;
	  labelMinY = labelOffsetY;
	  labelMaxY = labelMinY + labelHeight - 1;
	}

      updateOutput = update;
      if (sourceMapName[0] != '\0')
	{
	  sourceMap = NULL;
	  sourceMapWidth = (oWidth + sourceMapFactor - 1) / sourceMapFactor;
	  sourceMapHeight = (oHeight + sourceMapFactor - 1) / sourceMapFactor;
	  sourceMapSize = sourceMapWidth * sourceMapHeight;
#if 0
	  if (updateOutput)
	    {
	      // in order to update, we must successfully read in
	      //   the old source map
	      if (!ReadMap(sourceMapName, &sourceMap, &oldSourceMapLevel,
			   &oldSourceMapWidth, &oldSourceMapHeight,
			   &oldSourceMapMinX, &oldSourceMapMinY,
			   imName0, imName1,
			   msg) ||
		  oldSourceMapLevel != sourceMapLevel ||
		  oldSourceMapWidth != sourceMapWidth ||
		  oldSourceMapHeight != sourceMapHeight ||
		  oldSourceMapMinX != 0 ||
		  oldSourceMapMinY != 0)
		updateOutput = 0;
	    }
	  if (!updateOutput)
	    {
#endif
	      //	      printf("malloc %zu bytes for sourceMap\n",
	      //		     sourceMapSize * sizeof(MapElement));
	      sourceMap = (MapElement *) malloc(sourceMapSize * 
						sizeof(MapElement));
	      if (sourceMap == 0)
		Error("malloc of sourceMap failed; errno = %d\n", errno);

	      for (y = 0; y < sourceMapHeight; ++y)
		for (x = 0; x < sourceMapWidth; ++x)
		  {
		    sourceMap[y*sourceMapWidth + x].x = 0.0;
		    sourceMap[y*sourceMapWidth + x].y = 0.0;
		    sourceMap[y*sourceMapWidth + x].c = -1.0;
		  }
#if 0
	    }
#endif
	}
      else
	{
	  sourceMapWidth = 0;
	  sourceMapHeight = 0;
	  sourceMap = NULL;
	}

      for (hs = 1; ; hs *= 2)
	{
	  endY = oMinY - 1;
	  for (hi = 0; hi < hs; ++hi)
	    {
	      startY = endY + 1;
	      if (cols == 1)
		{
		  endY = oMinY + (hi + 1) * oHeight / hs - 1;
		  endY = ((endY - startY + 1) / reductionFactor) *
		    reductionFactor + startY - 1;
		  outHeight = rows * th;
		  outWidth = tw;
		}
	      else
		{
		  outHeight = ((oHeight / hs + th * reductionFactor - 1) /
			       (th * reductionFactor)) *
		    (th * reductionFactor);
		  if (startY + ((int) outHeight) - 1 > oMaxY)
		    outHeight = oMaxY - startY + 1;
		  outWidth = tw;
		  if (outHeight == 0)
		    Error("Insufficent memory to construct a single row of tiles.\n");
		  endY = startY + ((int) outHeight) - 1;
		}
	      canvasHeight = endY - startY + 1;
	      if (canvasHeight < reductionFactor)
		Error("Insufficent memory to construct a horizontal row\n");

	      memset(increase, 0, oWidth * sizeof(size_t));
	      memset(decrease, 0, oWidth * sizeof(size_t));

	      for (i = startImage; i <= endImage; ++i)
		{
		  iMinX = (int) floor(images[i].minX);
		  iMaxX = (int) ceil(images[i].maxX);
		  iMinY = (int) floor(images[i].minY);
		  iMaxY = (int) ceil(images[i].maxY);
		  if (iMinY > endY || iMaxY < startY)
		    continue;
		  memoryRequired = images[i].width * images[i].height *
		    (sizeof(unsigned char) + sizeof(unsigned char)) +
		    ((images[i].width + 7 ) / 8) * images[i].height +
		    images[i].mapBytes;
		  if (iMinX <= oMinX && iMaxX >= oMinX)
		    increase[0] += memoryRequired;
		  else if (iMinX <= oMaxX && iMaxX >= oMinX)
		    increase[iMinX - oMinX] += memoryRequired;
		  if (iMaxX >= oMinX && iMaxX < oMaxX-1)
		    decrease[iMaxX - oMinX + 1] += memoryRequired;
		}
	      
	      memoryRequired = 0;
	      maxMemoryRequired = 0;
	      for (i = 0; i < oWidth; ++i)
		{
		  memoryRequired += increase[i];
		  memoryRequired -= decrease[i];
		  if (memoryRequired > maxMemoryRequired)
		    maxMemoryRequired = memoryRequired;
		}
	      maxMemoryRequired += sourceMapWidth * sourceMapHeight *
		sizeof(MapElement);
	      maxMemoryRequired += outHeight * outWidth * sizeof(unsigned char);
	      maxMemoryRequired += 2 * oWidth * sizeof(size_t);
	      maxMemoryRequired += 2 * reductionFactor * canvasHeight * sizeof(unsigned char);
	      if (maxMemoryRequired > ((size_t) memoryLimit) * 1000000)
		break;
	    }
	  if (hi >= hs)
	    break;
	}
      if (hs > 1)
	printf("Splitting rendering into %d horizontal strips\n", hs);

      // do horizontal strips one-at-a-time
      endY = oMinY - 1;
      for (hi = 0; hi < hs; ++hi)
	{
	  //	  PrintUsage();

	  startY = endY + 1;
	  if (cols == 1)
	    {
	      endY = oMinY + (int) ((hi + 1) * oHeight / hs) - 1;
	      endY = ((endY - startY + 1) / reductionFactor) *
		reductionFactor + startY - 1;
	      outHeight = rows * th;
	      outWidth = tw;
	      if (hi == 0)
		{
		  //	printf("malloc %zu bytes for out  (cols == 1)\n",
		  //		 outHeight * outWidth *
		  //		 sizeof(unsigned char));
		  out = (unsigned char *) malloc(outHeight * outWidth *
						 sizeof(unsigned char));
		  if (out == 0)
		    Error("malloc of out failed; errno = %d\n", errno);
		}
	      outMinY = oMinY;
	      startRow = 0;
	      endRow = rows-1;
	    }
	  else
	    {
	      outHeight = ((oHeight / hs + th * reductionFactor - 1) /
			   (th * reductionFactor)) *
		(th * reductionFactor);
	      if (startY + ((int) outHeight) - 1 > oMaxY)
		outHeight = oMaxY - startY + 1;
	      outWidth = tw;
	      if (outHeight == 0)
		Error("Insufficent memory to construct a single row of tiles.\n");
	      endY = startY + outHeight - 1;
	      //	      printf("malloc %zu bytes for out (cols != 1)\n",
	      //		     outHeight * outWidth *
	      //		     sizeof(unsigned char));
	      out = (unsigned char *) malloc(outHeight * outWidth *
					     sizeof(unsigned char));
	      if (out == 0)
		Error("malloc of out (2) failed; errno = %d\n", errno);

	      outMinY = startY;
	      startRow = (startY - oMinY) / (th * reductionFactor);
	      endRow = (endY - oMinY) / (th * reductionFactor);
	    }

	  memset(increase, 0, oWidth * sizeof(size_t));
	  memset(decrease, 0, oWidth * sizeof(size_t));
	  for (i = startImage; i <= endImage; ++i)
	    {
	      iMinX = (int) floor(images[i].minX);
	      iMaxX = (int) ceil(images[i].maxX);
	      iMinY = (int) floor(images[i].minY);
	      iMaxY = (int) ceil(images[i].maxY);
	      if (iMinY > endY || iMaxY < startY)
		continue;
	      memoryRequired = images[i].width * images[i].height *
		(sizeof(unsigned char) + sizeof(unsigned char)) +
		((images[i].width + 7 ) / 8) * images[i].height +
		images[i].mapBytes;
	      if (iMinX <= oMinX && iMaxX >= oMinX)
		increase[0] += memoryRequired;
	      else if (iMinX <= oMaxX && iMaxX >= oMinX)
		increase[iMinX - oMinX] += memoryRequired;
	      if (iMaxX >= oMinX && iMaxX < oMaxX-1)
		decrease[iMaxX - oMinX + 1] += memoryRequired;
	    }

	  //	  PrintUsage();

	  imageMem = 0;
	  canvasWidth = 0;
	  canvasHeight = endY - startY + 1;
	  canvasMinX = oMinX;
	  canvasMinY = startY;
	  weightMinX = oMinX;
	  weightMinY = startY;
	  tx = 0;
	  ty = (startY - outMinY) / reductionFactor;
	  tCol = 0;

	  // render from left to right across the output
	  startX = oMinX;
	  do {
	    //	    PrintUsage();
	    // see how far we can increase endX while still satisfying
	    //   memory constraints
	    endX = startX-1;
	    memoryRequired = imageMem;
	    memoryRequired += sourceMapWidth * sourceMapHeight *
	      sizeof(MapElement);
	    memoryRequired += outHeight * outWidth * sizeof(unsigned char);
	    //	    printf("memreq1 = %zu\n", memoryRequired);
	    memoryRequired += 2 * oWidth * sizeof(size_t);
	    //	    printf("memreq2 = %zu\n", memoryRequired);
	    memoryRequired += (startX - canvasMinX) * canvasHeight *
	      sizeof(unsigned char);
	    //	    printf("memreq3 = %zu\n", memoryRequired);
	    if (!overlay)
	      memoryRequired += images[oi].width * images[oi].height *
		(sizeof(unsigned char) + sizeof(unsigned char)) +
		((images[oi].width + 7) / 8) * images[oi].height +
		images[oi].mapBytes;
	    //	    printf("memreq4 = %zu   (%d %d %d %d %d %d)\n", memoryRequired,
	    //		   oi, images[oi].width, images[oi].height, images[oi].mapBytes,
	    //		   startImage, endImage);
	    while (memoryRequired < ((size_t) memoryLimit) * 1000000 &&
		   endX < oMaxX)
	      {
		++endX;
		memoryRequired += canvasHeight *
		  (sizeof(unsigned char) + sizeof(unsigned short));
		if (overlay)
		  {
		    memoryRequired += increase[endX - oMinX];
		    memoryRequired -= decrease[endX - oMinX];
		  }
	      }
	    if (memoryRequired >= ((size_t) memoryLimit) * 1000000)
	      --endX;
	    if (endX < startX)
	      Error("Internal error: not enough memory to construct a single pixel column\n%d %d %d %d %zu",
		    startX, endX, startY, endY, memoryRequired);
	    if (endX < oMaxX &&
		endX - canvasMinX + 1 < reductionFactor)
	      Error("Internal error: not enough memory to construct a column of the output image\n");
	    printf("Rendering area from x=%d to x=%d, y=%d to y=%d\n",
		   startX, endX, startY, endY);

	    //	    PrintUsage();


	    // resize the canvas
	    oldCanvasWidth = canvasWidth;
	    canvasWidth = endX - canvasMinX + 1;
	    //	    printf("realloc canvas from %zu bytes to %zu bytes\n",
	    //		   oldCanvasWidth * canvasHeight * sizeof(unsigned char),
	    //		   canvasWidth * canvasHeight * sizeof(unsigned char));
	    canvas = (unsigned char *) realloc(canvas, canvasWidth * canvasHeight * sizeof(unsigned char));
	    weightWidth = endX - startX + 1;
	    weightHeight = canvasHeight;
	    //	    printf("malloc %zu bytes for weight\n",
	    //		   weightWidth * weightHeight * sizeof(unsigned short));
	    weight = (unsigned short *) malloc(weightWidth * weightHeight * sizeof(unsigned short));
	    if (weight == 0)
 	      Error("malloc of weight failed; errno = %d\n", errno);
	    memset(weight, 0xff, weightWidth * weightHeight * sizeof(unsigned short));

	    // shift any pixels present in the canvas to their new positions
	    for (y = canvasHeight-1; y >= 0; --y)
	      {
		if (oldCanvasWidth > 0)
		  memmove(&canvas[y*canvasWidth], &canvas[y*oldCanvasWidth],
			  oldCanvasWidth * sizeof(unsigned char));
		memset(&canvas[y*canvasWidth+oldCanvasWidth], 0, canvasWidth-oldCanvasWidth);
	      }
	    oldCanvasWidth = canvasWidth;

	    //	    PrintUsage();

	    // paint the old images
	    printf("Applying loaded image maps");
	    fflush(stdout);
	    nProcessed = 0;
	    for (i = startImage; i <= endImage; ++i)
	      if (images[i].minX < startX &&
		  images[i].maxX >= startX &&
		  images[i].minY <= endY &&
		  images[i].maxY >= startY)
		{
		  PaintImage(i, startX, endX, startY, endY);
		  if ((nProcessed % 50) == 0 && nProcessed != 0)
		    printf(" %d\n    ", nProcessed);
		  printf(".");
		  fflush(stdout);
		  ++nProcessed;
		}

	    //	    PrintUsage();

	    // paint the new images
	    printf("\nApplying new image maps");
	    fflush(stdout);
	    nProcessed = 0;
	    for (i = startImage; i <= endImage; ++i)
	      if (images[i].minX >= startX &&
		  images[i].minX <= endX &&
		  images[i].minY <= endY &&
		  images[i].maxY >= startY)
		{
		  PaintImage(i, startX, endX, startY, endY);
		  if ((nProcessed % 50) == 0 && nProcessed != 0)
		    printf(" %d\n    ", nProcessed);
		  printf(".");
		  fflush(stdout);
		  ++nProcessed;
		}

	    //	    PrintUsage();

	    // paint the label if requested
	    if (labelWidth > 0)
	      {
		printf("\nPainting label");
		fflush(stdout);

		if (labelMinX <= endX &&
		    labelMaxX >= startX &&
		    labelMinY <= endY &&
		    labelMaxY >= startY)
		  {
		    len = strlen(label);
		    for (i = 0; i < len; ++i)
		      {
			ci = label[i] - 32;
			if (ci < 0 || ci > 94)
			  continue;
			charMinX = labelMinX + i * labelWidth;
			charMaxX = charMinX + labelWidth - 1;
			if (charMinX > endX ||
			    charMaxX < startX)
			  continue;
			for (dy = 0; dy < labelHeight; ++dy)
			  {
			    y = labelMinY + dy;
			    if (y < startY || y > endY)
			      continue;
			    sy = ((float) dy) / labelHeight * fontHeight;
			    ey = ((float) (dy+1)) / labelHeight * fontHeight;
			    isy = (int) floor(sy);
			    iey = (int) floor(ey);
			    for (dx = 0; dx < labelWidth; ++dx)
			      {
				x = charMinX + dx;
				if (x < startX || x > endX)
				  continue;
				sx = ((float) dx) / labelWidth * fontWidth;
				ex = ((float) (dx+1)) / labelWidth * fontWidth;
				isx = (int) floor(sx);
				iex = (int) floor(ex);
				if (iex >= fontWidth)
				  iex = fontWidth-1;
				ws = 0.0;
				v = 0.0;
				for (iy = isy; iy <= iey; ++iy)
				  {
				    wy = 1.0;
				    if (iy == isy)
				      wy -= sy - isy;
				    if (iy == iey)
				      wy -= iey + 1.0 - ey;
				    for (ix = isx; ix <= iex; ++ix)
				      {
					wx = 1.0;
					if (ix == isx)
					  wx -= sx - isx;
					if (ix == iex)
					  wx -= iex + 1.0 - ex;
					ws += wx * wy;
					v += wx * wy * font[(ci * fontHeight + iy) * fontWidth + ix];
				      }
				  }
				if (ws == 0.0)
				  Error("Internal error: label weight is 0\n");
				ilv = (int) floor(v / ws + 0.5);
				if (ilv != 255)
				  canvas[(y - startY) * canvasWidth + x - canvasMinX] = 255-ilv;
			      }
			  }
		      }
		  }
	      }

	    // transfer to output tile(s)
	    printf("\nSaving rendered area");
	    fflush(stdout);
	    nProcessed = 0;
	    cx = 0;
            a = reductionFactor * reductionFactor;
            offset = a / 2;
	    ny = canvasHeight / reductionFactor;
	    for (;;)
	      {
		nx = (((int) canvasWidth) - cx) / reductionFactor;
		if (tw - tx < nx)
		  nx = tw - tx;
		if (nx <= 0)
		  break;

		if (reductionFactor > 1)
		  for (y = 0; y < ny; ++y)
		    {
		      iy = y * reductionFactor;
		      for (x = 0; x < nx; ++x)
			{
			  ix = x * reductionFactor + cx;
			  sum = 0;
			  valid = 1;
			  for (dy = 0; dy < reductionFactor; ++dy)
			    for (dx = 0; dx < reductionFactor; ++dx)
			      {
				iv = canvas[(iy+dy)*canvasWidth + (ix+dx)];
				if (iv != 0)
				  sum += iv;
				else
				  valid = 0;
			      }
			  if (valid)
			    out[(ty+y)*tw+tx+x] = (sum + offset) / a;
			  else
			    out[(ty+y)*tw+tx+x] = 0;
			}
		    }
		else
		  for (y = 0; y < ny; ++y)
		    memcpy(&out[(ty+y)*tw+tx],
			   &canvas[y*canvasWidth + cx],
			   nx);

		if (hi == hs-1)
		  {
		    // last stripe, so blacken tile bottom
		    for (y = ty+ny; y < outHeight; ++y)
		      memset(&out[y*tw+tx], 0, nx);
		  }
	    
		cx += nx * reductionFactor;
		tx += nx;
		if (endX == oMaxX)
		  {
		    // blacken tile right side
		    if (tw > tx)
		      for (y = 0; y < outHeight; ++y)
			memset(&out[y*tw+tx], 0, tw - tx);
		    tx = tw;
		  }

		if (tx == tw &&
		    (cols == 1 && hi == hs-1 ||
		     cols > 1))
		  {
		    WriteTiles(tCol, startRow, endRow,
			       images[startImage].name);
		    tx = 0;
		    ++tCol;
		  }
	      }
	    printf("\n");

	    //	    PrintUsage();

	    // save extra pixel columns
	    canvasWidth -= cx;
	    canvasMinX += cx;
	    if (canvasWidth > 0)
	      for (y = 0; y < canvasHeight; ++y)
		memmove(&canvas[y*canvasWidth],
			&canvas[y*oldCanvasWidth + cx],
			canvasWidth * sizeof(unsigned char));
	    //	    printf("realloc canvas from %zu bytes to %zu bytes\n",
	    //		   oldCanvasWidth * canvasHeight * sizeof(unsigned char),
	    //		   canvasWidth * canvasHeight * sizeof(unsigned char));
	    canvas = (unsigned char *)
	      realloc(canvas,
		      canvasWidth * canvasHeight * sizeof(unsigned char));
	    weightMinX += (int) weightWidth;

	    //	    printf("freeing %zu bytes from weight\n",
	    //		   weightHeight * weightWidth * sizeof(unsigned short));
	    free(weight);

	    startX = endX + 1;

	    //	    PrintUsage();
	  } while (startX <= oMaxX);

	  if (cols > 1)
	    {
	      //	      printf("freeing %zu bytes from out (cols > 1)\n",
	      //		     outHeight * outWidth *
	      //		     sizeof(unsigned char));
	      free(out);
	    }
	}

      //      PrintUsage();

      if (cols == 1)
	{
	  //	  printf("freeing %zu bytes from out (cols == 1)\n",
	  //		 outHeight * outWidth *
	  //		 sizeof(unsigned char));
	  free(out);
	}

      if (sourceMap != NULL)
	{
	  if (!CreateDirectories(sourceMapName))
	    Error("Could not create directories for source map file %s\n",
		  sourceMapName);
	  if (!WriteMap(sourceMapName, sourceMap, sourceMapLevel,
			sourceMapWidth, sourceMapHeight,
			0, 0,
			mapsName, imageListName,
			UncompressedMap, msg))
	    Error("Could not write source map %s:\n%s\n", sourceMapName, msg);
	  //	  printf("freeing %zu bytes from sourceMap\n",
	  //		     sourceMapSize * sizeof(MapElement));
	  free(sourceMap);
	  sourceMap = NULL;
	}
    }

  //  PrintUsage();

  printf("\nWriting size file... ");
  fflush(stdout);

  if (strlen(outputName) == 0)
    strcpy(fn, "size");
  else if (outputName[strlen(outputName)-1] == '/')
    sprintf(fn, "%ssize", outputName);
  else
    sprintf(fn, "%s.size", outputName);
  f = fopen(fn, "w");
  if (f == NULL)
    Error("Could not open size file %s for writing.\n", fn);
  fprintf(f, "%d %d\n%zdx%zd%+d%+d\n%zdx%zd\n",
	  rows, cols,
	  oWidth, oHeight, oMinX, oMinY,
	  oWidth / reductionFactor, oHeight / reductionFactor);
  fclose(f);
  printf("done.\n");
  fflush(stdout);
  return(0);
}

void
PaintImage (int i, int minX, int maxX, int minY, int maxY)
{
  int j;
  int iw, ih;
  int x, y;
  InverseMap *invMap;
  float xv, yv;
  int ixv, iyv;
  float rx, ry;
  float rrx, rry;
  float rv, dv;
  size_t maskBytes;
  size_t mbpl;
  unsigned char *mask;
  unsigned char *dist;
  float *distance;
  float dst;
  int idst;
  float cx, cy;
  float offset;
  float r00, r01, r10, r11;
  float d00, d01, d10, d11;
  int nx, ny;
  int pMinX, pMaxX, pMinY, pMaxY;
  unsigned char *image;
  float w;
  float d, d2;
  int row, col;
  int ix, iy;
  float mFactor;
  int offsetX, offsetY;
  int v;
  size_t newCanvasSize;
  int mxMin, myMin;
  char imName0[PATH_MAX];
  char imName1[PATH_MAX];
  char msg[PATH_MAX+256];
  char fn[PATH_MAX];
  struct stat statBuf;
  int complete;
  int iixv, iiyv;
  float rb00, rb01, rb10, rb11;
  float rw00, rw01, rw10, rw11;
  float rb, rw;
  MapElement *imap;
  float imapFactor;
  int imapw, imaph;
  int imapXMin, imapYMin;
  char imapName0[PATH_MAX], imapName1[PATH_MAX];
  float xvi, yvi;
  int sourceMapMask;
  int smi;
  int targetMapSize;
  MapElement *targetMap;
  int targetMapWidth, targetMapHeight;
  int targetMapMask;
  int tmi;
  int testX, testY;
  int mw, mh;
  MapElement *map;
  float spacing;
  float rxp, ryp;
  unsigned char *imask;
  size_t imbpl;
  int warned;

  /* read in map if necessary */
  if (images[i].map == NULL)
    {
      sprintf(fn, "%s%s.map", mapsName, images[i].name);
      if (!ReadMap(fn, &(images[i].map), &(images[i].mLevel),
		   &(images[i].mw), &(images[i].mh),
		   &(images[i].mxMin), &(images[i].myMin),
		   imName0, imName1,
		   msg))
	Error("Could not read map %s:\n  error: %s\n",
	      fn, msg);
      if (rotation != 0.0)
	{
	  map = images[i].map;
	  spacing = (1 << images[i].mLevel) * mapScale;
	  mw = images[i].mw;
	  mh = images[i].mh;
	  for (y = 0; y < mh; ++y)
	    for (x = 0; x < mw; ++x)
	      {
		rx = map[y*mw+x].x * spacing;
		ry = map[y*mw+x].y * spacing;
		rxp = rx - rotationX;
		ryp = ry - rotationY;
		rx = cosRot * rxp + sinRot * ryp + rotationX;
		ry = -sinRot * rxp + cosRot * ryp + rotationY;
		map[y*mw+x].x = rx / spacing;
		map[y*mw+x].y = ry / spacing;
	      }
	}

      imageMem += images[i].mapBytes;
    }
  if (images[i].invMap == NULL)
    images[i].invMap = InvertMap(images[i].map,
				 images[i].mw,
				 images[i].mh);
	
  /* read in image if necessary */
  if (images[i].image == NULL)
    {
      sprintf(fn, "%s%s", imagesName, images[i].name);
      if (!ReadImage(fn, &(images[i].image),
		     &iw, &ih,
		     -1, -1, -1, -1,
		     msg))
	Error("Could not read image %s:\n  error: %s\n",
	      fn, msg);
      if (iw != images[i].width ||
	  ih != images[i].height)
	Error("Dimensions of image %s do not match those in images list.\n",
	      fn);
      imageMem += images[i].width * images[i].height;

      if (targetMapsName[0] != '\0')
	{
	  targetMapWidth = (iw + targetMapsFactor - 1) / targetMapsFactor;
	  targetMapHeight = (ih + targetMapsFactor - 1) / targetMapsFactor;
	  targetMapSize = targetMapWidth * targetMapHeight;
	  //	  printf("malloc %zu bytes for targetMap\n",
	  //		 targetMapSize * sizeof(MapElement));
	  images[i].targetMap =
	    (MapElement *) malloc(targetMapSize * sizeof(MapElement));
 	  if (images[i].targetMap == 0)
 	    Error("malloc of targetMap failed; errno = %d\n", errno);
	  memset(images[i].targetMap, 0,
		 targetMapSize * sizeof(MapElement));
	}
      else
	images[i].targetMap = NULL;
    }
  targetMapWidth = (images[i].width + targetMapsFactor - 1) / targetMapsFactor;
  targetMapHeight = (images[i].height + targetMapsFactor - 1) / targetMapsFactor;
  targetMapSize = targetMapWidth * targetMapHeight;
  targetMap = images[i].targetMap;
  targetMapMask = targetMapsFactor - 1;

  /* read in mask if necessary */
  if (images[i].mask == NULL)
    if (masksName[0] != '\0')
      {
	sprintf(fn, "%s%s", masksName, images[i].name);
	if (!ReadBitmap(fn, &(images[i].mask),
			&iw, &ih,
			-1, -1, -1, -1,
			msg))
	  Error("Could not read image %s:\n  error: %s\n",
		fn, msg);
	if (maskScale == 1.0 &&
	    (iw != images[i].width ||
	     ih != images[i].height) ||
	    maskScale < 1.0 &&
	    ((int) floor(maskScale * iw + 0.5) != images[i].width ||
	     (int) floor(maskScale * ih + 0.5) != images[i].height) ||
	    maskScale > 1.0 &&
	    ((int) floor(images[i].width / maskScale) != iw ||
	     (int) floor(images[i].height / maskScale) != ih))
	  Error("Dimensions of mask %s do not match those in images list (%f %d %d %d %d)\n",
		fn, maskScale, iw, ih, images[i].width, images[i].height);
	mbpl = (images[i].width + 7) / 8;
	maskBytes = mbpl * images[i].height;
	if (maskScale != 1.0)
	  {
	    mask = (unsigned char *) malloc(maskBytes);
	    memset(mask, 0, maskBytes);
	    imask = images[i].mask;
	    imbpl = (iw + 7) / 8;
	    for (y = 0; y < images[i].height; ++y)
	      {
		iy = (int) floor(y / maskScale + 0.001);
		for (x = 0; x < images[i].width; ++x)
		  {
		    ix = (int) floor(x / maskScale + 0.001);
		    if (imask[iy*imbpl + (ix >> 3)] & (0x80 >> (ix & 7)))
		      mask[y*mbpl + (x >> 3)] |= 0x80 >> (x & 7);
		  }
	      }
	    free(images[i].mask);
	    images[i].mask = mask;
	  }
	
	mask = images[i].mask;
	for (j = 0; j < maskBytes; ++j)
	  mask[j] ^= 0xff;
	imageMem += maskBytes;
      }
    else
      {
	maskBytes = ((images[i].width + 7) / 8) * images[i].height;
	//	printf("malloc %zu bytes for mask\n", maskBytes);
	mask = images[i].mask = (unsigned char *) malloc(maskBytes);
 	if (mask == 0)
 	  Error("malloc of mask failed; errno = %d\n", errno);
	memset(images[i].mask, 0, maskBytes);
	imageMem += maskBytes;
      }
		
  /* read in intensity map if necessary */
  if (images[i].imap == NULL)
    if (imapsName[0] != '\0')
      {
	sprintf(fn, "%s%s.map", imapsName, images[i].name);
	if (!ReadMap(fn, &(images[i].imap), &(images[i].imapLevel),
		     &(images[i].imapw), &(images[i].imaph),
		     &imapXMin, &imapYMin,
		     imapName0, imapName1,
		     msg))
	  Error("Could not read map %s:\n  error: %s\n",
		fn, msg);
	if (imapXMin != 0 || imapYMin != 0)
	  Error("Can not handle partial intensity map: %s\n", fn);
	imap = images[i].imap;
	imageMem += images[i].imapw * images[i].imaph * sizeof(MapElement);
      }
		
  /* compute distance table if necessary */
  if (images[i].dist == NULL)
    {
      iw = images[i].width;
      ih = images[i].height;
      //      printf("malloc %zu bytes for dist\n", iw * ih * sizeof(unsigned char));
      images[i].dist = (unsigned char*) malloc(iw * ih * sizeof(unsigned char));
      if (images[i].dist == 0)
 	Error("malloc of images[i].dist failed; errno = %d\n", errno);
      distance = (float*) malloc(iw * ih * sizeof(float));
      if (distance == 0)
 	Error("malloc of distance failed; errno = %d\n", errno);
      computeDistance(EUCLIDEAN_DISTANCE,
		      iw, ih, images[i].mask,
		      distance);

      /* if distance from edge is less, use that */
      dist = images[i].dist;
      for (y = 0; y < ih; ++y)
	for (x = 0; x < iw; ++x)
	  {
	    dst = distance[y*iw+x];
	    d = x + 1;
	    if (d < dst)
	      dst = d;
	    d = iw - x;
	    if (d < dst)
	      dst = d;
	    d = y + 1;
	    if (d < dst)
	      dst = d;
	    d = ih - y;
	    if (d < dst)
	      dst = d;
	    idst = (int) floor(64.0 * dst);
	    if (idst > 255)
	      idst = 255;
	    else if (idst < 0)
	      idst = 0;
	    dist[y*iw + x] = idst;
	  }
      free(distance);
      imageMem += iw * ih;
    }

  /* paint the image on the canvas */
  pMinX = minX;
  pMaxX = maxX;
  pMinY = minY;
  pMaxY = maxY;
  if (images[i].minX > pMinX)
    pMinX = images[i].minX;
  if (images[i].maxX < pMaxX)
    pMaxX = images[i].maxX;
  if (images[i].minY > pMinY)
    pMinY = images[i].minY;
  if (images[i].maxY < pMaxY)
    pMaxY = images[i].maxY;

  mFactor = (1 << images[i].mLevel) * mapScale;
  invMap = images[i].invMap;
  image = images[i].image;
  mask = images[i].mask;
  imap = images[i].imap;
  if (imap != NULL)
    {
      imapFactor = (1 << images[i].imapLevel) * imapScale;
      imapw = images[i].imapw;
      imaph = images[i].imaph;
    }
  dist = images[i].dist;
  iw = images[i].width;
  ih = images[i].height;
  mxMin = images[i].mxMin;
  myMin = images[i].myMin;
  //  printf("mxMin = %d myMin = %d\n", mxMin, myMin);
  cx = (iw - 1) / 2.0;
  cy = (ih - 1) / 2.0;
  mbpl = (iw + 7) / 8;
  offset = -1000000.0;
  printf("Rendering %s from x=%d to %d y=%d to %d\n",
	 images[i].name, pMinX, pMaxX, pMinY, pMaxY);
  testX = (pMinX + pMaxX) / 2;
  testY = (pMinY + pMaxY) / 2;
  warned = 0;
  for (y = pMinY; y <= pMaxY; ++y)
    for (x = pMinX; x <= pMaxX; ++x)
      {
	//	if (y == testY && x == testX)
	//	  printf("TEST STARTED\n");
	if (!Invert(invMap, &xv, &yv, (x + 0.5) / mFactor, (y + 0.5) / mFactor))
	  {
	    //	    if (y == testY && x == testX)
	    //	      printf("TEST INVERT FAILED %f %f %f %f\n", (x + 0.5)/mFactor,
	    //		     (y+0.5)/mFactor, xv, yv);
	    continue;
	  }
	xv = (xv + mxMin) * mFactor - 0.5;
	yv = (yv + myMin) * mFactor - 0.5;
	//	if (y == testY && x == testX)
	//	  printf("TEST INVERT OK %f %f %f %f\n", (x + 0.5)/mFactor,
	//		 (y+0.5)/mFactor, xv, yv);
	ixv = (int) floor(xv);
	iyv = (int) floor(yv);
	if (ixv < -1 || ixv >= iw ||
	    iyv < -1 || iyv >= ih)
	  continue;
	rrx = xv - ixv;
	rry = yv - iyv;
	if (ixv >= 0 && iyv >= 0 &&
	    !(mask[iyv * mbpl + (ixv >> 3)] & (0x80 >> (ixv & 7))))
	  {
	    r00 = image[iyv * iw + ixv];
	    d00 = dist[iyv * iw + ixv];
	    if (targetMap != NULL && ((ixv | iyv) & targetMapMask) == 0)
	      {
		tmi = (iyv >> targetMapsLevel) * targetMapWidth +
		  (ixv >> targetMapsLevel);
		if (tmi >= targetMapSize)
		  Error("tmi out-of-bounds: %d %d\n%d %d\n",
			tmi, targetMapSize,
			targetMapWidth, targetMapHeight);
		targetMap[tmi].x = (float) x;
		targetMap[tmi].y = (float) y;
		targetMap[tmi].c = 1.0;
	      }
	  }
	else
	  {
	    r00 = 0.0;
	    d00 = 0.0;
	  }
	if (ixv >= 0 && iyv+1 < ih &&
	    !(mask[(iyv + 1) * mbpl + (ixv >> 3)] & (0x80 >> (ixv & 7))))
	  {
	    r01 = image[(iyv + 1) * iw + ixv];
	    d01 = dist[(iyv + 1) * iw + ixv];
	    if (targetMap != NULL && ((ixv | (iyv+1)) & targetMapMask) == 0)
	      {
		tmi = ((iyv + 1) >> targetMapsLevel) * targetMapWidth +
		  (ixv >> targetMapsLevel);
		if (tmi >= targetMapSize)
		  Error("tmi out-of-bounds: %d %d\n%d %d\n",
			tmi, targetMapSize,
			targetMapWidth, targetMapHeight);
		targetMap[tmi].x = (float) x;
		targetMap[tmi].y = (float) y;
		targetMap[tmi].c = 1.0;
	      }
	  }
	else
	  {
	    r01 = 0.0;
	    d01 = 0.0;
	  }
	if (ixv+1 < iw && iyv >= 0 &&
	    !(mask[iyv * mbpl + ((ixv+1) >> 3)] & (0x80 >> ((ixv+1) & 7))))
	  {
	    r10 = image[iyv * iw + ixv + 1];
	    d10 = dist[iyv * iw + ixv + 1];
	    if (targetMap != NULL && (((ixv + 1) | iyv) & targetMapMask) == 0)
	      {
		tmi = (iyv >> targetMapsLevel) * targetMapWidth +
		  ((ixv + 1) >> targetMapsLevel);
		if (tmi >= targetMapSize)
		  Error("tmi out-of-bounds: %d %d\n%d %d\n",
			tmi, targetMapSize,
			targetMapWidth, targetMapHeight);
		targetMap[tmi].x = (float) x;
		targetMap[tmi].y = (float) y;
		targetMap[tmi].c = 1.0;
	      }
	  }
	else
	  {
	    r10 = 0.0;
	    d10 = 0.0;
	  }
	if (ixv+1 < iw && iyv+1 < ih &&
	    !(mask[(iyv+1) * mbpl + ((ixv+1) >> 3)] & (0x80 >> ((ixv+1) & 7))))
	  {
	    r11 = image[(iyv + 1) * iw + ixv + 1];
	    d11 = dist[(iyv + 1) * iw + ixv + 1];
	    if (targetMap != NULL && (((ixv + 1) | (iyv+1)) & targetMapMask) == 0)
	      {
		tmi = ((iyv + 1) >> targetMapsLevel) * targetMapWidth +
		  ((ixv + 1) >> targetMapsLevel);
		if (tmi >= targetMapSize)
		  Error("tmi out-of-bounds: %d %d\n%d %d\n",
			tmi, targetMapSize,
			targetMapWidth, targetMapHeight);
		targetMap[tmi].x = (float) x;
		targetMap[tmi].y = (float) y;
		targetMap[tmi].c = 1.0;
	      }
	  }
	else
	  {
	    r11 = 0.0;
	    d11 = 0.0;
	  }
	rv = r00 * (rrx - 1.0) * (rry - 1.0)
	  - r10 * rrx * (rry - 1.0) 
	  - r01 * (rrx - 1.0) * rry
	  + r11 * rrx * rry;
	dv = d00 * (rrx - 1.0) * (rry - 1.0)
	  - d10 * rrx * (rry - 1.0) 
	  - d01 * (rrx - 1.0) * rry
	  + d11 * rrx * rry;
	if (dv <= 0.0)
	  continue;

	if (imap != NULL)
	  {
	    /* lookup xv,yv in intensity map, if present */
	    xvi = (xv + 0.5) / imapFactor;
	    yvi = (yv + 0.5) / imapFactor;
	    iixv = (int) floor(xvi);
	    iiyv = (int) floor(yvi);
	    rrx = xvi - iixv;
	    rry = yvi - iiyv;
	    while (iixv < 0)
	      {
		++iixv;
		rrx -= 1.0;
	      }
	    while (iixv >= imapw-1)
	      {
		--iixv;
		rrx += 1.0;
	      }
	    while (iiyv < 0)
	      {
		++iiyv;
		rry -= 1.0;
	      }
	    while (iiyv >= imaph-1)
	      {
		--iiyv;
		rry += 1.0;
	      }
	    rb00 = imap[iiyv*imapw+iixv].x;
	    rw00 = imap[iiyv*imapw+iixv].y;
	    rb01 = imap[(iiyv+1)*imapw+iixv].x;
	    rw01 = imap[(iiyv+1)*imapw+iixv].y;
	    rb10 = imap[iiyv*imapw+iixv+1].x;
	    rw10 = imap[iiyv*imapw+iixv+1].y;
	    rb11 = imap[(iiyv+1)*imapw+iixv+1].x;
	    rw11 = imap[(iiyv+1)*imapw+iixv+1].y;
	    rb = rb00 * (rrx - 1.0) * (rry - 1.0)
	      - rb10 * rrx * (rry - 1.0) 
	      - rb01 * (rrx - 1.0) * rry
	      + rb11 * rrx * rry;
	    rw = rw00 * (rrx - 1.0) * (rry - 1.0)
	      - rw10 * rrx * (rry - 1.0) 
	      - rw01 * (rrx - 1.0) * rry
	      + rw11 * rrx * rry;
	    if (rw <= rb)
	      {
		if (!warned)
		  {
		    printf("Warning: rw level (%f) is less than rb level (%f) for image %s at (%d %d)\n",
			   rw, rb, images[i].name, ixv, iyv);
		    warned = 1;
		  }
		rw = 255.0;
		rb = 0.0;
	      }
	    rv = (rv - rb) / (rw - rb);
	  }

	if (dv < 64.0)
	  w = 65535 - (int) floor(4.0 * dv);
	else
	  w = (int) floor(hypotf(xv - cx, yv - cy));
	if (w < weight[(y - minY) * weightWidth + x - weightMinX])
	  {
	    weight[(y - minY) * weightWidth + x - weightMinX] = w;
	    if (imap != NULL)
	      v = (int) floor(255.0 * rv + 0.5);
	    else
	      v = (int) floor(255.0 * (rv - blackValue) / range + 0.5);
	    if (v <= 0)
	      {
		if (w == 65535)
		  v = 0;
		else
		  v = 1;
	      }
	    else if (v > 255)
	      v = 255;
	    canvas[(y - minY) * canvasWidth + x - canvasMinX] = v;

	    if (sourceMap != NULL &&
		(((x - oMinX) | (y - oMinY)) & sourceMapMask) == 0)
	      {
		smi = ((y - oMinY) >> sourceMapLevel) * sourceMapWidth +
		  ((x - oMinX) >> sourceMapLevel);
		if (smi >= sourceMapSize)
		  Error("smi out-of-bounds\n");
		sourceMap[smi].x = (float) ixv;
		sourceMap[smi].y = (float) iyv;
		sourceMap[smi].c = (float) i;
	      }
	  }
      }

  /* free up if no longer required */
  if (images[i].maxX <= maxX || maxX >= oMaxX)
    {
      if (images[i].targetMap != NULL)
	{
	  sprintf(fn, "%s%s.map", targetMapsName, images[i].name);
	  if (!CreateDirectories(fn))
	    Error("Could not create directories for target map file %s\n",
		  fn);
	  if (!WriteMap(fn, images[i].targetMap, targetMapsLevel,
			targetMapWidth, targetMapHeight,
			0, 0,
			images[i].name, outputName,
			UncompressedMap, msg))
	    Error("Could not write target map %s:\n%s\n", fn, msg);
	  //	  printf("freeing %zu bytes from targetMap\n",
	  //		 targetMapHeight * targetMapWidth * sizeof(MapElement));
	  free(images[i].targetMap);
	  images[i].targetMap = NULL;
	}
      if (images[i].invMap != NULL)
	{
	  FreeInverseMap(images[i].invMap);
	  images[i].invMap = NULL;
	}
      if (images[i].map != NULL)
	{
	  free(images[i].map);
	  images[i].map = NULL;
	  imageMem -= images[i].mapBytes; // this accounts for the inverse
	                                  // map and target map as well
	}
      if (images[i].imap != NULL)
	{
	  free(images[i].imap);
	  images[i].imap = NULL;
	  imageMem -= images[i].imapw * images[i].imaph * sizeof(MapElement);
	}
      if (images[i].image != NULL)
	{
	  //	  printf("freeing %zu bytes from image\n",
	  //		 images[i].height * images[i].width * sizeof(unsigned char));
	  free(images[i].image);
	  images[i].image = NULL;
	  imageMem -= images[i].width * images[i].height;
	}
      if (images[i].mask != NULL)
	{
	  //	  printf("freeing %zu bytes from mask\n",
	  //		 (size_t) ((images[i].width + 7) / 8) * images[i].height);
	  free(images[i].mask);
	  images[i].mask = NULL;
	  imageMem -= ((images[i].width + 7) / 8) * images[i].height;
	}
      if (images[i].dist != NULL)
	{
	  //	  printf("freeing %zu bytes from dist\n",
	  //		 (size_t) (images[i].width  * images[i].height));
	  free(images[i].dist);
	  images[i].dist = NULL;
	  imageMem -= images[i].width * images[i].height;
	}
    }
}

void
WriteTiles (int col, int startRow, int endRow, char *iName)
{
  int row;
  char fn[PATH_MAX];
  char msg[PATH_MAX+256];

  for (row = startRow; row <= endRow; ++row)
    {
      if (overlay)
	{
	  if (tileWidth < 0 && tileHeight < 0)
	    sprintf(fn, "%s.tif", outputName);
	  else
	    sprintf(fn, "%sc%.2d%sr%.2d.tif", outputName, col+1,
		    tree ? "/" : "", row+1);
	}
      else
	{
	  if (tileWidth < 0 && tileHeight < 0)
	    sprintf(fn, "%s%s.tif",
		    outputName, iName);
	  else
	    sprintf(fn, "%s%s/c%.2d%sr%.2d.tif",
		    outputName, iName,
		    col+1, tree ? "/" : "", row+1);
	}
      if (!CreateDirectories(fn))
	Error("Could not create directories for output file %s\n", fn);
      if (!WriteImage(fn, &out[(row - startRow) * th * tw], (int) tw, (int) th,
		      compress ? HDiffDeflateImage : UncompressedImage,
		      msg))
	Error("Could not write output file %s:\n  error: %s\n",
	      fn, msg);
      if ((nProcessed % 50) == 0 && nProcessed != 0)
	printf(" %d\n    ", nProcessed);
      printf(".");
      fflush(stdout);
      ++nProcessed;
    }
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

unsigned int
Hash (char *s)
{
  unsigned int v = 0;
  char *p = s;
  while (*p != '\0')
    v = 37 * v + *p++;
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
	    Error("Directory hash table is full!\n");
	}
      if (dirHash[i] != NULL)
	continue;
      dirHash[i] = (char *) malloc(strlen(dn)+1);
      if (dirHash[i] == 0)
	Error("malloc of dirHash[i] failed; errno = %d\n", errno);
      strcpy(dirHash[i], dn);

      if (stat(dn, &statBuf) == 0)
	{
	  if (S_ISDIR(statBuf.st_mode) ||
	      S_ISLNK(statBuf.st_mode))
	    continue;
	  Error("Output path component %s is not a directory\n", dn);
	}
      if (errno != ENOENT)
	Error("Could not stat directory %s\n", dn);
      
      if (mkdir(dn, 0777) != 0)
	Error("Could not create directory %s\n", dn);
    }
  return(1);
}

void
PrintUsage ()
{
  struct rusage usage;

  printf("imageMem = %zu\n", imageMem);
  getrusage(RUSAGE_SELF, &usage);
  printf("getrusage_maxrss = %ld\n", usage.ru_maxrss);
}
