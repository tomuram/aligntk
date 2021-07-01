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
#include "imio.h"
#include "invert.h"
#include "dt.h"

#define LINE_LENGTH	255

/* TYPES */
typedef struct Image
{
  int next;		/* index of next image in hash bucket */
  char *name;   	/* name of this image */
  int width, height;    /* width and height in pixels */
  unsigned char *image; /* the image bytes (NULL if not loaded) */
  unsigned char *mask;  /* the image mask, if present */
  float *dist;          /* the distance array -- each element holds the distance
			   from the corresponding image pixel to the closest
			   masked pixel */
  float minX, maxX;     /* bounds of the area this image covers in the */
  float minY, maxY;	/*   final image */
  int needed;		/* true if this image is needed for the current set of
			   tiles */
  int mapBytes;         /* size of maps in bytes */
  int mLevel;		/* map level */
  int mw, mh;		/* map width and height */
  MapElement *map;      /* map of this image into the final image */
  InverseMap *imap;     /* inverse map that translates points in the final image
			   into points in this image */
} Image;

typedef struct CanvasElement
{
  float value;
  float weight;
} CanvasElement;


/* GLOBAL VARIABLES */
FILE *logFile = 0;

char imageListName[PATH_MAX];
char imagesName[PATH_MAX];
char masksName[PATH_MAX];
char mapsName[PATH_MAX];
char outputName[PATH_MAX];
int blend = 1;
int margin = -1;
int tileWidth = -1;
int tileHeight = -1;
int memoryLimit = 1024;   /* 1 GB */
int tree = 0;
float blackValue = 0.0;
float whiteValue = 255.0;
int compress = 0;

int nImages = 0;
Image *images = 0;
int *imageHashTable = 0;
CanvasElement *canvas = 0;
size_t canvasSize = 0;
size_t tw, th;
int oMinX, oMinY, oMaxX, oMaxY;
size_t oWidth, oHeight;
float range;
unsigned char *out;

#define DIR_HASH_SIZE	8192
char *dirHash[DIR_HASH_SIZE];

#define TOP	0
#define BOTTOM	1
#define LEFT	2
#define RIGHT	3

/* FORWARD DECLARATIONS */
void PartitionTiles (int startCol, int endCol, int startRow, int endRow,
		     int entryDir, int exitDir);
void GenerateTiles (int startCol, int endCol, int startRow, int endRow);
void Error (char *fmt, ...);
void Log (char *fmt, ...);
unsigned int Hash (char *s);
int CreateDirectories (char *fn);

int
main (int argc, char **argv, char **envp)
{
  int i, j;
  int n;
  int error;
  char fn[PATH_MAX];
  MapElement *map;
  int x, y;
  float iMinX, iMinY, iMaxX, iMaxY;
  float xv, yv;
  int spacing;
  float ax, bx, ay, by;
  int mLevel;
  int mw, mh;
  int mxMin, myMin;
  char imageName[PATH_MAX];
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
  int nProcessed;
  int ixv, iyv;
  float rx, ry;
  float rrx, rry;
  float rx00, rx01, rx10, rx11, ry00, ry01, ry10, ry11;
  float rv;
  cpu_set_t cpumask;

  // FIX -- TEMPORARY HACK
  CPU_ZERO(&cpumask);
  CPU_SET(7, &cpumask);
  sched_setaffinity(0, sizeof(cpu_set_t), &cpumask);

  sprintf(fn, "logs/applymap3.log");
  logFile = fopen(fn, "w");
  error = 0;
  imageListName[0] = '\0';
  imagesName[0] = '\0';
  masksName[0] = '\0';
  mapsName[0] = '\0';
  outputName[0] = '\0';
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
    else if (strcmp(argv[i], "-images") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(imagesName, argv[i]);
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
    else if (strcmp(argv[i], "-maps") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(mapsName, argv[i]);
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
      {
	if (++i == argc || sscanf(argv[i], "%d", &compress) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else error = 1;

  if (error)
    {
      fprintf(stderr, "Usage: applymap3 -image_list list_file -images image_prefix\n");
      fprintf(stderr, "              -maps map_prefix -output file_prefix\n");
      fprintf(stderr, "              [-masks mask_prefix]\n");
      fprintf(stderr, "              [-blend]\n");
      fprintf(stderr, "              [-mosaic]\n");
      fprintf(stderr, "              [-margin margin_thickness_in_pixels]\n");
      fprintf(stderr, "              [-tile tilewidthxtileheight]\n");
      fprintf(stderr, "              [-memory memory_limit_in_MB]\n");
      fprintf(stderr, "              [-tree]\n");
      fprintf(stderr, "              [-black black_value]\n");
      fprintf(stderr, "              [-white white_value]\n");
      fprintf(stderr, "              [-compress]\n");
      exit(1);
    }

  /* check that at least minimal parameters were supplied */
  if (imageListName[0] == '\0' || imagesName[0] == '\0' ||
      mapsName[0] == '\0' || outputName[0] == '\0')
    {
      fprintf(stderr, "-image_list, -images, -maps, and -output parameters must be specified.\n");
      exit(1);
    }
  range = whiteValue - blackValue;
  if (range == 0.0)
    Error("White value cannot be same as black value\n");

  /* read the image list */
  f = fopen(imageListName, "r");
  if (f == NULL)
    Error("Could not open file %s for reading\n", imageListName);
  nImages = 0;
  imagesSize = 0;
  while (fgets(line, LINE_LENGTH, f) != NULL)
    {
      nItems = sscanf(line, "%s%d%d",
		      imageName, &width, &height);
      if (nItems != 3)
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
      images[nImages].imap = NULL;

      ++nImages;
    }
  fclose(f);
  images = (Image *) realloc(images, nImages * sizeof(Image));
  printf("nImages = %d\n", nImages);

  printf("Previewing maps: ");
  fflush(stdout);

  oMinX = 1000000000;
  oMaxX = -1000000000;
  oMinY = 1000000000;
  oMaxY = -1000000000;
  map = NULL;
  nProcessed = 0;
  for (i = 0; i < nImages; ++i)
    {
      sprintf(fn, "%s%s.map", mapsName, images[i].name);
      if (!ReadMap(fn, &map, &mLevel,
		   &mw, &mh, &mxMin, &myMin,
		   imName0, imName1,
		   msg))
	Error("Could not read map %s:\n  error: %s\n",
	      fn, msg);

      spacing = 1 << mLevel;
      iMinX = 1000000000;
      iMaxX = -1000000000;
      iMinY = 1000000000;
      iMaxY = -1000000000;
      for (dir = 0; dir < 4; ++dir)
	{
	  switch (dir)
	    {
	    case 0:
	      n = images[i].width+1;
	      ax = 1.0;
	      bx = 0.0;
	      ay = 0.0;
	      by = 0.0;
	      break;

	    case 1:
	      n = images[i].width+1;
	      ax = 1.0;
	      bx = 0.0;
	      ay = 0.0;
	      by = images[i].height;
	      break;

	    case 2:
	      n = images[i].height-1;
	      ax = 0.0;
	      bx = 0.0;
	      ay = 1.0;
	      by = 1.0;
	      break;
	      
	    case 3:
	      n = images[i].height-1;
	      ax = 0.0;
	      bx = images[i].width;
	      ay = 1.0;
	      by = 1.0;
	      break;
	    }

	  for (j = 0; j < n; ++j)
	    {
	      xv = (ax * j + bx) / spacing;
	      yv = (ay * j + by) / spacing;
	      ixv = (int) floor(xv);
	      iyv = (int) floor(yv);
	      rrx = xv - ixv;
	      rry = yv - iyv;
	      ixv -= mxMin;
	      iyv -= myMin;
	      if (ixv < 0 || ixv >= mw-1 ||
		  iyv < 0 || iyv >= mh-1 ||
		  map[iyv*mw+ixv].c == 0.0 ||
		  map[(iyv+1)*mw+ixv].c == 0.0 ||
		  map[iyv*mw+ixv+1].c == 0.0 ||
		  map[(iyv+1)*mw+ixv+1].c == 0.0)
		continue;
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
	      rx = rx * spacing;
	      ry = ry * spacing;

	      if (rx < iMinX)
		iMinX = rx;
	      if (rx > iMaxX)
		iMaxX = rx;
	      if (ry < iMinY)
		iMinY = ry;
	      if (ry > iMaxY)
		iMaxY = ry;
	    }
	}
      images[i].mapBytes = mw * mh * (sizeof(MapElement) + sizeof(InverseMapElement) +
				      sizeof(unsigned char)) +
	5 * 5 * sizeof(long long) +
	sizeof(InverseMap) +
	256;
      images[i].minX = iMinX;
      images[i].maxX = iMaxX;
      images[i].minY = iMinY;
      images[i].maxY = iMaxY;
      
      if (iMinX < oMinX)
	oMinX = iMinX;
      if (iMaxX > oMaxX)
	oMaxX = iMaxX;
      if (iMinY < oMinY)
	oMinY = iMinY;
      if (iMaxY > oMaxY)
	oMaxY = iMaxY;

      if ((nProcessed % 50) == 0 && nProcessed != 0)
	printf(" %d\n    ", nProcessed);
      printf(".");
      fflush(stdout);
      ++nProcessed;
    }
  if (map != NULL)
    free(map);
  map = NULL;

  printf("\nAll images previewed.\n");
  oWidth = oMaxX - oMinX + 1;
  oHeight = oMaxY - oMinY + 1;
  printf("output width = %lu (%d to %d) output height = %lu (%d to %d)\n\n",
	 oWidth, oMinX, oMaxX,
	 oHeight, oMinY, oMaxY);

  if (tileWidth > 0)
    {
      cols = (oWidth + tileWidth - 1) / tileWidth;
      tw = tileWidth;
    }
  else
    {
      cols = 1;
      tw = oWidth;
    }
  if (tileHeight > 0)
    {
      rows = (oHeight + tileHeight - 1) / tileHeight;
      th = tileHeight;
    }
  else
    {
      rows = 1;
      th = oHeight;
    }
  if (tw * th >= ((size_t) memoryLimit) * 1000000)
    Error("Output image too large for memory: use smaller tiles\n");
  if (tileWidth > 0 || tileHeight > 0)
    printf("Tiling output into %d rows and %d columns\n",
	   rows, cols);
  out = (unsigned char *) malloc(tw * th);

  PartitionTiles(0, cols-1, 0, rows-1, TOP, RIGHT);

  printf("\nAll tiles completed.\n");

  printf("\nWriting size file... ");
  fflush(stdout);

  sprintf(fn, "%ssize", outputName);
  f = fopen(fn, "w");
  if (f == NULL)
    Error("Could not open size file %s for writing.\n", fn);
  fprintf(f, "%d %d\n", rows, cols);
  fclose(f);
  printf("done.\n");
  fflush(stdout);

  Log("DONE!\n");
  fclose(logFile);
  return(0);
}

void
PartitionTiles (int startCol, int endCol, int startRow, int endRow,
		int entryDir, int exitDir)
{
  int i;
  size_t spaceRequired;
  int minX, maxX, minY, maxY;
  size_t regionWidth, regionHeight;
  int breakLeftRight;
  int n;

  /* check how much space would be required to compute
     these columns (startCol-endCol, inclusive) */
  spaceRequired = (endCol - startCol + 1) * (endRow - startRow + 1) *
    ((size_t) tw) * ((size_t) th) *
    (sizeof(CanvasElement) + 1);
  minX = startCol * tw + oMinX;
  maxX = (endCol + 1) * tw - 1 + oMinX;
  minY = startRow * th + oMinY;
  maxY = (endRow + 1) * th - 1 + oMinY;
  for (i = 0; i < nImages; ++i)
    {
      if (images[i].minX > maxX ||
	  images[i].maxX < minX ||
	  images[i].minY > maxY ||
	  images[i].maxY < minY)
	continue;
      spaceRequired += images[i].mapBytes;
      spaceRequired += images[i].width * images[i].height + 32;
      spaceRequired += ((images[i].width + 7) / 8) * images[i].height + 32;
      spaceRequired += images[i].width * images[i].height * sizeof(float) + 32;
      spaceRequired += tw * th;
    }
  //  printf("Space required for startCol %d endCol %d = %lu\n",
  //	 startCol, endCol, spaceRequired);
  if (spaceRequired > ((size_t) memoryLimit) * 1000000)
    {
      if (startCol == endCol && startRow == endRow)
	{
	  if (tileWidth < 0)
	    Error("Output image too large for memory: use -tile option\n");
	  else
	    Error("Insufficient memory to construct tile c%0.2dr%0.2d\n",
		  startCol, startRow);
	}

      /* break into two parts */
      regionWidth = (endCol - startCol + 1) * tw;
      regionHeight = (endRow - startRow + 1) * th;
      
      if (regionWidth >= regionHeight)
	breakLeftRight = endCol > startCol;
      else
	breakLeftRight = endRow == startRow;
      
      if (breakLeftRight)
	{
	  /* break into left and right subparts */
	  n = (endCol - startCol + 1) / 2;
	  if (entryDir != RIGHT)
	    {
	      PartitionTiles(startCol, startCol + n - 1, startRow, endRow,
			    entryDir, RIGHT);
	      PartitionTiles(startCol + n, endCol, startRow, endRow,
			    LEFT, exitDir);
	    }
	  else
	    {
	      PartitionTiles(startCol + n, endCol, startRow, endRow,
			    RIGHT, LEFT);
	      PartitionTiles(startCol, startCol + n - 1, startRow, endRow,
			    RIGHT, exitDir);
	    }
	}
      else
	{
	  /* break into top and bottom subparts */
	  n = (endRow - startRow + 1) / 2;
	  if (entryDir != BOTTOM)
	    {
	      PartitionTiles(startCol, endCol, startRow, startRow + n - 1,
			    entryDir, BOTTOM);
	      PartitionTiles(startCol, endCol, startRow + n, endRow,
			    TOP, exitDir);
	    }
	  else
	    {
	      PartitionTiles(startCol, endCol, startRow + n, endRow,
			    BOTTOM, TOP);
	      PartitionTiles(startCol, endCol, startRow, startRow + n - 1,
			    BOTTOM, exitDir);
	    }
	}
      return;
    }

  GenerateTiles(startCol, endCol, startRow, endRow);
}

void
GenerateTiles (int startCol, int endCol, int startRow, int endRow)
{
  int minX, maxX, minY, maxY;
  int nProcessed;
  int i, j;
  int iw, ih;
  int x, y;
  InverseMap *imap;
  float xv, yv;
  int ixv, iyv;
  float rx, ry;
  float rrx, rry;
  float rv;
  size_t maskBytes;
  int mbpl;
  unsigned char *mask;
  float *dist;
  float dst;
  float cx, cy;
  float offset;
  float r00, r01, r10, r11;
  int nx, ny;
  int pMinX, pMaxX, pMinY, pMaxY;
  unsigned char *image;
  float weight;
  float d, d2;
  CanvasElement *ce;
  int row, col;
  int ix, iy;
  int mFactor;
  int offsetX, offsetY;
  int v;
  size_t newCanvasSize;
  int mxMin, myMin;
  char imName0[PATH_MAX];
  char imName1[PATH_MAX];
  char msg[PATH_MAX+256];
  char fn[PATH_MAX];

  nProcessed = 0;
  printf("Generating columns %d - %d, rows %d - %d:\n    ",
	 startCol, endCol, startRow, endRow);
  fflush(stdout);

  /* free up any maps and images not used in this area */
  minX = startCol * tw + oMinX;
  maxX = (endCol + 1) * tw - 1 + oMinX;
  minY = startRow * th + oMinY;
  maxY = (endRow + 1) * th - 1 + oMinY;
  for (i = 0; i < nImages; ++i)
    if (images[i].minX > maxX ||
	images[i].maxX < minX ||
	images[i].minY > maxY ||
	images[i].maxY < minY)
      {
	if (images[i].imap != NULL)
	  {
	    FreeInverseMap(images[i].imap);
	    images[i].imap = NULL;
	      }
	    if (images[i].map != NULL)
	      {
		free(images[i].map);
		images[i].map = NULL;
	      }
	    if (images[i].image != NULL)
	      {
		free(images[i].image);
		images[i].image = NULL;
	      }
	    if (images[i].mask != NULL)
	      {
		free(images[i].mask);
		images[i].mask = NULL;
	      }
	    if (images[i].dist != NULL)
	      {
		free(images[i].dist);
		images[i].dist = NULL;
	      }
	  }

  /* allocate the canvas */
      nx = maxX - minX + 1;
      ny = maxY - minY + 1;
      newCanvasSize = nx * ny * sizeof(CanvasElement);
      if (newCanvasSize != canvasSize)
	{
	  free(canvas);
	  canvasSize = newCanvasSize;
	  canvas = (CanvasElement*) malloc(canvasSize);
	}
      memset(canvas, 0, canvasSize);
      
      for (i = 0; i < nImages; ++i)
	if (images[i].minX <= maxX &&
	    images[i].maxX >= minX &&
	    images[i].minY <= maxY &&
	    images[i].maxY >= minY)
	  {
	    /* read in map if necessary */
	    if (images[i].map == NULL)
	      {
		sprintf(fn, "%s%s.map", mapsName, images[i].name);
		if (!ReadMap(fn, &(images[i].map), &(images[i].mLevel),
			     &(images[i].mw), &(images[i].mh),
			     &mxMin, &myMin,
			     imName0, imName1,
			     msg))
		  Error("Could not read map %s:\n  error: %s\n",
			fn, msg);
		if (mxMin != 0 || myMin != 0)
		  Error("Partial maps (%s) not currently supported by applymap3\n",
			fn);
	      }
	    if (images[i].imap == NULL)
	      images[i].imap = InvertMap(images[i].map,
					 images[i].mw,
					 images[i].mh);

	    /* read in image if necessary */
	    if (images[i].image == NULL)
	      {
		sprintf(fn, "%s%s.tif", imagesName, images[i].name);
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
	      }

	    /* read in mask if necessary */
	    if (images[i].mask == NULL)
	      if (masksName[0] != '\0')
		{
		  sprintf(fn, "%s%s.pbm", masksName, images[i].name);
		  if (!ReadBitmap(fn, &(images[i].mask),
				  &iw, &ih,
				  -1, -1, -1, -1,
				  msg))
		    Error("Could not read image %s:\n  error: %s\n",
			  fn, msg);
		  if (iw != images[i].width ||
		      ih != images[i].height)
		    Error("Dimensions of mask %s do not match those in images list.\n",
			  fn);
		  mask = images[i].mask;
		  maskBytes = ((iw + 7) / 8) * ih;
		  for (j = 0; j < maskBytes; ++j)
		    mask[j] ^= 0xff;
		}
	      else
		{
		  maskBytes = ((images[i].width + 7) / 8) * images[i].height;
		  mask = images[i].mask = (unsigned char *) malloc(maskBytes);
		  memset(images[i].mask, 0, maskBytes);
		}
		
	    if (images[i].dist == NULL)
	      {
		iw = images[i].width;
		ih = images[i].height;
		images[i].dist = (float*) malloc(iw * ih * sizeof(float));
		computeDistance(EUCLIDEAN_DISTANCE,
				iw, ih, images[i].mask,
				images[i].dist);

		/* if distance from edge is less, use that */
		dist = images[i].dist;
		for (y = 0; y < ih; ++y)
		  for (x = 0; x < iw; ++x)
		    {
		      dst = dist[y*iw+x];
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
		      dist[y*iw + x] = dst;
		    }
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

	    mFactor = 1 << images[i].mLevel;
	    imap = images[i].imap;
	    image = images[i].image;
	    mask = images[i].mask;
	    dist = images[i].dist;
	    iw = images[i].width;
	    ih = images[i].height;
	    cx = images[i].width / 2;
	    cy = images[i].height / 2;
	    mbpl = (iw + 7) / 8;
	    offset = -1000000.0;
	    for (y = pMinY; y <= pMaxY; ++y)
	      for (x = pMinX; x <= pMaxX; ++x)
		{
		  if (!Invert(imap, &xv, &yv, (x + 0.5) / mFactor, (y + 0.5) / mFactor))
		    continue;
		  xv *= mFactor;
		  yv *= mFactor;
		  ixv = (int) floor(xv);
		  iyv = (int) floor(yv);
		  rrx = xv - ixv;
		  rry = yv - iyv;
		  if (ixv < 0 || ixv >= iw-1 ||
		      iyv < 0 || iyv >= ih-1)
		    continue;
		  r00 = image[iyv * iw + ixv];
		  r01 = image[(iyv + 1) * iw + ixv];
		  r10 = image[iyv * iw + ixv + 1];
		  r11 = image[(iyv + 1) * iw + ixv + 1];
		  rv = r00 * (rrx - 1.0) * (rry - 1.0)
		    - r10 * rrx * (rry - 1.0) 
		    - r01 * (rrx - 1.0) * rry
		    + r11 * rrx * rry;

		  ce = &(canvas[(y - minY) * nx + x - minX]);
		  if (blend)
		    {
		      //		      if (ixv == 10 && iyv == 10 && dist[iyv * iw + ixv] != 0.0 &&
		      //			  strstr(images[i].name, "q0") != NULL)
		      //			abort();
		      if (dist[iyv * iw + ixv] >= margin)
			weight = 1.0;
		      else
			weight = dist[iyv * iw + ixv] / margin;
		      ce->value += weight * rv;
		      ce->weight += weight;
		    }
		  else
		    {
		      if (mask[iyv * mbpl + (ixv >> 3)] & (0x80 >> (ixv & 7)))
			continue;
		      d2 = sqrt((ixv - cx) * (ixv - cx) + (iyv - cy) * (iyv - cy));
		      if (d2 + offset < ce->weight)
			{
			  ce->value = rv;
			  ce->weight = d2 + offset;
			}
		    }
		}

	    if ((nProcessed % 50) == 0 && nProcessed != 0)
	      printf(" %d\n    ", nProcessed);
	    printf(".");
	    fflush(stdout);
	    ++nProcessed;
	  }
      
      printf("\nWriting output files... ");
      fflush(stdout);

      for (col = startCol; col <= endCol; ++col)
	for (row = startRow; row <= endRow; ++row)
	  {
	    offsetX = (col - startCol) * tw;
	    offsetY = (row - startRow) * th;
	    for (y = 0; y < th; ++y)
	      {
		iy = y + offsetY;
		for (x = 0; x < tw; ++x)
		  {
		    ix = x + offsetX;
		    if (ix >= oWidth || iy >= oHeight)
		      {
			out[y*tw+x] = 0;
			continue;
		      }
		    ce = &(canvas[iy * nx + ix]);
		    if (blend)
		      {
			if (ce->weight == 0.0)
			  v = 0;
			else
			  v = (int) floor(255.999 *
					  (ce->value / ce->weight - blackValue) /
					  range);
		      }
		    else
			v = (int) floor(255.999 * (ce->value - blackValue) / range);
		    if (v < 0)
		      v = 0;
		    else if (v > 255)
		      v = 255;
		    out[y*tw+x] = v;
		  }
	      }
	    if (tileWidth < 0 && tileHeight < 0)
	      strcpy(fn, outputName);
	    else
	      sprintf(fn, "%sc%0.2d%sr%0.2d.tif", outputName, col,
		      tree ? "/" : "", row);
	    if (!CreateDirectories(fn))
	      Error("Could not create directories for output file %s\n", fn);
	    if (!WriteImage(fn, out, (int) tw, (int) th,
			    compress ? HDiffDeflateImage : UncompressedImage,
			    msg))
	      Error("Could not write output file %s:\n  error: %s\n",
		    fn, msg);
	  }

      printf("done.\n\n");
      fflush(stdout);
    }



void Log (char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);
  vfprintf(logFile, fmt, args);
  va_end(args);
  fflush(logFile);
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

