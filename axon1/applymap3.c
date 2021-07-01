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
  InverseMap *invMap;     /* inverse map that translates points in the final image
			   into points in this image */
  /* intensity map is optional */
  int imapLevel;	/* intensity map level */
  int imapw, imaph;	/* intensity map width and height */
  MapElement *imap;	/* intensity map */
} Image;

typedef struct CanvasElement
{
  float value;
  float weight;
} CanvasElement;


/* GLOBAL VARIABLES */
FILE *logFile = 0;

int resume = 0;
char imageListName[PATH_MAX];
char imagesName[PATH_MAX];
char extension[PATH_MAX];
char masksName[PATH_MAX];
char mapsName[PATH_MAX];
char imapsName[PATH_MAX];
char outputName[PATH_MAX];
char outputExtension[PATH_MAX];
int overlay = 0;
int blend = 0;
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
void PartitionTiles (int startImage, int endImage,
		     int startCol, int endCol,
		     int startRow, int endRow,
		     int entryDir, int exitDir);
void GenerateTiles (int startImage, int endImage,
		    int startCol, int endCol,
		    int startRow, int endRow);
void Cleanup (int startImage, int endImage);
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
  int irx, iry;
  cpu_set_t cpumask;
  int nOutputImages;
  int startImage, endImage;
  int imapLevel;
  int imapw, imaph;
  int imapXMin, imapYMin;
  char imapName0[PATH_MAX], imapName1[PATH_MAX];
  MapElement *imap;
  int regionWidth, regionHeight, regionOffsetX, regionOffsetY;

  // FIX -- TEMPORARY HACK
  //  CPU_ZERO(&cpumask);
  //  CPU_SET(7, &cpumask);
  //  sched_setaffinity(0, sizeof(cpu_set_t), &cpumask);

  sprintf(fn, "logs/applymap3.log");
  logFile = fopen(fn, "w");
  error = 0;
  imageListName[0] = '\0';
  imagesName[0] = '\0';
  masksName[0] = '\0';
  mapsName[0] = '\0';
  imapsName[0] = '\0';
  outputName[0] = '\0';
  strcpy(extension, "tif");
  strcpy(outputExtension, "tif");
  regionWidth = regionHeight = regionOffsetX = regionOffsetY = -1;
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
    else if (strcmp(argv[i], "-pgm") == 0)
      strcpy(extension, "pgm");
    else if (strcmp(argv[i], "-output_pgm") == 0)
      strcpy(outputExtension, "pgm");
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
    else if (strcmp(argv[i], "-imaps") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(imapsName, argv[i]);
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
    else error = 1;

  if (error)
    {
      fprintf(stderr, "Usage: applymap3 -image_list list_file -images image_prefix\n");
      fprintf(stderr, "              -maps map_prefix -output file_prefix\n");
      fprintf(stderr, "              [-masks mask_prefix]\n");
      fprintf(stderr, "              [-imaps imaps_prefix]\n");
      fprintf(stderr, "              [-blend]\n");
      fprintf(stderr, "              [-mosaic]\n");
      fprintf(stderr, "              [-margin margin_thickness_in_pixels]\n");
      fprintf(stderr, "              [-tile tilewidthxtileheight]\n");
      fprintf(stderr, "              [-memory memory_limit_in_MB]\n");
      fprintf(stderr, "              [-tree]\n");
      fprintf(stderr, "              [-black black_value]\n");
      fprintf(stderr, "              [-white white_value]\n");
      fprintf(stderr, "              [-compress]\n");
      fprintf(stderr, "              [-resume]\n");
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
      images[nImages].invMap = NULL;
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
  imap = NULL;
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

	  //	  printf("iw = %d ih = %d n = %d\n",
	  //		images[i].width, images[i].height, n);
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
	      /* extrapolate as necessary */
	      while (ixv < 0)
		{
		  ++ixv;
		  rrx -= 1.0;
		}
	      while (ixv >= mw-1)
		{
		  --ixv;
		  rrx += 1.0;
		}
	      while (iyv < 0)
		{
		  ++iyv;
		  rry -= 1.0;
		}
	      while (iyv >= mh-1)
		{
		  --iyv;
		  rry += 1.0;
		}
	      if (map[iyv*mw+ixv].c == 0.0 ||
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

	      irx = (int) floor(rx + 0.01);
	      iry = (int) floor(ry + 0.01);
	      if (irx < iMinX)
		iMinX = irx;
	      if (iry < iMinY)
		iMinY = iry;
	      irx = (int) floor(rx - 0.01);
	      iry = (int) floor(ry - 0.01);
	      if (irx > iMaxX)
		iMaxX = irx;
	      if (iry > iMaxY)
		iMaxY = iry;
	      //	      printf("dir = %d j = %d xv = %f yv = %f rx = %f ry = %f irx = %d iry = %d\n",
	      //		     dir, j, xv, yv, rx, ry, irx, iry);
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
	}

      if ((nProcessed % 50) == 0 && nProcessed != 0)
	printf(" %d\n    ", nProcessed);
      printf(".");
      fflush(stdout);
      ++nProcessed;
    }
  if (map != NULL)
    free(map);
  map = NULL;
  if (imap != NULL)
    free(imap);
  imap = NULL;

  printf("\nAll images previewed.\n");
  oWidth = oMaxX - oMinX + 1;
  oHeight = oMaxY - oMinY + 1;
  if (regionWidth > 0)
    {
      oWidth = regionWidth;
      oHeight = regionHeight;
      oMinX = regionOffsetX;
      oMinY = regionOffsetY;
      oMaxX = oMinX + regionWidth - 1;
      oMaxY = oMinY + regionHeight - 1;
    }
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

  if (overlay)
    nOutputImages = 1;
  else
    nOutputImages = nImages;

  for (i = 0; i < nOutputImages; ++i)
    {
      if (overlay)
	{
	  startImage = 0;
	  endImage = nImages-1;
	}
      else
	{
	  startImage = i;
	  endImage = i;
	}
      PartitionTiles(startImage, endImage, 0, cols-1, 0, rows-1, TOP, RIGHT);
      Cleanup(startImage, endImage);
      printf("\nAll tiles of output image %d completed.\n", i);
    }

  printf("\nWriting size file... ");
  fflush(stdout);


  sprintf(fn, "%ssize", outputName);
  f = fopen(fn, "w");
  if (f == NULL)
    Error("Could not open size file %s for writing.\n", fn);
  fprintf(f, "%d %d\n%dx%d%+d%+d\n",
	  rows, cols,
	  oWidth, oHeight, oMinX, oMinY);
  fclose(f);
  printf("done.\n");
  fflush(stdout);

  Log("DONE!\n");
  fclose(logFile);
  return(0);
}

void
PartitionTiles (int startImage, int endImage,
		int startCol, int endCol,
		int startRow, int endRow,
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
  for (i = startImage; i <= endImage; ++i)
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
	      PartitionTiles(startImage, endImage,
			     startCol, startCol + n - 1,
			     startRow, endRow,
			     entryDir, RIGHT);
	      PartitionTiles(startImage, endImage,
			     startCol + n, endCol,
			     startRow, endRow,
			     LEFT, exitDir);
	    }
	  else
	    {
	      PartitionTiles(startImage, endImage,
			     startCol + n, endCol,
			     startRow, endRow,
			     RIGHT, LEFT);
	      PartitionTiles(startImage, endImage,
			     startCol, startCol + n - 1,
			     startRow, endRow,
			     RIGHT, exitDir);
	    }
	}
      else
	{
	  /* break into top and bottom subparts */
	  n = (endRow - startRow + 1) / 2;
	  if (entryDir != BOTTOM)
	    {
	      PartitionTiles(startImage, endImage,
			     startCol, endCol,
			     startRow, startRow + n - 1,
			    entryDir, BOTTOM);
	      PartitionTiles(startImage, endImage,
			     startCol, endCol,
			     startRow + n, endRow,
			    TOP, exitDir);
	    }
	  else
	    {
	      PartitionTiles(startImage, endImage,
			     startCol, endCol,
			     startRow + n, endRow,
			    BOTTOM, TOP);
	      PartitionTiles(startImage, endImage,
			     startCol, endCol,
			     startRow, startRow + n - 1,
			    BOTTOM, exitDir);
	    }
	}
      return;
    }

  GenerateTiles(startImage, endImage, startCol, endCol, startRow, endRow);
}

void
GenerateTiles (int startImage, int endImage,
	       int startCol, int endCol,
	       int startRow, int endRow)
{
  int minX, maxX, minY, maxY;
  int nProcessed;
  int i, j;
  int iw, ih;
  int x, y;
  InverseMap *invMap;
  float xv, yv;
  int ixv, iyv;
  float rx, ry;
  float rrx, rry;
  float rv, dv;
  size_t maskBytes;
  int mbpl;
  unsigned char *mask;
  float *dist;
  float dst;
  float cx, cy;
  float offset;
  float r00, r01, r10, r11;
  float d00, d01, d10, d11;
  int nx, ny;
  int pMinX, pMaxX, pMinY, pMaxY;
  unsigned char *image;
  float weight;
  float d;
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
  struct stat statBuf;
  int complete;
  float xvi, yvi;
  int iixv, iiyv;
  float rb00, rb01, rb10, rb11;
  float rw00, rw01, rw10, rw11;
  float rb, rw;
  MapElement *imap;
  int imapFactor;
  int imapw, imaph;
  int imapXMin, imapYMin;
  char imapName0[PATH_MAX], imapName1[PATH_MAX];

  if (resume)
    {
      complete = 1;
      for (col = startCol; col <= endCol && complete; ++col)
	for (row = startRow; row <= endRow && complete; ++row)
	  {
	    if (tileWidth < 0 && tileHeight < 0)
	      strcpy(fn, outputName);
	    else
	      sprintf(fn, "%sc%0.2d%sr%0.2d.tif", outputName, col+1,
		      tree ? "/" : "", row+1);
	    if (stat(fn, &statBuf) != 0 || !S_ISREG(statBuf.st_mode))
	      complete = 0;
	  }
      if (complete)
	return;
    }

  nProcessed = 0;
  printf("Generating columns %d - %d, rows %d - %d:\n    ",
	 startCol+1, endCol+1, startRow+1, endRow+1);
  fflush(stdout);

  /* free up any maps and images not used in this area */
  minX = startCol * tw + oMinX;
  maxX = (endCol + 1) * tw - 1 + oMinX;
  minY = startRow * th + oMinY;
  maxY = (endRow + 1) * th - 1 + oMinY;
  for (i = startImage; i <= endImage; ++i)
    if (images[i].minX > maxX ||
	images[i].maxX < minX ||
	images[i].minY > maxY ||
	images[i].maxY < minY)
      {
	if (images[i].invMap != NULL)
	  {
	    FreeInverseMap(images[i].invMap);
	    images[i].invMap = NULL;
	  }
	if (images[i].map != NULL)
	  {
	    free(images[i].map);
	    images[i].map = NULL;
	  }
	if (images[i].imap != NULL)
	  {
	    free(images[i].imap);
	    images[i].imap = NULL;
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
      if (canvas != NULL)
	free(canvas);
      canvasSize = newCanvasSize;
      canvas = (CanvasElement*) malloc(canvasSize);
    }
  memset(canvas, 0, canvasSize);
      
  for (i = startImage; i <= endImage; ++i)
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
	if (images[i].invMap == NULL)
	  images[i].invMap = InvertMap(images[i].map,
				       images[i].mw,
				       images[i].mh);
	
	/* read in image if necessary */
	if (images[i].image == NULL)
	  {
	    sprintf(fn, "%s%s.%s", imagesName, images[i].name, extension);
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
	      imap = images[i].imap;
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
	invMap = images[i].invMap;
	image = images[i].image;
	mask = images[i].mask;
	imap = images[i].imap;
	if (imap != NULL)
	  {
	    imapFactor = 1 << images[i].imapLevel;
	    imapw = images[i].imapw;
	    imaph = images[i].imaph;
	  }
	dist = images[i].dist;
	iw = images[i].width;
	ih = images[i].height;
	cx = (iw - 1) / 2.0;
	cy = (ih - 1) / 2.0;
	mbpl = (iw + 7) / 8;
	offset = -1000000.0;
	for (y = pMinY; y <= pMaxY; ++y)
	  for (x = pMinX; x <= pMaxX; ++x)
	    {
	      if (!Invert(invMap, &xv, &yv, (x + 0.5) / mFactor, (y + 0.5) / mFactor))
		continue;
	      xv = xv * mFactor - 0.5;
	      yv = yv * mFactor - 0.5;
	      ixv = (int) floor(xv);
	      iyv = (int) floor(yv);
	      rrx = xv - ixv;
	      rry = yv - iyv;
	      //	      if (y == pMinY || x == pMinX)
	      //		printf("Sampling: x = %d y = %d xv = %f yv = %f ixv = %d iyv = %d rrx = %f rry = %f\n",
	      //		       x, y, xv, yv, ixv, iyv, rrx, rry);
	      if (ixv < -1 || ixv >= iw ||
		  iyv < -1 || iyv >= ih)
		continue;
	      if (ixv >= 0 && iyv >= 0 &&
		  !(mask[iyv * mbpl + (ixv >> 3)] & (0x80 >> (ixv & 7))))
		{
		  r00 = image[iyv * iw + ixv];
		  d00 = dist[iyv * iw + ixv];
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
		    Error("rw level (%f) is less than rb level (%f) for image %s at (%d %d)\n",
			  rw, rb, images[i].name, ixv, iyv);
		  rv = (rv - rb) / (rw - rb);
		}

	      ce = &(canvas[(y - minY) * nx + x - minX]);
	      if (blend)
		{
		  if (dv >= margin)
		    weight = 1.0;
		  else
		    weight = dv / margin;
		  ce->value += weight * rv;
		  ce->weight += weight;
		}
	      else
		{
		  if (dv < 1.0)
		    weight = -dv;
		  else
		    weight = hypotf(xv - cx, yv - cy) + offset;
		  if (weight < ce->weight)
		    {
		      ce->value = rv;
		      ce->weight = weight;
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
		if (imapsName[0] != '\0')
		  {
		    if (blend)
		      {
			if (ce->weight == 0.0)
			  v = 0;
			else
			  v = (int) floor(255.0 * ce->value / ce->weight + 0.5);
		      }
		    else
		      v = (int) floor(255.0 * ce->value + 0.5);
		  }
		else
		  {
		    if (blend)
		      {
			if (ce->weight == 0.0)
			  v = 0;
			else
			  v = (int) floor(255.0 *
					  (ce->value / ce->weight - blackValue) /
					  range + 0.5);
		      }
		    else
		      v = (int) floor(255.0 * (ce->value - blackValue) / range + 0.5);
		  }
		if (v <= 0)
		  {
		    if (ce->weight == 0.0)
		      v = 0;
		    else
		      v = 1;
		  }
		else if (v > 255)
		  v = 255;
		out[y*tw+x] = v;
	      }
	  }
	if (overlay)
	  {
	    if (tileWidth < 0 && tileHeight < 0)
	      sprintf(fn, "%s.%s", outputName, outputExtension);
	    else
	      sprintf(fn, "%sc%0.2d%sr%0.2d.%s", outputName, col+1,
		      tree ? "/" : "", row+1,
		      outputExtension);
	  }
	else
	  {
	    if (tileWidth < 0 && tileHeight < 0)
	      sprintf(fn, "%s%s.%s",
		      outputName, images[startImage].name,
		      outputExtension);
	    else
	      sprintf(fn, "%s%s/c%0.2d%sr%0.2d.%s",
		      outputName, images[startImage].name,
		      col+1, tree ? "/" : "", row+1,
		      outputExtension);
	  }
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

void
Cleanup (int startImage, int endImage)
{
  int i;

  /* free up any maps and images */
  for (i = startImage; i <= endImage; ++i)
    {
      if (images[i].invMap != NULL)
	{
	  FreeInverseMap(images[i].invMap);
	  images[i].invMap = NULL;
	}
      if (images[i].map != NULL)
	{
	  free(images[i].map);
	  images[i].map = NULL;
	}
      if (images[i].imap != NULL)
	{
	  free(images[i].imap);
	  images[i].imap = NULL;
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

