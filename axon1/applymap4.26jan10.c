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
  MapElement *map;      /* map of this image into the final image */
  InverseMap *invMap;     /* inverse map that translates points in the final image
			   into points in this image */
  /* intensity map is optional */
  int imapLevel;	/* intensity map level */
  int imapw, imaph;	/* intensity map width and height */
  MapElement *imap;	/* intensity map */
} Image;

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
int overlay = 0;
int blend = 1;
int margin = -1;
int tileWidth = -1;
int tileHeight = -1;
int memoryLimit = 1024;   /* 1 GB */
int tree = 0;
float blackValue = 0.0;
float whiteValue = 255.0;
int compress = 0;
int reductionFactor = 1;

int nImages = 0;
Image *images = 0;
int *imageHashTable = 0;
unsigned char *canvas = 0;
unsigned short *weight = 0;
size_t canvasSize = 0;
int canvasWidth, canvasHeight;
int oldCanvasWidth;
int canvasMinX, canvasMinY;
int weightWidth, weightHeight;
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


#define DIR_HASH_SIZE	8192
char *dirHash[DIR_HASH_SIZE];

/* FORWARD DECLARATIONS */
void PaintImage (int i, int minX, int maxX, int minY, int maxY);
void WriteTiles (int col, int startRow, int endRow, char *iName);
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
  float minX, minY, maxX, maxY;
  int iMinX, iMinY, iMaxX, iMaxY;
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

  // FIX -- TEMPORARY HACK
  //  CPU_ZERO(&cpumask);
  //  CPU_SET(7, &cpumask);
  //  sched_setaffinity(0, sizeof(cpu_set_t), &cpumask);

  sprintf(fn, "logs/applymap4.log");
  logFile = fopen(fn, "w");
  error = 0;
  imageListName[0] = '\0';
  imagesName[0] = '\0';
  masksName[0] = '\0';
  mapsName[0] = '\0';
  imapsName[0] = '\0';
  outputName[0] = '\0';
  strcpy(extension, "tif");
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
    else if (strcmp(argv[i], "-reduction") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d", &reductionFactor) != 1)
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
      fprintf(stderr, "              [-region WxH+X+Y]\n");
      fprintf(stderr, "              [-reduction reduction_factor]\n");
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
      minX = 1000000000;
      maxX = -1000000000;
      minY = 1000000000;
      maxY = -1000000000;
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

      if ((nProcessed % 50) == 0 && nProcessed != 0)
	printf(" %d\n    ", nProcessed);
      printf(".");
      fflush(stdout);
      ++nProcessed;
    }

  printf("\nAll images previewed.\n");
  oWidth = oMaxX - oMinX + 1;
  oHeight = oMaxY - oMinY + 1;
  printf("output width = %lu (%d to %d) output height = %lu (%d to %d)\n\n",
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
  oMaxX = oMinX + oWidth - 1;
  oMaxY = oMinY + oHeight - 1;

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

  for (oi = 0; oi < nOutputImages; ++oi)
    {
      if (overlay)
	{
	  startImage = 0;
	  endImage = nImages-1;
	}
      else
	{
	  startImage = oi;
	  endImage = oi;
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
		  outHeight = ((oHeight / hs) / (th * reductionFactor)) *
		    (th * reductionFactor);
		  outWidth = tw;
		  if (outHeight == 0)
		    Error("Insufficent memory to construct a single row of tiles.\n");
		  endY = startY + outHeight - 1;
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
	      if (hi == 0)
		out = (unsigned char *) malloc(outHeight * outWidth *
					       sizeof(unsigned char));
	      outMinY = oMinY;
	      startRow = 0;
	      endRow = rows-1;
	    }
	  else
	    {
	      outHeight = ((oHeight / hs) / (th * reductionFactor)) *
		(th * reductionFactor);
	      outWidth = tw;
	      if (outHeight == 0)
		Error("Insufficent memory to construct a single row of tiles.\n");
	      endY = startY + outHeight - 1;
	      out = (unsigned char *) malloc(outHeight * outWidth *
					     sizeof(unsigned char));
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

	  startX = oMinX;
	  do {
	    // see how far we can increase endX while still satisfying
	    //   memory constraints
	    endX = startX-1;
	    memoryRequired = imageMem;
	    memoryRequired += outHeight * outWidth * sizeof(unsigned char);
	    printf("memreq1 = %llu\n", memoryRequired);
	    memoryRequired += 2 * oWidth * sizeof(size_t);
	    printf("memreq2 = %llu\n", memoryRequired);
	    memoryRequired += (startX - canvasMinX) * canvasHeight *
	      sizeof(unsigned char);
	    printf("memreq3 = %llu\n", memoryRequired);
	    if (!overlay)
	      memoryRequired += images[oi].width * images[oi].height *
		(sizeof(unsigned char) + sizeof(unsigned char)) +
		((images[oi].width + 7) / 8) * images[oi].height +
		images[oi].mapBytes;
	    printf("memreq4 = %llu   (%d %d %d %d %d %d)\n", memoryRequired,
		   oi, images[oi].width, images[oi].height, images[oi].mapBytes,
		   startImage, endImage);
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
	      Error("Internal error: not enough memory to construct a single pixel column\n%d %d %d %d %llu",
		    startX, endX, startY, endY, memoryRequired);
	    if (endX < oMaxX &&
		endX - canvasMinX + 1 < reductionFactor)
	      Error("Internal error: not enough memory to construct a column of the output image\n");
	    printf("Rendering area from x=%d to x=%d, y=%d to y=%d\n",
		   startX, endX, startY, endY);

	    // resize the canvas
	    oldCanvasWidth = canvasWidth;
	    canvasWidth = endX - canvasMinX + 1;
	    canvas = (unsigned char *) realloc(canvas, canvasWidth * canvasHeight * sizeof(unsigned char));
	    weightWidth = endX - startX + 1;
	    weightHeight = canvasHeight;
	    weight = (unsigned short *) malloc(weightWidth * weightHeight * sizeof(unsigned short));
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
		nx = (canvasWidth - cx) / reductionFactor;
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

	    // save extra pixel columns
	    canvasWidth -= cx;
	    canvasMinX += cx;
	    if (canvasWidth > 0)
	      for (y = 0; y < canvasHeight; ++y)
		memmove(&canvas[y*canvasWidth],
			&canvas[y*oldCanvasWidth + cx],
			canvasWidth * sizeof(unsigned char));
	    canvas = (unsigned char *)
	      realloc(canvas,
		      canvasWidth * canvasHeight * sizeof(unsigned char));
	    weightMinX += weightWidth;

	    free(weight);

	    startX = endX + 1;
	  } while (startX <= oMaxX);

	  if (cols > 1)
	    free(out);
	}
      if (cols == 1)
	free(out);
    }

  printf("\nWriting size file... ");
  fflush(stdout);

  if (strlen(outputName) != 0)
    sprintf(fn, "%s.size", outputName);
  else
    strcpy(fn, "size");
  f = fopen(fn, "w");
  if (f == NULL)
    Error("Could not open size file %s for writing.\n", fn);
  fprintf(f, "%d %d\n%dx%d%+d%+d\n%dx%d\n",
	  rows, cols,
	  oWidth, oHeight, oMinX, oMinY,
	  oWidth / reductionFactor, oHeight / reductionFactor);
  fclose(f);
  printf("done.\n");
  fflush(stdout);

  Log("DONE!\n");
  fclose(logFile);
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
  int mbpl;
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
  int iixv, iiyv;
  float rb00, rb01, rb10, rb11;
  float rw00, rw01, rw10, rw11;
  float rb, rw;
  MapElement *imap;
  int imapFactor;
  int imapw, imaph;
  int imapXMin, imapYMin;
  char imapName0[PATH_MAX], imapName1[PATH_MAX];
  float xvi, yvi;

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
      imageMem += images[i].mapBytes;
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
      imageMem += images[i].width * images[i].height;
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
	imageMem += maskBytes;
      }
    else
      {
	maskBytes = ((images[i].width + 7) / 8) * images[i].height;
	mask = images[i].mask = (unsigned char *) malloc(maskBytes);
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
	imap = images[i].imap;
	imageMem += images[i].imapw * images[i].imaph * sizeof(MapElement);
      }
		
  /* compute distance table if necessary */
  if (images[i].dist == NULL)
    {
      iw = images[i].width;
      ih = images[i].height;
      images[i].dist = (unsigned char*) malloc(iw * ih * sizeof(unsigned char));
      distance = (float*) malloc(iw * ih * sizeof(float));
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
	  }
      }

  /* free up if no longer required */
  if (images[i].maxX <= maxX || maxX >= oMaxX)
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
	  imageMem -= images[i].mapBytes; // this accounts for the inverse
	                                  // map as well
	}
      if (images[i].imap != NULL)
	{
	  free(images[i].imap);
	  images[i].imap = NULL;
	  imageMem -= images[i].imapw * images[i].imaph * sizeof(MapElement);
	}
      if (images[i].image != NULL)
	{
	  free(images[i].image);
	  images[i].image = NULL;
	  imageMem -= images[i].width * images[i].height;
	}
      if (images[i].mask != NULL)
	{
	  free(images[i].mask);
	  images[i].mask = NULL;
	  imageMem -= ((images[i].width + 7) / 8) * images[i].height;
	}
      if (images[i].dist != NULL)
	{
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
	    sprintf(fn, "%sc%0.2d%sr%0.2d.tif", outputName, col+1,
		    tree ? "/" : "", row+1);
	}
      else
	{
	  if (tileWidth < 0 && tileHeight < 0)
	    sprintf(fn, "%s%s.tif",
		    outputName, iName);
	  else
	    sprintf(fn, "%s%s/c%0.2d%sr%0.2d.tif",
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
