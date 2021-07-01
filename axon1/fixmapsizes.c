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

#define LINE_LENGTH	255

char imageListName[PATH_MAX];
char imagesName[PATH_MAX];
char mapsName[PATH_MAX];
char outputName[PATH_MAX];
int tileWidth = -1;
int tileHeight = -1;

int oMinX, oMinY, oMaxX, oMaxY;
size_t oWidth, oHeight;

#define TOP	0
#define BOTTOM	1
#define LEFT	2
#define RIGHT	3

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

  error = 0;
  imageListName[0] = '\0';
  imagesName[0] = '\0';
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
    else if (strcmp(argv[i], "-tile") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%dx%d",
				  &tileWidth, &tileHeight) != 2)
	  {
	    error = 1;
	    break;
	  }
      }
    else error = 1;

  if (error)
    {
      fprintf(stderr, "Usage: fixmapsizes -image_list list_file -images image_prefix\n");
      fprintf(stderr, "              -maps map_prefix -output file_prefix\n");
      fprintf(stderr, "              [-tile tilewidthxtileheight]\n");
      exit(1);
    }

  /* check that at least minimal parameters were supplied */
  if (imageListName[0] == '\0' || imagesName[0] == '\0' ||
      mapsName[0] == '\0' || outputName[0] == '\0')
    {
      fprintf(stderr, "-image_list, -images, -maps, and -output parameters must be specified.\n");
      exit(1);
    }

  printf("Previewing maps: ");
  fflush(stdout);
  oMinX = 1000000000;
  oMaxX = -1000000000;
  oMinY = 1000000000;
  oMaxY = -1000000000;
  map = NULL;
  nProcessed = 0;

  /* read the image list */
  f = fopen(imageListName, "r");
  if (f == NULL)
    {
      fprintf(stderr, "Could not open file %s for reading\n", imageListName);
      exit(1);
    }
  while (fgets(line, LINE_LENGTH, f) != NULL)
    {
      nItems = sscanf(line, "%s%d%d",
		      imageName, &width, &height);
      if (nItems != 3)
	{
	  fprintf(stderr, "Malformed line in %s:\n%s\n", imageListName, line);
	  exit(1);
	}

      sprintf(fn, "%s%s.map", mapsName, imageName);
      if (!ReadMap(fn, &map, &mLevel,
		   &mw, &mh, &mxMin, &myMin,
		   imName0, imName1,
		   msg))
	{
	  fprintf(stderr, "Could not read map %s:\n  error: %s\n",
		  fn, msg);
	  exit(1);
	}

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
	      n = width+1;
	      ax = 1.0;
	      bx = 0.0;
	      ay = 0.0;
	      by = 0.0;
	      break;

	    case 1:
	      n = width+1;
	      ax = 1.0;
	      bx = 0.0;
	      ay = 0.0;
	      by = height;
	      break;

	    case 2:
	      n = height-1;
	      ax = 0.0;
	      bx = 0.0;
	      ay = 1.0;
	      by = 1.0;
	      break;
	      
	    case 3:
	      n = height-1;
	      ax = 0.0;
	      bx = width;
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
  fclose(f);

  printf("\nAll images previewed.\n");
  oWidth = oMaxX - oMinX + 1;
  oHeight = oMaxY - oMinY + 1;
  printf("output width = %lu (%d to %d) output height = %lu (%d to %d)\n\n",
	 oWidth, oMinX, oMaxX,
	 oHeight, oMinY, oMaxY);

  if (tileWidth > 0)
    cols = (oWidth + tileWidth - 1) / tileWidth;
  else
    cols = 1;
  if (tileHeight > 0)
    rows = (oHeight + tileHeight - 1) / tileHeight;
  else
    rows = 1;
  if (tileWidth > 0 || tileHeight > 0)
    printf("Tiling output into %d rows and %d columns\n",
	   rows, cols);

  printf("\nWriting size file... ");
  fflush(stdout);

  sprintf(fn, "%ssize", outputName);
  f = fopen(fn, "w");
  if (f == NULL)
    {
      fprintf(stderr, "Could not open size file %s for writing.\n", fn);
      exit(1);
    }
  fprintf(f, "%d %d\n%dx%d%+d%+d\n",
	  rows, cols,
	  oWidth, oHeight, oMinX, oMinY);
  fclose(f);
  printf("done.\n");
  fflush(stdout);
  return(0);
}
