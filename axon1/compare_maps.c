#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <errno.h>

#include "imio.h"
#include "invert.h"

void Error (char *fmt, ...);

int
main (int argc, char **argv)
{
  int invert;
  float thresholdC;

  char map1Name[PATH_MAX];
  MapElement *map1;
  int mLevel;
  int mw, mh;
  int mxMin, myMin;
  char imgName[PATH_MAX], refName[PATH_MAX];
  
  char map2Name[PATH_MAX];
  MapElement *map2;
  int mLevel2;
  int mw2, mh2;
  int mxMin2, myMin2;
  char imgName2[PATH_MAX], refName2[PATH_MAX];
  
  int x, y;
  int i;
  int error;
  char msg[PATH_MAX+256];

  int n;
  float x1, y1, c1;
  float x2, y2, c2;
  float d, xd, yd, cd;
  float maxD, maxXD, maxYD, maxCD;
  float meanD, meanXD, meanYD, meanCD;


  error = 0;
  invert = 0;
  thresholdC = 0.0;
  for (i = 1; i < argc; ++i)
    if (strcmp(argv[i], "-map1") == 0)
      {
	if (++i == argc)
          {
            error = 1;
            break;
          }
        strcpy(map1Name, argv[i]);
      }
    else if (strcmp(argv[i], "-map2") == 0)
      {
        if (++i == argc)
          {
            error = 1;
            break;
          }
	invert = 0;
        strcpy(map2Name, argv[i]);
      }
    else
      {
	error = 1;
	break;
      }

  if (error)
    {
      if (i < argc)
        fprintf(stderr, "Invalid option: %s\n", argv[i]);
      else
        fprintf(stderr, "Incomplete option: %s\n", argv[i-1]);
      fprintf(stderr, "\n");

      fprintf(stderr, "Usage: compare_maps -map1 <map_name>\n");
      fprintf(stderr, "            [-map2 <map_name>]\n");
      exit(1);
    }

  /* check that at least minimal parameters were supplied */
  if (map1Name[0] == '\0')
    Error("-map1 parameter must be specified.\n");
  if (map2Name[0] == '\0')
    Error("-map2 or -inverse_map2 parameter must be specified.\n");

  if (!ReadMap(map1Name, &map1, &mLevel, &mw, &mh,
	       &mxMin, &myMin, imgName, refName, msg))
    Error("Could not read map %s:\n%s\n", map1Name, msg);

  if (!ReadMap(map2Name, &map2, &mLevel2, &mw2, &mh2,
	       &mxMin2, &myMin2, imgName2, refName2, msg))
    Error("Could not read map %s:\n%s\n", map2Name, msg);

  if (mLevel != mLevel2)
    {
      printf("map levels differ: %d vs. %d\n", mLevel, mLevel2);
      exit(0);
    }

  if (mw != mw2)
    {
      printf("map widths differ: %d vs %d\n", mw, mw2);
      exit(0);
    }

  if (mh != mh2)
    {
      printf("map heights differ: %d vs %d\n", mh, mh2);
      exit(0);
    }

  maxXD = 0.0;
  maxYD = 0.0;
  maxD = 0.0;
  maxCD = 0.0;
  meanXD = 0.0;
  meanYD = 0.0;
  meanD = 0.0;
  meanCD = 0.0;
  n = 0;
  for (y = 0; y < mh; ++y)
    for (x = 0; x < mw; ++x)
      {
	x1 = map1[y*mw+x].x * (1 << mLevel);
	y1 = map1[y*mw+x].y * (1 << mLevel);
	c1 = map1[y*mw+x].c;

	x2 = map2[y*mw+x].x * (1 << mLevel);
	y2 = map2[y*mw+x].y * (1 << mLevel);
	c2 = map2[y*mw+x].c;

	if (c1 == 0.0 && c2 == 0.0)
	  continue;

	++n;

	cd = c2 - c1;
	if (fabs(cd) > fabs(maxCD))
	  maxCD = cd;
	meanCD += cd;

	xd = x2 - x1;
	if (fabs(xd) > fabs(maxXD))
	  maxXD = xd;
	meanXD += xd;

	yd = y2 - y1;
	if (fabs(yd) > fabs(maxYD))
	  maxYD = yd;
	meanYD += yd;

	d = sqrt(xd * xd + yd * yd);
	if (d > maxD)
	  maxD = d;
	meanD += d;
      }

  if (n == 0)
    {
      printf("No valid map points found.\n");
      exit(0);
    }

  meanCD = meanCD / n;
  meanXD = meanXD / n;
  meanYD = meanYD / n;
  meanD = meanD / n;

  printf("maximum distance: %f\n", maxD);
  printf("average distance: %f\n", meanD);
  printf("maximum absolute delta: x: %f  y: %f  c: %f\n", maxXD, maxYD, maxCD);
  printf("average delta: x: %f  y: %f  c: %f\n", meanXD, meanYD, meanCD);

  free(map1);
  free(map2);

  return(0);
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
