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

void Error (char *fmt, ...);

int
main (int argc, char **argv)
{
  double weight1, weight2;
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
  
  char outputMapName[PATH_MAX];
  MapElement *omap;

  int x, y;
  double x1, y1, c1;
  double xv, yv;
  double xp, yp;
  int ix, iy;
  double rx, ry, rc;
  double rx00, rx01, rx10, rx11;
  double ry00, ry01, ry10, ry11;
  double rc00, rc01, rc10, rc11;
  double rrx, rry;

  int dx, dy;
  double weight;
  double totalWeight;
  double trx, try;
  double x2, y2;

  int i;
  int error;
  char msg[PATH_MAX+256];


  error = 0;
  weight1 = -1.0;
  weight2 = -1.0;
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
        strcpy(map2Name, argv[i]);
      }
    else if (strcmp(argv[i], "-weight1") == 0)
      {
        if (weight1 >= 0.0 ||
	    ++i == argc ||
	    sscanf(argv[i], "%lf", &weight1) != 1 ||
	    weight1 < 0.0)
          {
            error = 1;
            break;
          }
      }
    else if (strcmp(argv[i], "-weight2") == 0)
      {
        if (weight2 >= 0.0 ||
	    ++i == argc ||
	    sscanf(argv[i], "%lf", &weight2) != 1 ||
	    weight2 < 0.0)
          {
            error = 1;
            break;
          }
      }
    else if (strcmp(argv[i], "-inverse_weight1") == 0)
      {
        if (weight1 >= 0.0 ||
	    ++i == argc ||
	    sscanf(argv[i], "%lf", &weight1) != 1 ||
	    weight1 <= 0.0)
          {
            error = 1;
            break;
          }
	weight1 = 1.0 / weight1;
      }
    else if (strcmp(argv[i], "-inverse_weight2") == 0)
      {
        if (weight2 >= 0.0 ||
	    ++i == argc ||
	    sscanf(argv[i], "%lf", &weight2) != 1 ||
	    weight2 <= 0.0)
          {
            error = 1;
            break;
          }
	weight2 = 1.0 / weight2;
      }
    else if (strcmp(argv[i], "-output") == 0)
      {
        if (++i == argc)
          {
            error = 1;
            break;
          }
        strcpy(outputMapName, argv[i]);
      }
    else if (strcmp(argv[i], "-threshold_c") == 0)
      {
        if (++i == argc || sscanf(argv[i], "%f", &thresholdC) != 1)
          {
            error = 1;
            break;
          }
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

      fprintf(stderr, "Usage: blend_maps -map1 <map_name>\n");
      fprintf(stderr, "            -map2 <map_name>\n");
      fprintf(stderr, "            -output <output_map_name>\n");
      fprintf(stderr, "            [-weight1 <double>]\n");
      fprintf(stderr, "            [-weight2 <double>]\n");
      fprintf(stderr, "            [-inverse_weight1 <double>]\n");
      fprintf(stderr, "            [-inverse_weight2 <double>]\n");
      exit(1);
    }

  /* check that at least minimal parameters were supplied */
  if (map1Name[0] == '\0')
    Error("-map1 parameter must be specified.\n");
  if (map2Name[0] == '\0')
    Error("-map2 parameter must be specified.\n");
  if (outputMapName[0] == '\0')
    Error("-output parameter must be specified.\n");
  if (weight1 < 0.0)
    weight1 = 1.0;
  if (weight2 < 0.0)
    weight2 = 1.0;

  if (!ReadMap(map1Name, &map1, &mLevel, &mw, &mh,
	       &mxMin, &myMin, imgName, refName, msg))
    Error("Could not read map %s:\n%s\n", map1Name, msg);

  if (!ReadMap(map2Name, &map2, &mLevel2, &mw2, &mh2,
	       &mxMin2, &myMin2, imgName2, refName2, msg))
    Error("Could not read map %s:\n%s\n", map2Name, msg);

  omap = (MapElement*) malloc(mw*mh*sizeof(MapElement));

  totalWeight = weight1 + weight2;
  for (y = 0; y < mh; ++y)
    for (x = 0; x < mw; ++x)
      {
	x1 = map1[y*mw+x].x * (1 << mLevel);
	y1 = map1[y*mw+x].y * (1 << mLevel);
	c1 = map1[y*mw+x].c;
	xv = (x + mxMin) * (1 << mLevel) / (1 << mLevel2);
	yv = (y + myMin) * (1 << mLevel) / (1 << mLevel2);
	ix = ((int) floor(xv)) - mxMin2;
	iy = ((int) floor(yv)) - myMin2;
	rrx = xv - (mxMin2 + ix);
	rry = yv - (myMin2 + iy);

	if (ix < -1 || ix == -1 && rrx < 0.999 ||
	    ix == mw2-1 && rrx > 0.001 || ix >= mw2 ||
	    iy < -1 || iy == -1 && rry < 0.999 ||
	    iy == mh2-1 && rry > 0.001 || iy >= mh2)
	  {
	    omap[y*mw+x].x = 0.0;
	    omap[y*mw+x].y = 0.0;
	    omap[y*mw+x].c = 0.0;
	    continue;
	  }

	while (ix < 0)
	  {
	    ++ix;
	    rrx -= 1.0;
	  }
	while (iy < 0)
	  {
	    ++iy;
	    rry -= 1.0;
	  }
	while (ix >= mw2-1)
	  {
	    --ix;
	    rrx += 1.0;
	  }
	while (iy >= mh2-1)
	  {
	    --iy;
	    rry += 1.0;
	  }
	    
	rx00 = map2[iy*mw2+ix].x;
	ry00 = map2[iy*mw2+ix].y;
	rc00 = map2[iy*mw2+ix].c;
	rx01 = map2[(iy+1)*mw2+ix].x;
	ry01 = map2[(iy+1)*mw2+ix].y;
	rc01 = map2[(iy+1)*mw2+ix].c;
	rx10 = map2[iy*mw2+ix+1].x;
	ry10 = map2[iy*mw2+ix+1].y;
	rc10 = map2[iy*mw2+ix+1].c;
	rx11 = map2[(iy+1)*mw2+ix+1].x;
	ry11 = map2[(iy+1)*mw2+ix+1].y;
	rc11 = map2[(iy+1)*mw2+ix+1].c;
	rx = rx00 * (rrx - 1.0) * (rry - 1.0)
	  - rx10 * rrx * (rry - 1.0) 
	  - rx01 * (rrx - 1.0) * rry
	  + rx11 * rrx * rry;
	ry = ry00 * (rrx - 1.0) * (rry - 1.0)
	  - ry10 * rrx * (rry - 1.0) 
	  - ry01 * (rrx - 1.0) * rry
	  + ry11 * rrx * rry;
	rc = rc00;
	if (rc01 < rc)
	  rc = rc01;
	if (rc10 < rc)
	  rc = rc10;
	if (rc11 < rc)
	  rc = rc11;
	if (rc < 0.0)
	  Error("rc is negative! %f\n", rc);
	x2 = rx * (1 << mLevel2);
	y2 = ry * (1 << mLevel2);

	omap[y*mw+x].x = (weight1 * x1 + weight2 * x2) / totalWeight /
	  (1 << mLevel);
	omap[y*mw+x].y = (weight1 * y1 + weight2 * y2) / totalWeight /
	  (1 << mLevel);
	omap[y*mw+x].c = (c1 < rc) ? c1 : rc;
      }

  if (thresholdC != 0.0)
    for (y = 0; y < mh; ++y)
      for (x = 0; x < mw; ++x)
	if (omap[y*mw+x].c >= thresholdC)
	  omap[y*mw+x].c = 1.0;

  if (!WriteMap(outputMapName, omap, mLevel, mw, mh,
		mxMin, myMin,
		imgName, refName,
		UncompressedMap, msg))
    Error("Could not write map %s:\n%s\n",
	  outputMapName, msg);

  free(map1);
  free(map2);
  free(omap);

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
