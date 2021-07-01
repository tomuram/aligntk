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
  
  InverseMap *invMap2;

  char outputMapName[PATH_MAX];
  MapElement *omap;

  int x, y;
  float x1, y1, c1;
  float xv, yv;
  float xp, yp;
  int ix, iy;
  float rx, ry, rc;
  float rx00, rx01, rx10, rx11;
  float ry00, ry01, ry10, ry11;
  float rc00, rc01, rc10, rc11;
  float rrx, rry;

  int dx, dy;
  double weight;
  double totalWeight;
  double trx, try;

  int i;
  int error;
  char msg[PATH_MAX+256];

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
    else if (strcmp(argv[i], "-inverse_map2") == 0)
      {
        if (++i == argc)
          {
            error = 1;
            break;
          }
        strcpy(map2Name, argv[i]);
	invert = 1;
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

      fprintf(stderr, "Usage: compose_maps -map1 <map_name>\n");
      fprintf(stderr, "            [-map2 <map_name>]\n");
      fprintf(stderr, "            [-inverse_map2 <map_name>]\n");
      fprintf(stderr, "            -output <output_map_name>]\n");
      fprintf(stderr, "            [-threshold_c <float>]\n");
      exit(1);
    }

  /* check that at least minimal parameters were supplied */
  if (map1Name[0] == '\0')
    Error("-map1 parameter must be specified.\n");
  if (map2Name[0] == '\0')
    Error("-map2 or -inverse_map2 parameter must be specified.\n");
  if (outputMapName[0] == '\0')
    Error("-output parameter must be specified.\n");

  if (!ReadMap(map1Name, &map1, &mLevel, &mw, &mh,
	       &mxMin, &myMin, imgName, refName, msg))
    Error("Could not read map %s:\n%s\n", map1Name, msg);

  if (!ReadMap(map2Name, &map2, &mLevel2, &mw2, &mh2,
	       &mxMin2, &myMin2, imgName2, refName2, msg))
    Error("Could not read map %s:\n%s\n", map2Name, msg);

  omap = (MapElement*) malloc(mw*mh*sizeof(MapElement));

  if (invert)
    {
      invMap2 = InvertMap(map2, mw2, mh2);

      for (y = 0; y < mh; ++y)
	for (x = 0; x < mw; ++x)
	  {
	    x1 = map1[y*mw+x].x * (1 << mLevel) / (1 << mLevel2);
	    y1 = map1[y*mw+x].y * (1 << mLevel) / (1 << mLevel2);
	    c1 = map1[y*mw+x].c;

	    // next line is temporary to correct for bug in warp3
	    if (c1 != 0.0 && c1 < 0.001)
	      c1 = 0.001;

	    if (c1 == 0.0)
	      {
		omap[y*mw+x].x = 0.0;
		omap[y*mw+x].y = 0.0;
		omap[y*mw+x].c = 0.0;
		continue;
	      }
	    xp = 0.0;
	    yp = 0.0;
	    if (Invert(invMap2, &xp, &yp, x1, y1))
	      {
		//		printf("Invert of (%f %f) produced (%f %f)\n",
		//		       x1, y1, xp, yp);
		omap[y*mw+x].x = xp * (1 << mLevel2) / (1 << mLevel);
		omap[y*mw+x].y = yp * (1 << mLevel2) / (1 << mLevel);
		omap[y*mw+x].c = c1;
	      }
	    else
	      {
		//		printf("Invert of (%f %f) failed\n", x1, y1);
		omap[y*mw+x].x = 0.0;
		omap[y*mw+x].y = 0.0;
		omap[y*mw+x].c = 0.0;
	      }
	  }
    }
  else
    {
      for (y = 0; y < mh; ++y)
	for (x = 0; x < mw; ++x)
	  {
	    x1 = map1[y*mw+x].x * (1 << mLevel);
	    y1 = map1[y*mw+x].y * (1 << mLevel);
	    c1 = map1[y*mw+x].c;

	    // next line is temporary to correct for bug in warp3
	    if (c1 != 0.0 && c1 < 0.001)
	      c1 = 0.001;

	    if (c1 == 0.0)
	      {
		omap[y*mw+x].x = 0.0;
		omap[y*mw+x].y = 0.0;
		omap[y*mw+x].c = 0.0;
		continue;
	      }
	    if (c1 < 0.0)
	      Error("c1 is negative! %f\n", c1);
	    xv = x1 / (1 << mLevel2);
	    yv = y1 / (1 << mLevel2);
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

	    rc = rc00;
	    if (rc01 < rc)
	      rc = rc01;
	    if (rc10 < rc)
	      rc = rc10;
	    if (rc11 < rc)
	      rc = rc11;
	    if (rc < 0.0)
	      Error("rc is negative! %f\n", rc);
	    
	    if (rc == 0.0)
	      {
		omap[y*mw+x].x = 0.0;
		omap[y*mw+x].y = 0.0;
		omap[y*mw+x].c = 0.0;
		continue;
	      }

	    rx = rx00 * (rrx - 1.0) * (rry - 1.0)
	      - rx10 * rrx * (rry - 1.0) 
	      - rx01 * (rrx - 1.0) * rry
	      + rx11 * rrx * rry;
	    ry = ry00 * (rrx - 1.0) * (rry - 1.0)
	      - ry10 * rrx * (rry - 1.0) 
	      - ry01 * (rrx - 1.0) * rry
	      + ry11 * rrx * rry;
	    xp = rx * (1 << mLevel2);
	    yp = ry * (1 << mLevel2);

	    if (isnan(xp) || isnan(yp))
	      Error("nan in xp or yp at (%d,%d).\n", x, y);
	    if (isnan(rc) || isnan(c1))
	      Error("nan in input maps at (%d,%d).\n", x, y);
	    omap[y*mw+x].x = xp / (1 << mLevel);
	    omap[y*mw+x].y = yp / (1 << mLevel);
	    if (c1 < rc)
	      rc = c1;
	    omap[y*mw+x].c = rc;
	  }
    }

  if (thresholdC != 0.0)
    for (y = 0; y < mh; ++y)
      for (x = 0; x < mw; ++x)
	if (omap[y*mw+x].c >= thresholdC)
	  omap[y*mw+x].c = 1.0;
	else
	  omap[y*mw+x].c = 0.0;

  if (!WriteMap(outputMapName, omap, mLevel, mw, mh,
		mxMin, myMin,
		imgName, invert ? imgName2 : refName2,
		UncompressedMap, msg))
    Error("Could not write map %s:\n%s\n",
	  outputMapName, msg);

  free(map1);
  free(map2);
  if (invert)
    FreeInverseMap(invMap2);
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
