#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <errno.h>

#include "imio.h"

#define LINE_LENGTH	255

int
main (int argc, char **argv)
{
  FILE *f;
  int mw, mh;
  int ix, iy;
  int x, y;
  int i;
  char msg[PATH_MAX + 1024];
  int mLevel;
  int mxMin, myMin;
  char imgName[PATH_MAX], refName[PATH_MAX];
  int ix0, iy0;
  double v;
  MapElement *map;
  int mw1, mh1;
  int n;
  float xv, yv;
  int mFactor;
  MapElement rmap[4];
  double angleSum;
  double theta;
  double ct, st;
  double stx, sty;
  double tx, ty;
  double xp, yp;
  double xMax, yMax;
  int rLevel;
  int rFactor;
  
  if (argc != 3)
    {
      fprintf(stderr, "Usage: bestrst input.map output.map\n");
      exit(1);
    }

  if (!ReadMap(argv[1], &map, &mLevel, &mw, &mh,
	       &mxMin, &myMin,
	       imgName, refName,
	       msg))
    {
      fprintf(stderr, "Could not read map %s:\n  error: %s\n",
	      argv[1], msg);
      exit(1);
    }
  mFactor = 1 << mLevel;

  /* go through all map points */
  mw1 = mw - 1;
  mh1 = mh - 1;
  n = 0;
  angleSum = 0.0;
  for (y = 0; y < mh; ++y)
    for (x = 0; x < mw; ++x)
      {
	if (map[y*mw+x].c == 0.0)
	  continue;
	if (x < mw1 && map[y*mw+x+1].c != 0.0)
	  {
	    angleSum += atan2(map[y*mw+x+1].y - map[y*mw+x].y,
			      map[y*mw+x+1].x - map[y*mw+x].x);
	    ++n;
	  }
	if (y < mh1 && map[(y+1)*mw+x].c != 0.0)
	  {
	    angleSum += atan2(map[(y+1)*mw+x].y - map[y*mw+x].y,
			      map[(y+1)*mw+x].x - map[y*mw+x].x) - M_PI/2.0;
	    ++n;
	  }
      }	
  if (n == 0)
    {
      fprintf(stderr, "No valid map points found.\n");
      exit(1);
    }
  theta = angleSum / n;
  ct = cos(theta);
  st = sin(theta);

  /* now determine translation */
  stx = 0.0;
  sty = 0.0;
  n = 0;
  for (y = 0; y < mh; ++y)
    for (x = 0; x < mw; ++x)
      {
	if (map[y*mw+x].c == 0.0)
	  continue;
	xv = mFactor * (x + mxMin);
	yv = mFactor * (y + myMin);
	xp = ct * xv - st * yv;
	yp = st * xv + ct * yv;
	stx += mFactor * map[y*mw+x].x - xp;
	sty += mFactor * map[y*mw+x].y - yp;
	++n;
      }
  tx = stx / n;
  ty = sty / n;

  xMax = mFactor * (mw + mxMin);
  yMax = mFactor * (mh + myMin);
  rLevel = 0;
  rFactor = 1;
  while (xMax > rFactor || yMax > rFactor)
    rFactor = 1 << (++rLevel);
  xMax = rFactor;
  yMax = rFactor;

  rmap[0].x = tx / rFactor;
  rmap[0].y = ty / rFactor;
  rmap[0].c = 1.0;
  rmap[1].x = (ct * xMax + tx) / rFactor;
  rmap[1].y = (st * xMax + ty) / rFactor;
  rmap[1].c = 1.0;
  rmap[2].x = (-st * yMax + tx) / rFactor;
  rmap[2].y = (ct * yMax + ty) / rFactor;
  rmap[2].c = 1.0;
  rmap[3].x = (ct * xMax - st * yMax + tx) / rFactor;
  rmap[3].y = (st * xMax + ct * yMax + ty) / rFactor;
  rmap[3].c = 1.0;

  printf("rigid transformation is theta = %f degrees  tx = %f  ty = %f\n",
	 theta * 180.0 / M_PI, tx, ty);
  
  /* write new map out */
  if (!WriteMap(argv[2], rmap, rLevel, 2, 2, 0, 0,
		imgName, refName,
		UncompressedMap, msg))
    {
      fprintf(stderr, "Could not write rigid map %s:\n%s\n",
	      argv[2], msg);
      exit(1);
    }

  /* deallocate all map data structures */
  free(map);
}
