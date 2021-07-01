/*
 * best_rigid - find best rigid transform approximating the given nonlinear map
 *
 *  Copyright (c) 2010-2011 National Resource for Biomedical
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
 *       NIH NCRR grant 5P41RR006009
 *
 *  HISTORY
 *    2010  Written by Greg Hood (ghood@psc.edu)
 */

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
  int error;
  char inputName[PATH_MAX];
  char outputName[PATH_MAX];
  
  error = 0;
  inputName[0] = '\0';
  outputName[0] = '\0';
  for (i = 0; i < argc; ++i)
    if (strcmp(argv[i], "-input") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-input error\n");
	    break;
	  }
	strcpy(inputName, argv[i]);
      }
    else if (strcmp(argv[i], "-output") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-output error\n");
	    break;
	  }
	strcpy(outputName, argv[i]);
      }
    else
      {
	fprintf(stderr, "Unknown option: %s\n", argv[i]);
	error = 1;
      }

  if (error)
    {
      if (i >= argc)
	fprintf(stderr, "Incomplete option: %s\n\n", argv[i-1]);

      fprintf(stderr, "Usage: best_rigid -input in.map -output out.map\n");
      exit(1);
    }
  if (inputName[0] == '\0' || outputName[0] == '\0')
    {
      fprintf(stderr, "Both -input and -output map files must be specified.\n");
      exit(1);
    }

  if (!ReadMap(inputName, &map, &mLevel, &mw, &mh,
	       &mxMin, &myMin,
	       imgName, refName,
	       msg))
    {
      fprintf(stderr, "Could not read map %s:\n  error: %s\n",
	      inputName, msg);
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
  if (!WriteMap(outputName, rmap, rLevel, 2, 2, 0, 0,
		imgName, refName,
		UncompressedMap, msg))
    {
      fprintf(stderr, "Could not write rigid map %s:\n%s\n",
	      outputName, msg);
      exit(1);
    }

  /* deallocate all map data structures */
  free(map);
}
