/*
 * graphmap.c  - produce a pnm graph of a map
 *
 *   Copyright (c) 2010 Pittsburgh Supercomputing Center
 * 
 *  This program is distributed in the hope that it will    *
 *  be useful, but WITHOUT ANY WARRANTY; without even the   *
 *  implied warranty of MERCHANTABILITY or FITNESS FOR A    *
 *  PARTICULAR PURPOSE.  Neither Carnegie Mellon University *
 *  nor any of the authors assume any liability for         *
 *  damages, incidental or otherwise, caused by the         *
 *  installation or use of this software.                   *
 *                                                          *
 *  CLINICAL APPLICATIONS ARE NOT RECOMMENDED, AND THIS     *
 *  SOFTWARE HAS NOT BEEN EVALUATED BY THE UNITED STATES    *
 *  FDA FOR ANY CLINICAL USE.                               *
 *
 *  HISTORY
 *    2010  Written by Greg Hood (ghood@psc.edu)
 */

#include <stdio.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>

#include "imio.h"

void hsv_to_rgb (unsigned char* r, unsigned char *g, unsigned char *b,
		 float h, float s, float v);
void DrawLine(unsigned char *img, float px0, float py0, float px1, float py1,
	      int imageWidth, int imageHeight,
	      unsigned char r, unsigned char g, unsigned char b);

int
main (int argc, char **argv)
{
  int x, y;
  int i, j, k;
  FILE *f;
  float rx, ry, rc;
  int nx, ny;
  int ix, iy;
  int dx, dy;
  int xp, yp;
  unsigned char *grid;
  int factor;
  MapElement *map;
  char msg[PATH_MAX+256];
  float hue;
  unsigned char r, g, b;
  int n;
  float xv, yv;
  float minX, maxX, minY, maxY;
  float scale;
  float offsetX, offsetY;
  int outputGridWidth, outputGridHeight;
  int mLevel;
  int mw, mh;
  int mxMin, myMin;
  char imName0[PATH_MAX], imName1[PATH_MAX];

  outputGridWidth = 1800;
  outputGridHeight = 1400;

  if (argc != 3)
    {
      fprintf(stderr, "Usage: graphmap mapname.map output.pnm\n");
      exit(1);
    }

  if (!ReadMap(argv[1], &map, &mLevel,
	       &mw, &mh, &mxMin, &myMin,
	       imName0, imName1,
	       msg))
    {
      fprintf(stderr, "Could not read map %s:\n  error: %s\n",
	      argv[1], msg);
      exit(1);
    }

  factor = 1 << mLevel;
  /* determine the scaling and offset for the output grid */
  minX = 1.0e30;
  maxX = -1.0e30;
  minY = 1.0e30;
  maxY = -1.0e30;
  for (y = 0; y < mh; ++y)
    for (x = 0; x < mw; ++x)
      {
	rx = factor * map[y*mw + x].x;
	ry = factor * map[y*mw + x].y;
	rc = map[y*mw + x].c;
	if (rx < minX)
	  minX = rx;
	if (rx > maxX)
	  maxX = rx;
	if (ry < minY)
	  minY = ry;
	if (ry > maxY)
	  maxY = ry;
      }
      
  scale = 0.95 * outputGridWidth / (maxX - minX);
  if (0.95 * outputGridHeight / (maxY - minY) < scale)
    scale = 0.95 * outputGridHeight / (maxY - minY);
  offsetX = 0.5 * outputGridWidth - scale * 0.5 * (maxX + minX);
  offsetY = 0.5 * outputGridHeight - scale * 0.5 * (maxY + minY);
  printf("ogw = %d ogh = %d scale = %f offsetX = %f offsetY = %f\n",
	 outputGridWidth, outputGridHeight, scale, offsetX, offsetY);

  n = outputGridWidth * outputGridHeight * 3;
  grid = (unsigned char *) malloc(n);
  memset(grid, 255, n);

  nx = mw;
  ny = mh;
  for (y = 0; y < ny; ++y)
    for (x = 0; x < nx; ++x)
      for (dy = -1; dy <= 0; ++dy)
	{
	  yp = y + dy;
	  if (yp < 0)
	    continue;
	  for (dx = -1; dx <= 1; ++dx)
	    {
	      xp = x + dx;
	      if (xp < 0)
		continue;
	      if (xp >= nx || dy == 0 && dx >= 0)
		break;
	      
	      rc = map[y*nx+x].c;
	      if (map[yp*nx+xp].c < rc)
		rc = map[yp*nx+xp].c;
	      hue = 4.0 * rc;
	      hsv_to_rgb(&r, &g, &b, hue, 1.0, 1.0);

	      DrawLine(grid,
		       scale * factor * map[y*nx+x].x + offsetX,
		       scale * factor * map[y*nx+x].y + offsetY,
		       scale * factor * map[yp*nx+xp].x + offsetX,
		       scale * factor * map[yp*nx+xp].y + offsetY,
		       outputGridWidth, outputGridHeight,
		       r, g, b);
	      }
	    }

  for (y = 0; y < ny; ++y)
    for (x = 0; x < nx; ++x)
      {
	if (map[y*nx+x].c > 0.0)
	  {
	    ix = (int) (scale * factor * map[y*nx+x].x + offsetX);
	    iy = (int) (scale * factor * map[y*nx+x].y + offsetY);
	    for (dy = -1; dy <= 1; ++dy)
	      for (dx = -1; dx <= 1; ++dx)
		{
		  xp = ix + dx;
		  yp = iy + dy;
		  if (xp >= 0 && xp < outputGridWidth &&
		      yp >= 0 && yp < outputGridHeight)
		    {
		      grid[3*(yp * outputGridWidth + xp)] = 0;
		      grid[3*(yp * outputGridWidth + xp) + 1] = 0;
		      grid[3*(yp * outputGridWidth + xp) + 2] = 0;
		    }
		}
	  }
      }

  f = fopen(argv[2], "w");
  if (f == NULL)
    {
      fprintf(stderr, "Could not open output grid file %s\n", argv[2]);
      exit(1);
    }
  fprintf(f, "P6\n%d %d\n%d\n", outputGridWidth, outputGridHeight, 255);
  if (fwrite(grid, sizeof(unsigned char), n, f) != n)
    {
      fprintf(stderr, "Could not write to grid file %s\n", argv[2]);
      exit(1);
    }
  fclose(f);
}

void hsv_to_rgb (unsigned char* r, unsigned char *g, unsigned char *b,
		 float h, float s, float v)
{ 
  /* H is given on [0, 6] or UNDEFINED. S and V are given on [0, 1]. 
     RGB are each returned on [0, 255]. */
#define RETURN_RGB(rd,gr,bl)	{ *r = floor(255.999999 * rd); *g = floor(255.999999 * gr); *b = floor(255.999999 * bl); return; }
  
  float m, n, f;
  int i; 
  
  i = floor(h);
  f = h - i; 
  if ( !(i&1) )
    f = 1 - f; /* if i is even */
  m = v * (1 - s); 
  n = v * (1 - s * f); 
  
  switch (i)
    { 
    case 6: 
    case 0: RETURN_RGB(v, n, m); 
    case 1: RETURN_RGB(n, v, m); 
    case 2: RETURN_RGB(m, v, n);
    case 3: RETURN_RGB(m, n, v); 
    case 4: RETURN_RGB(n, m, v); 
    case 5: RETURN_RGB(v, m, n); 
    } 
} 

void
DrawLine(unsigned char *img, float px0, float py0, float px1, float py1,
	 int imageWidth, int imageHeight,
	 unsigned char r, unsigned char g, unsigned char b)
{
  int x0, y0, x1, y1;
  int steep;
  int tmp;
  int error;
  int yStep;
  int x, y;
  int dx, dy;
  float dist;
  float hue;
  unsigned char *p;

  x0 = (int) floor(px0);
  y0 = (int) floor(py0);
  x1 = (int) floor(px1);
  y1 = (int) floor(py1);

  // use the Bresenham Algorithm for lines
  //    (derived from the Wikipedia entry on
  //     "Bresenham's line algorithm")
  steep = abs(y1 - y0) > abs(x1 - x0);
  if (steep)
    {
      tmp = x0;
      x0 = y0;
      y0 = tmp;
      tmp = x1;
      x1 = y1;
      y1 = tmp;
    }
  if (x0 > x1)
    {
      tmp = x0;
      x0 = x1;
      x1 = tmp;
      tmp = y0;
      y0 = y1;
      y1 = tmp;
    }
  dx = x1 - x0;
  dy = abs(y1 - y0);
  error = -dx / 2;
  y = y0;
  if (y0 < y1)
    yStep = 1;
  else
    yStep = -1;
  for (x = x0; x <= x1; ++x)
    {
      if (steep)
	{
	  if (y >= 0 && y < imageWidth &&
	      x >= 0 && x < imageHeight)
	    {
	      p = &img[3*(x * imageWidth + y)];
	      *p++ = r;
	      *p++ = g;
	      *p++ = b;
	    }
	}
      else
	{
	  if (x >= 0 && x < imageWidth &&
	      y >= 0 && y < imageHeight)
	    {
	      p = &img[3*(y * imageWidth + x)];
	      *p++ = r;
	      *p++ = g;
	      *p++ = b;
	    }
	}
      error += dy;
      if (error > 0)
	{
	  y += yStep;
	  error -= dx;
	}
    }
}

