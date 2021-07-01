#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "invert.h"
#include "imio.h"

#define ONE_MINUS_EPSILON	0.9999994
#define NEGATIVE_EPSILON	-0.0000006

int
BilinearInvert (float *pAlpha, float *pBeta, float u1, float u2,
		float x1, float x2, float y1, float y2,
		float z1, float z2, float w1, float w2);

InverseMap*
InvertMap (MapElement *map, int nx, int ny)
{
  float xMin, xMax, yMin, yMax; 	/* bounds in destination space */
  InverseMap *inverseMap;
  InverseMapElement *e;
  int x, y;				/* point in source space */
  float xv, yv;				/* point in destination space */
  int minX, maxX, minY, maxY;           /* bounds of a single source tile in
					   destination space */
  float area;
  float scale;                          /* scale in destination space */
  int nxp, nyp;                         /* dimensions in destination space */
  int nxm1, nym1;
  int ix, iy;				/* point in destination space */
  float alpha, beta;
  int overlap;
  float pxv, pyv;
  float xvp, yvp;
  float pxvp, pyvp;
  int tx, ty;				/* tile in destination space */
  int corner;
  int side;
  float den;
  float ua, ub;
  int n;
  int maxElements;

  /* allocate the map */
  inverseMap = (InverseMap*) malloc(sizeof(InverseMap));
  inverseMap->map = map;
  inverseMap->nx = nx;
  inverseMap->ny = ny;
  
  /* find the bounding box */
  xMin = 1000000000.0;
  xMax = -1000000000.0;
  yMin = 1000000000.0;
  yMax = -1000000000.0;
  for (y = 0; y < ny; ++y)
    for (x = 0; x < nx; ++x)
      {
	xv = map[y*nx + x].x;
	yv = map[y*nx + x].y;
	if (xv < xMin)
	  xMin = xv;
	if (xv > xMax)
	  xMax = xv;
	if (yv < yMin)
	  yMin = yv;
	if (yv > yMax)
	  yMax = yv;
      }
  inverseMap->xMin = xMin;
  inverseMap->yMin = yMin;

  /* create the table and initialize entries*/
  area = (xMax - xMin) * (yMax - yMin);
  scale = sqrt(area / nx / ny);
  nxp = (int) ceil((xMax - xMin) / scale);
  nyp = (int) ceil((yMax - yMin) / scale);
  inverseMap->inverseMap = (InverseMapElement*) malloc(nxp * nyp * sizeof(InverseMapElement));
  inverseMap->scale = scale;
  inverseMap->nxp = nxp;
  inverseMap->nyp = nyp;
  e = inverseMap->inverseMap;
  for (iy = 0; iy < nyp; ++iy)
    for (ix = 0; ix < nxp; ++ix)
      {
	e->minX = 1000000000;
	e->maxX = -1000000000;
	e->minY = 1000000000;
	e->maxY = -1000000000;
	++e;
      }

  /* set table entries */
  nym1 = ny - 1;
  nxm1 = nx - 1;
  for (y = 0; y < nym1; ++y)
    for (x = 0; x < nxm1; ++x)
      {
	minX = 1000000000;
	minY = 1000000000;
	maxX = -1000000000;
	maxY = -1000000000;
	for (corner = 0; corner < 4; ++corner)
	  {
	    switch (corner)
	      {
	      case 0:
		xv = map[y*nx + x].x;
		yv = map[y*nx + x].y;
		break;
	      case 1:
		xv = map[y*nx + x + 1].x;
		yv = map[y*nx + x + 1].y;
		break;
	      case 2:
		xv = map[(y+1)*nx + x + 1].x;
		yv = map[(y+1)*nx + x + 1].y;
		break;
	      case 3:
		xv = map[(y+1)*nx + x].x;
		yv = map[(y+1)*nx + x].y;
		break;
	      }
	    ix = floor((xv - xMin) / scale);
	    iy = floor((yv - yMin) / scale);
	    if (ix < minX)
	      minX = ix;
	    if (ix > maxX)
	      maxX = ix;
	    if (iy < minY)
	      minY = iy;
	    if (iy > maxY)
	      maxY = iy;
	  }
	for (ty = minY; ty <= maxY; ++ty)
	  for (tx = minX; tx <= maxX; ++tx)
	    // for (ty = minY; ty < maxY; ++ty)
	    //	  for (tx = minX; tx < maxX; ++tx)
	    {
	      /* check for overlap between quadrilateral and square */
	      
	      /* first check if the center of the square lies within
		 the quadrilateral */
	      xv = scale * (tx + 0.5) + xMin;
	      yv = scale * (ty + 0.5) + yMin;
	      overlap = 0;
	      if (BilinearInvert(&alpha, &beta, xv, yv,
				 map[y*nx + x].x, map[y*nx + x].y,
				 map[(y+1)*nx + x].x, map[(y+1)*nx + x].y,
				 map[(y+1)*nx + (x+1)].x, map[(y+1)*nx + (x+1)].y,
				 map[y*nx + (x+1)].x, map[y*nx + (x+1)].y))
		overlap = 1;

	      pxv = map[(ty+1)*nx + tx].x;
	      pyv = map[(ty+1)*nx + tx].y;
	      for (corner = 0; corner < 4 && !overlap; ++corner)
		{
		  switch (corner)
		    {
		    case 0:
		      xv = map[y*nx + x].x;
		      yv = map[y*nx + x].y;
		      break;
		    case 1:
		      xv = map[y*nx + x + 1].x;
		      yv = map[y*nx + x + 1].y;
		      break;
		    case 2:
		      xv = map[(y+1)*nx + x + 1].x;
		      yv = map[(y+1)*nx + x + 1].y;
		      break;
		    case 3:
		      xv = map[(y+1)*nx + x].x;
		      yv = map[(y+1)*nx + x].y;
		      break;
		    }
		  pxvp = tx * scale + xMin;
		  pyvp = ty * scale + yMin;
		  for (side = 0; side < 4 && !overlap; ++side)
		    {
		      switch (side)
			{
			case 0:
			  xvp = (tx + 1) * scale + xMin;
			  yvp = pyvp;
			  break;
			case 1:
			  yvp = (ty + 1) * scale + yMin;
			  break;
			case 2:
			  xvp = tx * scale + xMin;
			  break;
			case 3:
			  yvp = ty * scale + yMin;
			  break;
			}
		      den = (yvp - pyvp) * (xv - pxv) - (xvp - pxvp) * (yv - pyv);
		      if (den < -1.0e-15 || den > 1.0e-15)
			{
			  ua = ((xvp - pxvp) * (pyv - pyvp) - (yvp - pyvp) * (pxv - pxvp)) / den;
			  ub = ((xv - pxv) * (pyv - pyvp) - (yv - pyv) * (pxv - pxvp)) / den;
			  overlap = (ua >= 0.0) && (ua <= 1.0) && (ub >= 0.0) && (ub <= 1.0);
			}
		      pxvp = xvp;
		      pyvp = yvp;
		    }
		  pxv = xv;
		  pyv = yv;
		}
	      if (overlap)
		{
		  e = &(inverseMap->inverseMap[ty * nxp + tx]);
		  if (x < e->minX)
		    e->minX = x;
		  if (x > e->maxX)
		    e->maxX = x;
		  if (y < e->minY)
		    e->minY = y;
		  if (y > e->maxY)
		    e->maxY = y;
		}
	    }
      }

  e = inverseMap->inverseMap;
  maxElements = 1;
  for (iy = 0; iy < nyp; ++iy)
    for (ix = 0; ix < nxp; ++ix)
      {
	if (e->minX <= e->maxX)
	  {
	    n = (e->maxX - e->minX + 1) * (e->maxY - e->minY + 1);
	    if (n > maxElements)
	      maxElements = n;
	  }
	++e;
      }
  inverseMap->tried = (long long*) malloc(maxElements * sizeof(long long));
  inverseMap->marked = (unsigned char *) malloc(nx * ny * sizeof(unsigned char));
  memset(inverseMap->marked, 0, nx * ny * sizeof(unsigned char));
  return(inverseMap);
}

int
Invert (InverseMap *inverseMap, float *xv, float *yv, float xvp, float yvp)
{
  int ix, iy;
  int minX, maxX, minY, maxY;
  int nxp, nyp;
  InverseMapElement *e;
  int nx;
  int nTried;
  int nNeighborsTried;
  int x, y;
  MapElement *me;
  int i;
  float alpha, beta;

  ix = (int) floor((xvp - inverseMap->xMin) / inverseMap->scale);
  iy = (int) floor((yvp - inverseMap->yMin) / inverseMap->scale);
  if (ix < 0 || iy < 0 || ix >= inverseMap->nxp || iy >= inverseMap->nyp)
    return(0);
  e = &(inverseMap->inverseMap[iy * inverseMap->nxp + ix]);
  minX = e->minX;
  maxX = e->maxX;
  if (minX > maxX)
    return(0);
  minY = e->minY;
  maxY = e->maxY;
  nx = inverseMap->nx;
  x = inverseMap->lastX;
  y = inverseMap->lastY;
  nTried = 0;
  nNeighborsTried = 0;
  for (;;)
    {
      if (x < minX)
	x = minX;
      else if (x > maxX)
	x = maxX;
      if (y < minY)
	y = minY;
      else if (y > maxY)
	y = maxY;
      /* check if we've already tried this MapElement */
      if (inverseMap->marked[y*nx+x])
	{
	  if (nTried >= (maxX - minX + 1) * (maxY - minY + 1))
	    {
	      fprintf(stderr, "Error: invert() could not locate the proper MapElement\n");
	      /* unmark all MapElements that were tried */
	      for (i = 0; i < nTried; ++i)
		inverseMap->marked[inverseMap->tried[i]] = 0;
	      return(0);
	    }
	  if (nNeighborsTried < 8)
	    {
	      x += (int) floor(3.0 * drand48() - 1.0);
	      y += (int) floor(3.0 * drand48() - 1.0);
	    }
	  else
	    {
	      nNeighborsTried = 0;
	      x = minX + (int) floor(drand48() * (maxX - minX + 1));
	      y = minY + (int) floor(drand48() * (maxY - minY + 1));
	    }
	  continue;
	}
      me = &(inverseMap->map[y * nx + x]);
      if (BilinearInvert(&alpha, &beta,
			 xvp, yvp,
			 me->x, me->y,
			 (me+nx)->x, (me+nx)->y,
			 (me+nx+1)->x, (me+nx+1)->y,
			 (me+1)->x, (me+1)->y))
	{
	  /* unmark all MapElements that were tried */
	  for (i = 0; i < nTried; ++i)
	    inverseMap->marked[inverseMap->tried[i]] = 0;
	  *xv = x + alpha;
	  *yv = y + beta;
	  return(1);
	}
      inverseMap->marked[y*nx + x] = 1;
      inverseMap->tried[nTried] = y*nx + x;
      ++nTried;
      ++nNeighborsTried;

      if (alpha != 0.0 || beta != 0.0)
	{
	  /* estimate the x,y indices of the MapElement to try */
	  x += (int) floor(alpha);
	  y += (int) floor(beta);
	}
      else
	{
	  x = minX + (int) floor(drand48() * (maxX - minX + 1));
	  y = minY + (int) floor(drand48() * (maxY - minY + 1));
	  nNeighborsTried = 0;
	}
    }
}

void
FreeInverseMap (InverseMap *inverseMap)
{
  free(inverseMap->inverseMap);
  free(inverseMap->tried);
  free(inverseMap->marked);
  free(inverseMap);
}

int
BilinearInvert (float *pu, float *pv, float xp, float yp,
		float x00, float y00, float x01, float y01,
		float x11, float y11, float x10, float y10)
{
  double xt, yt;
  double A, B, C;
  double disc, sdisc;
  double size;
  double den1, den2;
  double t;
  double dx, dy;
  double dxp00, dyp00;
  double dxp10, dyp10;
  double dxp01, dyp01;
  double xx, xy, yx, yy;
  double u1, v1, u2, v2;
  int sol1, sol2;
  double dist1, dist2;
  int valid1, valid2;
  double det;
  double dxh, dyh, dxv, dyv;
  double dxhp, dyhp, dxvp, dyvp;
  double d;

#if 0
  printf("xp = %f yp = %f\n", xp, yp);
  printf("x00 = %f y00 = %f\n", x00, y00);
  printf("x01 = %f y01 = %f\n", x01, y01);
  printf("x11 = %f y11 = %f\n", x11, y11);
  printf("x10 = %f y10 = %f\n", x10, y10);
#endif

  /* First determine if the point (xp, yp) is inside the quadrilateral; if so,
     we must return a valid value */

  /* We must be very careful here in regard to boundary conditions:
      the point (xp, yp) is considered inside the quadrilateral if it is equal
      to the upper left corner (x00, y00), lies on the line segment
      from (x00, y00) to (x10, y10), exclusive, or lies on the line
      segment from (x00, y00) to (x01, y01), exclusive.

      The tests for these boundary conditions should be computed in a way
      that yields consistent results for adjacent quadrilaterals. */

  dxp00 = xp - x00;
  dyp00 = yp - y00;
  dxp10 = xp - x10;
  dyp10 = yp - y10;
  dxp01 = xp - x01;
  dyp01 = yp - y01;
  
  dxh = x10 - x00;
  dyh = y10 - y00;
  dxv = x01 - x00;
  dyv = y01 - y00;
  dxvp = x11 - x10;
  dyvp = y11 - y10;
  dxhp = x11 - x01;
  dyhp = y11 - y01;

  size = hypot(dxh, dyh);
  d = hypot(dxv, dyv);
  if (d > size)
    size = d;
  d = hypot(dxvp, dyvp);
  if (d > size)
    size = d;
  d = hypot(dxhp, dyhp);
  if (d > size)
    size = d;

  if (dxh * dyp00 - dyh * dxp00 < 0.0 ||
      dxv * dyp00 - dyv * dxp00 > 0.0 ||
      dxvp * dyp10 - dyvp * dxp10 <= 0.0 ||
      dxhp * dyp01 - dyhp * dxp01 >= 0.0)
    {
      /* point lies outside the quadrilateral */

      det = dxh * dyv - dxv * dyh;
      if (det < 1.0e-10 * size)
	{
	  fprintf(stderr, "Error: determinant of the map element is zero\n");
	  *pu = 0.0;
	  *pv = 0.0;
	  return(0);
	}
      u1 = (dyv * dxp00 - dxv * dyp00) / det;
      v1 = (-dyh * dxp00 + dxh * dyp00) / det;

      /* make sure u, v are not in the [0,1)x[0,1) square */
      if (u1 >= 0.0 && u1 < 1.0 && v1 >= 0.0 && v1 < 1.0)
	{
	  /* coerce it outside */
	  if (fabs(u1 - 0.5) >= fabs(v1 - 0.5))
	    {
	      if (u1 < 0.5)
		u1 = NEGATIVE_EPSILON;
	      else
		u1 = 1.0;
	    }
	  else
	    {
	      if (v1 < 0.5)
		v1 = NEGATIVE_EPSILON;
	      else
		v1 = 1.0;
	    }
	}
      /*      printf("returning u = %f v = %f\n", u1, v1); */
      *pu = u1;
      *pv = v1;
      return(1);
    }

  xt = dxhp - dxh;
  yt = dyhp - dyh;

  A = - xt * dyh + yt * dxh;
  B = xt * dyp00 - dxv * dyh
    - yt * dxp00 + dyv * dxh;
  C = dxv * dyp00 - dyv * dxp00;

  if (fabs(A) < 1.0e-10 * size * size)
    {
      /* edges are close to parallel so we just have to solve a linear equation */
      if (fabs(B) < 1.0e-10 * size * size)
	{
	  fprintf(stderr, "BilinearInvert: both A and B in quadratic equation are effectively zero\n");
	  *pu = 0.0;
	  *pv = 0.0;
	  return(0);
	}
      u1 = - C / B;
      den1 = dxv + xt * u1;
      den2 = dyv + yt * u1;
      if (fabs(den1) < 1.0e-5 * size &&
	  fabs(den2) < 1.0e-5 * size)
	{
	  /* report error */
	  fprintf(stderr, "BilinearInvert: den1 (%f) and den2 (%f) effectively equal 0\n", den1, den2);
	  *pu = 0.0;
	  *pv = 0.0;
	  return(0);
	}
      if (fabs(den1) > fabs(den2))
	v1 = (-dxh * u1 + dxp00) / den1;
      else
	v1 = (-dyh * u1 + dyp00) / den2;

      /* coerce solution into range */
      if (u1 <= 0.0)
	u1 = 0.0;
      else if (u1 >= 1.0)
	u1 = ONE_MINUS_EPSILON;
      if (v1 <= 0.0)
	v1 = 0.0;
      else if (v1 >= 1.0)
	v1 = ONE_MINUS_EPSILON;
      /*      printf("u1 = %f  v1 = %f\n", u1, v1); */
      *pu = u1;
      *pv = v1;
      return(1);
    }

  sol1 = 0;
  sol2 = 0;
  disc = B*B - 4.0*A*C;
  if (disc < 0.0)
    {
      fprintf(stderr, "Error: discriminant in BilinearInvert is negative (%f %f %f %f)\n",
	      A, B, C, disc);
      *pu = 0.0;	
      *pv = 0.0;
      return(0);
    }
  sdisc = sqrt(disc);
  t = -0.5 * (B + sdisc);
  u1 = t / A;
  den1 = dxv + xt * u1;
  den2 = dyv + yt * u1;
  if (fabs(den1) >= 1.0e-5 * size ||
      fabs(den2) >= 1.0e-5 * size)
    {
      if (fabs(den1) > fabs(den2))
	v1 = (-dxh * u1 + dxp00) / den1;
      else
	v1 = (-dyh * u1 + dyp00) / den2;
      sol1 = 1;
      valid1 = 1;
      dist1 = 0.0;
      if (u1 < 0.0)
	{
	  valid1 = 0;
	  dist1 += fabs(u1);
	}
      else if (u1 >= 1.0)
	{
	  valid1 = 0;
	  dist1 += u1 - 1.0;
	}
      if (v1 < 0.0)
	{
	  valid1 = 0;
	  dist1 += fabs(v1);
	}
      else if (v1 >= 1.0)
	{
	  valid1 = 0;
	  dist1 += v1 - 1.0;
	}
    }

  t = -0.5 * (B - sdisc);
  u2 = t / A;
  den1 = dxv + xt * u2;
  den2 = dyv + yt * u2;
  if (fabs(den1) >= 1.0e-5 * size ||
      fabs(den2) >= 1.0e-5 * size)
    {
      if (fabs(den1) > fabs(den2))
	v2 = (-dxh * u2 + dxp00) / den1;
      else
	v2 = (-dyh * u2 + dyp00) / den2;
      sol2 = 1;
      valid2 = 1;
      dist2 = 0.0;
      if (u2 < 0.0)
	{
	  valid2 = 0;
	  dist2 += fabs(u2);
	}
      else if (u2 >= 1.0)
	{
	  valid2 = 0;
	  dist2 += u2 - 1.0;
	}
      if (v2 < 0.0)
	{
	  valid2 = 0;
	  dist2 += fabs(v2);
	}
      else if (v2 >= 1.0)
	{
	  valid2 = 0;
	  dist2 += v2 - 1.0;
	}
    }

  if (sol1 && valid1)
    {
      *pu = u1;
      *pv = v1;
      return(1);
    }
  if (sol2 && valid2)
    {
      *pu = u2;
      *pv = v2;
      return(1);
    }
  if (sol1 && (!sol2 || dist1 <= dist2))
    {
      if (dist1 > 1.0e-5 * size)
	{
	  /* report error */
	  fprintf(stderr, "BilinearInvert: solution1 out-of-bounds\n");
	  *pu = 0.0;	
	  *pv = 0.0;
	  return(0);
	}
	  
      /* coerce solution into range */
      if (u1 <= 0.0)
	u1 = 0.0;
      else if (u1 >= 1.0)
	u1 = ONE_MINUS_EPSILON;
      if (v1 <= 0.0)
	v1 = 0.0;
      else if (v1 >= 1.0)
	v1 = ONE_MINUS_EPSILON;
      *pu = u1;
      *pv = v1;
      return(1);
    }
  else if (sol2)
    {
      if (dist2 > 1.0e-5 * size)
	{
	  /* report error */
	  fprintf(stderr, "BilinearInvert: solution2 out-of-bounds\n");
	  *pu = 0.0;	
	  *pv = 0.0;
	  return(0);
	}

      /* coerce solution into range */
      if (u2 <= 0.0)
	u2 = 0.0;
      else if (u2 >= 1.0)
	u2 = ONE_MINUS_EPSILON;
      if (v2 <= 0.0)
	v2 = 0.0;
      else if (v2 >= 1.0)
	v2 = ONE_MINUS_EPSILON;
      *pu = u2;
      *pv = v2;
      return(1);
    }
  else
    {
      /* report error */
      fprintf(stderr, "BilinearInvert: no solution found\n");
      *pu = 0.0;
      *pv = 0.0;
      return(0);
    }
}
