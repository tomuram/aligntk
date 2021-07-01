/* function to compute the Euclidean distance transform of a binary image;
   algorithm adapted from A. Meijster, J.B.T.M. Roerdink and W.H. Hesselink,
   A general algorithm for computing distance transforms in linear time. In:
   Mathematical Morphology and its Apllications to Image and Signal Processing,
   Kluwer Acad. Publ., 2000, pp. 331-340. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "dt.h"


void
computeDistance (int type,
		 int nx, int ny,
		 unsigned char *mask,
		 float *dist)
{
  int *g;
  int mw;
  int q;
  int *s, *t;
  int u;
  int w, z;
  int x, y;

  // make sure there at least one marked pixel in mask; otherwise,
  //   all distances are infinity
  mw = (nx + 7) >> 3;
  for (y = 0; y < ny; ++y)
    for (x = 0; x < mw; ++x)
      if (mask[y*mw+x] != 0)
	goto nonempty;
  for (y = 0; y < ny; ++y)
    for (x = 0; x < nx; ++x)
      dist[y*nx+x] = 1.0e30;
  return;

 nonempty:
  g = (int *) malloc(nx * ny * sizeof(int));
  s = (int *) malloc(nx * sizeof(int));
  t = (int *) malloc(nx * sizeof(int));
  for (x = 0; x < nx; ++x)
    {
      if ((mask[x >> 3] & (0x80 >> (x & 7))) != 0)
	g[x] = 0;
      else
	g[x] = nx + ny;
      for (y = 1; y < ny; ++y)
	if ((mask[y * mw + (x >> 3)] & (0x80 >> (x & 7))) != 0)
	  g[y*nx + x] = 0;
	else
	  g[y*nx + x] = g[(y - 1) * nx + x] + 1;
      for (y = ny - 2; y >= 0; --y)
	if (g[(y + 1) * nx + x] < g[y*nx + x])
	  g[y*nx + x] = g[(y + 1) * nx + x] + 1;
    }

  for (y = 0; y < ny; ++y)
    {
      q = 0;
      s[0] = 0;
      t[0] = 0;
      for (u = 1; u < nx; ++u)
	{
	  switch (type)
	    {
	    case EUCLIDEAN_DISTANCE:
	    case EUCLIDEAN_DISTANCE_SQUARED:
	      while (q >= 0 &&
		     (t[q] - s[q]) * (t[q] - s[q]) +
		     g[y*nx + s[q]] * g[y*nx + s[q]] >
		     (t[q] - u) * (t[q] - u) + g[y*nx + u] * g[y*nx + u])
		--q;
	      break;
	    case MANHATTAN_DISTANCE:
	      while (q >= 0 &&
		     abs(t[q] - s[q]) + g[y*nx + s[q]] >
		     abs(t[q] - u) + g[y*nx + u])
		--q;
	      break;
	    case CHESSBOARD_DISTANCE:
	      while (q >= 0)
		{
		  w = abs(t[q] - s[q]);
		  if (g[y*nx + s[q]] > w)
		    w = g[y*nx + s[q]];
		  z = abs(t[q] - u);
		  if (g[y*nx + u] > z)
		    z = g[y*nx + u];
		  if (w <= z)
		    break;
		  --q;
		}
	      break;
	    }
	  if (q < 0)
	    {
	      q = 0;
	      s[0] = u;
	    }
	  else
	    {
	      switch (type)
		{
		case EUCLIDEAN_DISTANCE:
		case EUCLIDEAN_DISTANCE_SQUARED:
		  w = (u * u - s[q] * s[q] + g[y*nx+u] * g[y*nx+u] -
		       g[y*nx+s[q]] * g[y*nx+s[q]]) / (2 * (u - s[q])) + 1;
		  break;
		case MANHATTAN_DISTANCE:
		  if (g[y*nx + u] >= g[y*nx + s[q]] + u - s[q])
		    w = 1000000000;
		  else if (g[y*nx + s[q]] > g[y*nx + u] + u - s[q])
		    w = -1000000000;
		  else
		    w = (g[y*nx + u] - g[y*nx + s[q]] + u + s[q]) / 2 + 1;
		  break;
		case CHESSBOARD_DISTANCE:
		  if (g[y*nx + s[q]] <= g[y*nx + u])
		    {
		      w = s[q] + g[y*nx + u];
		      z = (s[q] + u) / 2;
		      if (z > w)
			w = z;
		    }
		  else
		    {
		      w = u - g[y*nx + s[q]];
		      z = (s[q] + u) / 2;
		      if (z < w)
			w = z;
		    }
		  ++w;
		  break;
		}
	      if (w < nx)
		{
		  ++q;
		  s[q] = u;
		  t[q] = w;
		}
	    }
	}
      for (u = nx-1; u >= 0; --u)
	{
	  switch (type)
	    {
	    case EUCLIDEAN_DISTANCE:
	      dist[y * nx + u] = sqrt((double) ((u - s[q]) * (u - s[q]) +
						g[y*nx + s[q]] * g[y*nx + s[q]]));
	      break;
	    case EUCLIDEAN_DISTANCE_SQUARED:
	      dist[y * nx + u] = (u - s[q]) * (u - s[q]) +
		g[y*nx + s[q]] * g[y*nx + s[q]];
	      break;
	    case MANHATTAN_DISTANCE:
	      dist[y * nx + u] = abs(u - s[q]) + g[y*nx + s[q]];
	      break;
	    case CHESSBOARD_DISTANCE:
	      w = abs(u - s[q]);
	      if (g[y*nx + s[q]] > w)
		w = g[y*nx + s[q]];
	      dist[y * nx + u] = w;
	      break;
	    }
	  if (u == t[q])
	    --q;
	}
    }
  free(g);
  free(s);
  free(t);
}
