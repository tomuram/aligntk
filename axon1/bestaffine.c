// ComputeMapping and svd adapted from nform.cpp and utility.cpp of
//   Reconstruct by John Fiala
//
//    Those portions are covered by the following:
//    Copyright (C) 1999-2007  John Fiala (fiala@bu.edu)
//
//    This is free software created with funding from the NIH. You may
//    redistribute it and/or modify it under the terms of the GNU General
//    Public License published by the Free Software Foundation (www.gnu.org).
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License version 2 for more details.



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

void ComputeMapping (double *a, double *b,
		     double *fx, double *fy, double *rx, double *ry,
		     int nopts, int d);
void svd (double **a, int m, int n, double *w, double **v);
double pythag (double a, double b);

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
  double theta;
  double ct, st;
  double stx, sty;
  double tx, ty;
  double xp, yp;
  double xMax, yMax;
  int rLevel;
  int rFactor;
  double a[6], b[6];
  double *fx, *fy, *rx, *ry;

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
  n = 0;
  for (y = 0; y < mh; ++y)
    for (x = 0; x < mw; ++x)
      if (map[y*mw+x].c != 0.0)
	++n;
  fx = (double *) malloc(n * sizeof(double));
  fy = (double *) malloc(n * sizeof(double));
  rx = (double *) malloc(n * sizeof(double));
  ry = (double *) malloc(n * sizeof(double));
  i = 0;
  for (y = 0; y < mh; ++y)
    for (x = 0; x < mw; ++x)
      if (map[y*mw+x].c != 0.0)
	{
	  fx[i] = mFactor * (x + mxMin);
	  fy[i] = mFactor * (y + myMin);
	  rx[i] = mFactor * map[y*mw+x].x;
	  ry[i] = mFactor * map[y*mw+x].y;
	  ++i;
	}
  if (i != n)
    {
      fprintf(stderr, "Internal error.\n");
      exit(1);
    }

  if (n < 4)
    {
      fprintf(stderr, "Not enough valid map points found.\n");
      exit(1);
    }

  /* now determine affine transformation */
  ComputeMapping(a, b, fx, fy, rx, ry, n, 3);

  printf("Best affine:  x = %gx %+gy %+g\n",
	 a[1], a[2], a[0]);
  printf("              y = %gx %+gy %+g\n",
	 b[1], b[2], b[0]);

  xMax = mFactor * (mw + mxMin);
  yMax = mFactor * (mh + myMin);
  rLevel = 0;
  rFactor = 1;
  while (xMax > rFactor || yMax > rFactor)
    rFactor = 1 << (++rLevel);
  xMax = rFactor;
  yMax = rFactor;

  rmap[0].x = a[0] / rFactor;
  rmap[0].y = b[0] / rFactor;
  rmap[0].c = 1.0;
  rmap[1].x = (a[1] * xMax + a[0]) / rFactor;
  rmap[1].y = (b[1] * xMax + b[0]) / rFactor;
  rmap[1].c = 1.0;
  rmap[2].x = (a[2] * yMax + a[0]) / rFactor;
  rmap[2].y = (b[2] * yMax + b[0]) / rFactor;
  rmap[2].c = 1.0;
  rmap[3].x = (a[1] * xMax + a[2] * yMax + a[0]) / rFactor;
  rmap[3].y = (b[1] * xMax + b[2] * yMax + b[0]) / rFactor;
  rmap[3].c = 1.0;

  /* write new map out */
  if (!WriteMap(argv[2], rmap, rLevel, 2, 2, 0, 0,
		imgName, refName,
		UncompressedMap, msg))
    {
      fprintf(stderr, "Could not write affine map %s:\n%s\n",
	      argv[2], msg);
      exit(1);
    }

  /* deallocate all map data structures */
  free(map);

  /* now determine betst quadratic transformation (for comparison) */
  ComputeMapping(a, b, fx, fy, rx, ry, n, 6);

  printf("Best quadratic:  x = %gx^2 %+gy^2 %+gxy %+gx %+gy %+g\n",
	 a[4], a[5], a[3], a[1], a[2], a[0]);
  printf("                 y = %gx^2 %+gy^2 %+gxy %+gx %+gy %+g\n",
	 b[4], b[5], b[3], b[1], b[2], b[0]);

  return(0);
}

void 
ComputeMapping (double *a, double *b,
		double *fx, double *fy, double *rx, double *ry,
		int nopts, int d)
{
  int i, j, k;
  double **U, *w, **V, *ux;
  double s, scale, wmin, wmax; // vectors have zero as first index

  U = (double**) malloc(nopts * sizeof(double*));
  for (k=0;k<nopts;k++) U[k] = (double*) malloc(d * sizeof(double));
  V = (double**) malloc(d * sizeof(double*));
  for (k=0;k<d;k++) V[k] = (double*) malloc(d * sizeof(double));
  w = (double*) malloc(d * sizeof(double));
  ux = (double*) malloc(d * sizeof(double));

  scale = 1.0;  // get scaling factor for matrix U
  for (i=0; i<nopts; i++)
    {
      if ( fx[i] > scale ) scale = fx[i];
      if ( fy[i] > scale ) scale = fy[i];
    }
  // compute U from points, d limits
  for (i=0; i<nopts; i++) {
    // which terms of a,b are computed
    if ( d > 0 ) U[i][0] = 1.0;
    if ( d > 1 ) U[i][1] = fx[i]/scale;
    if ( d > 2 ) U[i][2] = fy[i]/scale;
    if ( d > 3 ) U[i][3] = fx[i]*fy[i]/(scale*scale);
    if ( d > 4 ) U[i][4] = fx[i]*fx[i]/(scale*scale);
    if ( d > 5 ) U[i][5] = fy[i]*fy[i]/(scale*scale);
  }

  svd( U, nopts, d, w, V );   // compute SVD of U

  wmax = 0.0;                 // find max singular value
  for (i=0; i<d; i++) if (w[i] > wmax) wmax = w[i];
  wmin = wmax*1.0e-11;       // zero all small singular values
  for (i=0; i<d; i++) if (w[i] < wmin) w[i] = 0.0;
  //                  t
  for (i=0; i<d; i++) {      // solve U*diag(w)*V *a = rx
    s = 0.0;
    if ( w[i] ) {
      for (j=0; j<nopts; j++)      //      t
	s += U[j][i]*rx[j];        // s = U *rx
      s /= w[i];                   // ux = diag(1/w)*s
    }
    ux[i] = s;
  }
  for (i=0; i<d; i++) {
    // this->a = V*ux
    s = 0.0;
    for (j=0; j<d; j++) s += V[i][j]*ux[j];
    a[i] = s;
  }
  //                  t
  for (i=0; i<d; i++) {
    // solve U*diag(w)*V *b = ry
    s = 0.0;
    if ( w[i] ) {
      for (j=0; j<nopts; j++)
	//      t
	s += U[j][i]*ry[j];
      // s = U *ry
      s /= w[i];
      // ux = diag(1/w)*s
    }
    ux[i] = s;
  }
  for (i=0; i<d; i++) {
    // this->b = V*ux
    s = 0.0;
    for (j=0; j<d; j++) s += V[i][j]*ux[j];
    b[i] = s;
  }
  // unscale the resulting params
  if ( d > 1 ) { a[1] /= scale; b[1] /= scale; }
  if ( d > 2 ) { a[2] /= scale; b[2] /= scale; }
  if ( d > 3 ) { a[3] /= scale*scale; b[3] /= scale*scale; }
  if ( d > 4 ) { a[4] /= scale*scale; b[4] /= scale*scale; }
  if ( d > 5 ) { a[5] /= scale*scale; b[5] /= scale*scale; }

  for (i=d; i<6; i++) {
    // set remaining parameters to zero!
    a[i] = 0.0;
    b[i] = 0.0;
  }

  // free working memory
  free(ux);
  free(w);
  for (k=0;k<d;k++) free(V[k]);
  free(V);
  for (k=0;k<nopts;k++) free(U[k]);
  free(U);
}


void
svd (double **a, int m, int n, double *w, double **v)
// compute the singular value decomposition of A = U w V_transpose
{
  // where A is m x n real matrix using Householder bidiagonalization
  // and a variant of the QR algorithm, on return A is overwritten by U
  // to get A back compute R[i][j] += A[i][k]*V[j][k]*s[k] over k,j,i
  // based on Golub & Reinsch algorithm as written by Forsythe et al. (1977)
  // converted to C++ in the manner of NRinC but with zero-based indices
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
  int flag,i,its,j,jj,k,l,nm;
  double anorm,c,f,g,h,s,scale,x,y,z,*rv1,sum;

  rv1 = (double*) malloc(n*sizeof(double));
  g = 0.0; scale = 0.0; anorm = 0.0;
  for (i=0;i<n;i++)
    // householder reduction to bidiagonal form...
    {
      l = i+1;
      rv1[i] = scale*g;
      g = 0.0; s = 0.0; scale = 0.0;
      if ( i < m )
	{
	  for (k=i;k<m;k++) scale += fabs(a[k][i]);
	  if ( scale )
	    {
	      for (k=i;k<m;k++)
		{
		  a[k][i] /= scale;
		  s += a[k][i]*a[k][i];
		}
	      f = a[i][i];
	      g = -SIGN(sqrt(s),f);
	      h = f*g-s;
	      a[i][i] = f-g;
	      for (j=l;j<n;j++)
		{
		  for (s=0.0,k=i;k<m;k++) s += a[k][i]*a[k][j];
		  f = s/h;
		  for (k=i;k<m;k++) a[k][j] += f*a[k][i];
		}
	      for (k=i;k<m;k++) a[k][i] *= scale;
	    }
	}
      w[i] = scale*g;
      g = 0.0;
      s = 0.0;
      scale = 0.0;
      if ( (i < m) && (i != n-1) )
	{
	  for (k=l;k<n;k++) scale += fabs(a[i][k]);
	  if ( scale )
	    {
	      for (k=l;k<n;k++)
		{
		  a[i][k] /= scale;
		  s += a[i][k]*a[i][k];
		}
	      f = a[i][l];
	      g = -SIGN(sqrt(s),f);
	      h = f*g-s;
	      a[i][l] = f-g;
	      for (k=l;k<n;k++) rv1[k] = a[i][k]/h;
	      for (j=l;j<m;j++)
		{
		  for (s=0.0,k=l;k<n;k++) s += a[j][k]*a[i][k];
		  for (k=l;k<n;k++) a[j][k] += s*rv1[k];
		}
	      for (k=l;k<n;k++) a[i][k] *= scale;
	    }
	}
      sum = fabs(w[i])+fabs(rv1[i]);
      if (sum > anorm)
	anorm = sum;
    }

  for (i=n-1;i>=0;i--)
    // accumulation of right-hand transformations
    {
      if (i < n-1)
	{
	  if ( g )
	    {
	      for (j=l;j<n;j++)
		v[j][i] = (a[i][j]/a[i][l])/g;
	      // double division avoids underflow
	      for (j=l;j<n;j++)
		{
		  for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
		  for (k=l;k<n;k++) v[k][j] += s*v[k][i];
		}
	    }
	  for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
	}
      v[i][i] = 1.0;
      g = rv1[i];
      l = i;
    }

  for (i=((m < n) ? m-1 : n-1);i>=0;i--)
    // accumluation of left-hand transformations
    {
      l = i+1;
      g = w[i];
      for (j=l;j<n;j++) a[i][j] = 0.0;
      if ( g )
	{
	  g = 1.0/g;
	  // invert to avoid repeated division
	  for (j=l;j<n;j++)
	    {
	      for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
	      f = (s/a[i][i])*g;
	      // double division avoids underflow
	      for (k=i;k<m;k++) a[k][j] += f*a[k][i];
	    }
	  for (j=i;j<m;j++) a[j][i] *= g;
	}
      else for (j=i;j<m;j++) a[j][i] = 0.0;
      ++a[i][i];
    }

  for (k=n-1;k>=0;k--)
    // diagonalization of bidiagonal form
    {
      for (its=1;its<=30;its++)
	// give up after 30 iterations if not converging
	{
	  flag = 1;
	  for (l=k;l>=0;l--)
	    // test for splitting
	    {
	      nm = l-1;
	      if ((double)(fabs(rv1[l])+anorm) == anorm)
		{
		  flag = 0;
		  // rv1(1) is always zero, so there is no exit
		  break;
		}
	      if ((double)(fabs(w[nm])+anorm) == anorm) break;
	    }
	  if ( flag )
	    // cancellation of rv1[l] if l > 0
	    {
	      c = 0.0;
	      s = 1.0;
	      for (i=l;i<=k;i++)
		{
		  f = s*rv1[i];
		  rv1[i] = c*rv1[i];
		  if ((double)(fabs(f)+anorm) == anorm) break;
		  g = w[i];
		  h = pythag(f,g);
		  w[i] = h;
		  h = 1.0/h;
		  c = g*h;
		  s = -f*h;
		  for (j=0;j<m;j++)
		    {
		      y = a[j][nm];
		      z = a[j][i];
		      a[j][nm] = y*c+z*s;
		      a[j][i] = z*c-y*s;
		    }
		}
	    }
	  // test for convergence
	  z = w[k];
	  if ( l == k )
	    // shift from bottom 2x2 minor
	    {
	      if ( z < 0.0 )
		{
		  w[k] = -z;
		  for (j=0;j<n;j++) v[j][k] = -v[j][k];
		}
	      break;
	    }
	  if ( its == 30 )
	    {
	      fprintf(stderr, "Convergence failure in svd()\n");
	      break;
	    }
	  x = w[l];
	  nm = k-1;
	  y = w[nm];
	  g = rv1[nm];
	  h = rv1[k];
	  f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
	  g = pythag(f,1.0);
	  f = ((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
	  c = 1.0; s = 1.0;
	  // next QR transformation
	  for (j=l;j<=nm;j++)
	    {
	      i = j+1;
	      g = rv1[i];
	      y = w[i];
	      h = s*g;
	      g = c*g;
	      z = pythag(f,h);
	      rv1[j] = z;
	      c = f/z;
	      s = h/z;
	      f = x*c+g*s;
	      g = g*c-x*s;
	      h = y*s;
	      y *= c;
	      for (jj=0;jj<n;jj++)
		{
		  x = v[jj][j];
		  z = v[jj][i];
		  v[jj][j] = x*c+z*s;
		  v[jj][i] = z*c-x*s;
		}
	      z = pythag(f,h);
	      w[j] = z;
	      // rotation can be arbitrary of z is zero
	      if ( z )
		{
		  z = 1.0/z;
		  c = f*z;
		  s = h*z;
		}
	      f = c*g+s*y;
	      x = c*y-s*g;
	      for (jj=0;jj<m;jj++)
		{
		  y = a[jj][j];
		  z = a[jj][i];
		  a[jj][j] = y*c+z*s;
		  a[jj][i] = z*c-y*s;
		}
	    }
	  rv1[l] = 0.0;
	  rv1[k] = f;
	  w[k] = x;
	}
    }
  free(rv1);
#undef SIGN
}

double
pythag (double a, double b)
{
  double absa, absb, r;
  absa = fabs(a);
  absb = fabs(b);
  if ( absa > absb )
    {
      r = absb/absa;
      if ( r == 0.0 ) return( absa );
      else return( absa*sqrt(1.0+r*r) );
    }
  else
    {
      r = absa/absb;
      if ( r == 0.0 ) return( absb );
      else return( absb*sqrt(1.0+r*r) );
    }
}

