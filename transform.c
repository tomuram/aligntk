/*
 * transform.cc  - transform images and masks
 *
 *  Copyright (c) 2008-2013 National Resource for Biomedical
 *                          Supercomputing,
 *                          Pittsburgh Supercomputing Center,
 *                          Carnegie Mellon University
 *
 *  Distribution of this code is prohibited without the prior written
 *  consent of the National Resource for Biomedical Supercomputing.
 *
 *  Acknowledgements:
 *     Development of this code was supported in part by
 *       NIH NCRR grant 5P41RR006009 and
 *       NIH NIGMS grant P41GM103712
 *
 *  This program is distributed in the hope that it will
 *  be useful, but WITHOUT ANY WARRANTY; without even the
 *  implied warranty of MERCHANTABILITY or FITNESS FOR A
 *  PARTICULAR PURPOSE.  Neither Carnegie Mellon University
 *  nor any of the authors assume any liability for
 *  damages, incidental or otherwise, caused by the
 *  installation or use of this software.
 *
 *  CLINICAL APPLICATIONS ARE NOT RECOMMENDED, AND THIS
 *  SOFTWARE HAS NOT BEEN EVALUATED BY THE UNITED STATES
 *  FDA FOR ANY CLINICAL USE.
 *
 *  HISTORY
 *    2008-2013  Written by Greg Hood (ghood@psc.edu)
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <dirent.h>

#include "imio.h"

#define SHIFT   0
#define CROP    1
#define ROTATE  2
#define SCALE   3
#define STRETCH 4
#define SIZE    5
#define RESIZE  6

typedef struct Op
{
  int type;
  double p0;
  double p1;
  double p2;
  double p3;
} Op;

typedef struct Constraint
{
  /* a*x + b*y + c >= 0 */
  double a;
  double b;
  double c;
} Constraint;

/* FORWARD DECLARATIONS */
void ApplyTransform(double m[2][3], double t[2][3], int n, Constraint *con);
int Compare (const void *x, const void *y);

int
main (int argc, char **argv)
{
  int i, j;
  int n;
  int len;
  int error;
  int maskWidth, maskHeight;
  int w, h;
  double minX, maxX, minY, maxY;
  double xShift, yShift;
  double theta;
  double width, height;
  double scale;
  double xScale, yScale;
  char inputName[PATH_MAX];
  char maskName[PATH_MAX];
  char outputName[PATH_MAX];
  char outputMaskName[PATH_MAX];
  int oversamplingFactor;
  int nOps;
  Op* ops;
  int v;
  int x, y;
  int dx, dy;
  unsigned char *image;
  unsigned char *mask;
  unsigned char *output;
  unsigned char *outputMask;
  int nConstraints;
  Constraint* constraints;
  int mh, mw;
  double t[2][3];
  double m[2][3];
  int xSize, ySize;
  int maskPresent;
  int mbpl, ombpl;
  int outputWidth, outputHeight;
  int valid;
  int nx, ny;
  int sum;
  double px, py;
  int outside;
  int osf, osf2;
  char msg[PATH_MAX+256];

  error = 0;
  inputName[0] = '\0';
  maskName[0] = '\0';
  outputName[0] = '\0';
  outputMaskName[0] = '\0';
  oversamplingFactor = 16;
  nOps = 0;
  ops = NULL;
  for (i = 1; i < argc; ++i)
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
    else if (strcmp(argv[i], "-mask") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-mask error\n");
	    break;
	  }
	strcpy(maskName, argv[i]);
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
    else if (strcmp(argv[i], "-output_mask") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-output_mask error\n");
	    break;
	  }
	strcpy(outputMaskName, argv[i]);
      }
    else if (strcmp(argv[i], "-shift") == 0)
      {
	if (i+2 >= argc ||
	    sscanf(argv[i+1], "%lf", &xShift) != 1 ||
	    sscanf(argv[i+2], "%lf", &yShift) != 1)
	  {
	    error = 1;
	    fprintf(stderr, "-shift <x_shift> <y_shift> error\n");
	    break;
	  }
	++nOps;
	ops = (Op*) realloc(ops, nOps*sizeof(Op));
	ops[nOps-1].type = SHIFT;
	ops[nOps-1].p0 = xShift;
	ops[nOps-1].p1 = yShift;
	ops[nOps-1].p2 = 0.0;
	ops[nOps-1].p3 = 0.0;
	i += 2;
      }
    else if (strcmp(argv[i], "-crop") == 0)
      {
	if (i+4 >= argc ||
	    sscanf(argv[i+1], "%lf", &minX) != 1 ||
	    sscanf(argv[i+2], "%lf", &maxX) != 1 ||
	    sscanf(argv[i+3], "%lf", &minY) != 1 ||
	    sscanf(argv[i+4], "%lf", &maxY) != 1 ||
	    minX > maxX || minY > maxY)
	  {
	    error = 1;
	    fprintf(stderr, "-crop <x_min> <x_max> <y_min> <y_max> error\n");
	    break;
	  }
	++nOps;
	ops = (Op*) realloc(ops, nOps*sizeof(Op));
	ops[nOps-1].type = CROP;
	ops[nOps-1].p0 = minX;
	ops[nOps-1].p1 = maxX;
	ops[nOps-1].p2 = minY;
	ops[nOps-1].p3 = maxY;
	i += 4;
      }
    else if (strcmp(argv[i], "-rotate") == 0)
      {
	if (i+1 >= argc ||
	    sscanf(argv[i+1], "%lf", &theta) != 1)
	  {
	    error = 1;
	    fprintf(stderr, "-rotate <theta_degrees_ccw> error\n");
	    break;
	  }
	++nOps;
	ops = (Op*) realloc(ops, nOps*sizeof(Op));
	ops[nOps-1].type = ROTATE;
	ops[nOps-1].p0 = theta;
	ops[nOps-1].p1 = 0.0;
	ops[nOps-1].p2 = 0.0;
	ops[nOps-1].p3 = 0.0;
	i += 1;
      }
    else if (strcmp(argv[i], "-scale") == 0)
      {
	if (i+1 > argc ||
	    sscanf(argv[i+1], "%lf", &scale) != 1)
	  {
	    error = 1;
	    fprintf(stderr, "-scale <scale> error\n");
	    break;
	  }
	++nOps;
	ops = (Op*) realloc(ops, nOps*sizeof(Op));
	ops[nOps-1].type = SCALE;
	ops[nOps-1].p0 = scale;
	ops[nOps-1].p1 = 0.0;
	ops[nOps-1].p2 = 0.0;
	ops[nOps-1].p3 = 0.0;
	i += 1;
      }
    else if (strcmp(argv[i], "-stretch") == 0)
      {
	if (i+2 > argc ||
	    sscanf(argv[i+1], "%lf", &xScale) != 1 ||
	    sscanf(argv[i+2], "%lf", &yScale) != 1)
	  {
	    error = 1;
	    fprintf(stderr, "-stretch <x_scale> <y_scale> error\n");
	    break;
	  }
	++nOps;
	ops = (Op*) realloc(ops, nOps*sizeof(Op));
	ops[nOps-1].type = STRETCH;
	ops[nOps-1].p0 = xScale;
	ops[nOps-1].p1 = yScale;
	ops[nOps-1].p2 = 0.0;
	ops[nOps-1].p3 = 0.0;
	i += 2;
      }
    else if (strcmp(argv[i], "-size") == 0)
      {
	if (i+2 > argc ||
	    sscanf(argv[i+1], "%lf", &width) != 1 ||
	    sscanf(argv[i+2], "%lf", &height) != 1 ||
	    width <= 0 || height <= 0)
	  {
	    error = 1;
	    fprintf(stderr, "-size <width> <height> error\n");
	    break;
	  }
	++nOps;
	ops = (Op*) realloc(ops, nOps*sizeof(Op));
	ops[nOps-1].type = SIZE;
	ops[nOps-1].p0 = width;
	ops[nOps-1].p1 = height;
	ops[nOps-1].p2 = 0.0;
	ops[nOps-1].p3 = 0.0;
	i += 2;
      }
    else if (strcmp(argv[i], "-resize") == 0)
      {
	if (i+2 > argc ||
	    sscanf(argv[i+1], "%lf", &width) != 1 ||
	    sscanf(argv[i+2], "%lf", &height) != 1 ||
	    width <= 0 || height <= 0)
	  {
	    error = 1;
	    fprintf(stderr, "-resize <width> <height> error\n");
	    break;
	  }
	++nOps;
	ops = (Op*) realloc(ops, nOps*sizeof(Op));
	ops[nOps-1].type = RESIZE;
	ops[nOps-1].p0 = width;
	ops[nOps-1].p1 = height;
	ops[nOps-1].p2 = 0.0;
	ops[nOps-1].p3 = 0.0;
	i += 2;
      }
    else if (strcmp(argv[i], "-oversampling") == 0)
      {
	if (i+1 > argc ||
	    sscanf(argv[i+1], "%d", &oversamplingFactor) != 1)
	  {
	    error = 1;
	    fprintf(stderr, "-oversampling error\n");
	    break;
	  }
      }
    else
      {
	fprintf(stderr, "Unknown option: %s\n", argv[i]);
	error = 1;
      }
      

  if (error)
    {
      fprintf(stderr, "Usage: transform -image_list images.lst\n");
      fprintf(stderr, "              -input directory_prefix\n");
      fprintf(stderr, "             [-mask directory_prefix]\n");
      fprintf(stderr, "              -output file_prefix\n");
      fprintf(stderr, "             [-output_mask file_prefix]\n");
      fprintf(stderr, "             [-shift shift_x shift_y]\n");
      fprintf(stderr, "             [-crop min_x max_x min_y max_y]\n");
      fprintf(stderr, "             [-rotate theta_degrees_ccw]\n");
      fprintf(stderr, "             [-scale scale]\n");
      fprintf(stderr, "             [-stretch scale_x scale_y]\n");
      fprintf(stderr, "             [-size width height]");
      fprintf(stderr, "             [-resize width height]");
      exit(1);
    }

  /* check that at least minimal parameters were supplied */
  if (inputName[0] == '\0' ||
      outputName[0] == '\0')
    {
      fprintf(stderr, "-input and -output parameters must be specified.\n");
      exit(1);
    }

  /* read the image in */
  if (!ReadImage(inputName, &image, &w, &h,
		 -1, -1, -1, -1,
		 msg))
    {
      fprintf(stderr, "%s\n", msg);
      exit(1);
    }

  nConstraints = 4;
  constraints = (Constraint *) malloc(nConstraints * sizeof(Constraint));
  constraints[0].a = 1.0;
  constraints[0].b = 0.0;
  constraints[0].c = 0.0;
  constraints[1].a = -1.0;
  constraints[1].b = 0.0;
  constraints[1].c = w;
  constraints[2].a = 0.0;
  constraints[2].b = 1.0;
  constraints[2].c = 0.0;
  constraints[3].a = 0.0;
  constraints[3].b = -1.0;
  constraints[3].c = h;
  t[0][0] = 1.0;
  t[0][1] = 0.0;
  t[0][2] = 0.0;
  t[1][0] = 0.0;
  t[1][1] = 1.0;
  t[1][2] = 0.0;
  xSize = w;
  ySize = h;

  /* read in the mask, if one exists */
  maskPresent = 0;
  if (maskName[0] != '\0')
    {
      if (!ReadBitmap(maskName, &mask,
		      &mw, &mh,
		      -1, -1, -1, -1,
		      msg))
	{
	  fprintf(stderr, "%s\n", msg);
	  exit(1);
	}

      if (mw != w || mh != h)
	{
	  fprintf(stderr, "Inconsistent mask size.\n");
	  exit(1);
	}
      mbpl = (w + 7) >> 3;
      maskPresent = 1;
    }
  if (!maskPresent)
    {
      mbpl = (w + 7) >> 3;
      mask = (unsigned char *) malloc(h * mbpl);
      memset(mask, 0xff, h*mbpl);
    }

  osf = 1;
  for (i = 0; i < nOps; ++i)
    switch(ops[i].type)
      {
      case SHIFT:
	if (ops[i].p0 != floor(ops[i].p0) ||
	    ops[i].p1 != floor(ops[i].p1))
	  osf = oversamplingFactor;
	m[0][0] = 1.0;
	m[0][1] = 0.0;
	m[0][2] = ops[i].p0;
	m[1][0] = 0.0;
	m[1][1] = 1.0;
	m[1][2] = ops[i].p1;
	ApplyTransform(m, t, nConstraints, constraints);
	break;

      case CROP:
	if (ops[i].p0 != floor(ops[i].p0) ||
	    ops[i].p1 != floor(ops[i].p1) ||
	    ops[i].p2 != floor(ops[i].p2) ||
	    ops[i].p3 != floor(ops[i].p3))
	  osf = oversamplingFactor;
	nConstraints += 4;
	constraints = (Constraint *) realloc(constraints,
					     nConstraints * sizeof(Constraint));
	constraints[nConstraints-4].a = 1.0;
	constraints[nConstraints-4].b = 0.0;
	constraints[nConstraints-4].c = -ops[i].p0;
	constraints[nConstraints-3].a = -1.0;
	constraints[nConstraints-3].b = 0.0;
	constraints[nConstraints-3].c = ops[i].p1;
	constraints[nConstraints-2].a = 0.0;
	constraints[nConstraints-2].b = 1.0;
	constraints[nConstraints-2].c = -ops[i].p2;
	constraints[nConstraints-1].a = 0.0;
	constraints[nConstraints-1].b = -1.0;
	constraints[nConstraints-1].c = ops[i].p3;
	break;

      case ROTATE:
	if (fmod(ops[i].p0, 90.0) != 0.0)
	  osf = oversamplingFactor;
	m[0][0] = cos(ops[i].p0 * M_PI / 180.0);
	m[0][1] = sin(ops[i].p0 * M_PI / 180.0);
	m[0][2] = 0.0;
	m[1][0] = -sin(ops[i].p0 * M_PI / 180.0);
	m[1][1] = cos(ops[i].p0 * M_PI / 180.0);
	m[1][2] = 0.0;
	ApplyTransform(m, t, nConstraints, constraints);
	break;

      case SCALE:
	if (ops[i].p0 < 1.0 || ops[i].p0 != floor(ops[i].p0))
	  osf = oversamplingFactor;
	m[0][0] = ops[i].p0;
	m[0][1] = 0.0;
	m[0][2] = 0.0;
	m[1][0] = 0.0;
	m[1][1] = ops[i].p0;
	m[1][2] = 0.0;
	ApplyTransform(m, t, nConstraints, constraints);
	break;

      case STRETCH:
	if (ops[i].p0 < 1.0 || ops[i].p1 < 1.0 ||
	    ops[i].p0 != floor(ops[i].p0) || ops[i].p1 != floor(ops[i].p1))
	  osf = oversamplingFactor;
	m[0][0] = ops[i].p0;
	m[0][1] = 0.0;
	m[0][2] = 0.0;
	m[1][0] = 0.0;
	m[1][1] = ops[i].p1;
	m[1][2] = 0.0;
	ApplyTransform(m, t, nConstraints, constraints);
	break;

      case SIZE:
	if (ops[i].p0 != floor(ops[i].p0) ||
	    ops[i].p1 != floor(ops[i].p1))
	  osf = oversamplingFactor;
	nConstraints += 4;
	constraints = (Constraint *) realloc(constraints,
					     nConstraints * sizeof(Constraint));
	constraints[nConstraints-4].a = 1.0;
	constraints[nConstraints-4].b = 0.0;
	constraints[nConstraints-4].c = 0.0;
	constraints[nConstraints-3].a = -1.0;
	constraints[nConstraints-3].b = 0.0;
	constraints[nConstraints-3].c = ops[i].p0;
	constraints[nConstraints-2].a = 0.0;
	constraints[nConstraints-2].b = 1.0;
	constraints[nConstraints-2].c = 0.0;
	constraints[nConstraints-1].a = 0.0;
	constraints[nConstraints-1].b = -1.0;
	constraints[nConstraints-1].c = ops[i].p1;
	xSize = floor(ops[i].p0);
	ySize = floor(ops[i].p1);
	break;

      case RESIZE:
	if (ops[i].p0 != floor(ops[i].p0) ||
	    ops[i].p1 != floor(ops[i].p1) ||
	    ops[i].p0 < xSize ||
	    ops[i].p1 < ySize ||
	    fmod(ops[i].p0, (double) xSize) != 0.0 ||
	    fmod(ops[i].p1, (double) ySize) != 0.0)
	  osf = oversamplingFactor;
	m[0][0] = ops[i].p0 / xSize;
	m[0][1] = 0.0;
	m[0][2] = 0.0;
	m[1][0] = 0.0;
	m[1][1] = ops[i].p1 / ySize;
	m[1][2] = 0.0;
	ApplyTransform(m, t, nConstraints, constraints);
	nConstraints += 4;
	constraints = (Constraint *) realloc(constraints,
					     nConstraints * sizeof(Constraint));
	constraints[nConstraints-4].a = 1.0;
	constraints[nConstraints-4].b = 0.0;
	constraints[nConstraints-4].c = 0.0;
	constraints[nConstraints-3].a = -1.0;
	constraints[nConstraints-3].b = 0.0;
	constraints[nConstraints-3].c = ops[i].p0;
	constraints[nConstraints-2].a = 0.0;
	constraints[nConstraints-2].b = 1.0;
	constraints[nConstraints-2].c = 0.0;
	constraints[nConstraints-1].a = 0.0;
	constraints[nConstraints-1].b = -1.0;
	constraints[nConstraints-1].c = ops[i].p1;
	xSize = floor(ops[i].p0);
	ySize = floor(ops[i].p1);
	break;
      }
	  
  outputWidth = xSize;
  outputHeight = ySize;
  output = (unsigned char *) malloc(outputHeight*outputWidth);
  ombpl = (outputWidth + 7) / 8;
  outputMask = (unsigned char *) malloc(outputHeight*ombpl);
  memset(outputMask, 0, outputHeight*ombpl);

  /* compute the output image */
  osf2 = osf * osf;
  for (y = 0; y < outputHeight; ++y)
    for (x = 0; x < outputWidth; ++x)
      {
	sum = 0;
	valid = 1;
	for (dy = 0; dy < osf; ++dy)
	  {
	    py = y + (((double) dy) + 0.5) / osf;
	    for (dx = 0; dx < osf; ++dx)
	      {
		px = x + (((double) dy) + 0.5) / osf;
		outside = 0;
		for (i = 0; i < nConstraints; ++i)
		  if (constraints[i].a * px +
		      constraints[i].b * py +
		      constraints[i].c <= 0.0)
		    {
		      outside = 1;
		      valid = 0;
		      break;
		    }
		if (!outside)
		  {
		    nx = floor(t[0][0] * px + t[0][1] * py + t[0][2]);
		    ny = floor(t[1][0] * px + t[1][1] * py + t[1][2]);
		    if (nx >= 0 && nx < w &&
			ny >= 0 && ny < h)
		      {
			sum += image[ny*w + nx];
			if ((mask[ny*mbpl + (nx >> 3)] & (0x80 >> (nx & 7))) == 0)
			  valid = 0;
		      }
		    else
		      valid = 0;
		  }
	      }
	  }
	output[y*outputWidth+x] = (sum + (osf2 >> 1)) / osf2;
	if (valid)
	  outputMask[y*ombpl + (x >> 3)] |= 0x80 >> (x & 7);
      }

  /* write out the output image */
  if (!WriteImage(outputName, output,
		  outputWidth, outputHeight,
		  UncompressedImage,
		  msg))
    {
      fprintf(stderr, "%s\n", msg);
      exit(1);
    }
    
  /* write out the output mask */
  if (outputMaskName[0] != '\0')
    {
      if (!WriteBitmap(outputMaskName, outputMask,
		       outputWidth, outputHeight,
		       UncompressedBitmap,
		       msg))
	{
	  fprintf(stderr, "%s\n", msg);
	  exit(1);
	}
    }

  free(image);
  free(mask);
  free(output);
  free(outputMask);
  free(constraints);
  
  exit(0);
}

void
ApplyTransform(double m[2][3], double t[2][3], int n, Constraint *con)
{
  int i, j;
  double det;
  double tmp[2][3];
  double inv[2][3];

  det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
  /*  printf("det = %f\n", det); */
  inv[0][0] = m[1][1] / det;
  inv[0][1] = -m[0][1] / det;
  inv[0][2] = (m[0][1] * m[1][2] - m[1][1] * m[0][2]) / det;
  inv[1][0] = -m[1][0] / det;
  inv[1][1] = m[0][0] / det;
  inv[1][2] = (m[0][2] * m[1][0] - m[1][2] * m[0][0]) / det;

  tmp[0][0] = t[0][0] * inv[0][0] + t[0][1] * inv[1][0];
  tmp[0][1] = t[0][0] * inv[0][1] + t[0][1] * inv[1][1];
  tmp[0][2] = t[0][0] * inv[0][2] + t[0][1] * inv[1][2] + t[0][2];
  tmp[1][0] = t[1][0] * inv[0][0] + t[1][1] * inv[1][0];
  tmp[1][1] = t[1][0] * inv[0][1] + t[1][1] * inv[1][1];
  tmp[1][2] = t[1][0] * inv[0][2] + t[1][1] * inv[1][2] + t[1][2];
  for (i = 0; i < 2; ++i)
    for (j = 0; j < 3; ++j)
      t[i][j] = tmp[i][j];

  /*  printf("TRANSFORM:  %f %f %f\n", t[0][0], t[0][1], t[0][2]);
      printf("            %f %f %f\n", t[1][0], t[1][1], t[1][2]); */

  for (i = 0; i < n; ++i)
    {
      tmp[0][0] = con[i].a * inv[0][0] + con[i].b * inv[1][0];
      tmp[0][1] = con[i].a * inv[0][1] + con[i].b * inv[1][1];
      tmp[0][2] = con[i].a * inv[0][2] + con[i].b * inv[1][2] + con[i].c;
      con[i].a = tmp[0][0];
      con[i].b = tmp[0][1];
      con[i].c = tmp[0][2];
      /*      printf("CONSTRAINT %d:  %f * x + %f *y + %f > 0\n", i,
	      con[i].a, con[i].b, con[i].c); */
    }
}

