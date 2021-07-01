/*
 *  merge_images.c  -  merges multiple co-registered images into a single image
 *
 *  Copyright (c) 2012-2013 National Resource for Biomedical
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
 *       NIH NIGMS grant P41GM103712
 *
 *  HISTORY
 *    2012  Written by Greg Hood (ghood@psc.edu)
 */

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "imio.h"

#define MAX_INPUTS	32

int
main (int argc, char **argv)
{
  int i;
  char msg[PATH_MAX+256];
  char *inputName[MAX_INPUTS];
  unsigned char *input;
  int iw, ih;
  char outputName[PATH_MAX];
  unsigned char *output;
  int ow, oh;
  int x, y;
  int nInputs;
  int error;

  error = 0;
  nInputs = 0;
  for (i = 1; i < argc; ++i)
    if (strcmp(argv[i], "-input") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	if (nInputs >= MAX_INPUTS)
	  {
	    fprintf(stderr, "Maximum number of -input arguments exceeeded.\n");
	    exit(1);
	  }
	inputName[nInputs] = (char *) malloc(strlen(argv[i]) + 1);
	strcpy(inputName[nInputs], argv[i]);
      }
    else if (strcmp(argv[i], "-output") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(outputName, argv[i]);
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
      
      fprintf(stderr, "Usage: merge_images -input image0.tif\n");
      fprintf(stderr, "                    -input image1.tif\n");
      fprintf(stderr, "                    ...\n");
      fprintf(stderr, "                    -output out.tif\n");
      exit(1);
    }
	
  fprintf(stderr, "Merging images ");
  for (i = 0; i < nInputs; ++i)
    fprintf(stderr, "%s ", inputName[i]);
  fprintf(stderr, "into %s\n", outputName);
  fflush(stderr);

  for (i = 0; i < nInputs; ++i)
    {
      if (!ReadImage(inputName[i], &input, &iw, &ih, -1, -1, -1, -1, msg))
	{
	  fprintf(stderr, "Could not read image %s:\n%s\n", inputName[i], msg);
	  exit(1);
	}
      if (i == 0)
	{
	  ow = iw;
	  oh = ih;
	  output = (unsigned char *) malloc(oh*ow);
	  memset(output, 0, oh*ow);
	}
      else if (iw != ow || ih != oh)
	{
	  fprintf(stderr, "Image %s has a different size.\n", inputName[i]);
	  exit(1);
	}

      for (y = 0; y < ih; ++y)
	for (x = 0; x < iw; ++x)
	  if (output[y*ow+x] == 0 &&
	      input[y*iw+x] != 0)
	    output[y*ow+x] = input[y*iw+x];
      free(input);
      input = NULL;
    }

  if (!WriteImage(outputName, output, ow, oh,
		  UncompressedImage, msg))
    {
      fprintf(stderr, "Could not write image %s:\n%s\n",
	      outputName, msg);
      exit(1);
    }

  free(output);
  return(0);
}
