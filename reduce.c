#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <errno.h>
#include <math.h>

#include "imio.h"

int
main (int argc, char **argv)
{
  int iw, ih;
  unsigned char *img;
  int x, y;
  int dx, dy;
  int i;
  int iv;
  unsigned char inputName[PATH_MAX];
  unsigned char outputName[PATH_MAX];
  char msg[PATH_MAX + 256];
  unsigned char *out = NULL;
  float v;
  int factor;
  int ow, oh;
  int error;
  unsigned long long liw, lih;
  unsigned long long low, loh;
  unsigned long long n;

  error = 0;
  inputName[0] = '\0';
  outputName[0] = '\0';
  factor = 0;
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
    else if (strcmp(argv[i], "-factor") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d", &factor) != 1)
	  {
	    fprintf(stderr, "-factor error\n");
	    error = 1;
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
      if (i >= argc)
	fprintf(stderr, "Incomplete option: %s\n\n", argv[i-1]);

      fprintf(stderr, "Usage: reduce -input in.tif -output out.tif\n");
      fprintf(stderr, "              -factor int_value\n");
      exit(1);
    }

  if (factor <= 0)
    {
      fprintf(stderr, "reduction factor must be a positive integer\n");
      exit(1);
    }

  img = NULL;
  if (!ReadImage(inputName, &img,
		 &iw, &ih,
		 -1, -1, -1, -1,
		 msg))
    {
      fprintf(stderr, "%s", msg);
      exit(1);
    }
  ow = iw / factor;
  oh = ih / factor;
  liw = iw;
  lih = ih;
  low = ow;
  loh = oh;
  n = low * loh;
  if (n == 0)
    {
      fprintf(stderr, "Image size was too small to reduce.\n");
      exit(1);
    }
  out = (unsigned char *) malloc(n * sizeof(unsigned char));
  memset(out, 0, n);
	    
  for (y = 0; y < oh; ++y)
    for (x = 0; x < ow; ++x)
      {
	v = 0.0;
	for (dy = 0; dy < factor; ++dy)
	  for (dx = 0; dx < factor; ++dx)
	    v += img[(factor*y+dy)*liw + factor*x+dx];
	v /= factor * factor;
	iv = (int) floor(v + 0.5);
	if (iv < 0)
	  iv = 0;
	else if (iv > 255)
	  iv = 255;
	out[y * low + x] = iv;
      }
  free(img);
  img = NULL;

  if (!WriteImage(outputName, out, ow, oh, UncompressedImage, msg))
    {
      fprintf(stderr, "Could not open %s for writing\n", outputName);
      exit(1);
    }
  free(out);
}
