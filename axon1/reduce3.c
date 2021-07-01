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
  int ow, oh;
  unsigned long long liw, lih;
  unsigned long long low, loh;
  unsigned long long n;
  unsigned char *img;
  int x, y;
  int dx, dy;
  int i;
  int iv;
  char msg[PATH_MAX + 256];
  unsigned char *out = NULL;
  float v;
  int factor;

  if (argc != 4 ||
      sscanf(argv[1], "%d", &factor) != 1)
    {
      fprintf(stderr, "Usage: reduce3 reduction_factor input.tif output.tif\n");
      exit(1);
    }

  img = NULL;
  if (!ReadImage(argv[2], &img,
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

  if (!WriteImage(argv[3], out, ow, oh, UncompressedImage, msg))
    {
      fprintf(stderr, "Could not open %s for writing\n", argv[3]);
      exit(1);
    }
  free(out);
}
