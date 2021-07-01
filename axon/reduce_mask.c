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
  int n;
  int iw, ih;
  unsigned char *mask;
  int x, y;
  int dx, dy;
  int i;
  char msg[PATH_MAX + 256];
  unsigned char *out = NULL;
  int factor;
  int ow, oh;
  int imbpl, ombpl;

  if (argc != 4 ||
      sscanf(argv[1], "%d", &factor) != 1)
    {
      fprintf(stderr, "Usage: reduce_mask reduction_factor input.pbm output.pbm\n");
      exit(1);
    }

  mask = NULL;
  if (!ReadBitmap(argv[2], &mask,
		 &iw, &ih,
		 -1, -1, -1, -1,
		 msg))
    {
      fprintf(stderr, "%s", msg);
      exit(1);
    }
  imbpl = (iw + 7) >> 3;
  ow = iw / factor;
  ombpl = (ow + 7) >> 3;
  oh = ih / factor;
  n = ombpl * oh;
  out = (unsigned char *) malloc(n * sizeof(unsigned char));
  memset(out, 0, n);
	    
  for (y = 0; y < oh; ++y)
    for (x = 0; x < ow; ++x)
      {
	for (dy = 0; dy < factor; ++dy)
	  for (dx = 0; dx < factor; ++dx)
	    if (mask[(factor*y+dy)*imbpl + ((factor*x+dx) >> 3)] &
		(0x80 >> ((factor*x+dx) & 7)))
	      {
		out[y*ombpl+ (x>>3)] |= 0x80 >> (x & 7);
		goto outerLoop;
	      }
      outerLoop: ;
      }
  free(mask);
  mask = NULL;

  if (!WriteBitmap(argv[3], out, ow, oh, UncompressedBitmap, msg))
    {
      fprintf(stderr, "Could not open %s for writing\n", argv[3]);
      exit(1);
    }
  free(out);
}
