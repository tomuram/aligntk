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
  int i;
  char msg[PATH_MAX + 256];
  int mbpl;
  unsigned char *mask = NULL;

  if (argc != 3)
    {
      fprintf(stderr, "Usage: genmask7 input.tif output.pbm\n");
      exit(1);
    }

  img = NULL;
  if (!ReadImage(argv[1], &img,
		 &iw, &ih,
		 -1, -1, -1, -1,
		 msg))
    {
      fprintf(stderr, "%s", msg);
      exit(1);
    }
  mbpl = (iw + 7) / 8;
  mask = (unsigned char *) malloc(ih * mbpl);
  memset(mask, 0, ih * mbpl);
  for (y = 0; y < ih; ++y)
    for (x = 0; x < iw; ++x)
      if (img[y*iw+x] != 0)
	mask[y * mbpl + (x >> 3)] |= 0x80 >> (x & 7);
  free(img);

  if (!WriteBitmap(argv[2], mask, iw, ih, UncompressedBitmap, msg))
    {
      fprintf(stderr, "Could not open %s for writing\n", argv[2]);
      exit(1);
    }
  free(mask);
  return(0);
}
