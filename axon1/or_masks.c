#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "imio.h"

int
main (int argc, char **argv)
{
  int i;
  char msg[PATH_MAX+256];
  unsigned char *inputMask1;
  int iw, ih;
  int imbpl;
  unsigned char *inputMask2;
  int iw2, ih2;
  unsigned char *outputMask;
  int x, y;

  if (!ReadBitmap(argv[1], &inputMask1, &iw, &ih,
		  -1, -1, -1, -1,
		  msg))
    {
      fprintf(stderr, "Could not read input mask %s:\n  error: %s\n",
	      argv[1], msg);
      exit(1);
    }
  imbpl = (iw + 7) >> 3;
  if (!ReadBitmap(argv[2], &inputMask2, &iw2, &ih2,
		  -1, -1, -1, -1,
		  msg))
    {
      fprintf(stderr, "Could not read 2nd input mask %s:\n  error: %s\n",
	      argv[2], msg);
      exit(1);
    }
  if (iw2 != iw || ih2 != ih)
    {
      fprintf(stderr, "Sizes of input masks do not match.\n");
      exit(1);
    }
  outputMask = malloc(ih * imbpl * sizeof(unsigned char));
  memset(outputMask, 0, ih * imbpl);

  for (y = 0; y < ih; ++y)
    for (i = 0; i < imbpl; ++i)
      outputMask[y*imbpl + i] = inputMask1[y*imbpl + i] | inputMask2[y*imbpl + i];

  if (!WriteBitmap(argv[3], outputMask, iw, ih,
		   UncompressedBitmap, msg))
    {
      fprintf(stderr, "Could not write output mask %s:\n  error: %s\n",
	      argv[3], msg);
      exit(1);
    }

  return(0);
}

