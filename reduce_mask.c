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
  unsigned long long n;
  int iw, ih;
  unsigned char *mask;
  int x, y;
  int dx, dy;
  int i;
  unsigned char inputName[PATH_MAX];
  unsigned char outputName[PATH_MAX];
  char msg[PATH_MAX + 256];
  unsigned char *out = NULL;
  int factor;
  int ow, oh;
  unsigned long long imbpl, ombpl;
  int error;

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

      fprintf(stderr, "Usage: reduce_mask -input in.pbm -output out.pbm\n");
      fprintf(stderr, "                   -factor int_value\n");
      exit(1);
    }

  if (factor <= 0)
    {
      fprintf(stderr, "reduction factor must be a positive integer\n");
      exit(1);
    }

  mask = NULL;
  if (!ReadBitmap(inputName, &mask,
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

  if (!WriteBitmap(outputName, out, ow, oh, UncompressedBitmap, msg))
    {
      fprintf(stderr, "Could not open %s for writing\n", outputName);
      exit(1);
    }
  free(out);
}
