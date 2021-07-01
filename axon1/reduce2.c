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
  int startSection, endSection;
  char ddirName[PATH_MAX];
  char sdirName[PATH_MAX];
  int section;
  char fn[PATH_MAX];
  FILE *f;
  int n;
  int iw, ih;
  unsigned char *img;
  int x, y;
  int dx, dy;
  int i;
  int iv;
  char msg[PATH_MAX + 256];
  unsigned char *out = NULL;
  float v;
  int factor;
  int ow, oh;

  if (argc != 5 ||
      strlen(argv[1]) < 1 ||
      sscanf(argv[2], "%d", &startSection) != 1 ||
      sscanf(argv[3], "%d", &endSection) != 1 ||
      sscanf(argv[4], "%d", &factor) != 1)
    {
      fprintf(stderr, "Usage: reduce2 dataset_dir start_section end_section reduction_factor\n");
      exit(1);
    }

  if (chdir(argv[1]) != 0)
    {
      fprintf(stderr, "Could not change to directory %s\n", argv[1]);
      exit(1);
    }

  if (getcwd(ddirName, PATH_MAX) == NULL)
    {
      fprintf(stderr, "Could not get pathname of working directory\n");
      exit(1);
    }

  for (section = startSection; section <= endSection; ++section)
    {
      sprintf(fn, "%s/sections8/images/z%0.4d.tif", ddirName, section);
      img = NULL;
      if (!ReadImage(fn, &img,
		     &iw, &ih,
		     -1, -1, -1, -1,
		     msg))
	{
	  fprintf(stderr, "%s", msg);
	  exit(1);
	}
      ow = iw / factor;
      oh = ih / factor;
      n = ow * oh;
      out = (unsigned char *) malloc(n * sizeof(unsigned char));
      memset(out, 0, n);
	    
      for (y = 0; y < oh; ++y)
	for (x = 0; x < ow; ++x)
	  {
	    v = 0.0;
	    for (dy = 0; dy < factor; ++dy)
	      for (dx = 0; dx < factor; ++dx)
		v += img[(factor*y+dy)*iw + factor*x+dx];
	    v /= factor * factor;
	    iv = (int) floor(v + 0.5);
	    if (iv < 0)
	      iv = 0;
	    else if (iv > 255)
	      iv = 255;
	    out[y * ow + x] = iv;
	  }
      free(img);
      img = NULL;

      sprintf(fn, "%s/sections%d/z%0.4d.tif",
	      ddirName, 8*factor, section);
      if (!WriteImage(fn, out, ow, oh, UncompressedImage, msg))
	{
	  fprintf(stderr, "Could not open %s for writing\n", fn);
	  exit(1);
	}
      free(out);
    }
}
