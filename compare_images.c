#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <stdarg.h>

#include "imio.h"

int
main (int argc, char **argv)
{
  FILE *f;
  int iw1, ih1;
  int iw2, ih2;
  int i;
  int error;
  int x, y;
  float d;
  double sd2;
  double rmsd;
  char image1Name[PATH_MAX];
  char image2Name[PATH_MAX];
  char outputName[PATH_MAX];
  unsigned char *image1;
  unsigned char *image2;
  char errorMsg[PATH_MAX+256];

  error = 0;
  image1Name[0] = '\0';
  image2Name[0] = '\0';
  outputName[0] = '\0';
  for (i = 1; i < argc; ++i)
    if (strcmp(argv[i], "-image1") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-image1 error\n");
	    break;
	  }
	strcpy(image1Name, argv[i]);
      }
    else if (strcmp(argv[i], "-image2") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-image2 error\n");
	    break;
	  }
	strcpy(image2Name, argv[i]);
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
    else
      error = 1;

  if (argc == 1 || error)
    {
      fprintf(stderr, "Usage: compare_images -image1 image1.tif\n");
      fprintf(stderr, "                      -image2 image2.tif\n");
      fprintf(stderr, "                       [-output rms.out]\n");
      exit(1);
    }

  /* check that at least minimal parameters were supplied */
  if (image1Name[0] == '\0' || image2Name[0] == '\0')
    {
      fprintf(stderr, "-image1 and -image2 parameters must be specified.\n");
      exit(1);
    }

  if (!ReadImage(image1Name, &image1, &iw1, &ih1, -1, -1, -1, -1, errorMsg))
    {
      fprintf(stderr, "Error reading image %s:\n  %s\n", image1Name, errorMsg);
      exit(1);
    }
  if (!ReadImage(image2Name, &image2, &iw2, &ih2, -1, -1, -1, -1, errorMsg))
    {
      fprintf(stderr, "Error reading image %s:\n  %s\n", image2Name, errorMsg);
      exit(1);
    }

  /* make sure width and height match */
  if (iw2 != iw1 || ih2 != ih1)
    {
      fprintf(stderr, "Image dimensions are not consistent.\n");
      exit(1);
    }


  sd2 = 0.0;
  for (y = 0; y < ih1; ++y)
    for (x = 0; x < iw1; ++x)
      {
	d = image1[y*iw1+x] - image2[y*iw2+x];
	sd2 += d*d;
      }
  rmsd = sqrt(sd2 / (ih1 * iw1) );

  if (outputName[0] != '\0')
    {
      f = fopen(outputName, "w");
      fprintf(f, "%lf\n", rmsd);
      fclose(f);
    }
  printf("RMS difference = %lf\n", rmsd);
  return(0);
}
