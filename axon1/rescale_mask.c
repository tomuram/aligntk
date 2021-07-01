#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "imio.h"

int
main (int argc, char **argv)
{
  int error;
  int i;
  int and;
  char inputName[PATH_MAX];
  char modelName[PATH_MAX];
  char outputName[PATH_MAX];
  char msg[PATH_MAX+256];
  int factor;
  unsigned char *inputMask;
  int iw, ih;
  int imbpl;
  unsigned char *modelMask;
  int ow, oh;
  int ombpl;
  unsigned char *outputMask;
  int x, y;
  int ix, iy;

  inputName[0] = '\0';
  modelName[0] = '\0';
  outputName[0] = '\0';
  and = 0;
  factor = 0;
  for (i = 1; i < argc; ++i)
    if (strcmp(argv[i], "-factor") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d", &factor) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-input") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(inputName, argv[i]);
      }
    else if (strcmp(argv[i], "-model") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(modelName, argv[i]);
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
    else if (strcmp(argv[i], "-and") == 0)
      and = 1;
    else
      {
	fprintf(stderr, "Unrecognized argument: %s\n", argv[i]);
	error = 1;
	break;
      }

  if (error ||
      inputName[0] == '\0' ||
      modelName[0] == '\0' ||
      outputName[0] == '\0')
    {
      fprintf(stderr, "Usage: rescale_mask\n");
      fprintf(stderr, "         -factor upscale_factor\n");
      fprintf(stderr, "         -input mask_name\n");
      fprintf(stderr, "         -model mask_name\n");
      fprintf(stderr, "         -output mask_name\n");
      exit(1);
    }
  
  if (!ReadBitmap(inputName, &inputMask, &iw, &ih,
		  -1, -1, -1, -1,
		  msg))
    {
      fprintf(stderr, "Could not read input mask %s:\n  error: %s\n",
	      inputName, msg);
      exit(1);
    }
  imbpl = (iw + 7) >> 3;
  if (!ReadBitmap(modelName, &modelMask, &ow, &oh,
		  -1, -1, -1, -1,
		  msg))
    {
      fprintf(stderr, "Could not read model mask %s:\n  error: %s\n",
	      inputName, msg);
      exit(1);
    }
  ombpl = (ow + 7) >> 3;
  outputMask = malloc(oh * ombpl * sizeof(unsigned char));
  memset(outputMask, 0, oh * ombpl);

  for (y = 0; y < oh; ++y)
    for (x = 0; x < ow; ++x)
      {
	ix = x / factor;
	iy = y / factor;
	if (ix >= iw || iy >= ih)
	  continue;
	if (inputMask[iy * imbpl + (ix >> 3)] & (0x80 >> (ix & 7)))
	  outputMask[y*ombpl + (x >> 3)] |= 0x80 >> (x & 7);
      }

  if (and)
    for (y = 0; y < oh; ++y)
      for (i = 0; i < ombpl; ++i)
	outputMask[y*ombpl + i] &= modelMask[y*ombpl + i];

  if (!WriteBitmap(outputName, outputMask, ow, oh,
		   UncompressedBitmap, msg))
    {
      fprintf(stderr, "Could not write output mask %s:\n  error: %s\n",
	      outputName, msg);
      exit(1);
    }

  return(0);
}
