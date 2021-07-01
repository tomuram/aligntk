#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
  int error;
  char inputMaskName1[PATH_MAX];
  char inputMaskName2[PATH_MAX];
  char outputMaskName[PATH_MAX];
  int orOperation;
  int andOperation;

  error = 0;
  inputMaskName1[0] = '\0';
  inputMaskName2[0] = '\0';
  outputMaskName[0] = '\0';
  orOperation = 0;
  andOperation = 0;
  for (i = 1; i < argc; ++i)
    if (strcmp(argv[i], "-input1") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(inputMaskName1, argv[i]);
      }
    else if (strcmp(argv[i], "-input2") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(inputMaskName2, argv[i]);
      }
    else if (strcmp(argv[i], "-output") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(outputMaskName, argv[i]);
      }
    else if (strcmp(argv[i], "-or") == 0)
      orOperation = 1;
    else if (strcmp(argv[i], "-and") == 0)
      andOperation = 1;
    else error = 1;
	
  if (error)
    {
      fprintf(stderr, "Usage: combine_masks -input1 mask1.pbm\n");
      fprintf(stderr, "                     -input2 mask2.pbm\n");
      fprintf(stderr, "                     -output out.pbm\n");
      fprintf(stderr, "                     [-and]\n");
      fprintf(stderr, "                     [-or]\n");
      exit(1);
    }

  if ((andOperation ^ orOperation) == 0)
    {
      fprintf(stderr, "Either -and or -or operation must be specified.\n");
      exit(1);
    }
  if (inputMaskName1[0] == '\0' ||
      inputMaskName2[0] == '\0' ||
      outputMaskName[0] == '\0')
    {
      fprintf(stderr, "-input1, -input2, and -output options must be specified.\n");
      exit(1);
    }
  
  if (!ReadBitmap(inputMaskName1, &inputMask1, &iw, &ih,
		  -1, -1, -1, -1,
		  msg))
    {
      fprintf(stderr, "Could not read input mask %s:\n  error: %s\n",
	      argv[1], msg);
      exit(1);
    }
  imbpl = (iw + 7) >> 3;
  if (!ReadBitmap(inputMaskName2, &inputMask2, &iw2, &ih2,
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

  if (orOperation)
    {
      for (y = 0; y < ih; ++y)
	for (i = 0; i < imbpl; ++i)
	  outputMask[y*imbpl + i] = inputMask1[y*imbpl + i] | inputMask2[y*imbpl + i];
    }
  else
    {
      for (y = 0; y < ih; ++y)
	for (i = 0; i < imbpl; ++i)
	  outputMask[y*imbpl + i] = inputMask1[y*imbpl + i] & inputMask2[y*imbpl + i];
    }

  if (!WriteBitmap(outputMaskName, outputMask, iw, ih,
		   UncompressedBitmap, msg))
    {
      fprintf(stderr, "Could not write output mask %s:\n  error: %s\n",
	      argv[3], msg);
      exit(1);
    }

  return(0);
}

