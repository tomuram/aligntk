#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

#include "imio.h"
#include "invert.h"

int
main (int argc, char **argv)
{
  int error;
  int i;
  int and;
  int invert;
  char mapName[PATH_MAX];
  char inputName[PATH_MAX];
  char modelName[PATH_MAX];
  char outputName[PATH_MAX];
  MapElement *map;
  int mLevel;
  int mw, mh;
  int mxMin, myMin;
  char imName0[PATH_MAX], imName1[PATH_MAX];
  char msg[PATH_MAX+256];
  float factor;
  unsigned char *inputMask;
  int iw, ih;
  int imbpl;
  unsigned char *modelMask;
  int ow, oh;
  int ombpl;
  unsigned char *outputMask;
  int x, y;
  float xv, yv;
  int ixv, iyv;
  float rx, ry;
  int irx, iry;
  float rrx, rry;
  float rx00, ry00, rc00;
  float rx01, ry01, rc01;
  float rx10, ry10, rc10;
  float rx11, ry11, rc11;
  InverseMap *inverseMap;

  error = 0;
  mapName[0] = '\0';
  inputName[0] = '\0';
  modelName[0] = '\0';
  outputName[0] = '\0';
  and = 0;
  invert = 0;
  for (i = 1; i < argc; ++i)
    if (strcmp(argv[i], "-invert") == 0)
      invert = 1;
    else if (strcmp(argv[i], "-map") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(mapName, argv[i]);
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
      error = 1;

  if (error ||
      mapName[0] == '\0' ||
      inputName[0] == '\0' ||
      modelName[0] == '\0' ||
      outputName[0] == '\0')
    {
      fprintf(stderr, "Usage: applymap2mask\n");
      fprintf(stderr, "         -map map_name\n");
      fprintf(stderr, "         -input mask_name\n");
      fprintf(stderr, "         -model mask_name\n");
      fprintf(stderr, "         -output mask_name\n");
      fprintf(stderr, "         [-invert]\n");
      fprintf(stderr, "         [-and]\n");
      exit(1);
    }
  
  if (!ReadMap(mapName, &map, &mLevel, &mw, &mh,
	       &mxMin, &myMin, imName0, imName1, msg))
    {
      fprintf(stderr, "Could not read map %s:\n  error: %s\n",
	      mapName, msg);
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

  if (invert)
    inverseMap = InvertMap(map, mw, mh);

  factor = 1 << mLevel;
  for (y = 0; y < oh; ++y)
    for (x = 0; x < ow; ++x)
      {
	xv = (x + 0.5) / factor;
	yv = (y + 0.5) / factor;
	if (invert)
	  {
	    if (!Invert(inverseMap, &rx, &ry, xv, yv))
	      continue;
	    rx = (rx + mxMin) * factor;
	    ry = (ry + myMin) * factor;
	    irx = (int) floor(rx);
	    iry = (int) floor(ry);
	    if (irx < 0 || irx >= iw ||
		iry < 0 || iry >= ih)
	      continue;

	    if (inputMask[iry * imbpl + (irx >> 3)] & (0x80 >> (irx & 7)))
	      outputMask[y*ombpl + (x >> 3)] |= 0x80 >> (x & 7);
	  }
	else
	  {
	    ixv = (int) floor(xv);
	    iyv = (int) floor(yv);
	    rrx = xv - ixv;
	    rry = yv - iyv;
	    ixv -= mxMin;
	    iyv -= myMin;
	    if (ixv < 0 || ixv >= mw-1 ||
		iyv < 0 || iyv >= mh-1)
	      continue;

	    rx00 = map[iyv*mw+ixv].x;
	    ry00 = map[iyv*mw+ixv].y;
	    rc00 = map[iyv*mw+ixv].c;
	    rx01 = map[(iyv+1)*mw+ixv].x;
	    ry01 = map[(iyv+1)*mw+ixv].y;
	    rc01 = map[(iyv+1)*mw+ixv].c;
	    rx10 = map[iyv*mw+ixv+1].x;
	    ry10 = map[iyv*mw+ixv+1].y;
	    rc10 = map[iyv*mw+ixv+1].c;
	    rx11 = map[(iyv+1)*mw+ixv+1].x;
	    ry11 = map[(iyv+1)*mw+ixv+1].y;
	    rc11 = map[(iyv+1)*mw+ixv+1].c;

	    if (rc00 == 0.0 ||
		rc01 == 0.0 ||
		rc10 == 0.0 ||
		rc11 == 0.0)
	      continue;

	    rx = rx00 * (rrx - 1.0) * (rry - 1.0)
	      - rx10 * rrx * (rry - 1.0) 
	      - rx01 * (rrx - 1.0) * rry
	      + rx11 * rrx * rry;
	    ry = ry00 * (rrx - 1.0) * (rry - 1.0)
	      - ry10 * rrx * (rry - 1.0) 
	      - ry01 * (rrx - 1.0) * rry
	      + ry11 * rrx * rry;
	    irx = (int) floor(rx * factor);
	    iry = (int) floor(ry * factor);
	    if (irx < 0 || irx >= iw ||
		iry < 0 || iry >= ih)
	      continue;

	    if (inputMask[iry * imbpl + (irx >> 3)] & (0x80 >> (irx & 7)))
	      outputMask[y*ombpl + (x >> 3)] |= 0x80 >> (x & 7);
	  }
	
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
