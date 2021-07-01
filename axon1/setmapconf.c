#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include "imio.h"

int
main (int argc, char **argv)
{
  double angle, tx, ty;
  MapElement *map;
  char imName0[PATH_MAX], imName1[PATH_MAX];
  char msg[PATH_MAX+256];
  double xv, yv;
  double xvp, yvp;
  double factor;
  int level;
  int mw, mh;
  int mxMin, myMin;
  int x, y;
  double cost, sint;
  double conf;
    
  if (argc != 4 ||
      sscanf(argv[3], "%lf", &conf) != 1)
    {
      fprintf(stderr, "Usage: setmapconf input.map output.map confidence\n");
      exit(1);
    }

  if (!ReadMap(argv[1], &map, &level, &mw, &mh,
	       &mxMin, &myMin, imName0, imName1,
	       msg))
    {
      fprintf(stderr, "Could not read map %s:\n%s\n", argv[1], msg);
      exit(1);
    }

  factor = (double) (1 << level);
  for (y = 0; y < mh; ++y)
    for (x = 0; x < mw; ++x)
      map[y*mw+x].c = conf;

  if (!WriteMap(argv[2], map, level, mw, mh, mxMin, myMin, imName0, imName1,
		UncompressedMap, msg))
    {
      fprintf(stderr, "Could not write map %s:\n%s\n", argv[2], msg);
      exit(1);
    }

  return(0);
}
