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
    
  if (argc != 6 ||
      sscanf(argv[3], "%lf", &angle) != 1 ||
      sscanf(argv[4], "%lf", &tx) != 1 ||
      sscanf(argv[5], "%lf", &ty) != 1)
    {
      fprintf(stderr, "Usage: rotate_map input.map output.map rotation_in_degrees x_translation y_translation\n");
      exit(1);
    }

  if (!ReadMap(argv[1], &map, &level, &mw, &mh,
	       &mxMin, &myMin, imName0, imName1,
	       msg))
    {
      fprintf(stderr, "Could not read map %s:\n%s\n", argv[1], msg);
      exit(1);
    }

  if (angle == 0.0 || angle == -360.0 || angle == 360.0)
    {
      cost = 1.0;
      sint = 0.0;
    }
  else if (angle == 90.0 || angle == -270.0)
    {
      cost = 0.0;
      sint = 1.0;
    }
  else if (angle == 180.0 || angle == -180.0)
    {
      cost = -1.0;
      sint = 0.0;
    }
  else if (angle == 270.0 || angle == -90.0)
    {
      cost = 0.0;
      sint = -1.0;
    }
  else
    {
      cost = cos(angle * M_PI / 180.0);
      sint = sin(angle * M_PI / 180.0);
    }

  factor = (double) (1 << level);
  for (y = 0; y < mh; ++y)
    for (x = 0; x < mw; ++x)
      {
	if (map[y*mw+x].c == 0.0)
	  continue;

	xv = map[y*mw+x].x;
	yv = map[y*mw+x].y;

	xvp = cost * xv + sint * yv + tx / factor;
	yvp = -sint * xv + cost * yv + ty / factor;

	map[y*mw+x].x = xvp;
	map[y*mw+x].y = yvp;
      }

  if (!WriteMap(argv[2], map, level, mw, mh, mxMin, myMin, imName0, imName1,
		UncompressedMap, msg))
    {
      fprintf(stderr, "Could not write map %s:\n%s\n", argv[2], msg);
      exit(1);
    }

  return(0);
}
