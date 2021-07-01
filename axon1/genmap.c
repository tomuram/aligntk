#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include "imio.h"

int
main (int argc, char **argv)
{
  MapElement map[4];
  char name[PATH_MAX];
  char fn[PATH_MAX];
  char msg[PATH_MAX];
  int level;
  double angle;
  double scale;
  double tx;
  double ty;
  int x, y;
  double xv, yv;
  double xp, yp;
  int mFactor;
  double a;

  if (argc != 7 ||
      sscanf(argv[2], "%d", &level) != 1 ||
      sscanf(argv[3], "%lf", &angle) != 1 ||
      sscanf(argv[4], "%lf", &scale) != 1 ||
      sscanf(argv[5], "%lf", &tx) != 1 ||
      sscanf(argv[6], "%lf", &ty) != 1)
    {
      fprintf(stderr, "Usage: gen_map map_name level rotation_angle_in_degrees scale tx ty\n");
      exit(1);
    }

  a = angle * M_PI / 180.0;

  mFactor = 1 << level;
  for (y = 0; y < 2; ++y)
    for (x = 0; x < 2; ++x)
      {
	xv = mFactor * x;
	yv = mFactor * y;
	xp = scale * (cos(a) * xv - sin(a) * yv) + tx;
	yp = scale * (sin(a) * xv + cos(a) * yv) + ty;
	map[y*2+x].x = xp / mFactor;
	map[y*2+x].y = yp / mFactor;
	map[y*2+x].c = 1.0;
      }

  sprintf(name, "rot%+f\n", angle);
  if (!WriteMap(argv[1], map, level,
		2, 2, 0, 0,
		name, "unit_frame",
		UncompressedMap,
		msg))
    {
      fprintf(stderr, "WriteMap error: %s\n", msg);
      exit(1);
    }
  return(0);
}
