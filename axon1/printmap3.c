#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <limits.h>
#include <math.h>

#include "imio.h"

void Error (char *fmt, ...);

int
main (int argc, char **argv)
{
  int mw, mh;
  int x, y;
  MapElement *map;
  int level;
  int xMin, yMin;
  char imgn[PATH_MAX], refn[PATH_MAX];
  char msg[PATH_MAX+1024];
  double xv, yv;
  double rrx, rry;
  double rx00, ry00, rc00;
  double rx01, ry01, rc01;
  double rx10, ry10, rc10;
  double rx11, ry11, rc11;
  double rx, ry, rc;

  if (argc != 4)
    Error("Usage: printmap <mapfile> <x> <y>\n");
  if (sscanf(argv[2], "%lf", &xv) != 1 ||
      sscanf(argv[3], "%lf", &yv) != 1)
    Error("Usage: printmap <mapfile> <x> <y>\n");
  if (!ReadMap(argv[1], &map,
	      &level, &mw, &mh, &xMin, &yMin,
	      imgn, refn,
	       msg))
    {
      fprintf(stderr, "%s\n", msg);
      exit(1);
    }
	      
  x = (int) floor(xv);
  y = (int) floor(yv);
  if (x < 0 || x >= mw)
    Error("x value is out-of-range\n");
  if (y < 0 || y >= mh)
    Error("y value is out-of-range\n");
  printf("image = %s  ref = %s  map level = %d\n", imgn, refn, level);
  printf("mw = %d  mh = %d\n", mw, mh);
  printf("xMin = %d  yMin = %d\n", xMin, yMin);
  rrx = xv - x;
  rry = yv - y; 
  if (rrx == 0.0 && rry == 0.0)
    printf("[%d, %d] maps with confidence %f to [%f, %f]\n",
	   x, y,
	   map[y*mw+x].c,
	   map[y*mw+x].x,
	   map[y*mw+x].y);
  else
    {
      if (x < 0 || x >= mw-1)
	Error("x value is out-of-range\n");
      if (y < 0 || y >= mh-1)
	Error("y value is out-of-range\n");

      rx00 = map[y*mw+x].x;
      ry00 = map[y*mw+x].y;
      rc00 = map[y*mw+x].c;
      rx01 = map[(y+1)*mw+x].x;
      ry01 = map[(y+1)*mw+x].y;
      rc01 = map[(y+1)*mw+x].c;
      rx10 = map[y*mw + x + 1].x;
      ry10 = map[y*mw + x + 1].y;
      rc10 = map[y*mw + x + 1].c;
      rx11 = map[(y+1)*mw + x + 1].x;
      ry11 = map[(y+1)*mw + x + 1].y;
      rc11 = map[(y+1)*mw + x + 1].c;

      rx = rx00 * (rrx - 1.0) * (rry - 1.0)
	- rx10 * rrx * (rry - 1.0) 
	- rx01 * (rrx - 1.0) * rry
	+ rx11 * rrx * rry;
      ry = ry00 * (rrx - 1.0) * (rry - 1.0)
	- ry10 * rrx * (rry - 1.0) 
	- ry01 * (rrx - 1.0) * rry
	+ ry11 * rrx * rry;
      rc = rc00 * (rrx - 1.0) * (rry - 1.0)
	- rc10 * rrx * (rry - 1.0) 
	- rc01 * (rrx - 1.0) * rry
	+ rc11 * rrx * rry;

      printf("[%f, %f] maps with confidence %f to [%f, %f]\n",
	     xv, yv,
	     rc, rx, ry);
    }
  return(0);
}

void Error (char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
  exit(1);
}
