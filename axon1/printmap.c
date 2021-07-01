#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

void Error (char *fmt, ...);

int
main (int argc, char **argv)
{
  FILE *f;
  int mapHeader[3];
  int mw, mh;
  int mw2, mh2;
  float *xMap, *yMap, *cMap;
  int x, y;

  if (argc != 4)
    Error("Usage: printmap <mapfile> <x> <y>\n");
  f = fopen(argv[1], "r");
  if (f == NULL)
    Error("Could not open map file %s\n", argv[1]);
  if (fread(mapHeader, sizeof(int), 3, f) != 3)
    Error("Could not read map header from file %s", argv[1]);
  mw = mapHeader[0];
  mh = mapHeader[1];
  mw2 = mw + 2;
  mh2 = mh + 2;
  xMap = (float *) malloc(mw2 * mh2 * sizeof(float));
  yMap = (float *) malloc(mw2 * mh2 * sizeof(float));
  cMap = (float *) malloc(mw2 * mh2 * sizeof(float));
  if (fread(xMap, sizeof(float), mw2 * mh2, f) != mw2 * mh2 ||
      fread(yMap, sizeof(float), mw2 * mh2, f) != mw2 * mh2 ||
      fread(cMap, sizeof(float), mw2 * mh2, f) != mw2 * mh2)
    Error("Could not read map from file %s\n", argv[1]);
  fclose(f);
  
  if (sscanf(argv[2], "%d", &x) != 1 ||
      sscanf(argv[3], "%d", &y) != 1)
    Error("Usage: printmap <mapfile> <x> <y>\n");
  if (x < -1 || x > mw)
    Error("x value is out-of-range\n");
  if (y < -1 || y > mh)
    Error("y value is out-of-range\n");
  printf("[%d, %d] maps with confidence %f to [%f, %f] in slice %d\n",
	 x, y,
	 cMap[(y + 1) * mw2 + (x + 1)],
	 xMap[(y + 1) * mw2 + (x + 1)], yMap[(y + 1) * mw2 + (x + 1)],
	 mapHeader[2]);
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
