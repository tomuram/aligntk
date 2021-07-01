#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <errno.h>

#include "imio.h"
#include "invert.h"

void Error (char *fmt, ...);

int
main (int argc, char **argv)
{
  int width, height;
  int power;
  float extraC;

  char mapName[PATH_MAX];
  MapElement *map;
  int mLevel;
  int mw, mh;
  int mxMin, myMin;
  char imgName[PATH_MAX], refName[PATH_MAX];
  
  char outputMapName[PATH_MAX];
  int ow, oh;
  int oxMin, oyMin;
  MapElement *omap;

  int factor;
  int x, y;
  float x1, y1, c1;
  float xv, yv;
  float xp, yp;
  int ix, iy;
  float rx, ry, rc;
  float rx00, rx01, rx10, rx11;
  float ry00, ry01, ry10, ry11;
  float rc00, rc01, rc10, rc11;
  float rrx, rry;

  int dx, dy;
  double weight;
  double totalWeight;
  double trx, try;

  int i;
  int error;
  char msg[PATH_MAX+256];

  error = 0;
  mapName[0] = '\0';
  outputMapName[0] = '\0';
  width = -1;
  height = -1;
  power = 4;
  extraC = 0.0;
  for (i = 1; i < argc; ++i)
    if (strcmp(argv[i], "-map") == 0)
      {
	if (++i == argc)
          {
            error = 1;
            break;
          }
        strcpy(mapName, argv[i]);
      }
    else if (strcmp(argv[i], "-output") == 0)
      {
        if (++i == argc)
          {
            error = 1;
            break;
          }
        strcpy(outputMapName, argv[i]);
      }
    else if (strcmp(argv[i], "-width") == 0)
      {
        if (++i == argc || sscanf(argv[i], "%d", &width) != 1)
          {
            error = 1;
            break;
          }
      }
    else if (strcmp(argv[i], "-height") == 0)
      {
        if (++i == argc || sscanf(argv[i], "%d", &height) != 1)
          {
            error = 1;
            break;
          }
      }
    else if (strcmp(argv[i], "-power") == 0)
      {
        if (++i == argc || sscanf(argv[i], "%d", &power) != 1)
          {
            error = 1;
            break;
          }
      }
    else if (strcmp(argv[i], "-extra_c") == 0)
      {
        if (++i == argc || sscanf(argv[i], "%f", &extraC) != 1)
          {
            error = 1;
            break;
          }
      }
    else
      {
	error = 1;
	break;
      }

  if (error)
    {
      if (i < argc)
        fprintf(stderr, "Invalid option: %s\n", argv[i]);
      else
        fprintf(stderr, "Incomplete option: %s\n", argv[i-1]);
      fprintf(stderr, "\n");

      fprintf(stderr, "Usage: extrapolate_map -map <map_name>\n");
      fprintf(stderr, "            -output <output_map_name>]\n");
      fprintf(stderr, "            -extra_c <float>\n");
      fprintf(stderr, "            [-width <int>]\n");
      fprintf(stderr, "            [-height <int>]\n");
      fprintf(stderr, "            [-power <int (range [0-4])>]\n");
      exit(1);
    }

  /* check that at least minimal parameters were supplied */
  if (mapName[0] == '\0')
    Error("-map parameter must be specified.\n");
  if (outputMapName[0] == '\0')
    Error("-output parameter must be specified.\n");
  if (extraC <= 0.0)
    Error("-extra_c must be specified and greater than 0.\n");

  if (!ReadMap(mapName, &map, &mLevel, &mw, &mh,
	       &mxMin, &myMin, imgName, refName, msg))
    Error("Could not read map %s:\n%s\n", mapName, msg);

  factor = 1 << mLevel;
  ow = mw;
  oh = mh;
  oxMin = mxMin;
  oyMin = myMin;
  if (width >= 0)
    {
      ow = ((width + factor - 1) / factor) + 1;
      oxMin = 0;
    }
  if (height >= 0)
    {
      oh = ((height + factor - 1) / factor) + 1;
      oyMin = 0;
    }
  omap = (MapElement*) malloc(ow*oh*sizeof(MapElement));

  // extrapolate all entries with 0 confidence
  for (y = 0; y < oh; ++y)
    for (x = 0; x < ow; ++x)
      {
	if ((x + oxMin) >= mxMin && (x + oxMin) < mxMin+mw &&
	    (y + oyMin) >= myMin && (y + oyMin) < myMin+mh &&
	    map[(y+oyMin-myMin)*mw + (x+oxMin-mxMin)].c != 0.0)
	  {
	    // no need to extrapolate
	    omap[y*ow+x].x = map[(y+oyMin-myMin)*mw + (x+oxMin-mxMin)].x;
	    omap[y*ow+x].y = map[(y+oyMin-myMin)*mw + (x+oxMin-mxMin)].y;
	    omap[y*ow+x].c = map[(y+oyMin-myMin)*mw + (x+oxMin-mxMin)].c;
	    continue;
	  }

	// use a filter when doing the extrapolation
	trx = 0.0;
	try = 0.0;
	totalWeight = 0.0;
	for (iy = 0; iy < mh-1; ++iy)
	  {
	    if (iy >= (y + oyMin - myMin))
	      dy = iy - (y + oyMin - myMin);
	    else
	      dy = (y + oyMin - myMin) - (iy + 1);
	    rry = (y + oyMin - myMin) - iy;

	    for (ix = 0; ix < mw-1; ++ix)
	      {
		if (map[iy*mw+ix].c == 0.0 ||
		    map[(iy+1)*mw+ix].c == 0.0 ||
		    map[iy*mw+ix+1].c == 0.0 ||
		    map[(iy+1)*mw+ix+1].c == 0.0)
		  continue;
		    
		if (ix >= (x + oxMin - mxMin))
		  dx = ix - (x + oxMin - mxMin);
		else
		  dx = (x + oxMin - mxMin) - (ix + 1);
		rrx = (x + oxMin - mxMin) - ix;

		rx00 = map[iy*mw+ix].x;
		ry00 = map[iy*mw+ix].y;
		rx01 = map[(iy+1)*mw+ix].x;
		ry01 = map[(iy+1)*mw+ix].y;
		rx10 = map[iy*mw+ix+1].x;
		ry10 = map[iy*mw+ix+1].y;
		rx11 = map[(iy+1)*mw+ix+1].x;
		ry11 = map[(iy+1)*mw+ix+1].y;

		switch (power)
		  {
		  case 0:
		    weight = 1.0;
		    break;
		  case 1:
		    weight = 1.0 / sqrt((double) (dx*dx + dy*dy));
		    break;
		  case 2:
		    weight = 1.0 / ((double) (dx*dx + dy*dy));
		    break;
		  case 3:
		    weight = 1.0 / ((double) (dx*dx + dy*dy)) /
		      sqrt((double) (dx*dx + dy*dy));
		    break;
		  case 4:
		    weight = 1.0 / ((double) (dx*dx + dy*dy)) /
		      ((double) (dx*dx + dy*dy));
		    break;
		  default:
		    Error("power = %d is unsupported.\n", power);
		    break;
		  }
		if (isnan(weight))
		  Error("weight is nan at (%d,%d)\n", x, y);
		if (isinf(weight))
		  Error("weight is inf at (%d,%d)\n", x, y);
		    
		trx += weight * (rx00 * (rrx - 1.0) * (rry - 1.0)
				 - rx10 * rrx * (rry - 1.0)
				 - rx01 * (rrx - 1.0) * rry
				 + rx11 * rrx * rry);
		try += weight * (ry00 * (rrx - 1.0) * (rry - 1.0)
				 - ry10 * rrx * (rry - 1.0)
				 - ry01 * (rrx - 1.0) * rry
				 + ry11 * rrx * rry);
		
		totalWeight += weight;
	      }
	  }
	if (totalWeight == 0.0)
	  Error("Could not extrapolate map element (%d,%d).\n", x, y);

	omap[y*ow+x].x = trx / totalWeight;
	omap[y*ow+x].y = try / totalWeight;
	omap[y*ow+x].c = extraC;
      }

  if (!WriteMap(outputMapName, omap, mLevel, ow, oh,
		oxMin, oyMin,
		imgName, refName,
		UncompressedMap, msg))
    Error("Could not write map %s:\n%s\n",
	  outputMapName, msg);

  free(map);
  free(omap);
  return(0);
}

void Error (char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
  fflush(stderr);
  abort();
}
