#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdarg.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include "reader.h"

void hsv_to_rgb (unsigned char* r, unsigned char *g, unsigned char *b,
		 float h, float s, float v);
void DrawLine(unsigned char *img, float px0, float py0, float px1, float py1,
	      float nomD, int imageWidth, int imageHeight);

int
main (int argc, char **argv)
{
  FILE *f;
  int x, y, z;
  int i, j, k;
  char inputName[PATH_MAX];
  char outputName[PATH_MAX];
  int dx, dy;
  int xp, yp;
  unsigned char red, green, blue;
  MapElement *map;
  int error;
  int level;
  int mapWidth, mapHeight;
  int xMin, yMin;
  char imgName[PATH_MAX], refName[PATH_MAX];
  char errorMsg[PATH_MAX + 256];
  int imageWidth, imageHeight;
  int mapFactor;
  int mapOffsetX, mapOffsetY;
  unsigned char *grid;

  error = 0;
  inputName[0] = '\0';
  outputName[0] = '\0';

  for (i = 1; i < argc; ++i)
    if (strcmp(argv[i], "-input") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(inputName, argv[i]);
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
    else error = 1;

  if (error)
    {
      fprintf(stderr, "Usage: plotmap -input input.map -output grid.pnm\n");
      exit(1);
    }
      
  /* check that at least minimal parameters were supplied */
  if (inputName[0] == '\0')
    {
      fprintf(stderr, "-input must be specified\n");
      exit(1);
    }
  if (outputName[0] == '\0')
    {
      fprintf(stderr, "-output must be specified.\n");
      exit(1);
    }

  if (!ReadMap(inputName, &map,
	       &level, &mapWidth, &mapHeight,
	       &xMin, &yMin,
	       imgName, refName,
	       errorMsg))
    {
      fprintf(stderr, "Error reading map %s:\n  %s\n", inputName, errorMsg);
      exit(1);
    }

  imageWidth = 1024;
  imageHeight = 1024;
  mapOffsetX = 16;
  mapOffsetY = 16;
  mapFactor = (imageWidth - mapOffsetX) / mapWidth;
  if ((imageHeight - mapOffsetY) / mapHeight < mapFactor)
    mapFactor = (imageHeight - mapOffsetX) / mapHeight;

  grid = (unsigned char *) malloc(imageWidth * imageHeight * 3);
  memset(grid, 0, imageWidth * imageHeight * 3);
  for (y = 0; y < mapHeight; ++y)
    for (x = 0; x < mapWidth; ++x)
      for (dy = -1; dy <= 0; ++dy)
	{
	  yp = y + dy;
	  if (yp < 0)
	    continue;
	  for (dx = -1; dx <= 1; ++dx)
	    {
	      xp = x + dx;
	      if (xp < 0)
		continue;
	      if (xp >= mapWidth || dy == 0 && dx >= 0)
		break;
	      
	      DrawLine(grid,
		       mapFactor*map[y*mapWidth+x].x + mapOffsetX,
		       mapFactor*map[y*mapWidth+x].y + mapOffsetY,
		       mapFactor*map[yp*mapWidth+xp].x + mapOffsetX,
		       mapFactor*map[yp*mapWidth+xp].y + mapOffsetY,
		       ((dx * dy != 0) ? M_SQRT2 : 1.0),
		       imageWidth, imageHeight);
	    }
	}
  
  f = fopen(outputName, "w");
  if (f == NULL)
    {
      fprintf(stderr, "Could not open output grid file %s\n", outputName);
      exit(1);
    }
  fprintf(f, "P6\n%d %d\n%d\n", imageWidth, imageHeight, 255);
  if (fwrite(grid, sizeof(unsigned char),
	     3*imageWidth*imageHeight, f) !=
      3*imageWidth*imageHeight)
    {
      fprintf(stderr, "Could not write to grid file %s\n", outputName);
      exit(1);
    }
  fclose(f);
  free(grid);
  return(0);
}

void hsv_to_rgb (unsigned char* r, unsigned char *g, unsigned char *b,
		 float h, float s, float v)
{ 
  /* H is given on [0, 6] or UNDEFINED. S and V are given on [0, 1]. 
     RGB are each returned on [0, 255]. */
#define RETURN_RGB(rd,gr,bl)	{ *r = floor(255.999999 * rd); *g = floor(255.999999 * gr); *b = floor(255.999999 * bl); return; }
  
  float m, n, f;
  int i; 
  
  i = floor(h);
  f = h - i; 
  if ( !(i&1) )
    f = 1 - f; /* if i is even */
  m = v * (1 - s); 
  n = v * (1 - s * f); 
  
  switch (i)
    { 
    case 6: 
    case 0: RETURN_RGB(v, n, m); 
    case 1: RETURN_RGB(n, v, m); 
    case 2: RETURN_RGB(m, v, n);
    case 3: RETURN_RGB(m, n, v); 
    case 4: RETURN_RGB(n, m, v); 
    case 5: RETURN_RGB(v, m, n); 
    } 
} 

void
DrawLine(unsigned char *img, float px0, float py0, float px1, float py1,
	 float nomD, int imageWidth, int imageHeight)
{
  int x0, y0, x1, y1;
  int steep;
  int tmp;
  int error;
  int yStep;
  int x, y;
  int dx, dy;
  float dist;
  unsigned char r, g, b;
  float hue;
  unsigned char *p;

  x0 = (int) (px0 + 0.5);
  y0 = (int) (py0 + 0.5);
  x1 = (int) (px1 + 0.5);
  y1 = (int) (py1 + 0.5);

  dist = hypot(x1 - x0, y1 - y0);
  hue = 4.0 * (dist / nomD) - 2.0;
  if (hue < 0.0)
    hue = 0.0;
  if (hue > 5.0)
    hue = 5.0;
  hsv_to_rgb(&r, &g, &b, hue, 1.0, 1.0);

  // use the Bresenham Algorithm for lines
  //    (derived from the Wikipedia entry on
  //     "Bresenham's line algorithm")
  steep = abs(y1 - y0) > abs(x1 - x0);
  if (steep)
    {
      tmp = x0;
      x0 = y0;
      y0 = tmp;
      tmp = x1;
      x1 = y1;
      y1 = tmp;
    }
  if (x0 > x1)
    {
      tmp = x0;
      x0 = x1;
      x1 = tmp;
      tmp = y0;
      y0 = y1;
      y1 = tmp;
    }
  dx = x1 - x0;
  dy = abs(y1 - y0);
  error = -dx / 2;
  y = y0;
  if (y0 < y1)
    yStep = 1;
  else
    yStep = -1;
  for (x = x0; x <= x1; ++x)
    {
      if (steep)
	{
	  if (y >= 0 && y < imageWidth &&
	      x >= 0 && x < imageHeight)
	    {
	      p = &img[3*(x * imageWidth + y)];
	      *p++ = r;
	      *p++ = g;
	      *p++ = b;
	    }
	}
      else
	{
	  if (x >= 0 && x < imageWidth &&
	      y >= 0 && y < imageHeight)
	    {
	      p = &img[3*(y * imageWidth + x)];
	      *p++ = r;
	      *p++ = g;
	      *p++ = b;
	    }
	}
      error += dy;
      if (error > 0)
	{
	  y += yStep;
	  error -= dx;
	}
    }
}
