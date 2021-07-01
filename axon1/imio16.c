#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>
#include <tiffio.h>
#include <stdarg.h>
#include <limits.h>

#include "imio16.h"

int ReadHeader16 (FILE *f, char *tc, uint32 *w, uint32 *h, int *m);

int
ReadImages16 (char *filename, int nImages,
	      unsigned short **pixels,
	      int *width, int *height,
	      char *error)
{
  int len;
  uint32 iw, ih;
  unsigned short *buffer;
  unsigned long bufferSize, count;
  FILE *f;
  char tc;
  int m;
  int w, h;
  int y;
  int ni;
  int i;

  len = strlen(filename);
  if (len < 5)
    {
      sprintf(error, "Image filename is too short: %s\n", filename);
      return(0);
    }
  if (strcmp(&filename[len-4], ".pgm") == 0)
    {
      f = fopen(filename, "r");
      if (f == NULL)
	{
	  sprintf(error, "Could not open file %s for reading\n",
		  filename);
	  return(0);
	}
      for (ni = 0; ni < nImages; ++ni)
	{
	  if (!ReadHeader16(f, &tc, &iw, &ih, &m))
	    {
	      sprintf(error, "Image file %s not binary pgm.  (%c %d %d %d)\n", filename,
		      tc, iw, ih, m);
	      return(0);
	    }
	  if ((buffer = (unsigned short *) malloc(iw * ih * sizeof(unsigned short))) == NULL)
	    {
	      sprintf(error, "Could not allocate enough memory for image in file %s\n", filename);
	      return(0);
	    }
	  if (fread(buffer, sizeof(unsigned short), iw * ih, f) != iw * ih)
	    {
	      sprintf(error, "Image file %s apparently truncated.\n",
		      filename);
	      return(0);
	    }
	  for (i = 0; i < iw * ih; ++i)
	    buffer[i] = ((buffer[i] & 0xff) << 8) | ((buffer[i] >> 8) & 0xff);
	  width[ni] = iw;
	  height[ni] = ih;
	  pixels[ni] = buffer;
	}
      fclose(f);
      return(1);
    }
  else
    {
      sprintf(error, "Unrecognized file extension for image file %s\n", filename);
      return(0);
    }
}

int
ReadHeader16 (FILE *f, char *tc, unsigned int *w, unsigned int *h, int *m)
{
  int c;
  int v;

  c = fgetc(f);
  while (c == '#')
    {
      while ((c = fgetc(f)) != EOF && c != '\n') ;
    }
  if (c != 'P')
    {
      fprintf(stderr, "ERROR2\n");
      return(0);
    }
  c = fgetc(f);
  if (c != '4' && c != '5' && c != '6')
    {
      fprintf(stderr, "ERROR3\n");
      return(0);
    }
  *tc = c;
  c = fgetc(f);

  while (c == ' ' || c == '\t' || c == '\n' || c == '#')
    if (c == '#')
      {
	while ((c = fgetc(f)) != EOF && c != '\n') ;
      }
    else
      c = fgetc(f);
  if (c >= '0' && c <= '9')
    {
      v = 0;
      while (c != EOF && c >= '0' && c <= '9')
	{
	  v = 10 * v + (c - '0');
	  c = fgetc(f);
	}
      *w = v;
    }
  else
    {
      fprintf(stderr, "ERROR4\n");
      return(0);
    }

  while (c == ' ' || c == '\t' || c == '\n' || c == '#')
    if (c == '#')
      {
	while ((c = fgetc(f)) != EOF && c != '\n') ;
      }
    else
      c = fgetc(f);
  if (c >= '0' && c <= '9')
    {
      v = 0;
      while (c != EOF && c >= '0' && c <= '9')
	{
	  v = 10 * v + (c - '0');
	  c = fgetc(f);
	}
      *h = v;
    }
  else
    {
      fprintf(stderr, "ERROR5\n");
      return(0);
    }

  if (*tc == '4')
    return(1);

  while (c == ' ' || c == '\t' || c == '\n' || c == '#')
    if (c == '#')
      while ((c = fgetc(f)) != EOF && c != '\n') ;
    else
      c = fgetc(f);
  if (c >= '0' && c <= '9')
    {
      v = 0;
      while (c != EOF && c >= '0' && c <= '9')
	{
	  v = 10 * v + (c - '0');
	  c = fgetc(f);
	}
      *m = v;
    }
  else
    {
      fprintf(stderr, "ERROR1\n");
      return(0);
    }

  return(1);
}

