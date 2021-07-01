#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>
#include <tiffio.h>
#include <stdarg.h>
#include <limits.h>

#include "reader.h"

int ReadHeader (FILE *f, char *tc, uint32 *w, uint32 *h, int *m);

int
ReadImage (char *filename, unsigned char **pixels,
	   int *width, int *height,
	   int minX, int maxX, int minY, int maxY,
	   char *error)
{
  int len;
  TIFF *image;
  uint32 iw, ih;
  uint16 photo, bps, spp, fillorder;
  tsize_t stripSize;
  unsigned long imageOffset, result;
  int stripMax, stripCount;
  unsigned char *buffer;
  unsigned long bufferSize, count;
  FILE *f;
  char tc;
  int m;
  unsigned char *ptr;
  uint32 xmin, xmax, ymin, ymax;
  int w, h;
  int y;

  len = strlen(filename);
  if (len < 5)
    {
      sprintf(error, "Image filename is too short: %s\n", filename);
      return(0);
    }
  if (strcasecmp(&filename[len-4], ".tif") == 0 ||
      len > 5 && strcasecmp(&filename[len-5], ".tiff") == 0)
    {
      // Open the TIFF image
      if ((image = TIFFOpen(filename, "r")) == NULL)
	{
	  sprintf(error, "Could not open TIFF image: %s\n", filename);
	  return(0);
	}

      // Check that it is of a type that we support
      if (TIFFGetField(image, TIFFTAG_BITSPERSAMPLE, &bps) == 0 || bps != 8)
	{
	  sprintf(error, "Either undefined or unsupported number of bits per sample (bps = %d) in tiff image %s\n", bps, filename);
	  return(0);
	}
      
      if (TIFFGetField(image, TIFFTAG_SAMPLESPERPIXEL, &spp) == 0 || spp != 1)
	{
	  sprintf(error, "Either undefined or unsupported number of samples per pixel (spp = %d) in tiff image %s\n", spp, filename);
	  return(0);
	}

      if (TIFFGetField(image, TIFFTAG_IMAGEWIDTH, &iw) == 0)
	{
	  sprintf(error, "TIFF image %s does not define its width\n",
		  filename);
	  return(0);
	}
  
      if (TIFFGetField(image, TIFFTAG_IMAGELENGTH, &ih) == 0)
	{
	  sprintf(error, "Image does not define its height (length)\n",
		  filename);
	  return(0);
	}
  
      // Read in the possibly multiple strips
      stripSize = TIFFStripSize(image);
      stripMax = TIFFNumberOfStrips(image);
      imageOffset = 0;

      bufferSize = stripMax * stripSize;
      //      Log("TIFF: %ld %ld %ld %d %d %d %d\n",
      //	  (long) stripSize, (long) stripMax, (long) bufferSize,
      //	  (int) bps, (int) spp, (int) iw, (int) ih);

      if ((buffer = (unsigned char *) malloc(bufferSize)) == NULL)
        {
          sprintf(error, "Could not allocate enough memory for the uncompressed image (bytes = %d) of TIFF file %s\n", bufferSize, filename);
          return(0);
        }

      for (stripCount = 0; stripCount < stripMax; ++stripCount)
        {
          if ((result = TIFFReadEncodedStrip(image, stripCount,
                                             buffer + imageOffset,
                                             stripSize)) == -1)
            {
              sprintf(error, "Read error on input strip number %d in TIFF file %s\n",
                      stripCount, filename);
              return(0);
            }
	  //	  if (result != 5200)
	  //	    Log("TIFF result = %d\n", (int) result);
          imageOffset += result;
        }

      // Deal with photometric interpretations
      if (TIFFGetField(image, TIFFTAG_PHOTOMETRIC, &photo) == 0)
	{
	  sprintf(error, "TIFF file %s has an undefined photometric interpretation\n",
		filename);
	  return(0);
	}
  
      if (photo != PHOTOMETRIC_MINISBLACK)
	{
	  // flip bits
	  for (count = 0; count < bufferSize; ++count)
	    buffer[count] = ~buffer[count];
	}

      TIFFClose(image);
    }
  else if (strcmp(&filename[len-4], ".pgm") == 0)
    {
      f = fopen(filename, "r");
      if (f == NULL)
	{
	  sprintf(error, "Could not open file %s for reading\n",
		  filename);
	  return(0);
	}
      if (!ReadHeader(f, &tc, &iw, &ih, &m))
	{
	  sprintf(error, "Image file %s not binary pgm.\n", filename);
	  return(0);
	}
      if ((buffer = (unsigned char *) malloc(iw * ih)) == NULL)
	{
	  sprintf(error, "Could not allocate enough memory for image in file %s\n", filename);
	  return(0);
	}
      if (fread(buffer, 1, iw * ih, f) != iw * ih)
	{
	  sprintf(error, "Image file %s apparently truncated.\n",
		  filename);
	  return(0);
	}
      fclose(f);
    }
  else
    {
      sprintf(error, "Unrecognized file extension for image file %s\n", filename);
      return(0);
    }

  if (minX <= 0 && (maxX < 0 || maxX == iw-1) &&
      minY <= 0 && (maxY < 0 || maxY == ih-1))
    {
      *width = iw;
      *height = ih;
      *pixels = buffer;
      return(1);
    }

  if (minX < 0)
    xmin = 0;
  else
    xmin = minX;
  if (maxX < 0)
    xmax = iw-1;
  else
    xmax = maxX;
  if (minY < 0)
    ymin = 0;
  else
    ymin = minY;
  if (maxY < 0)
    ymax = ih-1;
  else
    ymax = maxY;
  w = xmax - xmin + 1;
  h = ymax - ymin + 1;
  *width = w;
  *height = h;

  //  Log("resizing from %d %d to %d %d (%d) %d %d %d %d\n", iw, ih, w, h, w * h, xmin, xmax,
  //      ymin, ymax);
  if ((ptr = (unsigned char *) malloc(w * h)) == NULL)
    {
      free(buffer);
      sprintf(error, "Could not allocate enough memory for requested image subregion of %s\n", filename);
      return(0);
    }
  *pixels = ptr;
  for (y = ymin; y <= ymax; ++y)
    {
      if (y >= ih || xmin >= iw)
	memset(ptr, 0, w);
      else if (xmax < iw)
	memcpy(ptr, buffer + y * iw + xmin, w);
      else if (w < iw - xmin)
	memcpy(ptr, buffer + y * iw + xmin, w);
      else
	{
	  memcpy(ptr, buffer + y * iw + xmin, iw - xmin);
	  memset(ptr+iw-xmin, 0, xmax - iw + 1);
	}
      ptr += w;
    }
  free(buffer);
  return(1);
}

int
ReadBitmap (char *filename, unsigned char **bitmap,
	    int *width, int *height,
	    int minX, int maxX, int minY, int maxY,
	    char *error)
{
  int len;
  gzFile gzf;
  FILE *f;
  uint32 iw, ih;
  char tc;
  int ibpl;
  unsigned char *buffer;
  unsigned char *ptr;
  int rbpl;
  int m;
  uint32 xmin, xmax, ymin, ymax;
  int w, h;
  int x, y;

  len = strlen(filename);
  if (len < 5)
    {
      sprintf(error, "Image filename is too short: %s\n", filename);
      return(0);
    }

  if (strcasecmp(&filename[len-4], ".pbm") == 0)
    {
      f = fopen(filename, "rb");
      if (f == NULL)
	{
	  sprintf(error, "Could not open file %s for reading\n",
		  filename);
	  return(0);
	}
      if (!ReadHeader(f, &tc, &iw, &ih, &m) || tc != '4')
	{
	  sprintf(error, "Mask file not binary pbm: %s\n", filename);
	  return(0);
	}

      ibpl = (iw + 7) >> 3;
      if ((buffer = (unsigned char *) malloc(ih * ibpl)) == NULL)
	{
	  sprintf(error, "Could not allocate enough memory for image in file %s\n", filename);
	  return(0);
	}
      if (fread(buffer, 1, ih*ibpl, f) != ih * ibpl)
	{
	  fclose(f);
	  sprintf(error, "Image file %s apparently truncated.\n", filename);
	  return(0);
	}
      fclose(f);
    }
  else if (len > 7 && strcasecmp(&filename[len-7], ".pbm.gz") == 0)
    {
      gzf = gzopen(filename, "rb");
      if (!ReadGZHeader(gzf, &tc, &iw, &ih, &m) || tc != '4')
	{
	  sprintf(error, "Mask file not binary pbm: %s\n", filename);
	  return(0);
	}

      ibpl = (iw + 7) >> 3;
      if ((buffer = (unsigned char *) malloc(ih * ibpl)) == NULL)
	{
	  sprintf(error, "Could not allocate enough memory for image in file %s\n", filename);
	  return(0);
	}
      if (gzread(gzf, buffer, ih*ibpl) != ih * ibpl)
	{
	  fclose(f);
	  sprintf(error, "Image file %s apparently truncated.\n", filename);
	  return(0);
	}
      gzclose(gzf);
    }
  else
    {
      sprintf(error, "Unrecognized file extension for bitmap file %s\n", filename);
      return(0);
    }

  if (minX <= 0 && (maxX < 0 || maxX == iw-1) &&
      minY <= 0 && (maxY < 0 || maxY == ih-1))
    {
      *bitmap = buffer;
      *width = iw;
      *height = ih;
      return(1);
    }

  if (minX < 0)
    xmin = 0;
  else
    xmin = minX;
  if (maxX < 0)
    xmax = iw-1;
  else
    xmax = maxX;
  if (minY < 0)
    ymin = 0;
  else
    ymin = minY;
  if (maxY < 0)
    ymax = ih-1;
  else
    ymax = maxY;
  w = xmax - xmin + 1;
  h = ymax - ymin + 1;
  *width = w;
  *height = h;
  rbpl = (w + 7) >> 3;

  if ((ptr = (unsigned char *) malloc(h * rbpl)) == NULL)
    {
      free(buffer);
      sprintf(error, "Could not allocate enough memory for requested bitmap subregion of %s\n", filename);
      return(0);
    }
  *bitmap = ptr;
  for (y = ymin; y < ymax; ++y)
    {
      if (y >= ih)
	memset(ptr, 0, rbpl);
      else if (xmax < iw && (xmin & 7) == 0)
	memcpy(ptr, buffer + y * ibpl + (xmin >> 3), rbpl);
      else
	{
	  memset(ptr, 0, rbpl);
	  for (x = xmin; x <= xmax; ++x)
	    if (x < iw && (buffer[y*ibpl+(x >> 3)] & (0x80 >> (x & 7))) != 0)
	      ptr[(x - xmin) >> 3] |= (0x80 >> ((x - xmin) & 7));
	}
      ptr += rbpl;
    }
  return(1);
}

int
ReadHeader (FILE *f, char *tc, unsigned int *w, unsigned int *h, int *m)
{
  int c;
  int v;

  c = fgetc(f);
  while (c == '#')
    {
      while ((c = fgetc(f)) != EOF && c != '\n') ;
    }
  if (c != 'P')
    return(0);
  c = fgetc(f);
  if (c != '4' && c != '5' && c != '6')
    return(0);
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
    return(0);

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
    return(0);

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
    return(0);

  return(1);
}

int
ReadGZHeader (gzFile f, char *tc, unsigned int *w, unsigned int *h, int *m)
{
  int c;
  int v;

  c = gzgetc(f);
  while (c == '#')
    {
      while ((c = gzgetc(f)) != EOF && c != '\n') ;
    }
  if (c != 'P')
    return(0);
  c = gzgetc(f);
  if (c != '4' && c != '5' && c != '6')
    return(0);
  *tc = c;
  c = gzgetc(f);

  while (c == ' ' || c == '\t' || c == '\n' || c == '#')
    if (c == '#')
      {
	while ((c = gzgetc(f)) != EOF && c != '\n') ;
      }
    else
      c = gzgetc(f);
  if (c >= '0' && c <= '9')
    {
      v = 0;
      while (c != EOF && c >= '0' && c <= '9')
	{
	  v = 10 * v + (c - '0');
	  c = gzgetc(f);
	}
      *w = v;
    }
  else
    return(0);

  while (c == ' ' || c == '\t' || c == '\n' || c == '#')
    if (c == '#')
      {
	while ((c = gzgetc(f)) != EOF && c != '\n') ;
      }
    else
      c = gzgetc(f);
  if (c >= '0' && c <= '9')
    {
      v = 0;
      while (c != EOF && c >= '0' && c <= '9')
	{
	  v = 10 * v + (c - '0');
	  c = gzgetc(f);
	}
      *h = v;
    }
  else
    return(0);

  if (*tc == '4')
    return(1);

  while (c == ' ' || c == '\t' || c == '\n' || c == '#')
    if (c == '#')
      while ((c = gzgetc(f)) != EOF && c != '\n') ;
    else
      c = gzgetc(f);
  if (c >= '0' && c <= '9')
    {
      v = 0;
      while (c != EOF && c >= '0' && c <= '9')
	{
	  v = 10 * v + (c - '0');
	  c = gzgetc(f);
	}
      *m = v;
    }
  else
    return(0);

  return(1);
}

int ReadMap (char *filename,
	     MapElement** map,
	     int *level,
	     int *width, int *height,
	     int *xMin, int *yMin,
	     char *imageName, char *referenceName,
	     char *error)
{
  char imgName[PATH_MAX], refName[PATH_MAX];
  int newWidth, newHeight;

  FILE *f = fopen(filename, "r");
  if (f == NULL)
    {
      sprintf(error, "Cannot open file %s\n", filename);
      return(0);
    }
  if (fgetc(f) != 'M' || fgetc(f) != '1' || fgetc(f) != '\n' ||
      fscanf(f, "%d%d%d%d%d%s%s",
	     level,
	     &newWidth, &newHeight,
	     xMin, yMin,
	     imgName, refName) != 7 ||
      fgetc(f) != '\n')
    {
      sprintf(error, "Cannot read header of map file %s\n", filename);
      return(0);
    }
  if (imageName != NULL)
    strcpy(imageName, imgName);
  if (referenceName != NULL)
    strcpy(referenceName, refName);
  if (*width <= 0 || *height <= 0)
    *map = NULL;
  if (*map != NULL && (newWidth != *width || newHeight != *height))
    {
      free(*map);
      *map = NULL;
    }
  if (*map == NULL)
    *map = (MapElement *) malloc(newWidth * newHeight * sizeof(MapElement));
  if (*map == NULL)
    {
      sprintf(error, "Could not allocate space for map %s\n", filename);
      return(0);
    }
  *width = newWidth;
  *height = newHeight;
  if (fread(*map, sizeof(MapElement), newWidth * newHeight, f) != newWidth * newHeight)
    {
      sprintf(error, "Could not read map from file %s\n", filename);
      return(0);
    }
  fclose(f);
  return(1);
}

