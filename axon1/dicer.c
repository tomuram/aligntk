#include <stdio.h>
#include <tiffio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <dirent.h>

int ReadHeader (FILE *f, char *tc, int *w, int *h, int *m);
void Error (char *fmt, ...);

int
main (int argc, char **argv)
{
  TIFF *image;
  uint16 bpp;
  uint32 width, height;
  unsigned char *tile;
  unsigned char *buffer;
  unsigned long bufferSize, count;
  FILE *f;
  int inputAndOutputFilesGiven;
  int i;
  char outputFileName[PATH_MAX];
  int len;
  unsigned long long nBytes;
  int w, h;
  char tc;
  int m;
  uint32 rowsperstrip = (uint32) -1;
  int row, col;
  int nRows, nCols;
  int scanLine;
  int size;
  unsigned char *img;
  int digits;
  int x, y;

  size = 2048;
  if (argc < 3 || argc > 4 ||
      strlen(argv[1]) <= 4 ||
      strcasecmp(&(argv[1][strlen(argv[1])-4]), ".pgm") != 0 ||
      argc == 4 && sscanf(argv[3], "%d", &size) != 1)
    Error("Usage: dicer input.pgm output_prefix [size]\n");
  
  // open the pgm image
  if ((f = fopen(argv[1], "r")) == NULL)
    Error("Could not open incoming image: %s\n", argv[i]);
      
  if (!ReadHeader(f, &tc, &w, &h, &m))
    Error("Could not read header of file %s\n", argv[i]);

  if (m == 255)
    bpp = 8;
  else if (m == 65535)
    bpp = 16;
  else
    Error("Invalid max value field in pgm file %s\n", argv[1]);
  nBytes = ((unsigned long long) w) * h * (bpp / 8);
  img = (unsigned char *) malloc(nBytes);
  if (fread(img, nBytes, 1, f) != 1)
    Error("Could not read from file %s\n", argv[i]);
  tile = (unsigned char *) malloc(size * size);

  nRows = (h + size - 1) / size;
  nCols = (w + size - 1) / size;
  if (nRows >= nCols)
    digits = (int) ceil(log((double) nRows) / log(10.0));
  else
    digits = (int) ceil(log((double) nCols) / log(10.0));
  if (digits < 1)
    digits = 1;
  printf("nRows = %d nCols = %d digits = %d\n",
	 nRows, nCols, digits);
  
  for (col = 0; col < nCols; ++col)
    for (row = 0; row < nRows; ++row)
      {
	memset(tile, 0, size*size);
	for (y = 0; y < size; ++y)
	  {
	    if (row * size + y >= h)
	      continue;
	    for (x = 0; x < size; ++x)
	      {
		if (col * size + x >= w)
		  continue;
		tile[y * size + x] = img[(row * size + y) *
					 ((unsigned long long) w) +
					 (col * size + x)];
	      }
	  }

	sprintf(outputFileName, "%s%sc%0.*dr%0.*d.tif",
		argv[2], (argv[2][strlen(argv[2])-1] == '/' ? "" : "."),
		digits, col, digits, row);

	// open the TIFF output file
	if ((image = TIFFOpen(outputFileName, "w")) == NULL)
	  Error("Could not open output file %s\n", outputFileName);

	TIFFSetField(image, TIFFTAG_IMAGEWIDTH, (uint32) size);
	TIFFSetField(image, TIFFTAG_IMAGELENGTH, (uint32) size);
	TIFFSetField(image, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
	TIFFSetField(image, TIFFTAG_SAMPLESPERPIXEL, 1);
	TIFFSetField(image, TIFFTAG_BITSPERSAMPLE, bpp);
	TIFFSetField(image, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	TIFFSetField(image, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	TIFFSetField(image, TIFFTAG_ROWSPERSTRIP,
		     TIFFDefaultStripSize(image, rowsperstrip));

	if (TIFFScanlineSize(image) > size * (bpp / 8))
	  buffer = (unsigned char *)_TIFFmalloc(TIFFScanlineSize(image));
	else
	  buffer = (unsigned char *)_TIFFmalloc(size * (bpp / 8));

	for (scanLine = 0; scanLine < size; ++scanLine)
	  {
	    memcpy(buffer, &tile[scanLine * size], size);
	    if (TIFFWriteScanline(image, buffer, scanLine, 0)  <  0)
	      Error("Could not write to tif file %s\n", outputFileName);
	  }

	(void) TIFFClose(image);
	_TIFFfree(buffer);
      }
}


int
ReadHeader (FILE *f, char *tc, int *w, int *h, int *m)
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

void Error (char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
  fflush(stderr);
  exit(1);
}
