#include <stdio.h>
#include <tiffio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

int ReadHeader (FILE *f, char *tc, int *w, int *h, int *m);
void Error (char *fmt, ...);

int
main (int argc, char **argv)
{
  TIFF *image;
  uint16 bpp;
  uint32 width, height;
  unsigned char *buffer;
  unsigned long bufferSize, count;
  FILE *f;
  int inputAndOutputFilesGiven;
  int i;
  char *outputFileName;
  int len;
  unsigned long long nBytes;
  int w, h;
  char tc;
  int m;
  uint32 rowsperstrip = (uint32) -1;
  int row;

  inputAndOutputFilesGiven = 0;
  if (argc == 3 && strlen(argv[1]) > 4 && strlen(argv[2]) > 4 &&
      strcasecmp(&(argv[1][strlen(argv[1])-4]), ".pgm") == 0 &&
      strcmp(&(argv[2][strlen(argv[2])-4]), ".tif") == 0)
    inputAndOutputFilesGiven = 1;
  else
    for (i = 1; i < argc; ++i)
      if (strlen(argv[i]) <= 4 ||
	  strcasecmp(&(argv[i][strlen(argv[i])-4]), ".pgm") != 0)
	Error("Usage: pgm2tif input.pgm output.tif\n\tor: pgm2tif input1.pgm input2.pgm ... inputN.pgm\n");
  
  for (i = 1; inputAndOutputFilesGiven ? i < 2 : i < argc; ++i)
    {
      // open the pgm image
      if ((f = fopen(argv[i], "r")) == NULL)
	Error("Could not open incoming image: %s\n", argv[i]);
      
      if (!ReadHeader(f, &tc, &w, &h, &m))
	Error("Could not read header of file %s\n", argv[i]);

      nBytes = ((unsigned long long) w) * h;
      if (nBytes >= (4ULL * 1024) * (1024 * 1024))
	Error("pgm file %s is too large to convert to tif (4GB limit)\n",
	      argv[i]);

      if (m == 255)
	bpp = 8;
      else if (m == 65535)
	bpp = 16;
      else
	Error("Invalid max value field in pgm file %s\n", argv[1]);

      if (inputAndOutputFilesGiven)
	outputFileName = argv[2];
      else
	{
	  outputFileName = (char *) malloc(strlen(argv[i]) + 1);
	  strcpy(outputFileName, argv[i]);
	  len = strlen(outputFileName);
	  if (len > 4 && outputFileName[len-4] == '.')
	    strcpy(&outputFileName[len-3], "tif");
	  else
	    Error("Cannot construct output file name from input file %s\n",
		  argv[i]);
	}

      // open the TIFF output file
      if ((image = TIFFOpen(outputFileName, "w")) == NULL)
	Error("Could not open output file %s\n", outputFileName);

      TIFFSetField(image, TIFFTAG_IMAGEWIDTH, (uint32) w);
      TIFFSetField(image, TIFFTAG_IMAGELENGTH, (uint32) h);
      TIFFSetField(image, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
      TIFFSetField(image, TIFFTAG_SAMPLESPERPIXEL, 1);
      TIFFSetField(image, TIFFTAG_BITSPERSAMPLE, bpp);
      TIFFSetField(image, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
      TIFFSetField(image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
      TIFFSetField(image, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
      TIFFSetField(image, TIFFTAG_ROWSPERSTRIP,
		   TIFFDefaultStripSize(image, rowsperstrip));

      if (TIFFScanlineSize(image) > w * (bpp / 8))
	buffer = (unsigned char *)_TIFFmalloc(TIFFScanlineSize(image));
      else
	buffer = (unsigned char *)_TIFFmalloc(w * (bpp / 8));

      for (row = 0; row < h; row++)
	{
	  if (fread(buffer, w * (bpp/8), 1, f) != 1)
	    Error("Could not read from file %s\n", argv[i]);

	  if (TIFFWriteScanline(image, buffer, row, 0)  <  0)
	    Error("Could not write to tif file %s\n", outputFileName);
	}

      (void) TIFFClose(image);
      _TIFFfree(buffer);

      if (!inputAndOutputFilesGiven)
	free(outputFileName);
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
