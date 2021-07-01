#include <stdio.h>
#include <tiffio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

void Error (char *fmt, ...);

int
main (int argc, char **argv)
{
  TIFF *image;
  uint16 photo, bps, spp, fillorder;
  uint32 width, height;
  tsize_t stripSize;
  unsigned long imageOffset, result;
  int stripMax, stripCount;
  unsigned char *buffer;
  unsigned long bufferSize, count;
  FILE *f;
  int inputAndOutputFilesGiven;
  int i;
  char *outputFileName;
  int len;

  inputAndOutputFilesGiven = 0;
  if (argc == 3 && strlen(argv[1]) > 4 && strlen(argv[2]) > 4 &&
      (strcasecmp(&(argv[1][strlen(argv[1])-4]), ".tif") == 0 ||
       strcasecmp(&(argv[1][strlen(argv[1])-5]), ".tiff") == 0) &&
      strcmp(&(argv[2][strlen(argv[2])-4]), ".pgm") == 0)
    inputAndOutputFilesGiven = 1;
  else
    for (i = 1; i < argc; ++i)
      if (strlen(argv[i]) <= 4 ||
	  strcasecmp(&(argv[i][strlen(argv[i])-4]), ".tif") != 0 &&
	  strcasecmp(&(argv[i][strlen(argv[i])-5]), ".tiff") != 0)
	Error("Usage: tif2pgm input.tif output.pgm\n\tor: tif2pgm input1.tif input2.tif ... inputN.tif\n");
  
  for (i = 1; inputAndOutputFilesGiven ? i < 2 : i < argc; ++i)
    {
      // Open the TIFF image
      if ((image = TIFFOpen(argv[i], "r")) == NULL)
	Error("Could not open incoming image: %s\n", argv[i]);

      // Check that it is of a type that we support
      if (TIFFGetField(image, TIFFTAG_BITSPERSAMPLE, &bps) == 0 || bps != 8)
	Error("Either undefined or unsupported number of bits per sample (bps = %d)\n", bps);
      
      if (TIFFGetField(image, TIFFTAG_SAMPLESPERPIXEL, &spp) == 0 || spp != 1)
	Error("Either undefined or unsupported number of samples per pixel (spp = %d)\n", spp);

      if (TIFFGetField(image, TIFFTAG_IMAGEWIDTH, &width) == 0)
	Error("Image does not define its width\n");
  
      if (TIFFGetField(image, TIFFTAG_IMAGELENGTH, &height) == 0)
	Error("Image does not define its height (length)\n");
  
      // Read in the possibly multiple strips
      stripSize = TIFFStripSize(image);
      stripMax = TIFFNumberOfStrips(image);
      imageOffset = 0;
  
      bufferSize = TIFFNumberOfStrips(image) * stripSize;
      if ((buffer = (unsigned char *) malloc(bufferSize)) == NULL)
	Error("Could not allocate enough memory for the uncompressed image (bytes = %d)\n", bufferSize);
      
      for (stripCount = 0; stripCount < stripMax; ++stripCount)
	{
	  if ((result = TIFFReadEncodedStrip(image, stripCount,
					     buffer + imageOffset,
					     stripSize)) == -1)
	    Error("Read error on input strip number %d\n", stripCount);
	  imageOffset += result;
	}

      // Deal with photometric interpretations
      if (TIFFGetField(image, TIFFTAG_PHOTOMETRIC, &photo) == 0)
	Error("Image has an undefined photometric interpretation\n");
  
      if (photo != PHOTOMETRIC_MINISBLACK)
	{
	  // flip bits
	  printf("Fixing the photometric interpretation\n");
	  for(count = 0; count < bufferSize; ++count)
	    buffer[count] = ~buffer[count];
	}

      TIFFClose(image);
     
      if (inputAndOutputFilesGiven)
	outputFileName = argv[2];
      else
	{
	  outputFileName = (char *) malloc(strlen(argv[i]) + 1);
	  strcpy(outputFileName, argv[i]);
	  len = strlen(outputFileName);
	  if (len > 4 && outputFileName[len-4] == '.')
	    strcpy(&outputFileName[len-3], "pgm");
	  else if (len > 5 && outputFileName[len-5] == '.')
	    strcpy(&outputFileName[len-4], "pgm");
	  else
	    Error("Cannot construct output file name from input file %s\n",
		  argv[i]);
	}
      f = fopen(outputFileName, "w");
      if (f == NULL)
	Error("Could not open output file %s\n", outputFileName);
      fprintf(f, "P5\n%d %d\n255\n", width, height);
      if (fwrite(buffer, 1, width*height, f) != width*height)
	Error("Could not write pixels to output file %s\n", outputFileName);
      if (!inputAndOutputFilesGiven)
	free(outputFileName);
    }
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
