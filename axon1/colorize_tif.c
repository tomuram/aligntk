#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <tiffio.h>
#include <stdarg.h>
#include <limits.h>
#include <math.h>

int
ReadTiffImage (char *filename, unsigned char **buffer,
	       int *width, int *height,
	       char *error);

int
main (int argc, char **argv)
{
  int w, h;
  int x, y;
  char msg[1024];
  unsigned char *pixels;
  TIFF *image;
  uint32 rowsperstrip = (uint32) -1;
  float red, green, blue;
  unsigned char *img;
  char filename[1024];

  if (argc < 6 ||
      sscanf(argv[3], "%f", &red) != 1 ||
      sscanf(argv[4], "%f", &green) != 1 ||
      sscanf(argv[5], "%f", &blue) != 1)
    {
      fprintf(stderr, "Usage: colorize_tif input_file output_file red green blue\n");
      fprintf(stderr, "   where red, green, blue are floats in the range of 0.0 to 1.0\n");
      exit(1);
    }

  printf("red = %f  green = %f  blue = %f\n", red, green, blue);

  if (!ReadTiffImage(argv[1], &img, &w, &h, msg))
    {
      fprintf(stderr, "Could not read image %s\n%s\n", argv[1], msg);
      exit(1);
    }
  printf("Read image of size %d x %d\n", w, h);

  pixels = (unsigned char *) malloc(w * h * 3);
  for (y = 0; y < h; ++y)
    for (x = 0; x < w; ++x) 
      {
	pixels[(y*w+x)*3] = (int) floor(red * img[y*w+x] + 0.5);
	pixels[(y*w+x)*3+1] = (int) floor(green * img[y*w+x] + 0.5);
	pixels[(y*w+x)*3+2] = (int) floor(blue * img[y*w+x] + 0.5);
      }

  printf("Finished constructing image.\n");

  // open the TIFF output file
  strcpy(filename, argv[2]);

  printf("Going to open %s for writing.\n", filename);

  if ((image = TIFFOpen(filename, "w")) == NULL)
    {
      fprintf(stderr, "Could not open file %s for writing\n", argv[2]);
      return(0);
    }

  TIFFSetField(image, TIFFTAG_IMAGEWIDTH, (uint32) w);
  TIFFSetField(image, TIFFTAG_IMAGELENGTH, (uint32) h);
  TIFFSetField(image, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
  TIFFSetField(image, TIFFTAG_SAMPLESPERPIXEL, 3);
  TIFFSetField(image, TIFFTAG_BITSPERSAMPLE, 8);
  TIFFSetField(image, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
  TIFFSetField(image, TIFFTAG_COMPRESSION, COMPRESSION_ADOBE_DEFLATE);
  TIFFSetField(image, TIFFTAG_PREDICTOR, PREDICTOR_HORIZONTAL);
  TIFFSetField(image, TIFFTAG_ROWSPERSTRIP,
	       TIFFDefaultStripSize(image, rowsperstrip));
  for (y = 0; y < h; ++y)
    if (TIFFWriteScanline(image, &pixels[y*w*3], y, 0)  <  0)
      {
	fprintf(stderr, "Could not write to tif file %s\n", filename);
	return(0);
      }
  TIFFClose(image);

  printf("Finished writing image.\n");
  free(img);
  free(pixels);
  return(0);
}  



int
ReadTiffImage (char *filename, unsigned char **buffer,
	       int *width, int *height,
	       char *error)
{
  TIFF *image;
  uint32 iw, ih;
  uint16 photo, bps, spp, fillorder;
  tsize_t stripSize;
  unsigned long imageOffset, result;
  int stripMax, stripCount;
  unsigned long bufferSize, count;
  unsigned char *b;

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
      sprintf(error, "Image %s does not define its height (length)\n",
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

  if ((b = (unsigned char *) malloc(bufferSize)) == NULL)
    {
      sprintf(error, "Could not allocate enough memory for the uncompressed image (bytes = %lu) of TIFF file %s\n", bufferSize, filename);
      return(0);
    }

  for (stripCount = 0; stripCount < stripMax; ++stripCount)
    {
      if ((result = TIFFReadEncodedStrip(image, stripCount,
					 b + imageOffset,
					 stripSize)) == -1)
	{
	  sprintf(error, "Read error on input strip number %d in TIFF file %s\n",
		  stripCount, filename);
	  return(0);
	}
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
	b[count] = ~b[count];
    }
  
  TIFFClose(image);

  *buffer = b;
  *width = iw;
  *height = ih;
  return(1);
}

