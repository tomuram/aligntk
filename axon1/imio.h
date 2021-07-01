//
// imio.h - functions to read and write images, bitmap masks, and maps
//  
#ifndef IMIO_H
#define IMIO_H

#ifdef __cplusplus
extern "C" {
#endif

  enum ImageCompression { UncompressedImage = 0,
				  HDiffDeflateImage = 1};

  enum BitmapCompression { UncompressedBitmap = 0,
			   GZBitmap = 1};

  enum MapCompression { UncompressedMap = 0 };
			
  typedef struct MapElement {
    float x;
    float y;
    float c;
  } MapElement;
  
  int ReadImageSize (char *filename,
		     int *width, int *height,
		     char *error);

  int ReadImage (char *filename, unsigned char **pixels,
	       int *width, int *height,
	       int minX, int maxX, int minY, int maxY,
	       char *error);

  int WriteImage (char *filename, unsigned char *pixels,
		  int width, int height,
		  enum ImageCompression compressionMode,
		  char *error);

  int ReadBitmapSize (char *filename,
		      int *width, int *height,
		      char *error);

  int ReadBitmap (char *filename, unsigned char **bitmap,
		  int *width, int *height,
		  int minX, int maxX, int minY, int maxY,
		  char *error);

  int WriteBitmap (char *filename, unsigned char *bitmap,
		   int width, int height,
		   enum BitmapCompression compressionMethod,
		   char *error);

  int ReadMap (char *filename, MapElement **map,
	       int *level,
	       int *width, int *height,
	       int *xMin, int *yMin,
	       char *imageName, char *referenceName,
	       char *error);

  int WriteMap (char *filename, MapElement *map,
		int level,
		int width, int height,
		int xMin, int yMin,
		char *imageName, char *referenceName,
		enum MapCompression compressionMethod,
		char *error);

#ifdef __cplusplus
}
#endif

#endif /* READER_H */
