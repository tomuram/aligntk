//
// reader.h - functions to read images and bitmap masks
//  
#ifndef READER_H
#define READER_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct MapElement {
  float x;
  float y;
  float c;
} MapElement;


int ReadImage (char *filename, unsigned char **pixels,
	       int *width, int *height,
	       int minX, int maxX, int minY, int maxY,
	       char *error);

int ReadBitmap (char *filename, unsigned char **bitmap,
		int *width, int *height,
		int minX, int maxX, int minY, int maxY,
		char *error);

int ReadMap (char *filename, MapElement **map,
	     int *level,
	     int *width, int *height,
	     int *xMin, int *yMin,
	     char *imageName, char *referenceName,
	     char *error);

#ifdef __cplusplus
}
#endif

#endif /* READER_H */
