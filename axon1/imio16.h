//
// imio16.h - functions to read and write 16-bit images
//  
#ifndef IMIO16_H
#define IMIO_H

#ifdef __cplusplus
extern "C" {
#endif

  int ReadImages16 (char *filename, int nImages, 
		    unsigned short **pixels,
		    int *width, int *height,
		    char *error);

#ifdef __cplusplus
}
#endif

#endif /* READER_H */
