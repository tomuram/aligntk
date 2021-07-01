//
// dt.h - function to compute distance tranformation of various types
//  
#ifndef EDT_H
#define EDT_H

#ifdef __cplusplus
extern "C" {
#endif

#define EUCLIDEAN_DISTANCE		0
#define EUCLIDEAN_DISTANCE_SQUARED	1
#define MANHATTAN_DISTANCE		2
#define CHESSBOARD_DISTANCE		3

void computeDistance (int type, int nx, int ny, unsigned char *mask, float *dist);

#ifdef __cplusplus
}
#endif

#endif /* EDT_H */
