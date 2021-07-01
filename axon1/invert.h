#ifndef INVERT_H
#define INVERT_H

#ifdef __cplusplus
extern "C" {
#endif

#include "imio.h"

  typedef struct InverseMapElement
  {
    int minX, maxX;
    int minY, maxY;
  } InverseMapElement;
  
  typedef struct InverseMap {
    MapElement *map;   /* the original map */
    int nx, ny;	       /* the number of elements in the original map */
    int nxp, nyp;      /* the number of elements in the inverse map table */
    InverseMapElement *inverseMap;  /* the bounding boxes for the inverse map */
    float xMin, yMin;  /* the upper left corner of the entire inverse map */
    float scale;       /* the length of each side of an InverseMapEntry square */
    int lastX, lastY;  /* the last used position in the original map */
    unsigned char *marked; /* the MapElements that have been marked */
    long long *tried;  /* the MapElements that have been tried */
  } InverseMap;

  InverseMap* InvertMap (MapElement *map, int nx, int ny);
  int Invert (InverseMap *inverseMap, float *x, float *y, float xp, float yp);
  void FreeInverseMap (InverseMap *inverseMap);

#ifdef __cplusplus
}
#endif

#endif /* INVERT_H */
