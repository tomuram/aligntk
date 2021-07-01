#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <errno.h>
#include <math.h>

#include "imio.h"

int
FloodFill (int w, int h,
	   int *cluster, unsigned char *mask,
	   int id, int ix, int iy);

int
main (int argc, char **argv)
{
  int iw, ih;
  unsigned char *img;
  int x, y;
  int i;
  char msg[PATH_MAX + 256];
  int mbpl;
  unsigned char *mask = NULL;
  unsigned char *newMask = NULL;
  unsigned char *tmpMask = NULL;
  int lowerThreshold, upperThreshold;
  int sp;
  int *count;
  int *cluster;
  int nClusters;
  int maxCluster;
  int maxClusterSize;
  unsigned char *bitMask = NULL;
  int erode;

  if (argc != 6 ||
      sscanf(argv[1], "%d", &lowerThreshold) != 1 ||
      sscanf(argv[2], "%d", &upperThreshold) != 1 ||
      sscanf(argv[3], "%d", &erode) != 1)
    {
      fprintf(stderr, "Usage: genmask8 threshold erode input.pgm output.pbm\n");
      exit(1);
    }

  img = NULL;
  if (!ReadImage(argv[4], &img,
		 &iw, &ih,
		 -1, -1, -1, -1,
		 msg))
    {
      fprintf(stderr, "%s", msg);
      exit(1);
    }
  mask = (unsigned char *) malloc(ih * iw);
  newMask = (unsigned char *) malloc(ih * iw);
  memset(mask, 0, ih * iw);
  for (y = 0; y < ih; ++y)
    for (x = 0; x < iw; ++x)
      if (img[y*iw+x] >= lowerThreshold &&
	  img[y*iw+x] <= upperThreshold)
	mask[y * iw + x] = 1;
  free(img);

  /* first erode */
  for (i = 0; i < erode; ++i)
    {
      for (y = 0; y < ih; ++y)
	for (x = 0; x < iw; ++x)
	  newMask[y*iw+x] = mask[y*iw+x] &&
	    (x == 0 || mask[y*iw+x-1]) &&
	    (x >= iw-1 || mask[y*iw+x+1]) &&
	    (y == 0 || mask[(y-1)*iw+x]) &&
	    (y >= ih-1 || mask[(y+1)*iw+x]);
      tmpMask = mask;
      mask = newMask;
      newMask = tmpMask;
    }
      free(newMask);




  // clear all the cluster bits
  cluster = (int *) malloc(ih * iw * sizeof(int));
  for (y = 0; y < ih; ++y)
    for (x = 0; x < iw; ++x)
      cluster[y * iw + x] = -1;

  count = (int *) malloc(ih * iw * sizeof(int));
  nClusters = 0;
  for (y = 0; y < ih; ++y)
    for (x = 0; x < iw; ++x)
      {
	if (!mask[y * iw + x] || cluster[y * iw + x] >= 0)
	  continue;
	count[nClusters] = FloodFill(iw, ih, cluster, mask, nClusters, x, y);
	printf("Cluster %d has %d elements\n", nClusters, count[nClusters]);
	++nClusters;
      }

  // pick the largest cluster
  maxCluster = -1;
  maxClusterSize = -1;
  for (i = 0; i < nClusters; ++i)
    if (count[i] > maxClusterSize)
      {
	maxCluster = i;
	maxClusterSize = count[i];
      }
  free(count);

  for (y = 0; y < ih; ++y)
    for (x = 0; x < iw; ++x)
      {
	mask[y*iw+x] = cluster[y*iw+x] != maxCluster;
	cluster[y*iw+x] = -1;
      }

  for (x = 0; x < iw; ++x)
    {
      FloodFill(iw, ih, cluster, mask, 0, x, 0);
      FloodFill(iw, ih, cluster, mask, 0, x, ih-1);
    }
  for (y = 1; y < ih-1; ++y)
    {
      FloodFill(iw, ih, cluster, mask, 0, 0, y);
      FloodFill(iw, ih, cluster, mask, 0, iw-1, y);
    }
  free(mask);

  mbpl = (iw + 7) / 8;
  bitMask = (unsigned char *) malloc(mbpl * ih);
  memset(bitMask, 0, ih * mbpl);
  for (y = 0; y < ih; ++y)
    for (x = 0; x < iw; ++x)
      if (cluster[y*iw+x] < 0)
	bitMask[y * mbpl + (x >> 3)] |= 0x80 >> (x & 7);
  free(cluster);

  if (!WriteBitmap(argv[5], bitMask, iw, ih, UncompressedBitmap, msg))
    {
      fprintf(stderr, "Could not open %s for writing\n", argv[5]);
      exit(1);
    }
  free(bitMask);
  return(0);
}

int
FloodFill (int w, int h,
	   int *cluster, unsigned char *mask,
	   int id, int ix, int iy)
{
  unsigned int *stack = (unsigned int *) malloc(h * w * sizeof(unsigned int));
  int sp = 0;
  int count = 0;
  int x, y;
  cluster[iy * w + ix] = id;
  stack[sp++] = iy*w+ix;
  while (sp > 0)
    {
      ++count;
      --sp;
      y = stack[sp] / w;
      x = stack[sp] - y * w;
      
      /* check up */
      if (y > 0 && mask[(y-1)*w+x] && cluster[(y-1)*w+x] < 0)
	{
	  cluster[(y-1)*w+x] = id;
	  stack[sp++] = (y-1)*w+x;
	}
      /* check right */
      if (x < w-2 && mask[y*w+x+1] && cluster[y*w+x+1] < 0)
	{
	  cluster[y*w+x+1] = id;
	  stack[sp++] = y*w+x+1;
	}
      /* check down */
      if (y < h-2 && mask[(y+1)*w+x] && cluster[(y+1)*w+x] < 0)
	{
	  cluster[(y+1)*w+x] = id;
	  stack[sp++] = (y+1)*w+x;
	}
      /* check left */
      if (x > 0 && mask[y*w+x-1] && cluster[y*w+x-1] < 0)
	{
	  cluster[y*w+x-1] = id;
	  stack[sp++] = y*w+x-1;
	}
    }
  free(stack);
  return(count);
}
