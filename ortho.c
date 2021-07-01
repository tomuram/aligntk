/*
 * ortho.c  - generate a set of images orthogonal to
 *            a given image stack
 *
 *  Copyright (c) 2010-2013 National Resource for Biomedical
 *                          Supercomputing,
 *                          Pittsburgh Supercomputing Center,
 *                          Carnegie Mellon University
 *
 *  This file is part of AlignTK.
 *
 *  AlignTK is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  AlignTK is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with AlignTK.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Acknowledgements:
 *     Development of this code was supported in part by
 *       NIH NCRR grant 5P41RR006009 and
 *       NIH NIGMS grant P41GM103712
 *
 *
 *  HISTORY
 *    2010  Written by Greg Hood (ghood@psc.edu)
 */

#include <stdio.h>
#include <mpi.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <sched.h>
#include <limits.h>
#include <sys/resource.h>

#include "imio.h"

#define LINE_LENGTH	256

int p;      /* rank of this process */
int np;     /* number of processes in this run */

char imageListName[PATH_MAX];
char inputName[PATH_MAX];
char outputName[PATH_MAX];
char format[256];
int minX = -1;
int maxX = -1;
int minY = -1;
int maxY = -1;
int minZ = -1;
int maxZ = -1;
int xzImages = 1;
int memoryLimit = 1024;
int every = 1;

FILE *logFile = 0;

void Error (char *fmt, ...);
void Log (char *fmt, ...);

int
main (int argc, char **argv)
{
  char fn[PATH_MAX];
  FILE *f;
  int x, y, z;
  int i, j, k;
  int error;
  int im;
  int n;
  int w, h;
  MPI_Status status;
  struct stat sb;
  char hostName[256];
  unsigned int hv;
  int imageNamesSize;
  int imageNamesPos;
  char *imageNames;
  char line[LINE_LENGTH];
  char imageName[PATH_MAX];
  int imageNameLen;
  int width, height;
  char msg[PATH_MAX+256];
  cpu_set_t cpumask;
  int digits;
  char outputNameFormat[64];
  size_t outputImageSize;
  size_t sizePerImage;
  size_t blockSize;
  size_t nBlocks;
  int nPasses;
  int maxBlocks;
  int minImageNumber, maxImageNumber;
  int *sendCounts, *sendOffsets;
  int *receiveCounts, *receiveOffsets;
  int *imageCounts, *imageOffsets;
  int imageNumber;
  int nImagesPerProcess;
  int pass;
  int block;
  int source, target;
  unsigned char *img;
  unsigned char *buffer;
  int iv;
  int minV, maxV, minW, maxW;
  unsigned char *dst;
  int nVirtualImages;
  int nImages;
  char **images;
  
  /* DECLS */

  /* initialize MPI */
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    Error("Could not do MPI_Init\n");
  if (MPI_Comm_size(MPI_COMM_WORLD, &np) != MPI_SUCCESS)
    Error("Could not do MPI_Comm_size\n");
  if (MPI_Comm_rank(MPI_COMM_WORLD, &p) != MPI_SUCCESS)
    Error("Could not do MPI_Comm_rank\n");

  printf("I am node %d of %d\n", p, np);
  fflush(stdout);
  
#if 0
  // TEMPORARY FOR DEBUGGING
  struct rlimit rlim;
  rlim.rlim_cur = 10000000000;
  rlim.rlim_max = 10000000000;
  setrlimit(RLIMIT_CORE, &rlim);
#endif

  sprintf(fn, "logs/%.2d.log", p);
  logFile = fopen(fn, "w");
  gethostname(hostName, 255);
  Log("Process %d is running on %s\n", p, hostName);
  if (p == 0)
    {
      error = 0;
      imageListName[0] = '\0';
      inputName[0] = '\0';
      outputName[0] = '\0';
      format[0] = '\0';
      
      for (i = 0; i < argc; ++i)
	Log("ARGV[%d] = %s\n", i, argv[i]);
      for (i = 1; i < argc; ++i)
	if (strcmp(argv[i], "-input") == 0)
	  {
	    if (++i == argc)
	      {
		error = 1;
		break;
	      }
	    strcpy(inputName, argv[i]);
	  }
	else if (strcmp(argv[i], "-image_list") == 0)
	  {
	    if (++i == argc)
	      {
		error = 1;
		break;
	      }
	    strcpy(imageListName, argv[i]);
	  }
	else if (strcmp(argv[i], "-output") == 0)
	  {
	    if (++i == argc)
	      {
		error = 1;
		break;
	      }
	    strcpy(outputName, argv[i]);
	  }	    
	else if (strcmp(argv[i], "-x") == 0)
	  {
	    if (++i == argc ||
		sscanf(argv[i], "%d-%d", &minX, &maxX) < 1)
	      {
		error = 1;
		break;
	      }
	  }
	else if (strcmp(argv[i], "-y") == 0)
	  {
	    if (++i == argc ||
		sscanf(argv[i], "%d-%d", &minY, &maxY) < 1)
	      {
		error = 1;
		break;
	      }
	  }
	else if (strcmp(argv[i], "-z") == 0)
	  {
	    if (++i == argc ||
		sscanf(argv[i], "%d-%d", &minZ, &maxZ) < 1)
	      {
		error = 1;
		break;
	      }
	  }
	else if (strcmp(argv[i], "-format") == 0)
	  {
	    if (++i == argc)
	      {
		error = 1;
		break;
	      }
	    strcpy(format, argv[i]);
	  }
	else if (strcmp(argv[i], "-memory") == 0)
	  {
	    if (++i == argc ||
		sscanf(argv[i], "%d", &memoryLimit) != 1)
	      {
		error = 1;
		break;
	      }
	  }
	else if (strcmp(argv[i], "-yz_images") == 0)
	  xzImages = 0;
	else if (strcmp(argv[i], "-xz_images") == 0)
	  xzImages = 1;
	else if (strcmp(argv[i], "-every") == 0)
	  {
	    if (++i == argc ||
		sscanf(argv[i], "%d", &every) != 1)
	      {
		error = 1;
		break;
	      }
	  }
	else error = 1;

      if (error)
	{
	  fprintf(stderr, "Usage: ortho\n");
	  fprintf(stderr, "              -input input_prefix\n");
	  fprintf(stderr, "              -output output_prefix\n");
	  fprintf(stderr, "              -image_list images_file\n");
	  fprintf(stderr, "              [-x minX-maxX]\n");
	  fprintf(stderr, "              [-y minY-maxY]\n");
	  fprintf(stderr, "              [-z minZ-maxZ]\n");
	  fprintf(stderr, "              [-format zFormat]\n");
	  fprintf(stderr, "              [-memory memory_limit_in_MB]\n");
	  fprintf(stderr, "              [-every multiple]\n");
	  exit(1);
	}
      
      /* images_file contains one line for each image:
            image_name
	    image_name
            ...
	
      */

      /* check that at least minimal parameters were supplied */
      if (imageListName[0] == '\0' ||
	  inputName[0] == '\0' ||
	  outputName[0] == '\0')
	Error("-image_list, -input, and -output must be specified\n");
    }
  
  if (p == 0)
    {
      /* read the images file */
      f = fopen(imageListName, "r");
      if (f == NULL)
	Error("Could not open file %s for reading\n", imageListName);
      nImages = 0;
      imageNamesSize = 0;
      imageNamesPos = 0;
      imageNames = 0;
      minImageNumber = 1000000000;
      maxImageNumber = -1000000000;
      width = -1;
      height = -1;
      while (fgets(line, LINE_LENGTH, f) != NULL)
	{
	  if (line[0] == '\0' || line[0] == '#')
	    continue;
	  if (sscanf(line, "%s", imageName) != 1)
	    Error("Malformed line in %s:\n%s\n", imageListName, line);
	  if (format[0] != '\0')
	    {
	      if (sscanf(imageName, format, &imageNumber) != 1)
		Error("Could not parse image name %s according to format %s\n",
		      imageName, format);
	      if (minZ >= 0 && maxZ >= 0 &&
		  (imageNumber < minZ || imageNumber > maxZ))
		continue;
	      if (imageNumber < minImageNumber)
		minImageNumber = imageNumber;
	      if (imageNumber > maxImageNumber)
		maxImageNumber = imageNumber;
	    }

	  imageNameLen = strlen(imageName);
	  while (imageNamesPos + imageNameLen + 1 >= imageNamesSize)
	    {
	      imageNamesSize = (imageNamesSize > 0) ? imageNamesSize * 2 : 1024;
	      imageNames = (char *) realloc(imageNames, imageNamesSize);
	    }
	  strcpy(&imageNames[imageNamesPos], imageName);
	  imageNamesPos += imageNameLen + 1;

	  if (nImages == 0)
	    {
	      sprintf(fn, "%s%s.tif", inputName, imageName);
	      if (!ReadImage(fn, &img, &w, &h,
			     -1, -1, -1, -1,
			     msg))
		Error("Could not read image %s:\n%s\n",
		      fn, msg);
	      width = w;
	      height = h;
	    }

	  ++nImages;
	}
      fclose(f);

      imageNamesSize = imageNamesPos;
      imageNames = (char *) realloc(imageNames, imageNamesSize);

      if (minX < 0)
	{
	  minX = 0;
	  maxX = width-1;
	}
      else if (maxX < 0)
	maxX = minX;
      if (minX >= width)
	Error("minimum x value is greater than image width\n");
      if (maxX >= width)
	maxX = width - 1;
      
      if (minY < 0)
	{
	  minY = 0;
	  maxY = height-1;
	}
      else if (maxY < 0)
	maxY = minY;
      if (minY >= height)
	Error("minimum y value is greater than image height\n");
      if (maxY >= height)
	maxY = height - 1;
	      
      if (minZ < 0)
	{
	  if (format[0] == '\0')
	    {
	      minZ = 0;
	      maxZ = nImages - 1;
	    }
	  else
	    {
	      minZ = minImageNumber;
	      maxZ = maxImageNumber;
	    }
	}
      if (maxZ < 0)
	{
	  if (format[0] == '\0')
	    maxZ = nImages - 1;
	  else
	    maxZ = maxImageNumber;
	}
    }

  /* broadcast the image info */
  if (MPI_Bcast(inputName, PATH_MAX, MPI_CHAR, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(outputName, PATH_MAX, MPI_CHAR, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(format, 256, MPI_CHAR, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&minX, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&maxX, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&minY, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&maxY, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&minZ, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&maxZ, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&xzImages, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&memoryLimit, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&every, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&width, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&height, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&nImages, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Bcast(&imageNamesSize, 1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("Broadcast of nImages or imageNamesSize failed.\n");
  if (p != 0)
    imageNames = (char *) malloc(imageNamesSize * sizeof(char));
  if (MPI_Bcast(imageNames, imageNamesSize, MPI_CHAR, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    Error("Broadcast of imageNames failed.\n");

  /* compute number of passes */
  nVirtualImages = maxZ - minZ + 1;
  nImagesPerProcess = (nVirtualImages + np - 1) / np;
  if (xzImages)
    {
      minV = minY;
      maxV = maxY;
      minW = minX;
      maxW = maxX;
    }
  else
    {
      minV = minX;
      maxV = maxX;
      minW = minY;
      maxW = maxY;
    }
  if (maxV > 0)
    digits = ((int) floor(log10((double) maxV))) + 1;
  else
    digits = 1;
  sprintf(outputNameFormat, "%%s%%0.%dd.tif", digits);
  outputImageSize = ((size_t) (maxZ - minZ + 1)) * (maxW - minW + 1);
  sizePerImage = ((size_t) (maxW - minW + 1)) * nImagesPerProcess;
  blockSize = np * sizePerImage;
  maxBlocks = (memoryLimit * 1024LL * 1024LL -
	       ((long long) width) * height - outputImageSize) / blockSize;
  //  nBlocks = ((maxV - minV + 1) + np - 1) / np;
  nBlocks = ((maxV - minV + every) / every + np - 1) / np;
  if (nBlocks > maxBlocks)
    nBlocks = maxBlocks;
  //  nPasses = (maxV - minV + 1 + nBlocks * np - 1) / (nBlocks * np);
  nPasses = ((maxV - minV + every) / every + nBlocks * np - 1) / (nBlocks * np);

  if (p == 0)
    Log("Going to make %d passes through the image stack\n",
	nPasses);
  buffer = (unsigned char *) malloc(nBlocks * blockSize);
  if (buffer == NULL)
    Error("Could not allocate %zd bytes of buffer space.\n",
	  nBlocks * blockSize);

  images = (char **) malloc(nVirtualImages * sizeof(char *));
  memset(images, 0, nVirtualImages * sizeof(char *));
  imageNamesPos = 0;
  for (i = 0; i < nImages; ++i)
    {
      if (format[0] == '\0')
	images[i] = &imageNames[imageNamesPos];
      else
	{
	  if (sscanf(&imageNames[imageNamesPos], format, &imageNumber) != 1)
	    Error("Could not parse image name %s according to format %s\n",
		  imageName, format);
	  images[imageNumber - minZ] = &imageNames[imageNamesPos];
	}
      imageNamesPos += strlen(&imageNames[imageNamesPos]) + 1;
    }

  sendCounts = (int *) malloc(np * sizeof(int));
  sendOffsets = (int *) malloc(np * sizeof(int));
  receiveCounts = (int *) malloc(np * sizeof(int));
  receiveOffsets = (int *) malloc(np * sizeof(int));

  for (pass = 0; pass < nPasses; ++pass)
    {
      for (i = 0; i < nImagesPerProcess; ++i)
	{
	  iv = p * nImagesPerProcess + i;
	  if (iv >= nVirtualImages)
	    break;
	  if (images[iv] != NULL)
	    {
	      sprintf(fn, "%s%s.tif", inputName, images[iv]);
	      if (!ReadImage(fn, &img, &w, &h,
			     -1, -1, -1, -1,
			     msg))
		Error("Could not read image %s:\n%s\n",
		      fn, msg);
	      if (w != width || h != height)
		Error("Size of image %s (%dx%d) is inconsistent with first image (%dx%d).\n",
		      fn, w, h, width, height);
	    }
	  else
	    {
	      img = (unsigned char *) malloc(width * height);
	      memset(img, 0, width * height);
	    }

	  /* save the part destined for output */
	  for (block = 0; block < nBlocks; ++block)
	    for (target = 0; target < np; ++target)
	      if (xzImages)
		{
		  //		  y = minY + (pass * nBlocks + block) * np + target;
		  y = minY + ((pass * nBlocks + block) * np + target) * every;
		  if (y <= maxY)
		    memcpy(&buffer[(((size_t) (block * np + target)) *
				    nImagesPerProcess + i) *
				   (maxX - minX + 1)],
			   &img[y*w + minX],
			   maxX - minX + 1);
		  else
		    memset(&buffer[(((size_t) (block * np + target)) *
				    nImagesPerProcess + i) *
				   (maxX - minX + 1)],
			   0,
			   maxX - minX + 1);
		}
	      else
		{
		  //		  x = minX + (pass * nBlocks + block) * np + target;
		  x = minX + ((pass * nBlocks + block) * np + target) * every;
		  if (x <= maxX)
		    {
		      dst = &buffer[(((size_t) (block * np + target)) *
				     nImagesPerProcess + i) *
				    (maxY - minY + 1)];
		      for (y = minY; y <= maxY; ++y)
			*dst++ = img[y*w+x];
		    }
		  else
		    memset(&buffer[(((size_t) (block * np + target)) *
				    nImagesPerProcess + i) *
				   (maxY - minY + 1)],
			   0,
			   maxY - minY + 1);
		}
	  free(img);
	}

      img = (unsigned char *) malloc(np * nImagesPerProcess * (maxW - minW + 1));
      for (block = 0; block < nBlocks; ++block)
	{
	  iv  = minV + (pass * nBlocks + block) * np * every;
	  if (iv > maxV)
	    break;

	  memset(img, 0, outputImageSize);

	  for (target = 0; target < np; ++target)
	    {
	      //	      iv = minV + (pass * nBlocks + block) * np + target;
	      iv = minV + ((pass * nBlocks + block) * np + target) * every;
	      if (iv <= maxV)
		{
		  sendCounts[target] = nImagesPerProcess * (maxW - minW + 1);
		  sendOffsets[target] = target * nImagesPerProcess * (maxW - minW + 1);
		}
	      else
		{
		  sendCounts[target] = 0;
		  sendOffsets[target] = 0;
		}
	    }

	  //	  iv = minV + (pass * nBlocks + block) * np + p;
	  iv = minV + ((pass * nBlocks + block) * np + p) * every;
	  for (source = 0; source < np; ++source)
	    if (iv <= maxV)
	      {
		receiveCounts[source] = nImagesPerProcess * (maxW - minW + 1);
		receiveOffsets[source] = source * nImagesPerProcess * (maxW - minW + 1);
	      }
	    else
	      {
		receiveCounts[source] = 0;
		receiveOffsets[source] = 0;
	      }

	  if (MPI_Alltoallv(&buffer[((size_t) block) * np * nImagesPerProcess *
				    (maxW - minW + 1)],
			    sendCounts, sendOffsets, MPI_UNSIGNED_CHAR,
			    img, receiveCounts, receiveOffsets, MPI_UNSIGNED_CHAR,
			    MPI_COMM_WORLD) != MPI_SUCCESS)
	    Error("MPI_Alltoallv() failed.\n");

	  if (iv <= maxV)
	    {
	      sprintf(fn, outputNameFormat,
		      outputName,
		      iv);
	      if (!WriteImage(fn, img,
			      maxW - minW + 1, maxZ - minZ + 1,
			      UncompressedImage, msg))
		Error("Could not write output image %s:\n%s\n", fn, msg);
	    }
	}
      free(img);
    }

  Log("FINALIZING\n");
  MPI_Finalize();
  fclose(logFile);
  return(0);
}

void Error (char *fmt, ...)
{
  va_list args;

  if (logFile != NULL)
    {
      va_start(args, fmt);
      fprintf(logFile, "%f: ERROR: ", MPI_Wtime());
      vfprintf(logFile, fmt, args);
      va_end(args);
      fflush(logFile);
    }
  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
  fflush(stderr);
  abort();
}

void Log (char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);
  fprintf(logFile, "%f: ", MPI_Wtime());
  vfprintf(logFile, fmt, args);
  va_end(args);
  fflush(logFile);
}
