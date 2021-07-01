/*
 *  gen_pyramid.c  -  generates an image pyramid from a large image;
 *                    this is suitable as input to CATMAID
 *
 *  Copyright (c) 2014 Pittsburgh Supercomputing Center,
 *                     Carnegie Mellon University
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
 *       NIH NIGMS grant P41GM103712
 *
 *  HISTORY
 *    2014  Written by Greg Hood (ghood@psc.edu)
 */


#include <dirent.h>
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "imio.h"

int
main (int argc, char **argv)
{
  unsigned char inputName[PATH_MAX];
  unsigned char outputName[PATH_MAX];
  unsigned char subdirName[PATH_MAX];
  int outputTileWidth, outputTileHeight;
  int rowCol;
  char format[PATH_MAX];
  int quality;
  enum ImageCompression compressionMode;
  int tm[4];

  DIR *dir;
  struct dirent *de;
  int len;
  char baseName[PATH_MAX];
  char fn[PATH_MAX];
  char msg[PATH_MAX + 256];
  int error;

  int zeroBasis;
  int nHorizInputTiles, nVertInputTiles;
  int inputTileWidth, inputTileHeight;
  int iw, ih;

  unsigned char *out = NULL;
  unsigned char *img;
  int imw, imh;
  int row, col;
  int x, y;
  int dx, dy;
  int i, j;
  int tx, ty;
  int bx, by;
  int nLevels;
  int *nHorizOutputTiles, *nVertOutputTiles;
  int ow, oh;
  unsigned char ***ucTiles;
  unsigned int ***uiTiles;
  int **req;
  int minBX, maxBX, minBY, maxBY;
  int lvl;
  int srcMinX, srcMaxX, srcMinY, srcMaxY;
  int dstMinX, dstMaxX, dstMinY, dstMaxY;
  unsigned char *ucTile;
  unsigned int *uiTile;
  unsigned int *sTile;
  int shift;
  int offset;
  int ix, iy;
  int sdx, ddx;
  float rotation;
  int n90;
  int tmp[4];
  int ttx, tty;
  unsigned char *tImg;
  struct stat sb;

  error = 0;
  inputName[0] = '\0';
  outputName[0] = '\0';
  subdirName[0] = '\0';
  outputTileWidth = 0;
  outputTileHeight = 0;
  rowCol = -1;
  format[0] = '\0';
  nHorizInputTiles = 0;
  nVertInputTiles = 0;
  compressionMode = JpegQuality90;
  tm[0] = 1;
  tm[1] = 0;
  tm[2] = 0;
  tm[3] = 1;

  for (i = 1; i < argc; ++i)
    if (strcmp(argv[i], "-input") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-input error\n");
	    break;
	  }
	strcpy(inputName, argv[i]);
      }
    else if (strcmp(argv[i], "-output") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-output error\n");
	    break;
	  }
	strcpy(outputName, argv[i]);
      }
    else if (strcmp(argv[i], "-subdir") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-subdir error\n");
	    break;
	  }
	strcpy(subdirName, argv[i]);
      }
    else if (strcmp(argv[i], "-output_tile_size") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%dx%d", &outputTileWidth, &outputTileHeight) != 2)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-colrow_format") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-colrow_format error\n");
	    break;
	  }
	rowCol = 0;
	strcpy(format, argv[i]);
      }
    else if (strcmp(argv[i], "-rowcol_format") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-rowcol_format error\n");
	    break;
	  }
	rowCol = 1;
	strcpy(format, argv[i]);
      }
    else if (strcmp(argv[i], "-cols") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d", &nHorizInputTiles) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-rows") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d", &nVertInputTiles) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-quality") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d", &quality) != 1)
	  {
	    error = 1;
	    break;
	  }
	switch (quality)
	  {
	  case 70:
	    compressionMode = JpegQuality70;
	    break;
	  case 75:
	    compressionMode = JpegQuality75;
	    break;
	  case 80:
	    compressionMode = JpegQuality80;
	    break;
	  case 85:
	    compressionMode = JpegQuality85;
	    break;
	  case 90:
	    compressionMode = JpegQuality90;
	    break;
	  case 95:
	    compressionMode = JpegQuality95;
	    break;
	  default:
	    fprintf(stderr, "-quality must be one of 70, 75, 80, 85, 90, or 95\n");
	    error = 1;
	    break;
	  }
	if (error)
	  break;
      }
    else if (strcmp(argv[i], "-rotate") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%f", &rotation) != 1)
	  {
	    error = 1;
	    break;
	  }
	n90 = ((int) rotation) / 90;
	if (n90 * 90.0 != rotation)
	  {
	    fprintf(stderr, "Rotation must be an integral multiple of 90.\n");
	    error = 1;
	    break;
	  }
	n90 &= 3;
	//  for each 90 degrees
	//      multiply by matrix: (x) =  (0 -1) (x')
	//                          (y)    (1  0) (y')
	//  this will take output coordinates (x', y') and map into
        //    original coordinates (x, y)

	for (j = 0; j < n90; ++j)
	  {
	    tmp[0] = -tm[2];
	    tmp[1] = -tm[3];
	    tmp[2] = tm[0];
	    tmp[3] = tm[1];
	    memcpy(tm, tmp, 4 * sizeof(int));
	  }
      }
    else if (strcmp(argv[i], "-flip_x") == 0)
      {
	tm[0] = -tm[0];
	tm[2] = -tm[2];
      }
    else if (strcmp(argv[i], "-flip_y") == 0)
      {
	tm[1] = -tm[1];
	tm[3] = -tm[3];
      }
    else
      {
	fprintf(stderr, "Unrecognized option: %s\n", argv[i]);
	error = 1;
      }

  if (error)
    {
      if (i >= argc)
	fprintf(stderr, "Incomplete option: %s\n\n", argv[argc-1]);
      
      fprintf(stderr, "Usage: make_pyramid\n");
      fprintf(stderr, "           -colrow_format format_string\n");
      fprintf(stderr, "           -rowcol_format format_string\n");
      fprintf(stderr, "           -input directory\n");
      fprintf(stderr, "           -rows integer\n");
      fprintf(stderr, "           -cols integer\n");
      fprintf(stderr, "           -quality integer\n");
      fprintf(stderr, "           -output root_directory\n");
      fprintf(stderr, "           -output_tile_size WxH\n");
      fprintf(stderr, "           -quality integer\n");
      fprintf(stderr, "           -rotate ccw_direction_in_degrees\n");
      fprintf(stderr, "           -flip_x\n");
      fprintf(stderr, "           -flip_y\n");
      exit(1);
    }

  if (inputName[0] == '\0' || outputName[0] == '\0')
    {
      fprintf(stderr, "-input filename/directory and -output directory required.\n");
      exit(1);
    }
  if (outputTileWidth == 0 || outputTileHeight == 0)
    {
      fprintf(stderr, "Must specify a tile width and height with -output_tile_size WxH\n");
      exit(1);
    }

  if (stat(inputName, &sb) != 0)
    {
      fprintf(stderr, "Could not stat %s\n", inputName);
      exit(1);
    }
  if (S_ISREG(sb.st_mode))
    {
      nHorizInputTiles = 1;
      nVertInputTiles = 1;
      rowCol = -1;
    }
  if (format[0] == '\0')
    strcpy(format, "c%0.2dr%0.2d");

  if (nHorizInputTiles == 0 || nVertInputTiles == 0)
    {
      // find the number of input rows and columns
      dir = opendir(inputName);
      if (dir == NULL)
	{
	  fprintf(stderr, "Could not open directory %s for input tiles\n",
		  inputName);
	  exit(1);
	}
      nHorizInputTiles = 0;
      nVertInputTiles = 0;
      zeroBasis = 0;
      while ((de = readdir(dir)) != NULL)
	{
	  len = strlen(de->d_name);
	  if (len > 4 && strcasecmp(&(de->d_name[len-4]), ".tif") == 0 &&
	      sscanf(de->d_name, "c%dr%d", &col, &row) == 2)
	    {
	      if (row == 0 || col == 0)
		zeroBasis = 1;
	      if (col > nHorizInputTiles)
		nHorizInputTiles = col;
	      if (row > nVertInputTiles)
		nVertInputTiles = row;
	    }
	}
      closedir(dir);
      if (zeroBasis)
	{
	  ++nHorizInputTiles;
	  ++nVertInputTiles;
	}
    }

  if (zeroBasis)
    row = col = 0;
  else
    row = col = 1;

  if (rowCol < 0)
    strcpy(fn, inputName);
  else
    {
      if (rowCol)
	sprintf(baseName, format, row, col);
      else
	sprintf(baseName, format, col, row);
      sprintf(fn, "%s/%s", inputName, baseName);
    }
  if (!ReadImageSize(fn, &inputTileWidth, &inputTileHeight, msg))
    {
      fprintf(stderr, "Could not read image size of %s\n", fn);
      exit(1);
    }

  // if an axis-swapping transformation, exchange variables
  if (tm[0] == 0)
    {
      tmp[0] = nHorizInputTiles;
      nHorizInputTiles = nVertInputTiles;
      nVertInputTiles = tmp[0];
      tmp[0] = inputTileWidth;
      inputTileWidth = inputTileHeight;
      inputTileHeight = tmp[0];
    }

  // compute how many levels we need
  nLevels = 0;
  iw = nHorizInputTiles * inputTileWidth;
  ih = nVertInputTiles * inputTileHeight;
  while (iw > (1 << nLevels) * outputTileWidth ||
	 ih > (1 << nLevels) * outputTileHeight)
    ++nLevels;
  // nLevels is now a sufficient level to encompass the entire image;
  //   increment this to yield a count of the number of valid levels
  ++nLevels;

  printf("Output image pyramid will have %d levels.\n", nLevels);
  if (mkdir(outputName, 0777) != 0 &&
      errno != EEXIST)
    {
      fprintf(stderr, "Could not create output directory %s\n", outputName);
      exit(1);
    }

  nHorizOutputTiles = (int *) malloc(nLevels * sizeof(int));
  nVertOutputTiles = (int *) malloc(nLevels * sizeof(int));
  req = (int **) malloc(nLevels * sizeof(int *));
  ucTiles = (unsigned char ***) malloc(nLevels * sizeof(unsigned char **));
  uiTiles = (unsigned int ***) malloc(nLevels * sizeof(unsigned int **));
  for (lvl = 0; lvl < nLevels; ++lvl)
    {
      sprintf(fn, "%s/%d", outputName, lvl);
      if (mkdir(fn, 0777) != 0 &&
	  errno != EEXIST)
	{
	  fprintf(stderr, "Could not create directory %s\n", fn);
	  exit(1);
	}
      if (subdirName[0] != '\0')
	{
	  sprintf(fn, "%s/%d/%s", outputName, lvl, subdirName);
	  if (mkdir(fn, 0777) != 0 &&
	      errno != EEXIST)
	    {
	      fprintf(stderr, "Could not create directory %s\n", fn);
	      exit(1);
	    }
	}
      ow = (1 << lvl) * outputTileWidth;
      oh = (1 << lvl) * outputTileHeight;
      nHorizOutputTiles[lvl] = (iw + ow - 1) / ow;
      nVertOutputTiles[lvl] = (ih + oh - 1) / oh;
      ucTiles[lvl] = (unsigned char **) malloc(nVertOutputTiles[lvl] *
					     nHorizOutputTiles[lvl] *
					     sizeof(unsigned char *));
      uiTiles[lvl] = (unsigned int **) malloc(nVertOutputTiles[lvl] *
					      nHorizOutputTiles[lvl] *
					      sizeof(unsigned int *));
      req[lvl] = (int *) malloc(nVertOutputTiles[lvl] *
				nHorizOutputTiles[lvl] *
				sizeof(int));
      memset(ucTiles[lvl], 0, 
	     nVertOutputTiles[lvl] * nHorizOutputTiles[lvl] *
	     sizeof(unsigned char *));
      memset(uiTiles[lvl], 0, 
	     nVertOutputTiles[lvl] * nHorizOutputTiles[lvl] *
	     sizeof(unsigned int *));
      memset(req[lvl], 0,
	     nVertOutputTiles[lvl] * nHorizOutputTiles[lvl] *
	     sizeof(int));

      for (by = 0; by < nVertOutputTiles[lvl]; ++by)
	{
	  sprintf(fn, "%s/%d/%s/%d", outputName, lvl, subdirName, by);
	  if (mkdir(fn, 0777) != 0 &&
	      errno != EEXIST)
	    {
	      fprintf(stderr, "Could not create directory %s\n", fn);
	      exit(1);
	    }
	}
    }


  // first go through and see how many items (subtiles or input tiles) are required
  //    to output each tile
  for (tx = 0; tx < nHorizInputTiles; ++tx)
    for (ty = 0; ty < nVertInputTiles; ++ty)
      {
	minBX = (tx * inputTileWidth) / outputTileWidth;
	maxBX = (tx * inputTileWidth + inputTileWidth - 1) / outputTileWidth;
	minBY = (ty * inputTileHeight) / outputTileHeight;
	maxBY = (ty * inputTileHeight + inputTileHeight - 1) / outputTileHeight;
	for (by = minBY; by <= maxBY; ++by)
	  for (bx = minBX; bx <= maxBX; ++bx)
	    ++req[0][by * nHorizOutputTiles[0] + bx];
      }
  for (lvl = 0; lvl < nLevels-1; ++lvl)
    for (by = 0; by < nVertOutputTiles[lvl]; ++by)
      for (bx = 0; bx < nHorizOutputTiles[lvl]; ++bx)
	++req[lvl+1][(by >> 1) * nHorizOutputTiles[lvl+1] + (bx >> 1)];

  for (tx = 0; tx < nHorizInputTiles; ++tx)
    for (ty = 0; ty < nVertInputTiles; ++ty)
      {
	// ttx, tty are the transformed tx and ty
	if (tm[0] == 1)
	  ttx = tx;
	else if (tm[0] == -1)
	  ttx = nHorizInputTiles - tx - 1;
	else if (tm[1] == 1)
	  ttx = ty;
	else
	  ttx = nVertInputTiles - ty - 1;
	if (tm[3] == 1)
	  tty = ty;
	else if (tm[3] == -1)
	  tty = nVertInputTiles - ty - 1;
	else if (tm[2] == 1)
	  tty = tx;
	else
	  tty = nHorizInputTiles - tx - 1;

	if (rowCol < 0)
	  strcpy(fn, inputName);
	else
	  {
	    if (rowCol)
	      sprintf(baseName, format,
		      tty + (zeroBasis ? 0 : 1),
		      ttx + (zeroBasis ? 0 : 1));
	    else
	      sprintf(baseName, format,
		      ttx + (zeroBasis ? 0 : 1),
		      tty + (zeroBasis ? 0 : 1));
	    sprintf(fn, "%s/%s", inputName, baseName);
	  }
	if (!ReadImage(fn, &img, &imw, &imh, -1, -1, -1, -1, msg))
	    {
	      fprintf(stderr, "Could not read image %s:\n  error: %s\n",
		      fn, msg);
	      exit(1);
	    }
	if (imw * imh != inputTileWidth * inputTileHeight)
	  {
	    fprintf(stderr, "Dimensions of image %s are inconsistent.\n", fn);
	    exit(1);
	  }

	// check if a reordering within the tile is needed
	if (tm[0] != 1 || tm[3] != 1)
	  {
	    tImg = (unsigned char *) malloc(imw * imh * sizeof(unsigned char));
	    if (tm[0] < 0)
	      dx = inputTileWidth - 1;
	    else if (tm[1] < 0)
	      dx = inputTileHeight - 1;
	    else
	      dx = 0;
	    if (tm[2] < 0)
	      dy = inputTileWidth - 1;
	    else if (tm[3] < 0)
	      dy = inputTileHeight - 1;
	    else
	      dy = 0;
	    for (y = 0; y < inputTileHeight; ++y)
	      for (x = 0; x < inputTileWidth; ++x)
		tImg[y*inputTileWidth + x] = 
		  img[(tm[2] * x + tm[3] * y + dy) * imw + tm[0] * x + tm[1] * y + dx];
	    free(img);
	    img = tImg;
	    if (tm[0] == 0)
	      {
		tmp[0] = imw;
		imw = imh;
		imh = tmp[0];
	      }
	  }

	if (imw != inputTileWidth || imh != inputTileHeight)
	  {
	    fprintf(stderr, "Dimensions of image %s are inconsistent.\n", fn);
	    exit(1);
	  }
	
	minBX = (tx * inputTileWidth) / outputTileWidth;
	maxBX = (tx * inputTileWidth + inputTileWidth - 1) / outputTileWidth;
	minBY = (ty * inputTileHeight) / outputTileHeight;
	maxBY = (ty * inputTileHeight + inputTileHeight - 1) / outputTileHeight;
	for (by = minBY; by <= maxBY; ++by)
	  for (bx = minBX; bx <= maxBX; ++bx)
	    {
	      srcMinX = bx * outputTileWidth - tx * inputTileWidth;
	      srcMaxX = srcMinX + outputTileWidth - 1;
	      dstMinX = 0;
	      dstMaxX = outputTileWidth-1;
	      if (srcMinX < 0)
		{
		  dstMinX -= srcMinX;
		  srcMinX = 0;
		}
	      if (srcMaxX >= inputTileWidth)
		{
		  dstMaxX -= srcMaxX - inputTileWidth + 1;
		  srcMaxX = inputTileWidth-1;
		}
	      srcMinY = by * outputTileHeight - ty * inputTileHeight;
	      srcMaxY = srcMinY + outputTileHeight - 1;
	      dstMinY = 0;
	      dstMaxY = outputTileHeight-1;
	      if (srcMinY < 0)
		{
		  dstMinY -= srcMinY;
		  srcMinY = 0;
		}
	      if (srcMaxY >= inputTileHeight)
		{
		  dstMaxY -= srcMaxY - inputTileHeight + 1;
		  srcMaxY = inputTileHeight-1;
		}

	      ucTile = ucTiles[0][by * nHorizOutputTiles[0] + bx];
	      if (ucTile == NULL)
		{
		  ucTile = (unsigned char *) malloc(outputTileHeight *
						    outputTileWidth *
						    sizeof(unsigned char));
		  memset(ucTile, 0,
			 outputTileHeight * outputTileWidth * sizeof(unsigned char));
		  ucTiles[0][by * nHorizOutputTiles[0] + bx] = ucTile;
		}

	      for (y = dstMinY; y <= dstMaxY; ++y)
		memcpy(&ucTile[y * outputTileWidth + dstMinX],
		       &img[(srcMinY + y - dstMinY) * inputTileWidth + srcMinX],
		       dstMaxX - dstMinX + 1);
	      lvl = 0;
	      ix = bx;
	      iy = by;
	      while (lvl < nLevels &&
		     --req[lvl][iy * nHorizOutputTiles[lvl] + ix] == 0)
		{
		  if (lvl > 0)
		    {
		      uiTile = uiTiles[lvl][iy * nHorizOutputTiles[lvl] + ix];
		      if (uiTile == NULL)
			{
			  fprintf(stderr, "Internal error: uiTile is NULL\n");
			  exit(1);
			}
		      ucTile = (unsigned char *) malloc(outputTileHeight *
							outputTileWidth *
							sizeof(unsigned char));
		      shift = 2 * lvl;
		      offset = 1 << (shift - 1);
		      for (y = 0; y < outputTileHeight; ++y)
			{
			  ddx = y * outputTileWidth;
			  for (x = 0; x < outputTileWidth; ++x)
			    ucTile[x + ddx] = (uiTile[x + ddx] + offset) >> shift;
			}
		    }

		  sprintf(fn, "%s/%d/%s/%d/%d.jpg",
			  outputName, lvl, subdirName, iy, ix);
		  printf("Writing image %s\n", fn);
		  if (!WriteImage(fn, ucTile,
				  outputTileWidth, outputTileHeight,
				  compressionMode, msg))
		    {
		      fprintf(stderr, "Could not write tile %s\n%s\n", fn, msg);
		      exit(1);
		    }

		  if (lvl+1 < nLevels)
		    {
		      // propagate upwards
		      sTile = uiTiles[lvl+1][(iy >> 1) * nHorizOutputTiles[lvl+1] +
					     (ix >> 1)];
		      if (sTile == NULL)
			{
			  sTile = (unsigned int *) malloc(outputTileHeight *
							  outputTileWidth *
							  sizeof(unsigned int));
			  memset(sTile, 0,
				 outputTileHeight *
				 outputTileWidth *
				 sizeof(unsigned int));
			  uiTiles[lvl+1][(iy >> 1) * nHorizOutputTiles[lvl+1] +
					 (ix >> 1)] = sTile;
			}

		      for (y = 0; y < outputTileHeight; ++y)
			{
			  sdx = y * outputTileWidth;
			  ddx = (((iy & 1) * outputTileHeight + y) >> 1) *
			    outputTileWidth +
			    (((ix & 1) * outputTileWidth) >> 1);
			  if (lvl == 0)
			    for (x = 0; x < outputTileWidth; ++x)
			      sTile[(x >> 1) + ddx] += ucTile[sdx + x];
			  else
			    for (x = 0; x < outputTileWidth; ++x)
			      sTile[(x >> 1) + ddx] += uiTile[sdx + x];
			}
		    }

		  free(ucTile);
		  if (lvl > 0)
		    {
		      free(uiTile);
		      uiTiles[lvl][iy * nHorizOutputTiles[lvl] + ix] = NULL;
		    }
		  else
		    ucTiles[0][iy * nHorizOutputTiles[0] + ix] = NULL;
		  ucTile = NULL;
		  uiTile = NULL;

		  ++lvl;
		  ix >>= 1;
		  iy >>= 1;
		}
	      if (lvl < nLevels && req[lvl][iy * nHorizOutputTiles[lvl] + ix] < 0)
		{
		  fprintf(stderr, "Internal error: req[%d][%d*%d + %d] = %d\n",
			  lvl, iy, nHorizOutputTiles[lvl], ix,
			  req[lvl][iy * nHorizOutputTiles[lvl] + ix]);
		  exit(1);
		}
	    }
	free(img);
	img = 0;
      }
  printf("Pyramid generation completed.\n");
  return(0);
}
