/*
 *  clean_maps.cc  -  allows a user to clean problematic regions from maps
 *
 *  Copyright (c) 2008-2013 National Resource for Biomedical
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
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <dirent.h>
#include <stdarg.h>
#include <FL/gl.h>
#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Slider.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/Fl_Value_Slider.H>
#include <FL/fl_draw.H>
#include <GL/glu.h>
#include <string>
#include <vector>
#include <algorithm>
#include "imio.h"
#include "invert.h"

using std::string;
using std::max;
using std::min;
using std::vector;

#define MAX_LINE_LENGTH		255

struct Pair
{
  char *imageName;
  int imageMinX, imageMaxX, imageMinY, imageMaxY; 
  char *refName;
  int refMinX, refMaxX, refMinY, refMaxY;
  char *pairName;
  bool modified;
};

struct MapElementAndDistortion
{
  int x, y;
  float distortion;
};

bool lesserDistortion (const MapElementAndDistortion &m0,
		       const MapElementAndDistortion &m1)
{ return(m0.distortion < m1.distortion); }

class MyWindow;

class ImageWindow: public Fl_Gl_Window
{
public:
  ImageWindow (MyWindow *p, int w, int h);
  void setTexture (unsigned char *image, int w, int h);
  void computeWarped ();

private:
  int handle (int event);
  void draw ();
  void loadImage ();

public:
  bool textureValid;  // true if the texture is up-to-date
  char name[PATH_MAX];// file name 

protected:
  GLuint texture;
  bool textureInitialized;
  bool circleInitialized;

private:
  void initTexture (int w, int h);
  void initCircle ();
  MyWindow *parent;         // the parent window
  int textureWidth;       // the width in pixels of the texture map we are usin\g
  int textureHeight;      // the height in pixels of the texture map we are usi\ng
  float maxx;             // the portion of the texture width we are actually
                          //   using (0 < maxx <= 1.0)
  float maxy;             // the portion of the texture height we are actually
                          //   using (0 < maxx <= 1.0)
  int width, height;

  unsigned char *warped;
  float scale;        // how many pixels on the screen correspond to
                      //    one pixel in the image
  float xOffset, yOffset;
  GLuint circle;
};

#define VALID_BIT	0x01
#define REJECT_BIT	0x02	// red
#define ACCEPT_BIT	0x04	// blue
#define THRESHOLD_BIT	0x08	// orange or green
#define CLUSTER_BIT	0x10	// green

class MyWindow: public Fl_Window
{
 public:
  MyWindow (int x, int y, int w, int h);
  void loadSection ();
  void draw ();
  void saveSection ();
  void create ();
  void reshape (int w, int h);
  void selectSeed ();
  void updateMask ();

public:  
  int width, height;

  int index;              // the index within the pair list of the current pair
  int lastIndex;          // the last "good" index
  int increment;          // the last increment to index (-1 or 1)
  bool sectionsLoaded;    // true if we have already loaded the currentPair

  unsigned char *image;
  int imageWidth, imageHeight;
  int imageMinX, imageMaxX;
  int imageMinY, imageMaxY;

  int mapWidth, mapHeight;
  int mapOffsetX, mapOffsetY;
  int mapLevel;
  int mapFactor;
  MapElement *map;
  InverseMap *inverseMap;

  unsigned char *mapMask;
  int *cluster;
  vector<int> seeds;       // cluster seeds

  float avgOrthogonal;
  float avgDiagonal;

  Fl_Value_Slider* orthogonalSlider;
  Fl_Value_Slider* diagonalSlider;
  ImageWindow *imageWindow;
};

char extension[8];
char pairsFile[PATH_MAX];
char imagesFile[PATH_MAX];
char inputName[PATH_MAX];
char mapsName[PATH_MAX];
char corrName[PATH_MAX];
char outputName[PATH_MAX];
int reductionFactor = 1;
MyWindow* window;
Pair* pairs;
int nPairs;
bool readOnly = false;
bool partial = true;

unsigned char colors[16][3] =
  {{255, 0, 0},     // red
   {255, 255, 0},   // yellow
   {0, 255, 0},     // green
   {0, 255, 255},   // turquoise
   {255, 153, 0},   // orange

   {0, 51, 255},    // bright blue
   {204, 0, 255},   // bright purple
   {255, 204, 102}, // gold
   {0, 204, 255},   // light blue
   {255, 102, 153}, // pink

   {153, 0, 0},     // dark red
   {204, 255, 0},   // light green
   {255, 204, 204}, // flesh tone
   {204, 0, 153},   // lavender
   {204, 153, 0},   // amber

   {204, 255, 204}  // pale green
  };

int sliderHeight = 25;

void OrthogonalSliderCallback (Fl_Widget *w, void *data);
void DiagonalSliderCallback (Fl_Widget *w, void *data);

int
main (int argc, char **argv)
{
  FILE *f;
  int w, h;
  int m;
  int v;
  int i;
  int n;
  char tc;
  int error;
  int pos;
  char fn[PATH_MAX];
  int len;
  int z;
  bool maximizeWindow = false;
  int windowFactor = 1;
  int windowWidth, windowHeight;
  char imgn[PATH_MAX], refn[PATH_MAX], pairn[PATH_MAX];
  int imgMinX, imgMaxX, imgMinY, imgMaxY;
  int refMinX, refMaxX, refMinY, refMaxY;

  error = 0;
  strcpy(extension, "tif");
  pairsFile[0] = '\0';
  imagesFile[0] = '\0';
  inputName[0] = '\0';
  mapsName[0] = '\0';
  corrName[0] = '\0';
  outputName[0] = '\0';
  for (i = 1; i < argc; ++i)
    if (strcmp(argv[i], "-big") == 0)
      maximizeWindow = true;
    else if (strcmp(argv[i], "-input") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-input error\n");
	    break;
	  }
	strcpy(inputName, argv[i]);
      }
    else if (strcmp(argv[i], "-maps") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-maps error\n");
	    break;
	  }
	strcpy(mapsName, argv[i]);
      }
    else if (strcmp(argv[i], "-corr") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-corr error\n");
	    break;
	  }
	strcpy(corrName, argv[i]);
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
    else if (strcmp(argv[i], "-pairs") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-pairs error\n");
	    break;
	  }
	strcpy(pairsFile, argv[i]);
      }
    else if (strcmp(argv[i], "-images") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-images error\n");
	    break;
	  }
	strcpy(imagesFile, argv[i]);
      }
    else if (strcmp(argv[i], "-reduction") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%d", &reductionFactor) != 1)
	  {
	    error = 1;
	    fprintf(stderr, "-reduction error\n");
	    break;
	  }
      }
    else if (strcmp(argv[i], "-readonly") == 0)
      readOnly = true;
    else if (strcmp(argv[i], "-pgm") == 0)
      strcpy(extension, "pgm");
    else
      error = 1;

  if (argc == 1 || error)
    {
      fprintf(stderr, "Usage: clean_maps -input images_file_prefix\n");
      fprintf(stderr, "                  -maps maps_file_prefix\n");
      fprintf(stderr, "                  -output maps_file_prefix\n");
      fprintf(stderr, "                  [-pairs pairfile]\n");
      fprintf(stderr, "                  [-images imagesfile]\n");
      fprintf(stderr, "                  [-reduction image_reduction_factor]\n");
      fprintf(stderr, "                  [-readonly]\n");
      exit(1);
    }

  /* check that at least minimal parameters were supplied */
  if (inputName[0] == '\0' || mapsName[0] == '\0' ||
      pairsFile[0] == '\0' && imagesFile[0] == '\0')
    {
      fprintf(stderr, "-input, -maps, and -pairs (or -images) parameters must be specified.\n");
      exit(1);
    }
  if (!readOnly && (corrName[0] == '\0' || outputName[0] == '\0'))
    {
      fprintf(stderr, "-corr and -output parameters must be specified if not readonly.\n");
      exit(1);
    }

  pairs = 0;
  if (pairsFile[0] != '\0')
    {
      char line[MAX_LINE_LENGTH+1];

      f = fopen(pairsFile, "r");
      if (f == NULL)
	{
	  fprintf(stderr, "Could not open pairs file %s\n", pairsFile);
	  exit(1);
	}
      
      
      while (fgets(line, MAX_LINE_LENGTH, f) != NULL)
	{
	  if (line[0] == '\0' || line[0] == '#')
	    continue;
	  if (sscanf(line, "%s %d %d %d %d %s %d %d %d %d %s",
		     imgn, &imgMinX, &imgMaxX, &imgMinY, &imgMaxY,
		     refn, &refMinX, &refMaxX, &refMinY, &refMaxY,
		     pairn) != 11)
	    {
	      if (sscanf(line, "%s %s %s", imgn, refn, pairn) != 3)
		{
		  fprintf(stderr, "Invalid line in pairs file %s:\n%s\n", pairsFile, line);
		  exit(1);
		}
	      imgMinX = -1;
	      imgMaxX = -1;
	      imgMinY = -1;
	      imgMaxY = -1;
	      refMinX = -1;
	      refMaxX = -1;
	      refMinY = -1;
	      refMaxY = -1;
	    }

	  if ((nPairs & 1023) == 0)
	    pairs = (Pair *) realloc(pairs, (nPairs + 1024) * sizeof(Pair));
	  pairs[nPairs].imageName = (char *) malloc(strlen(imgn)+1);
	  strcpy(pairs[nPairs].imageName, imgn);
	  pairs[nPairs].imageMinX = imgMinX;
	  pairs[nPairs].imageMaxX = imgMaxX;
	  pairs[nPairs].imageMinY = imgMinY;
	  pairs[nPairs].imageMaxY = imgMaxY;
	  pairs[nPairs].refName = (char *) malloc(strlen(refn)+1);
	  strcpy(pairs[nPairs].refName, refn);
	  pairs[nPairs].refMinX = refMinX;
	  pairs[nPairs].refMaxX = refMaxX;
	  pairs[nPairs].refMinY = refMinY;
	  pairs[nPairs].refMaxY = refMaxY;
	  pairs[nPairs].pairName = (char *) malloc(strlen(pairn)+1);
	  strcpy(pairs[nPairs].pairName, pairn);
	  ++nPairs;
	}
      fclose(f);
      printf("%d pairs listed in pairs file.\n", nPairs);
    }
  else
    {
      f = fopen(imagesFile, "r");
      if (f == NULL)
	{
	  fprintf(stderr, "Could not open images file %s\n", imagesFile);
	  exit(1);
	}
      char line[MAX_LINE_LENGTH+1];
      while (fgets(line, MAX_LINE_LENGTH, f) != NULL)
	{
	  n = sscanf(line, "%s%d%d%d%d",
		     imgn, &imgMinX, &imgMaxX, &imgMinY, &imgMaxY);
	  if (n != 1 && n != 5)
	    {
	      fprintf(stderr, "Invalid line in images file %s\n", imagesFile);
	      exit(1);
	    }
	  if (n == 1)
	    {
	      imgMinX = -1;
	      imgMaxX = -1;
	      imgMinY = -1;
	      imgMaxY = -1;
	    }
	  if ((nPairs & 1023) == 0)
	    pairs = (Pair *) realloc(pairs, (nPairs + 1024) * sizeof(Pair));
	  pairs[nPairs].imageName = (char *) malloc(strlen(imgn)+1);
	  strcpy(pairs[nPairs].imageName, imgn);
	  pairs[nPairs].imageMinX = imgMinX;
	  pairs[nPairs].imageMaxX = imgMaxX;
	  pairs[nPairs].imageMinY = imgMinY;
	  pairs[nPairs].imageMaxY = imgMaxY;
	  pairs[nPairs].refName = (char *) malloc(1);
	  pairs[nPairs].refName[0] = '\0';
	  pairs[nPairs].refMinX = -1;
	  pairs[nPairs].refMaxX = -1;
	  pairs[nPairs].refMinY = -1;
	  pairs[nPairs].refMaxY = -1;
	  pairs[nPairs].pairName = (char *) malloc(1);
	  pairs[nPairs].pairName[0] = '\0';
	  ++nPairs;
	}
      fclose(f);
      printf("%d pairs listed in pairs file.\n", nPairs);
    }

  /*  Fl::lock(); */
  Fl::visual(FL_RGB);
  Fl::gl_visual(FL_RGB);
  window = new MyWindow(0, 0, 800, 600);
  window->create();
  window->loadSection();
  window->show();
  printf("Windows created\n");
  return(Fl::run());
}

MyWindow::MyWindow (int x, int y, int w, int h) : Fl_Window (x, y, w, h)
{
  index = 0;
  lastIndex = -1;
  increment = 1;
  sectionsLoaded = false;

  mapWidth = 0;
  mapHeight = 0;

  image = 0;

  mapFactor = 0;
  map = 0;
  inverseMap = 0;
  mapMask = 0;
  cluster = 0;

  width = -1;
  height = -1;

  orthogonalSlider = 0;
  diagonalSlider = 0;
  imageWindow = 0;
}

void
MyWindow::create ()
{
  seeds.push_back(456);
  orthogonalSlider = new Fl_Value_Slider(0, h() - 2*sliderHeight,
					 w(), sliderHeight,
					 "Orthogonal Threshold");
  orthogonalSlider->bounds(0.0, 1.0);
  orthogonalSlider->type(FL_HOR_NICE_SLIDER);
  orthogonalSlider->textsize(12);
  orthogonalSlider->selection_color(FL_GREEN);
  orthogonalSlider->precision(3);
  orthogonalSlider->callback(OrthogonalSliderCallback);
  orthogonalSlider->value(0.1);

  diagonalSlider = new Fl_Value_Slider(0, h()-sliderHeight,
				       w(), sliderHeight,
				       "Diagonal Threshold");
  diagonalSlider->bounds(0.0, 1.0);
  diagonalSlider->type(FL_HOR_NICE_SLIDER);
  diagonalSlider->textsize(12);
  diagonalSlider->selection_color(FL_GREEN);
  diagonalSlider->precision(3);
  diagonalSlider->callback(DiagonalSliderCallback);
  diagonalSlider->value(0.1);

  imageWindow = new ImageWindow(this, w(), h() - 2 * sliderHeight);
  size_range(100,100);
  seeds.push_back(789);
}

void
MyWindow::loadSection ()
{
  int n;
  int i;
  int x, y;
  float *warp;
  FILE *f;
  int refSection;
  char fn[PATH_MAX];
  char errorMsg[PATH_MAX + 256];
  char imgn[PATH_MAX], refn[PATH_MAX];
  double sumOrthogonal, sumDiagonal;
  double distX, distY;
  double distA, distB;

 retry:
  if (index < 0)
    {
      fprintf(stderr, "Could not find a valid image pair to display.\n");
      exit(1);
    }
  printf("Loading section %s\n", pairs[index].imageName);

  printf("image = %p\n", image);
  if (image != 0)
    {
      free(image);
      image = 0;
    }

  /* read in the image */
  if (pairs[index].imageMinX >= 0)
    imageMinX = pairs[index].imageMinX / reductionFactor;
  else
    imageMinX = -1;
  if (pairs[index].imageMaxX >= 0)
    imageMaxX = pairs[index].imageMaxX / reductionFactor;
  else
    imageMaxX = -1;
  if (pairs[index].imageMinY >= 0)
    imageMinY = pairs[index].imageMinY / reductionFactor;
  else
    imageMinY = -1;
  if (pairs[index].imageMaxY >= 0)
    imageMaxY = pairs[index].imageMaxY / reductionFactor;
  else
    imageMaxY = -1;

  sprintf(fn, "%s%s.%s", inputName, pairs[index].imageName, extension);
  if (!ReadImage(fn, &image, &imageWidth, &imageHeight,
		 imageMinX, imageMaxX,
		 imageMinY, imageMaxY,
		 errorMsg))
    {
      if (partial)
	{
	  index += increment;
	  if (index < 0 || index >= nPairs)
	    index = lastIndex;
	  goto retry;
	}
      fprintf(stderr, "Could not read in image %s -- error: %s\n", fn, errorMsg);
      exit(1);
    }
  if (imageMinX < 0)
    imageMinX = 0;
  if (imageMaxX < 0)
    imageMaxX = imageWidth - 1;
  if (imageMinY < 0)
    imageMinY = 0;
  if (imageMaxY < 0)
    imageMaxY = imageHeight - 1;

  /* read in the map */
  sprintf(fn, "%s%s.map", mapsName,
	  (pairs[index].pairName[0] != '\0' ?
	   pairs[index].pairName : pairs[index].imageName));
  if (!ReadMap(fn, &map, &mapLevel, &mapWidth, &mapHeight,
	       &mapOffsetX, &mapOffsetY, imgn, refn, errorMsg))
    {
      if (partial)
	{
	  index += increment;
	  if (index < 0 || index >= nPairs)
	    index = lastIndex;
	  goto retry;
	}	  
      fprintf(stderr, "Error reading map %s:\n  %s\n", fn, errorMsg);
      exit(1);
    }

  /* construct the inverse map */
  mapFactor = (1 << mapLevel);
  inverseMap = InvertMap(map, mapWidth, mapHeight);

  /* compute the average ratios */
  sumOrthogonal = 0.0;
  sumDiagonal = 0.0;
  n = 0;
  for (y = 0; y < mapHeight-1; ++y)
    for (x = 0; x < mapWidth-1; ++x)
      {
	if (map[y * mapWidth + x].c == 0.0 || map[y * mapWidth + x + 1].c == 0.0 ||
	    map[(y + 1) * mapWidth + x].c == 0.0 || map[(y + 1) * mapWidth + x + 1].c == 0.0)
	  continue;
	++n;
	distX = 0.5 * hypot(map[y * mapWidth + x].x - map[y * mapWidth + x + 1].x,
		      map[y * mapWidth + x].y - map[y * mapWidth + x + 1].y) +
	  0.5 * hypot(map[(y + 1) * mapWidth + x].x - map[(y + 1) * mapWidth + x + 1].x,
		      map[(y + 1) * mapWidth + x].y - map[(y + 1) * mapWidth + x + 1].y);
	distY = 0.5 * hypot(map[y * mapWidth + x].x - map[(y + 1) * mapWidth + x].x,
			    map[y * mapWidth + x].y - map[(y + 1) * mapWidth + x].y) +
	  0.5 * hypot(map[y * mapWidth + x + 1].x - map[(y + 1) * mapWidth + x + 1].x,
		      map[y * mapWidth + x + 1].y - map[(y + 1) * mapWidth + x + 1].y);
	sumOrthogonal += distY / distX;
	distA = hypot(map[y * mapWidth + x].x - map[(y + 1) * mapWidth + x + 1].x,
		      map[y * mapWidth + x].y - map[(y + 1) * mapWidth + x + 1].y);
	distB = hypot(map[y * mapWidth + x + 1].x - map[(y + 1) * mapWidth + x].x,
		      map[y * mapWidth + x + 1].y - map[(y + 1) * mapWidth + x].y);
	sumDiagonal += distB/distA;
      }
  avgOrthogonal = sumOrthogonal / n;
  avgDiagonal = sumDiagonal / n;
  printf("avg orthogonal ratio = %f\n", avgOrthogonal);
  printf("avg diagonal ratio = %f\n", avgDiagonal);

  /* set up the map mask */
  printf("mapMask = %p\n", image);
  if (mapMask != 0)
    free(mapMask);
  mapMask = (unsigned char *) malloc(mapHeight * mapWidth * sizeof(unsigned char));
  memset(mapMask, 0, mapHeight * mapWidth * sizeof(unsigned char));
  printf("cluster = %p\n", image);
  if (cluster != 0)
    free(cluster);
  cluster = (int *) malloc(mapHeight * mapWidth * sizeof(int));
  for (int y = 0; y < mapHeight; ++y)
    for (int x = 0; x < mapWidth; ++x)
      cluster[y * mapWidth + x] = -1;

  seeds.clear();

  sprintf(fn, "%s%s.corr", corrName,
	  (pairs[index].pairName[0] != '\0' ?
	   pairs[index].pairName : pairs[index].imageName));
  f = fopen(fn, "r");
  if (f != NULL)
    {
      double oThreshold, dThreshold;
      if (fscanf(f, "%lf %lf", &oThreshold, &dThreshold) != 2)
	{
	  fprintf(stderr, "Invalid first line in %s\n", fn);
	  exit(1);
	}

      orthogonalSlider->value(oThreshold);
      diagonalSlider->value(dThreshold);
      while (fscanf(f, "%d %d", &x, &y) == 2 && x >= 0 && y >= 0)
	seeds.push_back(y*mapWidth+x);
      while (fscanf(f, "%d %d", &x, &y) == 2 && x >= 0 && y >= 0)
	mapMask[y*mapWidth+x] |= REJECT_BIT;
      while (fscanf(f, "%d %d", &x, &y) == 2 && x >= 0 && y >= 0)
	mapMask[y*mapWidth+x] |= ACCEPT_BIT;
      fclose(f);
    }
  else
    {
      selectSeed();
    }
  updateMask();

  printf("Set mapFactor to %d  (%d %d)\n", mapFactor, imageWidth, mapWidth);
  printf("image dim = (%d %d) imageMin = (%d %d)  mapOffset = (%d %d)\n",
	 imageWidth, imageHeight, imageMinX, imageMinY,
	 mapOffsetX, mapOffsetY);

  if (imageWindow != 0)
    imageWindow->computeWarped();

  lastIndex = index;
  sectionsLoaded = true;
}

void
MyWindow::saveSection ()
{
  if (readOnly)
    return;

  char fn[PATH_MAX];
  sprintf(fn, "%s%s.corr", corrName,
	  (pairs[index].pairName[0] != '\0' ?
	   pairs[index].pairName : pairs[index].imageName));
  FILE *f;
  while ((f = fopen(fn, "w")) == NULL)
    {
      fprintf(stderr, "Could not open file %s for writing.\n", fn);
      fprintf(stderr, "   Will retry operation in 60 seconds.\n");
      sleep(60);
    }
  fprintf(f, "%f %f\n",
	  orthogonalSlider->value(), diagonalSlider->value());
  for (int i = 0; i < seeds.size(); ++i)
    fprintf(f, "%d %d\n", seeds[i] % mapWidth, seeds[i] / mapWidth);
  fprintf(f, "-1 -1\n");
  for (int y = 0; y < mapHeight-1; ++y)
    for (int x = 0; x < mapWidth-1; ++x)
      if (mapMask[y*mapWidth+x] & REJECT_BIT)
	fprintf(f, "%d %d\n", x, y);
  fprintf(f, "-1 -1\n");
  for (int y = 0; y < mapHeight-1; ++y)
    for (int x = 0; x < mapWidth-1; ++x)
      if (mapMask[y*mapWidth+x] & ACCEPT_BIT)
	fprintf(f, "%d %d\n", x, y);
  fprintf(f, "-1 -1\n");
  fclose(f);

  MapElement *outMap = (MapElement*) malloc(mapHeight * mapWidth *
					    sizeof(MapElement));
  for (int y = 0; y < mapHeight; ++y)
    for (int x = 0; x < mapWidth; ++x)
      {
	outMap[y*mapWidth + x].x = map[y*mapWidth + x].x;
	outMap[y*mapWidth + x].y = map[y*mapWidth + x].y;
	if (x > 0 && y > 0 && mapMask[(y-1)*mapWidth + x-1] & CLUSTER_BIT ||
	    x > 0 && mapMask[y*mapWidth + x-1] & CLUSTER_BIT ||
	    y > 0 && mapMask[(y-1)*mapWidth + x] & CLUSTER_BIT ||
	    mapMask[y*mapWidth + x] & CLUSTER_BIT)
	  outMap[y*mapWidth + x].c = map[y*mapWidth + x].c;
	else
	  outMap[y*mapWidth + x].c = 0.0;
      }
  sprintf(fn, "%s%s.map", outputName,
	  (pairs[index].pairName[0] != '\0' ?
	   pairs[index].pairName : pairs[index].imageName));
  char errorMsg[PATH_MAX + 256];
  while (!WriteMap(fn, outMap, mapLevel, mapWidth, mapHeight,
		   mapOffsetX, mapOffsetY,
		   pairs[index].imageName,
		   (pairs[index].refName[0] != '\0' ?
		 pairs[index].refName : pairs[index].imageName),
		UncompressedMap, errorMsg))
    {
      fprintf(stderr, "Could not write map %s\n", fn);
      fprintf(stderr, "   Will retry operation in 60 seconds.\n");
      sleep(60);
    }
  free(outMap);
}

void
ImageWindow::computeWarped ()
{
  float xMin, xMax, yMin, yMax;
  float xRange, yRange;
  float scaleX, scaleY;
  float xp, yp;
  int x, y;
  float rx, ry, rc;
  int irx, iry;
  float rrx, rry;
  float r00, r01, r10, r11;
  float rv;
  float xv, yv;
  int b;

  int mapWidth = parent->mapWidth;
  int mapHeight = parent->mapHeight;
  MapElement* map = parent->map;
  int mapFactor = parent->mapFactor;
  
  xMin = 1.0e+30;
  xMax = -1.0e+30;
  yMin = 1.0e+30;
  yMax = -1.0e+30;
  for (y = 0; y < mapHeight; ++y)
    for (x = 0; x < mapWidth; ++x)
      {
	if (map[y*mapWidth + x].c == 0.0)
	  continue;
	xv = map[y*mapWidth + x].x * mapFactor / reductionFactor;
	yv = map[y*mapWidth + x].y * mapFactor / reductionFactor;
	if (xv < xMin)
	  xMin = xv;
	if (xv > xMax)
	  xMax = xv;
	if (yv < yMin)
	  yMin = yv;
	if (yv > yMax)
	  yMax = yv;
      }
  printf("xMin=%f xMax=%f yMin=%f yMax=%f\n",
	 xMin, xMax, yMin, yMax);
  xRange = xMax - xMin;
  xMin -= 0.05 * xRange;
  xMax += 0.05 * xRange;
  yRange = yMax - yMin;
  yMin -= 0.05 * yRange;
  yMax += 0.05 * yRange;
  scaleX = width / (xMax - xMin);
  scaleY = height / (yMax - yMin);
  scale = scaleX;
  if (scaleY < scale)
    scale = scaleY;
  xOffset = 0.5 * (width - scale * (xMin + xMax));
  yOffset = 0.5 * (height - scale * (yMin + yMax));
  printf("cw: scale=%f xo=%f yo=%f w=%d h=%d\n", scale, xOffset, yOffset,
	 width, height);

  if (warped != 0)
    free(warped);
  warped = (unsigned char *) malloc(width * height * sizeof(unsigned char));
  memset(warped, 0, width * height * sizeof(unsigned char));

  unsigned char* image = parent->image;
  if (image == 0)
    return;
  int imageWidth = parent->imageWidth;
  int imageHeight = parent->imageHeight;
  int imageMinX = parent->imageMinX;
  int imageMinY = parent->imageMinY;
  InverseMap* inverseMap = parent->inverseMap;
  float mapOffsetX = parent->mapOffsetX;
  float mapOffsetY = parent->mapOffsetY;
    
  for (y = 0; y < height; ++y)
    for (x = 0; x < width; ++x)
      {
	xp = ((x + 0.5 - xOffset) / scale) * reductionFactor / mapFactor;
	yp = ((y + 0.5 - yOffset) / scale) * reductionFactor / mapFactor;
	if (!Invert(inverseMap, &xv, &yv, xp, yp))
	  {
	    warped[y*width + x] = 0;
	    continue;
	  }
	rx = (xv + mapOffsetX) * mapFactor / reductionFactor - 0.5 - imageMinX;
	ry = (yv + mapOffsetY) * mapFactor / reductionFactor - 0.5 - imageMinY;
	irx = (int) floor(rx);
	iry = (int) floor(ry);
	if (irx < 0 || irx >= imageWidth-1 ||
	    iry < 0 || iry >= imageHeight-1)
	  {
	    warped[y*width + x] = 0;
	    continue;
	  }
	rrx = rx - irx;
	rry = ry - iry;
	r00 = image[iry*imageWidth + irx];
	r01 = image[(iry+1)*imageWidth + irx];
	r10 = image[iry*imageWidth + irx + 1];
	r11 = image[(iry+1)*imageWidth + irx + 1];
	rv = r00 * (rrx - 1.0) * (rry - 1.0)
	  - r10 * rrx * (rry - 1.0) 
	  - r01 * (rrx - 1.0) * rry
	  + r11 * rrx * rry;
	if (rv < -10.0 || rv > 1000.0)
	  {
	    printf("rv = %f  x = %d y = %d irx = %d iry = %d\n",
		   rv, x, y, irx, iry);
	    exit(1);
	  }
	b = (int) (rv);
	if (b < 0)
	  b = 0;
	else if (b > 255)
	  b = 255;
	warped[y * width + x] = b;
      }
  textureValid = false;

  printf("warped set: w=%d h=%d\n", width, height);
}

void
ImageWindow::initTexture (int w, int h)
{
  printf("initTexture called with w=%d h=%d\n", w, h);

  // initialize the texture
  if (texture > 0)
    glDeleteTextures(1, &texture);
  glGenTextures(1, &texture); 
  //  printf("textures = %d %d %d\n", texture[0], texture[1], texture[2]);
  textureWidth = w;
  textureHeight = h;
  unsigned char* black = (unsigned char *) malloc(textureWidth * textureHeight);
  memset(black, 0, textureWidth*textureHeight);
  glEnable(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D, texture);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_REPEAT); 
  glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE); 
  glBindTexture(GL_TEXTURE_2D, 0);
  glDisable(GL_TEXTURE_2D);
  free(black);
}

void
ImageWindow::setTexture (unsigned char *image, int w, int h)
{
  int x, y;
  unsigned char b;

  printf("ImageWindow::setTexture(%p, %d, %d)\n",
	 image, w, h);
  unsigned char *img = (unsigned char*) malloc(3 * textureWidth * textureHeight);
  unsigned char *p = img;
  for (y = 0; y < textureHeight; ++y)
    for (x = 0; x < textureWidth; ++x)
      {
	if (x < w && y < h)
	  b = image[y * w + x];
	else
	  b = 0;
	*p++ = b;
	*p++ = b;
	*p++ = b;
      }
  glBindTexture(GL_TEXTURE_2D, texture);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, textureWidth, textureHeight, 0,
	       GL_RGB, GL_UNSIGNED_BYTE, img);
  glBindTexture(GL_TEXTURE_2D, 0);
  free(img);
}

ImageWindow::ImageWindow(MyWindow *p, int w, int h) :
  Fl_Gl_Window (0, 0, w, h)
{
  parent = p;
  textureValid = false;
  textureInitialized = false;
  textureWidth = -1;
  textureHeight = -1;
  circleInitialized = false;
  warped = 0;
  width = -1;
  height = -1;
  name[0] = '\0';
}

void
MyWindow::draw ()
{
  if (w() != width || h() != height)
    reshape(w(), h());
  char title[128];
  sprintf(title, "map %s->%s", pairs[index].imageName, pairs[index].refName);
  label(title);
  Fl_Group::draw();
}

void
MyWindow::reshape (int w, int h)
{
  printf("MyWindow::reshape(%d, %d)\n", w, h);
  width = w;
  height = h;
  orthogonalSlider->resize(0, height - 2*sliderHeight,
			   width, sliderHeight);
  orthogonalSlider->redraw();
  diagonalSlider->resize(0, height - sliderHeight,
			 width, sliderHeight);
  diagonalSlider->redraw();
  imageWindow->resize(0, 0, width, height - 2 * sliderHeight);
  imageWindow->redraw();
}

void
MyWindow::selectSeed ()
{
  // find the top 10% of map elements in terms of low distortion score;
  vector<MapElementAndDistortion> m;
  MapElementAndDistortion me;
  int x, y;
  double distX, distY;
  double distA, distB;
  float diagonal, orthogonal;
  for (y = 0; y < mapHeight-1; ++y)
    for (x = 0; x < mapWidth-1; ++x)
      {
	if (map[y * mapWidth + x].c == 0.0 || map[y * mapWidth + x + 1].c == 0.0 ||
	    map[(y + 1) * mapWidth + x].c == 0.0 || map[(y + 1) * mapWidth + x + 1].c == 0.0)
	  continue;
	distX = 0.5 * hypot(map[y * mapWidth + x].x - map[y * mapWidth + x + 1].x,
			    map[y * mapWidth + x].y - map[y * mapWidth + x + 1].y) +
	  0.5 * hypot(map[(y + 1) * mapWidth + x].x - map[(y + 1) * mapWidth + x + 1].x,
		      map[(y + 1) * mapWidth + x].y - map[(y + 1) * mapWidth + x + 1].y);
	distY = 0.5 * hypot(map[y * mapWidth + x].x - map[(y + 1) * mapWidth + x].x,
			    map[y * mapWidth + x].y - map[(y + 1) * mapWidth + x].y) +
	  0.5 * hypot(map[y * mapWidth + x + 1].x - map[(y + 1) * mapWidth + x + 1].x,
		      map[y * mapWidth + x + 1].y - map[(y + 1) * mapWidth + x + 1].y);
	orthogonal = distY / distX;
	distA = hypot(map[y * mapWidth + x].x - map[(y + 1) * mapWidth + x + 1].x,
		      map[y * mapWidth + x].y - map[(y + 1) * mapWidth + x + 1].y);
	distB = hypot(map[y * mapWidth + x + 1].x - map[(y + 1) * mapWidth + x].x,
		      map[y * mapWidth + x + 1].y - map[(y + 1) * mapWidth + x].y);
	diagonal = distB / distA;
	me.x = x;
	me.y = y;
	me.distortion = fabs(orthogonal - avgOrthogonal) / avgOrthogonal +
	  fabs(diagonal - avgDiagonal) / avgDiagonal;
	m.push_back(me);
      }
  sort(m.begin(), m.end(), lesserDistortion);

  //   find the average location of those
  double sumX, sumY;
  double avgX, avgY;
  int n = (m.size() + 9) / 10;
  if (n == 0)
    return;
  for (int i = 0; i < n; ++i)
    {
      sumX += m[i].x;
      sumY += m[i].y;
    }
  avgX = sumX / n;
  avgY = sumY / n;

  //   pick the element that is located closest to the average
  double closestDist = 1.0e30;
  int closest = -1;
  double dist;
  for (int i = 0; i < n; ++i)
    {
      dist = hypot(m[i].x - avgX, m[i].y - avgY);
      if (dist < closestDist)
	{
	  closest = i;
	  closestDist = dist;
	}
    }
  if (closest >= 0)
    seeds.push_back(m[closest].y * mapWidth + m[closest].x);
}

void
MyWindow::updateMask ()
{
  int x, y;
  double distX, distY;
  double distA, distB;
  bool overThreshold;
  float diagonal, orthogonal;

  float oThreshold = orthogonalSlider->value();
  float dThreshold = diagonalSlider->value();
  for (y = 0; y < mapHeight-1; ++y)
    for (x = 0; x < mapWidth-1; ++x)
      {
	if (map[y*mapWidth+x].c == 0.0 || map[y*mapWidth+x+1].c == 0.0 ||
	    map[(y+1)*mapWidth+x].c == 0.0 || map[(y+1)*mapWidth+x+1].c == 0.0)
	  {
	    mapMask[y * mapWidth + x] = 0;
	    continue;
	  }
	mapMask[y * mapWidth + x] |= VALID_BIT;

	if (mapMask[y*mapWidth + x] & REJECT_BIT)
	  {
	    mapMask[y*mapWidth + x] &= ~THRESHOLD_BIT;
	    continue;
	  }
	if (mapMask[y*mapWidth + x] & ACCEPT_BIT)
	  {
	    mapMask[y*mapWidth + x] |= THRESHOLD_BIT;
	    continue;
	  }

	overThreshold = false;

	distX = 0.5 * hypot(map[y * mapWidth + x].x - map[y * mapWidth + x + 1].x,
			    map[y * mapWidth + x].y - map[y * mapWidth + x + 1].y) +
	  0.5 * hypot(map[(y + 1) * mapWidth + x].x - map[(y + 1) * mapWidth + x + 1].x,
		      map[(y + 1) * mapWidth + x].y - map[(y + 1) * mapWidth + x + 1].y);
	distY = 0.5 * hypot(map[y * mapWidth + x].x - map[(y + 1) * mapWidth + x].x,
			    map[y * mapWidth + x].y - map[(y + 1) * mapWidth + x].y) +
	  0.5 * hypot(map[y * mapWidth + x + 1].x - map[(y + 1) * mapWidth + x + 1].x,
		      map[y * mapWidth + x + 1].y - map[(y + 1) * mapWidth + x + 1].y);
	orthogonal = distY / distX;
	if (fabs(orthogonal - avgOrthogonal) / avgOrthogonal >= oThreshold)
	  overThreshold = true;
	distA = hypot(map[y * mapWidth + x].x - map[(y + 1) * mapWidth + x + 1].x,
		      map[y * mapWidth + x].y - map[(y + 1) * mapWidth + x + 1].y);
	distB = hypot(map[y * mapWidth + x + 1].x - map[(y + 1) * mapWidth + x].x,
		      map[y * mapWidth + x + 1].y - map[(y + 1) * mapWidth + x].y);
	diagonal = distB / distA;
	if (fabs(diagonal - avgDiagonal) / avgDiagonal >= dThreshold)
	  overThreshold = true;

	/* update mask */
	if (overThreshold)
	  mapMask[y * mapWidth + x] &= ~THRESHOLD_BIT;
	else
	  mapMask[y * mapWidth + x] |= THRESHOLD_BIT;
      }

  // clear all the cluster bits
  for (y = 0; y < mapHeight-1; ++y)
    for (x = 0; x < mapWidth-1; ++x)
      {
	mapMask[y * mapWidth + x] &= ~CLUSTER_BIT;
	cluster[y * mapWidth + x] = -1;
      }

  unsigned int *stack = (unsigned int *) malloc(mapHeight * mapWidth * sizeof(int));
  printf("updateMask: seeds.size = %zd\n", seeds.size());
  for (int i = 0; i < seeds.size(); ++i)
    {
      y = seeds[i] / mapWidth;
      x = seeds[i] - y * mapWidth;
      if (cluster[y * mapWidth + x] >= 0)
	continue;
      cluster[y * mapWidth + x] = i;
      int sp = 0;
      stack[sp++] = y*mapWidth+x;
      int count = 0;
      while (sp > 0)
	{
	  ++count;
	  --sp;
	  y = stack[sp] / mapWidth;
	  x = stack[sp] - y * mapWidth;
	  mapMask[y * mapWidth + x] |= CLUSTER_BIT;
	  
	  /* check up */
	  if (y > 0 && (mapMask[(y-1)*mapWidth+x] & THRESHOLD_BIT) &&
	      cluster[(y-1)*mapWidth+x] < 0)
	    {
	      cluster[(y-1)*mapWidth+x] = i;
	      stack[sp++] = (y-1)*mapWidth+x;
	    }
	  /* check right */
	  if (x < mapWidth-2 && (mapMask[y*mapWidth+x+1] & THRESHOLD_BIT) &&
	      cluster[y*mapWidth+x+1] < 0)
	    {
	      cluster[y*mapWidth+x+1] = i;
	      stack[sp++] = y*mapWidth+x+1;
	    }
	  /* check down */
	  if (y < mapHeight-2 && (mapMask[(y+1)*mapWidth+x] & THRESHOLD_BIT) &&
	      cluster[(y+1)*mapWidth+x] < 0)
	    {
	      cluster[(y+1)*mapWidth+x] = i;
	      stack[sp++] = (y+1)*mapWidth+x;
	    }
	  /* check left */
	  if (x > 0 && (mapMask[y*mapWidth+x-1] & THRESHOLD_BIT) &&
	      cluster[y*mapWidth+x-1] < 0)
	    {
	      cluster[y*mapWidth+x-1] = i;
	      stack[sp++] = y*mapWidth+x-1;
	    }
	}
      printf("flood fill marked %d elements\n", count);
    }
  free(stack);
}

void
ImageWindow::draw ()
{
  int x, y;
  float xv, yv;
  float xv0, yv0, xv1, yv1;
  int ix0, iy0, ix1, iy1;
  double distX, distY;
  double distA, distB;
  bool overThreshold;
  float diagonal, orthogonal;
  unsigned char m;
  int color0, color1;
  int color;
  int code;

  if (!valid())
    {
      printf("vid viewport = %d %d\n", w(), h());
      glViewport(0, 0, w(), h());

      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      glOrtho(0.0, w(), 0.0, h(), 0.0, 1.0);

      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
      gluLookAt(0.0, 0.0, 0.0,
                0.0, 0.0, -1.0,
                0.0, 1.0, 0.0);

      textureValid = false;
      if (!circleInitialized)
	initCircle();
    }

  if (w() != width || h() != height)
    {
      width = w();
      height = h();
      computeWarped();
    }

  // set the textures if necessary
  //  printf("draw called: tv=%d tw=%d th=%d\n",
  //	 textureValid ? 1 : 0, textureWidth, textureHeight);
  if (!textureValid)
    {
      // resize texture if necessary
      int newTextureWidth = 64;
      while (newTextureWidth < width)
	newTextureWidth <<= 1;
      int newTextureHeight = 64;
      while (newTextureHeight < height)
	newTextureHeight <<= 1;
      printf("ntw=%d nth=%d\n", newTextureWidth, newTextureHeight);
      if (newTextureWidth != textureWidth ||
	  newTextureHeight != textureHeight)
	initTexture(newTextureWidth, newTextureHeight);

      maxx = ((float) (width)) / textureWidth;
      maxy = ((float) (height)) / textureHeight;

      //      printf("setting textures %d %d %d %d %d %d\n", textureWidth, textureHeight, imageWidth, imageHeight, refWidth, refHeight);
      setTexture(warped, width, height);
      textureValid = true;
    }

  char name[128];
  char title[128];
  sprintf(name, "map %s", (pairs[parent->index].pairName[0] != '\0' ?
			   pairs[parent->index].pairName :
			   pairs[parent->index].imageName));
  label(title);

  //  printf("ImageWindow::draw  aw=%d ah=%d\n", width, height);
  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  glEnable(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D, texture);
  glBegin(GL_QUADS);
  glTexCoord2f(0.0, 0.0);
  glVertex2f(0.0, (float) (height));
  glTexCoord2f(maxx, 0.0);
  glVertex2f((float) (width), (float) (height));
  glTexCoord2f(maxx, maxy);
  glVertex2f((float) (width), 0.0);
  glTexCoord2f(0.0, maxy);
  glVertex2f(0.0, 0.0);
  glEnd();
  glDisable(GL_TEXTURE_2D);

  glBegin(GL_LINES);
  int mapWidth = parent->mapWidth;
  int mapHeight = parent->mapHeight;
  MapElement* map = parent->map;
  unsigned char *mapMask = parent->mapMask;
  int mapFactor = parent->mapFactor;
  for (y = 0; y < mapHeight; ++y)
    for (x = 0; x < mapWidth; ++x)
      {
	xv0 = scale * map[y*mapWidth+x].x * mapFactor / reductionFactor + xOffset;
	yv0 = scale * map[y*mapWidth+x].y * mapFactor / reductionFactor + yOffset;
	
	if (x < mapWidth-1 && map[y*mapWidth+x+1].c > 0.0)
	  {
	    /* draw horizontal line segment */
	    if (y > 0)
	      {
		m = mapMask[(y-1)*mapWidth+x];
		if (!(m & VALID_BIT))
		  color0 = -1;
		else if (m & REJECT_BIT)
		  color0 = 0; // red
		else if (m & ACCEPT_BIT)
		  color0 = 5; // blue
		else if (m & CLUSTER_BIT)
		  color0 = 2; // green
		else if (m & THRESHOLD_BIT)
		  color0 = 1; // yellow
		else
		  color0 = 4; // orange
	      }
	    else
	      color0 = -1;

	    if (y < mapHeight-1)
	      {
		m = mapMask[y*mapWidth+x];
		if (!(m & VALID_BIT))
		  color1 = -1;
		else if (m & REJECT_BIT)
		  color1 = 0; // red
		else if (m & ACCEPT_BIT)
		  color1 = 5; // blue
		else if (m & CLUSTER_BIT)
		  color1 = 2; // green
		else if (m & THRESHOLD_BIT)
		  color1 = 1; // yellow
		else
		  color1 = 4; // orange
	      }
	    else
	      color1 = -1;

	    if (color0 == 5 || color1 == 5)
	      color = 5;
	    else if (color0 == 0 || color1 == 0)
	      color = 0;
	    else if (color0 == 2 || color1 == 2)
	      color = 2;
	    else if (color0 == 4 || color1 == 4)
	      color = 4;
	    else if (color0 < 0 && color1 < 0)
	      color = -1;
	    else
	      color = 1;
	    if (color >= 0)
	      {
		glColor3f(colors[color][0] / 255.0,
			  colors[color][1] / 255.0,
			  colors[color][2] / 255.0);
		xv1 = scale * map[y*mapWidth+x+1].x * mapFactor / reductionFactor + xOffset;
		yv1 = scale * map[y*mapWidth+x+1].y * mapFactor / reductionFactor + yOffset;
		glVertex3f(xv0, height - yv0, -0.5);
		glVertex3f(xv1, height - yv1, -0.5);
	      }
	  }

	if (y < mapHeight-1)
	  {
	    /* draw vertical line segment */
	    if (x > 0)
	      {
		m = mapMask[y*mapWidth+x-1];
		if (!(m & VALID_BIT))
		  color0 = -1;
		else if (m & REJECT_BIT)
		  color0 = 0; // red
		else if (m & ACCEPT_BIT)
		  color0 = 5; // blue
		else if (m & CLUSTER_BIT)
		  color0 = 2; // green
		else if (m & THRESHOLD_BIT)
		  color0 = 1; // yellow
		else
		  color0 = 4; // orange
	      }
	    else
	      color0 = -1;

	    if (x < mapWidth-1)
	      {
		m = mapMask[y*mapWidth+x];
		if (!(m & VALID_BIT))
		  color1 = -1;
		else if (m & REJECT_BIT)
		  color1 = 0; // red
		else if (m & ACCEPT_BIT)
		  color1 = 5; // blue
		else if (m & CLUSTER_BIT)
		  color1 = 2; // green
		else if (m & THRESHOLD_BIT)
		  color1 = 1; // yellow
		else
		  color1 = 4; // orange
	      }
	    else color1 = -1;

	    if (color0 == 5 || color1 == 5)
	      color = 5;
	    else if (color0 == 0 || color1 == 0)
	      color = 0;
	    else if (color0 == 2 || color1 == 2)
	      color = 2;
	    else if (color0 == 4 || color1 == 4)
	      color = 4;
	    else if (color0 < 0 && color1 < 0)
	      color = -1;
	    else
	      color = 1;
	    if (color >= 0)
	      {
		glColor3f(colors[color][0] / 255.0,
			  colors[color][1] / 255.0,
			  colors[color][2] / 255.0);
		xv1 = scale * map[(y+1)*mapWidth+x].x * mapFactor / reductionFactor + xOffset;
		yv1 = scale * map[(y+1)*mapWidth+x].y * mapFactor / reductionFactor + yOffset;
		glVertex3f(xv0, height - yv0, -0.5);
		glVertex3f(xv1, height - yv1, -0.5);
	      }
	  }
      }
  glEnd();

  for (y = 0; y < mapHeight; ++y)
    for (x = 0; x < mapWidth; ++x)
      {
	if (map[y*mapWidth+x].c > 1.0)
	  {
	    code = (int) map[y*mapWidth+x].c;
	    for (int i = 0; i < 7; ++i)
	      {
		color = (code >> (3*i+1)) & 7;
		if (color == 0)
		  continue;
		xv = scale * map[y*mapWidth+x].x * mapFactor / reductionFactor + xOffset;
		yv = scale * map[y*mapWidth+x].y * mapFactor / reductionFactor + yOffset;
		if (i < 6)
		  {
		    xv += 5.0 * cos(i * 60.0 * M_PI / 180);
		    yv += 5.0 * sin(i * 60.0 * M_PI / 180);
		  }
		glPushMatrix();
		glTranslatef(xv, height - yv, -1.0);
		glScalef(3.0, 3.0, 3.0);
		glColor3f(colors[color][0]/255.0, colors[color][1]/255.0, colors[color][2]/255.0);
		glCallList(circle);
		glPopMatrix();
	      }
	  }
      }

#if 0
  glBegin(GL_LINES);
  glColor3f(1.0, 0.0, 0.0);
  glVertex3f(100.0, 100.0, -1.0);
  glVertex3f(200.0, 100.0, -1.0);
  glColor3f(0.0, 1.0, 0.0);
  glVertex3f(200.0, 200.0, -1.0);
  glVertex3f(300.0, 200.0, -1.0);
  glEnd();
#endif
  
  glFlush();
}

int
ImageWindow::handle (int event)
{
  int ex, ey;
  float xp, yp;
  float xv, yv;
  int x, y;

  switch (event)
    {
    case FL_PUSH:
    case FL_DRAG:
      switch (Fl::event_button())
	{
	case FL_LEFT_MOUSE:
	  ex = Fl::event_x();
	  ey = Fl::event_y();
	  xp = ((ex + 0.5 - xOffset) / scale) * reductionFactor /
	    parent->mapFactor;
	  yp = ((ey + 0.5 - yOffset) / scale) * reductionFactor /
	    parent->mapFactor;
	  if (!Invert(parent->inverseMap, &xv, &yv, xp, yp))
	    return(1);
	  x = (int) floor(xv);
	  y = (int) floor(yv);
	  if (x < 0 || x >= parent->mapWidth - 1 ||
	      y < 0 || y >= parent->mapHeight - 1)
	    return(1);

	  if (Fl::event_state() & FL_SHIFT)
	    {
	      if (Fl::event_state() & FL_CTRL)
		/* unmark map element for inclusion */
		parent->mapMask[y * (parent->mapWidth) + x] &= ~ACCEPT_BIT;
	      else
		/* mark map element for inclusion */
		parent->mapMask[y * (parent->mapWidth) + x] |= ACCEPT_BIT;
	    }
	  else
	    {
	      if (Fl::event_state() & FL_CTRL)
		/* unmark map element for exclusion */
		parent->mapMask[y * (parent->mapWidth) + x] &= ~REJECT_BIT;
	      else
		/* mark map element for exclusion */
		parent->mapMask[y * (parent->mapWidth) + x] |= REJECT_BIT;
	    }
	  parent->updateMask();
	  redraw();
	  break;

	case FL_MIDDLE_MOUSE:
	  if (event == FL_DRAG)
	    return(1);
	  ex = Fl::event_x();
	  ey = Fl::event_y();
	  xp = ((ex + 0.5 - xOffset) / scale) * reductionFactor /
	    parent->mapFactor;
	  yp = ((ey + 0.5 - yOffset) / scale) * reductionFactor /
	    parent->mapFactor;
	  if (!Invert(parent->inverseMap, &xv, &yv, xp, yp))
	    return(1);
	  x = (int) floor(xv);
	  y = (int) floor(yv);
	  if (x < 0 || x >= parent->mapWidth - 1 ||
	      y < 0 || y >= parent->mapHeight - 1)
	    return(1);

	  if (Fl::event_state() & FL_SHIFT)
	    {
	      /* remove seed */
	      int cluster = parent->cluster[y * parent->mapWidth + x];
	      if (cluster < 0)
		return(1);
	      vector<int> newSeeds;
	      for (int i = 0; i < parent->seeds.size(); ++i)
		if (parent->cluster[parent->seeds[i]] != cluster)
		  newSeeds.push_back(parent->seeds[i]);
	      parent->seeds = newSeeds;
	    }
	  else
	    {
	      /* create seed */
	      int seed = y * parent->mapWidth + x;
	      for (int i = 0; i < parent->seeds.size(); ++i)
		if (parent->seeds[i] == seed)
		  return(1);
	      printf("seeds size before = %zd\n", parent->seeds.size());
	      parent->seeds.push_back(seed);
	      printf("seeds size after = %zd\n", parent->seeds.size());
	    }
	  parent->updateMask();
	  redraw();
	  break;

	case FL_RIGHT_MOUSE:
	  if (event == FL_DRAG)
	    return(1);
	  if (Fl::event_state() & FL_SHIFT)
	    {
	      /* back up to the previous section */
	      parent->saveSection();
	      if (Fl::event_state() & FL_CTRL &&
		  Fl::event_state() & FL_ALT)
		parent->index = max(parent->index-100, 0);
	      else if (Fl::event_state() & FL_CTRL)
		parent->index = max(parent->index-10, 0);
	      else
		parent->index = max(parent->index-1, 0);
	      parent->increment = -1;
	      parent->loadSection();
	      parent->imageWindow->textureValid = false;
	      parent->imageWindow->redraw();
	      parent->redraw();
	    }
	  else
	    {
	      /* go to the next section */
	      parent->saveSection();
	      if (Fl::event_state() & FL_CTRL &&
		  Fl::event_state() & FL_ALT)
		parent->index = min(parent->index+100, nPairs-1);
	      else if (Fl::event_state() & FL_CTRL)
		parent->index = min(parent->index+10, nPairs-1);
	      else
		parent->index = min(parent->index+1, nPairs-1);
	      parent->increment = 1;
	      parent->loadSection();
	      parent->imageWindow->textureValid = false;
	      parent->imageWindow->redraw();
	      parent->redraw();
	    }
	  break;
	default:
	  break;
	}
      break;
    default:
      return(Fl_Gl_Window::handle(event));
    }
  return(1);
}

void
OrthogonalSliderCallback (Fl_Widget *w, void *data)
{
  Fl_Value_Slider *s = (Fl_Value_Slider*) w;

  s->value();
  window->updateMask();
  window->imageWindow->redraw();
}

void
DiagonalSliderCallback (Fl_Widget *w, void *data)
{
  Fl_Value_Slider *s = (Fl_Value_Slider*) w;

  s->value();
  window->updateMask();
  window->imageWindow->redraw();
}

void
ImageWindow::initCircle ()
{
  int i;
  float x, y;

  circle = 1;
  glNewList(circle, GL_COMPILE);
  glBegin(GL_TRIANGLE_FAN);
  glVertex2f(0.0, 0.0);
  for (i = 0; i <= 16; ++i)
    {
      x = cos(i * 2 * M_PI / 16.0);
      y = sin(i * 2 * M_PI / 16.0);
      glVertex2f(x, y);
    }
  glEnd();
  glEndList();
}
