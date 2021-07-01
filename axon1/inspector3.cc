/*
 * inspector3.cc  - display pairs of images so that
 *                  alignment maps can be visually checked
 *
 *   Copyright (c) 2008-2010 Pittsburgh Supercomputing Center
 * 
 *  This program is distributed in the hope that it will    *
 *  be useful, but WITHOUT ANY WARRANTY; without even the   *
 *  implied warranty of MERCHANTABILITY or FITNESS FOR A    *
 *  PARTICULAR PURPOSE.  Neither Carnegie Mellon University *
 *  nor any of the authors assume any liability for         *
 *  damages, incidental or otherwise, caused by the         *
 *  installation or use of this software.                   *
 *                                                          *
 *  CLINICAL APPLICATIONS ARE NOT RECOMMENDED, AND THIS     *
 *  SOFTWARE HAS NOT BEEN EVALUATED BY THE UNITED STATES    *
 *  FDA FOR ANY CLINICAL USE.                               *
 *
 *  HISTORY
 *    2008-2009  Written by Greg Hood (ghood@psc.edu)
 *    2010  Fixes to correctly display pairs involving
 *           subimages with nonzero offsets (ghood@psc.edu)
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
#include <FL/fl_draw.H>
#include <GL/glu.h>
#include <string>
#include <algorithm>
#include "imio.h"
#include "invert.h"

using std::string;
using std::max;
using std::min;

int fmDebug = 0;

#define N_LEVELS	16       /* max image size (w or h) is 2^N_LEVELS */s
#define MAP(map,x,y)	map[(y)*mapWidth + x]
#define IMAGE(i,x,y)    (x >= 0 && x < w && y >= 0 && y < h ? i[(y)*w + (x)] : background)
#define RIMAGE(i,x,y)   (x >= 0 && x < rw && y >= 0 && y < rh ? i[(y)*rw + (x)] : background)
#define LINE_LENGTH	255

struct CPoint
{
  float ix, iy;    // point in image
  float rx, ry;	   // point in reference
  int color;       // what color to use for these points
  float wrx, wry;  // point in image that maps to the reference point
};

struct Pair
{
  char *imageName;
  int imageMinX, imageMaxX, imageMinY, imageMaxY; 
  char *refName;
  int refMinX, refMaxX, refMinY, refMaxY;
  char *pairName;
  bool modified;
};

class Window: public Fl_Gl_Window
{
 public:
  Window (int x, int y, int w, int h);
  void draw ();
  int handle (int);
  void loadSection ();
  void saveSection ();
  void setTexture (unsigned int type, unsigned char *image,
		   int w, int h);
  bool backmap (float *rx, float *ry, float x, float y);
  bool forwardmap (float *px, float *py, float *pc, float x, float y);
  void reshape ();

 protected:
  GLuint texture[3];
  bool circleInitialized;

 private:
  void initTextures(int w, int h);
  void initCircle();

  bool setTextures;       // if true, the textures need to be set before drawing
  int textureWidth;       // the width in pixels of the texture map we are using
  int textureHeight;      // the height in pixels of the texture map we are using
  float maxx;             // the portion of the texture width we are actually
                          //   using (0 < maxx <= 1.0)
  float maxy;             // the portion of the texture height we are actually
                          //   using (0 < maxy <= 1.0)

  int type;               // the type of the current image displayed in the
                          //     window:
                          //   0: the original image
                          //   1: the warped reference
                          //   2: the original reference
  int index;              // the index within the pair list of the current pair
  int lastIndex;          // the last "good" index
  int increment;          // the last increment to index (-1 or 1)
  bool sectionsLoaded;    // true if we have already loaded the currentPair
  bool displayOriginal;   // true if we should display the original reference
                          //   rather than the warped reference
  bool displayMap;        // true if we should display the map superimposed
                          //   on the reference

  bool maxSize;       // true if the window should be kept to the maximum
                      //   size consistent with the image's aspect ratio
  float scale;        // how many pixels on the screen correspond to one pixel in
                      //   the image
  float offsetX;      // the x pixel position in the window that the left side of the
                      //   0 pixel position in the image maps to 
  float offsetY;      // the y pixel position in the window that the upper side of the
                      //   0 pixel position in the image maps to
  int dragX;	      // cursor X position when drag started
  int dragY;	      // cursor Y position when drag started
  float dragScale;    // scale value when drag started
  float dragOffsetX;  // offsetX value when drag started
  float dragOffsetY;  // offsetY value when drag started
  float dragRadius;   // distance from cursor to center of window when drag started
  GLuint circle;      // display list for a circle

  int nCpts;          // number of corresponding points
  CPoint *cpts;       // an array of the corresponding points

  int imageWidth, imageHeight;
  int imageMinX, imageMaxX, imageMinY, imageMaxY;
  int refWidth, refHeight;
  int refMinX, refMaxX, refMinY, refMaxY;
  int displayWidth, displayHeight;
  int windowWidth, windowHeight;  // desired window width and height
  int reductionFactor;

  unsigned char *image;
  unsigned char *ref;
  unsigned char *warped;

  int mapWidth, mapHeight;
  int mapOffsetX, mapOffsetY;
  int mapFactor;
  MapElement *map;
  InverseMap *inverseMap;
};

char pairsFile[PATH_MAX];
char inputName[PATH_MAX];
char mapName[PATH_MAX];
char cptsName[PATH_MAX];
float background = 0.0;
Window* window;
Pair* pairs;
int nPairs;
bool readOnly = false;
char extension[8] = "pgm";
int baseReductionFactor = 1;
bool partial = false;
bool gray = true;

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

int availableHeight;
int availableWidth;

int ParseRange (char *s, int *pos, int *minValue, int *maxValue);
int ParseValue (char *s, int *pos, int *value);
int ReadHeader (FILE *f, char *tc, int *w, int *h, int *m);
void Beep();

#if 0
extern "C" {
  int ReadImage (char *filename, unsigned char **pixels,
		 int *width, int *height,
		 int minX, int maxX, int minY, int maxY,
		 char *error);

  void Log(char *fmt, ...)
  {
    va_list args;
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
  }
}
#endif


int MyReadMap (char *filename,
	     MapElement** map,
	     int *level,
	     int *width, int *height,
	     int *xMin, int *yMin,
	     char *imageName, char *referenceName,
	     char *error)
{
  char imgName[PATH_MAX], refName[PATH_MAX];
  int newWidth, newHeight;

  FILE *f = fopen(filename, "r");
  if (f == NULL)
    {
      sprintf(error, "Cannot open file %s\n", filename);
      return(0);
    }
  if (fgetc(f) != 'M' || fgetc(f) != '1' || fgetc(f) != '\n' ||
      fscanf(f, "%d%d%d%d%d%s%s",
	     level,
	     &newWidth, &newHeight,
	     xMin, yMin,
	     imgName, refName) != 7 ||
      fgetc(f) != '\n')
    {
      sprintf(error, "Cannot read header of map file %s\n", filename);
      return(0);
    }
  if (imageName != NULL)
    strcpy(imageName, imgName);
  if (referenceName != NULL)
    strcpy(referenceName, refName);
  if (*width <= 0 || *height <= 0)
    *map = NULL;
  if (*map != NULL && (newWidth != *width || newHeight != *height))
    {
      free(*map);
      *map = NULL;
    }
  if (*map == NULL)
    *map = (MapElement *) malloc(newWidth * newHeight * sizeof(MapElement));
  if (*map == NULL)
    {
      sprintf(error, "Could not allocate space for map %s\n", filename);
      return(0);
    }
  *width = newWidth;
  *height = newHeight;
  if (fread(*map, sizeof(MapElement), newWidth * newHeight, f) != newWidth * newHeight)
    {
      sprintf(error, "Could not read map from file %s\n", filename);
      return(0);
    }
  fclose(f);
  return(1);
}

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
  int minS, maxS;
  int minSection, maxSection;
  DIR *dir;
  struct dirent *de;
  char fn[PATH_MAX];
  int len;
  int z;
  bool maximizeWindow = false;
  int windowFactor = 1;
  int windowWidth, windowHeight;
  char imgn[PATH_MAX], refn[PATH_MAX], pairn[PATH_MAX];
  int imgMinX, imgMaxX, imgMinY, imgMaxY;
  int refMinX, refMaxX, refMinY, refMaxY;
  char line[LINE_LENGTH+1];

  error = 0;
  pairsFile[0] = '\0';
  inputName[0] = '\0';
  mapName[0] = '\0';
  cptsName[0] = '\0';
  minSection = 0;
  maxSection = 1000000000;
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
    else if (strcmp(argv[i], "-map") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-map error\n");
	    break;
	  }
	strcpy(mapName, argv[i]);
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
    else if (strcmp(argv[i], "-cpts") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-cpts error\n");
	    break;
	  }
	strcpy(cptsName, argv[i]);
      }
    else if (strcmp(argv[i], "-reduction") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d", &baseReductionFactor) != 1)
	  {
	    error = 1;
	    fprintf(stderr, "-reduction error\n");
	    break;
	  }
      }
    else if (strcmp(argv[i], "-readonly") == 0)
      readOnly = true;
    else if (strcmp(argv[i], "-tif") == 0)
      strcpy(extension, "tif");
    else if (strcmp(argv[i], "-partial") == 0)
      partial = true;
    else if (strcmp(argv[i], "-no_gray") == 0)
      gray = false;
    else
      error = 1;

  if (error)
    {
      fprintf(stderr, "Usage: inspector3 -input file_prefix -map file_prefix\n");
      fprintf(stderr, "                  -pairs pairfile\n");
      fprintf(stderr, "                  [-cpts cpts_prefix] [-readonly]\n");
      fprintf(stderr, "                  [-reduction factor]\n");
      fprintf(stderr, "                  [-tif]\n");
      fprintf(stderr, "                  [-partial]\n");
      fprintf(stderr, "   where ranges are expressed as: integer\n");
      fprintf(stderr, "                              or: integer-integer\n");
      fprintf(stderr, "                              or: integer-\n");
      fprintf(stderr, "                              or: -integer\n");
      exit(1);
    }

  /* check that at least minimal parameters were supplied */
  if (inputName[0] == '\0' || mapName[0] == '\0' || pairsFile[0] == '\0')
    {
      fprintf(stderr, "-input, -map, and -pairs parameters must be specified.\n");
      exit(1);
    }

  pairs = 0;
  f = fopen(pairsFile, "r");
  if (f == NULL)
    {
      fprintf(stderr, "Could not open pairs file %s\n", pairsFile);
      exit(1);
    }
  while (fgets(line, LINE_LENGTH, f) != NULL)
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

  /*  Fl::lock(); */
  Fl::visual(FL_RGB);
  Fl::gl_visual(FL_RGB);
  availableHeight = Fl::h() - 20;
  availableWidth = Fl::w() - 10;
  //  availableHeight = 925;
  //  availableWidth = 1270;
  window = new Window(0, 0, 256, 256);
  window->loadSection();
  window->reshape();
  window->show();
  printf("Windows created\n");
  return(Fl::run());
}

Window::Window (int x, int y, int w, int h) : Fl_Gl_Window (x, y, w, h)
{
  texture[0] = texture[1] = texture[2] = 0;
  textureWidth = -1;
  textureHeight = -1;
  setTextures = true;
  circleInitialized = false;
  type = 0;
  index = 0;
  lastIndex = -1;
  increment = 1;
  sectionsLoaded = false;

  displayOriginal = (mapName[0] == '\0');
  displayMap = 0;
  dragX = 0;
  dragY = 0;
  dragScale = scale;
  dragOffsetX = 0.0;
  dragOffsetY = 0.0;
  dragRadius = 10.0;

  nCpts = 0;
  cpts = 0;

  mapWidth = 0;
  mapHeight = 0;

  image = 0;
  ref = 0;
  warped = 0;

  mapFactor = 0;
  map = 0;
}

void
Window::reshape ()
{
  if (windowWidth != this->w() || windowHeight != this->h())
    size(windowWidth, windowHeight);
}

void
Window::draw ()
{
  int x, y;
  float xv, yv;
  float xv0, yv0, xv1, yv1;
  int ix0, iy0, ix1, iy1;

  //  printf("Entering draw... nCpts = %d\n", nCpts);
  if (!valid())
    {
      if (!circleInitialized)
	initCircle();

      //      printf("vid viewport = %d %d\n", w(), h());
      glViewport(0, 0, w(), h());

      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      glOrtho(0.0, w(), 0.0, h(), 0.0, 1.0);

      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
      gluLookAt(0.0, 0.0, 0.0,
                0.0, 0.0, -1.0,
                0.0, 1.0, 0.0);

      setTextures = true;
    }

  // set the textures if necessary
  if (setTextures)
    {
      //      printf("Entering setTextures... nCpts = %d\n", nCpts);

      // resize texture if necessary
      int newTextureWidth = 64;
      while (newTextureWidth < displayWidth)
	newTextureWidth <<= 1;
      int newTextureHeight = 64;
      while (newTextureHeight < displayHeight)
	newTextureHeight <<= 1;
      if (newTextureWidth != textureWidth ||
	  newTextureHeight != textureHeight)
	initTextures(newTextureWidth, newTextureHeight);

      maxx = ((float) displayWidth) / textureWidth;
      maxy = ((float) displayHeight) / textureHeight;

      //      printf("setting textures %d %d %d %d %d %d\n", textureWidth, textureHeight, imageWidth, imageHeight, refWidth, refHeight);
      setTexture(0, image, imageWidth, imageHeight);
      setTexture(1, warped, imageWidth, imageHeight);
      setTexture(2, ref, refWidth, refHeight);
      setTextures = false;

      //      printf("Leaving setTextures... nCpts = %d\n", nCpts);
    }

  char name[128];
  char title[128];
  if (mapName[0] != '\0')
    {
      switch (type)
	{
	case 0:
	  sprintf(name, "image %s", pairs[index].imageName);
	  break;
	case 1:
	  sprintf(name, "image %s warped to %s",
		  pairs[index].refName,
		  pairs[index].imageName);
	  break;
	case 2:
	  sprintf(name, "image %s", pairs[index].refName);
	  break;
	}
      sprintf(title, "%s Pair %s: %s   \n",
	      pairs[index].modified ? "**" : "  ",
	      pairs[index].pairName,
	      name);
    }
  else
    {
      switch (type)
	{
	case 0:
	  sprintf(title, "image %s", pairs[index].imageName);
	  break;
	case 2:
	  sprintf(title, "image %s", pairs[index].refName);
	  break;
	}
    }
  label(title);

  glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  glEnable(GL_TEXTURE_2D);
  //  printf("type = %d, texture[type] = %d\n", type, (int) texture[type]);
  glBindTexture(GL_TEXTURE_2D, texture[type]);
  glBegin(GL_QUADS);
  glTexCoord2f(0.0, 0.0);
  //  glTexCoord2f(0.0, maxy);
  glVertex2f(offsetX, h() - offsetY);
  glTexCoord2f(maxx, 0.0);
  //  glTexCoord2f(maxx, 1.0 - maxy);
  glVertex2f(scale * displayWidth + offsetX,
	     h() - offsetY);
  glTexCoord2f(maxx, maxy);
  //  glTexCoord2f(maxx, 1.0);
  glVertex2f(scale * displayWidth + offsetX,
	     h() - (scale * displayHeight + offsetY));
  glTexCoord2f(0.0, maxy);
  //  glTexCoord2f(0.0, 1.0);
  glVertex2f(offsetX,
	     h() - (scale * displayHeight + offsetY));
  //  printf("quad: (%f %f) (%f %f) (%f %f) (%f %f)\noffsetX = %f offsetY = %f\nmaxx = %f maxy = %f\ntw = %d th = %d\n",
  //	 offsetX, h() - offsetY,
  //	 scale * displayWidth + offsetX, h() - offsetY,
  //	 scale * displayWidth + offsetX, h() - (scale * displayHeight + offsetY),
  //	 offsetX, h() - (scale * displayHeight + offsetY),
  //	 offsetX, offsetY,
  //	 maxx, maxy,
  //	 textureWidth, textureHeight);
  glEnd();
  glDisable(GL_TEXTURE_2D);

  if (displayMap && type == 2)
    {
      glColor3f(0.0, 0.0, 0.5);
      glBegin(GL_LINES);
      for (y = 0; y < mapHeight; ++y)
	for (x = 0; x < mapWidth; ++x)
	  {
	    if (MAP(map, x, y).c == 0.0)
	      continue;
	    xv0 = scale * MAP(map, x, y).x * mapFactor / reductionFactor + offsetX;
	    yv0 = scale * MAP(map, x, y).y * mapFactor / reductionFactor + offsetY;
	    
	    if (x < mapWidth-1 && MAP(map, x+1, y).c != 0.0)
	      {
		/* draw horizontal line segment */
		xv1 = scale * MAP(map, x+1, y).x * mapFactor / reductionFactor + offsetX;
		yv1 = scale * MAP(map, x+1, y).y * mapFactor / reductionFactor + offsetY;
		glVertex3f(xv0, h() - yv0, -0.5);
		glVertex3f(xv1, h() - yv1, -0.5);
	      }

	    if (y < mapHeight-1 && MAP(map, x, y+1).c != 0.0)
	      {
		/* draw vertical line segment */
		xv1 = scale * MAP(map, x, y + 1).x * mapFactor / reductionFactor + offsetX;
		yv1 = scale * MAP(map, x, y + 1).y * mapFactor / reductionFactor + offsetY;
		glVertex3f(xv0, h() - yv0, -0.5);
		glVertex3f(xv1, h() - yv1, -0.5);
	      }
	  }
      glEnd();
    }

#if 0
  glBegin(GL_LINES);
  glColor3f(1.0, 0.0, 0.0);
  glVertex3f(100.0, 100.0, -1.0);
  glVertex3f(200.0, 100.0, -1.0);
  glEnd();
#endif
  
#if 1
//  printf("nCpts = %d\n", nCpts);
  for (int i = 0; i < nCpts; ++i)
    {
      bool valid;
      switch (type)
	{
	case 0:
	  valid = (cpts[i].ix >= 0);
	  xv = scale * ((cpts[i].ix - imageMinX) / reductionFactor + 0.5) + offsetX;
	  yv = scale * ((cpts[i].iy - imageMinY) / reductionFactor + 0.5) + offsetY;
	  break;
	case 1:
	  valid = (cpts[i].wrx >= 0);
	  xv = scale * ((cpts[i].wrx - imageMinX) / reductionFactor + 0.5) + offsetX;
	  yv = scale * ((cpts[i].wry - imageMinY) / reductionFactor + 0.5) + offsetY;
	  //	  printf("valid %d x %f y %f\n", valid, x, y);
	  break;
	case 2:
	  valid = (cpts[i].rx >= 0);
	  xv = scale * ((cpts[i].rx - refMinX) / reductionFactor + 0.5) + offsetX;
	  yv = scale * ((cpts[i].ry - refMinY) / reductionFactor + 0.5) + offsetY;
	  break;
	}
      if (valid)
	{
	  glPushMatrix();
	  glTranslatef(xv, h() - yv, -1.0);
	  glScalef(5.0, 5.0, 1.0);
	  int ci = cpts[i].color & 15;
	  glColor3f(colors[ci][0]/255.0, colors[ci][1]/255.0, colors[ci][2]/255.0);
	  glCallList(circle);
	  glPopMatrix();
	}
    }
#endif
  glFlush();
}

/*              LEFT           MIDDLE             RIGHT
NONE           toggle          place_point        next image
SHIFT          toggle/orig     remove_point       prev image   
CTRL           move
SHIFT+CTRL    reduce/enlarge   
*/

int
Window::handle (int event)
{
  //  printf("Entering handle... nCpts = %d\n", nCpts);
  switch (event)
    {
    case FL_PUSH:
      switch (Fl::event_button())
	{
	case FL_LEFT_MOUSE:
	  if (Fl::event_state() & FL_CTRL)
	    {
	      dragX = Fl::event_x();
	      dragY = Fl::event_y();
	      if (Fl::event_state() & FL_SHIFT)
		{
		  dragRadius = sqrt((double) ((dragX - w()/2) * (dragX - w()/2)) +
				    (double) ((dragY - h()/2) * (dragY - h()/2)));
		  if (dragRadius < 10.0)
		    dragRadius = 10.0;
		  dragScale = scale;
		}
	      dragOffsetX = offsetX;
	      dragOffsetY = offsetY;
	    }
	  else if (Fl::event_state() & FL_SHIFT)
	    {
	      if (mapName[0] != '\0')
		{
		  if (displayOriginal)
		    displayMap = !displayMap;
		  displayOriginal = !displayOriginal;
		  type = displayOriginal ? 2 : 1;
		}
	      redraw();
	    }
	  else
	    {
	      if (type == 0)
		type = (displayOriginal || mapName[0] == '\0') ? 2 : 1;
	      else
		type = 0;
	      redraw();
	    }
	  break;
	case FL_MIDDLE_MOUSE:
	  if (mapName[0] == '\0')
	    break;
	  if (Fl::event_state() & FL_SHIFT)
	    {
	      /* delete correspondence point */
	      int ex = Fl::event_x();
	      int ey = Fl::event_y();
	      float rx, ry, rc;
	      float x, y;

	      /* find the point that is closest */
	      float minDist = 1000000000.0;
	      float dist;
	      int minPt;
	      switch (type)
		{
		case 0:
		  for (int i = 0; i < nCpts; ++i)
		    {
		      x = scale * cpts[i].ix / reductionFactor + offsetX;
		      y = scale * cpts[i].iy / reductionFactor + offsetY;
		      dist = hypot(x - (ex + 0.5), y - (ey + 0.5));
		      if (dist < minDist)
			{
			  minDist = dist;
			  minPt = i;
			}
		    }
		  break;
		case 1:
		  for (int i = 0; i < nCpts; ++i)
		    {
		      if (!forwardmap(&rx, &ry, &rc, cpts[i].rx, cpts[i].ry))
			continue;
		      x = scale * rx / reductionFactor + offsetX;
		      y = scale * ry / reductionFactor + offsetY;
		      dist = hypot(x - (ex + 0.5), y - (ey + 0.5));
		      if (dist < minDist)
			{
			  minDist = dist;
			  minPt = i;
			}
		    }
		  break;
		case 2:
		  for (int i = 0; i < nCpts; ++i)
		    {
		      x = scale * cpts[i].rx / reductionFactor + offsetX;
		      y = scale * cpts[i].ry / reductionFactor + offsetY;
		      dist = hypot(x - (ex + 0.5), y - (ey + 0.5));
		      if (dist < minDist)
			{
			  minDist = dist;
			  minPt = i;
			}
		    }
		  break;
		}
	      if (minDist > 5.0)
		break;
	      for (int i = minPt + 1; i < nCpts; ++i)
		cpts[i-1] = cpts[i];
	      --nCpts;
	      cpts = (CPoint*) realloc(cpts, nCpts * sizeof(CPoint));
	    }
	  else
	    {
	      /* create correspondence point */
	      int ex = Fl::event_x();
	      int ey = Fl::event_y();
	      printf("Adding corr pt at window pixel %d %d\n", ex, ey);
	      float x, y;
	      int ix, iy;

	      /* check if we should create a new pair */
	      if (nCpts == 0 || (cpts[nCpts-1].ix >= 0 && cpts[nCpts-1].rx >= 0))
		{
		  cpts = (CPoint*) realloc(cpts, (nCpts+1) * sizeof(CPoint));
		  cpts[nCpts].ix = -1;
		  cpts[nCpts].iy = -1;
		  cpts[nCpts].rx = -1;
		  cpts[nCpts].ry = -1;
		  cpts[nCpts].wrx = -1;
		  cpts[nCpts].wry = -1;
		  int minUnusedColor;
		  for (minUnusedColor = 0;
		       ;
		       ++minUnusedColor)
		    {
		      int i;
		      for (i = 0; i < nCpts; ++i)
			if (cpts[i].color == minUnusedColor)
			  break;
		      if (i == nCpts)
			break;
		    }
		  cpts[nCpts].color = minUnusedColor;
		  ++nCpts;
		  //		  printf("Incremented nCpts to %d\n", nCpts);
		}
	      switch (type)
		{
		case 0:
		  cpts[nCpts-1].ix = (int) floor(((ex - offsetX) / scale) * reductionFactor + imageMinX + 0.5);
		  cpts[nCpts-1].iy = (int) floor(((ey - offsetY) / scale) * reductionFactor + imageMinY + 0.5);
		  printf("Marking point (%f %f) in image\n",
			 cpts[nCpts-1].ix, cpts[nCpts-1].iy); 
		  break;

		case 1:
		  {
		    float rx, ry, rc;
		    x = (ex - offsetX) / scale;
		    y = (ey - offsetY) / scale;
		    cpts[nCpts-1].wrx = x * reductionFactor + imageMinX;
		    cpts[nCpts-1].wry = y * reductionFactor + imageMinY;
		    if (forwardmap(&rx, &ry, &rc,
				   cpts[nCpts-1].wrx, cpts[nCpts-1].wry))
		      {
			cpts[nCpts-1].rx = (int) floor(rx + 0.5);
			cpts[nCpts-1].ry = (int) floor(ry + 0.5);
			printf("Marking point (%f %f) in reference\n",
			       cpts[nCpts-1].rx, cpts[nCpts-1].ry);
		      }
		    else
		      {
			printf("failed to set rx,ry\n");
			cpts[nCpts-1].rx = -1;
			cpts[nCpts-1].ry = -1;
		      }
		  }
		  break;

		case 2:
		  x = (ex - offsetX) / scale;
		  y = (ey - offsetY) / scale;
		  cpts[nCpts-1].rx = x * reductionFactor + refMinX;
		  cpts[nCpts-1].ry = y * reductionFactor + refMinY;
		  printf("Marking point (%f %f) in reference\n",
			 cpts[nCpts-1].rx,
			 cpts[nCpts-1].ry);
		  if (!backmap(&cpts[nCpts-1].wrx, &cpts[nCpts-1].wry,
			       (float) cpts[nCpts-1].rx,
			       (float) cpts[nCpts-1].ry))
		    {
		      printf("failed to set wrx,wry because backmap failed\n");
		      cpts[nCpts-1].wrx = -1.0;
		      cpts[nCpts-1].wry = -1.0;
		    }
		  break;
		}
	    }
	  redraw();
	  break;
	case FL_RIGHT_MOUSE:
	  if (mapName[0] != '\0')
	    displayOriginal = false;
	  if (Fl::event_state() & FL_SHIFT)
	    {
	      /* back up to the previous section */
	      saveSection();
	      if (Fl::event_state() & FL_CTRL &&
		  Fl::event_state() & FL_ALT)
		index = max(index-100, 0);
	      else if (Fl::event_state() & FL_CTRL)
		index = max(index-10, 0);
	      else
		index = max(index-1, 0);
	      increment = -1;
	      loadSection();
	      offsetX = 0.0;
	      offsetY = 0.0;
	      reshape();
	      redraw();
	    }
	  else
	    {
	      /* go to the next section */
	      saveSection();
	      if (Fl::event_state() & FL_CTRL &&
		  Fl::event_state() & FL_ALT)
		index = min(index+100, nPairs-1);
	      else if (Fl::event_state() & FL_CTRL)
		index = min(index+10, nPairs-1);
	      else
		index = min(index+1, nPairs-1);
	      increment = 1;
	      loadSection();
	      offsetX = 0.0;
	      offsetY = 0.0;
	      reshape();
	      redraw();
	    }
	  break;
	default:
	  break;
	}
      break;
    case FL_RELEASE:
      break;
    case FL_DRAG:
      switch (Fl::event_button())
	{
	case FL_LEFT_MOUSE:
	  if (Fl::event_state() & FL_CTRL)
	    {
	      int ex = Fl::event_x();
	      int ey = Fl::event_y();
	      if (Fl::event_state() & FL_SHIFT)
		{
		  float radius = hypot((double) (ex - w()/2),
				       (double) (ey - h()/2));
		  if (radius < 10.0)
		    radius = 10.0;
		  scale = dragScale * radius / dragRadius;
		  offsetX = w()/2 - scale / dragScale * (w()/2 - dragOffsetX);
		  offsetY = h()/2 - scale / dragScale * (h()/2 - dragOffsetY);
		}
	      else
		{
		  offsetX = dragOffsetX + (ex - dragX);
		  offsetY = dragOffsetY + (ey - dragY);
		}
	      redraw();
	    }
	  break;
	default:
	  break;
	}
      break;
    default:
      return(Fl_Gl_Window::handle(event));
    }
  //  printf("Leaving handle... nCpts = %d\n", nCpts);
  return(1);
}

void
Window::loadSection ()
{
  int level;
  int w, h;
  int n;
  int i;
  int x, y;
  float *warp;
  FILE *f;
  int refSection;
  float rx, ry, rc;
  int irx, iry;
  float rrx, rry;
  float r00, r01, r10, r11;
  float rv;
  char fn[PATH_MAX];
  char tc;
  int m;
  int b;
  float xv, yv;
  int ixv, iyv;
  float rx00, rx01, rx10, rx11, ry00, ry01, ry10, ry11;
  int rw, rh;
  char errorMsg[PATH_MAX + 256];
  char imgn[PATH_MAX], refn[PATH_MAX];
  int iv;
  int nw, nh;
  int count;
  int newImageMinX, newImageMaxX, newImageMinY, newImageMaxY;
  int newRefMinX, newRefMaxX, newRefMinY, newRefMaxY;
  int deltaX, deltaY, delta;

 retry:
  if (index < 0)
    {
      fprintf(stderr, "Could not find a valid image pair to display.\n");
      exit(1);
    }
  if (pairs[index].refName[0] != '\0')
    printf("Loading sections %s and %s  (pair %s)  (%p %p %p %p)\n", pairs[index].imageName, pairs[index].refName,
	   pairs[index].pairName,
	   image, ref, warped,
	   map);
  else
    printf("Loading section %s  (%p %p %p)\n", pairs[index].imageName,
	   image, ref, warped);

  if (image != 0)
    {
      free(image);
      image = 0;
    }
  if (ref != 0)
    {
      free(ref);
      ref = 0;
    }
  if (warped != 0)
    {
      free(warped);
      warped = 0;
    }

  /* read in the image */
  if (pairs[index].imageMinX >= 0)
    imageMinX = pairs[index].imageMinX / baseReductionFactor;
  else
    imageMinX = -1;
  if (pairs[index].imageMaxX >= 0)
    imageMaxX = pairs[index].imageMaxX / baseReductionFactor;
  else
    imageMaxX = -1;
  if (pairs[index].imageMinY >= 0)
    imageMinY = pairs[index].imageMinY / baseReductionFactor;
  else
    imageMinY = -1;
  if (pairs[index].imageMaxY >= 0)
    imageMaxY = pairs[index].imageMaxY / baseReductionFactor;
  else
    imageMaxY = -1;

  sprintf(fn, "%s%s.%s", inputName, pairs[index].imageName, extension);
#if 1
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

  count = 0;
  printf("Image limits: %d %d %d %d\n",
	 pairs[index].imageMinX, pairs[index].imageMaxX,
	 pairs[index].imageMinY, pairs[index].imageMaxY);
  for (y = 0; y < imageHeight; ++y)
    for (x = 0; x < imageWidth; ++x)
      if (image[y * imageWidth + x] > 0)
	++count;
  printf("Image is %f%% non-black.\n", ((double) count) * 100.0 / imageHeight / imageWidth);

#endif

  if (pairs[index].refName[0] == '\0')
    {
      sectionsLoaded = true;
      type = 0;
      return;
    }
    
  /* read in the reference image */
  if (pairs[index].refMinX >= 0)
    refMinX = pairs[index].refMinX / baseReductionFactor;
  else
    refMinX = -1;
  if (pairs[index].refMaxX >= 0)
    refMaxX = pairs[index].refMaxX / baseReductionFactor;
  else
    refMaxX = -1;
  if (pairs[index].refMinY >= 0)
    refMinY = pairs[index].refMinY / baseReductionFactor;
  else
    refMinY = -1;
  if (pairs[index].refMaxY >= 0)
    refMaxY = pairs[index].refMaxY / baseReductionFactor;
  else
    refMaxY = -1;

  sprintf(fn, "%s%s.%s", inputName, pairs[index].refName, extension);
  if (!ReadImage(fn, &ref, &refWidth, &refHeight,
		 refMinX, refMaxX,
		 refMinY, refMaxY,
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
  if (refMinX < 0)
    refMinX = 0;
  if (refMaxX < 0)
    refMaxX = refWidth - 1;
  if (refMinY < 0)
    refMinY = 0;
  if (refMaxY < 0)
    refMaxY = refHeight - 1;

  if (mapName[0] == '\0')
    {
      nCpts = 0;
      if (cpts != 0)
	{
	  free(cpts);
	  cpts = 0;
	}
      lastIndex = index;
      sectionsLoaded = true;
      return;
    }

  /* read in the map */
  sprintf(fn, "%s%s.map", mapName, pairs[index].pairName);
  if (!MyReadMap(fn, &map, &level, &mapWidth, &mapHeight,
	       &mapOffsetX, &mapOffsetY, imgn, refn, errorMsg))
    {
      if (partial)
	{
	  index += increment;
	  if (index < 0 || index >= nPairs)
	    index = lastIndex;
	  goto retry;
	}	  
      fprintf(stderr, "Error reading map %s:\n  %s\n", mapName, errorMsg);
      exit(1);
    }

  /* construct the inverse map */
  inverseMap = InvertMap(map, mapWidth, mapHeight);

  /* compute the warped image */
  mapFactor = (1 << level);
  printf("Set mapFactor to %d  (%d %d)\n", mapFactor, imageWidth, mapWidth);
  h = imageHeight;
  w = imageWidth;
  rh = refHeight;
  rw = refWidth;
  printf("image dim = (%d %d) imageMinX = %d  imageMinY = %d refMinX = %d refMinY = %d mapOffset = (%d %d)\n",
	 imageWidth, imageHeight, imageMinX, imageMinY, refMinX, refMinY,
	 mapOffsetX, mapOffsetY);

  int npp = 0;
  warped = (unsigned char *) malloc(imageWidth * imageHeight * sizeof(unsigned char));
  memset(warped, 0, imageWidth * imageHeight * sizeof(unsigned char));
  for (y = 0; y < h; ++y)
    for (x = 0; x < w; ++x)
      {
	//	fmDebug = (y == 1000 && x == 3995);
	if (!forwardmap(&rx, &ry, &rc,
			(float) (x + imageMinX + 0.5) * baseReductionFactor,
			(float) (y + imageMinY + 0.5) * baseReductionFactor))
	  {
	    if (fmDebug)
	      printf("FMDEBUG pt skipped\n");
	    continue;
	  }
	if (fmDebug)
	  printf("FMDEBUG pt %d %d %f %f\n",
		 x, y, rx, ry);
	if (x == 0 && y == 0)
	  {
	    printf("computed fwd map of (0, 0) to be %f %f\n",
		   rx, ry);
	  }
	rx = rx / baseReductionFactor - refMinX - 0.5;
	ry = ry / baseReductionFactor - refMinY - 0.5;
	irx = (int) floor(rx);
	iry = (int) floor(ry);
	rrx = rx - irx;
	rry = ry - iry;
	if (x == 0 && y == 0)
	  {
	    printf("irx,y (%d %d) rrx,y (%f %f)\n",
		   irx, iry, rrx, rry);
	  }
	r00 = RIMAGE(ref, irx, iry);
	r01 = RIMAGE(ref, irx, iry + 1);
	r10 = RIMAGE(ref, irx + 1, iry);
	r11 = RIMAGE(ref, irx + 1, iry + 1);
	rv = r00 * (rrx - 1.0) * (rry - 1.0)
	  - r10 * rrx * (rry - 1.0) 
	  - r01 * (rrx - 1.0) * rry
	  + r11 * rrx * rry;
	if (gray)
	  {
	    if (rc == 0.0)
	      b = 0;
	    else
	      b = (int) ((0.5 + 0.5 * rc) * rv);
	  }
	else
	  {
	    if (rc == 0.0)
	      b = 0;
	    else
	      b = (int) rv;
	  }
	if (b < 0)
	  b = 0;
	else if (b > 255)
	  b = 255;
	warped[y * w + x] = b;
      }
  //  printf("warped[70, 200] set to %d\n", warped[200*w+70]);

  // read in the correspondence points
  nCpts = 0;
  if (cpts != 0)
    {
      free(cpts);
      cpts = 0;
    }
  sprintf(fn, "%s%s.pts", cptsName, pairs[index].pairName);
  f = fopen(fn, "r");
  char line[LINE_LENGTH+1];
  float imx, imy, refx, refy;
  if (f != NULL)
    {
      //      printf("Opened cpts file %s\n", fn);
      while (fgets(line, LINE_LENGTH, f) != NULL)
	{
	  sscanf(line, "%f%f%f%f", &imx, &imy, &refx, &refy);
	  cpts = (CPoint*) realloc(cpts, (nCpts+1)*sizeof(CPoint));
	  cpts[nCpts].ix = imx;
	  cpts[nCpts].iy = imy;
	  cpts[nCpts].rx = refx;
	  cpts[nCpts].ry = refy;
	  cpts[nCpts].color = nCpts;
	  //	  printf("read in corres pt %d %d %d %d\n", imx, imy, refx, refy);
	  if (!backmap(&cpts[nCpts].wrx, &cpts[nCpts].wry, refx, refy))
	    {
	      printf("backmap failed\n");
	      cpts[nCpts].wrx = -1;
	      cpts[nCpts].wry = -1;
	    }
	  printf("LOADED CPTS[%d]: (%f %f) (%f %f) (%f %f)\n", nCpts,
		 cpts[nCpts].ix, cpts[nCpts].iy,
		 cpts[nCpts].rx, cpts[nCpts].ry,
		 cpts[nCpts].wrx, cpts[nCpts].wry);
	  ++nCpts;
	  //	  printf("Incremented nCpts to %d\n", nCpts);
	}
      fclose(f);
    }

  displayWidth = max(imageWidth, refWidth);
  displayHeight = max(imageHeight, refHeight);

  scale = ((float) availableWidth) / displayWidth;
  if (((float) availableHeight) / displayHeight < scale)
    scale = ((float) availableHeight) / displayHeight;
  reductionFactor = baseReductionFactor;
  while (scale <= 0.5)
    {
      // reduce the resolution by a factor of 2
      reductionFactor *= 2;

      newImageMinX = (imageMinX + 1) / 2;
      newImageMaxX = (imageMaxX - 1) / 2;
      newImageMinY = (imageMinY + 1) / 2;
      newImageMaxY = (imageMaxY - 1) / 2;

      nw = newImageMaxX - newImageMinX + 1;
      nh = newImageMaxY - newImageMinY + 1;
      
      deltaX = 2 * newImageMinX - imageMinX;
      deltaY = 2 * newImageMinY - imageMinY;
      delta = deltaY * imageWidth + deltaX;
      for (y = 0; y < nh; ++y)
	for (x = 0; x < nw; ++x)
	  {
	    iv = image[2*y*imageWidth + 2*x + delta];
	    iv += image[2*y*imageWidth + 2*x + 1 + delta];
	    iv += image[(2*y+1)*imageWidth + 2*x + delta];
	    iv += image[(2*y+1)*imageWidth + 2*x + 1 + delta];
	    image[y*nw+x] = (iv + 2) >> 2;

	    iv = warped[2*y*imageWidth + 2*x + delta];
	    iv += warped[2*y*imageWidth + 2*x + 1 + delta];
	    iv += warped[(2*y+1)*imageWidth + 2*x + delta];
	    iv += warped[(2*y+1)*imageWidth + 2*x + 1 + delta];
	    warped[y*nw+x] = (iv + 2) >> 2;
	  }
      imageWidth = nw;
      imageHeight = nh;
      imageMinX = newImageMinX;
      imageMaxX = newImageMaxX;
      imageMinY = newImageMinY;
      imageMaxY = newImageMaxY;
      image = (unsigned char *) realloc(image, imageHeight*imageWidth);
      warped = (unsigned char *) realloc(warped, imageHeight*imageWidth);

      newRefMinX = (refMinX + 1) / 2;
      newRefMaxX = (refMaxX - 1) / 2;
      newRefMinY = (refMinY + 1) / 2;
      newRefMaxY = (refMaxY - 1) / 2;

      nw = newRefMaxX - newRefMinX + 1;
      nh = newRefMaxY - newRefMinY + 1;
      deltaX = 2 * newRefMinX - refMinX;
      deltaY = 2 * newRefMinY - refMinY;
      delta = deltaY * refWidth + deltaX;
      for (y = 0; y < nh; ++y)
	for (x = 0; x < nw; ++x)
	  {
	    iv = ref[2*y*refWidth + 2*x + delta];
	    iv += ref[2*y*refWidth + 2*x + 1 + delta];
	    iv += ref[(2*y+1)*refWidth + 2*x + delta];
	    iv += ref[(2*y+1)*refWidth + 2*x + 1 + delta];
	    ref[y*nw+x] = (iv + 2) >> 2;
	  }
      refWidth = nw;
      refHeight = nh;
      refMinX = newRefMinX;
      refMaxX = newRefMaxX;
      refMinY = newRefMinY;
      refMaxY = newRefMaxY;
      ref = (unsigned char *) realloc(ref, refHeight*refWidth);

      displayWidth = max(imageWidth, refWidth);
      displayHeight = max(imageHeight, refHeight);

      scale = ((float) availableWidth) / displayWidth;
      if (((float) availableHeight) / displayHeight < scale)
	scale = ((float) availableHeight) / displayHeight;
    }
  offsetX = 0.0;
  offsetY = 0.0;
  windowWidth = (int) (scale * displayWidth);
  windowHeight = (int) (scale * displayHeight);

  lastIndex = index;
  sectionsLoaded = true;
  setTextures = true;
  //  printf("Leaving loadSection... nCpts = %d\n", nCpts);
}

void
Window::saveSection ()
{
  if (!readOnly)
    {
      char fn[PATH_MAX];
      sprintf(fn, "%s%s.pts", cptsName, pairs[index].pairName);
      bool fileExists = false;
      FILE *f = fopen(fn, "r");
      if (f != NULL)
	{
	  fileExists = true;
	  fclose(f);
	}
      bool first = true;
      for (int i = 0; i < nCpts; ++i)
	{
	  if (cpts[i].ix < 0 || cpts[i].rx < 0)
	    continue;
	  if (first)
	    {
	      first = false;
	      f = fopen(fn, "w");
	      if (f == NULL)
		{
		  fprintf(stderr, "Unable to write correspondence points file %s\n", fn);
		  Beep();
		  break;
		}
	    }
	  fprintf(f, "%f %f %f %f\n",
		  cpts[i].ix, cpts[i].iy,
		  cpts[i].rx, cpts[i].ry);
	}
      if (first)
	{
	  if (fileExists)
	    unlink(fn);
	}
      else
	{
	  if (f != NULL)
	    fclose(f);
	}
      printf("Saved correspondence points in %s\n", fn);
    }
  if (cpts != 0)
    {
      free(cpts);
      cpts = 0;
    }
  nCpts = 0;
}

void
Window::initTextures (int w, int h)
{
  // initialize the textures
  if (texture[0] > 0)
    glDeleteTextures(3, texture);
  glGenTextures(3, texture); 
  //  printf("textures = %d %d %d\n", texture[0], texture[1], texture[2]);
  textureWidth = w;
  textureHeight = h;
  unsigned char* black = (unsigned char *) malloc(3 * textureWidth * textureHeight);
  memset(black, 0, 3*textureWidth*textureHeight);
  glEnable(GL_TEXTURE_2D);
  for (int i = 0; i < 3; ++i)
    {
      glBindTexture(GL_TEXTURE_2D, texture[i]);
        
      glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_REPEAT);
      glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_REPEAT); 
      glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE); 
      //      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, textureWidth, textureHeight, 0,
      //      		   GL_RGB, GL_UNSIGNED_BYTE, black);
      glBindTexture(GL_TEXTURE_2D, 0);
    }
  glDisable(GL_TEXTURE_2D);
  free(black);
}

void
Window::initCircle ()
{
  int i;
  float x, y;

  circle = 1;
  glNewList(circle, GL_COMPILE);
  glBegin(GL_LINE_LOOP);
  for (i = 0; i < 100; ++i)
    {
      x = cos(i * 2 * M_PI / 100.0);
      y = sin(i * 2 * M_PI / 100.0);
      glVertex2f(x, y);
    }
  glEnd();
  glEndList();
}

void
Window::setTexture (unsigned int type, unsigned char *image, int w, int h)
{
  int x, y;
  unsigned char b;

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
  //  printf("Binding texture %d\n", type);
  glBindTexture(GL_TEXTURE_2D, texture[type]);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, textureWidth, textureHeight, 0,
	       GL_RGB, GL_UNSIGNED_BYTE, img);
  //  glTexSubImage2D(GL_TEXTURE_2D, 0,
  //  		  0, 0, textureWidth, textureHeight,
  //  		  GL_RGB, GL_UNSIGNED_BYTE, img);
  //  gluBuild2DMipmaps(GL_TEXTURE_2D, 3, textureWidth, textureHeight,
  //		   GL_RGB, GL_UNSIGNED_BYTE, img);
  glBindTexture(GL_TEXTURE_2D, 0);
  free(img);
}

bool
Window::backmap (float *wrx, float *wry, float x, float y)
{
  float rx, ry;
  if (!Invert(inverseMap, &rx, &ry, x / mapFactor, y / mapFactor))
    return(false);
  *wrx = (rx + mapOffsetX) * mapFactor;
  *wry = (ry + mapOffsetY) * mapFactor;
  return(true);
}

int
ReadHeader (FILE *f, char *tc, int *w, int *h, int *m)
{
  int c;
  int v;

  c = fgetc(f);
  while (c == '#')
    {
      while ((c = fgetc(f)) != EOF && c != '\n') ;
    }
  if (c != 'P')
    return(0);
  c = fgetc(f);
  if (c != '4' && c != '5' && c != '6')
    return(0);
  *tc = c;
  c = fgetc(f);

  while (c == ' ' || c == '\t' || c == '\n' || c == '#')
    if (c == '#')
      {
	while ((c = fgetc(f)) != EOF && c != '\n') ;
      }
    else
      c = fgetc(f);
  if (c >= '0' && c <= '9')
    {
      v = 0;
      while (c != EOF && c >= '0' && c <= '9')
	{
	  v = 10 * v + (c - '0');
	  c = fgetc(f);
	}
      *w = v;
    }
  else
    return(0);

  while (c == ' ' || c == '\t' || c == '\n' || c == '#')
    if (c == '#')
      {
	while ((c = fgetc(f)) != EOF && c != '\n') ;
      }
    else
      c = fgetc(f);
  if (c >= '0' && c <= '9')
    {
      v = 0;
      while (c != EOF && c >= '0' && c <= '9')
	{
	  v = 10 * v + (c - '0');
	  c = fgetc(f);
	}
      *h = v;
    }
  else
    return(0);

  if (*tc == '4')
    return(1);

  while (c == ' ' || c == '\t' || c == '\n' || c == '#')
    if (c == '#')
      while ((c = fgetc(f)) != EOF && c != '\n') ;
    else
      c = fgetc(f);
  if (c >= '0' && c <= '9')
    {
      v = 0;
      while (c != EOF && c >= '0' && c <= '9')
	{
	  v = 10 * v + (c - '0');
	  c = fgetc(f);
	}
      *m = v;
    }
  else
    return(0);

  return(1);
}

int
ParseRange (char *s, int *pos, int *minValue, int *maxValue)
{
  int i, j;
  int v;
  int minSpecified, maxSpecified;

  i = *pos;
  v = 0;
  while (s[i] != '\0' && s[i] >= '0' && s[i] <= '9')
    v = 10 * v + (s[i++] - '0');
  minSpecified = (i != *pos);
  *minValue = v;
  if (s[i] == '-')
    {
      j = ++i;
      v = 0;
      while (s[i] != '\0' && s[i] >= '0' && s[i] <= '9')
	v = 10 * v + (s[i++] - '0');
      maxSpecified = (i != j);
      if (maxSpecified)
	*maxValue = v;
      else
	*maxValue = 0x7fffffff;
    }
  else
    if (minSpecified)
      *maxValue = v;
    else
      return(0);
  *pos = i;
  return(1);
}

void
Beep ()
{
  printf("\007");
  fflush(stdout);
}

bool
Window::forwardmap (float *px, float *py, float *pc, float x, float y)
{
  float xv, yv;
  int ixv, iyv;
  int irx, iry;
  float rrx, rry;
  float r00, r01, r10, r11;
  float rv;
  float rx00, rx01, rx10, rx11;
  float ry00, ry01, ry10, ry11;
  float rc00, rc01, rc10, rc11;
  float rx, ry, rc;

  // use bilinear interpolation to find value
  xv = x / mapFactor - mapOffsetX;
  yv = y / mapFactor - mapOffsetY;
  ixv = ((int) (xv + 2.0)) - 2;
  iyv = ((int) (yv + 2.0)) - 2;
  rrx = xv - ixv;
  rry = yv - iyv;
  if (ixv < 0 || ixv >= mapWidth-1 ||
      iyv < 0 || iyv >= mapHeight-1)
    return(false);

  rx00 = MAP(map, ixv, iyv).x;
  ry00 = MAP(map, ixv, iyv).y;
  rc00 = MAP(map, ixv, iyv).c;
  rx01 = MAP(map, ixv, iyv + 1).x;
  ry01 = MAP(map, ixv, iyv + 1).y;
  rc01 = MAP(map, ixv, iyv + 1).c;
  rx10 = MAP(map, ixv + 1, iyv).x;
  ry10 = MAP(map, ixv + 1, iyv).y;
  rc10 = MAP(map, ixv + 1, iyv).c;
  rx11 = MAP(map, ixv + 1, iyv + 1).x;
  ry11 = MAP(map, ixv + 1, iyv + 1).y;
  rc11 = MAP(map, ixv + 1, iyv + 1).c;

  rx = rx00 * (rrx - 1.0) * (rry - 1.0)
    - rx10 * rrx * (rry - 1.0) 
    - rx01 * (rrx - 1.0) * rry
    + rx11 * rrx * rry;
  ry = ry00 * (rrx - 1.0) * (rry - 1.0)
    - ry10 * rrx * (rry - 1.0) 
    - ry01 * (rrx - 1.0) * rry
    + ry11 * rrx * rry;
  if (rc00 == 0.0 || rc01 == 0.0 ||
      rc10 == 0.0 || rc11 == 0.0)
    rc = 0.0;
  else
    rc = rc00 * (rrx - 1.0) * (rry - 1.0)
      - rc10 * rrx * (rry - 1.0) 
      - rc01 * (rrx - 1.0) * rry
      + rc11 * rrx * rry;
  *px = mapFactor * rx;
  *py = mapFactor * ry;
  *pc = rc;
  return(true);
}
