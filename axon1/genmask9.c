#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <dirent.h>

#include "par.h"
#include "imio.h"

#define MASK_BLACK_METHOD		0
#define MASK_WHITE_METHOD		1
#define FILL_FROM_BOUNDARY_METHOD	2
#define FILL_BLACK_FROM_BOUNDARY_METHOD	3
#define FILL_WHITE_FROM_BOUNDARY_METHOD	4
#define SELECT_NONBLACK_CENTER_METHOD	5

#define LINE_LENGTH	255

typedef struct Context {
  /* NOTE: any new fields added to this struct should
     be also added to PackContext and UnpackContext */ 
  char inputName[PATH_MAX];
  char outputName[PATH_MAX];
  int maskMethod;
  float centerWidthPercentage;
  float centerHeightPercentage;
} Context;

typedef struct Task {
  /* NOTE: any new fields added to this struct should
     be also added to PackTask and UnpackTask */ 
  char imageName[PATH_MAX];
} Task;

typedef struct Result {
  char message[PATH_MAX + 1024];
} Result;

/* GLOBAL VARIABLES FOR MASTER */
int nProcessed = 0;

/* GLOBAL VARIABLES FOR MASTER & WORKER */
Context c;
Task t;
Result r;

/* GLOBAL VARIABLES FOR WORKER */
int background;

/* FORWARD DECLARATIONS */
void MasterTask (int argc, char **argv, char **envp);
void MasterResult ();
void WorkerContext ();
void WorkerTask ();
void PackContext ();
void UnpackContext ();
void PackTask ();
void UnpackTask ();
void PackResult ();
void UnpackResult ();
int Black (int v);
int NotBlack (int v);
int Background (int v);
void FloodFill (unsigned char *image, unsigned char *mask,
		size_t w, size_t h, int x, int y, int (*predicate)());
void Error (char *fmt, ...);

int
main (int argc, char **argv, char **envp)
{
  par_process(argc, argv, envp,
              (void (*)()) MasterTask, MasterResult,
              WorkerContext, WorkerTask, NULL,
              PackContext, UnpackContext,
              PackTask, UnpackTask,
              PackResult, UnpackResult);
}

/* MASTER PROCEDURES */

void
MasterTask (int argc,
            char **argv,
            char **envp)
{
  int i;
  int n;
  int len;
  char imageListName[PATH_MAX];
  int error;
  FILE *f;
  char fn[PATH_MAX];
  char line[LINE_LENGTH+1];
  int nItems;

  error = 0;
  c.inputName[0] = '\0';
  c.outputName[0] = '\0';
  c.maskMethod = FILL_FROM_BOUNDARY_METHOD;
  c.centerWidthPercentage = 25.0;
  c.centerHeightPercentage = 25.0;
  imageListName[0] = '\0';
  for (i = 1; i < argc; ++i)
    if (strcmp(argv[i], "-input") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-input error\n");
	    break;
	  }
	strcpy(c.inputName, argv[i]);
      }
    else if (strcmp(argv[i], "-output") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-output error\n");
	    break;
	  }
	strcpy(c.outputName, argv[i]);
      }
    else if (strcmp(argv[i], "-image_list") == 0)
      {
	if (++i == argc)
	  {
	    fprintf(stderr, "-image_list error\n");
	    error = 1;
	    break;
	  }
	strcpy(imageListName, argv[i]);
      }
    else if (strcmp(argv[i], "-center_width") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%f", &c.centerWidthPercentage) != 1)
	  {
	    error = 1;
	    break;
	  }
      }	
    else if (strcmp(argv[i], "-center_height") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%f", &c.centerHeightPercentage) != 1)
	  {
	    error = 1;
	    break;
	  }
      }	
    else if (strcmp(argv[i], "-black") == 0)
      c.maskMethod = MASK_BLACK_METHOD;
    else if (strcmp(argv[i], "-white") == 0)
      c.maskMethod = MASK_WHITE_METHOD;
    else if (strcmp(argv[i], "-fill_from_boundary") == 0)
      c.maskMethod = FILL_FROM_BOUNDARY_METHOD;
    else if (strcmp(argv[i], "-fill_black_from_boundary") == 0)
      c.maskMethod = FILL_BLACK_FROM_BOUNDARY_METHOD;
    else if (strcmp(argv[i], "-fill_white_from_boundary") == 0)
      c.maskMethod = FILL_WHITE_FROM_BOUNDARY_METHOD;
    else if (strcmp(argv[i], "-select_nonblack_center") == 0)
      c.maskMethod = SELECT_NONBLACK_CENTER_METHOD;
    else
      {
	fprintf(stderr, "Unknown option: %s\n", argv[i]);
	error = 1;
      }
      

  if (error)
    {
      fprintf(stderr, "Usage: genmask9 -input file_prefix -output file_prefix\n");
      fprintf(stderr, "                -image_list image_list_file\n");
      fprintf(stderr, "                [-black]\n");
      fprintf(stderr, "                [-white]\n");
      fprintf(stderr, "                [-fill_from_boundary]\n");
      fprintf(stderr, "                [-fill_black_from_boundary]\n");
      fprintf(stderr, "                [-fill_white_from_boundary]\n");
      fprintf(stderr, "                [-select_nonblack_center]\n");
      fprintf(stderr, "                [-center_width percentage]\n");
      fprintf(stderr, "                [-center_height percentage]\n");
      exit(1);
    }

  /* check that at least minimal parameters were supplied */
  if (c.inputName[0] == '\0' || c.outputName[0] == '\0' ||
      imageListName[0] == '\0')
    {
      fprintf(stderr, "-input, -output, and -image_list parameters must be specified.\n");
      exit(1);
    }

  par_set_context();

  /* for all slices */
  printf("Processing images: ");
  fflush(stdout);

  /* read the images file */
  f = fopen(imageListName, "r");
  if (f == NULL)
    Error("Could not open file %s for reading\n", imageListName);
  while (fgets(line, LINE_LENGTH, f) != NULL)
    {
      len = strlen(line);
      if (len > 0 && line[len-1] == '\n')
	line[len-1] = '\0';
      if (line[0] == '\0' || line[0] == '#')
	continue;
      nItems = sscanf(line, "%s", t.imageName);
      if (nItems != 1)
	Error("Malformed line in %s:\n%s\n", imageListName, line);
      par_delegate_task();
    }
  fclose(f);
  par_finish();

  printf(" %d\nAll images completed.\n", nProcessed);
}

void
MasterResult ()
{
  if (r.message[0] != '\0')
    {
      fprintf(stderr, "\nThe following error was encountered by one of the worker processes:\n%s\n", r.message);
      exit(1);
    }
  if ((nProcessed % 50) == 0 && nProcessed != 0)
    printf(" %d \n                   ", nProcessed);
  printf(".");
  fflush(stdout);
  ++nProcessed;
}


/* WORKER PROCEDURES */

void
WorkerContext ()
{
  /* nothing necessary here */
}

void
WorkerTask ()
{
  size_t w, h;
  char fn[PATH_MAX];
  int i, j;
  int n;
  int b;
  int x, y;
  int dx, dy;
  size_t bpl;
  int side;
  int count[256];
  unsigned char *image;
  int imageWidth, imageHeight;
  unsigned char *mask;
  char msg[PATH_MAX+256];
  int hwx, hwy;

  /* get the image */
  sprintf(fn, "%s%s", c.inputName, t.imageName);
  if (!ReadImage(fn, &image, &imageWidth, &imageHeight, -1, -1, -1, -1, msg))
    Error("Could not read image file %s:\n%s\n", fn, msg);
  w = imageWidth;
  h = imageHeight;
  bpl = (w+7) >> 3;
  /* set the result bitmap to make it black */
  mask = (unsigned char *) malloc(h * bpl);
  memset(mask, 0xff, h*bpl);

  if (c.maskMethod == MASK_BLACK_METHOD || c.maskMethod == MASK_WHITE_METHOD)
    {
      background = (c.maskMethod == MASK_BLACK_METHOD) ? 0x00 : 0xff;
      for (y = 0; y < h; ++y)
	for (x = 0; x < w; ++x)
	  if (image[y*w+x] == background)
	    mask[y*bpl + (x >> 3)] &= ~(0x80 >> (x & 7));
    }
  else
    {
      if (c.maskMethod == SELECT_NONBLACK_CENTER_METHOD)
	{
	  // flood fill out from center
	  hwx = (int) (c.centerWidthPercentage * 0.01 * w * 0.5);
	        // half-width of center in x direction
	  hwy = (int) (c.centerHeightPercentage * 0.01 * h * 0.5);
	        // half-width of center in y direction
	  for (dy = -hwy; dy < hwy; ++dy)
	    for (dx = -hwx; dx < hwx; ++dx)
	      {
		x = w/2 + dx;
		y = h/2 + dy;
		if (x >= 0 && x < w &&
		    y >= 0 && y < h &&
		    image[y*w+x] != 0 &&
		    (mask[y*bpl + (x >> 3)] & (0x80 >> (x & 7))) != 0)
		  FloodFill(image, mask, w, h, x, y, NotBlack);
	      }
	  // invert mask and copy back to image
	  for (y = 0; y < h; ++y)
	    for (x = 0; x < w; ++x)
	      if ((mask[y*bpl + (x >> 3)] & (0x80 >> (x & 7))) != 0)
		image[y*w + x] = 0;
	      else
		image[y*w + x] = 255;
	  memset(mask, 0xff, h*bpl);
	}

      switch (c.maskMethod)
	{
	case FILL_FROM_BOUNDARY_METHOD:
	  /* make a histogram of the perimeter pixels */
	  memset(count, 0, 256*sizeof(int));
	  for (x = 0; x < imageWidth; ++x)
	    {
	      y = 0;
	      ++count[image[y*w + x]];
	      y = h-1;
	      ++count[image[y*w + x]];
	    }
	  for (y = 1; y < imageHeight-1; ++y)
	    {
	      x = 0;
	      ++count[image[y*w + x]];
	      x = w-1;
	      ++count[image[y*w + x]];
	    }
	  
	  /* find the most common value and consider that the background */
	  background = 0;
	  for (i = 1; i < 256; ++i)
	    if (count[i] > count[background])
	      background = i;
	  break;

	case FILL_BLACK_FROM_BOUNDARY_METHOD:
	  background = 0;
	  break;

	case FILL_WHITE_FROM_BOUNDARY_METHOD:
	  background = 0xff;
	  break;

	case SELECT_NONBLACK_CENTER_METHOD:
	  background = 0;
	  break;

	default:
          Error("maskMethod is out-of-range\n");
	  break;
	}

      /* flood fill in from all perimeter pixels with that value */
      for (x = 0; x < imageWidth; ++x)
	{
	  y = 0;
	  if (image[y*w + x] == background &&
	      (mask[y*bpl + (x >> 3)] & (0x80 >> (x & 7))) != 0)
	    FloodFill(image, mask, w, h, x, y, Background);
	  y = h-1;
	  if (image[y*w + x] == background &&
	      (mask[y*bpl + (x >> 3)] & (0x80 >> (x & 7))) != 0)
	    FloodFill(image, mask, w, h, x, y, Background);
	}
      for (y = 1; y < imageHeight-1; ++y)
	{
	  x = 0;
	  if (image[y*w + x] == background &&
	      (mask[y*bpl + (x >> 3)] & (0x80 >> (x & 7))) != 0)
	    FloodFill(image, mask, w, h, x, y, Background);
	  x = w-1;
	  if (image[y*w + x] == background &&
	      (mask[y*bpl + (x >> 3)] & (0x80 >> (x & 7))) != 0)
	    FloodFill(image, mask, w, h, x, y, Background);
	}
    }
  free(image);
      
  /* write out the mask */
  sprintf(fn, "%s%s.pbm", c.outputName, t.imageName);
  if (!WriteBitmap(fn, mask, w, h, UncompressedBitmap, msg))
    Error("Could not write mask %s:\n%s\n", fn, msg);
  free(mask);

  r.message[0] = '\0';
}

int
Black (int v)
{
  return(v == 0);
}

int
NotBlack (int v)
{
  return(v != 0);
}

int
Background (int v)
{
  return(v == background);
}

// the following is an adaptation of floodFillScanlineStack found on
//   "Lode's Computer Graphics Tutorial"
//   (http://lodev.org/cgtutor/floodfill.html)
void
FloodFill (unsigned char *image, unsigned char *mask,
	   size_t w, size_t h, int x, int y, int (*predicate)())
{
  int y1;
  int left, right;
  size_t sp;
  size_t *stack;
  size_t bpl;

  bpl = (w+7) >> 3;
  sp = 0;
  stack = (size_t *) malloc(h * w * sizeof(int));
  stack[sp++] = y*w+x;
  while (sp > 0)
    {
      --sp;
      y = stack[sp] / w;
      x = stack[sp] - y * w;
      y1 = y;
      while (y1 >= 0 &&
	     (*predicate)(image[y1*w + x]) &&
	     (mask[y1*bpl + (x >> 3)] & (0x80 >> (x & 7))) != 0)
	--y1;
      ++y1;
      left = 0;
      right = 0;
      while (y1 < h &&
	     (*predicate)(image[y1*w + x]) &&
	     (mask[y1*bpl + (x >> 3)] & (0x80 >> (x & 7))) != 0)
        {
	  mask[y1*bpl + (x >> 3)] &= ~(0x80 >> (x & 7));
	  if (!left &&
	      x > 0 &&
	      (*predicate)(image[y1*w + (x-1)]) &&
	     (mask[y1*bpl + ((x-1) >> 3)] & (0x80 >> ((x-1) & 7))) != 0)
            {
	      stack[sp++] = y1 * w + (x-1);
	      left = 1;
            }
	  else if (left &&
		   x > 0 &&
		   (!(*predicate)(image[y1*w + (x-1)]) ||
		    (mask[y1*bpl + ((x-1) >> 3)] & (0x80 >> ((x-1) & 7))) == 0))
	    left = 0;

	  if (!right &&
	      x < w-1 &&
	      (*predicate)(image[y1*w + (x+1)]) &&
	      (mask[y1*bpl + ((x+1) >> 3)] & (0x80 >> ((x+1) & 7))) != 0)
            {
	      stack[sp++] = y1 * w + (x+1);
	      right = 1;
            }
	  else if (right &&
		   x < w-1 &&
		   (!(*predicate)(image[y1*w + (x+1)]) ||
		    (mask[y1*bpl + ((x+1) >> 3)] & (0x80 >> ((x+1) & 7))) == 0))
	    right = 0;
	  ++y1;
	}
    }
  free(stack);
}

void
PackContext ()
{
  par_pkstr(c.inputName);
  par_pkstr(c.outputName);
  par_pkint(c.maskMethod);
  par_pkfloat(c.centerWidthPercentage);
  par_pkfloat(c.centerHeightPercentage);
}

void
UnpackContext ()
{
  par_upkstr(c.inputName);
  par_upkstr(c.outputName);
  c.maskMethod = par_upkint();
  c.centerWidthPercentage = par_upkfloat();
  c.centerHeightPercentage = par_upkfloat();
}

void
PackTask ()
{
  par_pkstr(t.imageName);
}

void
UnpackTask ()
{
  par_upkstr(t.imageName);
}

void
PackResult ()
{
  par_pkstr(r.message);
}

void
UnpackResult ()
{
  par_upkstr(r.message);
}

void Error (char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
  fflush(stderr);
  abort();
}
