#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <stdarg.h>

#include "par.h"
#include "imio.h"

#define LINE_LENGTH	255

typedef struct Context {
  /* NOTE: any new fields added to this struct should
     be also added to PackContext and UnpackContext */ 
  char inputName[PATH_MAX];
  char outputName[PATH_MAX];
  int factor;
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
  int error;
  FILE *f;
  char inputDirName[PATH_MAX];
  char inputPrefix[PATH_MAX];
  int inputPrefixLen;
  char fn[PATH_MAX];
  int w, h;
  char imageListName[PATH_MAX];
  char line[LINE_LENGTH+1];
  int nItems;

  error = 0;
  c.factor = 2;
  c.inputName[0] = '\0';
  c.outputName[0] = '\0';
  imageListName[0] = '\0';
  for (i = 1; i < argc; ++i)
    if (strcmp(argv[i], "-factor") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d", &c.factor) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-input") == 0)
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
    else
      {
	fprintf(stderr, "Unknown option: %s\n", argv[i]);
	error = 1;
      }
      

  if (error)
    {
      fprintf(stderr, "Usage: reduce -input file_prefix -output file_prefix\n");
      fprintf(stderr, "              -image_list image_list_file\n");
      fprintf(stderr, "              [-factor reduction_factor]\n");
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
  printf("Processing slices: ");
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

  printf(" %d\nAll slices completed.\n", nProcessed);
  exit(0);
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
  size_t w, ow;
  char fn[PATH_MAX];
  int x, y;
  int dx, dy;
  long a;
  long v;
  int offset;
  int ix, iy;
  int imageWidth, imageHeight;
  int outputWidth, outputHeight;
  unsigned char *image;
  unsigned char *output;
  char msg[PATH_MAX+256];


  /* get the image */
  sprintf(fn, "%s%s", c.inputName, t.imageName);
  if (!ReadImage(fn, &image, &imageWidth, &imageHeight, -1, -1, -1, -1, msg))
    Error("Could not read image file %s:\n%s\n", fn, msg);
  outputWidth = imageWidth / c.factor;
  outputHeight = imageHeight / c.factor;
  w = imageWidth;
  ow = outputWidth;

  output = (unsigned char *) malloc(outputHeight * ow);
  a = c.factor * c.factor;
  offset = a / 2;
  for (y = 0; y < outputHeight; ++y)
    {
      iy = y * c.factor;
      for (x = 0; x < outputWidth; ++x)
	{
	  ix = x * c.factor;
	  v = 0;
	  for (dy = 0; dy < c.factor; ++dy)
	    for (dx = 0; dx < c.factor; ++dx)
	      v += image[(iy + dy)*w + ix + dx];
	  output[y * ow + x] = (v + offset) / a;
	}	
    }
  free(image);

  /* write out the single image */
  sprintf(fn, "%s%s.tif", c.outputName, t.imageName);
  if (!WriteImage(fn, output, outputWidth, outputHeight,
		  UncompressedImage, msg))
    Error("Could not write output image %s:\n%s\n", fn, msg);
  free(output);

  r.message[0] = '\0';
}

void
PackContext ()
{
  par_pkstr(c.inputName);
  par_pkstr(c.outputName);
  par_pkint(c.factor);
}

void
UnpackContext ()
{
  par_upkstr(c.inputName);
  par_upkstr(c.outputName);
  c.factor = par_upkint();
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
