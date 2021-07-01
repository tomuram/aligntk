#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <dirent.h>

// genmask5 is for low res dataset in emmons_project

#include "par.h"

typedef struct Context {
  /* NOTE: any new fields added to this struct should
     be also added to PackContext and UnpackContext */ 
  int nDigits;
  char inputName[PATH_MAX];
  char outputName[PATH_MAX];
  char outputGrayscaleName[PATH_MAX];
} Context;

typedef struct Task {
  /* NOTE: any new fields added to this struct should
     be also added to PackTask and UnpackTask */ 
  int slice;
} Task;

typedef struct Result {
  char message[PATH_MAX + 1024];
} Result;

/* GLOBAL VARIABLES FOR MASTER */
int nProcessed = 0;

/* GLOBAL VARIABLES FOR MASTER & WORKER */
Context c;
Task t;
Result res;

/* GLOBAL VARIABLES FOR WORKER */
int imageWidth = 0;
int imageHeight = 0;
unsigned char *image = 0;
unsigned char *mask = 0;
unsigned char *grayscale = 0;
unsigned char *background = 0;

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
void FloodFill (int x, int y, int w, int h);
int Compare (const void *x, const void *y);
int ParseRange (char *s, int *pos, int *minValue, int *maxValue);
int ReadHeader (FILE *f, char *tc, int *w, int *h, int *m);
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
  int nSkipSlices;
  int maxSkipSlices;
  int *skipSlices;
  int len;
  char skipFile[PATH_MAX];
  int minSlice, maxSlice;
  int error;
  int pos;
  FILE *f;
  char inputDirName[PATH_MAX];
  char inputPrefix[PATH_MAX];
  int inputPrefixLen;
  char fn[PATH_MAX];
  int minN, maxN;
  int digitLen;
  DIR *dir;
  struct dirent *de;
  char tc;
  int w, h;
  int m;
  int skipIndex;
  int minS, maxS;

  error = 0;
  c.inputName[0] = '\0';
  c.outputName[0] = '\0';
  c.outputGrayscaleName[0] = '\0';
  skipFile[0] = '\0';
  minSlice = 0;
  maxSlice = 1000000000;
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
    else if (strcmp(argv[i], "-output_grayscale") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-output_grayscale error\n");
	    break;
	  }
	strcpy(c.outputGrayscaleName, argv[i]);
      }
    else if (strcmp(argv[i], "-skip") == 0)
      {
	if (++i == argc)
	  {
	    fprintf(stderr, "-skip error\n");
	    error = 1;
	    break;
	  }
	strcpy(skipFile, argv[i]);
      }
    else if (strcmp(argv[i], "-slice") == 0)
      {
	pos = 0;
	if (++i == argc ||
	    !ParseRange(argv[i], &pos, &minS, &maxS) ||
	    argv[i][pos] != '\0')
	  {
	    error = 1;
	    fprintf(stderr, "-slice error 1\n");
	    break;
	  }
	if (minS <= maxS)
	  {
	    minSlice = minS;
	    maxSlice = maxS;
	  }
	else
	  {
	    minSlice = maxS;
	    maxSlice = minS;
	  }
      }
  else
    {
      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      error = 1;
    }
      

  if (error)
    {
      fprintf(stderr, "Usage: genmask4 -input file_prefix -output file_prefix\n");
      fprintf(stderr, "              [-output_grayscale file_prefix]\n");
      fprintf(stderr, "              [-skip skip_file] [-slice slice_range]\n");
      fprintf(stderr, "   where slice_range is expressed as: integer\n");
      fprintf(stderr, "                              or: integer-integer\n");
      fprintf(stderr, "                              or: integer-\n");
      fprintf(stderr, "                              or: -integer\n");
      exit(1);
    }

  /* check that at least minimal parameters were supplied */
  if (c.inputName[0] == '\0' || c.outputName[0] == '\0')
    {
      fprintf(stderr, "Both -input and -output parameters must be specified.\n");
      exit(1);
    }

  /* get a sorted list of the slices to skip */
  nSkipSlices = 0;
  skipSlices = (int *) malloc(sizeof(int));
  maxSkipSlices = 1;
  if (skipFile[0] != '\0')
    {
      f = fopen(skipFile, "r");
      if (f == NULL)
	{
	  fprintf(stderr, "Could not open skip file %s\n", skipFile);
	  exit(1);
	}
      while (fscanf(f, "%d", &n) == 1)
	{
	  if (++nSkipSlices > maxSkipSlices)
	    {
	      maxSkipSlices *= 2;
	      skipSlices = realloc(skipSlices, maxSkipSlices * sizeof(int));
	    }
	  skipSlices[nSkipSlices-1] = n;
	}
      if (!feof(f))
	{
	  fprintf(stderr, "Read error while reading skip file %s\n", skipFile);
	  exit(1);
	}
      
      fclose(f);
      qsort(skipSlices, nSkipSlices, sizeof(int), Compare);
    }    

  /* check what images are present */
  for (i = strlen(c.inputName)-1; i >= 0 && c.inputName[i] != '/'; --i) ;
  if (i >= 0)
    {
      strncpy(inputDirName, c.inputName, i+1);
      inputDirName[i+1] = '\0';
      strcpy(inputPrefix, &c.inputName[i+1]);
    }
  else
    {
      strcpy(inputDirName, "./");
      strcpy(inputPrefix, c.inputName);
    }

  /* read the directory to look for files of the form nameNNNN.ppm */
  dir = opendir(inputDirName);
  if (dir == NULL)
    {
      fprintf(stderr, "Could not open directory %s for input files\n");
      exit(1);
    }
  inputPrefixLen = strlen(inputPrefix);
  minN = 1000000000;
  maxN = -1;
  c.nDigits = -1;
  while ((de = readdir(dir)) != NULL)
    {
      if (strncmp(de->d_name, inputPrefix, inputPrefixLen) == 0 &&
	  (len = strlen(de->d_name)) > 4 &&
	  strcmp(&(de->d_name[len-4]), ".ppm") == 0)
	{
	  if (sscanf(&(de->d_name[inputPrefixLen]), "%d%n",
		     &n, &digitLen) < 1)
	    continue;
	  if (c.nDigits < 0)
	    c.nDigits = digitLen;
	  else if (digitLen != c.nDigits)
	    {
	      fprintf(stderr, "Inconsistent number of digits in input file names: %s\n", de->d_name);
	      exit(1);
	    }
	  if (n < minN)
	    {
	      minN = n;
	      
	      /* read the header to determine image size */
	      sprintf(fn, "%s%s", inputDirName, de->d_name);
	      f = fopen(fn, "r");
	      if (f == NULL)
		{
		  fprintf(stderr, "Could not open image file %s", de->d_name);
		  exit(1);
		}
	      if (!ReadHeader(f, &tc, &w, &h, &m))
		{
		  fclose(f);
		  fprintf(stderr, "Invalid image file header: %s\n", de->d_name);
		  exit(1);
		}
	      if (tc != '6')
		{
		  fclose(f);
		  fprintf(stderr, "Image file not color ppm: %s\n", de->d_name);
		  exit(1);
		}
	      if (m != 255)
		{
		  fclose(f);
		  fprintf(stderr, "Image file %s does not have byte-size image components. (%d)\n", de->d_name, m);
		  exit(1);
		}
	      fclose(f);
	    }
	  if (n > maxN)
	    maxN = n;
	}
    }
  closedir(dir);
  if (maxN < 0)
    {
      fprintf(stderr, "No input image files found.\n");
      exit(1);
    }
  par_set_context();

  /* for all slices */
  printf("Processing slices: ");
  fflush(stdout);
  skipIndex = 0;
  for (t.slice = minN; t.slice <= maxN; ++t.slice)
    {
      if (t.slice < minSlice || t.slice > maxSlice)
	continue;
      while (skipIndex < nSkipSlices && skipSlices[skipIndex] < t.slice)
	++skipIndex;
      if (skipIndex < nSkipSlices && skipSlices[skipIndex] == t.slice)
	continue;
      par_delegate_task();
    }
  par_finish();

  printf(" %d\nAll slices completed.\n", nProcessed);
  exit(0);
}

void
MasterResult ()
{
  if (res.message[0] != '\0')
    {
      fprintf(stderr, "\nThe following error was encountered by one of the worker processes:\n%s\n", res.message);
      exit(1);
    }
  if ((nProcessed % 50) == 0 && nProcessed != 0)
    printf(" %d \n                   ", nProcessed);
  printf(".");
  ++nProcessed;
}

int
Compare (const void *x, const void *y)
{
  if (*((int *) x) < *((int *) y))
    return(-1);
  else if (*((int *) x) == *((int *) y))
    return(0);
  else
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


/* WORKER PROCEDURES */

void
WorkerContext ()
{
  /* nothing necessary here */
}

void
WorkerTask ()
{
  FILE *f;
  int w, h;
  char fn[PATH_MAX];
  int i, j;
  int n;
  char tc;
  int m;
  int x, y;
  int bpl;
  float cx, cy, cz;
  float min_radius_x, max_radius_x, min_radius_y, max_radius_y;
  float rz;
  float cutoff_radius_x, cutoff_radius_y;
  float r, g, b;
  float factor;

  /* get the image */
  sprintf(fn, "%s%0.*d.ppm", c.inputName, c.nDigits, t.slice);
  f = fopen(fn, "r");
  if (f == NULL)
    {
      sprintf(res.message, "Worker could not open image file %s\n", fn);
      return;
    }

  if (!ReadHeader(f, &tc, &w, &h, &m))
    {
      sprintf(res.message, "Image file not binary ppm: %s\n", fn);
      fclose(f);
      return;
    }
  if (tc != '6')
    {
      sprintf(res.message, "Image file not color ppm: %s\n", fn);
      fclose(f);
      return;
    }
  bpl = (w+7) >> 3;
  
  if (w != imageWidth || h != imageHeight)
    {
      image = (unsigned char *) realloc(image, 3*w*h);
      mask = (unsigned char *) realloc(mask, h*bpl);
      grayscale = (unsigned char *) realloc(grayscale, w*h);
      background = (unsigned char *) realloc(background, w*h);
      imageWidth = w;
      imageHeight = h;
    }

  if (fread(image, 1, 3*w*h, f) != 3*w*h)
    {
      fclose(f);
      sprintf(res.message, "Image file %s apparently truncated\n", fn);
      return;
    }
  fclose(f);

  /* go through all pixels */
  factor = 1.3;
  min_radius_x = 0.051 * factor;
  max_radius_x = 0.068 * factor;
  min_radius_y = 0.081 * factor;
  max_radius_y = 0.129 * factor;
  for (y = 0; y < h; ++y)
    for (x = 0; x < w; ++x)
      {
	r = image[3*(y*w + x)] / 255.0;
	g = image[3*(y*w + x) + 1] / 255.0;
	b = image[3*(y*w + x) + 2] / 255.0;

#if 0
	cx = -0.746*r + 0.662*g + 0.084*b;
	cy = 0.334*r + 0.480*g - 0.814*b;
	cz = -0.576*r - 0.577*g - 0.579*b;

	rz = - cz / (0.576 + 0.577 + 0.579);
	cutoff_radius_x = min_radius_x + rz * (max_radius_x - min_radius_x);
	cutoff_radius_y = min_radius_y + rz * (max_radius_y - min_radius_y);

	if (cx * cx / (cutoff_radius_x * cutoff_radius_x) +
	    cy * cy / (cutoff_radius_y * cutoff_radius_y) > 1.0)
	  background[y*w + x] = 1;
	else
	  background[y*w + x] = 0;
	if (x == 400 && y == 400)
	  {
	    printf("\n%d (%f %f %f) (%f %f %f) %f (%f %f) %d\n", t.slice, r, g, b, cx, cy, cz, rz,
		   cutoff_radius_x, cutoff_radius_y, (int) background[y*w+x]);
	  }
#endif

	grayscale[y*w + x] = (int) floor(255.0 * (r + g + b) / 3.0 + 0.5);
      }

  /* set the result bitmap to make it black */
  memset(mask, 0xff, h*bpl);

#if 0
  /* flood fill background in from all perimeter pixels */
  for (x = 0; x < w; ++x)
    {
      y = 0;
      if (background[y*w + x] &&
	  (mask[y*bpl + (x >> 3)] & (0x80 >> (x & 7))) != 0)
	FloodFill(x, y, w, h);
      y = h-1;
      if (background[y*w + x] &&
	  (mask[y*bpl + (x >> 3)] & (0x80 >> (x & 7))) != 0)
	FloodFill(x, y, w, h);
    }
  for (y = 1; y < h-1; ++y)
    {
      x = 0;
      if (background[y*w + x] &&
	  (mask[y*bpl + (x >> 3)] & (0x80 >> (x & 7))) != 0)
	FloodFill(x, y, w, h);
      x = w-1;
      if (background[y*w + x] &&
	  (mask[y*bpl + (x >> 3)] & (0x80 >> (x & 7))) != 0)
	FloodFill(x, y, w, h);
    }
#endif

  /* set all strongly color pixels to mask out writing */
  for (y = 0; y < h; ++y)
    for (x = 0; x < w; ++x)
      {
	r = image[3*(y*w + x)] / 255.0;
	g = image[3*(y*w + x) + 1] / 255.0;
	b = image[3*(y*w + x) + 2] / 255.0;

	cx = -0.746*r + 0.662*g + 0.084*b;
	cy = 0.334*r + 0.480*g - 0.814*b;
	cz = -0.576*r - 0.577*g - 0.579*b;

	rz = - cz / (0.576 + 0.577 + 0.579);
	cutoff_radius_x = min_radius_x + rz * (max_radius_x - min_radius_x);
	cutoff_radius_y = min_radius_y + rz * (max_radius_y - min_radius_y);

	if (cx * cx / (cutoff_radius_x * cutoff_radius_x) +
	    cy * cy / (cutoff_radius_y * cutoff_radius_y) > 1.0)
	  mask[y*bpl + (x >> 3)] &= ~(0x80 >> (x & 7));
      }

  /* write out the single image */
  sprintf(fn, "%s%0.*d.pbm", c.outputName, c.nDigits, t.slice);
  f = fopen(fn, "w");
  if (f == NULL)
    {
      sprintf(res.message, "Could not open output file %s\n", fn);
      return;
    }

  fprintf(f, "P4\n%d %d\n", w, h);
  if (fwrite(mask, 1, h*bpl, f) != h*bpl)
    {
      sprintf(res.message, "Could not write output file %s\n", fn);
      return;
    }
  fclose(f);

  /* write out the grayscale image */
  if (c.outputGrayscaleName[0] != '\0')
    {
      sprintf(fn, "%s%0.*d.pgm", c.outputGrayscaleName, c.nDigits, t.slice);
      f = fopen(fn, "w");
      if (f == NULL)
	{
	  sprintf(res.message, "Could not open output file %s\n", fn);
	  return;
	}

      fprintf(f, "P5\n%d %d\n255\n", w, h);
      if (fwrite(grayscale, 1, w*h, f) != w*h)
	{
	  sprintf(res.message, "Could not write grayscale file %s\n", fn);
	  return;
	}
      fclose(f);
    }

  res.message[0] = '\0';
}

void FloodFill (int x, int y, int w, int h)
{
  int y1;
  int left, right;
  int sp;
  unsigned int *stack;
  int bpl;

  bpl = (w+7) >> 3;
  sp = 0;
  stack = (unsigned int *) malloc(h * w * sizeof(int));
  stack[sp++] = y*w+x;
  while (sp > 0)
    {
      --sp;
      y = stack[sp] / w;
      x = stack[sp] - y * w;
      y1 = y;
      while (y1 >= 0 &&
	     background[y1*w + x] &&
	     mask[y1*bpl + (x >> 3)] & (0x80 >> (x & 7)) != 0)
	--y1;
      ++y1;
      left = 0;
      right = 0;
      while (y1 < h &&
	     background[y1*w + x] &&
	     mask[y1*bpl + (x >> 3)] & (0x80 >> (x & 7)) != 0)
        {
	  mask[y1*bpl + (x >> 3)] &= ~(0x80 >> (x & 7));
	  if (!left &&
	      x > 0 &&
	     background[y1*w + (x-1)] &&
	     (mask[y1*bpl + ((x-1) >> 3)] & (0x80 >> ((x-1) & 7))) != 0)
            {
	      stack[sp++] = y1 * w + (x-1);
	      left = 1;
            }
	  else if (left &&
		   x > 0 &&
		   (!background[y1*w + (x-1)] ||
		    (mask[y1*bpl + ((x-1) >> 3)] & (0x80 >> ((x-1) & 7))) == 0))
	    left = 0;

	  if (!right &&
	      x < w-1 &&
	      background[y1*w + (x+1)] &&
	      (mask[y1*bpl + ((x+1) >> 3)] & (0x80 >> ((x+1) & 7))) != 0)
            {
	      stack[sp++] = y1 * w + (x+1);
	      right = 1;
            }
	  else if (right &&
		   x < w-1 &&
		   (!background[y1*w + (x+1)] ||
		    (mask[y1*bpl + ((x+1) >> 3)] & (0x80 >> ((x+1) & 7))) == 0))
	    right = 0;
	  ++y1;
	}
    }
  free(stack);
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


void
PackContext ()
{
  int i;

  par_pkint(c.nDigits);
  par_pkstr(c.inputName);
  par_pkstr(c.outputName);
  par_pkstr(c.outputGrayscaleName);
}

void
UnpackContext ()
{
  int i;

  c.nDigits = par_upkint();
  par_upkstr(c.inputName);
  par_upkstr(c.outputName);
  par_upkstr(c.outputGrayscaleName);
}

void
PackTask ()
{
  par_pkint(t.slice);
}

void
UnpackTask ()
{
  t.slice = par_upkint();
}

void
PackResult ()
{
  par_pkstr(res.message);
}

void
UnpackResult ()
{
  par_upkstr(res.message);
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
