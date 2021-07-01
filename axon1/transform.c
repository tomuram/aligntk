#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <dirent.h>

#include "par.h"

#define SHIFT   0
#define CROP    1
#define ROTATE  2
#define SCALE   3
#define STRETCH 4
#define SIZE    5
#define RESIZE  6

typedef struct Op
{
  int type;
  double p0;
  double p1;
  double p2;
  double p3;
} Op;

typedef struct Constraint
{
  /* a*x + b*y + c >= 0 */
  double a;
  double b;
  double c;
} Constraint;

typedef struct Context {
  /* NOTE: any new fields added to this struct should
     be also added to PackContext and UnpackContext */ 
  int nDigits;
  char inputName[PATH_MAX];
  int maskDigits;
  char maskName[PATH_MAX];
  char outputName[PATH_MAX];
  char outputMaskName[PATH_MAX];
  int oversamplingFactor;
  int nOps;
  Op* ops;
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
Task task;
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
void ApplyTransform(double m[2][3], double t[2][3], int n, Constraint *con);
int Compare (const void *x, const void *y);
int ParseRange (char *s, int *pos, int *minValue, int *maxValue);
int ReadHeader (FILE *f, char *tc, int *w, int *h, int *m);

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
  char maskDirName[PATH_MAX];
  char maskPrefix[PATH_MAX];
  int maskPrefixLen;
  char fn[PATH_MAX];
  int maskWidth, maskHeight;
  int minZ, maxZ;
  int digitLen;
  DIR *dir;
  struct dirent *de;
  char tc;
  int w, h;
  int m;
  int skipIndex;
  int minS, maxS;
  double minX, maxX, minY, maxY;
  double xShift, yShift;
  double theta;
  double width, height;
  double scale;
  double xScale, yScale;

  error = 0;
  c.nDigits = 0;
  c.inputName[0] = '\0';
  c.maskDigits = 0;
  c.maskName[0] = '\0';
  c.outputName[0] = '\0';
  c.outputMaskName[0] = '\0';
  c.oversamplingFactor = 0;
  c.nOps = 0;
  c.ops = NULL;
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
    else if (strcmp(argv[i], "-mask") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-mask error\n");
	    break;
	  }
	strcpy(c.maskName, argv[i]);
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
    else if (strcmp(argv[i], "-output_mask") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-output_mask error\n");
	    break;
	  }
	strcpy(c.outputMaskName, argv[i]);
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
    else if (strcmp(argv[i], "-shift") == 0)
      {
	if (i+2 >= argc ||
	    sscanf(argv[i+1], "%lf", &xShift) != 1 ||
	    sscanf(argv[i+2], "%lf", &yShift) != 1)
	  {
	    error = 1;
	    fprintf(stderr, "-shift <x_shift> <y_shift> error\n");
	    break;
	  }
	++c.nOps;
	c.ops = (Op*) realloc(c.ops, c.nOps*sizeof(Op));
	c.ops[c.nOps-1].type = SHIFT;
	c.ops[c.nOps-1].p0 = xShift;
	c.ops[c.nOps-1].p1 = yShift;
	c.ops[c.nOps-1].p2 = 0.0;
	c.ops[c.nOps-1].p3 = 0.0;
	i += 2;
      }
    else if (strcmp(argv[i], "-crop") == 0)
      {
	if (i+4 >= argc ||
	    sscanf(argv[i+1], "%lf", &minX) != 1 ||
	    sscanf(argv[i+2], "%lf", &maxX) != 1 ||
	    sscanf(argv[i+3], "%lf", &minY) != 1 ||
	    sscanf(argv[i+4], "%lf", &maxY) != 1 ||
	    minX > maxX || minY > maxY)
	  {
	    error = 1;
	    fprintf(stderr, "-crop <x_min> <x_max> <y_min> <y_max> error\n");
	    break;
	  }
	++c.nOps;
	c.ops = (Op*) realloc(c.ops, c.nOps*sizeof(Op));
	c.ops[c.nOps-1].type = CROP;
	c.ops[c.nOps-1].p0 = minX;
	c.ops[c.nOps-1].p1 = maxX;
	c.ops[c.nOps-1].p2 = minY;
	c.ops[c.nOps-1].p3 = maxY;
	i += 4;
      }
    else if (strcmp(argv[i], "-rotate") == 0)
      {
	if (i+1 >= argc ||
	    sscanf(argv[i+1], "%lf", &theta) != 1)
	  {
	    error = 1;
	    fprintf(stderr, "-rotate <theta_degrees_ccw> error\n");
	    break;
	  }
	++c.nOps;
	c.ops = (Op*) realloc(c.ops, c.nOps*sizeof(Op));
	c.ops[c.nOps-1].type = ROTATE;
	c.ops[c.nOps-1].p0 = theta;
	c.ops[c.nOps-1].p1 = 0.0;
	c.ops[c.nOps-1].p2 = 0.0;
	c.ops[c.nOps-1].p3 = 0.0;
	i += 1;
      }
    else if (strcmp(argv[i], "-scale") == 0)
      {
	if (i+1 > argc ||
	    sscanf(argv[i+1], "%lf", &scale) != 1)
	  {
	    error = 1;
	    fprintf(stderr, "-scale <scale> error\n");
	    break;
	  }
	++c.nOps;
	c.ops = (Op*) realloc(c.ops, c.nOps*sizeof(Op));
	c.ops[c.nOps-1].type = SCALE;
	c.ops[c.nOps-1].p0 = scale;
	c.ops[c.nOps-1].p1 = 0.0;
	c.ops[c.nOps-1].p2 = 0.0;
	c.ops[c.nOps-1].p3 = 0.0;
	i += 1;
      }
    else if (strcmp(argv[i], "-stretch") == 0)
      {
	if (i+2 > argc ||
	    sscanf(argv[i+1], "%lf", &xScale) != 1 ||
	    sscanf(argv[i+2], "%lf", &yScale) != 1)
	  {
	    error = 1;
	    fprintf(stderr, "-stretch <x_scale> <y_scale> error\n");
	    break;
	  }
	++c.nOps;
	c.ops = (Op*) realloc(c.ops, c.nOps*sizeof(Op));
	c.ops[c.nOps-1].type = STRETCH;
	c.ops[c.nOps-1].p0 = xScale;
	c.ops[c.nOps-1].p1 = yScale;
	c.ops[c.nOps-1].p2 = 0.0;
	c.ops[c.nOps-1].p3 = 0.0;
	i += 2;
      }
    else if (strcmp(argv[i], "-size") == 0)
      {
	if (i+2 > argc ||
	    sscanf(argv[i+1], "%lf", &width) != 1 ||
	    sscanf(argv[i+2], "%lf", &height) != 1 ||
	    width <= 0 || height <= 0)
	  {
	    error = 1;
	    fprintf(stderr, "-size <width> <height> error\n");
	    break;
	  }
	++c.nOps;
	c.ops = (Op*) realloc(c.ops, c.nOps*sizeof(Op));
	c.ops[c.nOps-1].type = SIZE;
	c.ops[c.nOps-1].p0 = width;
	c.ops[c.nOps-1].p1 = height;
	c.ops[c.nOps-1].p2 = 0.0;
	c.ops[c.nOps-1].p3 = 0.0;
	i += 2;
      }
    else if (strcmp(argv[i], "-resize") == 0)
      {
	if (i+2 > argc ||
	    sscanf(argv[i+1], "%lf", &width) != 1 ||
	    sscanf(argv[i+2], "%lf", &height) != 1 ||
	    width <= 0 || height <= 0)
	  {
	    error = 1;
	    fprintf(stderr, "-resize <width> <height> error\n");
	    break;
	  }
	++c.nOps;
	c.ops = (Op*) realloc(c.ops, c.nOps*sizeof(Op));
	c.ops[c.nOps-1].type = RESIZE;
	c.ops[c.nOps-1].p0 = width;
	c.ops[c.nOps-1].p1 = height;
	c.ops[c.nOps-1].p2 = 0.0;
	c.ops[c.nOps-1].p3 = 0.0;
	i += 2;
      }
    else
      {
	fprintf(stderr, "Unknown option: %s\n", argv[i]);
	error = 1;
      }
      

  if (error)
    {
      fprintf(stderr, "Usage: transform -input file_prefix [-mask file_prefix]\n");
      fprintf(stderr, "              -output file_prefix [-output_mask file_prefix]\n");
      fprintf(stderr, "              [-skip skip_file] [-slice slice_range]\n");
      fprintf(stderr, "              [-shift shift_x shift_y]\n");
      fprintf(stderr, "              [-crop min_x max_x min_y max_y]\n");
      fprintf(stderr, "              [-rotate theta_degrees_ccw]\n");
      fprintf(stderr, "              [-scale scale]\n");
      fprintf(stderr, "              [-stretch scale_x scale_y]\n");
      fprintf(stderr, "              [-size width height]");
      fprintf(stderr, "              [-resize width height]");
      fprintf(stderr, "   where the slice range is expressed as: integer\n");
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

  /* read the directory to look for files of the form nameNNNN.ppm  (or .pgm) */
  dir = opendir(inputDirName);
  if (dir == NULL)
    {
      fprintf(stderr, "Could not open directory %s for input files\n");
      exit(1);
    }
  inputPrefixLen = strlen(inputPrefix);
  minZ = 1000000000;
  maxZ = -1;
  c.nDigits = -1;
  while ((de = readdir(dir)) != NULL)
    {
      if (strncmp(de->d_name, inputPrefix, inputPrefixLen) == 0 &&
	  (len = strlen(de->d_name)) > 4 &&
	  strcmp(&(de->d_name[len-4]), ".pgm") == 0)
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
	  if (n < minZ)
	    {
	      minZ = n;
	      
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
	      if (m != 255)
		{
		  fclose(f);
		  fprintf(stderr, "Image file %s does not have byte-size image components. (%d)\n", de->d_name, m);
		  exit(1);
		}
	      fclose(f);
	    }
	  if (n > maxZ)
	    maxZ = n;
	}
    }
  closedir(dir);
  if (maxZ < 0)
    {
      fprintf(stderr, "No input image files found.\n");
      exit(1);
    }

  if (c.maskName[0] != '\0')
    {
      /* check what masks are present */
      for (i = strlen(c.maskName)-1; i >= 0 && c.maskName[i] != '/'; --i) ;
      if (i >= 0)
	{
	  strncpy(maskDirName, c.maskName, i+1);
	  maskDirName[i+1] = '\0';
	  strcpy(maskPrefix, &c.maskName[i+1]);
	}
      else
	{
	  strcpy(maskDirName, "./");
	  strcpy(maskPrefix, c.maskName);
	}
      
      /* read the directory to look for files of the form nameNNNN.pbm */
      dir = opendir(maskDirName);
      if (dir == NULL)
	{
	  fprintf(stderr, "Could not open directory %s for mask files\n");
	  exit(1);
	}
      maskPrefixLen = strlen(maskPrefix);
      c.maskDigits = -1;
      while ((de = readdir(dir)) != NULL)
	{
	  if (strncmp(de->d_name, maskPrefix, maskPrefixLen) == 0 &&
	      (len = strlen(de->d_name)) > 4 &&
	      strcmp(&(de->d_name[len-4]), ".pbm") == 0)
	    {
	      if (sscanf(&(de->d_name[maskPrefixLen]), "%d%n",
			 &n, &digitLen) < 1)
		continue;
	      if (c.maskDigits < 0)
		{
		  c.maskDigits = digitLen;

		  /* read the header to determine image size */
		  sprintf(fn, "%s%s", maskDirName, de->d_name);
		  f = fopen(fn, "r");
		  if (f == NULL)
		    {
		      fprintf(stderr, "Could not open mask file %s", fn);
		      exit(1);
		    }
		  if (!ReadHeader(f, &tc, &maskWidth, &maskHeight, &m))
		    {
		      fclose(f);
		      fprintf(stderr, "Invalid image file header: %s\n", fn);
		      exit(1);
		    }
		  if (tc != '4')
		    {
		      fclose(f);
		      fprintf(stderr, "Bad .pbm header in mask file: %s\n", fn);
		      exit(1);
		    }
		  if (maskWidth != w || maskHeight != h)
		    {
		      fclose(f);
		      fprintf(stderr, "Mask size does not match image size: %s\n",
			      fn);
		      exit(1);
		    }
		  fclose(f);
		}
	      else if (digitLen != c.maskDigits)
		{
		  fprintf(stderr, "Inconsistent number of digits in mask file names: %s\n", de->d_name);
		  exit(1);
		}
	    }
	}
      closedir(dir);
    }

  par_set_context();

  /* for all slices */
  printf("Processing slices: ");
  fflush(stdout);
  skipIndex = 0;
  for (task.slice = minZ; task.slice <= maxZ; ++task.slice)
    {
      if (task.slice < minSlice || task.slice > maxSlice)
	continue;
      while (skipIndex < nSkipSlices && skipSlices[skipIndex] < task.slice)
	++skipIndex;
      if (skipIndex < nSkipSlices && skipSlices[skipIndex] == task.slice)
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
  if (r.message[0] != '\0')
    {
      fprintf(stderr, "\nThe following error was encountered by one of the worker processes:\n%s\n", r.message);
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
  int v;
  int x, y;
  int dx, dy;
  unsigned char *image;
  unsigned char *mask;
  unsigned char *output;
  unsigned char *outputMask;
  int nConstraints;
  Constraint* constraints;
  int oversamplingFactor;
  int of2;
  int mh, mw;
  double t[2][3];
  double m[2][3];
  int xSize, ySize;
  int maskPresent;
  int mbpl, ombpl;
  int outputWidth, outputHeight;
  int valid;
  int nx, ny;
  int sum;
  double px, py;
  int outside;

  /* get the image */
  sprintf(fn, "%s%0.*d.pgm", c.inputName, c.nDigits, task.slice);
  f = fopen(fn, "r");
  if (f == NULL)
    {
      sprintf(r.message, "Worker could not open image file %s\n", fn);
      return;
    }

  if (!ReadHeader(f, &tc, &w, &h, &v))
    {
      sprintf(r.message, "Image file not binary pgm or ppm: %s\n", fn);
      fclose(f);
      return;
    }
  if (tc != '5')
    {
      sprintf(r.message, "Image file not grayscale pgm: %s\n", fn);
      fclose(f);
      return;
    }

  image = (unsigned char *) malloc(w*h);
  if (fread(image, 1, w*h, f) != w*h)
    {
      fclose(f);
      sprintf(r.message, "Image file %s apparently truncated\n", fn);
      return;
    }
  fclose(f);
  nConstraints = 4;
  constraints = (Constraint *) malloc(nConstraints * sizeof(Constraint));
  constraints[0].a = 1.0;
  constraints[0].b = 0.0;
  constraints[0].c = 0.0;
  constraints[1].a = -1.0;
  constraints[1].b = 0.0;
  constraints[1].c = w;
  constraints[2].a = 0.0;
  constraints[2].b = 1.0;
  constraints[2].c = 0.0;
  constraints[3].a = 0.0;
  constraints[3].b = -1.0;
  constraints[3].c = h;
  t[0][0] = 1.0;
  t[0][1] = 0.0;
  t[0][2] = 0.0;
  t[1][0] = 0.0;
  t[1][1] = 1.0;
  t[1][2] = 0.0;
  xSize = w;
  ySize = h;

  /* read in the mask, if one exists */
  maskPresent = 0;
  if (c.maskName[0] != '\0')
    {
      sprintf(fn, "%s%0.*d.pbm", c.maskName, c.maskDigits, task.slice);
      f = fopen(fn, "r");
      if (f != NULL)
	{
	  if (!ReadHeader(f, &tc, &mw, &mh, &v) || tc != '4')
	    {
	      sprintf(r.message, "Mask file not binary pbm: %s\n", fn);
	      fclose(f);
	      return;
	    }
	  if (mw != w || mh != h)
	    {
	      sprintf(r.message, "Inconsistent mask sizes: %s\n", fn);
	      fclose(f);
	      return;
	    }
	  mbpl = (w + 7) >> 3;
	  mask = (unsigned char *) malloc(h*mbpl);
	  if (fread(mask, 1, h*mbpl, f) != h*mbpl)
	    {
	      fclose(f);
	      sprintf(r.message, "Mask file %s apparently truncated\n", fn);
	      return;
	    }
	  maskPresent = 1;
	  fclose(f);
	}
    }
  if (!maskPresent)
    {
      mbpl = (w + 7) >> 3;
      mask = (unsigned char *) malloc(h * mbpl);
      memset(mask, 0xff, h*mbpl);
    }

  oversamplingFactor = 1;
  for (i = 0; i < c.nOps; ++i)
    switch(c.ops[i].type)
      {
      case SHIFT:
	if (c.ops[i].p0 != floor(c.ops[i].p0) ||
	    c.ops[i].p1 != floor(c.ops[i].p1))
	  oversamplingFactor = 16;
	m[0][0] = 1.0;
	m[0][1] = 0.0;
	m[0][2] = c.ops[i].p0;
	m[1][0] = 0.0;
	m[1][1] = 1.0;
	m[1][2] = c.ops[i].p1;
	ApplyTransform(m, t, nConstraints, constraints);
	break;

      case CROP:
	if (c.ops[i].p0 != floor(c.ops[i].p0) ||
	    c.ops[i].p1 != floor(c.ops[i].p1) ||
	    c.ops[i].p2 != floor(c.ops[i].p2) ||
	    c.ops[i].p3 != floor(c.ops[i].p3))
	  oversamplingFactor = 16;
	nConstraints += 4;
	constraints = (Constraint *) realloc(constraints,
					     nConstraints * sizeof(Constraint));
	constraints[nConstraints-4].a = 1.0;
	constraints[nConstraints-4].b = 0.0;
	constraints[nConstraints-4].c = -c.ops[i].p0;
	constraints[nConstraints-3].a = -1.0;
	constraints[nConstraints-3].b = 0.0;
	constraints[nConstraints-3].c = c.ops[i].p1;
	constraints[nConstraints-2].a = 0.0;
	constraints[nConstraints-2].b = 1.0;
	constraints[nConstraints-2].c = -c.ops[i].p2;
	constraints[nConstraints-1].a = 0.0;
	constraints[nConstraints-1].b = -1.0;
	constraints[nConstraints-1].c = c.ops[i].p3;
	break;

      case ROTATE:
	if (fmod(c.ops[i].p0, 90.0) != 0.0)
	  oversamplingFactor = 16;
	m[0][0] = cos(c.ops[i].p0 * M_PI / 180.0);
	m[0][1] = sin(c.ops[i].p0 * M_PI / 180.0);
	m[0][2] = 0.0;
	m[1][0] = -sin(c.ops[i].p0 * M_PI / 180.0);
	m[1][1] = cos(c.ops[i].p0 * M_PI / 180.0);
	m[1][2] = 0.0;
	ApplyTransform(m, t, nConstraints, constraints);
	break;

      case SCALE:
	if (c.ops[i].p0 < 1.0 || c.ops[i].p0 != floor(c.ops[i].p0))
	  oversamplingFactor = 16;
	m[0][0] = c.ops[i].p0;
	m[0][1] = 0.0;
	m[0][2] = 0.0;
	m[1][0] = 0.0;
	m[1][1] = c.ops[i].p0;
	m[1][2] = 0.0;
	ApplyTransform(m, t, nConstraints, constraints);
	break;

      case STRETCH:
	if (c.ops[i].p0 < 1.0 || c.ops[i].p1 < 1.0 ||
	    c.ops[i].p0 != floor(c.ops[i].p0) || c.ops[i].p1 != floor(c.ops[i].p1))
	  oversamplingFactor = 16;
	m[0][0] = c.ops[i].p0;
	m[0][1] = 0.0;
	m[0][2] = 0.0;
	m[1][0] = 0.0;
	m[1][1] = c.ops[i].p1;
	m[1][2] = 0.0;
	ApplyTransform(m, t, nConstraints, constraints);
	break;

      case SIZE:
	if (c.ops[i].p0 != floor(c.ops[i].p0) ||
	    c.ops[i].p1 != floor(c.ops[i].p1))
	  oversamplingFactor = 16;
	nConstraints += 4;
	constraints = (Constraint *) realloc(constraints,
					     nConstraints * sizeof(Constraint));
	constraints[nConstraints-4].a = 1.0;
	constraints[nConstraints-4].b = 0.0;
	constraints[nConstraints-4].c = 0.0;
	constraints[nConstraints-3].a = -1.0;
	constraints[nConstraints-3].b = 0.0;
	constraints[nConstraints-3].c = c.ops[i].p0;
	constraints[nConstraints-2].a = 0.0;
	constraints[nConstraints-2].b = 1.0;
	constraints[nConstraints-2].c = 0.0;
	constraints[nConstraints-1].a = 0.0;
	constraints[nConstraints-1].b = -1.0;
	constraints[nConstraints-1].c = c.ops[i].p1;
	xSize = floor(c.ops[i].p0);
	ySize = floor(c.ops[i].p1);
	break;

      case RESIZE:
	if (c.ops[i].p0 != floor(c.ops[i].p0) ||
	    c.ops[i].p1 != floor(c.ops[i].p1) ||
	    c.ops[i].p0 < xSize ||
	    c.ops[i].p1 < ySize ||
	    fmod(c.ops[i].p0, (double) xSize) != 0.0 ||
	    fmod(c.ops[i].p1, (double) ySize) != 0.0)
	  oversamplingFactor = 16;
	m[0][0] = c.ops[i].p0 / xSize;
	m[0][1] = 0.0;
	m[0][2] = 0.0;
	m[1][0] = 0.0;
	m[1][1] = c.ops[i].p1 / ySize;
	m[1][2] = 0.0;
	ApplyTransform(m, t, nConstraints, constraints);
	nConstraints += 4;
	constraints = (Constraint *) realloc(constraints,
					     nConstraints * sizeof(Constraint));
	constraints[nConstraints-4].a = 1.0;
	constraints[nConstraints-4].b = 0.0;
	constraints[nConstraints-4].c = 0.0;
	constraints[nConstraints-3].a = -1.0;
	constraints[nConstraints-3].b = 0.0;
	constraints[nConstraints-3].c = c.ops[i].p0;
	constraints[nConstraints-2].a = 0.0;
	constraints[nConstraints-2].b = 1.0;
	constraints[nConstraints-2].c = 0.0;
	constraints[nConstraints-1].a = 0.0;
	constraints[nConstraints-1].b = -1.0;
	constraints[nConstraints-1].c = c.ops[i].p1;
	xSize = floor(c.ops[i].p0);
	ySize = floor(c.ops[i].p1);
	break;
      }
	  
  outputWidth = xSize;
  outputHeight = ySize;
  output = (unsigned char *) malloc(outputHeight*outputWidth);
  ombpl = (outputWidth + 7) / 8;
  outputMask = (unsigned char *) malloc(outputHeight*ombpl);
  memset(outputMask, 0, outputHeight*ombpl);

  /* compute the output image */
  if (c.oversamplingFactor != 0)
    oversamplingFactor = c.oversamplingFactor;

  of2 = oversamplingFactor * oversamplingFactor;
  for (y = 0; y < outputHeight; ++y)
    for (x = 0; x < outputWidth; ++x)
      {
	sum = 0;
	valid = 1;
	for (dy = 0; dy < oversamplingFactor; ++dy)
	  {
	    py = y + (((double) dy) + 0.5) / oversamplingFactor;
	    for (dx = 0; dx < oversamplingFactor; ++dx)
	      {
		px = x + (((double) dy) + 0.5) / oversamplingFactor;
		outside = 0;
		for (i = 0; i < nConstraints; ++i)
		  if (constraints[i].a * px +
		      constraints[i].b * py +
		      constraints[i].c <= 0.0)
		    {
		      outside = 1;
		      valid = 0;
		      break;
		    }
		if (!outside)
		  {
		    nx = floor(t[0][0] * px + t[0][1] * py + t[0][2]);
		    ny = floor(t[1][0] * px + t[1][1] * py + t[1][2]);
		    if (nx >= 0 && nx < w &&
			ny >= 0 && ny < h)
		      {
			sum += image[ny*w + nx];
			if ((mask[ny*mbpl + (nx >> 3)] & (0x80 >> (nx & 7))) == 0)
			  valid = 0;
		      }
		    else
		      valid = 0;
		  }
	      }
	  }
	output[y*outputWidth+x] = (sum + (of2 >> 1)) / of2;
	if (valid)
	  outputMask[y*ombpl + (x >> 3)] |= 0x80 >> (x & 7);
      }

  /* write out the output image */
  sprintf(fn, "%s%0.*d.pgm", c.outputName, c.nDigits, task.slice);
  f = fopen(fn, "w");
  if (f == NULL)
    {
      sprintf(r.message, "Could not open output file %s\n", fn);
      return;
    }
  fprintf(f, "P5\n%d %d\n255\n", outputWidth, outputHeight);
  if (fwrite(output, 1, outputHeight*outputWidth, f) != outputHeight*outputWidth)
    {
      sprintf(r.message, "Could not write output file %s\n", fn);
      return;
    }
  fclose(f);

  /* write out the output mask */
  if (c.outputMaskName[0] != '\0')
    {
      sprintf(fn, "%s%0.*d.pbm", c.outputMaskName, c.nDigits, task.slice);
      f = fopen(fn, "w");
      if (f == NULL)
	{
	  sprintf(r.message, "Could not open output mask file %s\n", fn);
	  return;
	}
      fprintf(f, "P4\n%d %d\n", outputWidth, outputHeight);
      if (fwrite(outputMask, 1, outputHeight*ombpl, f) != outputHeight*ombpl)
	{
	  sprintf(r.message, "Could not write output mask file %s\n", fn);
	  return;
	}
      fclose(f);
    }

  free(image);
  free(mask);
  free(output);
  free(outputMask);
  free(constraints);
  
  r.message[0] = '\0';
}

void
ApplyTransform(double m[2][3], double t[2][3], int n, Constraint *con)
{
  int i, j;
  double det;
  double tmp[2][3];
  double inv[2][3];

  det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
  /*  printf("det = %f\n", det); */
  inv[0][0] = m[1][1] / det;
  inv[0][1] = -m[0][1] / det;
  inv[0][2] = (m[0][1] * m[1][2] - m[1][1] * m[0][2]) / det;
  inv[1][0] = -m[1][0] / det;
  inv[1][1] = m[0][0] / det;
  inv[1][2] = (m[0][2] * m[1][0] - m[1][2] * m[0][0]) / det;

  tmp[0][0] = t[0][0] * inv[0][0] + t[0][1] * inv[1][0];
  tmp[0][1] = t[0][0] * inv[0][1] + t[0][1] * inv[1][1];
  tmp[0][2] = t[0][0] * inv[0][2] + t[0][1] * inv[1][2] + t[0][2];
  tmp[1][0] = t[1][0] * inv[0][0] + t[1][1] * inv[1][0];
  tmp[1][1] = t[1][0] * inv[0][1] + t[1][1] * inv[1][1];
  tmp[1][2] = t[1][0] * inv[0][2] + t[1][1] * inv[1][2] + t[1][2];
  for (i = 0; i < 2; ++i)
    for (j = 0; j < 3; ++j)
      t[i][j] = tmp[i][j];

  /*  printf("TRANSFORM:  %f %f %f\n", t[0][0], t[0][1], t[0][2]);
      printf("            %f %f %f\n", t[1][0], t[1][1], t[1][2]); */

  for (i = 0; i < n; ++i)
    {
      tmp[0][0] = con[i].a * inv[0][0] + con[i].b * inv[1][0];
      tmp[0][1] = con[i].a * inv[0][1] + con[i].b * inv[1][1];
      tmp[0][2] = con[i].a * inv[0][2] + con[i].b * inv[1][2] + con[i].c;
      con[i].a = tmp[0][0];
      con[i].b = tmp[0][1];
      con[i].c = tmp[0][2];
      /*      printf("CONSTRAINT %d:  %f * x + %f *y + %f > 0\n", i,
	      con[i].a, con[i].b, con[i].c); */
    }
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
  par_pkint(c.maskDigits);
  par_pkstr(c.maskName);
  par_pkstr(c.outputName);
  par_pkstr(c.outputMaskName);
  par_pkint(c.oversamplingFactor);
  par_pkint(c.nOps);
  for (i = 0; i < c.nOps; ++i)
    {
      par_pkint(c.ops[i].type);
      par_pkdouble(c.ops[i].p0);
      par_pkdouble(c.ops[i].p1);
      par_pkdouble(c.ops[i].p2);
      par_pkdouble(c.ops[i].p3);
    }
}

void
UnpackContext ()
{
  int i;

  c.nDigits = par_upkint();
  par_upkstr(c.inputName);
  c.maskDigits = par_upkint();
  par_upkstr(c.maskName);
  par_upkstr(c.outputName);
  par_upkstr(c.outputMaskName);
  c.oversamplingFactor = par_upkint();
  c.nOps = par_upkint();
  c.ops = (Op*) realloc(c.ops, c.nOps * sizeof(Op));
  for (i = 0; i < c.nOps; ++i)
    {
      c.ops[i].type = par_upkint();
      c.ops[i].p0 = par_upkdouble();
      c.ops[i].p1 = par_upkdouble();
      c.ops[i].p2 = par_upkdouble();
      c.ops[i].p3 = par_upkdouble();
    }
}

void
PackTask ()
{
  par_pkint(task.slice);
}

void
UnpackTask ()
{
  task.slice = par_upkint();
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
