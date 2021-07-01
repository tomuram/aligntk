#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <dirent.h>

#include "par.h"

typedef struct Context {
  /* NOTE: any new fields added to this struct should
     be also added to PackContext and UnpackContext */ 
  int nDigits;
  char inputName[PATH_MAX];
  char outputName[PATH_MAX];
  double radius;
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
Result r;

/* GLOBAL VARIABLES FOR WORKER */
int imageWidth = 0;
int imageHeight = 0;
unsigned char *image = 0;
unsigned char *output = 0;
int kernelOffset;
int kernelHeight;
int kernelWidth;
double *kernel = 0;


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
  int minS, maxS;
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

  error = 0;
  c.inputName[0] = '\0';
  c.outputName[0] = '\0';
  skipFile[0] = '\0';
  c.radius = 5.0;
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
    else if (strcmp(argv[i], "-radius") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%f", &c.radius) != 1)
	  {
	    fprintf(stderr, "-radius error\n");
	    error = 1;
	    break;
	  }
      }
  else
    {
      fprintf(stderr, "Unknown option: %s\n", argv[i]);
      error = 1;
    }
      

  if (error)
    {
      fprintf(stderr, "Usage: filter -input file_prefix -output file_prefix\n");
      fprintf(stderr, "              [-skip skip_file] [-slice slice_range]\n");
      fprintf(stderr, "              [-radius filter_radius_in_pixels]\n");
      fprintf(stderr, "   where ranges are expressed as: integer\n");
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

  /* read the directory to look for files of the form nameNNNN.pgm */
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
	      if (tc != '5')
		{
		  fclose(f);
		  fprintf(stderr, "Image file of wrong type: %s  %c\n", de->d_name, tc);
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
  double xv, yv, rv;
  int x, y;
  double v;
  int kx, ky;
  int ix, iy;
  double minV, maxV;
  int dx, dy;

  kernelOffset = ceil(5.0 * c.radius);
  kernelWidth = 2 * kernelOffset + 1;
  kernelHeight = 2 * kernelOffset + 1;
  kernel = (double*) malloc(kernelHeight * kernelWidth * sizeof(double));
  memset(kernel, 0, kernelHeight*kernelWidth*sizeof(double));
  for (y = 0; y < kernelHeight; ++y)
    for (dy = 0; dy < 100; ++dy)
      {
	yv = (y - kernelOffset) + (dy - 49.5) * 0.01;
	for (x = 0; x < kernelWidth; ++x)
	  for (dx = 0; dx < 100; ++dx)
	    {
	      xv = (x - kernelOffset) + (dx - 49.5) * 0.01;
	      rv = sqrt(xv * xv + yv * yv);
	      kernel[y * kernelWidth + x] +=
		0.0001 * (1.0 - rv * rv / c.radius / c.radius) *
		exp(- rv * rv / c.radius / c.radius);
	    }
      }
  minV = 0.0;
  maxV = 0.0;
  for (y = 0; y < kernelHeight; ++y)
    for (x = 0; x < kernelWidth; ++x)
      {
	if (kernel[y * kernelWidth + x] >= 0.0)
	  maxV += 255.0 * kernel[y * kernelWidth + x];
	else
	  minV += 255.0 * kernel[y * kernelWidth + x];
      }
  printf("KERNEL minV = %f maxV = %f\n", minV, maxV);
}

void
WorkerTask ()
{
  FILE *f;
  int w, h;
  char fn[PATH_MAX];
  int n;
  char tc;
  int m;
  int b;
  int x, y;
  int dx, dy;
  double v;
  int kx, ky;
  int ix, iy;
  double minV, maxV;

  /* get the image */
  sprintf(fn, "%s%0.*d.pgm", c.inputName, c.nDigits, t.slice);
  f = fopen(fn, "r");
  if (f == NULL)
    {
      sprintf(r.message, "Worker could not open image file %s\n", fn);
      return;
    }

  if (!ReadHeader(f, &tc, &w, &h, &m))
    {
      sprintf(r.message, "Image file not binary pgm: %s\n", fn);
      fclose(f);
      return;
    }
  if (tc != '5')
    {
      sprintf(r.message, "Image file not binary pgm: %s\n", fn);
      fclose(f);
      return;
    }
  
  if (w != imageWidth || h != imageHeight)
    {
      image = (unsigned char *) realloc(image, w*h);
      output = (unsigned char *) realloc(output, w*h);
      imageWidth = w;
      imageHeight = h;
    }

  if (fread(image, 1, w*h, f) != w*h)
    {
      fclose(f);
      sprintf(r.message, "Image file %s apparently truncated (%d %d)\n", fn,
	      w, h);
      return;
    }
  fclose(f);

  /* go through the image pixels and compute the value of the
     filtered image */
  for (y = 0; y < h; ++y)
    for (x = 0; x < w; ++x)
      {
	v = 0.0;
	minV = 0.0;
	maxV = 0.0;
	for (ky = 0; ky < kernelHeight; ++ky)
	  {
	    iy = y + ky - kernelOffset;
	    if (iy < 0 || iy >= h)
	      continue;
	    for (kx = 0; kx < kernelWidth; ++kx)
	      {
		ix = x + kx - kernelOffset;
		if (ix < 0 || ix >= w)
		  continue;
		if (kernel[ky * kernelWidth + kx] >= 0.0)
		  maxV += 255.0 * kernel[ky * kernelWidth + kx];
		else
		  minV += 255.0 * kernel[ky * kernelWidth + kx];
		v += kernel[ky * kernelWidth + kx] * image[iy*w + ix];
	      }
	  }
	output[y * w + x] = floor((v - minV) / (maxV - minV) * 255.999999);
      }

  /* write out the single image */
  sprintf(fn, "%s%0.*d.pgm", c.outputName, c.nDigits, t.slice);
  f = fopen(fn, "w");
  if (f == NULL)
    {
      sprintf(r.message, "Could not open output file %s\n", fn);
      return;
    }

  fprintf(f, "P5\n%d %d\n255\n", w, h);
  if (fwrite(output, 1, h*w, f) != h*w)
    {
      sprintf(r.message, "Could not write output file %s\n", fn);
      return;
    }
  fclose(f);

  r.message[0] = '\0';
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
  par_pkdouble(c.radius);
}

void
UnpackContext ()
{
  int i;

  c.nDigits = par_upkint();
  par_upkstr(c.inputName);
  par_upkstr(c.outputName);
  c.radius = par_upkdouble();
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
  par_pkstr(r.message);
}

void
UnpackResult ()
{
  par_upkstr(r.message);
}
