#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <dirent.h>
#include <mpi.h>

#include "par.h"

typedef struct Context {
  /* NOTE: any new fields added to this struct should
     be also added to PackContext and UnpackContext */ 
  int nDigits;
  char inputName[PATH_MAX];
  char outputName[PATH_MAX];
  int factor;
  int ppm;
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
long imageWidth = 0;
long imageHeight = 0;
long outputWidth = 0;
long outputHeight = 0;
unsigned char *image = 0;
unsigned char *output = 0;

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
  c.factor = 2;
  c.inputName[0] = '\0';
  c.outputName[0] = '\0';
  skipFile[0] = '\0';
  minSlice = 0;
  maxSlice = 1000000000;
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
      fprintf(stderr, "Usage: reduce -input file_prefix -output file_prefix\n");
      fprintf(stderr, "              [-factor reduction_factor]\n");
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

  /* read the directory to look for files of the form nameNNNN.ppm  (or .pgm) */
  dir = opendir(inputDirName);
  if (dir == NULL)
    {
      fprintf(stderr, "Could not open directory %s for input files\n",
	      inputDirName);
      exit(1);
    }
  inputPrefixLen = strlen(inputPrefix);
  minN = 1000000000;
  maxN = -1;
  c.nDigits = -1;
  c.ppm = -1;
  while ((de = readdir(dir)) != NULL)
    {
      if (strncmp(de->d_name, inputPrefix, inputPrefixLen) == 0 &&
	  (len = strlen(de->d_name)) > 4 &&
	  (strcmp(&(de->d_name[len-4]), ".pgm") == 0 ||
	   strcmp(&(de->d_name[len-4]), ".ppm") == 0))
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
	  if (c.ppm < 0)
	    c.ppm = (strcmp(&(de->d_name[len-4]), ".ppm") == 0);
	  else
	    if ((c.ppm ^ (strcmp(&(de->d_name[len-4]), ".ppm") == 0)) != 0)
	      {
		fprintf(stderr, "Mix of .pgm and .ppm files not allowed.\n");
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
  /* nothing necessary here */
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
  int x, y;
  int dx, dy;
  long a;
  long v;
  long vr, vg, vb;
  int offset;
  int ix, iy;
  int bpp;

  /* get the image */
  sprintf(fn, "%s%0.*d.%s", c.inputName, c.nDigits, t.slice,
	  c.ppm ? "ppm" : "pgm");
  f = fopen(fn, "r");
  if (f == NULL)
    {
      sprintf(r.message, "Worker could not open image file %s\n", fn);
      return;
    }
  printf("%f: Opened the input file\n", MPI_Wtime());

  if (!ReadHeader(f, &tc, &w, &h, &m))
    {
      sprintf(r.message, "Image file not binary pgm or ppm: %s\n", fn);
      fclose(f);
      return;
    }
  if (m != 255 || (tc != '5' && tc != '6'))
    {
      sprintf(r.message, "Image file not 8-bit pgm or ppm: %s\n", fn);
      fclose(f);
      return;
    }
  
  bpp = c.ppm ? 3 : 1;
  if (w != imageWidth || h != imageHeight)
    {
      imageWidth = w;
      imageHeight = h;
      image = (unsigned char *) realloc(image, imageHeight * imageWidth * bpp);
      outputWidth = w / c.factor;
      outputHeight = h / c.factor;
      output = (unsigned char *) realloc(output, outputHeight * outputWidth * bpp);
    }

  printf("%f: Reading the image\n", MPI_Wtime());
  if (fread(image, 1, imageWidth*imageHeight*bpp, f) !=
      imageWidth*imageHeight*bpp)
    {
      fclose(f);
      sprintf(r.message, "Image file %s apparently truncated\n", fn);
      return;
    }
  printf("%f: Finished reading the image (size %ld x %ld)\n", MPI_Wtime(),
	 imageWidth, imageHeight);
  fclose(f);
  printf("%f: Closed the image file\n", MPI_Wtime());
  

  printf("%d\n", c.factor);
  printf("%ld %ld\n", outputWidth, outputHeight);

  a = c.factor * c.factor;
  offset = a / 2;
  for (y = 0; y < outputHeight; ++y)
    {
      iy = y * c.factor;
      for (x = 0; x < outputWidth; ++x)
	{
	  ix = x * c.factor;
	  if (bpp == 1)
	    {
	      v = 0;
	      for (dy = 0; dy < c.factor; ++dy)
		for (dx = 0; dx < c.factor; ++dx)
		  v += image[(iy + dy)*imageWidth + ix + dx];
	      output[y * outputWidth + x] = (v + offset) / a;
	    }
	  else
	    {
	      vr = 0;
	      vg = 0;
	      vb = 0;
	      for (dy = 0; dy < c.factor; ++dy)
		for (dx = 0; dx < c.factor; ++dx)
		  {
		    vr += image[((iy + dy)*imageWidth + ix + dx)*bpp];
		    vg += image[((iy + dy)*imageWidth + ix + dx)*bpp + 1];
		    vb += image[((iy + dy)*imageWidth + ix + dx)*bpp + 2];
		  }
	      output[(y * outputWidth + x)*bpp] = (vr + offset) / a;
	      output[(y * outputWidth + x)*bpp+1] = (vg + offset) / a;
	      output[(y * outputWidth + x)*bpp+2] = (vb + offset) / a;
	    }
	}	
    }
  printf("%f: Finished the computation\n", MPI_Wtime());


  /* write out the single image */
  sprintf(fn, "%s%0.*d.%s", c.outputName, c.nDigits, t.slice,
	  c.ppm ? "ppm" : "pgm");
  f = fopen(fn, "w");
  if (f == NULL)
    {
      sprintf(r.message, "Could not open output file %s\n", fn);
      return;
    }
  printf("%f: Opened the output file\n", MPI_Wtime());


  fprintf(f, "P%c\n%ld %ld\n255\n", c.ppm ? '6' : '5',
	  outputWidth, outputHeight);
  if (fwrite(output, outputHeight*outputWidth*bpp, 1, f) != 1)
    {
      sprintf(r.message, "Could not write output file %s\n", fn);
      return;
    }
  printf("%f: Finished writing the output file\n", MPI_Wtime());
  fclose(f);
  printf("%f: Closed the output file\n", MPI_Wtime());

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
  par_pkint(c.factor);
  par_pkint(c.ppm);
}

void
UnpackContext ()
{
  int i;

  c.nDigits = par_upkint();
  par_upkstr(c.inputName);
  par_upkstr(c.outputName);
  c.factor = par_upkint();
  c.ppm = par_upkint();
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
