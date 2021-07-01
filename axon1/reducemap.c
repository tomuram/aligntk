#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdarg.h>
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
long mpw = 0;
long mph = 0;
long mpwp;
long mphp;
long ompw = 0;
long omph = 0;
long ompwp;
long omphp;
float *mapX = 0;
float *mapY = 0;
float *mapC = 0;
float *omapX = 0;
float *omapY = 0;
float *omapC = 0;

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
      fprintf(stderr, "Usage: reducemap -input file_prefix -output file_prefix\n");
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
    Error("Both -input and -output parameters must be specified.\n");

  /* get a sorted list of the slices to skip */
  nSkipSlices = 0;
  skipSlices = (int *) malloc(sizeof(int));
  maxSkipSlices = 1;
  if (skipFile[0] != '\0')
    {
      f = fopen(skipFile, "r");
      if (f == NULL)
	Error("Could not open skip file %s\n", skipFile);
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
	Error("Read error while reading skip file %s\n", skipFile);
      
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
    Error("Could not open directory %s for input files\n",
	  inputDirName);
  inputPrefixLen = strlen(inputPrefix);
  minN = 1000000000;
  maxN = -1;
  c.nDigits = -1;
  while ((de = readdir(dir)) != NULL)
    {
      if (strncmp(de->d_name, inputPrefix, inputPrefixLen) == 0 &&
	  (len = strlen(de->d_name)) > 4 &&
	  strcmp(&(de->d_name[len-4]), ".map") == 0)
	{
	  if (sscanf(&(de->d_name[inputPrefixLen]), "%d%n",
		     &n, &digitLen) < 1)
	    continue;
	  if (c.nDigits < 0)
	    c.nDigits = digitLen;
	  else if (digitLen != c.nDigits)
	    Error("Inconsistent number of digits in input file names: %s\n", de->d_name);
	  if (n < minN)
	    minN = n;
	  if (n > maxN)
	    maxN = n;
	}
    }
  closedir(dir);
  if (maxN < 0)
    Error("No input map files found.\n");
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
  int offset;
  int ix, iy;
  int mapHeader[3];

  /* get the map */
  sprintf(fn, "%s%0.*d.map", c.inputName, c.nDigits, t.slice);
  f = fopen(fn, "r");
  if (f == NULL)
    {
      sprintf(r.message, "Worker could not open map file %s\n", fn);
      return;
    }

  if (fread(mapHeader, sizeof(int), 3, f) != 3)
    {
      sprintf(r.message, "Could not read map file header: %s\n", fn);
      fclose(f);
      return;
    }
  if (mapHeader[0] != mpw || mapHeader[1] != mph)
    {
      mpw = mapHeader[0];
      mph = mapHeader[1];
      mpwp = mpw + 2;
      mphp = mph + 2;
      mapX = (float *) realloc(mapX, mpwp * mphp * sizeof(float));
      mapY = (float *) realloc(mapY, mpwp * mphp * sizeof(float));
      mapC = (float *) realloc(mapC, mpwp * mphp * sizeof(float));
      ompw = mpw / c.factor;
      omph = mph / c.factor;
      ompwp = ompw + 2;
      omphp = omph + 2;
      omapX = (float *) realloc(omapX, ompwp * omphp * sizeof(float));
      omapY = (float *) realloc(omapY, ompwp * omphp * sizeof(float));
      omapC = (float *) realloc(omapC, ompwp * omphp * sizeof(float));
    }

  if (fread(mapX, sizeof(float), mpwp*mphp, f) != mpwp*mphp ||
      fread(mapY, sizeof(float), mpwp*mphp, f) != mpwp*mphp ||
      fread(mapC, sizeof(float), mpwp*mphp, f) != mpwp*mphp)
    {
      fclose(f);
      sprintf(r.message, "Map file %s apparently truncated\n", fn);
      return;
    }
  fclose(f);

  for (y = 0; y < omph; ++y)
    {
      iy = (int) floor(c.factor * y + (c.factor - 1.0) / 2.0);
      for (x = 0; x < ompw; ++x)
	{
	  ix = (int) floor(c.factor * x + (c.factor - 1.0) / 2.0);
	  omapX[(y + 1) * ompwp + (x + 1)] =
	    (mapX[(iy + 1) * mpwp + (ix + 1)] +
	     mapX[(iy + 1) * mpwp + (ix + 2)] +
	     mapX[(iy + 2) * mpwp + (ix + 1)] +
	     mapX[(iy + 2) * mpwp + (ix + 2)]) * 0.25;
	  omapY[(y + 1) * ompwp + (x + 1)] =
	    (mapY[(iy + 1) * mpwp + (ix + 1)] +
	     mapY[(iy + 1) * mpwp + (ix + 2)] +
	     mapY[(iy + 2) * mpwp + (ix + 1)] +
	     mapY[(iy + 2) * mpwp + (ix + 2)]) * 0.25;
	  omapC[(y + 1) * ompwp + (x + 1)] =
	    (mapC[(iy + 1) * mpwp + (ix + 1)] +
	     mapC[(iy + 1) * mpwp + (ix + 2)] +
	     mapC[(iy + 2) * mpwp + (ix + 1)] +
	     mapC[(iy + 2) * mpwp + (ix + 2)]) * 0.25;
	}	
    }
  /* extrapolate periphery */
  for (x = 0; x < ompw; ++x)
    {
      omapX[x + 1] = 2.0 * omapX[ompwp + x + 1]  - omapX[2*ompwp + x + 1];
      omapY[x + 1] = 2.0 * omapY[ompwp + x + 1]  - omapY[2*ompwp + x + 1];
      omapC[x + 1] = omapC[ompwp + x + 1]  * omapC[2*ompwp + x + 1];

      omapX[(omph + 1) * ompwp + x + 1] =
	2.0 * omapX[omph * ompwp + x + 1] - omapX[(omph - 1) * ompwp + x + 1];
      omapY[(omph + 1) * ompwp + x + 1] =
	2.0 * omapY[omph * ompwp + x + 1] - omapY[(omph - 1) * ompwp + x + 1];
      omapC[(omph + 1) * ompwp + x + 1] =
	omapC[omph * ompwp + x + 1] * omapC[(omph - 1) * ompwp + x + 1];
    }
  for (y = 0; y < omph; ++y)
    {
      omapX[(y + 1) * ompwp] =
	2.0 * omapX[(y + 1) * ompwp + 1]  - omapX[(y + 1) * ompwp + 2];
      omapY[(y + 1) * ompwp] =
	2.0 * omapY[(y + 1) * ompwp + 1]  - omapY[(y + 1) * ompwp + 2];
      omapC[(y + 1) * ompwp] =
	omapC[(y + 1) * ompwp + 1]  * omapC[(y + 1) * ompwp + 2];

      omapX[(y + 1) * ompwp + ompw + 1] =
	2.0 * omapX[(y + 1) * ompwp + ompw] -
	omapX[(y + 1) * ompwp + ompw - 1];
      omapY[(y + 1) * ompwp + ompw + 1] =
	2.0 * omapY[(y + 1) * ompwp + ompw] -
	omapY[(y + 1) * ompwp + ompw - 1];
      omapC[(y + 1) * ompwp + ompw + 1] =
	omapC[(y + 1) * ompwp + ompw] *
	omapC[(y + 1) * ompwp + ompw - 1];
    }
  omapX[0] = 0.5 * (2.0 * omapX[1] -
		    omapX[2] +
		    2.0 * omapX[ompwp] -
		    omapX[2*ompwp]);
  omapY[0] = 0.5 * (2.0 * omapY[1] -
		    omapY[2] +
		    2.0 * omapY[ompwp] -
		    omapY[2*ompwp]);
  omapC[0] = 0.5 * (omapC[1] *
		    omapC[2] +
		    omapC[ompwp] *
		    omapC[2*ompwp]);
  omapX[ompw + 1] = 0.5 * (2.0 * omapX[ompw] -
			   omapX[ompw-1] +
			   2.0 * omapX[ompwp + ompw + 1] -
			   omapX[2*ompwp + ompw + 1]);
  omapY[ompw + 1] = 0.5 * (2.0 * omapY[ompw] -
			   omapY[ompw-1] +
			   2.0 * omapY[ompwp + ompw + 1] -
			   omapY[2*ompwp + ompw + 1]);
  omapC[ompw + 1] = 0.5 * (omapC[ompw] *
			   omapC[ompw-1] +
			   omapC[ompwp + ompw + 1] *
			   omapC[2*ompwp + ompw + 1]);
  omapX[(omph + 1) * ompwp] = 0.5 * (2.0 * omapX[omph * ompwp] -
				     omapX[(omph - 1) * ompwp] +
				     2.0 * omapX[(omph + 1) * ompwp + 1] -
				     omapX[(omph + 1) * ompwp + 2]);
  omapY[(omph + 1) * ompwp] = 0.5 * (2.0 * omapY[omph * ompwp] -
				     omapY[(omph - 1) * ompwp] +
				     2.0 * omapY[(omph + 1) * ompwp + 1] -
				     omapX[(omph + 1) * ompwp + 2]);
  omapC[(omph + 1) * ompwp] = 0.5 * (omapC[omph * ompwp] *
				     omapC[(omph - 1) * ompwp] +
				     omapC[(omph + 1) * ompwp + 1] *
				     omapC[(omph + 1) * ompwp + 2]);
  omapX[(omph + 1) * ompwp + ompw + 1] =
    0.5 * (2.0 * omapX[omph * ompwp + ompw + 1] -
	   omapX[(omph - 1) * ompwp + ompw + 1] +
	   2.0 * omapX[(omph + 1) * ompwp + ompw] -
	   omapX[(omph + 1) * ompwp + ompw - 1]);
  omapY[(omph + 1) * ompwp + ompw + 1] =
    0.5 * (2.0 * omapY[omph * ompwp + ompw + 1] -
	   omapY[(omph - 1) * ompwp + ompw + 1] +
	   2.0 * omapY[(omph + 1) * ompwp + ompw] -
	   omapY[(omph + 1) * ompwp + ompw - 1]);
  omapC[(omph + 1) * ompwp + ompw + 1] =
    0.5 * (omapC[omph * ompwp + ompw + 1] *
	   omapC[(omph - 1) * ompwp + ompw + 1] +
	   omapC[(omph + 1) * ompwp + ompw] *
	   omapC[(omph + 1) * ompwp + ompw - 1]);
  
  /* scale values to output map size */
  for (y = -1; y <= omph; ++y)
    for (x = -1; x <= ompw; ++x)
      {
	omapX[(y + 1) * ompwp + (x + 1)] =
	  (2.0 * omapX[(y + 1) * ompwp + (x + 1)] - c.factor + 1.0) /
	  (2.0 * c.factor);
	omapY[(y + 1) * ompwp + (x + 1)] =
	  (2.0 * omapY[(y + 1) * ompwp + (x + 1)] - c.factor + 1.0) /
	  (2.0 * c.factor);
      }

  /* write out the map */
  sprintf(fn, "%s%0.*d.map", c.outputName, c.nDigits, t.slice);
  f = fopen(fn, "w");
  if (f == NULL)
    {
      sprintf(r.message, "Could not open output file %s\n", fn);
      return;
    }
  mapHeader[0] = ompw;
  mapHeader[1] = omph;
  if (fwrite(mapHeader, sizeof(int), 3, f) != 3 ||
      fwrite(omapX, sizeof(float), ompwp*omphp, f) != ompwp * omphp ||
      fwrite(omapY, sizeof(float), ompwp*omphp, f) != ompwp * omphp ||
      fwrite(omapC, sizeof(float), ompwp*omphp, f) != ompwp * omphp)
    {
      fclose(f);
      sprintf(r.message, "Could not write output file %s\n", fn);
      return;
    }
  fclose(f);

  r.message[0] = '\0';
}

void
PackContext ()
{
  int i;

  par_pkint(c.nDigits);
  par_pkstr(c.inputName);
  par_pkstr(c.outputName);
  par_pkint(c.factor);
}

void
UnpackContext ()
{
  int i;

  c.nDigits = par_upkint();
  par_upkstr(c.inputName);
  par_upkstr(c.outputName);
  c.factor = par_upkint();
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

void Error (char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
  fflush(stderr);
  abort();
}
