#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <dirent.h>

int imageWidth = 0;
int imageHeight = 0;
unsigned char *image = 0;
unsigned char *mask = 0;
int nDigits;
char inputName[PATH_MAX];
char outputName[PATH_MAX];
int nProcessed = 0;
int section = 10000;

int Compare (const void *x, const void *y);
int ParseRange (char *s, int *pos, int *minValue, int *maxValue);
int ReadHeader (FILE *f, char *tc, int *w, int *h, int *m);

int
main (int argc, char **argv, char **envp)
{
  int i;
  int n;
  int nSkipSlices;
  int maxSkipSlices;
  int *skipSlices;
  int len;
  char skipFile[PATH_MAX];
  int minSlice, maxSlice;
  int pos;
  FILE *f;
  FILE *of;
  char inputDirName[PATH_MAX];
  char inputPrefix[PATH_MAX];
  int inputPrefixLen;
  char fn[PATH_MAX];
  char ofn[PATH_MAX];
  int minN, maxN;
  int digitLen;
  DIR *dir;
  struct dirent *de;
  char tc;
  int w, h;
  int m;
  int skipIndex;
  int minS, maxS;
  unsigned char *black;
  unsigned char *line;
  int offset;
  unsigned char *c;
  int error;
  int s, t;

  error = 0;
  inputName[0] = '\0';
  outputName[0] = '\0';
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
      fprintf(stderr, "Usage: cross -input file_prefix -output file_prefix\n");
      fprintf(stderr, "              [-skip skip_file] [-slice slice_range]\n");
      fprintf(stderr, "   where slice_range is expressed as: integer\n");
      fprintf(stderr, "                              or: integer-integer\n");
      fprintf(stderr, "                              or: integer-\n");
      fprintf(stderr, "                              or: -integer\n");
      exit(1);
    }

  /* check that at least minimal parameters were supplied */
  if (inputName[0] == '\0' || outputName[0] == '\0')
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
  for (i = strlen(inputName)-1; i >= 0 && inputName[i] != '/'; --i) ;
  if (i >= 0)
    {
      strncpy(inputDirName, inputName, i+1);
      inputDirName[i+1] = '\0';
      strcpy(inputPrefix, &inputName[i+1]);
    }
  else
    {
      strcpy(inputDirName, "./");
      strcpy(inputPrefix, inputName);
    }

  /* read the directory to look for files of the form nameNNNN.ppm  (or .pgm) */
  dir = opendir(inputDirName);
  if (dir == NULL)
    {
      fprintf(stderr, "Could not open directory %s for input files\n");
      exit(1);
    }
  inputPrefixLen = strlen(inputPrefix);
  minN = 1000000000;
  maxN = -1;
  nDigits = -1;
  while ((de = readdir(dir)) != NULL)
    {
      if (strncmp(de->d_name, inputPrefix, inputPrefixLen) == 0 &&
	  (len = strlen(de->d_name)) > 4 &&
	  strcmp(&(de->d_name[len-4]), ".pgm") == 0)
	{
	  if (sscanf(&(de->d_name[inputPrefixLen]), "%d%n",
		     &n, &digitLen) < 1)
	    continue;
	  if (nDigits < 0)
	    nDigits = digitLen;
	  else if (digitLen != nDigits)
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

  h = maxN - minN + 1;
  printf("Size of output image (WxH) is %dx%d\n", w, h);
  printf("Processing slices: ");
  fflush(stdout);
  
  black = (unsigned char *) malloc(w);
  memset(black, 0, w);
  line = (unsigned char *) malloc(w);
  /* for all slices */
  for (s = minSlice; s <= maxSlice; ++s)
    {
      sprintf(ofn, "%s%0.*d.pgm", outputName, nDigits, s);
      of = fopen(ofn, "w");
      fprintf(of, "P5\n%d %d\n255\n", w, h);

      skipIndex = 0;
      for (t = minN; t <= maxN; ++t)
	{
	  while (skipIndex < nSkipSlices && skipSlices[skipIndex] < t)
	    ++skipIndex;
	  if (skipIndex < nSkipSlices && skipSlices[skipIndex] == t)
	    {
	      fwrite(black, 1, w, of);
	      continue;
	    }
	  sprintf(fn, "%s%0.*d.pgm", inputName, nDigits, t);
	  if ((f = fopen(fn, "r")) == NULL)
	    {
	      fprintf(stderr, "Cannot open input file %s for reading\n",
		      fn);
	      exit(1);
	    }
	  offset = s*w + 19;
	  fseek(f, (long) offset, SEEK_SET);
	  if (fread(line, 1, w, f) != w)
	    {
	      fprintf(stderr, "Cannot read scan line %d from file %s\n", t, fn);
	      exit(1);
	    }
	  fclose(f);
	  fwrite(line, 1, w, of);
	  if ((nProcessed % 50) == 0 && nProcessed != 0)
	    printf(" %d \n                   ", nProcessed);
	  printf(".");
	  fflush(stdout);
	  ++nProcessed;
	}
      fclose(of);
    }

  printf(" %d\nAll slices completed.\n", nProcessed);
  exit(0);
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
