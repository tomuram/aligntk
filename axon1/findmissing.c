#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <sys/types.h>
#include <dirent.h>
#include <stdlib.h>

int Compare (const void *x, const void *y);

int
main (int argc, char **argv, char **envp)
{
  int i;
  DIR *dir;
  struct dirent *de;
  int nSlices;
  int maxSlices;
  int *slices;
  int minN, maxN;
  int len;
  int digitLen;
  int nDigits;
  int n;
  int pgm;
  int ppm;
  int inputPrefixLen;
  int pos;
  char inputName[PATH_MAX];
  char inputDirName[PATH_MAX];
  char inputPrefix[PATH_MAX];

  if (argc != 2)
    {
      fprintf(stderr, "Usage: findmissing file_prefix\n");
      exit(1);
    }
  strcpy(inputName, argv[1]);

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
      strcpy(inputDirName, ".");
      strcpy(inputPrefix, inputName);
    }

  /* read the directory to look for files of the form nameNNNN.pgm or
     nameNNNN.ppm */
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
  nDigits = -1;
  ppm = 0;
  pgm = 0;
  nSlices = 0;
  slices = (int *) malloc(sizeof(int));
  maxSlices = 1;
  while ((de = readdir(dir)) != NULL)
    {
      if (strncmp(de->d_name, inputPrefix, inputPrefixLen) != 0 ||
	  (len = strlen(de->d_name)) <= 4 ||
	  (strcmp(&(de->d_name[len-4]), ".pgm") != 0 &&
	   strcmp(&(de->d_name[len-4]), ".ppm") != 0))
	continue;

      len = -1;
      if (sscanf(&(de->d_name[inputPrefixLen]), "%d%n.pgm%n",
		 &n, &digitLen, &len) >= 1 &&
	  inputPrefixLen + len == strlen(de->d_name))
	{
	  if (ppm)
	    {
	      fprintf(stderr, "Error: both .pgm and .ppm files found\n");
	      closedir(dir);
	      exit(1);
	    }
	  pgm = 1;
	}
      else
	{
	  len = -1;
	  if (sscanf(&(de->d_name[inputPrefixLen]), "%d%n.ppm%n",
		     &n, &digitLen, &len) >= 1 &&
	      inputPrefixLen + len == strlen(de->d_name))
	    {
	      if (pgm)
		{
		  fprintf(stderr, "Error: both .pgm and .ppm files found\n");
		  closedir(dir);
		  exit(1);
		}
	      ppm = 1;
	    }
	  else
	    continue;
	}
	  
      if (nDigits < 0)
	nDigits = digitLen;
      else if (digitLen != nDigits)
	{
	  fprintf(stderr, "Inconsistent number of digits in input file names: %s\n", de->d_name);
	  closedir(dir);
	  exit(1);
	}
      if (n < minN)
	minN = n;
      if (n > maxN)
	maxN = n;
      if (++nSlices > maxSlices)
	{
	  maxSlices *= 2;
	  slices = (int *) realloc(slices, maxSlices * sizeof(int));
	}
      slices[nSlices-1] = n;
    }
  closedir(dir);

  if (maxN < 0)
    {
      fprintf(stderr, "No input image files found.\n");
      exit(1);
    }
  qsort(slices, nSlices, sizeof(int), Compare);
  pos = 0;
  for (i = minN; i <= maxN; ++i)
    {
      if (pos < nSlices && slices[pos] < i)
	{
	  fprintf(stderr, "Error: more than one slice numbered %d found.\n",
		  slices[pos]);
	  exit(1);
	}
      if (pos < nSlices && slices[pos] == i)
	{
	  ++pos;
	  continue;
	}
      printf("%d\n", i);
    }
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
