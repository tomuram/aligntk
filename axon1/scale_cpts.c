#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <sys/types.h>
#include <dirent.h>
#include <stdlib.h>
#include <math.h>

int
main (int argc, char **argv, char **envp)
{
  int i;
  DIR *dir;
  struct dirent *de;
  int len;
  int digitLen;
  int inputDigits;
  int n;
  int inputPrefixLen;
  int pos;
  char inputName[PATH_MAX];
  char inputDirName[PATH_MAX];
  char inputPrefix[PATH_MAX];
  char outputName[PATH_MAX];
  int deltaLevel;
  char fn[PATH_MAX];
  char ofn[PATH_MAX];
  FILE *f;
  FILE *of;
  float factor;
  float x, y, rx, ry;
  int r;
  int error;
  char line[256];

  error = 0;
  factor = -1;
  inputName[0] = '\0';
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
    else if (strcmp(argv[i], "-delta_level") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d", &deltaLevel) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else
      error = 1;

  if (error)
    {
      fprintf(stderr, "Usage: scale_cpts -delta_level int\n");
      fprintf(stderr, "                  -input file_prefix\n");
      fprintf(stderr, "                  -output file_prefix\n");
      exit(1);
    }
	
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

  /* read the directory to look for files of the form nameNNNN.pts */
  dir = opendir(inputDirName);
  if (dir == NULL)
    {
      fprintf(stderr, "Could not open directory %s for input files\n",
	      inputDirName);
      exit(1);
    }
  inputDigits = -1;
  inputPrefixLen = strlen(inputPrefix);
  while ((de = readdir(dir)) != NULL)
    {
      if (strncmp(de->d_name, inputPrefix, inputPrefixLen) != 0 ||
	  (len = strlen(de->d_name)) <= 4 ||
	  strcmp(&(de->d_name[len-4]), ".pts") != 0)
	continue;
      
      len = -1;
      if (sscanf(&(de->d_name[inputPrefixLen]), "%d%n.pts%n",
		 &n, &digitLen, &len) < 1 ||
	  inputPrefixLen + len != strlen(de->d_name))
	continue;

      if (inputDigits < 0)
	inputDigits = digitLen;
      else if (digitLen != inputDigits)
	{
	  fprintf(stderr, "Inconsistent number of digits in input file names: %s\n", de->d_name);
	  closedir(dir);
	  exit(1);
	}

      sprintf(fn, "%s%0.*d.pts", inputName, inputDigits, n);
      sprintf(ofn, "%s%0.*d.pts", outputName, inputDigits, n);
      f = fopen(fn, "r");
      if (f == NULL)
	{
	  fprintf(stderr, "Could not open input file %s for reading.\n", fn);
	  exit(1);
	}
      of = fopen(ofn, "w");
      if (of == NULL)
	{
	  fprintf(stderr, "Could not open output file %s for writing.\n", ofn);
	  exit(1);
	}
      while (fgets(line, 255, f) != NULL)
	{
	  if (sscanf(line, "%f%f%d%f%f",
		     &x, &y, &r, &rx, &ry) != 5)
	    {
	      fprintf(stderr, "Bad line found in file %s\n", fn);
	      exit(1);
	    }
	  if (deltaLevel >= 0)
	    {
	      factor = 1 << deltaLevel;
	      x = (2.0 * x - factor + 1.0) / (2.0 * factor);
	      y = (2.0 * x - factor + 1.0) / (2.0 * factor);
	      rx = (2.0 * rx - factor + 1.0) / (2.0 * factor);
	      ry = (2.0 * ry - factor + 1.0) / (2.0 * factor);
	    }
	  else
	    {
	      factor = 1 << (-deltaLevel);
	      x = factor * x + (factor - 1.0) / 2.0;
	      y = factor * y + (factor - 1.0) / 2.0;
	      rx = factor * rx + (factor - 1.0) / 2.0;
	      ry = factor * ry + (factor - 1.0) / 2.0;
	    }
	  fprintf(of, "%d %d %d %d %d\n",
		  (int) floor(x + 0.5),
		  (int) floor(y + 0.5),
		  r,
		  (int) floor(rx + 0.5),
		  (int) floor(ry + 0.5));
	}
      fclose(f);
      fclose(of);
    }
  closedir(dir);
}
