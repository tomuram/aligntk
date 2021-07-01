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
  char inputDirName[PATH_MAX];
  char outputDirName[PATH_MAX];
  int deltaLevel;
  char fn[PATH_MAX];
  char ofn[PATH_MAX];
  FILE *f;
  FILE *of;
  float factor;
  float x, y, rx, ry;
  int error;
  char line[256];

  error = 0;
  factor = -1;
  inputDirName[0] = '\0';
  outputDirName[0] = '\0';
  for (i = 1; i < argc; ++i)
    if (strcmp(argv[i], "-input") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-input error\n");
	    break;
	  }
	strcpy(inputDirName, argv[i]);
      }
    else if (strcmp(argv[i], "-output") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-output error\n");
	    break;
	  }
	strcpy(outputDirName, argv[i]);
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
      fprintf(stderr, "Usage: scale_cpts3 -delta_level int\n");
      fprintf(stderr, "                  -input directory\n");
      fprintf(stderr, "                  -output directory\n");
      exit(1);
    }
	
  /* read the directory to look for files of the form name.pts */
  dir = opendir(inputDirName);
  if (dir == NULL)
    {
      fprintf(stderr, "Could not open directory %s for input files\n",
	      inputDirName);
      exit(1);
    }
  while ((de = readdir(dir)) != NULL)
    {
      if ((len = strlen(de->d_name)) <= 4 ||
	  strcmp(&(de->d_name[len-4]), ".pts") != 0)
	continue;
      
      sprintf(fn, "%s/%s", inputDirName, de->d_name);
      sprintf(ofn, "%s/%s", outputDirName, de->d_name);
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
	  if (sscanf(line, "%f%f%f%f",
		     &x, &y, &rx, &ry) != 4)
	    {
	      fprintf(stderr, "Bad line found in file %s\n", fn);
	      exit(1);
	    }
	  if (deltaLevel >= 0)
	    {
	      factor = 1 << deltaLevel;
	      x = x / factor;
	      y = y / factor;
	      rx = rx / factor;
	      ry = ry / factor;
	    }
	  else
	    {
	      factor = 1 << (-deltaLevel);
	      x = factor * x;
	      y = factor * y;
	      rx = factor * rx;
	      ry = factor * ry;
	    }
	  fprintf(of, "%d %d %d %d\n",
		  (int) floor(x + 0.5),
		  (int) floor(y + 0.5),
		  (int) floor(rx + 0.5),
		  (int) floor(ry + 0.5));
	}
      fclose(f);
      fclose(of);
    }
  closedir(dir);
}
