/*
 * construct_mapslist.c - filter a set of pairwise maps
 *                        according to a threshold and
 *                        print a maplist to stdout
 *
 *   Copyright (c) 2010 Pittsburgh Supercomputing Center
 * 
 *  This program is distributed in the hope that it will    *
 *  be useful, but WITHOUT ANY WARRANTY; without even the   *
 *  implied warranty of MERCHANTABILITY or FITNESS FOR A    *
 *  PARTICULAR PURPOSE.  Neither Carnegie Mellon University *
 *  nor any of the authors assume any liability for         *
 *  damages, incidental or otherwise, caused by the         *
 *  installation or use of this software.                   *
 *                                                          *
 *  CLINICAL APPLICATIONS ARE NOT RECOMMENDED, AND THIS     *
 *  SOFTWARE HAS NOT BEEN EVALUATED BY THE UNITED STATES    *
 *  FDA FOR ANY CLINICAL USE.                               *
 *
 *  HISTORY
 *    2010  Written by Greg Hood (ghood@psc.edu)
 */

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>

#define LINE_LENGTH	256

char pairListName[PATH_MAX];
char mapsName[PATH_MAX];
double minCorrelation = -1.0;
double maxCorrelation = 1.01;
double minDistortion = 0.0;
double maxDistortion = 1.0e30;
double weight = 1.0;

int
main (int argc, char **argv)
{
  char line[LINE_LENGTH];
  int i, j, k;
  int error;
  FILE *f, *sf;
  int nItems;
  char imageName[PATH_MAX];
  char scoreName[PATH_MAX];
  char imageName0[PATH_MAX], imageName1[PATH_MAX];
  char pairName[PATH_MAX];
  int minX0, maxX0, minY0, maxY0;
  int minX1, maxX1, minY1, maxY1;
  double energy, correlation, distortion, correspondence, constraining;

  error = 0;
  pairListName[0] = '\0';
  mapsName[0] = '\0';

  for (i = 1; i < argc; ++i)
    if (strcmp(argv[i], "-pairs") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-pairs error\n");
	    break;
	  }
	strcpy(pairListName, argv[i]);
      }
    else if (strcmp(argv[i], "-maps") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-maps error\n");
	    break;
	  }
	strcpy(mapsName, argv[i]);
      }
    else if (strcmp(argv[i], "-min_correlation") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%lf", &minCorrelation) != 1)
	  {
	    error = 1;
	    fprintf(stderr, "-min_correlation error\n");
	    break;
	  }
      }
    else if (strcmp(argv[i], "-max_correlation") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%lf", &maxCorrelation) != 1)
	  {
	    error = 1;
	    fprintf(stderr, "-max_correlation error\n");
	    break;
	  }
      }
    else if (strcmp(argv[i], "-min_distortion") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%lf", &minDistortion) != 1)
	  {
	    error = 1;
	    fprintf(stderr, "-min_distortion error\n");
	    break;
	  }
      }
    else if (strcmp(argv[i], "-max_distortion") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%lf", &maxDistortion) != 1)
	  {
	    error = 1;
	    fprintf(stderr, "-max_distortion error\n");
	    break;
	  }
      }
    else if (strcmp(argv[i], "-weight") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%lf", &weight) != 1)
	  {
	    error = 1;
	    fprintf(stderr, "-weight error\n");
	    break;
	  }
      }
    else
      {
	fprintf(stderr, "Unrecognized option: %s\n", argv[i]);
	error = 1;
      }

  if (argc == 1 || error)
    {
      fprintf(stderr, "Usage: construct_maplist -pairs pairlist_file\n");
      fprintf(stderr, "                   -maps maps_prefix\n");
      fprintf(stderr, "                  [-min_correlation score]\n");
      fprintf(stderr, "                  [-max_correlation score]\n");
      fprintf(stderr, "                  [-min_distortion score]\n");
      fprintf(stderr, "                  [-max_distortion score]\n");
      fprintf(stderr, "                  [-weight weight]\n");
      exit(1);
    }

  /* check that at least minimal parameters were supplied */
  if (pairListName[0] == '\0')
    {
      fprintf(stderr, "-pairs must be specified\n");
      exit(1);
    }
  if (mapsName[0] == '\0')
    {
      fprintf(stderr, "-maps must be specified\n");
      exit(1);
    }

  /* read the pairs list */
  f = fopen(pairListName, "r");
  if (f == NULL)
    {
      fprintf(stderr, "Could not open pairs list %s\n", pairListName);
      exit(1);
    }
  while (fgets(line, LINE_LENGTH, f) != NULL)
    {
      nItems = sscanf(line, "%s %d %d %d %d %s %d %d %d %d %s",
		      imageName0, &minX0, &maxX0, &minY0, &maxY0,
		      imageName1, &minX1, &maxX1, &minY1, &maxY1,
		      pairName);
      if (nItems != 11)
	{
	  fprintf(stderr, "Malformed line in %s:\n%s\n", pairListName, line);
	  exit(1);
	}

      sprintf(scoreName, "%s%s.score", mapsName, pairName);
      sf = fopen(scoreName, "r");
      if (sf == NULL)
	{
	  fprintf(stderr, "Could not open file %s\n", scoreName);
	  exit(1);
	}
      if (fscanf(sf, "%lf%lf%lf%lf%lf", &energy, &correlation, &distortion,
		 &correspondence, &constraining) != 5)
	{
	  fprintf(stderr, "Could not read scores from file %s\n", scoreName);
	  exit(1);
	}
      fclose(sf);

      if (correlation < minCorrelation || correlation >= maxCorrelation ||
	  distortion < minDistortion || distortion >= maxDistortion)
	continue;
      
      printf("%s %s %s %f\n", imageName0, imageName1, pairName, weight);
    }
  fclose(f);
  return(0);
}
