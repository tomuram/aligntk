#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <stdarg.h>

#include "imio.h"

int
main (int argc, char **argv)
{
  FILE *f;
  int i;
  int n;
  int error;
  char map1Name[PATH_MAX];
  char map2Name[PATH_MAX];
  char outputName[PATH_MAX];
  MapElement *map1, *map2;
  int mLevel1, mLevel2;
  int mw1, mh1;
  int mw2, mh2;
  int mox1, moy1;
  int mox2, moy2;
  char imgn[PATH_MAX];
  char refn[PATH_MAX];
  char errorMsg[PATH_MAX+256];
  float mFactor1, mFactor2;
  int x, y;
  float xv1, yv1;
  float xv2, yv2;
  float rx00, ry00, rc00;
  float rx01, ry01, rc01;
  float rx10, ry10, rc10;
  float rx11, ry11, rc11;
  int ixv, iyv;
  float rrx, rry;
  double sd2;
  double rmsd;
  float rx1, ry1;
  float rx2, ry2;

  error = 0;
  map1Name[0] = '\0';
  map2Name[0] = '\0';
  outputName[0] = '\0';
  for (i = 1; i < argc; ++i)
    if (strcmp(argv[i], "-map1") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-map1 error\n");
	    break;
	  }
	strcpy(map1Name, argv[i]);
      }
    else if (strcmp(argv[i], "-map2") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    fprintf(stderr, "-map2 error\n");
	    break;
	  }
	strcpy(map2Name, argv[i]);
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
    else
      error = 1;

  if (argc == 1 || error)
    {
      fprintf(stderr, "Usage: rms_compare_maps -map1 map1.map\n");
      fprintf(stderr, "                        -map2 map2.map\n");
      fprintf(stderr, "                       [-output rms.out]\n");
      exit(1);
    }

  /* check that at least minimal parameters were supplied */
  if (map1Name[0] == '\0' || map2Name[0] == '\0')
    {
      fprintf(stderr, "-map1 and -map2 parameters must be specified.\n");
      exit(1);
    }

  if (!ReadMap(map1Name, &map1, &mLevel1, &mw1, &mh1,
	       &mox1, &moy1, imgn, refn, errorMsg))
    {
      fprintf(stderr, "Error reading map %s:\n  %s\n", map1Name, errorMsg);
      exit(1);
    }
  if (!ReadMap(map2Name, &map2, &mLevel2, &mw2, &mh2,
	       &mox2, &moy2, imgn, refn, errorMsg))
    {
      fprintf(stderr, "Error reading map %s:\n  %s\n", map2Name, errorMsg);
      exit(1);
    }

  /* go through the points of the first map */
  mFactor1 = 1 << mLevel1;
  mFactor2 = 1 << mLevel2;
  sd2 = 0.0;
  n = 0;
  for (y = 0; y < mh1; ++y)
    {
      yv1 =  (y + moy1) * mFactor1;
      yv2 = yv1 / mFactor2 - moy2;
      iyv = (int) floor(yv2);
      rry = yv2 - iyv;
      for (x = 0; x < mw1; ++x)
	{
	  if (map1[y*mw1+x].c <= 0.0)
	    continue;
	  rx1 = map1[y*mw1+x].x;
	  ry1 = map1[y*mw1+x].y;
	  ++n;
	  xv1 =  (x + mox1) * mFactor1;
	  xv2 = xv1 / mFactor2 - mox2;
	  ixv = (int) floor(xv2);
	  rrx = xv2 - ixv;
	  if (ixv < 0 || ixv >= mw2 ||
	      iyv < 0 || iyv >= mh2)
	    {
	      //	      printf("PENALTY: %f %f %d %d\n", xv2, yv2, ixv, iyv);
	      sd2 += mFactor1*mFactor1*(mw1+mh1)*(mw1+mh1);
	      continue;
	    }

	  while (ixv < 0)
	    {
	      ++ixv;
	      rrx -= 1.0;
	    }
	  while (ixv >= mw2-1)
	    {
	      --ixv;
	      rrx += 1.0;
	    }
	  while (iyv < 0)
	    {
	      ++iyv;
	      rry -= 1.0;
	    }
	  while (iyv >= mh2-1)
	    {
	      --iyv;
	      rry += 1.0;
	    }
	  rx00 = map2[iyv*mw2+ixv].x;
	  ry00 = map2[iyv*mw2+ixv].y;
	  rc00 = map2[iyv*mw2+ixv].c;
	  rx01 = map2[(iyv+1)*mw2+ixv].x;
	  ry01 = map2[(iyv+1)*mw2+ixv].y;
	  rc01 = map2[(iyv+1)*mw2+ixv].c;
	  rx10 = map2[iyv*mw2+ixv+1].x;
	  ry10 = map2[iyv*mw2+ixv+1].y;
	  rc10 = map2[iyv*mw2+ixv+1].c;
	  rx11 = map2[(iyv+1)*mw2+ixv+1].x;
	  ry11 = map2[(iyv+1)*mw2+ixv+1].y;
	  rc11 = map2[(iyv+1)*mw2+ixv+1].c;

	  if (rc00 == 0.0 ||
	      rc01 == 0.0 ||
	      rc10 == 0.0 ||
	      rc11 == 0.0)
	    {
	      //	      printf("PENALTY2: %f %f %d %d %f %f %f %f\n", xv2, yv2, ixv, iyv,
	      //		     rc00, rc01, rc10, rc11);
	      sd2 += mFactor1*mFactor1*(mw1+mh1)*(mw1+mh1);
	      continue;
	    }

	  rx2 = rx00 * (rrx - 1.0) * (rry - 1.0)
	    - rx10 * rrx * (rry - 1.0) 
	    - rx01 * (rrx - 1.0) * rry
	    + rx11 * rrx * rry;
	  ry2 = ry00 * (rrx - 1.0) * (rry - 1.0)
	    - ry10 * rrx * (rry - 1.0) 
	    - ry01 * (rrx - 1.0) * rry
	    + ry11 * rrx * rry;
	  sd2 += (rx2 - rx1) * (rx2 - rx1) + (ry2 - ry1) * (ry2 - ry1);
	}
    }
  if (n == 0)
    {
      fprintf(stderr, "No points of overlap found in maps\n");
      exit(1);
    }
  rmsd = sqrt(sd2 / n);

  if (outputName[0] != '\0')
    {
      f = fopen(outputName, "w");
      fprintf(f, "%lf\n", rmsd);
      fclose(f);
    }
  printf("RMS difference = %lf\n", rmsd);
  return(0);
}
