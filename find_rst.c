/* 
 * find_rst.c  - find image rotation, scale, and translation
 *
 *  Copyright (c) 2009-2015 Pittsburgh Supercomputing Center,
 *                          Carnegie Mellon University
 *
 *  This file is part of AlignTK.
 *
 *  AlignTK is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  AlignTK is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with AlignTK.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Acknowledgements:
 *     Development of this code was supported in part by
 *       NIH NCRR grant 5P41RR006009 and
 *       NIH NIGMS grant P41GM103712
 *
 *  HISTORY
 *    7/2009  Written by Greg Hood (ghood@psc.edu) based on the paper:
 *              W.Pan, K. Qin, and Y. Chen.
 *              An adaptable-multilayer fractional fourier
 *              transform approach for image registration.
 *              IEEE Trans. Patt. Anal. Mach. Intell., 31(3):400-413, 2009.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>
#include <unistd.h>
#include <limits.h>
#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <time.h>
#include <fftw3.h>
#include <mpi.h>

#include "imio.h"
#include "dt.h"
#include "par.h"

#define MAX_FRAC_FT_RES_LEVELS	4
#define LINE_LENGTH		255

typedef struct Context {
  char type;                        /* 't' for TIFF, 'p' for PGM */ 
  char imageBasename[PATH_MAX];
  char imageMaskBasename[PATH_MAX];
  char referenceBasename[PATH_MAX];
  char referenceMaskBasename[PATH_MAX];
  char outputBasename[PATH_MAX];
  char logBasename[PATH_MAX];
  int mapLevel;			    /* desired output map level */
  int outputImages;		    /* if 1, output all FT images */

  int maxRes;                       /* maximum FFT resolution; images will be reduced to fit into
				       half this size */
  float distortion;                 /* distortion that can be tolerated, expressed as a pixel
				       displacement */
  float margin;                     /* half-width of Blackman filter, expressed in pixels */
  float minFeatureSize;             /* expressed in pixels, used for RS */
  float maxFeatureSize;		    /* expressed in pixels, used for RS */
  float minTransFeatureSize;        /* expressed in pixels, used for T */
  float maxTransFeatureSize;        /* expressed in pixels, used for T */
  float minRotationalSeparation;    /* expressed in degrees */
  float minScaleSeparation;         /* expressed as a percentage increase */
  float minTranslationalSeparation; /* expressed as a percentage of image width or height */
  float fracRes[MAX_FRAC_FT_RES_LEVELS];    /* fractional FT resolutions (0.0 <= res < 1.0) */

  int maxRSCandidates;              /* number of candidates to try for rotation & scale */
  int maxCandidates;		    /* maximum number of final (rotation,scale,translation) candidates
				       to output */

  float minTheta, maxTheta;         /* expressed as the rotation in degrees to get the image
				       to match up with the reference */
  float minScale, maxScale;         /* expressed as the scale factor to be applied to the image
				       to match up with the reference */
  float minTX, maxTX;		    /* expressed as the translation in pixels to get the image
				       to match up with reference (after application of rotation
				       and scale) */
  float minTY, maxTY;		    /* expressed as the translation in pixels to get the image
				       to match up with reference (after application of rotation
				       and scale) */
  int update;			    /* if 1, just update the maps for
				       the frames that have changed */
  int partial;			    /* if 1, don't abort if input image files
				       are missing */
  int nWorkers;			    /* number of worker processes */
} Context;


typedef struct Pair {
  char *imageName;
  int imageMinX, imageMaxX, imageMinY, imageMaxY;
  char *refName;
  int refMinX, refMaxX, refMinY, refMaxY;
  char *pairName;
} Pair;

typedef struct Task {
  /* NOTE: any new fields added to this struct should
     be also added to PackTask and UnpackTask */ 
  Pair pair;
} Task;

typedef struct Result {
  Pair pair;
  int updated;
  float quality;
  float rotation;
  float scale;
  float tx;
  float ty;
  char *message;
} Result;

typedef struct PositionValue
{
  int x;
  int y;
  float val;
  float radius;
} PositionValue;

typedef struct Transformation
{
  /* transformatino to be applied to image to match it with reference */
  float rotation;  /* apply rotation first (positive is clockwise in degrees) */
  float scale;     /* apply scaling next about the center of the image */
  float tx, ty;    /* apply translation last (positive tx moves to the right,
		       positive ty moves down) */
  float radius;    /* what filtering radius yielded this transformation */
  float quality;   /* estimated quality of this transformation (higher is better) */
} Transformation;

/* GLOBAL VARIABLES FOR MASTER */
int nResults = 0;
Result* results = 0;
#define DIR_HASH_SIZE	8192
char *dirHash[DIR_HASH_SIZE];
char summaryName[PATH_MAX] = "";

/* GLOBAL VARIABLES FOR MASTER & WORKER */
Context c;
Task t;
Result r;
FILE *logFile = NULL;

/* GLOBAL VARIABLES FOR WORKER */
int last_iw[2] = {-1, -1};
int last_ih[2] = {-1, -1};
int last_n = -1;
int last_blk_w = -1;

float *dist[2];
float *img, *window, *mag, *lp, *cc, *windowed[2];
float *cos_t, *sin_t, *gaussian;
fftwf_complex *fft_img, *fft_orig[2];
fftwf_complex *Y, *Z, *W;
fftwf_complex *fft_Y, *fft_Z, *ifft_W;
fftwf_complex *prod, *fft_cc, *fft_lp, *flp[2];
fftwf_complex *rot;
fftwf_complex *ccfilter;
fftwf_plan plan_img, plan_Y, plan_Z, plan_W, plan_lp, plan_cc;
PositionValue *values;
unsigned char *marked;


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
int SortByName (const void *x, const void *y);
int SortByQuality (const void *x, const void *y);
int ParseRange (char *s, int *pos, int *minValue, int *maxValue);
int ParseValue (char *s, int *pos, int *value);
size_t CountBits (unsigned char *p, size_t n);
int CreateDirectories (char *fn);
void Error (char *fmt, ...);
void Log (char *fmt, ...);
int CompareValues (const void *v1, const void *v2);
void fft_shift (fftwf_complex *fft, int n);
void fft_expand (fftwf_complex *fft, int n);
void fft_compress (fftwf_complex *fft, int n);
char *GetTimestamp (char *timestamp, size_t size);

/* NOTES:

   CONVENTIONS:
      All image coordinates use a left-hand coordinate system with +x to the right,
      and +y going down.  The upper left pixel has coordinates (0,0).

      Positive rotation angles correspond to clockwise rotation; this allows the
      conventional formulas:
          x' = cos(theta) * x - sin(theta) * y
	  y' = sin(theta) * x + cos(theta) * y
      to be used for rotation by angle theta.

   IMAGES ARE ROTATED, SCALED, AND TRANSLATED TO FIT INTO A POWER-OF-TWO SIZE
     ARRAY FOR PROCESSING.  The transformation from internal coordinates (x',y')
     to image coordinates (x,y) is given by:

    x = (scale * (cos(theta) * (x' - n/2) - sin(theta) * (y' - n/2)) * rFactor + width/2
    y = (scale * (sin(theta) * (x' - n/2) + cos(theta) * (y' - n/2)) * rFactor + height/2

    and in the other direction:
	   
    x' = (1.0/scale) * (cos(theta) * (x - width/2) + sin(theta) * (y - height/2)) / rFactor + n/2
    y' = (1.0/scale) * (-sin(theta) * (x - width/2) + cos(theta) * (y - height/2)) / rFactor + n/2


    Now, assume scale0 = 1.0 and theta0 = 0.0, and we find that the translational difference between
    image0 and image1 in internal coordinates is (deltaX,deltaY).  Then,
      x1' = x0' + deltaX
      y1' = y0' + deltaY
    Substituting in equations above (canceling n/2 terms, and multiplying by rFactor):
      (1.0/scale1) * (cos(theta1) * (x1 - width1/2) + sin(theta1) * (y1 - height1/2)) = x0 - width0/2 + deltaX * rFactor
      (1.0/scale1) * (-sin(theta1) * (x1 - width1/2) + cos(theta1) * (y1 - height1/2)) = y0 - height0/2 + deltaY * rFactor
    x0,y0 can be expressed as a rotation, scaling, and translation of x1,y1:
      x0 = (1.0/scale1) * (cos(theta1) * (x1 - width1/2) + sin(theta1) * (y1 - height1/2)) + width0/2 - deltaX * rFactor
      y0 = (1.0/scale1) * (-sin(theta1) * (x1 - width1/2) + cos(theta1) * (y1 - height1/2)) + height0/2 - deltaY * rFactor

      x0 = (1.0/scale1) * (cos(-theta1) * x1 - sin(-theta1) * y1) +  (1.0/scale1) * (-cos(theta1) * width1/2 - sin(theta1) * height1/2) + width0/2 - deltaX * rFactor
      y0 = (1.0/scale1) * (sin(-theta1) * x1 + cos(-theta1) * y1) + (1.0/scale1) * (sin(theta1) * width1/2 - cos(theta1) * height1/2) + height0/2 - deltaY * rFactor

      i.e., transformation 1->0:
              rotation by -theta1
              scaling by 1.0/scale1
              translation by:
                ((1.0/scale1) * (-cos(theta1) * width1/2 - sin(theta1) * height1/2) + width0/2 - deltaX * rFactor,
                 (1.0/scale1) * (sin(theta1) * width1/2 - cos(theta1) * height1/2) + height0/2 - deltaY * rFactor)


    x1,y1 can also be expressed as a rotation, scaling, and translation of x0,y0:
      
      x0 - width0/2 + deltaX * rFactor = (1.0/scale1) * (cos(theta1) * (x1 - width1/2) + sin(theta1) * (y1 - height1/2))
      y0 - height0/2 + deltaY * rFactor = (1.0/scale1) * (-sin(theta1) * (x1 - width1/2) + cos(theta1) * (y1 - height1/2))

      x1 - width1/2 = scale1 * (cos(-theta1) * (x0 - width0/2 + deltaX * rFactor) + sin(-theta) * (y0 - height0/2 + deltaY * rFactor))
      y1 - height1/2 = scale1 * (-sin(-theta1) * (x0 - width0/2 + deltaX * rFactor) + cos(-theta1) * (y0 - height0/2 + deltaY * rFactor))

      x1 - width1/2 = scale1 * (cos(-theta1) * (x0 - width0/2 + deltaX * rFactor) + sin(-theta) * (y0 - height0/2 + deltaY * rFactor))
      y1 - height1/2 = scale1 * (-sin(-theta1) * (x0 - width0/2 + deltaX * rFactor) + cos(-theta1) * (y0 - height0/2 + deltaY * rFactor))

      x1 = scale1 * (cos(-theta1) * (x0 - width0/2 + deltaX * rFactor) + sin(-theta1) * (y0 - height0/2 + deltaY * rFactor)) + width1/2
      y1 = scale1 * (-sin(-theta1) * (x0 - width0/2 + deltaX * rFactor) + cos(-theta1) * (y0 - height0/2 + deltaY * rFactor)) + height1/2

      x1 = scale1 * (cos(theta1) * x0 - sin(theta) * y0) +
           scale1 * (cos(theta1) * (- width0/2 + deltaX * rFactor) - sin(theta1) * (- height0/2 + deltaY * rFactor)) + width1/2
      y1 = scale1 * (sin(theta1) * x0 + cos(theta1) * y0) +
           scale1 * (sin(theta1) * (- width0/2 + deltaX * rFactor) + cos(theta1) * (- height0/2 + deltaY * rFactor)) + height1/2

      i.e., transformation 0->1:
              rotation by theta1
              scaling by scale1
	      translation by:
                (scale1 * (cos(theta1) * (- width0/2 + deltaX * rFactor) - sin(theta1) * (- height0/2 + deltaY * rFactor)) + width1/2,
                 scale1 * (sin(theta1) * (- width0/2 + deltaX * rFactor) + cos(theta1) * (- height0/2 + deltaY * rFactor)) + height1/2)

*/


int
main (int argc, char **argv, char **envp)
{
  par_process(argc, argv, envp,
              (void (*)()) MasterTask, MasterResult,
              WorkerContext, WorkerTask, NULL,
              PackContext, UnpackContext,
              PackTask, UnpackTask,
              PackResult, UnpackResult);
  return(0);
}


/* MASTER PROCEDURES */

void
MasterTask (int argc,
            char **argv,
            char **envp)
{
  int i;
  int error;
  char pairsFile[PATH_MAX];
  int nPairs;
  Pair *pairs;
  FILE *f;
  char imgn[PATH_MAX], refn[PATH_MAX], pairn[PATH_MAX];
  int imgMinX, imgMaxX, imgMinY, imgMaxY;
  int refMinX, refMaxX, refMinY, refMaxY;
  char fn[PATH_MAX];
  char outputDirName[PATH_MAX];
  int pn;
  char line[LINE_LENGTH+1];

  error = 0;
  c.type = 'p';
  c.imageBasename[0] = '\0';
  c.imageMaskBasename[0] = '\0';
  c.referenceBasename[0] = '\0';
  c.referenceMaskBasename[0] = '\0';
  c.outputBasename[0] = '\0';
  c.logBasename[0] = '\0';
  c.mapLevel = -1;                            // output the coarsest map possible
  c.outputImages = 0;

  c.maxRes = -1;                              // use FFTs that accommodate the full-res images
  c.distortion = -1.0;                        // let the distortion be 1% of max(image width,
                                              //                                 image height)
  c.margin = -1.0;                            // let the feature size determine the margins
  c.minFeatureSize = -1.0;                    // use features larger than 1 pixels in size
  c.maxFeatureSize = -1.0;		      // use features up to 10% of min(image width,
                                              //                               image height)
  c.minTransFeatureSize = -1.0;               // use features larger than 1 pixels in size
  c.maxTransFeatureSize = -1.0;               // use features as large as image
  c.minRotationalSeparation = 2.0;            // 2 degrees
  c.minScaleSeparation = 2.0;                 // 2 percent
  c.minTranslationalSeparation = 1.0;         // 1% shift
  for (i = 0; i < MAX_FRAC_FT_RES_LEVELS; ++i)
    c.fracRes[i] = 0.0;
  c.fracRes[0] = 1.0;
  c.fracRes[1] = 0.6;
  c.fracRes[2] = 0.2;
  c.maxRSCandidates = 4;
  c.maxCandidates = 1;                        // output only the best final RST candidate

  c.minTheta = 0.0;                           // consider all angles
  c.maxTheta = 360.0;
  c.minScale = 0.5;                           // consider scales between 0.5 and 2.0
  c.maxScale = 2.0;
  c.minTX = -100.0;                           // consider translations up to the full size of
  c.maxTX = 100.0;                            //   the images
  c.minTY = -100.0;
  c.maxTY = 100.0;
  c.update = 0;
  c.partial = 0;
  c.nWorkers = par_workers();

  r.pair.imageName = NULL;
  r.pair.refName = NULL;
  r.pair.pairName = NULL;
  r.message = (char *) malloc(PATH_MAX + 1024);

  pairsFile[0] = '\0';
  nPairs = 0;
  pairs = 0;
  for (i = 1; i < argc; ++i)
    if (strcmp(argv[i], "-images") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.imageBasename, argv[i]);
      }
    else if (strcmp(argv[i], "-mask") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.imageMaskBasename, argv[i]);
      }
    else if (strcmp(argv[i], "-reference") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.referenceBasename, argv[i]);
      }
    else if (strcmp(argv[i], "-reference_mask") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.referenceMaskBasename, argv[i]);
      }
    else if (strcmp(argv[i], "-logs") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.logBasename, argv[i]);
      }
    else if (strcmp(argv[i], "-pairs") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(pairsFile, argv[i]);
      }
    else if (strcmp(argv[i], "-output") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.outputBasename, argv[i]);
      }
    else if (strcmp(argv[i], "-summary") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(summaryName, argv[i]);
      }
    else if (strcmp(argv[i], "-output_images") == 0)
      {
	c.outputImages = 1;
      }
    else if (strcmp(argv[i], "-margin") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%f", &c.margin) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-max_res") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d", &c.maxRes) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-candidates") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d", &c.maxCandidates) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-rs_candidates") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d", &c.maxRSCandidates) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-distortion") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%f", &c.distortion) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-feature") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	if (sscanf(argv[i], "%f-%f", &c.minFeatureSize, &c.maxFeatureSize) != 2)
	  {
	    if (sscanf(argv[i], "%f", &c.minFeatureSize) != 1)
	      {
		error = 1;
		break;
	      }
	  }
      }
    else if (strcmp(argv[i], "-trans_feature") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	if (sscanf(argv[i], "%f-%f", &c.minTransFeatureSize, &c.maxTransFeatureSize) != 2)
	  {
	    if (sscanf(argv[i], "%f", &c.minTransFeatureSize) != 1)
	      {
		error = 1;
		break;
	      }
	  }
      }
    else if (strcmp(argv[i], "-rotation") == 0)
      { 
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	if (sscanf(argv[i], "%f-%f", &c.minTheta, &c.maxTheta) != 2)
	  {
	    if (sscanf(argv[i], "%f", &c.minTheta) != 1)
	      {
		error = 1;
		break;
	      }
	    c.maxTheta = c.minTheta;
	  }
      }
    else if (strcmp(argv[i], "-scale") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	if (sscanf(argv[i], "%f-%f", &c.minScale, &c.maxScale) != 2)
	  {
	    if (sscanf(argv[i], "%f", &c.minScale) != 1)
	      {
		error = 1;
		break;
	      }
	    c.maxScale = c.minScale;
	  }
      }
    else if (strcmp(argv[i], "-tx") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	if (sscanf(argv[i], "%f-%f", &c.minTX, &c.maxTX) != 2)
	  {
	    if (sscanf(argv[i], "%f", &c.minTX) != 1)
	      {
		error = 1;
		break;
	      }
	    c.maxTX = c.minTX;
	  }
      }
    else if (strcmp(argv[i], "-ty") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	if (sscanf(argv[i], "%f-%f", &c.minTY, &c.maxTY) != 2)
	  {
	    if (sscanf(argv[i], "%f", &c.minTY) != 1)
	      {
		error = 1;
		break;
	      }
	    c.maxTY = c.minTY;
	  }
      }
    else if (strcmp(argv[i], "-frac_res") == 0)
      {
	if (sscanf(argv[i], "%f,%f,%f", &c.fracRes[1], &c.fracRes[2], &c.fracRes[3]) != 3)
	  {
	    if (sscanf(argv[i], "%f,%f", &c.fracRes[1], &c.fracRes[2]) != 2)
	      {
		if (sscanf(argv[i], "%f", &c.fracRes[1]) != 1)
		  {
		    error = 1;
		    break;
		  }
	      }
	  }
      }
    else if (strcmp(argv[i], "-tif") == 0)
      c.type = 't';
    else if (strcmp(argv[i], "-update") == 0)
      c.update = 1;
    else if (strcmp(argv[i], "-partial") == 0)
      c.partial = 1;
  
    else
      {
	error = 1;
	break;
      }
 
  if (error)
    {
      if (i < argc)
	fprintf(stderr, "Invalid option: %s\n", argv[i]);
      else
	fprintf(stderr, "Incomplete option: %s\n", argv[i-1]);
      fprintf(stderr, "\n");

      fprintf(stderr, "Usage: find_rst -images <file_prefix>\n");
      fprintf(stderr, "            [-reference <reference_prefix>]\n");
      fprintf(stderr, "            -pairs <pair_file>\n");
      fprintf(stderr, "            -output <output_prefix>]\n");

      fprintf(stderr, "            [-margin <margin_size_in_pixels>]\n");
      fprintf(stderr, "            [-max_res <maximum_pixel_resolution>]\n");
      fprintf(stderr, "            [-min_scale <min_scale>]\n");
      fprintf(stderr, "            [-max_scale <max_scale>]\n");
      fprintf(stderr, "            [-min_rotation <min_rotation>]\n");
      fprintf(stderr, "            [-max_rotation <max_rotation>]\n");

      fprintf(stderr, "            [-min_rotation_separation <degrees>]\n");
      fprintf(stderr, "            [-min_scale_separation <percent>]\n");
      fprintf(stderr, "            [-update]\n");
      fprintf(stderr, "            [-partial]\n");
      fprintf(stderr, "            [-logs <log_file_prefix>]\n");
      exit(1);
    }

  /* check that at least minimal parameters were supplied */
  if (c.imageBasename[0] == '\0' || c.outputBasename[0] == '\0' ||
      pairsFile[0] == '\0')
    Error("-input, -output, and -pairs parameters must be specified.\n");
  f = fopen(pairsFile, "r");
  if (f == NULL)
    Error("Could not open pairs file %s\n", pairsFile);
  while (fgets(line, LINE_LENGTH, f) != NULL)
    {
      if (line[0] == '\0' || line[0] == '#')
	continue;
      if (sscanf(line, "%s %d %d %d %d %s %d %d %d %d %s",
		 imgn, &imgMinX, &imgMaxX, &imgMinY, &imgMaxY,
		 refn, &refMinX, &refMaxX, &refMinY, &refMaxY,
		 pairn) != 11)
	{
	  if (sscanf(line, "%s %s %s", imgn, refn, pairn) != 3)
	    Error("Invalid line in pairs file %s:\n%s\n", pairsFile, line);
	  imgMinX = -1;
	  imgMaxX = -1;
	  imgMinY = -1;
	  imgMaxY = -1;
	  refMinX = -1;
	  refMaxX = -1;
	  refMinY = -1;
	  refMaxY = -1;
	}
      if ((nPairs & 1023) == 0)
	pairs = (Pair *) realloc(pairs, (nPairs + 1024) * sizeof(Pair));
      pairs[nPairs].imageName = (char *) malloc(strlen(imgn)+1);
      strcpy(pairs[nPairs].imageName, imgn);
      pairs[nPairs].imageMinX = imgMinX;
      pairs[nPairs].imageMaxX = imgMaxX;
      pairs[nPairs].imageMinY = imgMinY;
      pairs[nPairs].imageMaxY = imgMaxY;
      pairs[nPairs].refName = (char *) malloc(strlen(refn)+1);
      strcpy(pairs[nPairs].refName, refn);
      pairs[nPairs].refMinX = refMinX;
      pairs[nPairs].refMaxX = refMaxX;
      pairs[nPairs].refMinY = refMinY;
      pairs[nPairs].refMaxY = refMaxY;
      pairs[nPairs].pairName = (char *) malloc(strlen(pairn)+1);
      strcpy(pairs[nPairs].pairName, pairn);
      ++nPairs;
    }
  fclose(f);

  /* check that output directories are writeable */ 
  Log("MASTER checking output directories\n");
  sprintf(fn, "%sTEST.map", c.outputBasename);
  f = fopen(fn, "w");
  if (f == NULL)
    {
      for (i = strlen(c.outputBasename)-1; i >= 0 && c.outputBasename[i] != '/'; --i) ;
      if (i != 0)
	strncpy(outputDirName, c.outputBasename, i+1);
      else
	outputDirName[0] = '.';
      outputDirName[i+1] = '\0';
      Error("Could not open test output file %s --\n       does directory %s exist and is it writeable?\n", fn, outputDirName);
    }
  fclose(f);
  unlink(fn);

 Log("MASTER setting context\n");

  par_set_context();

  Log("MASTER set context\n");

  /* for all slices */
  printf("Processing image pairs: ");
  fflush(stdout);
  t.pair.imageMinX = t.pair.imageMaxX = t.pair.imageMinY = t.pair.imageMaxY = -1;
  t.pair.refMinX = t.pair.refMaxX = t.pair.refMinY = t.pair.refMaxY = -1;
  Log("nPairs = %d\n", nPairs);
  for (pn = 0; pn < nPairs; ++pn)
    {
      memcpy(&(t.pair), &(pairs[pn]), sizeof(Pair));

      // make sure that output directories exist
      sprintf(fn, "%s%s.rst", c.outputBasename, t.pair.pairName);
      if (!CreateDirectories(fn))
	continue;

      Log("Delegating pair %d\n", pn);
      par_delegate_task();
    }
  par_finish();

  if (summaryName[0] != '\0')
    {
      f = fopen(summaryName, "w");
      if (f == NULL)
	Error("Could not open summary output file %s\n", summaryName);

      qsort(results, nResults, sizeof(Result), SortByName);
      fprintf(f, "Sorted by slice:\n");
      fprintf(f, "         SOURCE  TARGET   QUALITY     ROTATION    SCALE       TX          TY\n");
      for (i = 0; i < nResults; ++i)
	fprintf(f, "%s  %s  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f\n",
		results[i].pair.imageName,
		results[i].pair.refName,
		results[i].quality,
		results[i].rotation,
		results[i].scale,
		results[i].tx,
		results[i].ty);
      fprintf(f, "\n\nSorted by quality:\n");
      fprintf(f, "         SOURCE  TARGET   QUALITY     ROTATION    SCALE       TX          TY\n");
      qsort(results, nResults, sizeof(Result), SortByQuality);
      for (i = 0; i < nResults; ++i)
	fprintf(f, "%s  %s  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f\n",
		results[i].pair.imageName,
		results[i].pair.refName,
		results[i].quality,
		results[i].rotation,
		results[i].scale,
		results[i].tx,
		results[i].ty);
      fclose(f);
    }

  printf(" %d\nAll image pairs completed.\n", nResults);
}

void
MasterResult ()
{
  if (r.message[0] != '\0')
    Error("\nThe following error was encountered by one of the worker processes:\n");

  results = (Result*) realloc(results, (nResults + 1) * sizeof(Result));
  results[nResults] = r;
  results[nResults].message = NULL;
  r.pair.imageName = NULL;
  r.pair.refName = NULL;
  r.pair.pairName = NULL;

  if ((nResults % 50) == 0 && nResults != 0)
    printf(" %d \n                   ", nResults);
  printf(".");
  fflush(stdout);
  ++nResults;
}


/* WORKER PROCEDURES */
void
WorkerContext ()
{
  Log("WorkerContext called\n");
  r.pair.imageName = NULL;
  r.pair.refName = NULL;
  r.pair.pairName = NULL;
  r.message = (char *) malloc(PATH_MAX + 1024);
  r.message[0] = '\0';
  t.pair.imageName = NULL;
  t.pair.refName = NULL;
  t.pair.pairName = NULL;
}

void
WorkerTask ()
{
  unsigned char *image_in[2];
  unsigned char *image_mask_in[2];
  int iw[2], ih[2];
  int imw[2], imh[2];
  int eiw[2], eih[2];
  int n, n2, n_over_2, n_over_2_plus_1, n_times_2, n2_partial;
  char msg[2*PATH_MAX];
  int i, j, k;
  int x, y;
  float xv, yv;
  float v;
  float theta, rho;
  int it, ir;
  float logrhobase, logrhooffset;
  float rhobase;
  float vi, vr;
  float max_v;
  float v_max, v_min;
  int max_x, max_y;
  float p_r, p_i;
  float d, dst;
  int count;
  int offset_x, offset_y;
  float alpha;
  float a0, a1, a2;
  float wy;
  int nRes;
  float resMin[MAX_FRAC_FT_RES_LEVELS];
  int res;
  int imi;
  float cos_theta, sin_theta;
  int ixv, iyv;
  float rrx, rry;
  float rv00, rv01, rv10, rv11;
  float hpfx, hpfy;
  float filter;
  char fn[PATH_MAX];
  FILE *f;
  float xvp, yvp;
  int nCandidates, nCandidatesAdded;
  int dx, dy;
  int nit, nir;
  PositionValue *candidates;
  float radius, sigma;
  float minDegreeRotationalSeparation, minPercentageScaleSeparation;
  int blk_w;
  float rhoFactor;
  int delta_it, delta_ir;
  int target_it, target_ir;
  int max_delta_it, max_delta_ir;
  float minITSeparation, minIRSeparation;
  char imageName[PATH_MAX], imageMaskName[PATH_MAX];
  char refName[PATH_MAX], refMaskName[PATH_MAX];
  char outputName[PATH_MAX], mapName[PATH_MAX], scoreName[PATH_MAX];
  double sum;
  float mean[2];
  int th;
  Transformation *finalCandidates;
  int nFinalCandidates;
  float baseRotation;
  float rotation, scale, tx, ty, quality;
  float otx, oty;
  int ix, iy;
  int nix, niy;
  int smaller, larger;
  float ct, st;
  int rFactor, sFactor;
  MapElement *map;
  int md;
  int offset2_x, offset2_y;
  int mMinX, mMaxX, mMinY, mMaxY;
  int mw, mh;
  int iMinX, iMaxX, iMinY, iMaxY;
  int level;
  int maxDim;
  float maxRadius;
  int ccFilterMin, ccFilterMax;
  unsigned char *mask;
  struct stat sb;
  int computeMap;
  double outputTime;
  float deltaX, deltaY;

  Log("Worker received task %s -> %s\n", t.pair.imageName, t.pair.refName);
  /* construct filenames */
  //  sprintf(imageName, "%s%s.%s", c.imageBasename, t.pair.imageName,
  //	  c.type == 't' ? "tif" : "pgm");
  sprintf(imageName, "%s%s", c.imageBasename, t.pair.imageName);
  if (c.imageMaskBasename[0] != '\0')
    sprintf(imageMaskName, "%s%s", c.imageMaskBasename, t.pair.imageName);
  else
    imageMaskName[0] = '\0';
  if (c.referenceBasename[0] != '\0')
    sprintf(refName, "%s%s", c.referenceBasename, t.pair.refName);
  else
    sprintf(refName, "%s%s", c.imageBasename, t.pair.refName);
  if (c.referenceMaskBasename[0] != '\0')
    sprintf(refMaskName, "%s%s", c.referenceMaskBasename, t.pair.refName);
  else if (c.referenceBasename[0] == '\0' && c.imageMaskBasename[0] != '\0')
    sprintf(refMaskName, "%s%s", c.imageMaskBasename, t.pair.refName);
  else
    refMaskName[0] = '\0';

  if (c.outputBasename[0] != '\0')
    sprintf(outputName, "%s%s", c.outputBasename, t.pair.pairName);
  else
    outputName[0] = '\0';

  /* check if we can skip this task because of the -update option */
  computeMap = 0;
  if (!c.update)
    computeMap = 1;
  /* check if output map file already exists */
  sprintf(scoreName, "%s.score", outputName);
  if (stat(scoreName, &sb) == 0)
    outputTime = (double) sb.st_mtime;
  else
    computeMap = 1;
  if (!computeMap)
    {
      if (stat(imageName, &sb) == 0 &&
	  (double) (sb.st_mtime) > outputTime)
	computeMap = 1;
    }
  if (!computeMap && imageMaskName[0] != '\0')
    {
      if (stat(imageMaskName, &sb) == 0 &&
	  (double) (sb.st_mtime) > outputTime)
	computeMap = 1;
    }
  if (!computeMap)
    {
      if (stat(refName, &sb) == 0 &&
	  (double) (sb.st_mtime) > outputTime)
	computeMap = 1;
    }
  if (!computeMap && refMaskName[0] != '\0')
    {
      if (stat(refMaskName, &sb) == 0 &&
	  (double) (sb.st_mtime) > outputTime)
	computeMap = 1;
    }
  if (!computeMap)
    {
      sprintf(scoreName, "%s.score", outputName);
      f = fopen(scoreName, "r");
      if (f == NULL ||
	  fscanf(f, "%f%f%f%f%f",
		 &(r.quality),
		 &(r.rotation),
		 &(r.scale),
		 &(r.tx),
		 &(r.ty)) != 5)
	computeMap = 1;
      else
	r.rotation = r.rotation * M_PI / 180.0;
      if (f != NULL)
	fclose(f);
    }
  if (!computeMap)
    {
      Log("WORKER TASK skipping %s to %s since up-to-date\n",
	  imageName, refName);
      memcpy(&(r.pair), &(t.pair), sizeof(Pair));
      r.updated = 0;
      r.message[0] = '\0';
      return;
    }
  Log("WORKER TASK opening files for %s to %f\n",
      imageName, refName);

  if (c.partial)
    {
      computeMap = 1;
      if ((f = fopen(imageName, "r")) != NULL)
	fclose(f);
      else
	computeMap = 0;
      if ((f = fopen(imageMaskName, "r")) != NULL)
	fclose(f);
      else
	computeMap = 0;
      if ((f = fopen(refName, "r")) != NULL)
	fclose(f);
      else
	computeMap = 0;
      if ((f = fopen(refMaskName, "r")) != NULL)
	fclose(f);
      else
	computeMap = 0;
      if (!computeMap)
	{
	  Log("WORKER TASK skipping %s to %s since input files missing\n",
	      imageName, refName);
	  memcpy(&(r.pair), &(t.pair), sizeof(Pair));
	  r.updated = 0;
	  r.message[0] = '\0';
	  return;
	}
    }

  Log("\nRegistering image %s against reference %s\n", imageName, refName);

  if (!ReadImage(imageName, &image_in[0],
		 &iw[0], &ih[0],
		 t.pair.imageMinX, t.pair.imageMaxX,
		 t.pair.imageMinY, t.pair.imageMaxY, msg))
    Error("%s\n", msg);
  Log("Image %s read in.\n", imageName);
  if (!ReadImage(refName, &image_in[1],
		 &iw[1], &ih[1],
		 t.pair.refMinX, t.pair.refMaxX,
		 t.pair.refMinY, t.pair.refMaxY, msg))
    Error("%s\n", msg);
  Log("Image %s read in.\n", refName);

  image_mask_in[0] = NULL;
  image_mask_in[1] = NULL;
  if (imageMaskName[0] != '\0')
    {
      if (!ReadBitmap(imageMaskName, &image_mask_in[0],
		      &imw[0], &imh[0],
		      t.pair.imageMinX, t.pair.imageMaxX,
		      t.pair.imageMinY, t.pair.imageMaxY, msg))
	Error("%s\n", msg);
      Log("Mask %s read in.\n", imageMaskName);
      if (imw[0] != iw[0] || imh[0] != ih[0])
	Error("image mask %s has different dimensions (%d x %d) than image %s (%d x %d)\n",
	      imageMaskName, imw[0], imh[0],
	      imageName, iw[0], ih[0]);
    }
  if (refMaskName[0] != '\0')
    {
      if (!ReadBitmap(refMaskName, &image_mask_in[1],
		      &imw[1], &imh[1],
		      t.pair.refMinX, t.pair.refMaxX,
		      t.pair.refMinY, t.pair.refMaxY, msg))
	Error("%s\n", msg);
      Log("Mask %s read in.\n", refMaskName);
      if (imw[1] != iw[1] || imh[1] != ih[1])
	Error("reference mask %s has different dimensions than reference %s\n",
	      refMaskName, refName);
    }

  /* determine resolutions */
  n = 1;
  while (n < iw[0] + iw[1] || n < ih[0] + ih[1])
    n <<= 1;
  //  printf("full res: n = %d\n", n);
  rFactor = 1;
  if (c.maxRes > 0)
    while (n > c.maxRes)
      {
	++rFactor;
	n = 1;
	while (n < (iw[0] + iw[1]) / rFactor ||
	       n < (ih[0] + ih[1]) / rFactor)
	  n <<= 1;
      }

  //  printf("effective res: n = %d\n", n);
  maxDim = -1;
  for (imi = 0; imi < 2; ++imi)
    {
      if (iw[imi] > maxDim)
	maxDim = iw[imi];
      if (ih[imi] > maxDim)
	maxDim = ih[imi];

      if (iw[imi] != last_iw[imi] || ih[imi] != last_ih[imi])
	{
	  if (last_iw[imi] >= 0)
	    {
	      fftwf_free(windowed[imi]);
	      fftwf_free(dist[imi]);
	    }
	  windowed[imi] = (float*) fftwf_malloc(ih[imi] * iw[imi] * sizeof(float));
	  dist[imi] = (float *) fftwf_malloc(iw[imi]*ih[imi]*sizeof(float));

	  last_iw[imi] = iw[imi];
	  last_ih[imi] = ih[imi];
	}
    }

  n2 = n * n;
  n_over_2 = n / 2;
  n_over_2_plus_1 = n_over_2 + 1;
  n2_partial = n_over_2_plus_1 * n;
  n_times_2 = n * 2;
  if (c.maxFeatureSize > 0.0)
    logrhooffset = log(c.maxFeatureSize / rFactor);
  else if (c.margin > 0.0)
    logrhooffset = log(c.margin / rFactor / 2.0);
  else if (iw[0] < ih[0])
    logrhooffset = log(0.1 * iw[0] / rFactor);
  else
    logrhooffset = log(0.1 * ih[0] / rFactor);
  if (c.minFeatureSize > 0.0 &&
      n / (2.0 * c.minFeatureSize / rFactor) < n/2)
    logrhobase = (log(n / (2.0 * c.minFeatureSize / rFactor)) - logrhooffset) / n;
  else
    logrhobase = (log(n/2.0) - logrhooffset) / n;
  minITSeparation = n * c.minRotationalSeparation / 180.0;
  minIRSeparation = log(1.0 + 0.01 * c.minScaleSeparation) / logrhobase;
  candidates = (PositionValue*) malloc(c.maxRSCandidates * sizeof(PositionValue));
  finalCandidates = (Transformation*) malloc(c.maxCandidates * sizeof(Transformation));

  if (c.distortion >= 0.0)
    radius = 0.01 * c.distortion * maxDim / rFactor / 2.0;
  else
    radius = 0.01 * maxDim / rFactor / 2.0;

  if (c.margin >= 0.0)
    blk_w = (int) floor(2.0*c.margin);
  else if (c.maxFeatureSize > 0.0)
    blk_w = (int) floor(4.0*c.maxFeatureSize);
  else
    {
      blk_w = iw[0];
      if (ih[0] < blk_w)
	blk_w = ih[0];
      if (iw[1] < blk_w)
	blk_w = iw[1];
      if (ih[1] < blk_w)
	blk_w = ih[1];
      blk_w = blk_w / 2;
    }

  if (c.minTransFeatureSize > 0.0)
    ccFilterMax = (int) (n / (2.0 * c.minTransFeatureSize / rFactor) + 0.5);
  else
    ccFilterMax = n_over_2;
  if (c.maxTransFeatureSize > 0.0)
    ccFilterMin = (int) (n / (2.0 * c.maxTransFeatureSize / rFactor) + 0.5);
  else
    ccFilterMin = 0;

  nRes = MAX_FRAC_FT_RES_LEVELS;
  for (i = 0; i < MAX_FRAC_FT_RES_LEVELS; ++i)
    if (c.fracRes[i] <= 0.0)
      {
	nRes = i;
	break;
      }

  Log("Using rFactor=%d radius=%f blk_w=%d ccFilterMin=%d ccFilterMax=%d\n",
      rFactor, radius, blk_w, ccFilterMin, ccFilterMax);

#if 0
  // old code
  logrhooffset = log(0.03 * n);
  rhoFactor = 0.2;
  logrhobase = (log(rhoFactor*n-1.0)  - logrhooffset) / n;
#endif

  // allocate data structures
  if (blk_w != last_blk_w)
    {
      if (last_blk_w >= 0 && window != NULL)
	fftwf_free(window);
      if (blk_w >= 2)
	window = (float *) fftwf_malloc(blk_w/2 * sizeof(float));
      else
	window = NULL;
      last_blk_w = blk_w;
    }
  if (n != last_n)
    {
      if (last_n > 0)
	{
	  fftwf_free(img);
	  fftwf_free(mag);
	  fftwf_free(lp);
	  fftwf_free(cc);
	  fftwf_free(cos_t);
	  fftwf_free(sin_t);
	  fftwf_free(fft_img);
	  fftwf_free(fft_orig[0]);
	  fftwf_free(fft_orig[1]);
	  fftwf_free(Y);
	  fftwf_free(Z);
	  fftwf_free(W);
	  fftwf_free(fft_Y);
	  fftwf_free(fft_Z);
	  fftwf_free(ifft_W);
	  fftwf_free(prod);
	  fftwf_free(gaussian);
	  fftwf_free(fft_cc);
	  fftwf_free(fft_lp);
	  fftwf_free(flp[0]);
	  fftwf_free(flp[1]);
	  fftwf_free(rot);
	  fftwf_free(values);
	  fftwf_free(marked);
	  fftwf_free(ccfilter);

	  fftwf_destroy_plan(plan_img);
	  fftwf_destroy_plan(plan_Y);
	  fftwf_destroy_plan(plan_Z);
	  fftwf_destroy_plan(plan_W);
	  fftwf_destroy_plan(plan_lp);
	  fftwf_destroy_plan(plan_cc);
	}
      img = (float*) fftwf_malloc(n2 * sizeof(float));
      mag = (float*) fftwf_malloc(nRes * n2 * sizeof(float));
      lp = (float *) fftwf_malloc(n2 * sizeof(float));
      cc = (float*) fftwf_malloc(n2 * sizeof(float));
      cos_t = (float *) fftwf_malloc(n * sizeof(float));
      sin_t = (float *) fftwf_malloc(n * sizeof(float));
      fft_img = (fftwf_complex*) fftwf_malloc(n2 * sizeof(fftwf_complex));
      fft_orig[0] = (fftwf_complex*) fftwf_malloc(n2_partial * sizeof(fftwf_complex));
      fft_orig[1] = (fftwf_complex*) fftwf_malloc(n2_partial * sizeof(fftwf_complex));
      Y = (fftwf_complex*) fftwf_malloc(2*n*sizeof(fftwf_complex));
      Z = (fftwf_complex*) fftwf_malloc(2*n*sizeof(fftwf_complex));
      W = (fftwf_complex*) fftwf_malloc(2*n*sizeof(fftwf_complex));
      fft_Y = (fftwf_complex*) fftwf_malloc(2*n*sizeof(fftwf_complex));
      fft_Z = (fftwf_complex*) fftwf_malloc(2*n*sizeof(fftwf_complex));
      ifft_W = (fftwf_complex*) fftwf_malloc(2*n*sizeof(fftwf_complex));
      prod = (fftwf_complex*) fftwf_malloc(n2_partial * sizeof(fftwf_complex));
      gaussian = (float *) fftwf_malloc(n * sizeof(float));
      fft_cc = (fftwf_complex*) fftwf_malloc(n2_partial * sizeof(fftwf_complex));
      fft_lp = (fftwf_complex*) fftwf_malloc(n2_partial * sizeof(fftwf_complex));
      flp[0] = (fftwf_complex*) fftwf_malloc(n2_partial * sizeof(fftwf_complex));
      flp[1] = (fftwf_complex*) fftwf_malloc(n2_partial * sizeof(fftwf_complex));
      rot = (fftwf_complex*) fftwf_malloc(2 * n * sizeof(fftwf_complex));
      values = (PositionValue*) fftwf_malloc(n2 *sizeof(PositionValue));
      marked = (unsigned char *) fftwf_malloc(n2 * sizeof(unsigned char));
      ccfilter = (fftwf_complex*) fftwf_malloc(n2 * sizeof(fftwf_complex));

      plan_img = fftwf_plan_dft_r2c_2d(n, n, img, fft_img, FFTW_ESTIMATE);
      plan_Y = fftwf_plan_dft_1d(2*n, Y, fft_Y, FFTW_FORWARD, FFTW_ESTIMATE);
      plan_Z = fftwf_plan_dft_1d(2*n, Z, fft_Z, FFTW_FORWARD, FFTW_ESTIMATE);
      plan_W = fftwf_plan_dft_1d(2*n, W, ifft_W, FFTW_BACKWARD, FFTW_ESTIMATE);
      plan_lp = fftwf_plan_dft_r2c_2d(n, n, lp, fft_lp, FFTW_ESTIMATE);
      plan_cc = fftwf_plan_dft_c2r_2d(n, n, fft_cc, cc, FFTW_ESTIMATE);

      last_n = n;
    }

  // construct Blackman window
  alpha = 0.16;
  a0 = 0.5 * (1.0 - alpha);
  a1 = 0.5;
  a2 = 0.5 * alpha;
  for (x = 0; x < blk_w/2; ++x)
    window[x] = a0 - a1 * cos(2.0 * M_PI * x / (blk_w - 1)) +
      a2 * cos(4.0 * M_PI * x / (blk_w - 1));

  for (imi = 0; imi < 2; ++imi)
    {
      eiw[imi] = iw[imi] / rFactor;
      eih[imi] = ih[imi] / rFactor;
      offset_x = (n - eiw[imi]) >> 1;
      offset_y = (n - eih[imi]) >> 1;
      //      printf("iw = %d ih = %d\n", iw[imi], ih[imi]);
      //      printf("offset_x = %d offset_y = %d\n", offset_x, offset_y);

      sum = 0.0;
      for (y = 0; y < ih[imi]; ++y)
	for (x = 0; x < iw[imi]; ++x)
	  sum += image_in[imi][y*iw[imi] + x];
      mean[imi] = sum / ih[imi] / iw[imi];
      //    printf("mean[%d] = %f\n", imi, mean[imi]);
      //      mean[imi] = 0.0;

      sFactor = rFactor;
      while (sFactor < 4)
	sFactor += rFactor;
      //      printf("rFactor = %d  sFactor = %d\n", rFactor, sFactor);

      if (blk_w >= 2)
	{
	  // compute distance from mask
	  if (image_mask_in[imi] != NULL)
	    {
	      count = ih[imi] * ((iw[imi] + 7) >> 3);
	      mask = image_mask_in[imi];
	      for (i = 0; i < count; ++i)
		mask[i] ^= 0xff;
	      computeDistance(EUCLIDEAN_DISTANCE,
			      iw[imi], ih[imi],
			      mask,
			      dist[imi]);
	    }
	  else
	    for (y = 0; y < ih[imi]; ++y)
	      for (x = 0; x < iw[imi]; ++x)
		dist[imi][y*iw[imi] + x] = 1000000000.0;

	  // take min with distance from edge, and
	  // apply a Blackman window to smooth out edges
	  for (y = 0; y < ih[imi]; ++y)
	    for (x = 0; x < iw[imi]; ++x)
	      {
		dst = dist[imi][y*iw[imi] + x];
		d = x + 1;
		if (d < dst)
		  dst = d;
		d = iw[imi] - x;
		if (d < dst)
		  dst = d;
		d = y + 1;
		if (d < dst)
		  dst = d;
		d = ih[imi] - y;
		if (d < dst)
		  dst = d;
		dist[imi][y*iw[imi] + x] = dst;
		ix = (int) (dst + 0.5);
		if (ix >= blk_w/2)
		  windowed[imi][y*iw[imi]+x] = image_in[imi][y*iw[imi] + x] - mean[imi];
		else
		  windowed[imi][y*iw[imi]+x] = window[ix] *
		    (image_in[imi][y*iw[imi] + x] - mean[imi]);
	      }
	}
      else
	for (y = 0; y < ih[imi]; ++y)
	  for (x = 0; x <iw[imi]; ++x)
	    windowed[imi][y*iw[imi]+x] = image_in[imi][y*iw[imi] + x] - mean[imi];

      if (rFactor == 1)
	for (y = 0; y < n; ++y)
	  for (x = 0; x < n; ++x)
	    {
	      if (x >= offset_x && x < offset_x + iw[imi] &&
		  y >= offset_y && y < offset_y + ih[imi])
		img[y*n+x] = windowed[imi][(y - offset_y)*iw[imi] + (x - offset_x)];
	      else
		img[y*n+x] = 0.0;
	    }
      else
	for (y = 0; y < n; ++y)
	  for (x = 0; x < n; ++x)
	    {
	      v = 0.0;
	      count = 0;
	      for (dy = 0; dy < sFactor; ++dy)
		for (dx = 0; dx < sFactor; ++dx)
		  {
		    xv = (x - offset_x + ((float) dx) / sFactor) * rFactor;
		    yv = (y - offset_y + ((float) dy) / sFactor) * rFactor;

		    ixv = (int) floor(xv);
		    iyv = (int) floor(yv);
		    if (ixv >= 0 && ixv < iw[imi]-1 &&
			iyv >= 0 && iyv < ih[imi]-1)
		      {
			rrx = xv - ixv;
			rry = yv - iyv;
			rv00 = windowed[imi][iyv * iw[imi] + ixv];
			rv01 = windowed[imi][(iyv+1) * iw[imi] + ixv];
			rv10 = windowed[imi][iyv * iw[imi] + (ixv+1)];
			rv11 = windowed[imi][(iyv+1) * iw[imi] + (ixv+1)];
			v += rv00 * (rrx - 1.0) * (rry - 1.0)
			  - rv10 * rrx * (rry - 1.0) 
			  - rv01 * (rrx - 1.0) * rry
			  + rv11 * rrx * rry;
			++count;
		      }
		  }
	      if (count > 0)
		img[y*n+x] = v / count;
	      else
		img[y*n+x] = 0.0;
	    }

      if (c.outputImages)
	{
	  v_max = 0.0;
	  for (y = 0; y < n; ++y)
	    for (x = 0; x < n; ++x)
	      {
		if (fabs(img[y*n+x]) > v_max)
		  v_max = fabs(img[y*n+x]);
	      }
	  sprintf(fn, "%s%s.filtered%d.pgm",
		  c.outputBasename, t.pair.pairName, imi);
	  f = fopen(fn, "w");
	  fprintf(f, "P5\n%d %d\n255\n", n, n);
	  for (y = 0; y < n; ++y)
	    for (x = 0; x < n; ++x)
	      {
		i = (int) (127.99 * (img[y*n + x] / v_max) + 128.0);
		if (i < 0)
		  i = 0;
		else if (i > 255)
		  i = 255;
		fputc(i, f);
	      }
	  fclose(f);
	}

      fftwf_execute(plan_img);
      //      printf("Completed FT of image\n");

      memcpy(fft_orig[imi], fft_img, n2_partial*sizeof(fftwf_complex));

      // no need to do fractional FFTs if RS is known
      if (c.minScale == c.maxScale && c.minTheta == c.maxTheta)
	continue;

      fft_expand(fft_img, n);
      fft_shift(fft_img, n);

      for (y = 0; y < n; ++y)
	for (x = 0; x < n; ++x)
	  //	  mag[y*n+x] = hypot(fft_img[y*n+x][0],
	  //				 fft_img[y*n+x][1]);
	  mag[y*n+x] = log(hypot(fft_img[y*n+x][0],
				 fft_img[y*n+x][1]));
      resMin[0] = -0.5 * n;

      for (res = 1; res < nRes; ++res)
	{
	  alpha = c.fracRes[res] / n;
	  resMin[res] = -0.5 * n * c.fracRes[res];
	  //	  printf("resMin[%d] = %f  c.fracRes[%d] = %f\n",
	  //		 res, resMin[res], res, c.fracRes[res]);
	  for (j = n; j < n_times_2; ++j)
	    {
	      Y[j][0] = 0.0;
	      Y[j][1] = 0.0;
	    }
	  for (j = 0; j < n; ++j)
	    {
	      theta = -M_PI*j*j*alpha;
	      cos_t[j] = cos(theta);
	      sin_t[j] = sin(theta);
	      theta = M_PI*(j-n_over_2)*(j-n_over_2)*alpha;
	      Z[j][0] = cos(theta);
	      Z[j][1] = sin(theta);
	    }
	  for (j = n; j < n_times_2; ++j)
	    {
	      theta = M_PI*(j-n_over_2-n_times_2)*(j-n_over_2-n_times_2)*alpha;
	      Z[j][0] = cos(theta);
	      Z[j][1] = sin(theta);
	    }
	  fftwf_execute(plan_Z);

	  for (y = 0; y < n; ++y)
	    {
	      for (j = 0; j < n; ++j)
		{
		  Y[j][0] = img[y*n+j] * cos_t[j];
		  Y[j][1] = img[y*n+j] * sin_t[j];
		}
	      fftwf_execute(plan_Y);
	      
	      for (j = 0; j < n_times_2; ++j)
		{
		  W[j][0] = fft_Y[j][0] * fft_Z[j][0] - fft_Y[j][1] * fft_Z[j][1];
		  W[j][1] = fft_Y[j][0] * fft_Z[j][1] + fft_Y[j][1] * fft_Z[j][0];
		}
	      fftwf_execute(plan_W);

	      for (j = 0; j < n; ++j)
		{
		  fft_img[y * n + j][0] =
		    (Z[j][0] * ifft_W[j][0] + Z[j][1] * ifft_W[j][1]) / (2*n);
		  fft_img[y * n + j][1] =
		    (Z[j][0] * ifft_W[j][1] - Z[j][1] * ifft_W[j][0]) / (2*n);
		}
	    }
      
	  for (x = 0; x < n; ++x)
	    {
	      for (j = 0; j < n; ++j)
		{
		  Y[j][0] = fft_img[j*n+x][0] * cos_t[j] -
		    fft_img[j*n+x][1] * sin_t[j];
		  Y[j][1] = fft_img[j*n+x][0] * sin_t[j] +
		    fft_img[j*n+x][1] * cos_t[j];
		}
	      fftwf_execute(plan_Y);
      
	      for (j = 0; j < n_times_2; ++j)
		{
		  W[j][0] = fft_Y[j][0] * fft_Z[j][0] - fft_Y[j][1] * fft_Z[j][1];
		  W[j][1] = fft_Y[j][0] * fft_Z[j][1] + fft_Y[j][1] * fft_Z[j][0];
		}
	      fftwf_execute(plan_W);
	      
	      for (j = 0; j < n; ++j)
		{
		  fft_img[j * n + x][0] =
		    (Z[j][0] * ifft_W[j][0] + Z[j][1] * ifft_W[j][1]) / (2*n);
		  fft_img[j * n + x][1] =
		    (Z[j][0] * ifft_W[j][1] - Z[j][1] * ifft_W[j][0]) / (2*n);
		}
	    }

	  for (y = 0; y < n; ++y)
	    for (x = 0; x < n; ++x)
	      //	      mag[res*n2+ y*n + x] = hypot(fft_img[y*n+x][0],
	      //					       fft_img[y*n+x][1]);
	      mag[res*n2+ y*n + x] = log(hypot(fft_img[y*n+x][0],
					       fft_img[y*n+x][1]));
	}

      if (c.outputImages)
	for (res = 0; res < nRes; ++res)
	  {
	    v = 0.0;
	    for (y = 0; y < n; ++y)
	      for (x = 0; x < n; ++x)
		v += mag[res*n2 + y*n + x];

	    // printf("[%d] avg mag[%d] = %f\n", imi, res, v / n2);
	    // printf("[%d] mag[%d][0,0] = %f\n", imi, res, mag[res * n2 +
	    //          n_over_2 * n +
	    //          n_over_2]);
	    // printf("[%d] mag[%d][4,4] = %f\n", imi, res, mag[res * n2 +
	    //   (n_over_2 + (1 << res)) * n +
	    //   n_over_2 + (1 << res)]);
	    // printf("[%d] mag[%d][1000,4] = %f\n", imi, res, mag[res * n2 +
	    //   (n_over_2 + (1 << res)) * n +
	    //   n_over_2 + 250 * (1 << res)]);
	    // printf("[%d] mag[%d][-1000,1000] = %f\n", imi, res, mag[res * n2 +
	    //   (n_over_2 + 250 * (1 << res)) * n +
	    //   n_over_2 - 250 * (1 << res)]);

	    v = 2.0 * v / n2;
	    sprintf(fn, "%s%s.dft%d.%d.pgm",
		    c.outputBasename, t.pair.pairName, imi, res);
	    f = fopen(fn, "w");
	    fprintf(f, "P5\n%d %d\n255\n", n, n);
	    for (y = 0; y < n; ++y)
	      for (x = 0; x < n; ++x)
		{
		  i = 255.99 * mag[res*n2 + y*n + x] / v;
		  if (i > 255)
		    i = 255;
		  fputc(i, f);
		}
	    fclose(f);
	  }

      for (it = 0; it < n; ++it)
	{
	  theta = it * M_PI / n;
	  cos_theta = cos(theta);
	  sin_theta = sin(theta);
	  for (ir = 0; ir < n; ++ir)
	    {
	      rho = exp(ir * logrhobase + logrhooffset);
	      xv = rho * cos_theta;
	      yv = rho * sin_theta;
	      for (res = nRes - 1; res >= 0; --res)
		{
		  xvp = (xv - resMin[res]) / c.fracRes[res];
		  yvp = (yv - resMin[res]) / c.fracRes[res];
		  ixv = (int) floor(xvp);
		  iyv = (int) floor(yvp);
		  if (ixv >= 0 && ixv < n-1 &&
		      iyv >= 0 && iyv < n-1)
		    {
		      rrx = xvp - ixv;
		      rry = yvp - iyv;
		      rv00 = mag[res*n2 + iyv * n + ixv];
		      rv01 = mag[res*n2 + (iyv+1) * n + ixv];
		      rv10 = mag[res*n2 + iyv * n + (ixv+1)];
		      rv11 = mag[res*n2 + (iyv+1) * n + (ixv+1)];
		      v = rv00 * (rrx - 1.0) * (rry - 1.0)
			- rv10 * rrx * (rry - 1.0) 
			- rv01 * (rrx - 1.0) * rry
			+ rv11 * rrx * rry;
		      break;
		    }
		}
	      if (res < 0)
		Error("Could not find point in any resolution level: %d %d %f %f  %f %f %f\n", it, ir, xv, yv, resMin[0], resMin[1], resMin[2]);

#if 0
	      hpfx = cos(xv * M_PI / n);
	      hpfy = cos(yv * M_PI / n);
	      filter = (1.0 - hpfx*hpfy) * (2.0 - hpfx*hpfy);
#else
	      filter = 1.0;
#endif
	      lp[ir*n+it] = filter * v;
	    }
	}

      if (c.outputImages)
	{
	  v = 0.0;
	  for (ir = 0; ir < n; ++ir)
	    for (it = 0; it < n; ++it)
	      if (lp[ir*n+it] > v)
		v = lp[ir*n+it];
	  sprintf(fn, "%s%s.logpolar%d.pgm",
		  c.outputBasename, t.pair.pairName, imi);
	  f = fopen(fn, "w");
	  fprintf(f, "P5\n%d %d\n255\n", n, n);
	  for (ir = 0; ir < n; ++ir)
	    for (it = 0; it < n; ++it)
	      fputc((int) (255.99 * lp[ir*n+it] / v), f);
	  fclose(f);
	}
      fftwf_execute(plan_lp);
      memcpy(flp[imi], fft_lp, n2_partial * sizeof(fftwf_complex));
    }

  if (c.minScale == c.maxScale && c.minTheta == c.maxTheta)
    {
      candidates[0].x = fmod(c.minTheta + 360.0, 180.0) * n / 180.0;
      if (c.minScale >= 1.0)
	candidates[0].y = (int) (log(c.minScale) / logrhobase + 0.5);
      else
	candidates[0].y = (n + (int) floor(log(c.minScale) / logrhobase + 0.5)) % n;
      candidates[0].val = 1.0;
      candidates[0].radius = 0.0;
      nCandidates = 1;
      goto findTranslation;
    }

  for (i = 0; i < n2_partial; ++i)
    {
      p_r = flp[0][i][0] * flp[1][i][0] + flp[0][i][1] * flp[1][i][1];
      p_i = flp[0][i][1] * flp[1][i][0] - flp[0][i][0] * flp[1][i][1];
      d = hypot(p_r, p_i);
      prod[i][0] = p_r / d;
      prod[i][1] = p_i / d;
    }

#if 0
  printf("Completed construction of prod\n");
  count = 0;
  for (i = 0; i < n2_partial; ++i)
    {
      if (fabs(prod[i][0]) > 0.001 ||
	  fabs(prod[i][1]) > 0.001)
	++count;
    }
  printf("%d entries in prod are nonzero\n", count);

  printf("max_x = %d max_y = %d\n", max_x, max_y);
  printf("   val at max = %f\n", cc[max_y*n+max_x]);
  printf("%f %f\n", cc[256], cc[512+256]);

  printf("theta = %f degrees\n", (180.0 * max_x) / (n - 1));
  printf("factor = %f\n", max_y < (n / 2) ? exp(max_y * logrhobase) :
	 exp((max_y - n) * logrhobase));
#endif

  nCandidates = 0;
  //  printf("minITSeparation = %f  minIRSeparation = %f\n",
  //	 minITSeparation, minIRSeparation);

  /* convolve with a Gaussian of the given radius */
  if (radius == 0.0)
    {
      for (x = 0; x < n; ++x)
	gaussian[x] = 1.0 / n;
    }
  else
    {
      // n = 64  radius = 1 sigma = 10.185
      sigma = n / (2.0 * M_PI * radius);
      memset(gaussian, 0, n * sizeof(float));
      for (i = 0; ; ++i)
	{
	  v = exp(- i * i / (2.0 * sigma * sigma)) /
	    (sigma * sqrt(2.0 * M_PI));
	  gaussian[i % n] += v;
	  if (v < 1.0e-10)
	    break;
	}
      // printf("final i = %d\n", i);
      for (i = -1; ; --i)
	{
	  v = exp(- i * i / (2.0 * sigma * sigma)) /
	    (sigma * sqrt(2.0 * M_PI));
	  gaussian[(i + n2) % n] += v;
	  if (v < 1.0e-10 || i == -n2)
	    break;
	}
      // printf("final i = %d\n", i);
    }
	  
  // apply gaussian directly to partial FT
  for (y = 0; y < n; ++y)
    for (x = 0; x < n_over_2_plus_1; ++x)
      {
	fft_cc[y*n_over_2_plus_1+x][0] =
	  prod[y*n_over_2_plus_1+x][0] * gaussian[x] * gaussian[y];
	fft_cc[y*n_over_2_plus_1+x][1] =
	  prod[y*n_over_2_plus_1+x][1] * gaussian[x] * gaussian[y];
      }

#if 0
  count = 0;
  for (i = 0; i < n2_partial; ++i)
    {
      if (fabs(fft_cc[i][0]) > 0.001 ||
	  fabs(fft_cc[i][1]) > 0.001)
	++count;
    }
  printf("%d entries in fft_cc are nonzero\n", count);
#endif

  // perform inverse FFT
  fftwf_execute(plan_cc);
  // printf("Completed IFT of fft_cc\n");

  if (c.outputImages)
    {
      v_max = 0.0;
      v_min = 0.0;
      for (y = 0; y < n; ++y)
	for (x = 0; x < n; ++x)
	  {
	    if (cc[y*n+x] > v_max)
	      v_max = cc[y*n+x];
	    if (cc[y*n+x] < v_min)
	      v_min = cc[y*n+x];
	  }
      sprintf(fn, "%s%s.rscc.pgm",
	      c.outputBasename, t.pair.pairName);
      f = fopen(fn, "w");
      fprintf(f, "P5\n%d %d\n255\n", n, n);
      for (y = 0; y < n; ++y)
	for (x = 0; x < n; ++x)
	  {
	    j = 255.99 * (cc[y*n+x] - v_min) / (v_max - v_min);
	    if (j > 255)
	      j = 255;
	    fputc(j, f);
	  }
      fclose(f);
    }

  // sort by value
  i = 0; 
  for (ir = 0; ir < n; ++ir)
    for (it = 0; it < n; ++it)
      {
	values[i].x = it;
	values[i].y = ir;
	values[i].val = cc[ir*n+it];
	values[i].radius = radius;
	++i;
      }
  //      printf("values[0].val = %f\n", values[0].val);
  qsort(values, n2, sizeof(PositionValue), CompareValues);
  //      printf("values[0].val = %f\n", values[0].val);
      
  // pick the top candidates
  memset(marked, 0, n2*sizeof(unsigned char));
  nCandidatesAdded = 0;

  target_it = fmod((c.minTheta + c.maxTheta) / 2.0, 180.0) * n / 180.0;
  target_ir = 0.5 * log(c.minScale * c.maxScale) / logrhobase;
  if (target_ir < 0)
    target_ir += n;
  max_delta_it = ((c.maxTheta - c.minTheta) / 2.0 / 180.0) * n;
  max_delta_ir = log(sqrt(c.maxScale / c.minScale)) / logrhobase;
  //     printf("target_it = %d   mdit = %d\n", target_it, max_delta_it);
  //     printf("target_ir = %d   mdir = %d\n", target_ir, max_delta_ir);

  for (i = 0; i < n2; ++i)
    {
      it = values[i].x;
      ir = values[i].y;

      /* check if in valid region */
      delta_it = abs(it - target_it);
      if (delta_it > n_over_2)
	delta_it = n - delta_it;
      delta_ir = abs(ir - target_ir);
      if (delta_ir > n_over_2)
	delta_ir = n - delta_ir;
      if (delta_it > max_delta_it ||
	  delta_ir > max_delta_ir)
	continue;

      /* check if any neighbors have been marked */
      for (dy = -1; dy <= 1; ++dy) 
	for (dx = -1; dx <= 1; ++dx)
	  {
	    nit = it + dx;
	    nir = ir + dy;
	    if (nit < 0 || nit >= n ||
		nir < 0 || nir >= n ||
		dx == 0 && dy == 0)
	      continue;
	    if (marked[nir*n+nit])
	      goto nextPosition;
	  }

      /* check if far enough from already chosen candidates */
      for (j = 0; j < nCandidates; ++j)
	{
	  delta_it = abs(it - candidates[j].x);
	  if (delta_it > n_over_2)
	    delta_it = n - delta_it;
	  delta_ir = abs(ir - candidates[j].y);
	  if (delta_ir > n_over_2)
	    delta_ir = n - delta_ir;
	  if (delta_it < minITSeparation &&
	      delta_ir < minIRSeparation)
	    goto nextPosition;
	}
      
      /* add to candidate list */ 
      candidates[nCandidates].x = it;
      candidates[nCandidates].y = ir;
      candidates[nCandidates].val = values[i].val;
      candidates[nCandidates].radius = values[i].radius;
      ++nCandidates;
      if (++nCandidatesAdded >= c.maxRSCandidates)
	break;
	  
    nextPosition:
      marked[ir*n+it] = 1;
    }

 findTranslation:
  /* print out all candidates */
  for (i = 0; i < nCandidates; ++i)
    Log("candidate %d:  rot = %f or %f scale = %f\n val = %f radius = %f\n",
	i,
	180.0 * candidates[i].x / n,
	180.0 * candidates[i].x / n + 180.0,
	candidates[i].y < n_over_2 ? exp(candidates[i].y * logrhobase) :
	exp((candidates[i].y - n) * logrhobase),
	candidates[i].val,
	candidates[i].radius);

  /* create the final candidate array */
  nFinalCandidates = 0;

  for (x = 0; x < n_times_2; ++x)
    {
      rot[x][0] = cos(2.0 * x * M_PI / n);
      rot[x][1] = sin(2.0 * x * M_PI / n);
    }

  /* find translations for all candidates */
  for (i = 0; i < nCandidates; ++i)
    {
      /* rotate and scale the larger image down to the smaller
         using subsampling */

      if (candidates[i].y < n_over_2)
	{
	  /* the reference is at larger magnification than the image;
	     we will transform and resample the reference image so that
	     it matches image 0 in scale and rotation */
	  scale = exp(candidates[i].y * logrhobase);
	  baseRotation = -candidates[i].x * M_PI / n;
	  smaller = 0;
	  larger = 1;
	}
      else
	{
	  /* the image is at larger magnification than the reference;
	   we will transform and resample image 0 so that it matches
	   the reference in scale and rotation */
	  scale = exp(-(candidates[i].y - n) * logrhobase);
	  baseRotation = candidates[i].x * M_PI / n;
	  smaller = 1;
	  larger = 0;
	}
      /* if rotation is very close to a multiple of 180 degrees,
	 make it exactly 0 degrees */
      if (fabs(baseRotation) < 0.001 || fabs(fabs(baseRotation) - M_PI) < 0.001)
	baseRotation = 0.0;

      ct = cos(baseRotation);
      st = sin(baseRotation);
      /* if scaling is very close to 1.0, make it exactly 1.0 */
      if (fabs(scale - 1.0) < 0.001)
	scale = 1.0;

      Log("Considering RS candidate %d: rotation = %f scale = %f smaller = %d larger = %d\n",
	  i, baseRotation * 180.0 / M_PI, scale, smaller, larger);

      offset_x = (n - eiw[larger]) >> 1;
      offset_y = (n - eih[larger]) >> 1;

      if (rFactor == 1 && scale == 1.0 && baseRotation == 0.0)
	for (y = 0; y < n; ++y)
	  for (x = 0; x < n; ++x)
	    {
	      if (x >= offset_x && x < offset_x + iw[larger] &&
		  y >= offset_y && y < offset_y + ih[larger])
		img[y*n+x] = windowed[larger][(y - offset_y)*iw[larger] + (x - offset_x)];
	      else
		img[y*n+x] = mean[larger];
	    }
      else
	for (y = 0; y < n; ++y)
	  for (x = 0; x < n; ++x)
	    {
	      v = 0.0;
	      count = 0;
	      for (dy = 0; dy < sFactor; ++dy)
		for (dx = 0; dx < sFactor; ++dx)
		  {
		    xv = (x - n_over_2) + dx / ((float) sFactor);
		    yv = (y - n_over_2) + dy / ((float) sFactor);
		    
		    xvp = scale * (ct * xv - st * yv) + n_over_2 - offset_x;
		    yvp = scale * (st * xv + ct * yv) + n_over_2 - offset_y;

		    xvp *= rFactor;
		    yvp *= rFactor;

		    ixv = (int) floor(xvp);
		    iyv = (int) floor(yvp);
		    if (ixv >= 0 && ixv < iw[larger]-1 &&
			iyv >= 0 && iyv < ih[larger]-1)
		      {
			rrx = xvp - ixv;
			rry = yvp - iyv;
			rv00 = image_in[larger][iyv * iw[larger] + ixv];
			rv01 = image_in[larger][(iyv+1) * iw[larger] + ixv];
			rv10 = image_in[larger][iyv * iw[larger] + (ixv+1)];
			rv11 = image_in[larger][(iyv+1) * iw[larger] + (ixv+1)];
			v += rv00 * (rrx - 1.0) * (rry - 1.0)
			  - rv10 * rrx * (rry - 1.0) 
			  - rv01 * (rrx - 1.0) * rry
			  + rv11 * rrx * rry;
			++count;

			dst = dist[larger][iyv * iw[larger] + ixv];
		      }
		  }
	      if (count > 0)
		{
		  ix = (int) (dst * scale + 0.5);
		  if (ix >= blk_w/2)
		    img[y*n+x] = (v / count) - mean[larger];
		  else
		    img[y*n+x] = window[ix] * (v / count - mean[larger]);
		}
	      else
		img[y*n+x] = 0.0;
	    }

      if (c.outputImages && i == 0)
	{
	  v_max = 0.0;
	  for (y = 0; y < n; ++y)
	    for (x = 0; x < n; ++x)
	      {
		if (fabs(img[y*n+x]) > v_max)
		  v_max = fabs(img[y*n+x]);
	      }
	  sprintf(fn, "%s%s.rsimg.pgm",
		  c.outputBasename, t.pair.pairName);
	  f = fopen(fn, "w");
	  fprintf(f, "P5\n%d %d\n255\n", n, n);
	  for (y = 0; y < n; ++y)
	    for (x = 0; x < n; ++x)
	      {
		j = (int) (127.99 * (img[y*n + x] / v_max) + 128.0);
		if (j < 0)
		  j = 0;
		else if (j > 255)
		  j = 255;
		fputc(j, f);
	      }
	  fclose(f);
	}

      /* perform FFT */
      fftwf_execute(plan_img);

      for (th = 0; th <= 180; th += 180)
	{
	  rotation = baseRotation + (th / 180) * M_PI;
	  /* if the rotation is not in the specified range, ignore */
	  if (fmod(rotation * 180.0 / M_PI - c.minTheta + 360.0, 360.0) >
	      c.maxTheta - c.minTheta + 0.001)
	    continue;
	  ct = cos(rotation);
	  st = sin(rotation);

	  if (th == 180)
	    {
	      /* obtain FT of image rotated by 180 degrees by
		 rotating the phases of the conjugate of the FT of
		 the unrotated image */
	      for (y = 0; y < n; ++y)
		for (x = 0; x < n_over_2_plus_1; ++x)
		  {
		    vr = fft_img[y*n_over_2_plus_1+x][0];
		    vi = fft_img[y*n_over_2_plus_1+x][1];
		    fft_img[y*n_over_2_plus_1+x][0] =
		      vr * rot[x+y][0] + vi * rot[x+y][1];
		    fft_img[y*n_over_2_plus_1+x][1] =
		      vr * rot[x+y][1] - vi * rot[x+y][0];
		  }
	    }

	  /* compute cross-power spectrum */
	  for (j = 0; j < n2_partial; ++j)
	    {
	      p_r = fft_orig[smaller][j][0] * fft_img[j][0] +
		fft_orig[smaller][j][1] * fft_img[j][1];
	      p_i = fft_orig[smaller][j][1] * fft_img[j][0] -
		fft_orig[smaller][j][0] * fft_img[j][1];
	      d = hypot(p_r, p_i);
	      prod[j][0] = p_r / d;
	      prod[j][1] = p_i / d;
	    }

	  /* convolve with a Gaussian of the given radius */
	  if (radius == 0.0)
	    {
	      for (x = 0; x < n; ++x)
		gaussian[x] = 1.0 / n;
	    }
	  else
	    {
	      // n = 64  radius = 1 sigma = 10.185
	      sigma = n / (2.0 * M_PI * radius);
	      memset(gaussian, 0, n * sizeof(float));
	      for (j = 0; ; ++j)
		{
		  v = exp(- j * j / (2.0 * sigma * sigma)) /
		    (sigma * sqrt(2.0 * M_PI));
		  gaussian[j % n] += v;
		  if (v < 1.0e-10)
		    break;
		}
	      //		  printf("final j = %d\n", j);
	      for (j = -1; ; --j)
		{
		  v = exp(- j * j / (2.0 * sigma * sigma)) /
		    (sigma * sqrt(2.0 * M_PI));
		  gaussian[(j + 1024 * n) % n] += v;
		  if (v < 1.0e-10)
		    break;
		}
	      //		  printf("final j = %d\n", j);
	    }
	  
	  // apply gaussian directly to partial FT
	  for (y = 0; y < n; ++y)
	    for (x = 0; x < n_over_2_plus_1; ++x)
	      {
		fft_cc[y*n_over_2_plus_1+x][0] =
		  prod[y*n_over_2_plus_1+x][0] * gaussian[x] * gaussian[y];
		fft_cc[y*n_over_2_plus_1+x][1] =
		  prod[y*n_over_2_plus_1+x][1] * gaussian[x] * gaussian[y];
	      }

	  if (ccFilterMin > 0 || ccFilterMax < n_over_2)
	    {
	      // zero all the components outside of the desired frequency range
	      memcpy(ccfilter, fft_cc, n2_partial * sizeof(fftwf_complex));
	      fft_expand(ccfilter, n);
	      fft_shift(ccfilter, n);
	      for (y = 0; y < n; ++y)
		for (x = 0; x < n; ++x)
		  {
		    dx = abs(x - n_over_2);
		    dy = abs(x - n_over_2);
		    if (dx < ccFilterMin || dx > ccFilterMax ||
			dy < ccFilterMin || dy > ccFilterMax)
		      {
			ccfilter[y*n+x][0] = 0.0;
			ccfilter[y*n+x][1] = 0.0;
		      }
		  }
	      fft_shift(ccfilter, n);
	      fft_compress(ccfilter, n);
	      memcpy(fft_cc, ccfilter, n2_partial * sizeof(fftwf_complex));
	    }

	  /* perform IFFT */
	  fftwf_execute(plan_cc);

	  /* generate a picture */
	  if (c.outputImages && i == 0)
	    {
	      v_max = 0.0;
	      v_min = 0.0;
	      for (y = 0; y < n; ++y)
		for (x = 0; x < n; ++x)
		  {
		    if (cc[y*n+x] > v_max)
		      v_max = cc[y*n+x];
		    if (cc[y*n+x] < v_min)
		      v_min = cc[y*n+x];
		  }
	      sprintf(fn, "%s%s.cc.pgm",
		      c.outputBasename, t.pair.pairName);
	      f = fopen(fn, "w");
	      fprintf(f, "P5\n%d %d\n255\n", n, n);
	      for (y = 0; y < n; ++y)
		for (x = 0; x < n; ++x)
		  {
		    j = 255.99 * (cc[y*n+x] - v_min) / (v_max - v_min);
		    if (j > 255)
		      j = 255;
		    fputc(j, f);
		  }
	      fclose(f);
	    }

	  /* sort the values */
	  j = 0; 
	  for (y = 0; y < n; ++y)
	    for (x = 0; x < n; ++x)
	      {
		values[j].x = x;
		values[j].y = y;
		values[j].val = cc[y*n+x];
		values[j].radius = radius;
		++j;
	      }
	  //	  Log("before values[0].val = %f\n", values[0].val);
	  qsort(values, n2, sizeof(PositionValue), CompareValues);
	  //	  Log("after values[0].val = %f  radius = %f\n", values[0].val,
	  //	      radius);

	  /* determine highest peaks in the normalized cross-correlation */
	  memset(marked, 0, n2*sizeof(unsigned char));
	  for (j = 0; j < n2; ++j)
	    {
	      quality = values[j].val;
	      ix = values[j].x;
	      iy = values[j].y;

	      if (quality < 0.0 ||
		  nFinalCandidates == c.maxCandidates &&
		  quality <= finalCandidates[nFinalCandidates-1].quality)
		break;
		      
	      if (ix <= n_over_2)
		deltaX = -(float) ix;
	      else
		deltaX = - (float) (ix - n);
	      if (iy <= n_over_2)
		deltaY = -(float) iy;
	      else
		deltaY = -(float) (iy - n);

	      /* check if in valid region */
	      if (deltaX * rFactor < c.minTX * iw[0] * 0.01 || deltaX * rFactor > c.maxTX * iw[0] * 0.01 ||
		  deltaY * rFactor < c.minTY * ih[0] * 0.01 || deltaY * rFactor > c.maxTY * ih[0] * 0.01)
		{
		  //		  Log("Rejecting qual %f (val %f) candidate %d %d (%f %f) because not in valid region\n",
		  //		      quality, values[j].val,
		  //		      ix, iy, tx, ty);
		  continue;
		}

	      if (larger == 1)
		{
		  // from transformation 0->1 in comments at beginning of program:
		  tx = scale * (ct * (- iw[0]/2 + deltaX * rFactor) - st * (- ih[0]/2 + deltaY * rFactor)) + iw[1]/2;
		  ty = scale * (st * (- iw[0]/2 + deltaX * rFactor) + ct * (- ih[0]/2 + deltaY * rFactor)) + ih[1]/2;
		}
	      else
		{
		  // from transformation 1->0 in comments at beginning of program, but with 0 and 1 reversed:
		  tx = (1.0/scale) * (-ct * iw[0] / 2  - st * ih[0] / 2.0) + iw[1] / 2.0 - deltaX * rFactor;
		  ty = (1.0/scale) * (st * iw[0] / 2 - ct * ih[0] / 2.0) + ih[1] / 2.0 - deltaY * rFactor;
		}
	      //	      printf("Considering candidate %d: rot=%f scale=%f otx=%f oty=%f tx=%f ty=%f qual=%f rad=%f\n",
	      //		     j, rotation * 180.0 / M_PI, scale, otx, oty, tx, ty, quality, radius);

	      /* check if any neighbors have been marked */
	      for (dy = -1; dy <= 1; ++dy) 
		for (dx = -1; dx <= 1; ++dx)
		  {
		    nix = ix + dx;
		    niy = iy + dy;
		    if (nix < 0 || nix >= n ||
			niy < 0 || niy >= n ||
			dx == 0 && dy == 0)
		      continue;
		    if (marked[niy*n+nix])
		      goto nextFinalPosition;
		  }

	      /* check if far enough from already chosen candidates */
	      for (k = 0; k < nFinalCandidates; ++k)
		{
		  if (finalCandidates[k].rotation != rotation ||
		      finalCandidates[k].scale != scale)
		    continue;
		  if (quality > finalCandidates[k].quality)
		    break;
		  if (fabs(finalCandidates[k].tx - tx) < c.minTranslationalSeparation * 0.01 * iw[0] &&
		      fabs(finalCandidates[k].ty - ty) < c.minTranslationalSeparation * 0.01 * ih[0])
		    goto nextFinalPosition;
		}

	      /* add to candidate list, keeping candidates in sorted order */ 
	      if (nFinalCandidates == c.maxCandidates)
		if (quality > finalCandidates[nFinalCandidates-1].quality)
		  --nFinalCandidates;
		else
		  break;
	      for (k = 0; k < nFinalCandidates; ++k)
		if (quality > finalCandidates[k].quality)
		  break;
	      if (k < nFinalCandidates)
		memmove(&finalCandidates[k+1], &finalCandidates[k],
			(nFinalCandidates - k) * sizeof(Transformation));
	      if (larger == 1)
		{
		  finalCandidates[k].rotation = rotation;
		  finalCandidates[k].scale = scale;
		}
	      else
		{
		  finalCandidates[k].rotation = -rotation;
		  finalCandidates[k].scale = 1.0 / scale;
		}
	      finalCandidates[k].tx = tx;
	      finalCandidates[k].ty = ty;
	      finalCandidates[k].quality = quality;
	      finalCandidates[k].radius = radius;
	      Log("Added final candidate %d: %f %f %f %f %f %f\n",
		  k, rotation * 180.0 / M_PI, scale, tx, ty, quality, radius);

	      ++nFinalCandidates;

	      /* eliminate any other candidates that are close and lower in quality */
	      ++k;
	      while (k < nFinalCandidates)
		{
		  if (finalCandidates[k].rotation != rotation ||
		      finalCandidates[k].scale != scale ||
		      fabs(finalCandidates[k].tx - tx) >= c.minTranslationalSeparation * 0.01 * iw[0] ||
		      fabs(finalCandidates[k].ty - ty) >= c.minTranslationalSeparation * 0.01 * ih[0])
		    {
		      ++k;
		      continue;
		    }
		  if (k < nFinalCandidates-1)
		    memmove(&finalCandidates[k], &finalCandidates[k+1],
			    (nFinalCandidates - 1 - k) *
			    sizeof(Transformation));
		  --nFinalCandidates;
		}

	    nextFinalPosition:
	      marked[iy*n+ix] = 1;
	    }
	}
    }

  Log("nFinalCandidates = %d\n", nFinalCandidates);
  for (i = 0; i < nFinalCandidates; ++i)
      Log("  final candidate %d:  rot = %f scale = %f\n tx = %f ty = %f quality = %f radius = %f\n",
	  i,
	  finalCandidates[i].rotation * 180.0 / M_PI,
	  finalCandidates[i].scale,
	  finalCandidates[i].tx,
	  finalCandidates[i].ty,
	  finalCandidates[i].quality,
	  finalCandidates[i].radius);

  if (outputName[0] != '\0')
    {
      /* turn the best candidate into a map */
      rotation = finalCandidates[0].rotation;
      scale = finalCandidates[0].scale;
      tx = finalCandidates[0].tx;
      ty = finalCandidates[0].ty;
      offset_x = (n - eiw[0]) >> 1;
      offset_y = (n - eih[0]) >> 1;
      offset2_x = (n - eiw[1]) >> 1;
      offset2_y = (n - eih[1]) >> 1;
      ct = cos(rotation);
      st = sin(rotation);
#if 0
      printf("rot = %f scale = %f tx = %f ty = %f\n", rotation, scale, tx, ty);
      printf("offset_x = %d offset_y = %d\n", offset_x, offset_y);
      printf("offset2_x = %d offset2_y = %d\n", offset2_x, offset2_y);
      printf("ct = %f  st = %f\n", ct, st);
#endif

      if (t.pair.imageMinX >= 0)
	iMinX = t.pair.imageMinX;
      else
	iMinX = 0;
      if (t.pair.imageMaxX >= 0)
	iMaxX = t.pair.imageMaxX;
      else
	iMaxX = iMinX + iw[0];
      if (t.pair.imageMinY >= 0)
	iMinY = t.pair.imageMinY;
      else
	iMinY = 0;
      if (t.pair.imageMaxY >= 0)
	iMaxY = t.pair.imageMaxY;
      else
	iMaxY = iMinY + ih[0];

      if (c.mapLevel >= 0)
	level = c.mapLevel;
      else
	for (level = 0; iMaxX > (1<<level) || iMaxY > (1<<level); ++level) ;
      md = 1 << level;

      mMinX = iMinX >> level;
      mMaxX = (iMaxX + md - 1) >> level;
      mMinY = iMinY >> level;
      mMaxY = (iMaxY + md - 1) >> level;
      mw = mMaxX - mMinX + 1;
      mh = mMaxY - mMinY + 1;
      map = (MapElement*) malloc(mh * mw * sizeof(MapElement));
      //      printf("mw = %d mh = %d md = %d\n", mw, mh, md);

      for (y = mMinY; y <= mMaxY; ++y)
	{
	  yv = (float) (y * md - iMinY);
	  for (x = mMinX; x <= mMaxX; ++x)
	    {
	      xv = (float) (x * md - iMinX);
	      xvp = scale * (ct * xv - st * yv) + tx;
	      yvp = scale * (st * xv + ct * yv) + ty;
	      i = (y - mMinY) * mw + (x - mMinX);
	      map[i].x = xvp / md;
	      map[i].y = yvp / md;
	      map[i].c = 1.0;
	      //	      printf("Computed map [%d %d]: (%f %f) (%f %f) %d (%f %f)\n",
	      //		     x, y, xv, yv, xvp, yvp, i, map[i].x, map[i].y);
	    }
	}

      /* write out the map file */
      sprintf(mapName, "%s.map", outputName);
      if (!WriteMap(mapName, map, level, mw, mh, mMinX, mMinY,
		    t.pair.imageName, t.pair.refName,
		    UncompressedMap, msg))
	Error("%s\n", msg);
      free(map);
    }

  // write out the score for this mapping
  if (c.outputBasename[0] != '\0')
    {
      sprintf(scoreName, "%s.score", outputName);
      f = fopen(scoreName, "w");
      fprintf(f, "%f %f %f %f %f\n",
	      finalCandidates[0].quality,
	      finalCandidates[0].rotation * 180.0 / M_PI,
	      finalCandidates[0].scale, 
	      finalCandidates[0].tx,
	      finalCandidates[0].ty);
      fclose(f);
    }

  memcpy(&(r.pair), &(t.pair), sizeof(Pair));
  r.updated = 1;
  r.quality = finalCandidates[0].quality;
  r.rotation = finalCandidates[0].rotation;
  r.scale = finalCandidates[0].scale;
  r.tx = finalCandidates[0].tx;
  r.ty = finalCandidates[0].ty;
  r.message[0] = '\0';

  for (imi = 0; imi < 2; ++imi)
    {
      free(image_in[imi]);
      image_in[imi] = NULL;
      if (image_mask_in[imi] != NULL)
	{
	  free(image_mask_in[imi]);
	  image_mask_in[imi] = NULL;
	}
    }
  free(candidates);
  free(finalCandidates);
}

int CompareValues (const void *v1, const void *v2)
{
  if (((PositionValue*) v1)->val > ((PositionValue*) v2)->val)
    return(-1);
  else if (((PositionValue*) v1)->val == ((PositionValue*) v2)->val)
    return(0);
  else
    return(1);
}

void
fft_shift (fftwf_complex *fft, int n)
{
  int n2;
  int x, y;
  float tmp_r, tmp_i;

  n2 = n / 2;

  // swap quadrants 1 and 3, and 2 and 4
  for (x = 0; x < n2; x++)
    for (y = 0; y < n2; y++)
      {
	tmp_r = fft[y*n+x][0];
	tmp_i = fft[y*n+x][1];
	fft[y*n+x][0] = fft[(y+n2)*n+(x+n2)][0];
	fft[y*n+x][1] = fft[(y+n2)*n+(x+n2)][1];
	fft[(y+n2)*n+(x+n2)][0] = tmp_r;
	fft[(y+n2)*n+(x+n2)][1] = tmp_i;
	
	tmp_r = fft[y*n+(x+n2)][0];
	tmp_i = fft[y*n+(x+n2)][1];
	fft[y*n+(x+n2)][0] = fft[(y+n2)*n+x][0];
	fft[y*n+(x+n2)][1] = fft[(y+n2)*n+x][1];
	fft[(y+n2)*n+x][0] = tmp_r;
	fft[(y+n2)*n+x][1] = tmp_i;
      }
}

void
fft_expand (fftwf_complex *fft, int n)
{
  int n2, n21;
  int x, y;
  float tmp_r, tmp_i;

  n2 = n / 2;
  n21 = n2 + 1;

  // shift into final position
  for (y = n-1; y > 0; --y)
    for (x = n2; x >= 0; --x)
      {
	fft[y*n+x][0] = fft[y*n21+x][0];
	fft[y*n+x][1] = fft[y*n21+x][1];
      }

  // fill in other entries
  for (y = 0; y < n; ++y)
    for (x = n21; x < n; ++x)
      {
	fft[y*n+x][0] = fft[((n-y) % n)*n + (n-x)][0];
	fft[y*n+x][1] = -fft[((n-y) % n)*n + (n-x)][1];
      }
}

void
fft_compress (fftwf_complex *fft, int n)
{
  int n2, n21;
  int x, y;
  float tmp_r, tmp_i;

  n2 = n / 2;
  n21 = n2 + 1;

  // shift into final position
  for (y = 0; y < n; ++y)
    for (x = 0; x < n21; ++x)
      {
	fft[y*n21+x][0] = fft[y*n+x][0];
	fft[y*n21+x][1] = fft[y*n+x][1];
      }

  // zero other entries
  memset(&fft[n*n21], 0, n*(n-n21)*sizeof(fftwf_complex));
}

int
SortByQuality (const void *x, const void *y)
{
  Result *rx, *ry;
  rx = (Result *) x;
  ry = (Result *) y;
  if (rx->quality < ry->quality)
    return(-1);
  else if (rx->quality == ry->quality)
    return(0);
  else
    return(1);
}

int
SortByName (const void *x, const void *y)
{
  Result *rx, *ry;
  rx = (Result *) x;
  ry = (Result *) y;
  return(strcmp(rx->pair.imageName, ry->pair.imageName));
}

void
PackContext ()
{
  int i;

  par_pkbyte((unsigned char) c.type);
  par_pkstr(c.imageBasename);
  par_pkstr(c.imageMaskBasename);
  par_pkstr(c.referenceBasename);
  par_pkstr(c.referenceMaskBasename);
  par_pkstr(c.outputBasename);
  par_pkstr(c.logBasename);
  par_pkint(c.mapLevel);
  par_pkint(c.outputImages);

  par_pkint(c.maxRes);
  par_pkfloat(c.distortion);
  par_pkfloat(c.margin);
  par_pkfloat(c.minFeatureSize);
  par_pkfloat(c.maxFeatureSize);
  par_pkfloat(c.minTransFeatureSize);
  par_pkfloat(c.maxTransFeatureSize);
  par_pkfloat(c.minRotationalSeparation);
  par_pkfloat(c.minScaleSeparation);
  par_pkfloat(c.minTranslationalSeparation);
  for (i = 0; i < MAX_FRAC_FT_RES_LEVELS; ++i)
    par_pkfloat(c.fracRes[i]);
  par_pkint(c.maxRSCandidates);
  par_pkint(c.maxCandidates);

  par_pkfloat(c.minTheta);
  par_pkfloat(c.maxTheta);
  par_pkfloat(c.minScale);
  par_pkfloat(c.maxScale);

  par_pkfloat(c.minTX);
  par_pkfloat(c.maxTX);
  par_pkfloat(c.minTY);
  par_pkfloat(c.maxTY);

  par_pkint(c.update);
  par_pkint(c.partial);
  par_pkint(c.nWorkers);
}

void
UnpackContext ()
{
  int i;

  c.type = (char) par_upkbyte();
  par_upkstr(c.imageBasename);
  par_upkstr(c.imageMaskBasename);
  par_upkstr(c.referenceBasename);
  par_upkstr(c.referenceMaskBasename);
  par_upkstr(c.outputBasename);
  par_upkstr(c.logBasename);
  c.mapLevel = par_upkint();
  c.outputImages = par_upkint();

  c.maxRes = par_upkint();
  c.distortion = par_upkfloat();
  c.margin = par_upkfloat();
  c.minFeatureSize = par_upkfloat();
  c.maxFeatureSize = par_upkfloat();
  c.minTransFeatureSize = par_upkfloat();
  c.maxTransFeatureSize = par_upkfloat();
  c.minRotationalSeparation = par_upkfloat();
  c.minScaleSeparation = par_upkfloat();
  c.minTranslationalSeparation = par_upkfloat();
  for (i = 0; i < MAX_FRAC_FT_RES_LEVELS; ++i)
    c.fracRes[i] = par_upkfloat();
  c.maxRSCandidates = par_upkint();
  c.maxCandidates = par_upkint();

  c.minTheta = par_upkfloat();
  c.maxTheta = par_upkfloat();
  c.minScale = par_upkfloat();
  c.maxScale = par_upkfloat();

  c.minTX = par_upkfloat();
  c.maxTX = par_upkfloat();
  c.minTY = par_upkfloat();
  c.maxTY = par_upkfloat();

  c.update = par_upkint();
  c.partial = par_upkint();
  c.nWorkers = par_upkint();
}

void
PackPair (Pair *p)
{
  par_pkstr(p->imageName);
  par_pkint(p->imageMinX);
  par_pkint(p->imageMaxX);
  par_pkint(p->imageMinY);
  par_pkint(p->imageMaxY);
  par_pkstr(p->refName);
  par_pkint(p->refMinX);
  par_pkint(p->refMaxX);
  par_pkint(p->refMinY);
  par_pkint(p->refMaxY);
  par_pkstr(p->pairName);
}

void
UnpackPair (Pair *p)
{
  char s[PATH_MAX];

  if (p->imageName != NULL)
    free(p->imageName);
  par_upkstr(s);
  p->imageName = (char*) malloc(strlen(s)+1);
  strcpy(p->imageName, s);

  p->imageMinX = par_upkint();
  p->imageMaxX = par_upkint();
  p->imageMinY = par_upkint();
  p->imageMaxY = par_upkint();

  if (p->refName != NULL)
    free(p->refName);
  par_upkstr(s);
  p->refName = (char*) malloc(strlen(s)+1);
  strcpy(p->refName, s);

  p->refMinX = par_upkint();
  p->refMaxX = par_upkint();
  p->refMinY = par_upkint();
  p->refMaxY = par_upkint();

  if (p->pairName != NULL)
    free(p->pairName);
  par_upkstr(s);
  p->pairName = (char*) malloc(strlen(s)+1);
  strcpy(p->pairName, s);
}

void
PackTask ()
{
  PackPair(&(t.pair));
}

void
UnpackTask ()
{
  UnpackPair(&(t.pair));
}

void
PackResult ()
{
  PackPair(&(r.pair));
  par_pkint(r.updated);
  par_pkfloat(r.quality);
  par_pkfloat(r.rotation);
  par_pkfloat(r.scale);
  par_pkfloat(r.tx);
  par_pkfloat(r.ty);
  par_pkstr(r.message);
}

void
UnpackResult ()
{
  UnpackPair(&(r.pair));
  r.updated = par_upkint();
  r.quality = par_upkfloat();
  r.rotation = par_upkfloat();
  r.scale = par_upkfloat();
  r.tx = par_upkfloat();
  r.ty = par_upkfloat();
  par_upkstr(r.message);
}

size_t
CountBits (unsigned char *p, size_t n)
{
  size_t ii;
  unsigned char b;
  size_t sum = 0;

  for (ii = 0; ii < n; ++ii)
    {
      b = *p++;
      while (b != 0)
	{
	  if (b & 1)
	    ++sum;
	  b >>= 1;
	}
    }
  return(sum);
}

int
CreateDirectories (char *fn)
{
  int pos;
  char dn[PATH_MAX];
  char *slash;
  int len;
  struct stat statBuf;
  int hv;
  int i;

  for (pos = 0;
       (slash = strchr(&fn[pos], '/')) != NULL;
       pos = len + 1)
    {
      len = slash - fn;
      if (len == 0)
	continue;
      strncpy(dn, fn, len);
      dn[len] = '\0';

      /* check if in hash table */
      hv = 0;
      for (i = 0; i < len; ++i)
	hv = 239*hv + dn[i];
      hv &= DIR_HASH_SIZE-1;
      i = hv;
      while (dirHash[i] != NULL)
	{
	  if (strcmp(dirHash[i], dn) == 0)
	    break;
	  i = (i + 1) & (DIR_HASH_SIZE-1);
	  if (i == hv)
	    {
	      Log("Directory hash table is full!\n");
	      return(0);
	    }
	}
      if (dirHash[i] != NULL)
	continue;
      dirHash[i] = (char *) malloc(strlen(dn)+1);
      strcpy(dirHash[i], dn);

      if (stat(dn, &statBuf) == 0)
	{
	  if (S_ISDIR(statBuf.st_mode) ||
	      S_ISLNK(statBuf.st_mode))
	    continue;
	  Log("Output path component %s is not a directory\n", dn);
	  return(0);
	}
      if (errno != ENOENT)
	{
	  Log("Could not stat directory %s\n", dn);
	  return(0);
	}
      
      if (mkdir(dn, 0777) != 0)
	{
	  Log("Could not create directory %s\n", dn);
	  return(0);
	}
    }
  return(1);
}

void Error (char *fmt, ...)
{
  va_list args;
  char timestamp[32];

  if (logFile != NULL)
    {
      va_start(args, fmt);
      fprintf(logFile, "%s: ERROR: ", GetTimestamp(timestamp, 32));
      vfprintf(logFile, fmt, args);
      va_end(args);
      fflush(logFile);
    }
  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
  fflush(stderr);
  abort();
}

void Log (char *fmt, ...)
{
  va_list args;
  char timestamp[32];

  if (c.logBasename[0] == '\0')
    return;

  if (logFile == NULL)
    {
      int i;
      char logName[128];
      i = par_instance();
      if (i < 0)
	sprintf(logName, "%s0.log", c.logBasename);
      else if (c.nWorkers < 10)
	sprintf(logName, "%s%.1d.log", c.logBasename, i);
      else if (c.nWorkers < 100)
	sprintf(logName, "%s%.2d.log", c.logBasename, i);
      else if (c.nWorkers < 1000)
	sprintf(logName, "%s%.3d.log", c.logBasename, i);
      else
	sprintf(logName, "%s%.6d.log", c.logBasename, i);
      logFile = fopen(logName, "w");
      if (logFile == NULL)
	Error("Could not open log file %s\n", logName);
    }

  va_start(args, fmt);
  fprintf(logFile, "%s: ", GetTimestamp(timestamp, 32));
  vfprintf(logFile, fmt, args);
  va_end(args);
  fflush(logFile);
}

char *
GetTimestamp (char *timestamp, size_t size)
{
  struct timeval current_time;
  struct tm lt;
  size_t len;

  if (size < 24 ||
      gettimeofday(&current_time, NULL) != 0 ||
      strftime(timestamp, size, "%y%m%d %H:%M:%S",
	       localtime_r(&(current_time.tv_sec), &lt)) == 0)
    return(NULL);
  len = strlen(timestamp);
  if (size - len > 8)
    sprintf(&timestamp[len], ".%.6d", (int) current_time.tv_usec);
  return(timestamp);
}

  //  exp(x*i) = cos(x) + i sin(x)
  //
  //  for j=0:7                         
  //	  y(j+1)=x(j+1)*exp(-pi*i*j*j/16);  
  //      z(j+1)=exp(pi*i*(j-4)^2/16);      
  //  end
  //  for j=8:15                        
  //      y(j+1)=0;                         
  //      z(j+1)=exp(pi*i*(j-4-16)^2/16);       
  //  end
  //  w=fft(y).*fft(z);                
  //  iw=ifft(w);                      
  //  for k=0:7                        
  //      g(k+1)=exp(-pi*i*(k-4)^2/16)*iw(k+1);  
  //  end

  // >> g
  // g =
  // 1.1186 - 0.4470i   0.0858 - 0.5544i  -0.3448 + 0.3534i   0.7392 + 1.4457i
  // 2.2715 - 0.0000i   0.7392 - 1.4457i  -0.3448 - 0.3534i   0.0858 + 0.5544i
  // >> fftshift(fft(x))
  // ans =
  // 1.0769            -0.1471 - 0.0302i   1.1186 - 0.4470i  -0.3448 + 0.3534i
  // 2.2715            -0.3448 - 0.3534i   1.1186 + 0.4470i  -0.1471 + 0.0302i




// for j=0:7                     
//   y(j+1)=x(j+1)*exp(-pi*i*j*j/8);
//   z(j+1)=exp(pi*i*(j-4)^2/8);
// end
// for j=8:15
//   y(j+1)=0;
//   z(j+1)=exp(pi*i*(j-4-16)^2/8);
// end
// w=fft(y).*fft(z);
// iw=ifft(w);
// for k=0:7
//   g(k+1)=exp(-pi*i*(k-4)^2/8)*iw(k+1);
// end
//
// g =
// 1.0769 - 0.0000i  -0.1471 - 0.0302i   1.1186 - 0.4470i  -0.3448 + 0.3534i
// 2.2715 + 0.0000i  -0.3448 - 0.3534i   1.1186 + 0.4470i  -0.1471 + 0.0302i
// fftshift(fft(x))
// 1.0769            -0.1471 - 0.0302i   1.1186 - 0.4470i  -0.3448 + 0.3534i
// 2.2715            -0.3448 - 0.3534i   1.1186 + 0.4470i  -0.1471 + 0.0302i
