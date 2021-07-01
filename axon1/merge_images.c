#include <stdio.h>
#include <mpi.h>
#include <math.h>

#include "imio.h"
#include "dt.h"
#include "par.h"

typedef struct Context {
  /* NOTE: any new fields added to this struct should
     be also added to PackContext and UnpackContext */ 
  char type;
  char imageBasename[PATH_MAX];
  char maskBasename[PATH_MAX];
  char outputBasename[PATH_MAX];
  char logBasename[PATH_MAX];
} Context;

typedef struct Task {
  /* NOTE: any new fields added to this struct should
     be also added to PackTask and UnpackTask */ 
  char imageName[PATH_MAX];
} Task;

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
void Error (char *fmt, ...);
void Log (char *fmt, ...);


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

void
MasterTask (int argc,
            char **argv,
            char **envp)
{
  int i;
  int n;
  int len;
  char pairsFile[PATH_MAX];
  char outputPairsFile[PATH_MAX];
  char outputSortedPairsFile[PATH_MAX];
  int nPairs;
  Pair *pairs;
  int minSlice, maxSlice;
  int error;
  int pos;
  DIR *dir;
  struct dirent *de;
  char fn[PATH_MAX];
  FILE *f;
  char inputDirName[PATH_MAX];
  char inputPrefix[PATH_MAX];
  int inputPrefixLen;
  int minN, maxN;
  int digitLen;
  char maskDirName[PATH_MAX];
  char maskPrefix[PATH_MAX];
  int maskPrefixLen;
  int maskMinN, maskMaxN;
  char outputDirName[PATH_MAX];
  char tc;
  unsigned int w, h;
  int m;
  int xRes, yRes;
  int minS, maxS;
  int forward;
  struct stat statBuf;
  int pn;
  char imgn[2][PATH_MAX];
  char pairn[PATH_MAX];
  int imgMinX[2], imgMaxX[2], imgMinY[2], imgMaxY[2];
  unsigned int iw, ih;
  int imageSlice, refSlice;
  int imi;
  char line[LINE_LENGTH+1];
  FILE *opf;

  error = 0;
  c.type = '\0';
  c.imageBasename[0] = '\0';
  c.maskBasename[0] = '\0';
  c.discontinuityBasename[0] = '\0';
  c.strictMasking = 0;
  c.inputMapName[0] = '\0';
  c.cptsName[0] = '\0';
  c.constrainingMapName[0] = '\0';
  c.outputMapBasename[0] = '\0';
  c.outputWarpedBasename[0] = '\0';
  c.outputCorrelationBasename[0] = '\0';
  c.logBasename[0] = '\0';
  c.startLevel = -1;
  c.outputLevel = 0;
  c.minResolution = 1;
  c.distortion = 1.0;
  c.correspondence = 1.0;
  c.correspondenceThreshold = 0.0;
  c.constraining = 1.0;
  c.constrainingThreshold = 0.0;
  c.constrainingConfidenceThreshold = 0.5;
  c.quality = 10.0;
  c.minOverlap = 20.0;
  c.correlationHalfWidth = 31;
  c.writeAllMaps = 0;
  c.depth = 0;
  c.cptsMethod = -1;
  c.update = 0;
  c.partial = 0;
  c.nWorkers = par_workers();
  r.pair.imageName[0]  = r.pair.imageName[1] = NULL;
  r.pair.pairName = NULL;
  r.message = NULL;
  forward = 1;
  pairsFile[0] = '\0';
  outputPairsFile[0] = '\0';
  outputSortedPairsFile[0] = '\0';
  nPairs = 0;
  pairs = 0;
  minSlice = 0;
  maxSlice = 1000000000;
  memset(dirHash, 0, DIR_HASH_SIZE * sizeof(char*));

  for (i = 1; i < argc; ++i)
    if (strcmp(argv[i], "-depth") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d", &c.depth) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-min_res") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d", &c.minResolution) != 1)
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
	strcpy(c.maskBasename, argv[i]);
      }
    else if (strcmp(argv[i], "-input_map") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.inputMapName, argv[i]);
      }
    else if (strcmp(argv[i], "-cpts") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.cptsName, argv[i]);
      }
    else if (strcmp(argv[i], "-constraining_map") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.constrainingMapName, argv[i]);
      }
    else if (strcmp(argv[i], "-output") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.outputMapBasename, argv[i]);
      }
    else if (strcmp(argv[i], "-output_mask") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.outputMaskBasename, argv[i]);
      }
    else if (strcmp(argv[i], "-output_pairwise_mask") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.outputPairwiseMaskBasename, argv[i]);
      }
    else if (strcmp(argv[i], "-start_level") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d", &c.startLevel) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-output_level") == 0)
      {
	if (++i == argc || sscanf(argv[i], "%d", &c.outputLevel) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-output_warped") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.outputWarpedBasename, argv[i]);
      }
    else if (strcmp(argv[i], "-output_correlation") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(c.outputCorrelationBasename, argv[i]);
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
    else if (strcmp(argv[i], "-summary") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(summaryName, argv[i]);
      }
    else if (strcmp(argv[i], "-distortion") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%lf", &c.distortion) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-correspondence") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%lf", &c.correspondence) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-correspondence_threshold") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%lf", &c.correspondenceThreshold) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-constraining") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%lf", &c.constraining) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-constraining_threshold") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%lf", &c.constrainingThreshold) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-constraining_confidence_threshold") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%lf", &c.constrainingConfidenceThreshold) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-quality") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%lf", &c.quality) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-correlation_half_width") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%d", &c.correlationHalfWidth) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-min_overlap") == 0)
      {
	if (++i == argc ||
	    sscanf(argv[i], "%lf", &c.minOverlap) != 1)
	  {
	    error = 1;
	    break;
	  }
      }
    else if (strcmp(argv[i], "-all_maps") == 0)
      c.writeAllMaps = 1;
    else if (strcmp(argv[i], "-forward") == 0)
      forward = 1;
    else if (strcmp(argv[i], "-reverse") == 0)
      forward = 0;
    else if (strcmp(argv[i], "-update") == 0)
      c.update = 1;
    else if (strcmp(argv[i], "-partial") == 0)
      c.partial = 1;
    else if (strcmp(argv[i], "-strict_masking") == 0)
      c.strictMasking = 1;
    else if (strcmp(argv[i], "-tif") == 0)
      c.type = 't';
    else if (strcmp(argv[i], "-vis") == 0)
      vis = 1;
    else if (strcmp(argv[i], "-output_pairs") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(outputPairsFile, argv[i]);
      }
    else if (strcmp(argv[i], "-output_sorted_pairs") == 0)
      {
	if (++i == argc)
	  {
	    error = 1;
	    break;
	  }
	strcpy(outputSortedPairsFile, argv[i]);
      }
    else if (strcmp(argv[i], "-translation") == 0)
      {
	if (c.cptsMethod >= 0)
	  {
	    error = 1;
	    break;
	  }
	c.cptsMethod = TRANSLATION_METHOD;
      }
    else if (strcmp(argv[i], "-rigid") == 0)
      {
	if (c.cptsMethod >= 0)
	  {
	    error = 1;
	    break;
	  }
	c.cptsMethod = RIGID_METHOD;
      }
    else if (strcmp(argv[i], "-affine") == 0)
      {
	if (c.cptsMethod >= 0)
	  {
	    error = 1;
	    break;
	  }
	c.cptsMethod = AFFINE_METHOD;
      }
    else if (strcmp(argv[i], "-quadratic") == 0)
      {
	if (c.cptsMethod >= 0)
	  {
	    error = 1;
	    break;
	  }
	c.cptsMethod = QUADRATIC_METHOD;
      }
    else error = 1;

  if (error)
    {
      fprintf(stderr, "Usage: warp3 -input file_prefix -output file_prefix\n");
      fprintf(stderr, "              [-resolution HORIZxVERT]\n");
      fprintf(stderr, "              [-distortion distortion_weight]\n");
      fprintf(stderr, "              [-correspondence correspondence_weight]\n");
      fprintf(stderr, "              [-correspondence_threshold pixels]\n");
      fprintf(stderr, "              [-translation]\n");
      fprintf(stderr, "              [-rigid]\n");
      fprintf(stderr, "              [-affine]\n");
      fprintf(stderr, "              [-quadratic]\n");
      fprintf(stderr, "              [-quality quality_factor]\n");
      fprintf(stderr, "              [-min_res minimum_resolution_in_pixels]\n");
      fprintf(stderr, "              [-depth delta_depth]\n");
      fprintf(stderr, "              [-all_maps]\n");
      fprintf(stderr, "              [-update]\n");
      fprintf(stderr, "              [-partial]\n");
      fprintf(stderr, "              [-pairs <pair_file>]\n");
      fprintf(stderr, "              [-input_map <input_map_prefix>]\n");
      fprintf(stderr, "              [-constraining_map <constraining_map_prefix>]\n");
      fprintf(stderr, "              [-constraining constraining_weight]\n");
      fprintf(stderr, "              [-constraining_threshold pixels]\n");
      fprintf(stderr, "              [-constraining_confidence_threshold score]\n");
      fprintf(stderr, "              [-strict_masking]\n");
      fprintf(stderr, "              [-min_overlap percent]\n");
      fprintf(stderr, "              [-output_pairs <output_pair_file>]\n");
      fprintf(stderr, "              [-output_sorted_pairs <output_pair_file>]\n");
      fprintf(stderr, "              [-logs <log_file_directory>]\n");
      fprintf(stderr, "   where ranges are expressed as: integer\n");
      fprintf(stderr, "                              or: integer-integer\n");
      fprintf(stderr, "                              or: integer-\n");
      fprintf(stderr, "                              or: -integer\n");
      exit(1);
    }

  Log("MASTER starting on node %d\n", par_instance());

  Log("argc = %d\n", argc);
  for (i = 0; i < argc; ++i)
    Log("ARGV[%d] = %s\n", i, argv[i]);

  /* check that at least minimal parameters were supplied */
  if (c.imageBasename[0] == '\0' || c.outputMapBasename[0] == '\0' ||
      pairsFile[0] == '\0')
    Error("-input, -output, and -pairs parameters must be specified.\n");
  if (c.cptsMethod < 0)
    c.cptsMethod = AFFINE_METHOD;

  f = fopen(pairsFile, "r");
  if (f == NULL)
    Error("Could not open pairs file %s\n", pairsFile);
  while (fgets(line, LINE_LENGTH, f) != NULL)
    {
      if (line[0] == '\0' || line[0] == '#')
	continue;
      if (sscanf(line, "%s %d %d %d %d %s %d %d %d %d %s",
		 imgn[0], &imgMinX[0], &imgMaxX[0], &imgMinY[0], &imgMaxY[0],
		 imgn[1], &imgMinX[1], &imgMaxX[1], &imgMinY[1], &imgMaxY[1],
		 pairn) != 11)
	{
	  if (sscanf(line, "%s %s %s", imgn[0], imgn[1], pairn) != 3)
	    Error("Invalid line in pairs file %s:\n%s\n", pairsFile, line);
	  imgMinX[0] = -1;
	  imgMaxX[0] = -1;
	  imgMinY[0] = -1;
	  imgMaxY[0] = -1;
	  imgMinX[1] = -1;
	  imgMaxX[1] = -1;
	  imgMinY[1] = -1;
	  imgMaxY[1] = -1;
	}
      if ((nPairs & 1023) == 0)
	pairs = (Pair *) realloc(pairs, (nPairs + 1024) * sizeof(Pair));
      for (imi = 0; imi < 2; ++imi)
	{
	  pairs[nPairs].imageName[imi] = NULL;
	  CopyString(&(pairs[nPairs].imageName[imi]), imgn[imi]);
	  pairs[nPairs].imageMinX[imi] = imgMinX[imi];
	  pairs[nPairs].imageMaxX[imi] = imgMaxX[imi];
	  pairs[nPairs].imageMinY[imi] = imgMinY[imi];
	  pairs[nPairs].imageMaxY[imi] = imgMaxY[imi];
	}
      pairs[nPairs].pairName = NULL;
      CopyString(&(pairs[nPairs].pairName), pairn);
      ++nPairs;
    }
  fclose(f);

  /* check that output directories are writeable */ 
  Log("MASTER checking output directories\n");
  sprintf(fn, "%sTEST.map", c.outputMapBasename);
  if (!CreateDirectories(fn))
    Error("Could not create directory for output maps: %s\n",
	  c.outputMapBasename);
  f = fopen(fn, "w");
  if (f == NULL)
    {
      for (i = strlen(c.outputMapBasename)-1; i >= 0 && c.outputMapBasename[i] != '/'; --i) ;
      if (i != 0)
	strncpy(outputDirName, c.outputMapBasename, i+1);
      else
	outputDirName[0] = '.';
      outputDirName[i+1] = '\0';
      Error("Could not open test output file %s --\n       does directory %s exist and is it writeable?\n", fn, outputDirName);
    }
  fclose(f);
  unlink(fn);

  if (c.outputWarpedBasename[0] != '\0')
    {
      sprintf(fn, "%sTEST.pgm", c.outputWarpedBasename);
      if (!CreateDirectories(fn))
	Error("Could not create directory for output warped images: %s\n",
	      c.outputWarpedBasename);
      f = fopen(fn, "w");
      if (f == NULL)
	{
	  for (i = strlen(c.outputWarpedBasename)-1; i >= 0 && c.outputWarpedBasename[i] != '/'; --i) ;
	  if (i != 0)
	    strncpy(outputDirName, c.outputWarpedBasename, i+1);
	  else
	    outputDirName[0] = '.';
	  outputDirName[i+1] = '\0';
	  Error("Could not open test output image file %s --\n      does directory %s exist and is it writeable?\n", fn, outputDirName);
	}
      fclose(f);
      unlink(fn);
    }

  if (c.outputCorrelationBasename[0] != '\0')
    {
      sprintf(fn, "%sTEST.pgm", c.outputCorrelationBasename);
      if (!CreateDirectories(fn))
	Error("Could not create directory for output correlation: %s\n",
	      c.outputCorrelationBasename);
      f = fopen(fn, "w");
      if (f == NULL)
	{
	  for (i = strlen(c.outputCorrelationBasename)-1; i >= 0 && c.outputCorrelationBasename[i] != '/'; --i) ;
	  if (i != 0)
	    strncpy(outputDirName, c.outputCorrelationBasename, i+1);
	  else
	    outputDirName[0] = '.';
	  outputDirName[i+1] = '\0';
	  Error("Could not open test output image file %s --\n      does directory %s exist and is it writeable?\n", fn, outputDirName);
	}
      fclose(f);
      unlink(fn);
    }

  if (outputPairsFile[0] != '\0')
    {
      if (!CreateDirectories(outputPairsFile))
	Error("Could not create directory for output pairs file: %s\n",
	      outputPairsFile);
      f = fopen(outputPairsFile, "w");
      if (f == NULL)
	Error("Could not open output pairs file %s\n", outputPairsFile);
      fclose(f);
      unlink(outputPairsFile);
    }

  if (outputSortedPairsFile[0] != '\0')
    {
      if (!CreateDirectories(outputSortedPairsFile))
	Error("Could not create directory for output sorted pairs file: %s\n",
	      outputSortedPairsFile);
      f = fopen(outputSortedPairsFile, "w");
      if (f == NULL)
	Error("Could not open output sorted pairs file %s\n",
	      outputSortedPairsFile);
      fclose(f);
      unlink(outputSortedPairsFile);
    }

  if (summaryName[0] != '\0')
    {
      if (!CreateDirectories(summaryName))
	Error("Could not create directory for summary file: %s\n",
	      summaryName);
      f = fopen(summaryName, "w");
      if (f == NULL)
	Error("Could not open summary file %s\n",
	      outputSortedPairsFile);
      fclose(f);
      unlink(summaryName);
    }

  if (c.logBasename[0] != '\0')
    {
      sprintf(fn, "%sTEST.log", c.logBasename);
      if (!CreateDirectories(fn))
	Error("Could not create directory for test logfile: %s\n", fn);
      f = fopen(fn, "w");
      if (f == NULL)
	Error("Could not open test log file %s\n", fn);
      fclose(f);
      unlink(fn);
    }

  Log("MASTER setting context\n");

  par_set_context();

  Log("MASTER set context\n");

  /* for all slices */
  printf("Processing slices: ");
  fflush(stdout);
  for (imi = 0; imi < 2; ++imi)
    t.pair.imageMinX[imi] = t.pair.imageMaxX[imi] = t.pair.imageMinY[imi] = t.pair.imageMaxY[imi] = -1;
  Log("nPairs = %d\n", nPairs);
  for (pn = 0; pn < nPairs; ++pn)
    {
      for (imi = 0; imi < 2; ++imi)
	{
	  CopyString(&(t.pair.imageName[imi]), pairs[pn].imageName[imi]);
	  t.pair.imageMinX[imi] = pairs[pn].imageMinX[imi];
	  t.pair.imageMaxX[imi] = pairs[pn].imageMaxX[imi];
	  t.pair.imageMinY[imi] = pairs[pn].imageMinY[imi];
	  t.pair.imageMaxY[imi] = pairs[pn].imageMaxY[imi];
	}
      CopyString(&(t.pair.pairName), pairs[pn].pairName);

      // make sure that output directories exist
      sprintf(fn, "%s%s.map", c.outputMapBasename, t.pair.pairName);
      if (!CreateDirectories(fn))
	continue;
      if (c.outputWarpedBasename[0] != '\0')
	{
	  sprintf(fn, "%s%s.pgm", c.outputWarpedBasename, t.pair.pairName);
	  if (!CreateDirectories(fn))
	    continue;
	}
      if (c.outputCorrelationBasename[0] != '\0')
	{
	  sprintf(fn, "%s%s.pgm", c.outputCorrelationBasename, t.pair.pairName);
	  if (!CreateDirectories(fn))
	    continue;
	}

      Log("Delegating pair %d\n", pn);
      par_delegate_task();
    }
  par_finish();

  if (outputPairsFile[0] != '\0')
    {
      qsort(results, nResults, sizeof(Result), SortBySlice);
      opf = fopen(outputPairsFile, "w");
      if (opf == NULL)
	Error("Could not open output pairs file %s for writing.\n",
	      outputPairsFile);
      for (i = 0; i < nResults; ++i)
	if (results[i].updated)
	  fprintf(opf, "%s %d %d %d %d %s %d %d %d %d %s\n",
		  results[i].pair.imageName[0],
		  results[i].pair.imageMinX[0],
		  results[i].pair.imageMaxX[0],
		  results[i].pair.imageMinY[0],
		  results[i].pair.imageMaxX[0],
		  results[i].pair.imageName[1],
		  results[i].pair.imageMinX[1],
		  results[i].pair.imageMaxX[1],
		  results[i].pair.imageMinY[1],
		  results[i].pair.imageMaxX[1],
		  results[i].pair.pairName);
      fclose(opf);
    }
  if (outputSortedPairsFile[0] != '\0')
    {
      qsort(results, nResults, sizeof(Result), SortByEnergy);
      opf = fopen(outputSortedPairsFile, "w");
      if (opf == NULL)
	Error("Could not open output sorted pairs file %s for writing.\n",
	      outputSortedPairsFile);
      for (i = 0; i < nResults; ++i)
	if (results[i].updated)
	  fprintf(opf, "%s %d %d %d %d %s %d %d %d %d %s\n",
		  results[i].pair.imageName[0],
		  results[i].pair.imageMinX[0],
		  results[i].pair.imageMaxX[0],
		  results[i].pair.imageMinY[0],
		  results[i].pair.imageMaxX[0],
		  results[i].pair.imageName[1],
		  results[i].pair.imageMinX[1],
		  results[i].pair.imageMaxX[1],
		  results[i].pair.imageMinY[1],
		  results[i].pair.imageMaxX[1],
		  results[i].pair.pairName);
      fclose(opf);
    }

  if (summaryName[0] != '\0')
    {
      qsort(results, nResults, sizeof(Result), SortBySlice);
      f = fopen(summaryName, "w");
      if (f == NULL)
	Error("Could not open summary output file %s\n", summaryName);

      fprintf(f, "Sorted by slice:\n");
      fprintf(f, "IMAGE    REFERENCE  CORRELATION   DISTORTION    CORRESPOND    CONSTRAIN     ENERGY\n");
      for (i = 0; i < nResults; ++i)
	fprintf(f, "%8s %8s %12.6f  %12.6f  %12.6f  %12.6f  %12.6f\n",
		results[i].pair.imageName[0],
		results[i].pair.imageName[1],
		results[i].correlation,
		results[i].distortion,
		results[i].correspondence,
		results[i].constraining,
		results[i].distortion * c.distortion - results[i].correlation +
		results[i].correspondence * c.correspondence + results[i].constraining * c.constraining);
      fprintf(f, "\n\nSorted by energy:\n");
      fprintf(f, "IMAGE    REFERENCE  CORRELATION   DISTORTION    CORRESPOND    CONSTRAIN     ENERGY\n");
      qsort(results, nResults, sizeof(Result), SortByEnergy);
      for (i = 0; i < nResults; ++i)
	fprintf(f, "%8s %8s  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f\n",
		results[i].pair.imageName[0],
		results[i].pair.imageName[1],
		results[i].correlation,
		results[i].distortion,
		results[i].correspondence,
		results[i].constraining,
		results[i].distortion * c.distortion - results[i].correlation +
		results[i].correspondence * c.correspondence + results[i].constraining * c.constraining);
      fclose(f);
    }

  printf(" %d\nAll slices completed.\n", nResults);
}

void
MasterResult ()
{
  int imi;

  if (r.message != NULL)
    Error("\nThe following error was encountered by one of the worker processes:%s\n", r.message);

  results = (Result*) realloc(results, (nResults + 1) * sizeof(Result));
  for (imi = 0; imi < 2; ++imi)
    {
      results[nResults].pair.imageName[imi] = NULL;
      CopyString(&(results[nResults].pair.imageName[imi]), r.pair.imageName[imi]);
      results[nResults].pair.imageMinX[imi] = r.pair.imageMinX[imi];
      results[nResults].pair.imageMaxX[imi] = r.pair.imageMaxX[imi];
      results[nResults].pair.imageMinY[imi] = r.pair.imageMinY[imi];
      results[nResults].pair.imageMaxY[imi] = r.pair.imageMaxY[imi];
    }
  results[nResults].pair.pairName = NULL;
  CopyString(&(results[nResults].pair.pairName), r.pair.pairName);
  results[nResults].updated = r.updated;
  results[nResults].distortion = r.distortion;
  results[nResults].correlation = r.correlation;
  results[nResults].correspondence = r.correspondence;
  results[nResults].constraining = r.constraining;
  results[nResults].message = NULL;
  CopyString(&(results[nResults].message), r.message);

  if ((nResults % 50) == 0 && nResults != 0)
    printf(" %d \n                   ", nResults);
  printf(".");
  ++nResults;
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
