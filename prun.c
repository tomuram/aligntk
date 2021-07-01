#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <limits.h>
#include <unistd.h>
#include <mpi.h>
#include <sys/types.h>
#include <sys/wait.h>

#define MAX_COMMAND_LENGTH	4095

FILE *logFile = NULL;

void Error (char *fmt, ...);
void Log (char *fmt, ...);

int
main (int argc, char **argv)
{
  int p, np;
  int i;
  char command[MAX_COMMAND_LENGTH+1];
  char line[MAX_COMMAND_LENGTH+1];
  char cmd[MAX_COMMAND_LENGTH+1];
  char commandName[PATH_MAX];
  char listName[PATH_MAX];
  char logName[PATH_MAX];
  char fn[PATH_MAX];
  char hostName[PATH_MAX];
  int error;
  int n;
  FILE *f;
  FILE *pf;
  int len;
  int status;

  /* initialize MPI */
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    {
      fprintf(stderr, "Could not do MPI_Init\n");
      exit(1);
    }
  if (MPI_Comm_size(MPI_COMM_WORLD, &np) != MPI_SUCCESS)
    {
      fprintf(stderr, "Could not do MPI_Comm_size\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
      exit(1);
    }
  if (MPI_Comm_rank(MPI_COMM_WORLD, &p) != MPI_SUCCESS)
    {
      fprintf(stderr, "Could not do MPI_Comm_rank\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
      exit(1);
    }
      
  if (p == 0)
    {
      error = 0;
      listName[0] = '\0';
      logName[0] = '\0';
      commandName[0] = '\0';
      for (i = 1; i < argc; ++i)
	if (strcmp(argv[i], "-command") == 0)
	  {
	    if (++i == argc)
	      {
		error = 1;
		break;
	      }
	    strcpy(commandName, argv[i]);
	  }
	else if (strcmp(argv[i], "-list") == 0)
	  {
	    if (++i == argc)
	      {
		error = 1;
		break;
	      }
	    strcpy(listName, argv[i]);
	  }
	else if (strcmp(argv[i], "-logs") == 0)
	  {
	    if (++i == argc)
	      {
		error = 1;
		break;
	      }
	    strcpy(logName, argv[i]);
	  }
	else
	  error = 1;

      if (error || listName[0] == '\0')
	{
	  fprintf(stderr, "Usage: prun -list list_file\n");
	  fprintf(stderr, "            [-command command_file]\n");
	  fprintf(stderr, "            [-logs log_file_prefix]\n");
	  MPI_Abort(MPI_COMM_WORLD, 1);
	  exit(1);
	}
    }

  if (logName[0] == '\0')
    strcpy(logName, "logs/");
  if (MPI_Bcast(logName, PATH_MAX, MPI_CHAR,
		0, MPI_COMM_WORLD) != MPI_SUCCESS)
    {
      fprintf(stderr, "MPI_Bcast failed.\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
      exit(1);
    }
  sprintf(fn, "%s%.2d.log", logName, p);
  logFile = fopen(fn, "w");
  if (logFile == NULL)
    {
      fprintf(stderr, "Could not open log file %s\n", fn);
      MPI_Abort(MPI_COMM_WORLD, 1);
      exit(1);
    }

  gethostname(hostName, PATH_MAX);
  if (p == 0)
    {
      Log("prun master starting on node %s\n", hostName);
      if (commandName[0] != '\0')
	{
	  f = fopen(commandName, "r");
	  if (f == NULL)
	    Error("Could not open command file %s\n", commandName);
	  if (fgets(command, MAX_COMMAND_LENGTH+1, f) == NULL)
	    Error("Could not read command from file %s\n", commandName);
	  len = strlen(command);
	  if (len > 0 && command[len-1] == '\n')
	    command[len-1] = '\0';
	  if (command[0] == '\0')
	    Error("No command found in file %s\n", commandName);
	  fclose(f);
	}
      else
	strcpy(command, "%s");

      f = fopen(listName, "r");
      if (f == NULL)
	Error("Could not open list file %s\n", listName);
      while (fgets(line, MAX_COMMAND_LENGTH+1, f) != NULL)
	{
	  if (line[0] == '#')
	    continue;

	  len = strlen(line);
	  if (len > 0 && line[len-1] == '\n')
	    line[len-1] = '\0';
	  if (line[0] == '\0')
	    continue;

	  /* construct command to execute */
	  snprintf(cmd, MAX_COMMAND_LENGTH+1, command,
		   line, line, line, line,
		   line, line, line, line,
		   line, line, line, line,
		   line, line, line, line);
	  cmd[MAX_COMMAND_LENGTH] = '\0';
	  if (cmd[0] == '\0')
	    continue;

	  /* wait for a ready worker */
	  if (MPI_Recv(&n, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
		       MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS)
	    Error("MPI_Recv failed\n");
	  if (n < 1 || n >= np)
	    Error("Received out-of-range worker index (%d)\n", n);
	  Log("Delegating to worker %d: %s\n", n, cmd);

	  /* send command to that worker */
	  if (MPI_Send(cmd, MAX_COMMAND_LENGTH+1, MPI_CHAR,
		       n, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
	    Error("MPI_Send to worker %d failed\n", n);
	}
      fclose(f);

      /* wait for all workers to finish */
      Log("All commands delegated; waiting for completion\n");
      cmd[0] = '\0';
      for (i = 1; i < np; ++i)
	{
	  if (MPI_Recv(&n, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
		       MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS)
	    Error("Final MPI_Recv failed\n");
	  if (MPI_Send(cmd, MAX_COMMAND_LENGTH+1, MPI_CHAR,
		       n, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
	    Error("Final MPI_Send failed\n");
	}
      Log("All workers finished.\n");
    }
  else
    {
      Log("prun worker starting on node %s\n", hostName);
      for (;;)
	{
	  /* send a ready message */
	  if (MPI_Send(&p, 1, MPI_INT, 0, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
	    Error("MPI_Send of ready worker failed.\n");

	  /* receive a command */
	  if (MPI_Recv(cmd, MAX_COMMAND_LENGTH+1, MPI_CHAR, 0, 0,
		       MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS)
	    Error("MPI_Recv in worker failed.\n");
	  if (cmd[0] == '\0')
	    {
	      Log("Worker %d received request to finalize.\n", p);
	      break;
	    }

	  /* execute the command, piping stdout to log file */
	  Log("Worker %d running command %s\n", p, cmd);
	  pf = popen(cmd, "r");
	  if (pf == NULL)
	    Error("popen(%s) failed.\n", cmd);
	  while (fgets(line, MAX_COMMAND_LENGTH+1, pf) != NULL)
	    Log("%s", line);
	  status = pclose(pf);
	  Log("Worker %d finished running command %s (exit code = %d).\n\n", p, cmd, WEXITSTATUS(status));
	  //	  if (WEXITSTATUS(status) != 0)
	  //	    {
	  //	      Log("prun worker terminating because of non-zero exit code.\n");
	  //	      break;
	  //	    }
	}
    }

  fclose(logFile);
  MPI_Finalize();
  return(0);
}

void Error (char *fmt, ...)
{
  va_list args;

  if (logFile != NULL)
    {
      va_start(args, fmt);
      fprintf(logFile, "%f: ERROR: ", MPI_Wtime());
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

  va_start(args, fmt);
  fprintf(logFile, "%f: ", MPI_Wtime());
  vfprintf(logFile, fmt, args);
  va_end(args);
  fflush(logFile);
}
