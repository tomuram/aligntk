/*
 *	libpar.c	Parallel task manager
 *
 *	This module implements a single-master/multiple-slave style
 *	of parallel processing.  Both the master and worker programs
 *	(which typically are instances of the same code) should call
 *	ParallelProcess upon startup.  In this call they supply the
 *	routines to be used for task execution as well as communication.
 *
 *	The user-supplied (*master_task) procedure farms out individual
 *	tasks to workers by setting up a Task record and then calling
 *	par_delegate_task.  When a worker has completed a task, it sets up
 *	a Result record before returning from (*worker_task).
 *	(*master_result) is then automatically called on the master
 *	to process the values in the Result record.
 *
 *  Copyright (c) 1995-2011 National Resource for Biomedical
 *                          Supercomputing,
 *                          Pittsburgh Supercomputing Center,
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
 *       NIH NCRR grant 5P41RR006009.
 *
 *  HISTORY
 *	11/95 - Written by Greg Hood
 *	6/96  - Added capability for workers to return results to master
 *              (ghood@psc.edu)
 *	12/96 - Added Wait call and made group names passed through environment
 *		variables rather than by arguments (ghood@psc.edu)
 *	11/99 - added par_broadcast_context() for sending large contexts
 *		(ghood@psc.edu)
 *	3/00  - fixed bug in HandleRestart (ghood@psc.edu)
 *      3/01  - Converted to use either PVM or MPI (ghood@psc.edu)
 *      2/04  - Fixed bug that caused ready/idle list corruption upon
 *              worker termination (SWang@psych.uic.edu, ghood@psc.edu,
 *              welling@stat.cmu.edu)
 *      11/09 - Removed PVM support; made MPI support more
 *              robust (ghood@psc.edu); removed FIASCO dependencies
 *
 */

/*
 * POSSIBLE FUTURE OPTIMIZATIONS:
 *	The scheduling could be made more sophisticated so that tasks having
 *	    the same context would be preferentially assigned to the same
 *	    worker, in order to reduce context switching.
 *	The scheduling could be made predictive so that a task would not
 *	    be assigned to a historically slow worker near the end of the
 *	    run when a faster worker will likely soon be available.
 */

#include <sys/types.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <netinet/in.h>
#include <netdb.h>
#include <time.h>
#include <stdarg.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "par.h"

#define WORKER_RECEIVE_TIMEOUT	600000	/* seconds before a worker will timeout
					   waiting for a message from the
					   master */
#define PENDING_WORKERS_TIMEOUT	300	/* seconds to wait for pending
					   workers to show up before exiting */
#define REPORT_INTERVAL		5	/* we will report every this
					   many seconds when waiting
					   for pending workers to show up */
#define MAX_HOSTNAME_LENGTH	128     /* maximum # of chars in host name */

#define REQUEST_MSG		1	/* worker -> master */
#define RESTART_MSG		2	/* worker -> master */
#define CONTEXT_MSG		3	/* master -> worker  */
#define TASK_MSG		4	/* master -> worker  */
#define TERMINATE_MSG		6	/* master -> worker  */
#define WORKER_EXIT_MSG		7	/* PVM daemon -> master */
#define BROADCAST_CONTEXT_MSG	10	/* master -> worker */
#define BROADCAST_ACK_MSG	11	/* worker -> master */

#ifndef FALSE
#define FALSE			0
#endif
#ifndef TRUE
#define TRUE			1
#endif

#undef MIN
#define MIN(x,y)	((x) < (y) ? (x) : (y))
#undef MAX
#define MAX(x,y)	((x) > (y) ? (x) : (y))

typedef int Boolean;

typedef char Hostname[MAX_HOSTNAME_LENGTH];

typedef struct Buffer {
  unsigned char *buffer;	/* actual data in the buffer */
  int size;			/* size of the buffer in bytes */
  int position;			/* reading or writing position */
} Buffer;

typedef struct Context {
  int use_count;		/* # of times this structure is referenced */
  Buffer buffer;		/* the context buffer */
} Context;

typedef struct Task {
  struct Task *next;		/* next on the task queue */
  struct Task *prev;		/* previous on the task queue */
  Context *context;	        /* context associated with the task */
  int number;			/* number of this task */
  Buffer buffer;		/* the task buffer */
} Task;

typedef struct WorkerState {
  int tid;			/* the TID of this worker; note that for
				   MPI-based compilations, this is equal
				   to the rank */
  Hostname host;		/* the name of the host this worker
				   is running on */
  Boolean idle;			/* TRUE if this worker is ready to accept
				   a task and has no tasks currently
				   assigned to it */
  int next_idle;		/* the next worker on the idle list
				   (-1 terminates list) */
  Context* last_context;	/* the context that was last sent to the
				   worker */
  Task *first_task;		/* the first task (of at most 2) assigned to
				   this worker */
} WorkerState;

static Context *current_context = NULL;

static Task *first_queued_task = NULL;
static Task *last_queued_task = NULL;

static WorkerState workers[PAR_MAX_WORKERS];

static Boolean par = TRUE;	/* TRUE if we are going to run in parallel */

static int prog_argc;		/* the # of arguments this program was
				   invoked with, less -PAR_ arguments */
static char **prog_argv;	/* the arguments this program was invoked
				   with, less -PAR_ arguments */
static int par_argc;            /* the # of -PAR_ arguments this program
				   was invoked with */
static char **par_argv;         /* the -PAR_ arguments this program was
				   invoked with */
static char prog_name[128];	/* our program name */

static int rank = -1;		/* the rank (instance #) of this process within
				   the process group (0 indicates this is the
				   master) */

static int my_tid;		/* the TID of this process */
static int master_tid;		/* the TID of the master process */
static Hostname my_host_name;	/* the name of the host on which this process
				   is running */

static Boolean context_changed = FALSE;	/* TRUE if the context has changed
					   since the last call to
					   par_delegate_task */

static Par_Task task_number = 0; /* the next task number to be assigned */
static int tasks_outstanding = 0;/* counts the # of tasks that have been
				    delegated but not yet finished by the
				    workers */

static int n_workers = 0;	/* actual number of workers that are running */
static int n_workers_pending = 0; /* number of workers that should be starting,
				     but we haven't heard from yet */
static int idle_workers = -1;	/* a list of idle workers (i.e., those workers
				   having 0 tasks currently assigned to them);
				   a -1 terminates this list */

static int par_verbose = FALSE;	/* if TRUE, report on normal events such as
				   task delegation and receiving worker
				   results */

static unsigned char *task_completed = NULL;
                                /* a bit array that tells for each task
				   whether it has been completed (a 1
				   indicated it has) */
static int task_completed_size = 0;
                                /* the number of bytes in task_completed */
static int task_completed_first = 0;
                                /* the number of the first task represented
				   in the table; all tasks with ID's lower
				   than this have already completed */
static int task_completed_offset = 0;
                                /* the index into task_completed at which
				   the task_completed_first task can be
				   found */
static int broadcast_count = 0;	/* number of broadcast messsages that
				   have been sent */
static int broadcast_first_ack_count = -1;
                                /* count of the last acknowledge received from
				   the first process that this process relayed
				   the broadcast to */
static int broadcast_second_ack_count = -1;
                                /* count of the last acknowledge received from
				   the second process that this process relayed
				   the the broadcast to */
static int broadcast_first_forward_tid = -1; /* TID of the first process to
						forward a broadcast to */
static int broadcast_second_forward_tid = -1; /* TID of the second process to
						 forward a broadcast to */
static int broadcast_backward_tid = -1; /* TID of the process from which this
					   process receives broadcasts */
static int n_processes = 1;	/* number of entries in process_tids */
static int process_tids[PAR_MAX_WORKERS+1];
                                /* stores all TIDs for processes
				   involved in a broadcast */
static int my_process_index = 0; /* index into process_tids for my process */

static int n_hosts_pending = 0;	/* number of hosts in the host table which
				   for which max_workers > 0 but which we
				   have not yet heard from */

static struct timeval longest_task= {0,0}; 
                                /* run time of longest task so far */

static unsigned char *out_buffer = NULL;  /* output buffer */
static int out_size = 0;	/* size of output buffer in bytes */
static int out_position = 0;	/* index into output buffer */

static unsigned char *in_buffer = NULL;
static int in_size = 0;		/* size of input buffer in bytes */
static int in_position = 0;	/* index into input buffer */

static int sizeof_byte;		/* size of packed length of individual types */
static int sizeof_short;
static int sizeof_int;
static int sizeof_long;
static int sizeof_float;
static int sizeof_double;


/* the functions which must be supplied by the user of the this module */
static void (*par_master_task)() = NULL;	/* required */
static void (*par_master_result)() = NULL;	/* optional */
static void (*par_worker_context)() = NULL;	/* optional */
static void (*par_worker_task)() = NULL;	/* required */
static void (*par_worker_finalize)() = NULL;	/* optional */
static void (*par_pack_context)() = NULL;	/* optional */
static void (*par_unpack_context)() = NULL;	/* optional */
static void (*par_pack_task)() = NULL;		/* optional */
static void (*par_unpack_task)() = NULL;	/* optional */
static void (*par_pack_result)() = NULL;	/* optional */
static void (*par_unpack_result)() = NULL;	/* optional */

/* the internal functions */
static void ScanEnvironment();
static void DispatchTasks();
static void DispatchTask();
static int  FindReadyWorker();
static void HandleMessage();
static void HandleRequest();
static void HandleBroadcastAck();
static void HandleWorkerExit();
static void QueueTask();
static void RequeueTask();
static Context *ReuseContext();
static void DisuseContext();
static void FreeBuffer ();
static void PutOnIdleList();
static void RemoveFromIdleList();
static void DeclareWorkerDead();
static void PerformWorkerTasks();
static void ComposeRequest();
static int  MasterReceiveMessage(float timeout);
static int  WorkerReceiveMessage();
static void PrepareToSend();
static void Send();

static void StartMPI();
static void ExpandInBuffer(int size);
static void ExpandOutBuffer(int size);

static Boolean StringCopy (char *to, const char *from, const int maxChars);

static void Report (char *fmt, ...);
static void Error (char *fmt, ...);
static void Abort (char *fmt, ...);

void
par_process (int argc,
	     char **argv,
	     char **envp,
	     void (*master_task)(),
	     void (*master_result)(),
	     void (*worker_context)(),
	     void (*worker_task)(),
	     void (*worker_finalize)(),
	     void (*pack_context)(),
	     void (*unpack_context)(),
	     void (*pack_task)(),
	     void (*unpack_task)(),
	     void (*pack_result)(),
	     void (*unpack_result)())
{
  struct hostent *e;

  /* save the user-supplied functions */
  if (master_task == NULL)
    Abort("libpar: master_task cannot be NULL\n");
  par_master_task = master_task;
  par_master_result = master_result;
  par_worker_context = worker_context;
  if (worker_task == NULL)
    Abort("libpar: worker_task cannot be NULL\n");
  par_worker_task = worker_task;
  par_worker_finalize = worker_finalize;
  par_pack_context = pack_context;
  par_unpack_context = unpack_context;
  par_pack_task = pack_task;
  par_unpack_task = unpack_task;
  par_pack_result = pack_result;
  par_unpack_result = unpack_result;

  if (gethostname(my_host_name, MAX_HOSTNAME_LENGTH) < 0)
    Abort("gethostname failed\n");
  e = gethostbyname(my_host_name);
  if (e != NULL)
    StringCopy(my_host_name, e->h_name, MAX_HOSTNAME_LENGTH);
  else
    strcpy(my_host_name, "unknown_host");

  /* we have to enroll with MPI here in order for the command-line
     arguments to be valid */

  //  printf("Going to call MPI_Init with argc = %d and argv[1] = %s\n",
  //	 argc, argv[1]);
  if (MPI_Init(&argc, &argv) == MPI_SUCCESS &&
      MPI_Comm_size(MPI_COMM_WORLD, &n_workers_pending) == MPI_SUCCESS &&
      n_workers_pending > 1)
    {
      printf("running in parallel with MPI\n");
      n_workers_pending = 0;
      par = TRUE;
    }
  else
    printf("not running in parallel\n");

  /* check the relevant environment variables */
  ScanEnvironment(argc, argv);

  if (par) 
    {
      /* the parallelism flag is on, so start running in an MPI mode */
      StartMPI(argc, argv, envp);
    }

  if (!par) /* par may possibly have been turned off by StartMPI */
    {
      /* the parallelism flag is off, so just run the master task */
      (*master_task)(prog_argc, prog_argv, envp);
      if (worker_finalize != NULL)
	(*worker_finalize)();
    }

  if (par_verbose) Report("Tid %d at MPI_Finalize\n",my_tid);
  if (MPI_Finalize() != MPI_SUCCESS) {
    Abort("MPI_Finalize failed\n");
  }
}

static void
ScanEnvironment (int argc, char **argv)
{
  char *p;
  int i;
  int v;
  char s[256];

  if ((p = strrchr(argv[0], '/')) == NULL)
    StringCopy(prog_name, argv[0], 128);
  else
    StringCopy(prog_name, p+1, 128);

  if ((p = getenv("PAR_VERBOSE")) != NULL &&
      sscanf(p, "%d", &v) == 1)
    par_verbose = (v != 0);
  prog_argc = 0;
  prog_argv = (char **) malloc(argc * sizeof(char *));
  par_argc = 0;
  par_argv = (char **) malloc(argc * sizeof(char *));
  for (i = 0; i < argc; ++i)
    {
      par_argv[par_argc++] = argv[i];
      if (strncmp(argv[i], "-PAR_", 5) == 0) {
	if (strncmp(argv[i], "-PAR_ENABLE=", 12) == 0)
	  {
	    if (sscanf(&argv[i][12], "%d", &v) == 1)
	      par = v != 0;
	  }
	else if (strncmp(argv[i], "-PAR_CWD=", 9) == 0)
	  {
	    if (sscanf(&argv[i][9], "%s", s) == 1)
	      {
		if (chdir(s) != 0)
		  Error("Could not change to working directory %s\n", s);
	      }
	  }
	else if (strncmp(argv[i], "-PAR_VERBOSE=", 13) == 0)
	  {
	    if (sscanf(&argv[i][13], "%d", &v) == 1)
	      par_verbose = (v != 0);
	  }
      }
      else {
	prog_argv[prog_argc++]= argv[i];
      }
    }
}

Par_Task
par_delegate_task ()
{
  int i;
  int new_task_completed_size;
  int index, bit;
  Task* task;

  if (par_verbose)
    Report("par_delegate_task called.\n");

  /* if we are not running in parallel, do the worker task ourself */
  if (!par)
    {
      if (par_verbose)
	Report("Performing task %d myself\n", task_number);
      (*par_worker_task)();
      if (par_master_result != NULL)
	(*par_master_result)(task_number);
      return(task_number++);
    }

  if (task_number - task_completed_first >= 8*task_completed_size)
    {
      /* double the size of the task_completed table */
      if (task_completed_size == 0)
	{
	  new_task_completed_size = 1024;
	  task_completed = (unsigned char *) malloc(new_task_completed_size *
						    sizeof(unsigned char));
        }
      else
	{
	  new_task_completed_size = 2*task_completed_size;
	  task_completed = realloc(task_completed,
				   new_task_completed_size *
				   sizeof(unsigned char));
	}
      for (i = task_completed_offset; i < task_completed_size; ++i)
	task_completed[task_completed_size + i] = task_completed[i];
      task_completed_offset += task_completed_size;
      task_completed_size = new_task_completed_size;
    }
  index = (((task_number - task_completed_first) >> 3) + task_completed_offset)
    % task_completed_size;
  bit = (task_number - task_completed_first) & 7;
  task_completed[index] &= ~(0xff << bit);

  ++tasks_outstanding;
  if (context_changed)
    {
      current_context = (Context *) malloc(sizeof(Context));
      if (par_verbose)
	Report("Context allocated at %x\n", current_context);
      current_context->use_count = 1;
      out_position = 0;
      if (par_pack_context != NULL)
	(*par_pack_context)();
      current_context->buffer.buffer = out_buffer;
      current_context->buffer.size = out_size;
      current_context->buffer.position = out_position;
      out_buffer = NULL;
      out_size = 0;
      out_position = 0;
      context_changed = FALSE;
    }
  task = (Task*) malloc(sizeof(Task));
  task->next = NULL;
  task->prev = NULL;
  task->context = ReuseContext(current_context);
  task->number = task_number;

  out_position = 0;
  par_pkint(task_number);
  if (par_pack_task != NULL)
    (*par_pack_task)();
  task->buffer.buffer = out_buffer;
  task->buffer.size = out_size;
  task->buffer.position = out_position;
  out_buffer = NULL;
  out_size = 0;
  out_position = 0;

  if (par_verbose)
    Report("par_delegate_task queueing task.\n");

  QueueTask(task);

  if (par_verbose)
    Report("par_delegate_task dispatching tasks.\n");

  DispatchTasks();

  if (par_verbose)
    Report("par_delegate_task returning %d.\n", task_number+1);
  return(task_number++);
}

void
par_set_context ()
{
  /* if we are not running in parallel, call WorkerContext
     directly */
  if (!par)
    {
      if (par_worker_context != NULL)
	(*par_worker_context)();
      return;
    }

  if (current_context != NULL)
    DisuseContext(&current_context);
  context_changed = TRUE;
}

void
par_broadcast_context ()
{
  int total_time;
  int bufid;
  int i;
  int oldid;
  struct timeval current_time;

  /* if we are not running in parallel, call WorkerContext directly */
  if (!par)
    {
      if (par_worker_context != NULL)
	(*par_worker_context)();
      return;
    }

  if (par_verbose)
    Report("Master going to broadcast context %d\n", broadcast_count);

  /* wait for all tasks to complete */
  while (tasks_outstanding > 0)
    (void) MasterReceiveMessage(PAR_FOREVER);

  /* wait for any pending workers to start up */
  total_time = 0;
  while ((n_hosts_pending > 0 || n_workers_pending > 0) &&
	 total_time < PENDING_WORKERS_TIMEOUT)
    {
      if (gettimeofday(&current_time, NULL) != 0)
	Abort("Master could not get time of day.\n");
      //      printf("pbc: %d %d %d %d %d\n",
      //	     n_hosts_pending, n_workers_pending,
      //	     total_time,
      //	     current_time.tv_sec, current_time.tv_usec);
      if (MasterReceiveMessage(REPORT_INTERVAL) == 0)
	total_time += REPORT_INTERVAL;
    }
  if (gettimeofday(&current_time, NULL) != 0)
    Abort("Master could not get time of day.\n");
  //  printf("pbc_cls: %d %d %d %d %d\n",
  //	 n_hosts_pending, n_workers_pending,
  //	 total_time,
  //	 current_time.tv_sec, current_time.tv_usec);
  fflush(stdout);
  /* give up on any that have not responded */
  n_hosts_pending = 0;
  n_workers_pending = 0;

  /* if there are many broadcasts already in progress, we have to
     wait here for acknowledgements */
  while (broadcast_count > MIN(broadcast_first_ack_count,
			       broadcast_second_ack_count) +
	                   PAR_MAX_OVERLAPPED_BROADCASTS)
    MasterReceiveMessage(PAR_FOREVER);

  if (broadcast_count == 0)
    {
      if (n_workers <= 0)
	Abort("Broadcast operation not possible without any workers!\n");
      process_tids[0] = master_tid;
      for (i = 0; i < n_workers; ++i)
	process_tids[i+1] = workers[i].tid;
      n_processes = n_workers + 1;
      if (n_processes > 1)
	broadcast_first_forward_tid = process_tids[1];
      if (n_processes > 2)
	broadcast_second_forward_tid = process_tids[2];
    }
  if (par_verbose)
    Report("Master broadcasting context %d\n", broadcast_count);
  if (broadcast_first_forward_tid >= 0)
    {
      /* prepare message */
      PrepareToSend();
      par_pkint(broadcast_count);
      if (broadcast_count == 0)
	{
	  par_pkint(n_processes);
	  par_pkintarray(process_tids, n_processes);
	}
      if (par_pack_context != NULL)
	(*par_pack_context)();
      Send(broadcast_first_forward_tid, BROADCAST_CONTEXT_MSG);
    }
  else 
    broadcast_first_ack_count = broadcast_count;

  if (broadcast_second_forward_tid >= 0)
    {
      /* prepare message */
      PrepareToSend();
      par_pkint(broadcast_count);
      if (broadcast_count == 0)
	{
	  par_pkint(n_processes);
	  par_pkintarray(process_tids, n_processes);
	}
      if (par_pack_context != NULL)
	(*par_pack_context)();
      Send(broadcast_second_forward_tid, BROADCAST_CONTEXT_MSG);
    }
  else 
    broadcast_second_ack_count = broadcast_count;
  ++broadcast_count;
}

int
par_wait (float timeout)
{
  /* if we are not running in parallel, there is nothing to do */
  if (!par)
    return(0);
  return(MasterReceiveMessage(timeout));
}

void
par_finish ()
{
  int i;
  int total_time;

  /* if we are not running in parallel, there is nothing to do */
  if (!par)
    return;

  while (tasks_outstanding > 0)
    (void) MasterReceiveMessage(PAR_FOREVER);

  /* wait for all workers to start up before shutting down */
  while (n_workers_pending > 0)
    MasterReceiveMessage(PAR_FOREVER);

  if (par_verbose)
    Report("Sending terminate messages to %d workers\n",n_workers);
  for (i = 0; i < n_workers; ++i)
    if (MPI_Send(out_buffer, 0, MPI_PACKED,
		   workers[i].tid, TERMINATE_MSG,
		   MPI_COMM_WORLD) != MPI_SUCCESS)
      Abort("Master cannot send TERMINATE message to worker.\n");

  n_workers = 0;   /* consider all workers dead */
}

int
par_finished (Par_Task id)
{
  int index, bit;

  if (!par)
    return(id < task_number);

  if (id < task_completed_first)
    return(TRUE);
  if (id >= task_number)
    return(FALSE);
  index = (((id - task_completed_first) >> 3) + task_completed_offset) %
    task_completed_size;
  bit = (id - task_completed_first) & 7;
  return((task_completed[index] >> bit) & 1);
}

int
par_finished_up_to (Par_Task id)
{
  int index, bit, test_bits;

  if (!par)
    return(id < task_number);

  if (id < task_completed_first)
    return(TRUE);
  if (id >= (task_completed_first+8) ||
      id >= task_number)
    return(FALSE);
  index = (((id - task_completed_first) >> 3) + task_completed_offset) %
    task_completed_size;
  bit = (id - task_completed_first) & 7;
  test_bits = 0xff >> (7-bit);
  return((task_completed[index] & test_bits) == test_bits);
}

int
par_tasks_outstanding ()
{
  if (!par)
    return(0);
  return(tasks_outstanding);
}

Par_Task par_next_task()
{
  return(task_number);
}

int
par_enabled ()
{
  return(par);
}

int
par_workers ()
{
  if (!par)
    return(1);
  return(n_workers);
}

int
par_instance ()
{
  if (!par)
    return(0);
  return(rank);
}

const char* 
par_arch ()
{
  static char* arch= NULL;
  FILE *f = popen("uname -s", "r");
  if (f != NULL)
    {
      char s[1024];
      if (fscanf(f, "%s", s) == 1)
	arch = strdup(s);
      else
	arch = "DUMMY";
      pclose(f);
    }
  else
    arch = "DUMMY";

  return arch;
}

/* DispatchTasks does not return until all queued tasks have been assigned to
   workers for processing */
static void
DispatchTasks ()
{
  /* make sure all broadcasts have completed */
  while (broadcast_first_ack_count < (broadcast_count - 1) ||
	 broadcast_second_ack_count < (broadcast_count - 1))
    MasterReceiveMessage(PAR_FOREVER);

  if (par_verbose)
    Report("DispatchTasks: all broadcasts completed\n");

  /* dispatch all tasks */
  while (first_queued_task != NULL)
      DispatchTask();
}

static void
DispatchTask ()
{
  int n;
  int res;
  Task *task;

  /* take the first task off the queue */
  task = first_queued_task;
  first_queued_task = task->next;
  if (first_queued_task != NULL)
    first_queued_task->prev = NULL;
  else
    last_queued_task = NULL;
  task->next = NULL;

  /* find a worker to run it on */
  n = FindReadyWorker();

  /* check if the context needs to be sent */
  if (task->context != NULL &&
      task->context != workers[n].last_context)
    {
      if (MPI_Send(task->context->buffer.buffer, task->context->buffer.position,
		   MPI_PACKED, workers[n].tid, CONTEXT_MSG,
		   MPI_COMM_WORLD) != MPI_SUCCESS)
	Abort("Error sending context to worker.\n");

      DisuseContext(&workers[n].last_context);
      workers[n].last_context = ReuseContext(task->context);
    }

  /* send the task */
  if (par_verbose)
    Report("Master sending task %d to worker %d on %s\n",
	   task->number, n, workers[n].host);

  if (MPI_Send(task->buffer.buffer, task->buffer.position, MPI_PACKED,
	       workers[n].tid, TASK_MSG, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Error sending task to worker.\n");

  /* record that worker n is performing the task */
  if (workers[n].first_task == NULL)
    {
      workers[n].first_task = task;
      RemoveFromIdleList(n);
    }
  else Abort("Worker %d requested more than 1 task at a time.\n", n);
}

static int
FindReadyWorker ()
{
  for (;;)
    {
      if (par_verbose)
	Report("FindReadyWorker: idle_workers = %d\n", idle_workers);

      /* if a worker is totally idle, choose it */
      if (idle_workers >= 0)
	return(idle_workers);

      /* wait for something to happen */
      (void) MasterReceiveMessage(PAR_FOREVER);
    }
}

static int
MasterReceiveMessage (float timeout)
{
  int flag;
  int len;
  MPI_Status status;
  struct timeval end_time, current_time;
  struct timespec duration;

  if (par_verbose)
    Report("Master entering MasterReceiveMessage.\n");
  if (timeout == 0.0)
    {
      /* check if any message has arrived */
      if (MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
		     &flag, &status) != MPI_SUCCESS)
	Abort("Could not probe for messages in MasterReceiveMessage()\n");

      if (!flag)
	/* no message has arrived */
	return(0);

      /* a message is there; now receive it */
      if (MPI_Get_count(&status, MPI_PACKED, &len) != MPI_SUCCESS)
	Abort("Could not obtain number of bytes in message\n");

      if (par_verbose)
	Report("Master going to receive message of %d bytes\n", len);

      if (len > in_size)
	ExpandInBuffer(len);

      if (MPI_Recv(in_buffer, len, MPI_PACKED, MPI_ANY_SOURCE, MPI_ANY_TAG,
		   MPI_COMM_WORLD, &status) != MPI_SUCCESS)
	Abort("Master could not receive message.\n");
    }
  else if (timeout < 0.0)
    {
      /* we can't spawn new workers, so we'll have
	 to be patient and wait for one to free up */
      if (MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
		    &status) != MPI_SUCCESS)
	Abort("Could not probe for messages in MasterReceiveMessage()\n");

      /* a message is there; now receive it */
      if (MPI_Get_count(&status, MPI_PACKED, &len) != MPI_SUCCESS)
	Abort("Could not obtain number of bytes in message\n");

      if (par_verbose)
	Report("Master going to receive message of %d bytes\n", len);

      if (len > in_size)
	ExpandInBuffer(len);

      if (MPI_Recv(in_buffer, len, MPI_PACKED, MPI_ANY_SOURCE, MPI_ANY_TAG,
		   MPI_COMM_WORLD, &status) != MPI_SUCCESS)
	Abort("Master could not receive message.\n");
    }
  else
    {
      /* wait until a message arrives or the timeout has expired */
      /* since MPI does not support timeouts, and we don't want to use
	 multiple threads, we have to poll (yuck) */
      end_time.tv_sec = 0;
      if (gettimeofday(&end_time, NULL) != 0)
	Abort("Master could not get time of day.\n");
      end_time.tv_sec += floor(timeout);
      end_time.tv_usec += floor(1000000.0*(timeout - floor(timeout)));
      if (end_time.tv_usec >= 1000000)
	{
	  end_time.tv_sec++;
	  end_time.tv_usec -= 1000000;
	}
      duration.tv_sec = 0;
      duration.tv_nsec = 10000000;	/* 10 milliseconds */
      do {
	/* check if any message has arrived */
	if (MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
		       &flag, &status) != MPI_SUCCESS)
	  Abort("Could not probe for messages in MasterReceiveMessage()\n");

	if (flag)
	  break;

	nanosleep(&duration, NULL);

	if (gettimeofday(&current_time, NULL) != 0)
	  Abort("Master could not get time of day.\n");
      } while (current_time.tv_sec < end_time.tv_sec ||
	       (current_time.tv_sec == end_time.tv_sec &&
		current_time.tv_usec < end_time.tv_usec));

      if (!flag)
	return(0);

      /* a message is there; now receive it */
      if (MPI_Get_count(&status, MPI_PACKED, &len) != MPI_SUCCESS)
	Abort("Could not obtain number of bytes in message\n");
      if (len > in_size)
	ExpandInBuffer(len);
      
      if (MPI_Recv(in_buffer, len, MPI_PACKED, MPI_ANY_SOURCE, MPI_ANY_TAG,
		   MPI_COMM_WORLD, &status) != MPI_SUCCESS)
	Abort("Master could not receive message.\n");
    }

  in_position = 0;
  HandleMessage(status.MPI_TAG, status.MPI_SOURCE);
  return(1);
}

static void
HandleMessage (int msg_tag, int tid)
{
  switch (msg_tag)
    {
    case REQUEST_MSG:
      HandleRequest(tid);
      break;
	      
    case BROADCAST_ACK_MSG:
      HandleBroadcastAck(tid);
      break;

    case WORKER_EXIT_MSG:
      HandleWorkerExit();
      break;

    default:
      Error("Master received unknown message: %d\n", msg_tag);
    }
}

static void
HandleRequest (int tid)
{
  int tc;
  int n;
  int result_appended;
  int index, bit;
  int current;
  Hostname worker_host_name;

  /* locate the worker in the worker table */
  tc = par_upkint();
  if (tc < 0)
    par_upkstr(worker_host_name);

  if (par_verbose)
    Report("Master: HandleRequest called (tid = %d tc = %d %s\n",
	   tid, tc, (tc < 0) ? worker_host_name : ""); 

  for (n = 0; n < n_workers; ++n)
    if (workers[n].tid == tid)
      break;	/* found it */

  if (n >= n_workers)
    {
      if (par_verbose)
        Report("Received msg from new worker: %d\n", tid);
      --n_workers_pending;
      /* check that integer value sent is -1 */
      if (tc != -1)
        Abort("Initial worker message was invalid.\n");
      /* add it to the table of workers */
      workers[n_workers].tid = tid;
      strcpy(workers[n_workers].host, worker_host_name);
      workers[n_workers].idle = FALSE;
      workers[n_workers].last_context = NULL;
      workers[n_workers].first_task = NULL;
      PutOnIdleList(n_workers);
      if (par_verbose)
        Report("Worker %d started on host %s\n", n_workers, workers[n_workers].host);
      ++n_workers;
      return;
    }

  /* this is an old worker */
  /* update the first_task field */
  if (tc >= 0 && workers[n].first_task->number == tc)
    {
      /* we are done with the first task */
      if (par_verbose)
	Report("Master received result of task %d from worker %d on host %s\n",
	       tc, n, workers[n].host);
      result_appended = par_upkint();
      if (result_appended && par_master_result != NULL)
	{
	  if (par_unpack_result != NULL)
	    (*par_unpack_result)();
	  (*par_master_result)(tc);
	}
      --tasks_outstanding;
      DisuseContext(&workers[n].first_task->context);

      /* free up the task buffer */
      FreeBuffer(&workers[n].first_task->buffer);

      /* mark as finished */
      index = (((workers[n].first_task->number - task_completed_first) >> 3) +
	       task_completed_offset) % task_completed_size;
      bit = (workers[n].first_task->number - task_completed_first) & 7;
      task_completed[index] |= 1 << bit;
      if (index == task_completed_offset)
	{
	  current = (((task_number - task_completed_first) >> 3) +
		     task_completed_offset) % task_completed_size;
	  while (task_completed_offset != current)
	    {
	      if (task_completed[task_completed_offset] != 0xff)
		break;
	      task_completed_first += 8;
	      if (++task_completed_offset >= task_completed_size)
		task_completed_offset = 0;
	    }
	}

      free(workers[n].first_task);
      workers[n].first_task = NULL;
    }
  else
    if (tc != -1)
      Error("Worker %d on host %s completed task %d instead of %d\n",
	    n, workers[n].host, tc, workers[n].first_task->number);

  if (workers[n].first_task == NULL)
    PutOnIdleList(n);
}

static void
HandleBroadcastAck (int from_tid)
{
  int count;

  count = par_upkint();
  if (par_verbose)
    Report("Master received BROADCAST_ACK to %d from %d\n",
	   count, from_tid);
  if (from_tid == broadcast_first_forward_tid)
    broadcast_first_ack_count = count;
  else if (from_tid == broadcast_second_forward_tid)
    broadcast_second_ack_count = count;
  else
    Abort("Received BROADCAST_ACK_MESSAGE from unexpected source.\n");
}


static void
HandleWorkerExit ()
{
  int stid;
  int n;

  stid = par_upkint();
  for (n = 0; n < n_workers; ++n)
    if (workers[n].tid == stid)
      {
	Error("Worker %d on host %s exited.\n", n, workers[n].host);
	DeclareWorkerDead(n);
	return;
      }
  /* if we get here, we didn't have a record of the worker anyway,
     so ignore the message */
  if (par_verbose)
    Report ("Unknown Worker tid %d.\n", stid);
}

static void
QueueTask (Task* task)
{
  /* put the task at the end of the queued task list so that
     it will wait its turn before getting dispatched */
  task->next = NULL;
  task->prev = last_queued_task;
  if (last_queued_task != NULL)
    last_queued_task->next = task;
  else
    first_queued_task = task;
  last_queued_task = task;
}
     
static void
RequeueTask (Task *task)
{
  /* put the task at the front of the queued task list so that
     it will get dispatched as soon as possible */
  task->next = first_queued_task;
  task->prev = NULL;
  if (first_queued_task != NULL)
    first_queued_task->prev = task;
  else
    last_queued_task = task;
  first_queued_task = task;
}
     
static Context*
ReuseContext (Context *context)
{
  if (context != NULL)
    context->use_count++;
  return(context);
}

static void
DisuseContext (Context **pContext)
{
  Context *context;
  int i;
  Task* t;

  context = *pContext;
  *pContext = NULL;
  if (context == NULL || --context->use_count > 0)
    return;
  if (context->use_count < 0)
    fprintf(stderr, "CONTEXT_ERROR: use_count of context %p is %d\n",
	    context, context->use_count);
  FreeBuffer(&context->buffer);
  free(context);
  if (par_verbose)
    Report("Context freed at %x\n", context);

  /* make sure that context is not used anywhere */
  if (current_context == context)
    fprintf(stderr,
	    "CONTEXT ERROR: current_context is a pointer to deleted context %p\n",
	    context);
  for (i = 0; i < n_workers; ++i)
    {
      if (workers[i].last_context == context)
	fprintf(stderr,
		"CONTEXT ERROR: worker %d has a pointer to deleted context %p\n",
		i, context);
      if (workers[i].first_task != NULL &&
	  workers[i].first_task->context == context)
	fprintf(stderr,
		"CONTEXT ERROR: the first_task of worker %d has a pointer to deleted context %p\n",
		i, context);
    }
  for (t = first_queued_task; t != NULL; t = t->next)
    if (t->context == context)
      fprintf(stderr,
	      "CONTEXT ERROR: task %p on task queue has a pointer to deleted context %p\n",
	      t, context);
}

static void
FreeBuffer (Buffer *buffer)
{
  if (buffer->buffer != NULL)
    {
      free(buffer->buffer);
      buffer->buffer = NULL;
    }
}

static void
PutOnIdleList (int n)
{
  if (workers[n].idle)
    /* worker was already on idle list */
    return;

  /* put on the idle list */
  workers[n].idle = TRUE;
  workers[n].next_idle = idle_workers;
  idle_workers = n;
}

static void
RemoveFromIdleList (int n)
{
  int pi;
  int i;

  if (!workers[n].idle)
    /* worker wasn't on idle list */
    return;

  /* take this worker off the idle list */
  workers[n].idle = FALSE;
  pi = -1;
  i = idle_workers;
  while (i >= 0)
    {
      if (i == n)
	{
	  if (pi >= 0)
	    workers[pi].next_idle = workers[n].next_idle;
	  else
	    idle_workers = workers[n].next_idle;
	  return;
	}
      pi = i;
      i = workers[i].next_idle;
    }
  Abort("Corrupt idle list.\n");
}

static void
DeclareWorkerDead (int n)
{
  if (workers[n].first_task != NULL)
    RequeueTask(workers[n].first_task);
  RemoveFromIdleList(n);

  --n_workers;
  if (n != n_workers)
    {
      /* move the last worker into the slot formerly occupied by the
	 dead worker */
      workers[n] = workers[n_workers];
      if (workers[n_workers].idle)
	{
	  RemoveFromIdleList(n_workers);
	  workers[n].idle = FALSE;
	  PutOnIdleList(n);
	}
    }
}

static void
PerformWorkerTasks ()
{
  int bufid;
  int msg_type;
  int task_num;
  int from_tid;
  int i;
  int dest;
  int oldid;
  struct timeval task_start;
  struct timeval task_end;
  long task_sec;
  long task_usec;

#ifdef PTHREADS
  /* launch another thread that will actually execute the tasks */
#endif

  ComposeRequest(-1, FALSE);
  Send(master_tid, REQUEST_MSG);

  for (;;)
    {
      msg_type = WorkerReceiveMessage(&from_tid);

      switch (msg_type)
	{
	case CONTEXT_MSG:
	  if (par_unpack_context != NULL)
	    (*par_unpack_context)();
	  if (par_worker_context != NULL)
	    (*par_worker_context)();
	  break;
	case TASK_MSG:
	  task_num = par_upkint();
	  if (par_verbose)
	    Report("Worker received task %d\n", task_num);
	  if (par_unpack_task != NULL)
	    (*par_unpack_task)();
	  if (par_worker_task != NULL) {
	    if (gettimeofday(&task_start, NULL) != 0)
	      Abort("Worker could not get time of day.\n");
	    (*par_worker_task)();
	    if (gettimeofday(&task_end, NULL) != 0)
	      Abort("Worker could not get time of day.\n");
	    task_sec= task_end.tv_sec-task_start.tv_sec;
	    task_usec= task_end.tv_usec-task_start.tv_usec;
	    if (task_usec<0) 
	      {
		task_sec -= 1;
		task_usec += 1000000;
	      }
	    if (task_sec > longest_task.tv_sec 
		||(task_sec==longest_task.tv_sec 
		   && task_usec>longest_task.tv_usec))
	      {
		longest_task.tv_sec= task_sec;
		longest_task.tv_usec= task_usec;
	      }
	  }
	  if (par_master_result != NULL)
	    {
	      ComposeRequest(task_num, TRUE);
	      if (par_pack_result != NULL)
		(*par_pack_result)();
	      Send(master_tid, REQUEST_MSG);
	    }
	  else
	    {
	      ComposeRequest(task_num, FALSE);
	      Send(master_tid, REQUEST_MSG);
	    }
	  break;
	case BROADCAST_CONTEXT_MSG:
	  broadcast_count = par_upkint();
	  if (broadcast_count == 0)
	    {
	      n_processes = par_upkint();
	      par_upkintarray(process_tids, n_processes);
	      /* find my_process_index */
	      for (i = 0; i < n_processes; ++i)
		if (process_tids[i] == my_tid)
		  {
		    my_process_index = i;
		    break;
		  }
	      if (i >= n_processes)
		Abort("Could not find my process tid in process array.\n");
	      dest = (my_process_index << 1) + 1;
	      if (dest < n_processes)
		broadcast_first_forward_tid = process_tids[dest];
	      dest = (my_process_index << 1) + 2;
	      if (dest < n_processes)
		broadcast_second_forward_tid = process_tids[dest];
	      dest = (my_process_index - 1) >> 1;
	      broadcast_backward_tid = process_tids[dest];
	    }
	  if (par_verbose)
	    Report("Process %d received broadcast context %d\n",
		   my_process_index, broadcast_count);
	  /* handle the broadcast message locally */
	  if (par_unpack_context != NULL)
	    (*par_unpack_context)();
	  if (par_worker_context != NULL)
	    (*par_worker_context)();

	  if (broadcast_first_forward_tid >= 0)
	    {
	      if (par_verbose)
		Report("Process %d relaying context %d onto tid %d\n",
		       my_process_index, broadcast_count,
		       broadcast_first_forward_tid);

	      /* prepare message */
	      PrepareToSend();
	      par_pkint(broadcast_count);
	      if (broadcast_count == 0)
		{
		  par_pkint(n_processes);
		  par_pkintarray(process_tids, n_processes);
		}
	      if (par_pack_context != NULL)
		(*par_pack_context)();
	      Send(broadcast_first_forward_tid, BROADCAST_CONTEXT_MSG);
	    }
	  else
	    broadcast_first_ack_count = broadcast_count;

	  if (broadcast_second_forward_tid >= 0)
	    {
	      if (par_verbose)
		Report("Process %d relaying context %d onto tid %d\n",
		       my_process_index, broadcast_count,
		       broadcast_second_forward_tid);

	      /* prepare message */
	      PrepareToSend();
	      par_pkint(broadcast_count);
	      if (broadcast_count == 0)
		{
		  par_pkint(n_processes);
		  par_pkintarray(process_tids, n_processes);
		}
	      if (par_pack_context != NULL)
		(*par_pack_context)();
	      Send(broadcast_second_forward_tid, BROADCAST_CONTEXT_MSG);
	    }
	  else 
	    broadcast_second_ack_count = broadcast_count;

	  if (broadcast_first_ack_count == broadcast_count &&
	      broadcast_second_ack_count == broadcast_count)
	    {
	      /* we can acknowledge this immediately */
	      if (par_verbose)
		Report("Process %d sending immediate acknowledge of %d\n",
		       my_process_index, broadcast_count);
	      PrepareToSend();
	      par_pkint(broadcast_count);
	      Send(broadcast_backward_tid, BROADCAST_ACK_MSG);
	    }
	  break;
	case BROADCAST_ACK_MSG:
	  broadcast_count = par_upkint();
	  if (par_verbose)
	    Report("Process %d received acknowledgement of %d from %d\n",
		   my_process_index, broadcast_count, from_tid);
	  if (from_tid == broadcast_first_forward_tid)
	    {
	      broadcast_first_ack_count = broadcast_count;
	      if (broadcast_first_ack_count <= broadcast_second_ack_count)
		{
		  /* relay acknowledgement */
		  if (par_verbose)
		    Report("Process %d relaying acknowledgement of %d\n",
			   my_process_index, broadcast_count);
		  PrepareToSend();
		  par_pkint(broadcast_count);
		  Send(broadcast_backward_tid, BROADCAST_ACK_MSG);
		}
	    }
	  else if (from_tid == broadcast_second_forward_tid)
	    {
	      broadcast_second_ack_count = broadcast_count;
	      if (broadcast_second_ack_count <= broadcast_first_ack_count)
		{
		  /* relay acknowledgement */
		  if (par_verbose)
		    Report("Process %d relaying acknowledgement of %d\n",
			   my_process_index, broadcast_count);
		  PrepareToSend();
		  par_pkint(broadcast_count);
		  Send(broadcast_backward_tid, BROADCAST_ACK_MSG);
		}
	    }
	  else
	    Abort("Received BROADCAST_ACK_MSG from unexpected source.\n");
	  break;
	case TERMINATE_MSG:
	  if (par_verbose)
	    Report("Worker received TERMINATE_MSG: tid=%d\n", my_tid);
	  return;
	default:
	  Abort("Unknown message type received.\n");
	  break;
	}
    }
}

static void
ComposeRequest (int last_task_completed, int result_appended)
{
  if (par_verbose)
    Report("Worker requesting task\n");
  PrepareToSend();
  par_pkint(last_task_completed);
  if (last_task_completed < 0)
    par_pkstr(my_host_name);
  par_pkint(result_appended);
}

static int
WorkerReceiveMessage (int *pfrom_tid)
{
  int flag;
  int len;
  MPI_Status status;
  struct timeval end_time, current_time;
  struct timespec duration;

  /* wait until a message arrives or the timeout has expired */
  /* since MPI does not support timeouts, and we don't want to use
     multiple threads, we have to poll (yuck) */
  end_time.tv_sec = 0;
  if (gettimeofday(&end_time, NULL) != 0)
    Abort("Worker could not get time of day.\n");
  end_time.tv_sec += MAX(WORKER_RECEIVE_TIMEOUT,10*(longest_task.tv_sec+1));
  duration.tv_sec = 0;
  duration.tv_nsec = 10000000;	/* 10 milliseconds */
  do {
    /* check if any message has arrived */
    if (MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
		   &flag, &status) != MPI_SUCCESS)
      Abort("Could not probe for messages in WorkerReceiveMessage()\n");

    if (flag)
      break;

    nanosleep(&duration, NULL);

    if (gettimeofday(&current_time, NULL) != 0)
      Abort("Worker could not get time of day.\n");
  } while (current_time.tv_sec < end_time.tv_sec ||
	   (current_time.tv_sec == end_time.tv_sec &&
	    current_time.tv_usec < end_time.tv_usec));

  if (!flag)
    {
      Report("Worker timed out waiting for message. Exiting...\n");
      PrepareToSend ();
      par_pkint (my_tid);
      Send (master_tid, WORKER_EXIT_MSG);
      /* Fake msg so finalize will be done before exit.*/
      return (TERMINATE_MSG); 
    }

  /* a message is there; now receive it */
  if (MPI_Get_count(&status, MPI_PACKED, &len) != MPI_SUCCESS)
    Abort("Could not obtain number of bytes in message\n");
  if (len > in_size)
    ExpandInBuffer(len);

  if (MPI_Recv(in_buffer, len, MPI_PACKED, MPI_ANY_SOURCE, MPI_ANY_TAG,
	       MPI_COMM_WORLD, &status) != MPI_SUCCESS)
    Abort("Worker could not receive message.\n");

  in_position = 0;

  *pfrom_tid = status.MPI_SOURCE;
  return(status.MPI_TAG);
}

static void
PrepareToSend ()
{
  out_position = 0;
}

static void
Send (int tid, int tag)
{
  if (par_verbose)
    Report("Worker %d sending message %d to %d.\n",
	   rank, tag, tid);
  if (MPI_Send(out_buffer, out_position, MPI_PACKED,
	       tid, tag, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Cannot send message (tid = %d tag = %d)\n", tid, tag);
}

void
par_pkbyte(unsigned char v)
{
  if (out_position + sizeof_byte > out_size)
    ExpandOutBuffer(out_position + sizeof_byte);
  if (MPI_Pack(&v, 1, MPI_BYTE, out_buffer, out_size, &out_position,
	       MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not pack byte into MPI buffer\n");
}

void
par_pkshort(short v)
{
  if (out_position + sizeof_short > out_size)
    ExpandOutBuffer(out_position + sizeof_short);
  if (MPI_Pack(&v, 1, MPI_SHORT, out_buffer, out_size, &out_position,
	       MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not pack short into MPI buffer\n");
}

void
par_pkint(int v)
{
  if (out_position + sizeof_int > out_size)
    ExpandOutBuffer(out_position + sizeof_int);
  if (MPI_Pack(&v, 1, MPI_INT, out_buffer, out_size, &out_position,
	       MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not pack int into MPI buffer\n");
}

void
par_pklong (long v)
{
  if (out_position + sizeof_long > out_size)
    ExpandOutBuffer(out_position + sizeof_long);
  if (MPI_Pack(&v, 1, MPI_LONG, out_buffer, out_size, &out_position,
	       MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not pack long into MPI buffer\n");
}

void
par_pkfloat (float v)
{
  if (out_position + sizeof_float > out_size)
    ExpandOutBuffer(out_position + sizeof_float);
  if (MPI_Pack(&v, 1, MPI_FLOAT, out_buffer, out_size, &out_position,
	       MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not pack float into MPI buffer\n");
}

void
par_pkdouble (double v)
{
  if (out_position + sizeof_double > out_size)
    ExpandOutBuffer(out_position + sizeof_double);
  if (MPI_Pack(&v, 1, MPI_DOUBLE, out_buffer, out_size, &out_position,
	       MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not pack double into MPI buffer\n");
}

void
par_pkstr (char *v)
{
  int size;
  int string_len = strlen(v);

  if (MPI_Pack_size(string_len, MPI_CHAR, MPI_COMM_WORLD,
		    &size) != MPI_SUCCESS)
    Abort("Could not determine packing size of string\n");
  if (out_position + sizeof_int + size > out_size)
    ExpandOutBuffer(out_position + sizeof_int + size);
  if (MPI_Pack(&string_len, 1, MPI_INT, out_buffer, out_size, &out_position,
	       MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Pack(v, string_len, MPI_CHAR, out_buffer, out_size, &out_position,
	       MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not pack string into MPI buffer\n");
}

void
par_pkbytearray (unsigned char *p, int n)
{
  int size;
  if (MPI_Pack_size(n, MPI_BYTE, MPI_COMM_WORLD,
		    &size) != MPI_SUCCESS)
    Abort("Could not determine packing size of byte array\n");
  if (out_position + size > out_size)
    ExpandOutBuffer(out_position + size);
  if (MPI_Pack(p, n, MPI_BYTE, out_buffer, out_size, &out_position,
	       MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not pack byte array into MPI buffer\n");
}

void
par_pkshortarray (short *p, int n)
{
  int size;
  if (MPI_Pack_size(n, MPI_SHORT, MPI_COMM_WORLD,
		    &size) != MPI_SUCCESS)
    Abort("Could not determine packing size of short array\n");
  if (out_position + size > out_size)
    ExpandOutBuffer(out_position + size);
  if (MPI_Pack(p, n, MPI_SHORT, out_buffer, out_size, &out_position,
	       MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not pack short array into MPI buffer\n");
}

void
par_pkintarray (int *p, int n)
{
  int size;
  if (MPI_Pack_size(n, MPI_INT, MPI_COMM_WORLD,
		    &size) != MPI_SUCCESS)
    Abort("Could not determine packing size of int array\n");
  if (out_position + size > out_size)
    ExpandOutBuffer(out_position + size);
  if (MPI_Pack(p, n, MPI_INT, out_buffer, out_size, &out_position,
	       MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not pack int array into MPI buffer\n");
}

void
par_pklongarray (long *p, int n)
{
  int size;
  if (MPI_Pack_size(n, MPI_LONG, MPI_COMM_WORLD,
		    &size) != MPI_SUCCESS)
    Abort("Could not determine packing size of long array\n");
  if (out_position + size > out_size)
    ExpandOutBuffer(out_position + size);
  if (MPI_Pack(p, n, MPI_LONG, out_buffer, out_size, &out_position,
	       MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not pack long array into MPI buffer\n");
}

void
par_pkfloatarray (float *p, int n)
{
  int size;
  if (MPI_Pack_size(n, MPI_FLOAT, MPI_COMM_WORLD,
		    &size) != MPI_SUCCESS)
    Abort("Could not determine packing size of float array\n");
  if (out_position + size > out_size)
    ExpandOutBuffer(out_position + size);
  if (MPI_Pack(p, n, MPI_FLOAT, out_buffer, out_size, &out_position,
	       MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not pack float array into MPI buffer\n");
}

void
par_pkdoublearray (double *p, int n)
{
  int size;
  if (MPI_Pack_size(n, MPI_DOUBLE, MPI_COMM_WORLD,
		    &size) != MPI_SUCCESS)
    Abort("Could not determine packing size of double array\n");
  if (out_position + size > out_size)
    ExpandOutBuffer(out_position + size);
  if (MPI_Pack(p, n, MPI_DOUBLE, out_buffer, out_size, &out_position,
	       MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not pack double array into MPI buffer\n");
}

unsigned char
par_upkbyte ()
{
  unsigned char v;
  if (MPI_Unpack(in_buffer, in_size, &in_position,
		 &v, 1, MPI_BYTE, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not unpack byte from MPI buffer\n");
  return(v);
}

short
par_upkshort ()
{
  short v;
  if (MPI_Unpack(in_buffer, in_size, &in_position,
		 &v, 1, MPI_SHORT, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not unpack short from MPI buffer\n");
  return(v);
}

int
par_upkint ()
{
  int v;
  if (MPI_Unpack(in_buffer, in_size, &in_position,
		 &v, 1, MPI_INT, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not unpack int from MPI buffer\n");
  return(v);
}

long
par_upklong ()
{
  long v;
  if (MPI_Unpack(in_buffer, in_size, &in_position,
		 &v, 1, MPI_LONG, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not unpack long from MPI buffer\n");
  return(v);
}

float
par_upkfloat ()
{
  float v;
  if (MPI_Unpack(in_buffer, in_size, &in_position,
		 &v, 1, MPI_FLOAT, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not unpack float from MPI buffer\n");
  return(v);
}

double
par_upkdouble ()
{
  double v;
  if (MPI_Unpack(in_buffer, in_size, &in_position,
		 &v, 1, MPI_DOUBLE, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not unpack double from MPI buffer\n");
  return(v);
}

void
par_upkstr (char *s)
{
  int string_len;
  if (MPI_Unpack(in_buffer, in_size, &in_position,
		 &string_len, 1, MPI_INT, MPI_COMM_WORLD) != MPI_SUCCESS ||
      MPI_Unpack(in_buffer, in_size, &in_position,
		 s, string_len, MPI_CHAR, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not unpack string from MPI buffer\n");
  s[string_len] = '\0';
}

void
par_upkbytearray (unsigned char *p, int n)
{
  if (MPI_Unpack(in_buffer, in_size, &in_position,
		 p, n, MPI_BYTE, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not unpack byte array from MPI buffer\n");
}

void
par_upkshortarray (short *p, int n)
{
  if (MPI_Unpack(in_buffer, in_size, &in_position,
		 p, n, MPI_SHORT, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not unpack short array from MPI buffer\n");
}

void
par_upkintarray (int *p, int n)
{
  if (MPI_Unpack(in_buffer, in_size, &in_position,
		 p, n, MPI_INT, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not unpack int array from MPI buffer\n");
}

void
par_upklongarray (long *p, int n)
{
  if (MPI_Unpack(in_buffer, in_size, &in_position,
		 p, n, MPI_LONG, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not unpack long array from MPI buffer\n");
}

void
par_upkfloatarray (float *p, int n)
{
  if (MPI_Unpack(in_buffer, in_size, &in_position,
		 p, n, MPI_FLOAT, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not unpack float array from MPI buffer\n");
}

void
par_upkdoublearray (double *p, int n)
{
  if (MPI_Unpack(in_buffer, in_size, &in_position,
		 p, n, MPI_DOUBLE, MPI_COMM_WORLD) != MPI_SUCCESS)
    Abort("Could not unpack double array from MPI buffer\n");
}


static void
StartMPI (int argc,
	  char **argv,
	  char **envp)
{	
  if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS)
    Abort("Could not obtain rank\n");
  my_tid = rank;

  if (MPI_Pack_size(1, MPI_BYTE, MPI_COMM_WORLD,
		    &sizeof_byte) != MPI_SUCCESS ||
      MPI_Pack_size(1, MPI_SHORT, MPI_COMM_WORLD,
		    &sizeof_short) != MPI_SUCCESS ||
      MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD,
		    &sizeof_int) != MPI_SUCCESS ||
      MPI_Pack_size(1, MPI_LONG, MPI_COMM_WORLD,
		    &sizeof_long) != MPI_SUCCESS ||
      MPI_Pack_size(1, MPI_FLOAT, MPI_COMM_WORLD,
		    &sizeof_float) != MPI_SUCCESS ||
      MPI_Pack_size(1, MPI_DOUBLE, MPI_COMM_WORLD,
		    &sizeof_double) != MPI_SUCCESS)
    Abort("Could not determine packing size of elementary types\n");

  if (rank == 0)
    {
      /* I am the master */
      master_tid = 0;
      
      if (MPI_Comm_size(MPI_COMM_WORLD, &n_workers_pending) != MPI_SUCCESS)
	Abort("Cannot obtain number of workers\n");
      
      --n_workers_pending;	/* because the master is not a worker */

      if (n_workers_pending == 0)
	{
	  par = FALSE;
	  return;
	}

      if (par_verbose)
	Report("I am the master with %d workers.\n", n_workers_pending);

      (*par_master_task)(prog_argc, prog_argv, envp);

      /* make sure everything is finished in case par_master_task does not
	 call par_finish itself */
      par_finish();
    }
  else
    {
      /* I am a worker */
      master_tid = 0;

      PerformWorkerTasks();
      if (par_worker_finalize != NULL)
	(*par_worker_finalize)();

#ifdef ACCT
      PrintAcct(argv[0], my_tid);
#endif
    }
}

static void
ExpandInBuffer (int size)
{
  while (in_size < size)
    if (in_size == 0)
      in_size = 1024;
    else
      in_size *= 2;
  in_buffer = realloc(in_buffer, in_size);
  if (in_buffer == NULL)
    Abort("Could not expand in_buffer to %d bytes\n", in_size);
}

static void
ExpandOutBuffer (int size)
{
  while (out_size < size)
    if (out_size == 0)
      out_size = 1024;
    else
      out_size *= 2;
  out_buffer = realloc(out_buffer, out_size);
  if (out_buffer == NULL)
    Abort("Could not expand out_buffer to %d bytes\n", out_size);
}

static Boolean
StringCopy (char *to,
	    const char *from,
	    const int maxChars)
{
  if (strlen(from) < maxChars)
    {
      strcpy(to, from);
      return(TRUE);
    }
  else
    {
      strncpy(to, from, maxChars-1);
      to[maxChars-1] = '\0';
      return(FALSE);
    }
}

static void
Report (char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);
  vprintf(fmt, args);
  va_end(args);
  fflush(stdout);
}

static void
Error (char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
}

static void
Abort (char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
  abort();
}

