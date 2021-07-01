#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>

#include "par.h"

typedef struct Context {
  /* NOTE: any new fields added to this struct should
     be also added to PackContext and UnpackContext */ 
  int context_num;
} Context;

typedef struct Task {
  /* NOTE: any new fields added to this struct should
     be also added to PackTask and UnpackTask */ 
  int task_num;
} Task;

typedef struct Result {
  int result_num;
} Result;

/* GLOBAL VARIABLES FOR MASTER & WORKER */
Context c;
Task t;
Result r;

/* GLOBAL VARIABLES FOR MASTER */
int total;

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

int
main (int argc, char **argv, char **envp)
{
  printf("Going to call par_process\n");
  par_process(argc, argv, envp,
              (void (*)()) MasterTask, MasterResult,
              WorkerContext, WorkerTask, NULL,
              PackContext, UnpackContext,
              PackTask, UnpackTask,
              PackResult, UnpackResult);
  printf("Going to exit\n");
  return(0);
}

/* MASTER PROCEDURES */

void
MasterTask (int argc,
            char **argv,
            char **envp)
{
  c.context_num = 137;
  printf("Before par_set_context\n");
  par_set_context();
  printf("After par_set_context\n");

  /* for all slices */
  total = 0;
  for (t.task_num = 0; t.task_num < 64; ++t.task_num)
    par_delegate_task();
  printf("Before par_finish\n");
  par_finish();
  printf("After par_finish\n");

  printf("Total = %d\n", total);
}

void
MasterResult ()
{
  printf("MasterResult called\n");
  total += r.result_num;
}


/* WORKER PROCEDURES */

void
WorkerContext ()
{
  printf("WorkerContext called\n");
  /* nothing necessary here */
}

void
WorkerTask ()
{
  printf("WorkerTask called\n");
  r.result_num = t.task_num;
}

void
PackContext ()
{
  par_pkint(c.context_num);
}

void
UnpackContext ()
{
  c.context_num = par_upkint();
}

void
PackTask ()
{
  par_pkint(t.task_num);
}

void
UnpackTask ()
{
  t.task_num = par_upkint();
}

void
PackResult ()
{
  par_pkint(r.result_num);
}

void
UnpackResult ()
{
  r.result_num = par_upkint();
}
