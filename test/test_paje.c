#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>
#include <math.h>

#include "TP_solver.h" 
#include "TP_matrix.h" 
#include "TP_dense.h"
#include "TP_enum.h"
#include "TP_pivot_list.h" 
#include "TP_auxiliary.h"


int
main(int argc, char **argv) 
{
  TP_solver self;
  
  self = TP_solver_create();
  TP_solver_parse_args(self, argc, argv);
  TP_solver_read_matrix(self);
  TP_solver_init(self);
  printf("%d !\n",omp_get_thread_num());

#pragma omp parallel num_threads(self->exe_parms->nb_threads)
  {
    TP_verbose_trace_start_event(self->verbose, 0);
    sleep(1);
    TP_verbose_trace_stop_event(self->verbose);
    sleep(2);
    TP_verbose_trace_start_event(self->verbose, 1);
    sleep(1);
    TP_verbose_trace_stop_event(self->verbose);
  }
  TP_solver_finalize(self);

  TP_solver_destroy(self);

  return 0;
}
