#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>
#include <math.h>

#include "ParSHUM_solver.h" 
#include "ParSHUM_matrix.h" 
#include "ParSHUM_dense.h"
#include "ParSHUM_enum.h"
#include "ParSHUM_pivot_list.h" 
#include "ParSHUM_auxiliary.h"


int
main(int argc, char **argv) 
{
  ParSHUM_solver self;
  
  self = ParSHUM_solver_create();
  ParSHUM_solver_parse_args(self, argc, argv, 1);
  ParSHUM_solver_read_matrix(self);
  ParSHUM_solver_init(self);
  printf("%d !\n",omp_get_thread_num());

#pragma omp parallel num_threads(self->exe_parms->nb_threads)
  {
    ParSHUM_verbose_trace_start_event(self->verbose, 0);
    sleep(1);
    ParSHUM_verbose_trace_stop_event(self->verbose);
    sleep(2);
    ParSHUM_verbose_trace_start_event(self->verbose, 1);
    sleep(1);
    ParSHUM_verbose_trace_stop_event(self->verbose);
  }
  ParSHUM_solver_finalize(self);

  ParSHUM_solver_destroy(self);

  return 0;
}
