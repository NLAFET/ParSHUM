#include "ParSHUM_SBBD.h"
#include <mpi.h>

int
main(int argc, char **argv)
{
  ParSHUM_SBBD self;
  ParSHUM_vector X, B, SOL;
  ParSHUM_MPI_info MPI_info;
  
  self = ParSHUM_SBBD_create(0);

  ParSHUM_SBBD_parse_args(self, argc, argv);
  ParSUHM_SBBD_read_matrix(self);
  MPI_info = self->MPI_info;
  /* if (! MPI_info->rank) { */
  /*   X    = ParSHUM_vector_create(self->input_A->n); */
  /*   SOL  = ParSHUM_vector_create(self->input_A->m); */
  /*   B    = ParSHUM_vector_create(self->input_A->m); */
    
  /*   ParSHUM_vector_read_file(X, self->solver->exe_parms->RHS_file); */
    
  /*   /\* ParSHUM_vector_print(X, "X"); *\/ */
  /*   ParSHUM_matrix_SpMV(self->input_A, X, SOL); */
  /*   /\* ParSHUM_vector_print(SOL, "RHS"); *\/ */

  /*   ParSHUM_vector_copy(SOL, B); */
  /*   /\* ParSHUM_vector_print(B, "B"); *\/ */
  /* } */
  ParSHUM_SBBD_partition(self);

  /* ParSHUM_SBBD_factorize(self); */

  /* ParSHUM_SBBD_solve(self, B); */

  /* if (! MPI_info->rank) { */
  /*   self->solver->A = self->input_A; */
  /*   ParSHUM_solver_compute_norms(self->solver, B, SOL); */
  /* } */
  /* if (! MPI_info->rank)  */
  /*   ParSHUM_vector_print(B, "after solve"); */

  /* ParSHUM_SBBD_destroy(solver); */

  MPI_Finalize();

  return 0;
}
