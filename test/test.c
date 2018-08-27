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
  ParSHUM_vector X, B;
  
  /* Create the solver */
  self = ParSHUM_solver_create();

  /* Parse the arguments */
  ParSHUM_solver_parse_args(self, argc, argv, 1);
  /* Read the matrix */
  ParSHUM_solver_read_matrix(self);

  
  /* Initialize the vectors */
  X = ParSHUM_vector_create(self->A->n);
  B = ParSHUM_vector_create(self->A->n);  
  ParSHUM_vector_read_file(X, self->exe_parms->RHS_file);

  /* copy the vector X in B */
  ParSHUM_vector_copy(X, B);
  
  /* Initialize the solver */
  ParSHUM_solver_init(self);

  /* Perform the factorization */
  ParSHUM_solver_factorize(self);
  
  /* Perform the solve operation */
  ParSHUM_solver_solve(self, B);

  /* Compute the norms */
  ParSHUM_solver_compute_norms(self, B, X);

  /* Finalize the solver */
  ParSHUM_solver_finalize(self);

  /* Free all the data */
  ParSHUM_vector_destroy(X);
  ParSHUM_vector_destroy(B);
  ParSHUM_solver_destroy(self);

  return 0;
}
