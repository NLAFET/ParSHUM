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
  ParSHUM_vector X, B, SOL;

  /* Create the solver */
  self = ParSHUM_solver_create();

  /* Parse the arguments */
  ParSHUM_solver_parse_args(self, argc, argv, 1);
  /* Read the matrix */
  ParSHUM_solver_read_matrix(self);
  
  /* Initialize the vectors */
  X    = ParSHUM_vector_create(self->A->n);
  SOL  = ParSHUM_vector_create(self->A->m);  
  B    = ParSHUM_vector_create(self->A->m);  
  ParSHUM_vector_read_file(X, self->exe_parms->RHS_file); 

  /* ParSHUM_vector_print(X, "the real X"); */
  ParSHUM_matrix_SpMV(self->A, X, SOL);

  /* copy the vector SOL in B */
  ParSHUM_vector_copy(SOL, B);
  
  /* Initialize the solver */
  ParSHUM_solver_init(self);

  /* Perform the factorization */
  ParSHUM_solver_factorize(self);
  
  /* Perform the solve operation */
  ParSHUM_solver_solve(self, B);
  /* ParSHUM_vector_print(B, "computed solution"); */

  /* compute the norms */
  ParSHUM_solver_compute_norms(self, B, SOL);

  /* Finalize the solver */
  ParSHUM_solver_finalize(self);

  /* Free all the data */
  ParSHUM_vector_destroy(X);
  ParSHUM_vector_destroy(B);
  ParSHUM_solver_destroy(self);


  return 0;
}
