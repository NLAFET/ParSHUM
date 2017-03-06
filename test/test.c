#include <stdio.h>
#include <string.h>
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
  int retval = 0; 
  TP_vector X, rhs;

  
  self = TP_solver_create();
  TP_solver_init(self, argc, argv);
  
  X   = TP_vector_create(self->A->n);
  rhs = TP_vector_create(self->A->n);
  
  TP_vector_memset(X, 100.0);
  TP_matrix_SpMV(self->A, X, rhs);
  TP_vector_memset(X, 10.0);

  TP_solver_factorize(self);
  TP_solver_solve(self, rhs, rhs);
  
  TP_solver_finalize(self);
  
  TP_vector_destroy(X);
  TP_vector_destroy(rhs);
  TP_solver_destroy(self);
  
  return retval;
}
