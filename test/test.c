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
  TP_vector X, rhs, sol;
  int i;
  
  self = TP_solver_create();
  TP_solver_parse_args(self, argc, argv);
  TP_solver_read_matrix(self);
  TP_solver_init(self);
  
  X   = TP_vector_create(self->A->n);
  rhs = TP_vector_create(self->A->n);
  sol = TP_vector_create(self->A->n);
  
  for(i=0; i<X->n; i++)
    X->vect[i] = (1+i) * 200;
  TP_vector_memset(rhs, 0.0);
  TP_matrix_SpMV(self->A, X, rhs);
  TP_vector_copy(rhs, sol);
  
  TP_solver_factorize(self);
  
  TP_solver_solve(self, sol);
  
  TP_solver_copmpute_norms(self, X, sol, rhs);

  TP_solver_finalize(self);

  TP_vector_destroy(X);
  TP_vector_destroy(rhs);
  TP_vector_destroy(sol);
  TP_solver_destroy(self);

  return 0;
}
