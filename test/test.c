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
  ParSHUM_vector X, rhs, sol;
  int i;
  
  self = ParSHUM_solver_create();
  ParSHUM_solver_parse_args(self, argc, argv);
  ParSHUM_solver_read_matrix(self);
  ParSHUM_solver_init(self);
  
  X   = ParSHUM_vector_create(self->A->n);
  rhs = ParSHUM_vector_create(self->A->n);
  sol = ParSHUM_vector_create(self->A->n);
  
  if (self->exe_parms->RHS_file)
    ParSHUM_vector_read_file(X, self->exe_parms->RHS_file);
  else 
    for(i=0; i<X->n; i++)
      X->vect[i] = (1+i) * 200;
  ParSHUM_vector_memset(rhs, 0.0);
  ParSHUM_matrix_SpMV(self->A, X, rhs);
  ParSHUM_vector_copy(rhs, sol);
  
  ParSHUM_solver_factorize(self);
  
  ParSHUM_solver_solve(self, sol);
  ParSHUM_solver_copmpute_norms(self, X, sol, rhs);

  ParSHUM_solver_finalize(self);

  ParSHUM_vector_destroy(X);
  ParSHUM_vector_destroy(rhs);
  ParSHUM_vector_destroy(sol);
  ParSHUM_solver_destroy(self);

  return 0;
}
