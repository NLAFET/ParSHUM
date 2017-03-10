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



double get_swtime()
{
  struct timeval tp;
  /* struct timezone tzp; */

  gettimeofday(&tp, NULL);
  return  ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

double get_uwtime() 
{
  struct timeval t;
  /* struct timezone tzp; */

  gettimeofday(&t,NULL);
  return (double) t.tv_sec*1e6+t.tv_usec;
}

int 
main(int argc, char **argv)
{
  TP_solver self;
  int retval = 0; 
  TP_vector X, rhs;

  double t_sa;
  int i;
  
  self = TP_solver_create();
  TP_solver_init(self, argc, argv);
  
  X   = TP_vector_create(self->A->n);
  rhs = TP_vector_create(self->A->n);
  
  TP_vector_memset(X, 1.0);
  TP_vector_memset(rhs, 0.0);
  TP_matrix_SpMV(self->A, X, rhs);
  /* TP_vector_memset(X, 10.0); */

  // Factor
  t_sa = get_uwtime();
  TP_solver_factorize(self);
  double t_facto = get_uwtime() - t_sa;

  // Solve
  TP_vector sol;
  sol = TP_vector_create(self->A->n);
  TP_vector_memset(sol, 0.0);
  // init sol = b
  for (i = 0; i < self->A->n; i++)
     sol->vect[i] = rhs->vect[i];

  TP_solver_solve(self, sol, sol);

  // Compute errors

  // Compute forward error in inf norm
  double fwd_err = 0.0;
  for (i = 0; i < self->A->n; i++) {
     /* printf("x-x0[i] %e\n", abs(X->vect[i] - sol->vect[i])); */
     fwd_err = fmax(fwd_err, fabs(X->vect[i] - sol->vect[i]));
  }

  // Compute backward error in inf norm

  // Compute residual and norm of b
  TP_vector res_vec;
  double res = 0.0, x_norm = 0.0, b_norm = 0.0;

  res_vec = TP_vector_create(self->A->m);
  
  TP_matrix_SpMV(self->A, sol, res_vec);
  for (i = 0; i < self->A->m; i++) {
     res_vec->vect[i] = res_vec->vect[i] - rhs->vect[i];
     res = fmax(res, res_vec->vect[i]);
     b_norm = fmax(b_norm, rhs->vect[i]);
  }

  // compute norm of x
  for (i = 0; i < self->A->n; i++)
     x_norm = fmax(x_norm, sol->vect[i]);

  // Compute inf norm of A
  TP_vector row_norm;
  row_norm = TP_vector_create(self->A->n);
  TP_vector_memset(row_norm, 0.0);

  for (i = 0; i < self->A->col_ptr[self->A->n]; i++)
     row_norm->vect[self->A->row[i]] = row_norm->vect[self->A->row[i]] + fabs(self->A->val[i]); 

  double A_norm = 0.0;
  for (i = 0; i < self->A->n; i++)
     A_norm = fmax(A_norm, row_norm->vect[i]);

  // Scaled residual
  double tmp = A_norm * x_norm + b_norm;
  double bwd_err = res / tmp;

  TP_vector_destroy(sol);

  printf("[tp_test] Time for factorize (s): %e\n", t_facto/1e6);

  printf("[tp_test] Forward error || ||_inf : %e\n", fwd_err);
  printf("[tp_test] Backward error || ||_inf : %e\n", bwd_err);
  
  TP_solver_finalize(self);
  
  TP_vector_destroy(X);
  TP_vector_destroy(rhs);
  TP_solver_destroy(self);
  
  return retval;
}
