#ifndef   _TP_SOLVER_H
#define   _TP_SOLVER_H

#include <stdlib.h>
#include "TP_matrix.h" 
#include "TP_schur_matrix.h" 
#include "TP_dense.h"
#include "TP_verbose.h"

typedef struct _TP_solver *TP_solver;


struct _TP_solver {
  TP_matrix A;

  TP_matrix L;
  TP_matrix D;
  TP_matrix U;
  TP_schur_matrix S;
  TP_dense_matrix S_dense;

  int *row_perm;
  int *col_perm;
  int *random_col;
  int found_pivots;
  int done_pivots;
  int nb_threads;
  int verbosity;
  int nb_init_blocks;
  
  char *matrix_file;
  
  double value_tol;
  double marko_tol;
  double extra_space;
  double extra_space_inbetween;
  double sparse_part;
  
  TP_verbose verbose;
};

void print_sthg(int i);
void print_sthg2(int i, int j);
TP_solver TP_solver_create();
void      TP_solver_init(TP_solver self, int argc, char **argv);
void      TP_solver_read_matrix(TP_solver self);
void      TP_solver_factorize(TP_solver self);
void      TP_solver_solve(TP_solver self, TP_vector X, TP_vector rhs);
void      TP_solver_finalize(TP_solver self);
void      TP_solver_destroy(TP_solver self);

#endif // _TP_SOLVER_H
