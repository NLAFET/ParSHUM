#ifndef   _TP_SOLVER_H
#define   _TP_SOLVER_H

#include <stdlib.h>
#include "TP_enum.h"
#include "TP_matrix.h" 
#include "TP_schur_matrix.h" 
#include "TP_pivot_list.h" 
#include "TP_dense.h"
#include "TP_verbose.h"

typedef struct _TP_solver *TP_solver;

struct _TP_solver {
  TP_matrix A;
  TP_dense_2D A_debug;

  TP_matrix L;
  TP_matrix D;
  TP_matrix U;
  TP_schur_matrix S;
  TP_dense_matrix S_dense;
  TP_pivot_candidates candidates;

  int *row_perm;
  int *col_perm;
  int *invr_row_perm;
  int *invr_col_perm;
  int *random_col;
  int *previous_pivots;

  int found_pivots;
  int done_pivots;
  int verbosity;
  
  int debug;

  TP_exe_parms exe_parms;
  TP_verbose verbose;
};

TP_solver TP_solver_create();
void      TP_solver_init(TP_solver self);
void      TP_solver_read_matrix(TP_solver self);
void      TP_solver_factorize(TP_solver self);
void      TP_solver_solve(TP_solver self, TP_vector rhs);
void      TP_solver_iterative_refinement(TP_solver self, TP_vector X, 
					 TP_vector rhs, double wanted_precision);
void      TP_solver_parse_args(TP_solver self, int argc, char **argv);
void      TP_solver_copmpute_norms(TP_solver self,       TP_vector X,
				   TP_vector X_computed, TP_vector rhs);
void      TP_solver_finalize(TP_solver self);
void      TP_solver_destroy(TP_solver self);

#endif // _TP_SOLVER_H
