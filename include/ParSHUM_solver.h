#ifndef   _ParSHUM_SOLVER_H
#define   _ParSHUM_SOLVER_H

#include <stdlib.h>
#include <pthread.h>
#include "ParSHUM_enum.h"
#include "ParSHUM_matrix.h" 
#include "ParSHUM_U_matrix.h" 
#include "ParSHUM_L_matrix.h" 
#include "ParSHUM_Luby.h" 
#include "ParSHUM_schur_matrix.h" 
#include "ParSHUM_dense.h"
#include "ParSHUM_verbose.h"

typedef struct _ParSHUM_solver *ParSHUM_solver;
typedef struct  _ParSHUM_pivot_candidates *ParSHUM_pivot_candidates;
typedef struct  _ParSHUM_counters *ParSHUM_counters;

struct _ParSHUM_pivot_candidates
{
  int *row;
  int *marko;
  int *best_marko;
};

struct _ParSHUM_counters
{
  int *array;
  int *used_counters;
  int nb_used_counters;
};

struct _ParSHUM_solver {
  ParSHUM_matrix   A;
  ParSHUM_dense_2D A_debug;

  ParSHUM_L_matrix L;
  ParSHUM_matrix   D;
  ParSHUM_U_matrix U;
  ParSHUM_schur_matrix S;
  ParSHUM_dense_matrix S_dense;
  
  ParSHUM_pivot_candidates candidates;
  ParSHUM_counters *counters;

  ParSHUM_Luby Luby;

  ParSHUM_U_struct *U_struct;
  ParSHUM_U_struct *L_struct;

  pthread_mutex_t counters_lock;

  int *row_perm;
  int *col_perm;
  int *invr_row_perm;
  int *invr_col_perm;
  int *rows_count;
  int *cols_count;
  int *random_col;
  int *previous_pivots;
  int *cols;
  int *rows;
  int *distributions;
  int *seeds;
  int *logical_cols;
  int *logical_rows;
  void **workspace;

  int nb_row_singletons;
  int nb_col_singletons;
  int found_pivots;
  int done_pivots;
  int dense_pivots;
  int verbosity;
  int nb_counters;
  int size_counters;
  int allocated_U_struct;  
  int n_U_structs;  
  int nnz_U_structs;  
  int allocated_L_struct;  
  int n_L_structs;  
  int nnz_L_structs;  
  int previous_step_pivots;

  int debug;
  int step;

  ParSHUM_exe_parms exe_parms;
  ParSHUM_verbose verbose;
};

ParSHUM_solver ParSHUM_solver_create();
void           ParSHUM_solver_init(ParSHUM_solver self);
void           ParSHUM_solver_read_matrix(ParSHUM_solver self);
void           ParSHUM_solver_factorize(ParSHUM_solver self);
void           ParSHUM_solver_solve(ParSHUM_solver self, ParSHUM_vector rhs);
/* void           ParSHUM_solver_iterative_refinement(ParSHUM_solver self, ParSHUM_vector X,  */
/*     						      ParSHUM_vector rhs, double wanted_precision); */
void           ParSHUM_solver_parse_args(ParSHUM_solver self, int argc, char **argv);
void           ParSHUM_solver_compute_norms(ParSHUM_solver self, ParSHUM_vector X, ParSHUM_vector rhs);
void           ParSHUM_solver_alloc_counters(ParSHUM_solver solver, int **col_count, int **row_count);
void           ParSHUM_solver_dealloc_counters(ParSHUM_solver solver, int *col_count, int *row_count);

void           ParSHUM_solver_finalize(ParSHUM_solver self);
void           ParSHUM_solver_destroy(ParSHUM_solver self);

#endif // _ParSHUM_SOLVER_H
