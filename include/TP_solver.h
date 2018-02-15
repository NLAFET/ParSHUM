#ifndef   _TP_SOLVER_H
#define   _TP_SOLVER_H

#include <stdlib.h>
#include <pthread.h>
#include "TP_enum.h"
#include "TP_matrix.h" 
#include "TP_U_matrix.h" 
#include "TP_Luby.h" 
#include "TP_schur_matrix.h" 
#include "TP_dense.h"
#include "TP_verbose.h"

typedef struct _TP_solver *TP_solver;
typedef struct  _TP_pivot_candidates *TP_pivot_candidates;
typedef struct  _TP_counters *TP_counters;

struct _TP_pivot_candidates
{
  int *row;
  int *marko;
  int *best_marko;
};

struct _TP_counters
{
  int *array;
  int *used_counters;
  int nb_used_counters;
};

struct _TP_solver {
  TP_matrix   A;
  TP_dense_2D A_debug;

  TP_matrix   L;
  TP_matrix   D;
  TP_U_matrix U;
  TP_schur_matrix S;
  TP_dense_matrix S_dense;
  
  TP_pivot_candidates candidates;
  TP_counters *counters;

  TP_Luby Luby;

  TP_U_struct *U_struct;
  TP_U_struct *L_struct;

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
void      TP_solver_alloc_counters(TP_solver solver, int **col_count, int **row_count);
void      TP_solver_dealloc_counters(TP_solver solver, int *col_count, int *row_count);

void      TP_solver_finalize(TP_solver self);
void      TP_solver_destroy(TP_solver self);

#endif // _TP_SOLVER_H
