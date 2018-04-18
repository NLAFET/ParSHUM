#ifndef   _TP_SCHUR_MATRIX_H 
#define   _TP_SCHUR_MATRIX_H 

#include <omp.h>

#include "TP_matrix.h"
#include "TP_U_matrix.h"
#include "TP_L_matrix.h"
#include "TP_dense.h"
#include "TP_verbose.h"
#include "TP_internal_mem.h"

typedef struct _TP_schur_matrix *TP_schur_matrix;

struct _TP_schur_matrix {
  CSC_struct *CSC;
  CSR_struct *CSR;

  omp_lock_t *row_locks;
  omp_lock_t *col_locks;

  TP_internal_mem internal_mem;
  TP_verbose verbose;

  int **data_struct;
  int *base;

  int nb_threads;
  
  int n;
  int m;
  long nnz;
  
  int debug;
  
  double extra_space;
  double extra_space_inbetween;
};

TP_schur_matrix TP_schur_matrix_create();

void TP_schur_matrix_allocate(TP_schur_matrix self, int n, int m, long nnz, int debug, TP_verbose verbose, 
			      int nb_threads, double extra_space, double extra_space_inbetween);

void TP_schur_matrix_copy(TP_matrix A, TP_schur_matrix self, double value_tol);

void TP_schur_matrix_print(TP_schur_matrix self, char *mess);

void TP_schur_get_singletons(TP_schur_matrix self, int done_pivots, int previous_step_pivots,
			     int *nb_col_singletons, int *nb_row_singletons,
			     int *cols, int *rows, int *distributions,
			     int nb_done_pivots, int *col_perm, int *row_perm,
			     int *invr_col_perm, int *invr_row_perm);

void TP_schur_matrix_update_LD(TP_schur_matrix S, TP_L_matrix L, TP_U_matrix U, TP_matrix D,
			       int *row_perm, int *col_perm,  int nb_pivots,
			       int *invr_row_perm, int nb_row_singeltons,
			       int nb_col_singeltons, void **workspace);

void TP_schur_matrix_update_U(TP_schur_matrix S, TP_U_matrix U, TP_matrix L, 
			      int nb_pivots, int *row_perm,
			      TP_U_struct *U_struct, int U_new_n, int U_new_nnz);

void TP_schur_matrix_update_S(TP_schur_matrix S, TP_L_matrix L, TP_U_matrix U,
			      int *U_struct, int U_new_n, int *invr_row_perm,
			      int nb_pivots, int *row_perms, void **workspace,
			      double value_tol);

void 
TP_schur_matrix_update_S_rows(TP_schur_matrix S, int *L_struct, int L_new_n,
			      int L_new_nnz, int *invr_col_perm, int nb_pivots,
			      int *row_perms, int done_pivots, void **workspace);
  
void TP_schur_matrix_destroy(TP_schur_matrix self);

TP_dense_matrix  TP_schur_matrix_convert(TP_schur_matrix self, int done_pivots);

void TP_schur_check_doubles(TP_schur_matrix sefl);

void  TP_schur_matrix_check_pivots(TP_schur_matrix self,
				   int *row_perms, int *col_perms,
				   int *invr_row_perms, int *invr_col_perms,
				   int nb_pivots);

void  TP_schur_matrix_memory_check(TP_schur_matrix self);

void  TP_schur_matrix_check_symetry(TP_schur_matrix self);
void  TP_check_current_counters(TP_schur_matrix self,
				int *col_perm, int *row_perm, int nb_perms, 
				int *col_count, int *row_count, int base);

#endif //  _TP_SCHUR_MATRIX_H 

