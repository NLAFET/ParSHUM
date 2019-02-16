#ifndef   _ParSHUM_SCHUR_MATRIX_H 
#define   _ParSHUM_SCHUR_MATRIX_H 

#include <omp.h>

#include "ParSHUM_matrix.h"
#include "ParSHUM_U_matrix.h"
#include "ParSHUM_dense.h"
#include "ParSHUM_verbose.h"
#include "ParSHUM_internal_mem.h"

typedef struct _ParSHUM_schur_matrix *ParSHUM_schur_matrix;

typedef struct _CSC_struct CSC_struct;
typedef struct _CSR_struct CSR_struct;

typedef struct _ParSHUM_L_matrix *ParSHUM_L_matrix;


struct _CSR_struct {
  int nb_elem;
  int nb_free;
  int *col;
};

struct _CSC_struct {
  double col_max;
  int    nb_elem;
  int    nb_free;
  int    nb_numerical_eligible;
  int    nb_eligible;

  double *val;
  int    *row;
};

struct _ParSHUM_schur_matrix {
  CSC_struct *CSC;
  CSR_struct *CSR;

  omp_lock_t *row_locks;
  omp_lock_t *col_locks;

  ParSHUM_internal_mem internal_mem;
  ParSHUM_verbose verbose;

  long alignment;
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


struct _ParSHUM_L_matrix {
  int n;
  long nnz;
  CSC_struct *col;
};

ParSHUM_schur_matrix ParSHUM_schur_matrix_create();

void ParSHUM_schur_matrix_allocate(ParSHUM_schur_matrix self, int n, int m, long nnz,
				   int debug, ParSHUM_verbose verbose, int nb_threads,
				   double extra_space, double extra_space_inbetween);

void ParSHUM_schur_matrix_copy(ParSHUM_matrix A, ParSHUM_schur_matrix self, double value_tol);

void ParSHUM_schur_matrix_print(ParSHUM_schur_matrix self, char *mess);

void ParSHUM_schur_get_singletons(ParSHUM_schur_matrix self, int done_pivots, int previous_step_pivots,
				  double val_tol, int *nb_col_singletons, int *nb_row_singletons,
				  int *cols, int *rows, int *distributions,
				  int nb_BB_cols, int *col_perm, int *row_perm,
				  int *invr_col_perm, int *invr_row_perm, void **workspace);

void ParSHUM_schur_matrix_update_LD(ParSHUM_schur_matrix S, ParSHUM_L_matrix L,
				    ParSHUM_U_matrix U, ParSHUM_matrix D,
				    int *row_perm, int *col_perm,  int nb_pivots,
				    int *invr_row_perm, int nb_row_singeltons,
				    int nb_col_singeltons, void **workspace);

void ParSHUM_schur_matrix_update_U(ParSHUM_schur_matrix S, ParSHUM_U_matrix U,
				   ParSHUM_matrix L, int nb_pivots, int *row_perm,
				   ParSHUM_U_struct *U_struct, int U_new_n, int U_new_nnz);

void
ParSHUM_schur_matrix_update_S(ParSHUM_schur_matrix S, ParSHUM_L_matrix L, ParSHUM_U_matrix U,
			      int *U_struct, int U_new_n, int *L_struct, int L_new_n,  
			      int *row_perms, int *invr_col_perm, int *invr_row_perm, 
			      int nb_pivots, int done_pivots, double value_tol, void **workspace);
  
void ParSHUM_schur_matrix_destroy(ParSHUM_schur_matrix self);

ParSHUM_dense_matrix 
ParSHUM_schur_matrix_convert(ParSHUM_schur_matrix S, int done_pivots, 
			     int *col_perm, int *invr_col_perm,
			     int *row_perm, int *invr_row_perm);


void ParSHUM_schur_check_doubles(ParSHUM_schur_matrix sefl);

void  ParSHUM_schur_matrix_check_pivots(ParSHUM_schur_matrix self,
					int *row_perms, int *col_perms,
					int *invr_row_perms, int *invr_col_perms,
					int nb_pivots);

void  ParSHUM_schur_matrix_memory_check(ParSHUM_schur_matrix self);

void  ParSHUM_schur_matrix_check_symetry(ParSHUM_schur_matrix self);
void  ParSHUM_check_current_counters(ParSHUM_schur_matrix self,
				     int *col_perm, int *row_perm, int nb_perms, 
				     int *col_count, int *row_count, int base);

#endif //  _ParSHUM_SCHUR_MATRIX_H 

