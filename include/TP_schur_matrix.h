#ifndef   _TP_SCHUR_MATRIX_H 
#define   _TP_SCHUR_MATRIX_H 

#include <pthread.h>

#include "TP_matrix.h"
#include "TP_U_matrix.h"
#include "TP_dense.h"

typedef struct _TP_schur_matrix *TP_schur_matrix;
typedef struct _free_space  *free_space;


struct CSR_struct {
  int nb_elem;
  int nb_free;
  long offset;
};


struct CSC_struct {
  double col_max;
  int    nb_elem;
  int    nb_free;
  long offset;
};


struct _TP_schur_matrix {
  struct CSC_struct *CSC;
  struct CSR_struct *CSR;
  free_space unused_CSC;
  free_space unused_CSR;

  pthread_mutex_t *row_locks;
  pthread_mutex_t *col_locks;

  double *val;
  int    *row;
  int    *col;
  long **row_struct;

  int nb_threads;
  
  int n;
  int m;
  long nnz;
  
  int debug;

  long allocated_CSC;
  long allocated_CSR;
  
  
  double extra_space;
  double extra_space_inbetween;
};


TP_schur_matrix TP_schur_matrix_create();

void TP_schur_matrix_allocate(TP_schur_matrix self, int n, int m, long nnz, int debug,
			      int nb_threads, double extra_space, double extra_space_inbetween);

void TP_schur_matrix_copy(TP_matrix A, TP_schur_matrix self);

void TP_schur_matrix_print(TP_schur_matrix self, char *mess);

void TP_schur_matrix_update_LD(TP_schur_matrix S, TP_matrix L, TP_matrix D,
			       int *row_perm, int *col_perm, int nb_pivots);

void TP_schur_matrix_update_U(TP_schur_matrix S, TP_U_matrix U,
			      int nb_pivots, int *row_perm,
			      TP_U_struct *U_struct, int U_new_n, int U_new_nnz);

void TP_schur_matrix_update_S(TP_schur_matrix S, TP_matrix L, TP_U_matrix U,
			      TP_U_struct *U_struct, int U_new_n, int *invr_row_perm);

void TP_schur_matrix_destroy(TP_schur_matrix self);

TP_dense_matrix  TP_schur_matrix_convert(TP_schur_matrix self, int done_pivots);

void TP_schur_check_doubles(TP_schur_matrix sefl);

void  TP_schur_matrix_check_pivots(TP_schur_matrix self,
				   int *row_perms, int *col_perms,
				   int *invr_row_perms, int *invr_col_perms,
				   int nb_pivots);

void  TP_schur_matrix_memory_check(TP_schur_matrix self);

void  TP_schur_matrix_check_symetry(TP_schur_matrix self);
void  TP_print_GB(TP_schur_matrix self, char *mess);
void  TP_print_single_GB(free_space self, char *mess);
void  TP_check_current_counters(TP_schur_matrix self,
				int *col_perm, int *row_perm, int nb_perms, 
				int *col_count, int *row_count, int base);


#endif //  _TP_SCHUR_MATRIX_H 

