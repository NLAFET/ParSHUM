#ifndef   _TP_SCHUR_MATRIX_H 
#define   _TP_SCHUR_MATRIX_H 

#include "TP_matrix.h"
#include "TP_dense.h"


typedef struct _TP_schur_matrix *TP_schur_matrix;
typedef struct _free_space_CSC  *free_space_CSC;
typedef struct _free_space_CSR  *free_space_CSR;


struct CSR_struct {
  int nb_elem;
  int nb_free;
  int *col;
};


struct CSC_struct {
  double col_max;
  int    nb_elem;
  int    nb_free;
  int    *row;
  double *val;
};


struct _TP_schur_matrix {
  struct CSC_struct *CSC;
  struct CSR_struct *CSR;
  free_space_CSC unused_CSC;
  free_space_CSR unused_CSR;

  int n;
  int m;
  double **val;
  int    **row;
  int    **col;
  
  int nb_CSC_memories;
  int nb_CSR_memories;
  double extra_space;
  double extra_space_inbetween;
};


TP_schur_matrix TP_schur_matrix_create();

void TP_schur_matrix_allocate(TP_schur_matrix self, int n, int m, long nnz, 
			      double extra_space, double extra_space_inbetween);

void TP_schur_matrix_copy(TP_matrix A, TP_schur_matrix self);

void TP_schur_matrix_print(TP_schur_matrix self, char *mess);

void TP_schur_matrix_update_LD(TP_schur_matrix S, TP_matrix L, TP_matrix D,
			       int *row_perm, int *col_perm, int nb_pivots);

void TP_schur_matrix_update_U(TP_schur_matrix S, TP_matrix U,
			      int *row_perm, int *col_perm, int nb_pivots);

void TP_schur_matrix_update_S(TP_schur_matrix S, TP_matrix L, TP_matrix U,
			      int start, int end);

void TP_schur_matrix_destroy(TP_schur_matrix self);

TP_dense_matrix  TP_schur_matrix_convert(TP_schur_matrix self, int done_pivots);

void   TP_schur_matrix_check_perms(TP_schur_matrix self, int *row_perms, 
				   int *col_perms, int nb_pivots);

#endif //  _TP_SCHUR_MATRIX_H 
