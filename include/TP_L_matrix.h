#ifndef _TP_L_MATRIX
#define _TP_L_MATRIX

#include "TP_dense.h"
#include "TP_internal_mem.h"

typedef struct _TP_L_matrix *TP_L_matrix;

struct _TP_L_matrix {
  int n;
  long nnz;
  CSC_struct *col;
};

TP_L_matrix  TP_L_matrix_create(int n);
void         TP_L_matrix_solve(TP_L_matrix L, TP_vector RHS, int *perms);
void         TP_L_matrix_print(TP_L_matrix self, char *mess);
void         TP_L_matrix_destroy(TP_L_matrix self);

#endif //_TP_L_MATRIX
