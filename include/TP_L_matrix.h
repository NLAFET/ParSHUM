#ifndef _TP_L_MATRIX
#define _TP_L_MATRIX

#include "TP_dense.h"

typedef struct _TP_L_matrix *TP_L_matrix;

struct _TP_L_matrix {
  long *start;
  long *end;
  int *row;
  double *val;

  int n;
  long nnz;
  long allocated;
};

TP_L_matrix  TP_L_matrix_create(int n, long ini_nnz);
void         TP_L_matrix_realloc(TP_L_matrix self);
void         TP_L_matrix_solve(TP_L_matrix L, TP_vector RHS, int *perms);
void         TP_L_matrix_print(TP_L_matrix self, char *mess);
void         TP_L_matrix_destroy(TP_L_matrix self);

#endif //_TP_L_MATRIX
