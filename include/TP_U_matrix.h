#ifndef _TP_U_MATRIX
#define _TP_U_MATRIX

#include <pthread.h>

#include "TP_matrix.h"

typedef struct _U_col U_col;
typedef struct _U_matrix *TP_U_matrix;
typedef struct _TP_U_struct TP_U_struct;

struct _U_col {
  int nb_elem;
  int allocated;
  int cost;

  double *val;
  int *row;
};


struct _TP_U_struct
{
  int col;
  int nb_elem;
};

struct _U_matrix {
  int n;
  long nnz;

  U_col *col;
};

TP_U_matrix TP_U_matrix_create(TP_matrix A, double extra_space);
void        TP_U_matrix_solve(TP_U_matrix U, TP_matrix D, TP_vector rhs,
			      int *col_perms, int *row_perms, int nb_dense_pivots);
void        TP_U_col_realloc(U_col *self);
void        TP_U_matrix_print(TP_U_matrix self, char *mess);
void        TP_U_matrix_destroy(TP_U_matrix self);

#endif // _TP_U_MATRIX
