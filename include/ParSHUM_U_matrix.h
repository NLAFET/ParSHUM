#ifndef _ParSHUM_U_MATRIX
#define _ParSHUM_U_MATRIX

#include <omp.h>
#include "ParSHUM_matrix.h"

typedef struct _U_col U_col;
typedef struct _U_matrix *ParSHUM_U_matrix;
typedef struct _ParSHUM_U_struct ParSHUM_U_struct;

struct _U_col {
  int nb_elem;
  int allocated;
  int cost;

  double *val;
  int *row;
  omp_lock_t lock;
};


struct _ParSHUM_U_struct
{
  int col;
  int nb_elem;
};

struct _U_matrix {
  int n;
  long nnz;

  U_col *col;
};

ParSHUM_U_matrix ParSHUM_U_matrix_create(ParSHUM_matrix A, double extra_space);
void        ParSHUM_U_matrix_solve(ParSHUM_U_matrix U, ParSHUM_matrix D, ParSHUM_vector rhs,
			      int *col_perms, int *row_perms, int nb_dense_pivots);
void        ParSHUM_U_col_realloc(U_col *self);
void        ParSHUM_U_matrix_print(ParSHUM_U_matrix self, char *mess);
void        ParSHUM_U_matrix_destroy(ParSHUM_U_matrix self);

#endif // _ParSHUM_U_MATRIX
