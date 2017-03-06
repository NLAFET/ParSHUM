#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <plasma.h>

#include "TP_matrix.h"
#include "TP_dense.h"

TP_vector
TP_vector_create(int n)
{
  TP_vector self = malloc((size_t) sizeof(*self));
  self->vect     = malloc((size_t) n*sizeof(*self->vect));
  self->n        = n;
  return self;
}

void
TP_vector_permute(TP_vector self, int *perms)
{
  int     i;
  int     n    = self->n;
  double *vect = self->vect, tmp[n];
  
  memcpy(tmp, vect, n*sizeof(*vect));
  for(i = 0; i < n; i++) 
    vect[i] = tmp[perms[i]];
}

void 
    TP_vector_memset(TP_vector self, double val)
{
  int     n    = self->n, i;
  double *vals = self->vect;
  for( i=0; i < n; i++)
    vals[i] = val;
}

void
TP_vector_print(TP_vector self, char *mess)
{
  int n = self->n, i;
  double *vals = self->vect;
  printf("%s\n", mess);
  for( i=0; i < n; i++)
    printf("%d:(%f)  ", i, vals[i]);
  printf("\n");
}

void 
TP_vector_destroy(TP_vector self)
{
  free(self->vect);
  free(self);
}

TP_dense_matrix 
TP_dense_matrix_create(int n, int m)
{
  TP_dense_matrix self = calloc((size_t) 1,  sizeof(*self));
  self->val = calloc((size_t) n*m, sizeof(*self->val));
  self->original_rows = malloc((size_t) m * sizeof(*self->original_rows));
  self->original_cols = malloc((size_t) n * sizeof(*self->original_cols));
  self->pivots        = malloc((size_t) m * sizeof(*self->pivots));
  self->n = n;
  self->m = m;

  return self;
}


void
TP_dense_matrix_factorize(TP_dense_matrix self)
{
  int ret = plasma_dgetrf(self->m, self->n, self->val, self->m, self->pivots);  
  if (ret) 
    printf("The factorization of the dense schur is completed, but the entry U(%d,%d) has a zero on it.\n", ret, ret);
}

void
TP_dense_matrix_print(TP_dense_matrix self, char *mess)
{
  int n = self->n, m = self->m, row, col;
 
  printf("%s\n", mess);
  for(col = 0; col < n; col++)
    {
      printf("================%d======================\n", col);
      for(row = 0; row < m; row++) 
	printf("%d:(%f)  ", row, self->val[col*n + row]);
      printf("\n");
    }
  printf("\n");
}


void 
TP_dense_matrix_destroy(TP_dense_matrix self)
{
  free(self->original_rows);
  free(self->original_cols);
  free(self->pivots);
  free(self->val);
  free(self);
}
