#include <stdio.h>
#include <stdlib.h>
#include "TP_U_matrix.h"


TP_U_matrix 
TP_U_matrix_create(TP_matrix A, double extra_space)
{
  int i, n = A->n;
  TP_U_matrix self = calloc(1, sizeof(*self));
  self->col     = calloc(n, sizeof(*self->col));
  self->n       = n;
  
  for ( i = 0; i < n; i++)
    {
      int col_size = A->col_ptr[i+1] - A->col_ptr[i]; 
      self->col[i].val = malloc(col_size * sizeof(*self->col[i].val));
      self->col[i].row = malloc(col_size * sizeof(*self->col[i].row));
      self->col[i].allocated = col_size;
    }

  return self;
}

void 
TP_U_matrix_solve(TP_U_matrix U, TP_matrix D, TP_vector rhs,
		  int *col_perms, int *row_perms, int nb_dense_pivots)
{
  int col, n = U->n, i;
  double *rhs_val = rhs->vect;
  
  for( col = n-1; col >= n - nb_dense_pivots; col--)
    {
      int permuted_col = col_perms[col];
      U_col *u_col = &U->col[permuted_col];
      int nb_elem = u_col->nb_elem;
      double *u_val = u_col->val;
      int *u_rows = u_col->row;

      for( i = 0; i < nb_elem; i++)
  	rhs_val[row_perms[u_rows[i]]] -= u_val[i] * rhs_val[col];
    } 

  for( ; col >= 0; col--)
    {
      int permuted_col = col_perms[col];
      U_col *u_col = &U->col[permuted_col];
      int nb_elem = u_col->nb_elem;
      double *u_val = u_col->val;
      int *u_rows = u_col->row;

      rhs_val[col] /= D->val[col];

      for( i = 0; i < nb_elem; i++)
  	rhs_val[row_perms[u_rows[i]]] -= u_val[i] * rhs_val[col];
    }
}

/* TODO: take an argument, and doulbe allocated as long as it is smaller */
void
TP_U_col_realloc(U_col *self)
{
  self->allocated *=  2;
  self->val        = realloc(self->val, self->allocated * sizeof(*self->val));
  self->row        = realloc(self->row, self->allocated * sizeof(*self->row));
}

void
TP_U_matrix_print(TP_U_matrix self, char *mess)
{
  int n = self->n, i, j;
  fprintf(stdout, "%s\n", mess);
  fprintf(stdout, "printing the U matrix\n");

  for( i = 0; i < n; i++)
    {
      U_col *col = &self->col[i];
      int nb_elem = col->nb_elem;
      printf("================%d======================\n ", i);
      for( j = 0; j < nb_elem; j++)
	printf("%d(%e)  ", col->row[j], col->val[j]);
      printf("\n");
    }
  printf("\n");
}

void
TP_U_matrix_destroy(TP_U_matrix self)
{
  int i, n = self->n;

  for( i = 0; i < n; i++)
    {
      free(self->col[i].val);
      free(self->col[i].row);
    }
  free(self->col);
  free(self);
}
