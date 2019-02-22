#include <stdio.h>
#include <stdlib.h>
#include "ParSHUM_U_matrix.h"


ParSHUM_U_matrix 
ParSHUM_U_matrix_create(ParSHUM_matrix A, double extra_space)
{
  int i, n = A->n;
  ParSHUM_U_matrix self = calloc(1, sizeof(*self));
  self->col     = calloc((size_t) n, sizeof(*self->col));
  self->n       = n;
  
  for ( i = 0; i < n; i++)
    {
      int col_size = (int) A->col_ptr[i+1] - A->col_ptr[i]; 
      self->col[i].val = malloc((size_t) col_size * sizeof(*self->col[i].val));
      self->col[i].row = malloc((size_t) col_size * sizeof(*self->col[i].row));
      self->col[i].allocated = col_size;
      omp_init_lock(&self->col[i].lock);
    }

  return self;
}

void 
ParSHUM_U_BB_matrix_solve(ParSHUM_U_matrix U, ParSHUM_matrix D, ParSHUM_vector rhs, double *BB_rhs,
			  int *col_perms, int *row_perms, int nb_dense_pivots, int BB_cols)
{
  int col, n = U->n, i, j;
  double *rhs_val = rhs->vect;
  
  for( col = n-1, j = BB_cols - 1; col >= n - BB_cols; col--, j--)
    {
      int permuted_col = col_perms[col];
      U_col *u_col = &U->col[permuted_col];
      int nb_elem = u_col->nb_elem;
      double *u_val = u_col->val;
      int *u_rows = u_col->row;
      double tmp = BB_rhs[j];
      for( i = 0; i < nb_elem; i++)
  	rhs_val[u_rows[i]] -= u_val[i] * tmp;
    } 


  for( ; col >= n - nb_dense_pivots; col--)
    {
      int permuted_col = col_perms[col];
      U_col *u_col = &U->col[permuted_col];
      int nb_elem = u_col->nb_elem;
      double *u_val = u_col->val;
      int *u_rows = u_col->row;
      double tmp = rhs_val[row_perms[col]];
      for( i = 0; i < nb_elem; i++)
  	rhs_val[u_rows[i]] -= u_val[i] * tmp;
    } 
  
  /* printf("entering with %d\n", col); */
  for( ; col >= 0; col--)
    {
      int permuted_col = col_perms[col];
      U_col *u_col = &U->col[permuted_col];
      int nb_elem = u_col->nb_elem;
      double *u_val = u_col->val;
      int *u_rows = u_col->row;

      rhs_val[row_perms[col]] /= D->val[col];
      double tmp = rhs_val[row_perms[col]];
      for( i = 0; i < nb_elem; i++)
  	rhs_val[u_rows[i]] -= u_val[i] * tmp;
    }
}


void 
ParSHUM_U_matrix_solve(ParSHUM_U_matrix U, ParSHUM_matrix D, ParSHUM_vector rhs,
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
      double tmp = rhs_val[row_perms[col]];
      for( i = 0; i < nb_elem; i++)
  	rhs_val[u_rows[i]] -= u_val[i] * tmp;
    } 

  for( ; col >= 0; col--)
    {
      int permuted_col = col_perms[col];
      U_col *u_col = &U->col[permuted_col];
      int nb_elem = u_col->nb_elem;
      double *u_val = u_col->val;
      int *u_rows = u_col->row;

      rhs_val[row_perms[col]] /= D->val[col];
      double tmp = rhs_val[row_perms[col]];
      for( i = 0; i < nb_elem; i++)
  	rhs_val[u_rows[i]] -= u_val[i] * tmp;
    }
}

/* TODO: take an argument, and doulbe allocated as long as it is smaller */
void
ParSHUM_U_col_realloc(U_col *self)
{
  self->allocated *=  2;
  self->val        = realloc(self->val, (size_t) self->allocated * sizeof(*self->val));
  self->row        = realloc(self->row, (size_t) self->allocated * sizeof(*self->row));
}

void
ParSHUM_U_matrix_print(ParSHUM_U_matrix self, char *mess)
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
ParSHUM_U_matrix_destroy(ParSHUM_U_matrix self)
{
  int i, n = self->n;

  for( i = 0; i < n; i++)
    {
      free(self->col[i].val);
      free(self->col[i].row);
      omp_destroy_lock(&self->col[i].lock);
    }
  free(self->col);
  free(self);
}
