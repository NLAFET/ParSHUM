#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <strings.h>
#include <stdbool.h>

#include "spral_rutherford_boeing.h"
#include "spral_matrix_util.h"

#include "TP_enum.h"
#include "TP_auxiliary.h"
#include "TP_solver.h"
#include "TP_schur_matrix.h"
#include "TP_dense.h"

#include "TP_matrix.h"


TP_matrix
TP_matrix_create()
{
  TP_matrix self = calloc(1, sizeof(*self));
  return self;
}

void
TP_matrix_destroy(TP_matrix self)
{
  int    *index;
  long   *ptr;
  double *val;
  switch (self->type)  
    {
    case TP_CSC_matrix:
      ptr   = self->col_ptr;
      index = self->row;
      val   = self->val;
      break;

    case TP_CSR_matrix:
      ptr   = self->row_ptr;
      index = self->col;
      val   = self->val;
      break;

    case TP_Diag_matrix:
      free(self->val);
      free(self);
      return;

    case TP_Rutherford_matrix:
      spral_rb_free_handle(&self->handle);
       free(self);
      return;
    }

  free(index);
  free(ptr);
  free(val);
    
  free(self);
}


int
read_rutherford_boeing(TP_matrix self, const char*filename)
{
  int retval = 0; 
  enum spral_matrix_type matrix_type;
  struct spral_rb_read_options options;

  spral_rb_default_read_options(&options);
  options.values = 4;

  spral_rb_read(filename, &self->handle, &matrix_type, &self->m, &self->n, &self->col_ptr,
		&self->row, &self->val, &options, NULL, NULL, NULL);
  self->type = TP_Rutherford_matrix;

  self->nnz       = self->col_ptr[self->n];  
  self->allocated = self->col_ptr[self->n];  

  return retval;
}


void
TP_matrix_allocate(TP_matrix self, int n, int m, long nnz,
		   double extra_space, enum TP_matrix_type type)
{
  int  **index;
  long **ptr;
  double **val;
  int size_index;
  long size_ptr;

  switch (type)  
    {
    case TP_CSC_matrix:
      ptr   = &self->col_ptr;
      index = &self->row;
      val   = &self->val;
      size_index = nnz * extra_space;
      size_ptr   = n + 1;
      self->type = TP_CSC_matrix;
      break;

    case TP_CSR_matrix:
      ptr   = &self->row_ptr;
      index = &self->col;
      val   = &self->val;
      size_index = nnz * extra_space;
      size_ptr   = m + 1;
      self->type = TP_CSR_matrix;
      break;

    case TP_Diag_matrix:
      self->allocated = n * extra_space;
      self->val = calloc((size_t) n, sizeof(*self->val));
      self->type = TP_Diag_matrix;
      return;

    case TP_Rutherford_matrix:
      TP_warning(__FUNCTION__, __FILE__, __LINE__,"No memory is allocated, Rutherford Boeing matrices are allocated by SPRAL.");
      return;
    }

  self->allocated = size_index;

  *ptr   = calloc((size_t) size_ptr  , sizeof(**ptr));
  *index = calloc((size_t) size_index, sizeof(**index));
  *val   = calloc((size_t) size_index, sizeof(**val));
}


double
TP_matrix_get_val(TP_matrix A, int row, int col)
{
  int i;
  int col_start = A->col_ptr[col];
  int col_end   = A->col_ptr[col+1];
 
  for(i = col_start; i < col_end; i++)
    if (A->row[i] == row) 
      return A->val[i];

  char *mess;
  asprintf(&mess, "accesing not exisiting A[%d, %d]", row, col);
  TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, mess);
  return 0.0;
}

void 
TP_matrix_realloc(TP_matrix self, double extra_space)
{
  double **vect;
  int    **index;
  long   *allocated;

  vect      = &self->val;
  allocated = &self->allocated;

  switch (self->type)  
    {
    case TP_CSC_matrix:
    case TP_Rutherford_matrix:
      index = &self->row;
      break;
    case TP_CSR_matrix:
      index = &self->col;
      break;
    default:
      TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "unsuported format");
    }
  
  *allocated = *allocated * extra_space;
  *vect  = realloc(*vect,  (size_t) *allocated * sizeof(**vect));
  *index = realloc(*index, (size_t) *allocated * sizeof(**index));
}

void 
TP_print_matrix(TP_matrix self, char *mess)
{
  int n, i, j;
  long    *ptr;
  double  *array;
  int     *index;
  
  printf("%s\n", mess);
    
  switch ( self->type ) {
  case TP_CSC_matrix:
  case TP_Rutherford_matrix:
    printf("printing CSC matrix\n");
    ptr = self->col_ptr; array = self->val; index = self->row; n = self->n;
    break;
    
  case TP_CSR_matrix :
    printf("printing CSR matrix\n");
    ptr = self->row_ptr; array = self->val; index = self->col; n = self->m;
    break;
    
  case TP_Diag_matrix :
    printf("printing Diag matrix\n");
    for(i=0; i<self->n; i++)
     printf("(%f)\t",self->val[i]);
    printf("\n\n");
    return;
  }

  for (i = 0; i < n; i++) {
    printf("================%d======================\n ", i);
    for (j = ptr[i]; j < ptr[i+1]; j++)
      printf("%d(%f)  ", index[j], array[j]);
    printf("\n");
  }
  printf("\n");
}


int *
TP_matrix_rows_sizes(TP_matrix self)
{
  int  i, j, n;
  int  *row_sizes, *row;
  long *col_ptr;

  if ( self->type !=  TP_CSC_matrix &&  self->type != TP_Rutherford_matrix) 
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"this function is implemented for CSC matrix only");

  col_ptr = self->col_ptr;
  row     = self->row;
  n       = self->n;
  row_sizes = calloc((size_t) self->m, sizeof(*row_sizes));

  for(j = 0;  j < n; j++)
      for(i = col_ptr[j]; i < col_ptr[j+1]; i++)
	  row_sizes[row[i]]++;
  
  return row_sizes;
}

void 
TP_matrix_solve_L(TP_matrix L, TP_vector X, int *row_perms)
{
  int    col, n = L->n, i;
  double *L_val = L->val, *X_val = X->vect;
  long   *L_col_ptr = L->col_ptr;
  int    *L_row     = L->row;
  int    invr_perm[X->n];

  for(i = 0; i < X->n; i++)
    invr_perm[row_perms[i]] = i;

  for( col = 0; col < n; col++)
    {
      long col_start = L_col_ptr[col], col_end = L_col_ptr[col+1];

      for( i = col_start; i < col_end; i++)
  	X_val[invr_perm[L_row[i]]] -= L_val[i] * X_val[col];
    }
}

void 
TP_matrix_solve_UD(TP_matrix U, TP_matrix D,  TP_vector rhs , int *col_perm)
{
  int row, m = U->m, i;
  double *U_val = U->val, *rhs_val = rhs->vect;
  long   *U_row_ptr = U->row_ptr;
  int    *U_col = U->col;
  int invr_perm[rhs->n];

  for(i = 0; i < rhs->n; i++)
    invr_perm[col_perm[i]] = i;

  for( row = m-1; row >= 0; row--)
    {
      long row_start = U_row_ptr[row], row_end = U_row_ptr[row+1];

      for( i = row_start; i < row_end; i++) 
  	rhs_val[row] -= U_val[i] * rhs_val[invr_perm[U_col[i]]];
      rhs_val[row] /= D->val[row];
    }
}

void
TP_matrix_SpMV(TP_matrix A, TP_vector x, TP_vector y)
{
  int  col, i;
  int  *rows = A->row, n = A->n;
  long *col_ptr = A->col_ptr;
  double *A_val = A->val, *x_val = x->vect, *y_val = y->vect;

  switch (A->type)  
    {
    case TP_CSR_matrix:
    case TP_Rutherford_matrix:
      break;
    case TP_CSC_matrix:
    case TP_Diag_matrix:
      TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "this function does not support this format");
    }

  for(i = 0; i < n; i++)
    y_val[i] = 0.0;

  for(col = 0; col < n; col++) 
    {
      long col_start = col_ptr[col], col_end = col_ptr[col+1];
      for( i = col_start; i < col_end; i++)
	y_val[rows[i]] += A_val[i] * x_val[col];
    }
}

TP_matrix 
TP_matrix_create_random_matrix(int m, int n)
{
  int row, col;
  TP_matrix self = TP_matrix_create();
 
  self->n = n;
  self->m = m;
  self->allocated = self->nnz = n*m;
  self->col_ptr = malloc((size_t) (n+1)*sizeof(*self->col_ptr));
  self->row     = malloc((size_t) n*m  *sizeof(*self->row));
  self->val     = malloc((size_t) n*m  *sizeof(*self->val));

  self->col_ptr[0] = 0;
  /* self->col_ptr[1] = 3; */
  /* self->col_ptr[2] = 5; */
  /* self->col_ptr[3] = 7; */

  /* self->val[0] = 1.0; */
  /* self->val[1] = 1.0; */
  /* self->val[2] = 1.0; */
  /* self->val[3] = 4.0; */
  /* self->val[4] = 1.0; */
  /* self->val[5] = 1.0; */
  /* self->val[6] = 1.0; */

  /* self->row[0] = 0; */
  /* self->row[1] = 1; */
  /* self->row[2] = 2; */
  /* self->row[3] = 0; */
  /* self->row[4] = 1; */
  /* self->row[5] = 1; */
  /* self->row[6] = 2; */


  for(col = 0; col < n; col++)
    {
      self->col_ptr[col+1] = (col + 1) * m;
      for(row = 0; row < m; row++)
  	{
  	  /* self->row[col*m + row] = row; */
  	  /* self->val[col*m + row] = 1.00; */
  	  /* if (row == 0 && col == 2) */
  	  /*   self->val[col*m + row] = 0.0; */
  	  /* if (row == 0 && col == 1) */
  	  /*   self->val[col*m + row] = 4.0; */
  	  /* if (row == 2 && col == 1) */
  	  /*   self->val[col*m + row] = 0.0; */
	  
  	  self->val[col*m + row] =  ( (double)rand() / (double)RAND_MAX ) * 10;
  	  if (col == row)
  	    self->val[col*m + row] *= 4;
  	}
    }

  return self;
}
