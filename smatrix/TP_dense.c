#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <plasma.h>

#include "TP_matrix.h"
#include "TP_auxiliary.h"
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
TP_vector_read_mtl_file(TP_vector self, char *filename)
{
  FILE *file;
  int i, unused, n = self->n;
  double *vect;
  
  vect = self->vect;

  file = fopen(filename, "r");
  if (!file)
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"error while opening the vector file");
  
  for ( i = 0; i < n; i++) {
    int index;
    if ( fscanf(file, "%d %d %lf\n", &index, &unused, &vect[i]) != 3)
      TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"error while reading the vector file");
    if (index != (i + 1) ) 
      TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"the indices are wrong in the vector file");
  }
}

void 
TP_vector_read_file(TP_vector self, char *filename)
{
  char *file_ext = strrchr(filename, '.');
  if (!strcmp(file_ext, ".mtl"))
    TP_vector_read_mtl_file(self, filename);
  else 
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"unsupported vector file");
}

void
TP_vector_permute(TP_vector self, int *perms)
{
  int     i;
  int     n    = self->n;
  double *vect = self->vect, tmp[n];
  
  memcpy(tmp, vect, (size_t) n*sizeof(*vect));
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

double 
TP_vector_2norm(TP_vector self)
{
  int i;
  double res = 0.00;
  for(i=0; i<self->n; i++)
    res += self->vect[i] * self->vect[i];
  return sqrt(res);
}

void
TP_vector_copy(TP_vector src, TP_vector dst)
{
  int n_src = src->n, n_dst = dst->n;

  if (n_src != n_dst)
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"Vectors are not the same size");
  
  memcpy(dst->vect, src->vect, (size_t) n_src*sizeof(*dst->vect));
}

void
TP_vector_add(TP_vector A, double alpha, TP_vector B, double beta, TP_vector C)
{
  int n  = C->n, i;
  for(i=0; i<n; i++)
    C->vect[i] = alpha * A->vect[i] + beta * B->vect[i];
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

/* ******************************************************************************** */
/* ****************************DENSE MATRIX************************************* */
/* ******************************************************************************** */

TP_dense_matrix 
TP_dense_matrix_create(int n, int m)
{
  TP_dense_matrix self = calloc((size_t) 1,  sizeof(*self));

  self->val = calloc((size_t) n * m, sizeof(*self->val));
  self->original_rows = malloc((size_t) m * sizeof(*self->original_rows));
  self->original_cols = malloc((size_t) n * sizeof(*self->original_cols));
  self->pivots        = malloc((size_t) m * sizeof(*self->pivots));
  self->n = n;
  self->m = m;

  return self;
}


void
TP_dense_matrix_factorize(TP_dense_matrix self, int nb_threads)
{
  omp_set_num_threads(nb_threads);
  int ret = plasma_dgetrf(self->m, self->n, self->val, self->m, self->pivots);  
  if (ret) 
    printf("The factorization of the dense schur is completed, but the entry U(%d,%d) has a zero on it.\n", ret, ret);
}


int *
TP_dense_get_row_perms(TP_dense_matrix self)
{
  int m = self->m, i;
  int *row_perms = malloc((size_t) m * sizeof(*row_perms));
  
  memcpy(row_perms, self->original_rows, (size_t) m * sizeof(*row_perms));
  for(i = 0; i < m; i++) {
    int tmp = row_perms[i]; 
    /*  plasma is in freaking fortran  */
    int perm  = self->pivots[i] - 1 ;
    row_perms[i]    = row_perms[perm];
    row_perms[perm] = tmp;
  }
  return row_perms;
}

void
TP_dense_matrix_print(TP_dense_matrix self, char *mess)
{
  int n = self->n, m = self->m, row, col;
 
  printf("%s\n", mess);
  for(row = 0; row < m; row++)
    {
      for(col = 0; col < n; col++) 
	printf("[%d,%d]:(%3f)\t", row, col, self->val[col*n + row]);
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





/* ******************************************************************************** */
/* ****************************DENSE 2D MATRIX************************************* */
/* ******************************************************************************** */

TP_dense_2D
TP_dense_2D_create(int n, int m)
{
  if ( n < 0 || m < 0)
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"acquiring negative amount of memory ");
  TP_dense_2D self = calloc(1, sizeof(*self));
  self->n = n;
  self->m = m;
  self->vals = malloc((size_t) m * sizeof(*self->vals));
  for(int i = 0; i < m; i++)
    self->vals[i] = calloc((size_t) n, sizeof(**self->vals));
  
  return self;
}


TP_dense_2D
TP_dense_2D_permute(TP_dense_2D A, int *row_perm, int *col_perm)
{
  TP_dense_2D self = TP_dense_2D_create(A->n, A->m);
  int i, j;
  for(i = 0; i < self->m; i++) 
    for(j = 0; j < self->n; j++) 
      self->vals[i][j] = A->vals[row_perm[i]][col_perm[j]];

  TP_dense_2D_destroy(A);
  return self;
}

void
TP_dense_2D_facto(TP_dense_2D self)
{
  int i, j, k, n = self->n;
  double **vals = self->vals;
  if ( self->n != self->m)  
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"matrix must be square");

  for( i = 0; i < n; i++) 
    {
      double pivot = vals[i][i];
      for( j = i+1 ; j < n; j++) {
	vals[j][i] /= pivot;
	double L_ij = vals[j][i];
	for( k = i+1; k < n; k++)
	  vals[j][k] =  vals[j][k] - L_ij * vals[i][k] ;
      }
    }
}


void
TP_dense_2D_solve(TP_dense_2D self, TP_vector RHS)
{
  int i, j, n = self->n;
  double **vals = self->vals;
  double  *rhs = RHS->vect;

  if ( self->n != self->m)  
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"matrix must be square");

  for( i = 0; i < n; i++)
    for( j = i + 1; j < n; j++)
      rhs[j] -= vals[j][i] * rhs[i];

  for( i = n-1; i >= 0; i--) 
    {
      for( j = i+1; j < n; j++)
	rhs[i] -= vals[i][j] * rhs[j];
      rhs[i] /= vals[i][i];
    }
}

void 
TP_dense_2D_print(TP_dense_2D self, char *mess)
{
  int i, j, n = self->n, m = self->m;
  double **vals = self->vals;
  fprintf(stdout, "%s\n", mess);
  for(i = 0; i < m; i++) 
    {
      for(j = 0; j < n; j++)
	fprintf(stdout,"[%d,%d]: (%f)\t", i, j, vals[i][j]);
      fprintf(stdout, "\n");
    }
}

void
TP_dense_2D_destroy(TP_dense_2D self)
{
  int i;
  for(i=0; i<self->m; i++)
    free(self->vals[i]);
  free(self->vals);
  free(self);
}
