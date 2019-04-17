#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <plasma.h>

#include "ParSHUM_matrix.h"
#include "ParSHUM_auxiliary.h"
#include "ParSHUM_dense.h"

ParSHUM_vector
ParSHUM_vector_create(int n)
{
  ParSHUM_vector self = malloc((size_t) sizeof(*self));
  self->vect     = malloc((size_t) n*sizeof(*self->vect));
  self->n        = n;
  return self;
}

void
ParSHUM_vector_read_mtl_file(ParSHUM_vector self, char *filename)
{
  FILE *file;
  int i, unused, n = self->n;
  double *vect;
  
  vect = self->vect;

  file = fopen(filename, "r");
  if (!file)
    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"error while opening the vector file");
  
  for ( i = 0; i < n; i++) {
    int index;
    if ( fscanf(file, "%d %d %lf\n", &index, &unused, &vect[i]) != 3)
      ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"error while reading the vector file");
    if (index != (i + 1) ) 
      ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"the indices are wrong in the vector file");
  }
}

/* TODO: rewrite this bullshit */
void 
ParSHUM_vector_read_file(ParSHUM_vector self, char *filename)
{
  char *file_ext;
  
  if (filename) {
    file_ext = strrchr(filename, '.');

    if (!strcmp(file_ext, ".mtl")) {
      ParSHUM_vector_read_mtl_file(self, filename);
    } else   {
      for(int i=0; i < self->n; i++)
	self->vect[i] = (double) (1+i);
    }
  }  else { 
    for(int i=0; i < self->n; i++)
      self->vect[i] = (double) (1+i);
  }
}

void
ParSHUM_vector_permute(ParSHUM_vector self, int *perms, int n)
{
  int     i;
  double *vect = self->vect, tmp[n];
  
  memcpy(tmp, vect, (size_t) n*sizeof(*vect));
  for(i = 0; i < n; i++) 
    vect[i] = tmp[perms[i]];
}

void 
ParSHUM_vector_memset(ParSHUM_vector self, double val)
{
  int     n    = self->n, i;
  double *vals = self->vect;
  for( i=0; i < n; i++)
    vals[i] = val;
}

double 
ParSHUM_vector_2norm(ParSHUM_vector self)
{
  int i;
  double res = 0.00;
  for(i=0; i<self->n; i++)
    res += self->vect[i] * self->vect[i];
  return sqrt(res);
}

void
ParSHUM_vector_copy(ParSHUM_vector src, ParSHUM_vector dst)
{
  int n_src = src->n, n_dst = dst->n;

  if (n_src != n_dst)
    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"Vectors are not the same size");
  
  memcpy(dst->vect, src->vect, (size_t) n_src*sizeof(*dst->vect));
}

void
ParSHUM_vector_add(ParSHUM_vector A, double alpha, ParSHUM_vector B, double beta, ParSHUM_vector C)
{
  int n  = C->n, i;
  for(i=0; i<n; i++)
    C->vect[i] = alpha * A->vect[i] + beta * B->vect[i];
}

void
ParSHUM_vector_print(ParSHUM_vector self, char *mess)
{
  int n = self->n, i;
  double *vals = self->vect;
  printf("%s\n", mess);
  for( i=0; i < n; i++)
    printf("%d:(%f)  ", i, vals[i]);
  printf("\n");
}

void 
ParSHUM_vector_destroy(ParSHUM_vector self)
{
  free(self->vect);
  free(self);
}

/* ******************************************************************************** */
/* ****************************DENSE MATRIX************************************* */
/* ******************************************************************************** */

ParSHUM_dense_matrix 
ParSHUM_dense_matrix_create(int n, int m)
{
  ParSHUM_dense_matrix self = calloc((size_t) 1,  sizeof(*self));
  int i;
  self->val = calloc((size_t) n * m, sizeof(*self->val));
  self->pivots        = calloc((size_t) m, sizeof(*self->pivots));
  for (i = 0; i < m; i++) 
    self->pivots[i] = i+1;

  self->n = n;
  self->m = m;

  return self;
}


void
ParSHUM_dense_matrix_factorize(ParSHUM_dense_matrix self, int BB_cols, int nb_threads)
{
  int ret = 0;
  char mess[2048];

  omp_set_num_threads(nb_threads);
  if (self->n && self->m) { 
    ret = plasma_dgetrf(self->m, self->n - BB_cols, self->val, self->m, self->pivots);  
  }
  if (ret)  {
    snprintf(mess, 2048,"The factorization of the dense schur is completed, but the entry U(%d,%d) has a zero on it.\n", ret, ret);
    ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
  }

  if (BB_cols) {
    plasma_dlaswp(PlasmaRowwise, self->m, BB_cols,  &self->val[(self->n-BB_cols) * self->m], self->m,
    		  self->pivots, 1);
    
    
    plasma_dtrsm(PlasmaLeft, PlasmaLower, PlasmaNoTrans, PlasmaUnit, self->n - BB_cols, BB_cols, 1.0,
    		 self->val, self->m, &self->val[(self->n - BB_cols) * self->m], self->m);

    plasma_dgemm(PlasmaNoTrans, PlasmaNoTrans, self->m - self->n + BB_cols, BB_cols, self->n - BB_cols, -1.0, 
		 &self->val[self->n - BB_cols], self->m, &self->val[(self->n - BB_cols) * self->m], self->m, 
 		 1,  &self->val[(self->n - BB_cols) * self->m + self->n - BB_cols], self->m);
  }
}


int *
ParSHUM_dense_get_row_perms(ParSHUM_dense_matrix self, int *row_perm)
{
  int m = self->m, i;
  int *_row_perm = malloc((size_t) m * sizeof(*_row_perm));
  
  memcpy(_row_perm, row_perm, (size_t) m * sizeof(*row_perm));
  for(i = 0; i < m; i++) {
    int tmp = _row_perm[i]; 
    /*  plasma is in freaking fortran  */
    int perm  = self->pivots[i] - 1 ;
    _row_perm[i]    = _row_perm[perm];
    _row_perm[perm] = tmp;
  }
  return _row_perm;
}

void
ParSHUM_dense_matrix_get_RHS(ParSHUM_dense_matrix self, double *dense_RHS,
			     int *row_perms, double *RHS, ParSHUM_perm_type perms_type)
{
  int i,  m = self->m, n = self->n; 
  int *plasma_perm = self->pivots;
  if (perms_type != ParSHUM_perm_none)  
    for(i = 0; i < m; i++) 
      dense_RHS[i] = RHS[row_perms[i]];
  
  if (perms_type == ParSHUM_perm_both)  
    for(i = 0; i < m; i++) { 
      double tmp = dense_RHS[i];
      /*  plasma is in freaking fortran  */
      int perm   = plasma_perm[i] - 1;
      dense_RHS[i] = dense_RHS[perm];
      dense_RHS[perm] = tmp;
    }
}

void
ParSHUM_dense_matrix_update_RHS(ParSHUM_dense_matrix self, double *dense_RHS,
				int *row_perms, double *RHS)
{
  int i,  m = self->m;//, n = self->n; 
  /* int *plasma_perm = self->pivots; */

  /* for(i = n - 1; i >= 0; i--) { */
  /*   double tmp = dense_RHS[i]; */
  /*   /\*  plasma is in freaking fortran  *\/ */
  /*   int perm   = plasma_perm[i] - 1; */
  /*   dense_RHS[i] = dense_RHS[perm]; */
  /*   dense_RHS[perm] = tmp; */
  /* } */

  for(i = 0; i < m; i++) 
    RHS[row_perms[i]] = dense_RHS[i];
}

void
ParSHUM_dense_matrix_print(ParSHUM_dense_matrix self, char *mess)
{
  int n = self->n, m = self->m, row, col;
 
  printf("%s\n", mess);
  for(row = 0; row < m; row++)
    {
      for(col = 0; col < n; col++) 
	printf("[%d,%d]:(%3f)\t", col, row, self->val[col*m + row]);
      printf("\n");
    }
  printf("\n");
}


void 
ParSHUM_dense_matrix_destroy(ParSHUM_dense_matrix self)
{
  free(self->pivots);
  free(self->val);
  free(self);
}





/* ******************************************************************************** */
/* ****************************DENSE 2D MATRIX************************************* */
/* ******************************************************************************** */

ParSHUM_dense_2D
ParSHUM_dense_2D_create(int n, int m)
{
  if ( n < 0 || m < 0)
    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"acquiring negative amount of memory ");
  ParSHUM_dense_2D self = calloc(1, sizeof(*self));
  self->n = n;
  self->m = m;
  self->vals = malloc((size_t) m * sizeof(*self->vals));
  for(int i = 0; i < m; i++)
    self->vals[i] = calloc((size_t) n, sizeof(**self->vals));
  
  return self;
}


ParSHUM_dense_2D
ParSHUM_dense_2D_permute(ParSHUM_dense_2D A, int *row_perm, int *col_perm)
{
  ParSHUM_dense_2D self = ParSHUM_dense_2D_create(A->n, A->m);
  int i, j;
  for(i = 0; i < self->m; i++) 
    for(j = 0; j < self->n; j++) 
      self->vals[i][j] = A->vals[row_perm[i]][col_perm[j]];

  ParSHUM_dense_2D_destroy(A);
  return self;
}

void
ParSHUM_dense_2D_facto(ParSHUM_dense_2D self)
{
  int i, j, k, n = self->n;
  double **vals = self->vals;
  if ( self->n != self->m)  
    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"matrix must be square");

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
ParSHUM_dense_2D_solve(ParSHUM_dense_2D self, ParSHUM_vector RHS)
{
  int i, j, n = self->n;
  double **vals = self->vals;
  double  *rhs = RHS->vect;

  if ( self->n != self->m)  
    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"matrix must be square");

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
ParSHUM_dense_2D_print(ParSHUM_dense_2D self, char *mess)
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
ParSHUM_dense_2D_destroy(ParSHUM_dense_2D self)
{
  int i;
  for(i=0; i<self->m; i++)
    free(self->vals[i]);
  free(self->vals);
  free(self);
}
