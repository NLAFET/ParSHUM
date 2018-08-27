#ifndef   _ParSHUM_MATRIX_H 
#define   _ParSHUM_MATRIX_H 

#include "ParSHUM_dense.h"
#include "ParSHUM_enum.h"

typedef struct _ParSHUM_matrix *ParSHUM_matrix;

struct _ParSHUM_matrix {
  ParSHUM_matrix_type type;
  int n; 
  int m; 

  long allocated;
  long nnz;

  int *row;
  long *col_ptr; 
  double *val;

  int *col;
  long *row_ptr;

  // needed because the driver for reading the matrix 
  // done in fotran. It is additional memory for  
  // for each "fortran"  array (http://www.test-numerical.rl.ac.uk/spral/doc/sphinx/C/rutherford_boeing.html) 
  void *handle;
};

ParSHUM_matrix ParSHUM_matrix_create();
ParSHUM_matrix ParSHUM_matrix_create_random_matrix(int m, int n);
int ParSHUM_read_rutherford_boeing(ParSHUM_matrix self, const char*filename);

/* int       read_rutherford_boeing(ParSHUM_matrix self, const char*filename); */
void      ParSHUM_read_mtl_file(ParSHUM_matrix self, const char*filename);
void      ParSHUM_matrix_allocate(ParSHUM_matrix self, int n, int m, long nnz, double extra_space, ParSHUM_matrix_type type);
void      ParSHUM_matrix_copy(ParSHUM_matrix src, ParSHUM_matrix dest);

double    ParSHUM_matrix_get_val(ParSHUM_matrix A, int row, int col);
int      *ParSHUM_matrix_rows_sizes(ParSHUM_matrix self);
void      ParSHUM_matrix_print(ParSHUM_matrix self, char *mess);
void      ParSHUM_matrix_print_as_dense(ParSHUM_matrix self, char *mess);
void      ParSHUM_matrix_realloc(ParSHUM_matrix self);

void      ParSHUM_matrix_solve_L(ParSHUM_matrix L, ParSHUM_vector X, int *perms);
void      ParSHUM_matrix_solve_UD(ParSHUM_matrix U, ParSHUM_matrix D,  ParSHUM_vector rhs, int *perms);

ParSHUM_matrix ParSHUM_matrix_permute(ParSHUM_matrix A, int *col_perm, int *invr_row_perm);
void           ParSHUM_matrix_convert(ParSHUM_matrix self, ParSHUM_matrix_type type);

void      ParSHUM_matrix_SpMV(ParSHUM_matrix A, ParSHUM_vector x, ParSHUM_vector y);
double    ParSHUM_matrix_get_norm(ParSHUM_matrix self);

void      ParSHUM_print_LDU(ParSHUM_matrix A, ParSHUM_matrix L, ParSHUM_matrix D, ParSHUM_matrix U, int *row_perm);
ParSHUM_dense_2D ParSHUM_dense_2D_convert_sparse(ParSHUM_matrix A);

void      ParSHUM_matrix_destroy(ParSHUM_matrix self);

#endif // _ParSHUM_MATRIX_H 
