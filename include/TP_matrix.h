#ifndef   _TP_MATRIX_H 
#define   _TP_MATRIX_H 

#include "TP_dense.h"
#include "TP_enum.h"

typedef struct _TP_matrix *TP_matrix;

struct _TP_matrix {
  enum TP_matrix_type type;
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

TP_matrix TP_matrix_create();
int       TP_matrix_init_internal_struct(TP_matrix self);
TP_matrix TP_matrix_create_random_matrix(int m, int n);

int       read_rutherford_boeing(TP_matrix self, const char*filename);
void      TP_matrix_allocate(TP_matrix self, int n, int m, long nnz, double extra_space, enum TP_matrix_type type);

double    TP_matrix_get_val(TP_matrix A, int row, int col);
int      *TP_matrix_rows_sizes(TP_matrix self);
void      TP_print_matrix(TP_matrix self, char *mess);
void      TP_matrix_realloc(TP_matrix self, double extra_space);

void      TP_matrix_solve_L(TP_matrix L, TP_vector X, int *row_perms);
void      TP_matrix_solve_UD(TP_matrix U, TP_matrix D,  TP_vector rhs , int *col_perm);

void      TP_matrix_SpMV(TP_matrix A, TP_vector x, TP_vector y);
void      TP_matrix_destroy(TP_matrix self);

#endif // _TP_MATRIX_H 
