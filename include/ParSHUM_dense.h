#ifndef _ParSHUM_DENSE_H
#define _ParSHUM_DENSE_H

typedef  struct _ParSHUM_vector *ParSHUM_vector;
typedef  struct _ParSHUM_dense_matrix *ParSHUM_dense_matrix;
typedef  struct _ParSHUM_dense_2D *ParSHUM_dense_2D;

struct _ParSHUM_vector {
  double *vect; 
  int n;
};

struct _ParSHUM_dense_matrix {
  int n;
  int m; 
  double *val;
  int *original_rows;
  int *original_cols;
  int *pivots;
};

struct _ParSHUM_dense_2D {
  int n;
  int m;
  double **vals;
};

ParSHUM_vector ParSHUM_vector_create(int n);
void           ParSHUM_vector_read_file(ParSHUM_vector self, char *filename);
void           ParSHUM_vector_permute(ParSHUM_vector self, int *perms);
void           ParSHUM_vector_memset(ParSHUM_vector self, double val);
void           ParSHUM_vector_print(ParSHUM_vector self, char *mess);
double         ParSHUM_vector_2norm(ParSHUM_vector self);
void           ParSHUM_vector_copy(ParSHUM_vector src, ParSHUM_vector dst);
void           ParSHUM_vector_destroy(ParSHUM_vector self);

/* C = alpha A + beta B  */
void           ParSHUM_vector_add(ParSHUM_vector A, double alpha,
				  ParSHUM_vector B, double beta, ParSHUM_vector C);


ParSHUM_dense_matrix ParSHUM_dense_matrix_create(int n, int m);
void                 ParSHUM_dense_matrix_factorize(ParSHUM_dense_matrix self, int nb_threads);
int *                ParSHUM_dense_get_row_perms(ParSHUM_dense_matrix self);
void                 ParSHUM_dense_matrix_print(ParSHUM_dense_matrix self, char *mess);
void                 ParSHUM_dense_matrix_destroy(ParSHUM_dense_matrix slef);


ParSHUM_dense_2D   ParSHUM_dense_2D_create(int n, int m);
ParSHUM_dense_2D   ParSHUM_dense_2D_permute(ParSHUM_dense_2D A, int *row_perm, int *col_perm);
void               ParSHUM_dense_2D_facto(ParSHUM_dense_2D self);
void               ParSHUM_dense_2D_solve(ParSHUM_dense_2D self, ParSHUM_vector RHS);
void               ParSHUM_dense_2D_print(ParSHUM_dense_2D self, char *mess);
void               ParSHUM_dense_2D_destroy(ParSHUM_dense_2D self);

#endif // _ParSHUM_DENSE_H
