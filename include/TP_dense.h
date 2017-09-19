#ifndef _TP_DENSE_H
#define _TP_DENSE_H

typedef  struct _TP_vector *TP_vector;
typedef  struct _TP_dense_matrix *TP_dense_matrix;
typedef  struct _TP_dense_2D *TP_dense_2D;

struct _TP_vector {
  double *vect; 
  int n;
};

struct _TP_dense_matrix {
  int n;
  int m; 
  double *val;
  int *original_rows;
  int *original_cols;
  int *pivots;
};

struct _TP_dense_2D {
  int n;
  int m;
  double **vals;
};

TP_vector TP_vector_create(int n);
void      TP_vector_read_file(TP_vector self, char *filename);
void      TP_vector_permute(TP_vector self, int *perms);
void      TP_vector_memset(TP_vector self, double val);
void      TP_vector_print(TP_vector self, char *mess);
double    TP_vector_2norm(TP_vector self);
void      TP_vector_copy(TP_vector src, TP_vector dst);
void      TP_vector_destroy(TP_vector self);

/* C = alpha A + beta B  */
void      TP_vector_add(TP_vector A, double alpha, TP_vector B, double beta, TP_vector C);


TP_dense_matrix TP_dense_matrix_create(int n, int m);
void            TP_dense_matrix_factorize(TP_dense_matrix self, int nb_threads);
int *           TP_dense_get_row_perms(TP_dense_matrix self);
void            TP_dense_matrix_print(TP_dense_matrix self, char *mess);
void            TP_dense_matrix_destroy(TP_dense_matrix slef);


TP_dense_2D TP_dense_2D_create(int n, int m);
TP_dense_2D TP_dense_2D_permute(TP_dense_2D A, int *row_perm, int *col_perm);
void        TP_dense_2D_facto(TP_dense_2D self);
void        TP_dense_2D_solve(TP_dense_2D self, TP_vector RHS);
void        TP_dense_2D_print(TP_dense_2D self, char *mess);
void        TP_dense_2D_destroy(TP_dense_2D self);

#endif // _TP_DENSE_H
