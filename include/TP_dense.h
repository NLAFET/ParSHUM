#ifndef _TP_DENSE_H
#define _TP_DENSE_H

typedef  struct _TP_vector *TP_vector;
typedef  struct _TP_dense_matrix *TP_dense_matrix;

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

TP_vector TP_vector_create(int n);
void      TP_vector_permute(TP_vector self, int *perms);
void      TP_vector_memset(TP_vector self, double val);
void      TP_vector_print(TP_vector self, char *mess);
void      TP_vector_destroy(TP_vector self);

TP_dense_matrix TP_dense_matrix_create(int n, int m);
void            TP_dense_matrix_factorize(TP_dense_matrix self);
void            TP_dense_matrix_print(TP_dense_matrix self, char *mess);
void            TP_dense_matrix_destroy(TP_dense_matrix slef);

#endif // _TP_DENSE_H
