#include <stdio.h>
#include <stdlib.h>
#include "TP_L_matrix.h"



TP_L_matrix 
TP_L_matrix_create(int n, long nnz)
{
  TP_L_matrix self = malloc(sizeof(*self));
  self->start = calloc((size_t) n + 1,   sizeof(*self->start));
  self->end   = calloc((size_t) n,   sizeof(*self->end));
  self->val   = calloc((size_t) nnz, sizeof(*self->val));
  self->row   = calloc((size_t) nnz, sizeof(*self->row));
  self->n = 0;
  self->nnz = 0;
  self->allocated = nnz;
  
  return self;
}


void
TP_L_matrix_realloc(TP_L_matrix self) 
{
  self->allocated *= 2;
  self->val = realloc(self->val, (size_t) self->allocated * sizeof(*self->val));
  self->row = realloc(self->row, (size_t) self->allocated * sizeof(*self->row));
}

void
TP_L_matrix_solve(TP_L_matrix L, TP_vector RHS, int *perms)
{
  int    col, n = L->n, i;
  double *L_val = L->val, *rhs_val = RHS->vect;
  int    *L_row     = L->row;
  long *start = L->start, *end = L->end;

  for( col = 0; col < n; col++)
    {
      long col_start = start[col], col_end = end[col];
      for( i = col_start; i < col_end; i++)
  	rhs_val[perms[L_row[i]]] -= L_val[i] * rhs_val[col];
    }
}

void
TP_L_matrix_print(TP_L_matrix self, char *mess)
{
  int n = self->n, i;
  long j;
  double *vals = self->val;
  int    *rows = self->row;
  long *start = self->start, *end = self->end;

  printf("%s\n", mess);
  for (i = 0; i < n; i++) {
    printf("================%d======================\n ", i);
    for (j = start[i]; j < end[i]; j++)
      printf("%d(%e)  ", rows[j], vals[j]);
    printf("\n");
  }
  printf("\n");
}


void
TP_L_matrix_destroy(TP_L_matrix self)
{
  free(self->start);
  free(self->end);
  free(self->val);
  free(self->row);
  free(self);
}
