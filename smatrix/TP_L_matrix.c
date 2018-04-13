#include <stdio.h>
#include <stdlib.h>
#include "TP_L_matrix.h"



TP_L_matrix 
TP_L_matrix_create(int n)
{
  TP_L_matrix self = malloc(sizeof(*self));

  self->n = 0;
  self->nnz = 0;
  self->col = calloc((size_t) n, sizeof(*self->col));
  
  return self;
}

void
TP_L_matrix_solve(TP_L_matrix L, TP_vector RHS, int *perms)
{
  int    col, n = L->n, i;
  double *rhs_val = RHS->vect;

  for( col = 0; col < n; col++)
    {
      CSC_struct *CSC = &L->col[col];
      double *L_val = CSC->val;
      int    *L_row = CSC->row;
      int nb_elem = CSC->nb_elem;
      for( i = 0; i < nb_elem; i++)
  	rhs_val[perms[L_row[i]]] -= L_val[i] * rhs_val[col];
    }
}

void
TP_L_matrix_print(TP_L_matrix self, char *mess)
{
  int n = self->n, i, j;

  printf("%s\n", mess);
  for (i = 0; i < n; i++) 
    {
      CSC_struct *CSC = &self->col[i];
      double *L_val = CSC->val;
      int    *L_row = CSC->row;
      int nb_elem = CSC->nb_elem;
      printf("================%d======================\n ", i);
      for (j = 0; j < nb_elem; j++)
	printf("%d(%e)  ", L_row[j], L_val[j]);
      printf("\n");
    }
  printf("\n");
}

void
TP_L_matrix_destroy(TP_L_matrix self)
{
  free(self->col);
  free(self);
}
