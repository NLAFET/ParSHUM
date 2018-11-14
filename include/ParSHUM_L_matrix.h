#ifndef _ParSHUM_L_MATRIX
#define _ParSHUM_L_MATRIX

#include "ParSHUM_schur_matrix.h"
#include "ParSHUM_dense.h"


ParSHUM_L_matrix  ParSHUM_L_matrix_create(int n);
void              ParSHUM_L_matrix_solve(ParSHUM_L_matrix L, ParSHUM_vector RHS, int *perms);
void              ParSHUM_L_matrix_print(ParSHUM_L_matrix self, char *mess);
void              ParSHUM_L_matrix_destroy(ParSHUM_L_matrix self);

#endif //_ParSHUM_L_MATRIX
