#ifndef __ParSHUM_LUBY_H
#define __ParSHUM_LUBY_H

#include "ParSHUM_schur_matrix.h"

typedef struct _ParSHUM_Luby_ *ParSHUM_Luby;


struct _ParSHUM_Luby_ {
  
  double *score;
  int    *chosen;

  int *position;
  
  int chosen_base;
  int n;
};
 
ParSHUM_Luby  ParSHUM_Luby_create(ParSHUM_schur_matrix matrix);

void ParSHUM_Luby_destroy(ParSHUM_Luby self);

long ParSHUM_Luby_get_eligible(ParSHUM_schur_matrix matrix, ParSHUM_Luby Luby,
			       double value_tol, int *global_invr_col_perms,
			       int *global_invr_row_perms, int *cols, int first_col,
			       int last_col, int max_col_length);

void ParSHUM_Luby_get_candidates(ParSHUM_schur_matrix matrix, ParSHUM_Luby Luby,
				 int allowed_marko, int *cols,
				 int first_col, int last_col);

int ParSHUM_Luby_assign_score(ParSHUM_Luby Luby, ParSHUM_schur_matrix matrix,
			      long allowed_marko, int *seed,
			      int *col_perm, int *row_perm, 
			      int *cols, int first_col, int last_col);

void ParSHUM_Luby_first_pass(ParSHUM_Luby Luby, ParSHUM_schur_matrix matrix,
			     int *col_perm, int *row_perm, int nb_candidates);

int ParSHUM_Luby_second_pass(ParSHUM_schur_matrix matrix, ParSHUM_Luby Luby, 
			     int *col_perm, int *row_perm, int nb_candidates);

void ParSHUM_Luby_check_pivots(ParSHUM_Luby Luby, ParSHUM_schur_matrix matrix,
			       int *col_perms, int *row_perm, int nb_pivots);


#endif // __ParSHUM_LUBY_H
