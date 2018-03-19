#ifndef __TP_LUBY_H
#define __TP_LUBY_H

#include "TP_schur_matrix.h"

typedef struct _TP_Luby_ *TP_Luby;


struct _TP_Luby_ {
  
  double *score;
  int    *chosen;

  int *position;
  
  int chosen_base;
  int n;
};
 
TP_Luby  TP_Luby_create(TP_schur_matrix matrix);

void TP_Luby_destroy(TP_Luby self);

int TP_Luby_get_eligible(TP_schur_matrix matrix, TP_Luby Luby,
 			 double value_tol, int *global_invr_col_perms,
			 int *global_invr_row_perms, int *cols, int first_col,
			 int last_col, int max_col_length);

void TP_Luby_get_candidates(TP_schur_matrix matrix, TP_Luby Luby,
			    int allowed_marko, int *cols,
			    int first_col, int last_col);

int TP_Luby_assign_score(TP_Luby Luby, TP_schur_matrix matrix,
			 int allowed_marko, int *seed,
			 int *col_perm, int *row_perm, 
			 int *cols, int first_col, int last_col);

void TP_Luby_first_pass(TP_Luby Luby, TP_schur_matrix matrix,
			int *col_perm, int *row_perm, int nb_candidates);

int TP_Luby_second_pass(TP_schur_matrix matrix, TP_Luby Luby, 
			int *col_perm, int *row_perm, int nb_candidates);

void TP_Luby_check_pivots(TP_Luby Luby, TP_schur_matrix matrix,
			  int *col_perms, int *row_perm, int nb_pivots);


#endif // __TP_LUBY_H
