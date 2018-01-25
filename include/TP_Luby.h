#ifndef __TP_LUBY_H
#define __TP_LUBY_H

#include "TP_schur_matrix.h"

typedef struct _TP_Luby_ *TP_Luby;


struct _TP_Luby_ {
  
  double *col_max_val;
  int    *col_max_row;
  double *row_max_val;
  int    *row_max_col;

  int *invr_col_perm;
  int *invr_row_perm;
  
  int m;
  int n;
};
 
TP_Luby  TP_Luby_create(TP_schur_matrix matrix);

void TP_Luby_destroy(TP_Luby self);

int TP_Luby_get_eligible(TP_schur_matrix matrix,
			 TP_Luby Luby, double value_tol,
			 int first_col, int last_col, 
			 int max_col_length);

void TP_Luby_get_candidates(TP_schur_matrix matrix, TP_Luby Luby,
			    int best_marko, int first_col, int last_col);

int TP_Luby_assign_score(TP_Luby Luby, TP_schur_matrix matrix,
			 int *col_perm, int *row_perm,
 			 int first_col, int last_col);
  
void TP_Luby_first_pass(TP_Luby Luby, TP_schur_matrix matrix,
			int *col_perm, int *row_perm, int nb_candidates);

int TP_Luby_second_pass(TP_schur_matrix matrix, TP_Luby Luby, 
			int *col_perm, int *row_perm, int nb_candidates);

/* void TP_Luby_get_max(TP_Luby Luby, TP_schur_matrix matrix, */
/* 		     int first_col, int last_col); */
  
/* int TP_Luby_first_pass(TP_Luby Luby, int *col_perm, int *row_perm); */

/* void TP_Luby_second_pass(TP_schur_matrix matrix, TP_Luby Luby,  */
/* 			int nb_eligeble_pivots, int *col_perm, int *row_perm); */

/* int TP_Luby_discard(TP_schur_matrix matrix, TP_Luby Luby,  */
/* 		    int nb_eligible_pivots, int *col_perm, int *row_perm); */

void TP_Luby_print_scores(TP_Luby Luby, char *mess);

void TP_Luby_check_pivots(TP_Luby Luby, TP_schur_matrix matrix,
			  int *col_perms, int *row_perm, int nb_pivots);


#endif // __TP_LUBY_H
