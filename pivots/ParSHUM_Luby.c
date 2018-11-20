#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include "ParSHUM_auxiliary.h"

#include "ParSHUM_Luby.h"
 
ParSHUM_Luby 
ParSHUM_Luby_create(ParSHUM_schur_matrix matrix)
{
  ParSHUM_Luby self = calloc(1, sizeof(*self));
  int n = matrix->n;

  self->n = n;
  
  self->score    = calloc((size_t) n, sizeof(*self->score));
  self->chosen   = calloc((size_t) n, sizeof(*self->chosen));
  self->position = calloc((size_t) n, sizeof(*self->position));
 
  return self;
}

void
ParSHUM_Luby_destroy(ParSHUM_Luby self)
{
  free(self->score);
  free(self->chosen);
  free(self->position);

  free(self);
}

long
ParSHUM_Luby_get_eligible(ParSHUM_schur_matrix matrix, ParSHUM_Luby Luby,
			  double value_tol, int *global_invr_col_perms,
			  int *global_invr_row_perms, int *cols,
			  int first_col, int last_col, int max_col_length)
{
  long best_marko = LONG_MAX;
  int i, j;
  int *chosen = Luby->chosen, *position = Luby->position;
  int yes = Luby->chosen_base + 1;
  
  for(i = first_col; i < last_col; i++)
    {
      int col = cols[i];
      CSC_struct *CSC = &matrix->CSC[col];

      if ( global_invr_col_perms[col] != ParSHUM_UNUSED_PIVOT ||
      	   CSC->nb_elem > max_col_length) {
      	continue;
      }

      int *rows    = CSC->row;
      int col_degree  = CSC->nb_elem - 1;
      int col_nb_elem = CSC->nb_numerical_eligible;
      int best_position  = -1;
      int col_best_row = INT_MAX;

      for ( j = 0; j < col_nb_elem; j++)
	{
	  int row = rows[j];
	  /* we have to test this because of singeltons */
	  if (global_invr_row_perms[row] != ParSHUM_UNUSED_PIVOT ) 
	    continue;
	  int row_degree = matrix->CSR[row].nb_elem - 1;
	  if ( col_best_row > row_degree)  {

	    col_best_row = row_degree;
	    best_position = j;
	  }
	}

      if (best_position != -1) {
	long col_best_marko =  (long) col_best_row * col_degree;
	if(col_best_marko < best_marko)
	  best_marko = col_best_marko;
	position[col] = best_position;
	chosen[col] = yes;
      }
    }
  
  return best_marko;
}

int
ParSHUM_Luby_assign_score(ParSHUM_Luby Luby, ParSHUM_schur_matrix matrix,
			  long allowed_marko, int *seed,
			  int *col_perm, int *row_perm, 
			  int *cols, int first_col, int last_col)
{
  int i, nb_candidates = 0;
  int *chosen = Luby->chosen, *positions = Luby->position;
  double *scores = Luby->score;
  int yes = Luby->chosen_base + 1, potential_pivot = Luby->chosen_base + 2;
  int my_seed = *seed;

  for ( i = first_col; i < last_col; i++)
    {
      int col = cols[i];
      CSC_struct *CSC = &matrix->CSC[col];

      if (chosen[col] != yes)
      	continue;
      int row  = CSC->row[positions[col]];
      long marko_cost = (matrix->CSR[row].nb_elem - 1) * (CSC->nb_elem - 1);
      if  (  marko_cost > allowed_marko)
      	continue;

      scores[col] = ((double) marko_cost / (double) allowed_marko)
      	          + (CSC->val[positions[col]] /CSC->col_max);

      chosen[col] = potential_pivot;
      
      col_perm[nb_candidates  ] = col;
      row_perm[nb_candidates++] = row;
    }
  
  *seed = my_seed;
  return nb_candidates;
}

void
ParSHUM_Luby_first_pass(ParSHUM_Luby Luby, ParSHUM_schur_matrix matrix,
			int *col_perm, int *row_perm, int nb_candidates)
{
  int i, j;
  int *chosen = Luby->chosen;
  double *scores = Luby->score;
  int potential_pivot = Luby->chosen_base + 2, discarded_pivot = Luby->chosen_base + 3;

  for( i = 0; i < nb_candidates; i++)
    {
      int piv_row = row_perm[i];
      int piv_col = col_perm[i];
      CSR_struct *CSR = &matrix->CSR[piv_row];
      int nb_elem = CSR->nb_elem;
      int *cols = CSR->col;
      double piv_score = scores[piv_col];
      
      for ( j = 0; j < nb_elem; j++) 
	{
	  int col = cols[j];
	  if ( chosen[col] == potential_pivot || chosen[col] == discarded_pivot  ) {
	    if ( piv_score > scores[col] ) { 
	      chosen[col] = discarded_pivot;
	    } else if ( piv_score < scores[col] )  {
	      chosen[piv_col] = discarded_pivot;
	    } else if  (piv_col > col)  {
	      chosen[piv_col] = discarded_pivot;
	    } else if  (piv_col < col)  {
	      chosen[col] = discarded_pivot;
	    }
	  }
	}
    }
}
 
int
ParSHUM_Luby_second_pass(ParSHUM_schur_matrix matrix, ParSHUM_Luby Luby, 
			 int *col_perm, int *row_perm, int nb_candidates)
{
  int i;
  int *chosen = Luby->chosen;
  int *positions = Luby->position;
  int discarded_pivot = Luby->chosen_base + 3;

  for( i = 0; i < nb_candidates; )
    {
      int col = col_perm[i];

      if (chosen[col] == discarded_pivot )  {
	col_perm[i] = col_perm[--nb_candidates];
	row_perm[i] = row_perm[  nb_candidates];
	col_perm[nb_candidates] = ParSHUM_UNUSED_PIVOT;
	row_perm[nb_candidates] = ParSHUM_UNUSED_PIVOT;
      } else  {
	int position1 = positions[col];
	CSC_struct *CSC = &matrix->CSC[col];
	double *val = CSC->val;
	int    *row = CSC->row;
	int position2 = CSC->nb_elem - 1;
	int tmp_row;
	double tmp_val;
	tmp_row = row[position1];
	tmp_val = val[position1];
	row[position1] = row[position2];
	val[position1] = val[position2];
	row[position2] = tmp_row;
	val[position2] = tmp_val;
	i++;
      }
    }

  return nb_candidates;
}

void 
ParSHUM_Luby_check_pivots(ParSHUM_Luby Luby, ParSHUM_schur_matrix matrix,
			  int *col_perms, int *row_perms, int nb_pivots)
{
  int i, j, n = Luby->n;
  int *col_counts = calloc((size_t) n, sizeof(*col_counts));
  int *row_counts = calloc((size_t) n, sizeof(*row_counts));
  char mess[2048];

  for ( i = 0; i < nb_pivots; i++)
    {
      CSC_struct *CSC = &matrix->CSC[col_perms[i]];
      int *rows = CSC->row; 
      int nb_elem = CSC->nb_elem;
      for ( j = 0; j < nb_elem; j++)
	row_counts[rows[j]]++;
      CSR_struct *CSR = &matrix->CSR[row_perms[i]];
      int *cols = CSR->col;
      nb_elem = CSR->nb_elem;
      for ( j = 0; j < nb_elem; j++)
	col_counts[cols[j]]++;
    }

  for ( i = 0; i < nb_pivots; i++)
    {
      if (col_counts[col_perms[i]] != 1) {
	snprintf(mess, 2048, "%d is suppose to be a pivot but row_count = %d\n",
		 col_perms[i], col_counts[col_perms[i]]);
	ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
      }
      if (row_counts[row_perms[i]] != 1)  {
	snprintf(mess, 2048, "%d is suppose to be a pivot but col_count = %d\n",
		 row_perms[i], row_counts[row_perms[i]]);
	ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
      }
    }

  free(col_counts);
  free(row_counts);
}
