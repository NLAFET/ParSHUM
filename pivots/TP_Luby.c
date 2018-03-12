#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include "TP_auxiliary.h"

#include "TP_Luby.h"
 
TP_Luby 
TP_Luby_create(TP_schur_matrix matrix)
{
  TP_Luby self = calloc(1, sizeof(*self));
  int n = matrix->n, m = matrix->m;

  self->n = n;
  self->m = m;
  
  self->col_max_val = malloc((size_t) n * sizeof(*self->col_max_val));
  self->col_max_row = malloc((size_t) n * sizeof(*self->col_max_row));
  self->row_max_val = malloc((size_t) m * sizeof(*self->row_max_val));
  self->row_max_col = malloc((size_t) m * sizeof(*self->row_max_col));

  self->invr_col_perm = malloc((size_t) n * sizeof(*self->invr_col_perm));
  self->invr_row_perm = malloc((size_t) n * sizeof(*self->invr_row_perm));

  /* srand(brrr); */ 
 
  return self;
}

void
TP_Luby_destroy(TP_Luby self)
{
  free(self->col_max_val);
  free(self->col_max_row);
  free(self->row_max_val);
  free(self->row_max_col);

  free(self->invr_col_perm);
  free(self->invr_row_perm);

  free(self);
}

int
TP_Luby_get_eligible(TP_schur_matrix matrix, TP_Luby Luby,
		     double value_tol, int *global_invr_col_perms,
		     int *cols, int first_col, int last_col, 
		     int max_col_length)
{
  int best_marko = INT_MAX;
  int i, j, unused = max_col_length;
  unused++;

  /* TODO: get rid of all init */
  bzero(Luby->col_max_val, (size_t) Luby->n * sizeof(*Luby->col_max_val));
  int_array_memset(Luby->invr_col_perm, -1, Luby->n);

  for(i = first_col; i < last_col; i++)
    {
      int col = cols[i];
      CSC_struct *CSC = &matrix->CSC[col];
      if ( global_invr_col_perms[col] != TP_UNUSED_PIVOT) {
	CSC->nb_numerical_eligible = 0;
	continue;
      }
      int  *rows    = CSC->row;
      int col_degree  = CSC->nb_elem - 1;
      int col_nb_elem = CSC->nb_numerical_eligible;
      int nb_eligible = 0;
      
      for ( j = 0; j < col_nb_elem; j++)
	{
	  int row = rows[j];
	  if (!matrix->CSR[row].nb_elem)  {
	    printf("KO empty line\n");
	    continue;
	  }
	  int current_marko = col_degree * (matrix->CSR[row].nb_elem - 1);
	  if ( best_marko > current_marko)
	    best_marko = current_marko;
	}
    }

  return best_marko;
}

void
TP_Luby_get_candidates(TP_schur_matrix matrix, TP_Luby Luby,
		       int allowed_marko, int *cols,
		       int first_col, int last_col)
{
  int i, j;

  for (i = first_col; i < last_col; i++)
    {
      CSC_struct *CSC = &matrix->CSC[cols[i]];
      double *vals    = CSC->val;
      int *rows       = CSC->row;
      int col_degree  = CSC->nb_elem - 1 ;
      int     nb_elem = CSC->nb_numerical_eligible;
      /* 
	 TODO: this loop is not needed, the potential pivot could be 
	 picked now randomly or chose the entry with best markowitz cost,
	 or a combination of the strategies.
      */

      for (j = 0; j < nb_elem;) {
	int row = rows[j];
	if (!matrix->CSR[row].nb_elem) {
	  printf("KO in get candidates with the row\n");
	  continue;
	}
	int current_marko = col_degree  * (matrix->CSR[row].nb_elem - 1);
	if ( current_marko > allowed_marko) {
	  double tmp = vals[j];
	  rows[j] = rows[--nb_elem];
	  rows[nb_elem] = row;
	  vals[j] = vals[  nb_elem];
	  vals[nb_elem] = tmp;
	} else {
	  j++;
	}
      }

      CSC->nb_eligible = nb_elem;
    }
}

int
TP_Luby_assign_score(TP_Luby Luby, TP_schur_matrix matrix,
		     int *global_row_perms,
		     int *seed, int *col_perm, int *row_perm, 
		     int *cols, int first_col, int last_col)
{
  int i, nb_candidates = 0;
  double *col_max_val = Luby->col_max_val;
  int    *invr_col_perm  = Luby->invr_col_perm;
  int my_seed = *seed;

  for ( i = first_col; i < last_col; i++)
    {
      int col = cols[i];
      CSC_struct *CSC = &matrix->CSC[col];
      int nb_elem     = CSC->nb_eligible;
      if (!nb_elem)
	continue;
      int row  = CSC->row[abs(TP_rand_int(&my_seed, nb_elem))];
      if (global_row_perms[row] != TP_UNUSED_PIVOT) 
	continue;
      double score =  TP_rand_double(&my_seed);

      col_max_val[col] = score;
      invr_col_perm[col] = 1;
      
      col_perm[nb_candidates  ] = col;
      row_perm[nb_candidates++] = row;
    }
  
  *seed = my_seed;
  return nb_candidates;
}

void
TP_Luby_first_pass(TP_Luby Luby, TP_schur_matrix matrix,
		   int *col_perm, int *row_perm, int nb_candidates)
{
  int i, j;
  /* TODO: rename col_max_val to score */
  double *col_max_val     = Luby->col_max_val;
  int    *invr_col_perm   = Luby->invr_col_perm;

  for( i = 0; i < nb_candidates; i++)
    {
      int piv_row = row_perm[i];
      int piv_col = col_perm[i];
      CSR_struct *CSR = &matrix->CSR[piv_row];
      int nb_elem = CSR->nb_elem;
      int *cols = CSR->col;
      double piv_score = col_max_val[piv_col];
      
      for ( j = 0; j < nb_elem; j++) 
	{
	  int col = cols[j];
	  if ( invr_col_perm[col] == 1 || invr_col_perm[col] == 2  ) {
	    if ( piv_score > col_max_val[col] ) { 
	      invr_col_perm[col] = 2;
	    } else if ( piv_score < col_max_val[col] )  {
	      invr_col_perm[piv_col] = 2;
	    } else if  (piv_col > col)  {
	      invr_col_perm[piv_col] = 2;
	    } else if  (piv_col < col)  {
	      invr_col_perm[col] = 2;
	    }
	  }
	}
    }
}

int
TP_Luby_second_pass(TP_schur_matrix matrix, TP_Luby Luby, 
		    int *col_perm, int *row_perm, int nb_candidates)
{
  int i;
  int *invr_col_perm   = Luby->invr_col_perm;

  for( i = 0; i < nb_candidates; )
    {
      int col = col_perm[i];

      if (invr_col_perm[col] == 2 )  {
	invr_col_perm[col] = -100;
	col_perm[i] = col_perm[--nb_candidates];
	row_perm[i] = row_perm[  nb_candidates];
	col_perm[nb_candidates] = TP_UNUSED_PIVOT;
	row_perm[nb_candidates] = TP_UNUSED_PIVOT;
      } else  {
	i++;
      }
    }
  return nb_candidates;
}

void 
TP_Luby_check_pivots(TP_Luby Luby, TP_schur_matrix matrix,
		     int *col_perms, int *row_perms, int nb_pivots)
{
  int i, j, n = Luby->n;
  int *col_counts = calloc((size_t) n, sizeof(*col_counts));
  int *row_counts = calloc((size_t) n, sizeof(*row_counts));

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
	printf("%d is suppose to be a pivot but row_count = %d\n",  col_perms[i], col_counts[col_perms[i]]);
	GDB_BREAK;
      }
      if (row_counts[row_perms[i]] != 1)  {
	printf("%d is suppose to be a pivot but col_count = %d\n",  row_perms[i], row_counts[row_perms[i]]);
	GDB_BREAK;
      }
    }

  free(col_counts);
  free(row_counts);
}
