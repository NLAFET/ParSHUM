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
  int brrr, n = matrix->n, m = matrix->m, i;

  brrr = time(NULL);

  self->n = n;
  self->m = m;
  
  self->col_max_val = malloc(n * sizeof(*self->col_max_val));
  self->col_max_row = malloc(n * sizeof(*self->col_max_row));
  self->row_max_val = malloc(m * sizeof(*self->row_max_val));
  self->row_max_col = malloc(m * sizeof(*self->row_max_col));

  self->invr_col_perm = malloc(n * sizeof(*self->invr_col_perm));
  self->invr_row_perm = malloc(n * sizeof(*self->invr_row_perm));

  /* srand(brrr); */ 
 
  return self;
}

void
TP_Luby_destroy(TP_Luby self)
{
  int i, n = self->n;

  free(self->col_max_val);
  free(self->col_max_row);
  free(self->row_max_val);
  free(self->row_max_col);

  free(self->invr_col_perm);
  free(self->invr_row_perm);

  free(self);
}

int
TP_Luby_get_eligible(TP_schur_matrix matrix,
		     TP_Luby Luby, double value_tol,
		     int first_col, int last_col, 
		     int max_col_length)
{
  int best_marko = INT_MAX;
  int i, j, unused = max_col_length;
  unused++;
  bzero(Luby->col_max_val, Luby->n * sizeof(*Luby->col_max_val));
  bzero(Luby->row_max_val, Luby->m * sizeof(*Luby->row_max_val));
  int_array_memset(Luby->col_max_row, -1, Luby->n);
  int_array_memset(Luby->row_max_col, -1, Luby->m);
  int_array_memset(Luby->invr_row_perm, -1, Luby->n);
  int_array_memset(Luby->invr_col_perm, -1, Luby->n);

  for(i = first_col; i < last_col; i++)
    {
      CSC_struct *CSC = &matrix->CSC[i];
      double *vals    = CSC->val;
      double col_max  = CSC->col_max;
      int *rows       = CSC->row;
      int col_degree  = CSC->nb_elem - 1;
      int col_nb_elem = CSC->nb_elem;
      int nb_eligible = 0;
      
      for ( j = 0; j < col_nb_elem; j++)
	{
	  int row = rows[j];
	  double val = vals[j];
	  int current_marko = col_degree * (matrix->CSR[row].nb_elem - 1);
	  if ( fabs(val) >= value_tol * col_max ) {
	    rows[j] = rows[nb_eligible];
	    vals[j] = vals[nb_eligible];
	    rows[nb_eligible  ] = row;
	    vals[nb_eligible++] = val;
	    if ( best_marko > current_marko) 
	      best_marko = current_marko;
	  }
	}
      CSC->nb_eligible = nb_eligible;
    }

  return best_marko;
}

void
TP_Luby_get_candidates(TP_schur_matrix matrix, TP_Luby Luby,
		       int best_marko, int first_col, int last_col)
{
  int i, j;

  for (i = first_col; i < last_col; i++)
    {
      CSC_struct *CSC = &matrix->CSC[i];
      double *vals    = CSC->val;
      int *rows       = CSC->row;
      int col_degree  = CSC->nb_elem - 1 ;
      int     nb_elem = CSC->nb_eligible;
      
      for (j = 0; j < nb_elem;) {
	int row = rows[j];
	int current_marko = col_degree  * (matrix->CSR[row].nb_elem - 1);
	if ( current_marko > best_marko) {
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
		     int *col_perm, int *row_perm,
		     int first_col, int last_col)
{
  int i, nb_candidates = 0;
  double *col_max_val = Luby->col_max_val;
  double *row_max_val = Luby->row_max_val;
  int    *invr_col_perm  = Luby->invr_col_perm;
  
  for ( i = first_col; i < last_col; i++)
    {
      CSC_struct *CSC = &matrix->CSC[i];
      int nb_elem     = CSC->nb_eligible;
      if (!nb_elem)
	continue;
      int row  = CSC->row[rand() % nb_elem];
      double tmp =  (double) rand() / (double) RAND_MAX ;

      col_max_val[i  ] = tmp;
      row_max_val[row] = tmp;

      invr_col_perm[i] = 1;
      
      col_perm[nb_candidates  ] = i;
      row_perm[nb_candidates++] = row;
    }

  return nb_candidates;
}

void
TP_Luby_first_pass(TP_Luby Luby, TP_schur_matrix matrix,
		   int *col_perm, int *row_perm, int nb_candidates)
{
  int i, n = Luby->n, nb_pivots = 0, j;
  double *col_max_val     = Luby->col_max_val;
  int    *row_max_col     = Luby->row_max_col;
  double *row_max_val     = Luby->row_max_val;
  int    *invr_col_perm   = Luby->invr_col_perm;

  for( i = 0; i < nb_candidates; i++)
    {
      int row = row_perm[i];
      CSR_struct *CSR = &matrix->CSR[row];
      int nb_elem = CSR->nb_elem;
      int *cols = CSR->col;
      int total = 0, killed = 0, my_col = col_perm[i];
      double score = col_max_val[my_col];

      for ( j = 0; j < nb_elem; j++) 
	{
	  int col = cols[j];
	  if ( invr_col_perm[col] == 1 ) {
	    total++;
	    if ( score > col_max_val[col] ) { 
	      invr_col_perm[col] = 2;
	      killed++;
	    }
	  }
	}
      if ( killed != total - 1 )
	invr_col_perm[my_col] = 2;
      
    }
}

int
TP_Luby_second_pass(TP_schur_matrix matrix, TP_Luby Luby, 
		    int *col_perm, int *row_perm, int nb_candidates)
{
  int i, j;
  double *col_max_val = Luby->col_max_val;
  double *row_max_val = Luby->row_max_val;
  int *invr_col_perm   = Luby->invr_col_perm;

  for( i = 0; i < nb_candidates; )
    {
      int col = col_perm[i];
      int row = row_perm[i];      

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

int
TP_Luby_discard(TP_schur_matrix matrix, TP_Luby Luby, 
		int nb_eligible_pivots, int *col_perm, int *row_perm)
{
  int i, j, n = Luby->n; 
  int *invr_row_perm  = Luby->invr_row_perm;
  int *invr_col_perm  = Luby->invr_col_perm;
  int old = nb_eligible_pivots;

  for(i = 0; i < nb_eligible_pivots; ) 
    {
      if ( col_perm[i] < 0) {
	invr_col_perm[-(col_perm[i] + 1)] = -100;
	invr_row_perm[-(row_perm[i] + 1)] = -100;
	col_perm[i] = col_perm[--nb_eligible_pivots];
	row_perm[i] = row_perm[  nb_eligible_pivots];
      } else {
	i++;
      }
    }

  for( i = nb_eligible_pivots; i < old; i++)
    {
      col_perm[i] = -100;
      row_perm[i] = -100;
    }    

  /* TODO: move this in the previous loop */
  /* for(i = 0; i < nb_eligible_pivots; i++)  */
  /*   { */
  /*     int row     = row_perm[i]; */
  /*     int *cols   = matrix->CSR[row].col; */
  /*     int nb_elem = matrix->CSR[row].nb_elem; */

  /*     for ( j = 0; j < nb_elem; j++)   */
  /* 	Luby->scores[cols[j]].nb_elem = 0; */
  /*   } */

  /* for(i = 0; i < nb_eligible_pivots; i++)  */
  /*   { */
  /*     int col     = col_perm[i]; */
  /*     int *rows   = matrix->CSC[col].row; */
  /*     int nb_elem = matrix->CSC[col].nb_elem; */

  /*     for ( j = 0; j < nb_elem; j++)   */
  /* 	invr_row_perm[rows[j]] = 1; */
  /*   } */

  /* for( i = 0 ; i < n; i++)  */
  /*   { */
  /*     TP_Luby_score *luby_score = &Luby->scores[i]; */
  /*     int nb_elem = luby_score->nb_elem; */
  /*     if (!nb_elem)  */
  /* 	continue; */
  /*     int *rows = luby_score->row; */
  /*     double *scores = luby_score->score; */

  /*     for (j = 0; j < nb_elem; ) */
  /* 	{ */
  /* 	  int row = rows[j]; */
  /* 	  if ( invr_row_perm[row] > -1 )  { */
  /* 	    scores[j] = scores[--nb_elem]; */
  /* 	    rows[j] = rows[nb_elem]; */
  /* 	  } else { */
  /* 	    j++; */
  /* 	  } */
  /* 	} */
  /*     luby_score->nb_elem = nb_elem; */
  /*   } */

  /* Luby->step++; */
  
  return nb_eligible_pivots;
}

void
TP_Luby_print_scores(TP_Luby Luby, char *mess)
{
  int n = Luby->n;
  
  /* printf("%s\n", mess); */
  /* printf("Luby's scores by column\n"); */
  /* for(int i = 0; i < n; i++) */
  /*   { */
  /*     TP_Luby_score *luby_score = &Luby->scores[i]; */
  /*     int *rows    = luby_score->row; */
  /*     double *scores = luby_score->score; */
  /*     int nb_elem = luby_score->nb_elem; */
  /*     printf("================%d======================\n", i); */
  /*     for(int j = 0; j < nb_elem; j++) */
  /* 	printf("%d:(%e)  ", rows[j], scores[j]); */
  /*     printf("\n"); */
  /*   } */

}

void 
TP_Luby_check_pivots(TP_Luby Luby, TP_schur_matrix matrix,
		     int *col_perms, int *row_perms, int nb_pivots)
{
  int i, j, n = Luby->n;
  int *col_counts = calloc(n, sizeof(*col_counts));
  int *row_counts = calloc(n, sizeof(*row_counts));

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
