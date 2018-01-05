#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <limits.h>
#include "TP_auxiliary.h"

#include "TP_Luby.h"

double
randfrom(double min, double max) 
{
  double range = (max - min); 
  double div = RAND_MAX / range;
  return min + (rand() / div);
}
 
TP_Luby 
TP_Luby_create(int n, int m)
{
  TP_Luby self = calloc(1, sizeof(*self));

  self->n = n;
  self->m = m;
  
  self->scores = calloc( n, sizeof(*self->scores));
  
  self->col_max_val = malloc(n * sizeof(*self->col_max_val));
  self->col_max_row = malloc(n * sizeof(*self->col_max_row));
  self->row_max_val = malloc(m * sizeof(*self->row_max_val));
  self->row_max_col = malloc(m * sizeof(*self->row_max_col));

  self->invr_col_perm = malloc(n * sizeof(*self->invr_col_perm));
  self->invr_row_perm = malloc(n * sizeof(*self->invr_row_perm));
  srand(time(NULL));
  return self;
}

void 
TP_Luby_destroy(TP_Luby self)
{
  int i, n = self->n;

  for(i = 0; i < n; i++) 
    free(self->scores[i].score);

  free(self->col_max_val);
  free(self->col_max_row);
  free(self->row_max_val);
  free(self->row_max_col);

  free(self->invr_col_perm);
  free(self->invr_row_perm);

  free(self);
}

int
TP_Luby_get_candidates(TP_schur_matrix matrix,
		       TP_Luby Luby, double value_tol,
		       int first_col, int last_col, 
		       int max_col_length)
{
  int best_marko = INT_MAX;
  int i, j, unused = max_col_length;
  unused++;

  for(i = first_col; i < last_col; i++)
    {
      CSC_struct *CSC = &matrix->CSC[i];
      double *vals    = CSC->val;
      double col_max  = CSC->col_max, tmp_double;
      int *rows       = CSC->row;
      int col_nb_elem = CSC->nb_elem;
      TP_Luby_score *luby_score = &Luby->scores[i];
      int nb_eligible = 0, tmp_int;
      luby_score->nb_elem = 0;

      for ( j = 0; j < col_nb_elem; j++)
	{
	  int row = rows[j];
	  double val = vals[j];
	  int current_marko = (col_nb_elem - 1) * (matrix->CSR[row].nb_elem - 1);

	  if ( fabs(val) >= value_tol * col_max ) {
	    tmp_double = vals[nb_eligible];
	    tmp_int    = rows[nb_eligible];
	    vals[nb_eligible] = vals[j];
	    vals[j] = tmp_double;
	    rows[nb_eligible++] = rows[j];
	    rows[j] = tmp_int;
	  }
	  if ( best_marko > current_marko) 
	    best_marko = current_marko;
	}
      luby_score->nb_elem = nb_eligible;
    }

  return best_marko;
}

void
TP_Luby_assign_score(TP_schur_matrix matrix, TP_Luby Luby,
		     double marko_tol, int best_marko,
		     int first_col, int last_col)
{
  int i, j;

  for (i = first_col; i < last_col; i++)
    {
      CSC_struct *CSC = &matrix->CSC[i];
      TP_Luby_score *luby_score = &Luby->scores[i];
      double *vals    = CSC->val;
      int    *rows    = CSC->row;
      int col_nb_elem = CSC->nb_elem;
      int     nb_elem = luby_score->nb_elem;
      
      if (luby_score->allocated)  {
	luby_score->score = malloc(nb_elem * sizeof(*luby_score->score));
	luby_score->allocated = nb_elem;
      } else if (luby_score->allocated < nb_elem ) {
	luby_score->score = realloc(luby_score->score, nb_elem * sizeof(*luby_score->score));
	luby_score->allocated = nb_elem;
      }
      double *scores = luby_score->score;

      for (j = 0; j < nb_elem;) {
	int row    = rows[j];
	int current_marko = (col_nb_elem - 1) * (matrix->CSR[row].nb_elem - 1);
	double tmp_double;
	if ( current_marko > best_marko * marko_tol) {
	  nb_elem--;
	  tmp_double = vals[j];
	  vals[j] = vals[nb_elem];
	  vals[nb_elem]  = tmp_double;
	  rows[j] = rows[nb_elem];
	  rows[nb_elem] = row;
	} else { 
	  scores[j++] = randfrom(1.0, 15.0);
	}
      }
      luby_score->nb_elem = nb_elem;
    }
}

void
TP_Luby_get_max(TP_schur_matrix matrix, TP_Luby Luby,
		int first_col, int last_col)
{
  int i, j;
  double *col_max_val = Luby->col_max_val;
  int    *col_max_row = Luby->col_max_row;
  double *row_max_val = Luby->row_max_val;
  int    *row_max_col = Luby->row_max_col;

  bzero(Luby->col_max_val, Luby->n * sizeof(*Luby->col_max_val));
  bzero(Luby->row_max_val, Luby->m * sizeof(*Luby->row_max_val));
  int_array_memset(col_max_row, -1, Luby->n);
  int_array_memset(row_max_col, -1, Luby->m);

  int_array_memset(Luby->invr_col_perm, -1, Luby->n);
  int_array_memset(Luby->invr_row_perm, -1, Luby->m);

  for ( i = first_col; i < last_col; i++)
    {
      TP_Luby_score *luby_score = &Luby->scores[i];
      int nb_elem = luby_score->nb_elem;
      int *rows = matrix->CSC[i].row;
      double *scores = luby_score->score;
      for ( j = 0; j < nb_elem; j++)  {
	int row = rows[j];
	if (col_max_val[i] < scores[j]) {
	  col_max_val[i] = scores[j];
	  col_max_row[i] = row;
	}
	if (row_max_val[row] < scores[j]) {
	  row_max_val[row] = scores[j];
	  row_max_col[row] = i;
	}
      }
    }
}

int
TP_Luby_first_pass(TP_Luby Luby, int *col_perm, int *row_perm)
{
  int i, n = Luby->n, nb_pivots = 0;
  int *col_max     = Luby->col_max_row;
  double *col_max_val = Luby->col_max_val;
  int *row_max     = Luby->row_max_col;
  int *invr_row_perm = Luby->invr_row_perm;
  int *invr_col_perm = Luby->invr_col_perm;

  for( i = 0; i < n; i++)
    {
      int row = col_max[i];
      if ( col_max_val[i] >= 1.0 && row_max[row] == i ) {
	col_perm[nb_pivots] = i;
	row_perm[nb_pivots] = row;
	invr_col_perm[i]    = nb_pivots;
	invr_row_perm[row]  = nb_pivots++;
      }
    }

  return nb_pivots;
}

void
TP_Luby_second_pass(TP_schur_matrix matrix, TP_Luby Luby, 
		    int nb_eligible_pivots, int *col_perm, int *row_perm)
{
  int i, j;
  double *col_max_val = Luby->col_max_val;
  double *row_max_val = Luby->row_max_val;
  int *invr_row_perm  = Luby->invr_row_perm;
  int *invr_col_perm  = Luby->invr_col_perm;

  for( i = 0; i < nb_eligible_pivots; i++)
    {
      int col = col_perm[i];
      int row = row_perm[i];
      CSC_struct *CSC = &matrix->CSC[col];
      CSR_struct *CSR = &matrix->CSR[row];
      int nb_elem = CSC->nb_elem;
      int *rows = CSC->row;
      int *cols = CSR->col;
      double score = col_max_val[col];
      int eligble = 1;

      for ( j = 0; j < nb_elem; j++)
	{
	  int current_row = rows[j];
	  if (invr_row_perm[current_row] > -1 && score < row_max_val[current_row] )  {
	    col_perm[i] = -col - 1;
	    row_perm[i] = -row - 1;
	    eligble = 0;
	    break;
	  }
	}
      if (!eligble) 
	continue;
      nb_elem = CSR->nb_elem;
      for ( j = 0; j < nb_elem; j++)
	{
	  int current_col = cols[j];
	  if (invr_col_perm[current_col] > -1 && score < col_max_val[current_col] )  {
	    col_perm[i] = -col - 1;
	    row_perm[i] = -row - 1;
	    eligble = 0;
	    break;
	  }
	} 
      if (!eligble) 
	continue;
    }
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
	nb_eligible_pivots--;
	invr_col_perm[-(col_perm[i] + 1)] = -100;
	invr_row_perm[-(row_perm[i] + 1)] = -100;
	col_perm[i] = col_perm[nb_eligible_pivots];
	row_perm[i] = row_perm[nb_eligible_pivots];
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
  for(i = 0; i < nb_eligible_pivots; i++) 
    {
      int row     = row_perm[i];
      int *cols   = matrix->CSR[row].col;
      int nb_elem = matrix->CSR[row].nb_elem;

      for ( j = 0; j < nb_elem; j++)  
	Luby->scores[cols[j]].nb_elem = 0;
    }

  for(i = 0; i < nb_eligible_pivots; i++) 
    {
      int col     = col_perm[i];
      int *rows   = matrix->CSC[col].row;
      int nb_elem = matrix->CSC[col].nb_elem;

      for ( j = 0; j < nb_elem; j++)  
	invr_row_perm[rows[j]] = 1;
    }

  for( i = 0 ; i < n; i++) 
    {
      TP_Luby_score *luby_score = &Luby->scores[i];
      int nb_elem = luby_score->nb_elem;
      if (!nb_elem) 
	continue;
      CSC_struct *CSC = &matrix->CSC[i]; 
      int *rows = CSC->row;
      double *vals = CSC->val;
      double *scores = luby_score->score;

      for (j = 0; j < nb_elem; )
	{
	  int row = rows[j];
	  if ( invr_row_perm[row] > -1 )  {
	    nb_elem--;
	    double tmp_val = vals[j];
	    double tmp_score = scores[j]; 
	    vals[j] = vals[nb_elem];
	    vals[nb_elem] = tmp_val;
	    scores[j] = scores[nb_elem];
	    rows[j] = rows[nb_elem];
	    rows[nb_elem] = row;
	  } else {
	    j++;
	  }
	}
      luby_score->nb_elem = nb_elem;
    }

  return nb_eligible_pivots;
}

void
TP_Luby_print_scores(TP_Luby Luby, TP_schur_matrix matrix, char *mess)
{
  int n = Luby->n;
  
  printf("%s\n", mess);
  printf("Luby's scores by column\n");
  for(int i = 0; i < n; i++)
    {
      TP_Luby_score *luby_score = &Luby->scores[i];
      CSC_struct *CSC = &matrix->CSC[i];
      int *rows    = CSC->row;
      int nb_elem = luby_score->nb_elem;
      double *scores = luby_score->score;
      printf("================%d======================\n", i);
      for(int j = 0; j < nb_elem; j++)
	printf("%d:(%e)  ", rows[j], scores[j]);
      printf("\n");
    }

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
      if (col_counts[col_perms[i]] != 1) 
	printf("%d is suppose to be a pivot but row_count = %d\n",  col_perms[i], col_counts[col_perms[i]]);
      if (row_counts[row_perms[i]] != 1) 
	printf("%d is suppose to be a pivot but col_count = %d\n",  row_perms[i], row_counts[row_perms[i]]);
    }
} 
