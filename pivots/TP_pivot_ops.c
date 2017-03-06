#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h> 

#include "TP_auxiliary.h"

#include "TP_pivot_list.h"

struct pivot_candidates
{
  int best_markov;
  int nb_found;
  int *row;
  int *markov;
};


struct pivot_candidates *
get_cadidates(TP_schur_matrix matrix, double value_tol,
	      int first_col, int last_col)
{
  struct pivot_candidates *self;

  int best_markov = INT_MAX;
  int i, j;
  
  self           = malloc((size_t) sizeof(*self));
  self->row      = malloc( (last_col - first_col) * sizeof(*self->row));
  self->markov   = malloc( (last_col - first_col) * sizeof(*self->markov));
  self->nb_found = 0;
 
  for(i = first_col; i < last_col; i++)
    self->row[i] = -1;
  for(i = first_col; i < last_col; i++)
    self->markov[i] = INT_MAX;

  for(i = first_col; i < last_col; i++)
    {
      struct CSC_struct *CSC = &matrix->CSC[i];
      int *rows       = CSC->row;
      double *vals    = CSC->val;
      int col_nb_elem = CSC->nb_elem;
      double col_max  = CSC->col_max;

      for(j = 0; j < col_nb_elem; j++)
  	{
	  int row    = rows[j];
	  double val = vals[j]; 
  	  int current_markov = col_nb_elem * matrix->CSR[row].nb_elem;
  	  if ( fabs(val) >= value_tol * col_max &&
  	       current_markov < self->markov[i] ) {
  	    self->row[i]  = row;
  	    self->markov[i]  = current_markov;
	    self->nb_found++;
  	  }
  	}
      best_markov = ( self->markov[i] < best_markov  ) ? self->markov[i] : best_markov;
    }

  self->best_markov = best_markov;
  return self;
}


TP_pivot_list
get_possible_pivots(TP_schur_matrix matrix, int *random_col, 
		    double value_tol, double markov_tol,
		    int nb_init_block)
{
  TP_pivot_list self;
  int i, elems_per_block, elems_per_set = 0, nb_blocks = 1;
  int n = matrix->n;
  struct pivot_candidates *candidates;

  self = TP_pivot_list_create();
  
  candidates = get_cadidates(matrix, value_tol, 0, n);
  if ( nb_init_block  > candidates->nb_found ) 
    nb_init_block = candidates->nb_found;

  elems_per_block = candidates->nb_found / nb_init_block;
  
  for(i = 0; i < n; i++) {
    int col = random_col[i];
    if (elems_per_set > elems_per_block && nb_blocks < nb_init_block) {
      elems_per_set = 0;
      nb_blocks++;
    }

    if (candidates->row[col] == -1 ) 
      continue;

    if ( matrix->CSC[col].nb_elem * matrix->CSR[candidates->row[col]].nb_elem >
	 markov_tol * candidates->best_markov) 
      continue;
    /* else */
    /*   TP_pivot_list_insert_new_set(self, matrix, candidates->row[i], i, candidates->markov[i]); */
    
    if ( elems_per_set == 0) {
      TP_pivot_list_insert_new_set(self, matrix, candidates->row[col], col, candidates->markov[col]);
    } else {
      TP_pivot_set set = self->last;
      if( !set->cols_count[col] && !set->rows_count[candidates->row[col]] ) {
    	TP_pivot_cell cell = TP_pivot_cell_create(candidates->row[col], col, candidates->markov[col]);
    	add_cell_to_sorted_set(set, cell, matrix);
      }
    }
    elems_per_set++;
  }
    
  free(candidates->row);
  free(candidates->markov);
  free(candidates);

  return self;
}




TP_pivot_list 
merge_pivot_sets(TP_pivot_list self, TP_schur_matrix matrix)
{
  while (self->nb_elem > 1)
    {
      TP_pivot_list merged_list = TP_pivot_list_create();
      int i;
      int nb_merges = (self->nb_elem + 1) / 2;

      for( i = 0; i < nb_merges; i++)
	{
	  TP_pivot_set merged_set;
	  
	  merged_set = get_next_merging_set(self);
	  merged_set = merge_to_larger_set(merged_set, matrix);

/* #ifdef MERGING_INTO_ONE_SET */
	  /* merged_set = merge_sorted_sets(merged_set); */
	  /* merged_set = get_independent_pivots(merged_set, matrix); */
/* #endif */
	  TP_pivot_list_insert_set(merged_list, merged_set);
	}

      TP_pivot_list_destroy(self);
      self = merged_list;
  }

  return self;
}


void
add_cell_to_sorted_set(TP_pivot_set set, TP_pivot_cell cell, TP_schur_matrix matrix)
{
  TP_pivot_cell merged = set->cells;
  
  if (set->nb_elem == 0) 
    set->cells = cell;
  else if (merged->marko > cell->marko) {
    cell->next = set->cells;
    set->cells = cell;
  } else { 
    while(merged->next)
      {
	if ( merged->next->marko > cell->marko ) 
	  break;
	else
	  merged = merged->next;
      }
    TP_pivot_cell tmp = merged->next;
    merged->next = cell;
    cell->next   = tmp;
  }
  set->nb_elem++;
  update_counter(set->cols_count, matrix->CSR[cell->row].col, matrix->CSR[cell->row].nb_elem);
  update_counter(set->rows_count, matrix->CSC[cell->col].row, matrix->CSC[cell->col].nb_elem);
}


TP_pivot_set
merge_to_larger_set(TP_pivot_set self, TP_schur_matrix matrix)
{
  TP_pivot_set large, small;
  TP_pivot_cell cells_to_merge;
  int *rows_count, *cols_count;

  if (!self || !self->next )
    return self;

  if (self->next->next) 
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"there are more than two sets to merge, the sets following the second set will be ignored");
    
  if (self->nb_elem >= self->next->nb_elem) {
    large = self;
    small = self->next;
  } else {
    large = self->next;
    small = self;
  }

  cells_to_merge = small->cells;
  rows_count     = large->rows_count;
  cols_count     = large->cols_count;
  small->cells   = NULL;

  while(cells_to_merge)
    {
      if( !cols_count[cells_to_merge->col] && !rows_count[cells_to_merge->row] ) {
	TP_pivot_cell tmp = cells_to_merge;
	cells_to_merge = cells_to_merge->next;
	tmp->next = NULL;
	add_cell_to_sorted_set(large, tmp, matrix);
      }  else { 
	TP_pivot_cell tmp = cells_to_merge;
	cells_to_merge = cells_to_merge->next;
	free(tmp);
      }
    }

  TP_pivot_set_destroy(small);
  self = large;
  large->next = NULL;
  return large;
}
