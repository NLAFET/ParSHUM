#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h> 
#include <omp.h>

#include "TP_auxiliary.h"

#include "TP_pivot_list.h"


void
get_candidates(TP_schur_matrix matrix,
	       TP_pivot_candidates candidates,
	       double value_tol,
	       int first_col, int last_col, 
	       int max_col_length)
{
  int best_marko = INT_MAX;
  int i, j;
  int me =  omp_get_thread_num();

  for(i = first_col; i < last_col; i++)
    candidates->row[i] = -1;
  for(i = first_col; i < last_col; i++)
    candidates->marko[i] = INT_MAX;

  for(i = first_col; i < last_col; i++)
    {
      CSC_struct *CSC  = &matrix->CSC[i];
      int *rows       = CSC->row;
      double *vals    = CSC->val;
      int col_nb_elem = CSC->nb_elem;
      double col_max  = CSC->col_max;
      if ( col_nb_elem > max_col_length)
	continue;

      for(j = 0; j < col_nb_elem; j++)
  	{
	  int row    = rows[j];
	  double val = vals[j]; 
  	  int current_marko = (col_nb_elem - 1) * (matrix->CSR[row].nb_elem - 1);
  	  if ( fabs(val) >= value_tol * col_max &&
  	       current_marko < candidates->marko[i] ) {
  	    candidates->row[i]  = row;
  	    candidates->marko[i]  = current_marko;
  	  }
  	}
      best_marko = ( candidates->marko[i] < best_marko  ) ? candidates->marko[i] : best_marko;
    }

  candidates->best_marko[me] = best_marko;
}


void
create_init_set(TP_solver solver, 
		TP_schur_matrix matrix,
		TP_pivot_list list,
		int *random_col, 
		TP_pivot_candidates candidates,
		int max_marko, int start, int end)
{
  int i;
  TP_pivot_set set = TP_pivot_set_create(solver, matrix->n, matrix->m);
  
  for( i = start; i < end; i++) 
    {
      int col = random_col[i];

      if ( candidates->marko[col]  <=  max_marko) {
	TP_pivot_cell cell = TP_pivot_cell_create(candidates->row[col], col, candidates->marko[col]);
	if ( (add_cell_to_sorted_set(set, cell, matrix)) ) 
	  TP_pivot_cell_destroy(cell);
      }
    }

  if( set->nb_elem) 
#pragma omp critical
    {
      TP_pivot_list_insert_set(list, set);
    }
  else 
    TP_pivot_set_destroy(set, solver);
}

TP_pivot_list
get_possible_pivots(TP_solver solver, TP_schur_matrix matrix, int *random_col, 
                    TP_pivot_candidates candidates, int nb_threads,
		    double value_tol, double marko_tol,
		    int nb_candidates_per_block)
{
  TP_pivot_list self = TP_pivot_list_create();
  int i, candidates_per_set = 0;
  int n = matrix->n, best_marko;
  int max_col_length = (int)  ( matrix->nnz / (matrix->n - solver->done_pivots) );
  max_col_length /= 1.2;

#pragma omp parallel num_threads(nb_threads)
  {
    int me =  omp_get_thread_num();
    int start = me * ( n / nb_threads);
    int  end  = (me+1) * ( n / nb_threads);
    if ( me == (nb_threads - 1) )
      end = n;
    
    get_candidates(matrix, candidates, value_tol, start, end, max_col_length);
  }  

  best_marko = candidates->best_marko[0];
  for(i = 1; i < nb_threads; i++) 
    best_marko = (best_marko > candidates->best_marko[i]) ? candidates->best_marko[i] : best_marko;
  best_marko = (best_marko == 0) ? 1 : best_marko;

#pragma omp parallel num_threads(nb_threads)
  {
#pragma omp single
    {
      int last_step;
      for(i = 0, last_step = 0; i < n; i++) {
	int col = random_col[i];
	
	if (candidates->row[col] == -1 ) 
	  continue;
	
	candidates_per_set++;
	
	if (candidates_per_set > nb_candidates_per_block) {
	  candidates_per_set = 0;
#pragma omp task shared (self, random_col, candidates, solver) 
	  {
	    int start = last_step, end = i;
	    int markov_tolerance = (int) (best_marko * marko_tol);
	    /* best_marko * marko_tol could overbuff, if marko tol is too large  */
	    create_init_set(solver, matrix, self, random_col, candidates, markov_tolerance, start, end);
	  }
	  last_step = i;
	}
      }
      if (candidates_per_set)
#pragma omp task shared (self, random_col, candidates)
	{
	  create_init_set(solver, matrix, self, random_col, candidates, (int) (best_marko * marko_tol), last_step, i);
	}
#pragma omp taskwait
    } //single
  } //parallel
  
  return self;
}




TP_pivot_list 
merge_pivot_sets(TP_pivot_list self, TP_schur_matrix matrix, TP_solver solver)
{
  TP_pivot_list merged_list;
#pragma omp parallel shared(self, merged_list) num_threads(solver->exe_parms->nb_threads)
  {
#pragma omp single
    {
      while (self->nb_elem > 1)
	{
	  merged_list = TP_pivot_list_create();
	  int i;
	  int nb_merges = (self->nb_elem + 1) / 2;
	  
	  for( i = 0; i < nb_merges; i++)
	    {
#pragma omp task shared(self, merged_list)
              {
		TP_pivot_set merged_set;
		
#pragma omp critical
		{
		  merged_set = get_next_merging_set(self);
		}
		merged_set = merge_to_larger_set(merged_set, matrix, solver);
		
#pragma omp critical
		{
		  TP_pivot_list_insert_set(merged_list, merged_set);
		}
	      }// task
	    }// for 
#pragma omp taskwait
	  TP_pivot_list_destroy(self, solver);
	  self = merged_list;
	}// while
    } // single
  } // parallel

  return self;
}


TP_pivot_cell
add_cell_to_sorted_set(TP_pivot_set set, TP_pivot_cell cell, TP_schur_matrix matrix)
{
  TP_pivot_cell merged = set->cells;
  int base = set->base;
  if( set->cols_count[cell->col] >= base ||
      set->rows_count[cell->row] >= base )
    return cell;
  
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

  update_counter(set->cols_count, matrix->CSR[cell->row].col, matrix->CSR[cell->row].nb_elem, set->base);
  update_counter(set->rows_count, matrix->CSC[cell->col].row, matrix->CSC[cell->col].nb_elem, set->base);
  set->nb_elem++;

  return NULL;
}

TP_pivot_cell
get_independent_cells(TP_pivot_cell cell, int *cols_count, int *rows_count, int base)
{
  TP_pivot_cell first = NULL, last;

  while (cell)
    {
      TP_pivot_cell tmp = cell->next;
      if( cols_count[cell->col] >= base ||
	  rows_count[cell->row] >= base )
	{
	  free(cell);
	}
      else 
	{ 
	  if ( !first )  {
	    first = cell;
	    last = cell;
	  } else {
	    last->next = cell;
	    last = cell;
	  }
	  last->next = NULL;
	}
      cell = tmp;
    }
  return first;
}

/* TP_pivot_cell */
/* TP_cells_merge(TP_pivot_cell first, TP_pivot_cell second) */
/* { */
/*   if (!first && !second)  */
/*     return NULL; */
/*   if (!first)  */
/*     return second; */
/*   if(!second) */
/*     return first; */

/*   while ( second )  */
/*     { */
/*       TP_pivot_cell cell = second; */
/*       TP_pivot_cell cell1, cell2; */
/*       second = second->next; */

/*       if ( first->marko > cell->marko) { */
/* 	cell->next = first; */
/* 	first = cell; */
/* 	continue; */
/*       } */
/*       cell1 = first; */
/*       cell2 = first->next; */
/*       while ( cell2 ) */
/* 	{ */
/* 	  if( cell2->marko > cell->marko ) */
/* 	    break; */
/* 	  cell1 = cell2; */
/* 	  cell2 = cell2->next; */
/* 	}  */
/*       cell1->next = cell; */
/*       cell->next = cell2; */
/*     } */

/*   return first; */
/* } */

TP_pivot_cell
TP_cells_merge(TP_pivot_cell first, TP_pivot_cell second)
{
  TP_pivot_cell res, cells;
  if (!first && !second) 
    return NULL;
  if (!first) 
    return second;
  if (!second) 
    return first;
  if (first->marko <= second->marko) {
    res = cells = first;
    first = first->next;
  } else {
    res = cells = second;
    second = second->next;
  }
  cells->next = NULL;
  
  while (first && second)
    {
      if (first->marko <= second->marko) {
	cells->next = first; 
	while (first->next && first->marko <= second->marko) 
	  first = first->next;
	cells = first;
	first = first->next;
      } else {
	cells->next = second;
	while (second->next && second->marko <= first->marko) 
	  second = second->next;
	cells = second;
	second = second->next;
      }
    }

  if (!first) 
    cells->next = second;
  if (!second) 
    cells->next = first;
  
  return res;
}

int
update_both(TP_pivot_cell cell, TP_schur_matrix matrix, int *cols_count, int *rows_count, int base)
{
  int cells = 0;
  while ( cell)
    {
      update_counter(cols_count, matrix->CSR[cell->row].col, matrix->CSR[cell->row].nb_elem, base);
      update_counter(rows_count, matrix->CSC[cell->col].row, matrix->CSC[cell->col].nb_elem, base);
      cells++;
      cell = cell->next;
    }
  return cells;
}

TP_pivot_set
merge_to_larger_set(TP_pivot_set self, TP_schur_matrix matrix, TP_solver solver)
{
  TP_pivot_set large, small;
  TP_pivot_cell cells_to_merge;

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
  small->cells   = NULL;

  cells_to_merge = get_independent_cells(cells_to_merge, large->cols_count, large->rows_count, large->base);
  large->nb_elem += update_both(cells_to_merge, matrix, large->cols_count, large->rows_count, large->base);
  large->cells = TP_cells_merge(large->cells, cells_to_merge);

  /* while(cells_to_merge) */
  /*   { */
  /*     TP_pivot_cell cell = cells_to_merge; */
  /*     cells_to_merge = cells_to_merge->next; */
  /*     cell->next = NULL; */

  /*     if ( (add_cell_to_sorted_set(large, cell, matrix)) ) */
  /* 	TP_pivot_cell_destroy(cell); */
  /*   } */

  TP_pivot_set_destroy(small, solver);
  large->next = NULL;
  return large;
}
