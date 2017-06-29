#include <stdlib.h>
#include <stdio.h>

#include "TP_matrix.h"
#include "TP_auxiliary.h"

#include "TP_pivot_list.h"

TP_pivot_list 
TP_pivot_list_create()
{
  TP_pivot_list self = calloc( (size_t) 1, (size_t) sizeof(*self));
   return self;
}


TP_pivot_set 
TP_pivot_set_create(TP_solver solver, int n, int m)
{
  TP_pivot_set self = calloc( (size_t) 1, (size_t) sizeof(*self));
  TP_solver_alloc_counters(solver, &self->cols_count, &self->rows_count);
  self->base = solver->done_pivots;

   return self;
}

TP_pivot_cell
TP_pivot_cell_create(int row, int col, int marko)
{
  TP_pivot_cell self = calloc((size_t) 1, sizeof(*self));
  self->row    = row;
  self->col    = col;
  self->marko  = marko;
  self->next   = NULL;

  return self;
}


TP_pivot_list 
TP_pivot_list_insert_set(TP_pivot_list self, TP_pivot_set set)
{
  TP_pivot_set  tmp;

  tmp = self->sets;
  self->sets = set;
  set->next = tmp;
  self->nb_elem++;
  
  return self;
}


TP_pivot_set
get_next_merging_set(TP_pivot_list self)
{
  TP_pivot_set set;

  if(self->nb_elem < 3) {
    self->nb_elem = 0;
    set = self->sets;
    self->sets = NULL;
    return set;
  } else {
    set = self->sets;
    self->sets = self->sets->next->next;
    set->next->next = NULL;
    self->nb_elem -= 2;
  }

  return set;
}

/* TODO: clean up */
/* it seams that this fucntion merge the two first 
   sets in the list of sets, but does not check if they are 
   indpenedent or no. After this the funtion the  
   get_independent_pivots should be called */
TP_pivot_set 
merge_sorted_sets(TP_pivot_set self)
{
  TP_pivot_cell cell1, cell2, merged;

  if (!self ) 
    return NULL;
  if (!self->next  ) 
    return self;

  
  cell1 = self->cells;
  cell2 = self->next->cells;
  
  // init self->cells to the start of the new merged set
  // and merged is used to add new elements
  if (cell1->marko <= cell2->marko) 
    {
      self->cells = merged = cell1;
      cell1 = cell1->next;
    }
  else 
    {
      self->cells = merged = cell2;
      cell2 = cell2->next;
    }
  
  while ( cell1 && cell2) 
    {
      TP_pivot_cell  new_cell;
      if (cell1->marko < cell2->marko) 
	{
	  new_cell = cell1;
	  cell1 = cell1->next;
	}
      else 
	{
	  new_cell = cell2;
	  cell2 = cell2->next;
	}
      merged->next = new_cell; 
      merged = new_cell; 
    }
  
  if (!cell1) 
    merged->next = cell2; 
  if (!cell2) 
    merged->next = cell1;

  free(self->next);
  self->next = NULL;

  return self;
}


/* TODO: clean up */
/* TP_pivot_set */
/* get_independent_pivots(TP_pivot_set candidates, TP_schur_matrix matrix) */
/* { */
/*   TP_pivot_cell independent_cell; */
/*   TP_pivot_cell candidates_cell; */
/*   struct CSC_struct *CSC = matrix->CSC; */
/*   struct CSR_struct *CSR = matrix->CSR; */

/*   // if candidate is empty, return  */
/*   if (!candidates)  */
/*     return NULL; */
/*   candidates_cell = candidates->cells; */

/*   independent_cell = candidates_cell; */
/*   candidates_cell = candidates_cell->next; */
/*   independent_cell->next = NULL; */

/*   while (candidates_cell)  */
/*     { */
/*       int i, j, independent = 1; */
/*       int row = candidates_cell->row; */
/*       int col = candidates_cell->col; */
/*       int *rows = matrix->row + CSC[col].offset; */
/*       int *cols = matrix->col + CSR[row].offset; */
/*       int nb_rows = CSC[col].nb_elem; */
/*       int nb_cols = CSR[row].nb_elem; */

/*       for(i = 0; i < nb_rows; i++) */
/* 	{ */
/* 	  TP_pivot_cell tmp = independent_cell; */
/* 	  while(tmp)  */
/* 	    { */
/* 	      if(tmp->row == rows[i]) { */
/* 		independent = 0; */
/* 		break; */
/* 	      } */
/* 	      tmp = tmp->next; */
/* 	    } */
/* 	  if (!independent)  */
/* 	    break; */
/* 	} */

/*       for(j = 0; j < nb_cols; j++) */
/* 	{ */
/* 	  TP_pivot_cell tmp = independent_cell; */
/* 	  while(tmp)  */
/* 	    { */
/* 	      if(tmp->col == cols[j]) { */
/* 		independent = 0; */
/* 		break; */
/* 	      } */
/* 	      tmp = tmp->next; */
/* 	    } */
/* 	  if (!independent)  */
/* 	    break; */
/* 	} */
      
/*       // if the candidate is independet, put it on the end of the cells */
/*       if(independent) { */
/* 	TP_pivot_cell tmp = independent_cell;  */
/* 	while(tmp->next) */
/* 	  tmp = tmp->next; */
/* 	tmp->next = candidates_cell; */
/* 	candidates_cell = candidates_cell->next; */
/* 	tmp->next->next = NULL; */
/*       } else {  */
/* 	TP_pivot_cell tmp = candidates_cell; */
/* 	candidates_cell = candidates_cell->next; */
/* 	free(tmp); */
/*       }	 */
/*     } */

/*   candidates->cells = independent_cell; */
/*   return candidates; */
/* } */



void
print_pivot_list(TP_pivot_list self, char *mess)
{
  TP_pivot_set print_set = self->sets;

  printf("%s\n", mess);
  printf("List contains %d elements\n", self->nb_elem);
  while(print_set)
    {
      TP_pivot_cell cells = print_set->cells;
      printf("(( ");
      while( cells ) 
	{
	  printf("[col=%d, row=%d, marko=%d]  ", cells->col, cells->row, cells->marko);
	  cells = cells->next;
	}
      printf("))\n");
      print_set = print_set->next;
    }
}
void 
TP_pivot_cell_destroy(TP_pivot_cell self)
{
  free(self);
}

void
TP_pivot_set_destroy(TP_pivot_set self, TP_solver solver)
{
  TP_solver_dealloc_counters(solver, self->cols_count, self->rows_count);

  while(self->cells) {
    TP_pivot_cell tmp = self->cells;
    self->cells = self->cells->next;
    free(tmp);
  }
  free(self);
}

void
TP_pivot_list_destroy(TP_pivot_list self, TP_solver solver)
{
  TP_pivot_set set = self->sets;

  /* while the set is not empty */
  while(set)
    {
      TP_pivot_set  tmp = set->next;

      TP_pivot_set_destroy(set, solver);
      set = tmp;
    }

  //destroy list
  free(self);
}

