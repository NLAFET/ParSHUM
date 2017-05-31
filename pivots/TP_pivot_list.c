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
TP_pivot_set_create(int n, int m)
{
  TP_pivot_set self = calloc( (size_t) 1, (size_t) sizeof(*self));
  self->rows_count  = calloc((size_t) m, (size_t) sizeof(*self->rows_count));
  self->cols_count  = calloc((size_t) n, (size_t) sizeof(*self->cols_count));

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


TP_pivot_list TP_pivot_list_insert_new_set(TP_pivot_list self, TP_schur_matrix matrix, 
					   int row, int col, int marko)
{
  int m = matrix->m, n = matrix->n;
  TP_pivot_set  set  = TP_pivot_set_create(n, m);
  TP_pivot_cell cell = malloc( (size_t) sizeof(*cell));

  set->cells   = cell;
  set->next    = NULL;
  set->nb_elem = 1;
  update_counter(set->cols_count, matrix->col + matrix->CSR[row].offset, matrix->CSR[row].nb_elem);
  update_counter(set->rows_count, matrix->row + matrix->CSC[col].offset, matrix->CSC[col].nb_elem);

  cell->row   = row;
  cell->col   = col;
  cell->marko = marko;
  cell->next  = NULL;
  
  if (self->first) {  
    TP_pivot_set  tmp = self->last; 
    tmp->next = set;
    self->last = set;
  }  else 
    self->first = self->midle = self->last = set;
  
  // the midle pointer should be moved only on even number of elements 
  if ( ! ( ++self->nb_elem % 2 ) )  
    self->midle = self->midle->next;
  
  return self;
}


TP_pivot_list 
TP_pivot_list_insert_set(TP_pivot_list self, TP_pivot_set set)
{
  if (self->first) {  
    TP_pivot_set  tmp = self->last; 
    tmp->next  = set;
    self->last = set;
    set->next  = NULL;
  }  else 
    self->first = self->midle = self->last = set;
  
  // the midle pointer should be moved only on even number of elements 
  if ( ! ( ++self->nb_elem % 2 ) )  
    self->midle = self->midle->next;
  
  return self;
}


TP_pivot_set
get_next_merging_set(TP_pivot_list self)
{
  TP_pivot_set first  = self->first;
  TP_pivot_set second = self->midle;

  if (self->nb_elem <= 2) { 
    self->first = self->midle = self->last = NULL;
    if (self->nb_elem == 0) {
      first = NULL;
    } else if (self->nb_elem == 1) {
      first->next = NULL;
    } else if (self->nb_elem == 2) {
      first->next  = second;
      second->next = NULL;
    }
    self->nb_elem = 0;

  } else {

    // handle the second
    TP_pivot_set tmp = self->first; 
    while( tmp->next != second) 
      tmp = tmp->next;
    tmp->next = second->next;
    self->midle = second->next;
    second->next = NULL;
    
    // handle the first
    self->first = first->next;
    first->next = second;
    self->nb_elem -= 2;
  }

  return first;
}


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


TP_pivot_set
get_independent_pivots(TP_pivot_set candidates, TP_schur_matrix matrix)
{
  TP_pivot_cell independent_cell;
  TP_pivot_cell candidates_cell;
  struct CSC_struct *CSC = matrix->CSC;
  struct CSR_struct *CSR = matrix->CSR;

  // if candidate is empty, return 
  if (!candidates) 
    return NULL;
  candidates_cell = candidates->cells;

  independent_cell = candidates_cell;
  candidates_cell = candidates_cell->next;
  independent_cell->next = NULL;

  while (candidates_cell) 
    {
      int i, j, independent = 1;
      int row = candidates_cell->row;
      int col = candidates_cell->col;
      int *rows = matrix->row + CSC[col].offset;
      int *cols = matrix->col + CSR[row].offset;
      int nb_rows = CSC[col].nb_elem;
      int nb_cols = CSR[row].nb_elem;

      for(i = 0; i < nb_rows; i++)
	{
	  TP_pivot_cell tmp = independent_cell;
	  while(tmp) 
	    {
	      if(tmp->row == rows[i]) {
		independent = 0;
		break;
	      }
	      tmp = tmp->next;
	    }
	  if (!independent) 
	    break;
	}

      for(j = 0; j < nb_cols; j++)
	{
	  TP_pivot_cell tmp = independent_cell;
	  while(tmp) 
	    {
	      if(tmp->col == cols[j]) {
		independent = 0;
		break;
	      }
	      tmp = tmp->next;
	    }
	  if (!independent) 
	    break;
	}
      
      // if the candidate is independet, put it on the end of the cells
      if(independent) {
	TP_pivot_cell tmp = independent_cell; 
	while(tmp->next)
	  tmp = tmp->next;
	tmp->next = candidates_cell;
	candidates_cell = candidates_cell->next;
	tmp->next->next = NULL;
      } else { 
	TP_pivot_cell tmp = candidates_cell;
	candidates_cell = candidates_cell->next;
	free(tmp);
      }	
    }

  candidates->cells = independent_cell;
  return candidates;
}



void
print_pivot_list(TP_pivot_list self, char *mess)
{
  TP_pivot_set print_set = self->first;

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
TP_pivot_set_destroy(TP_pivot_set self)
{
  free(self->rows_count);
  free(self->cols_count);
  while(self->cells) {
    TP_pivot_cell tmp = self->cells;
    self->cells = self->cells->next;
    free(tmp);
  }
  free(self);
}

void
TP_pivot_list_destroy(TP_pivot_list self)
{
  TP_pivot_set set = self->first;

  // while the set is not empty
  while(set) 
    {
      TP_pivot_set  tmp = set->next; 

      // destroy set
      TP_pivot_set_destroy(set);
      set = tmp;
    }

  //destroy list
  free(self);
}

