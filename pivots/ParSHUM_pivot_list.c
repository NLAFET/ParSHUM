#include <stdlib.h>
#include <stdio.h>

#include "ParSHUM_matrix.h"
#include "ParSHUM_auxiliary.h"

#include "ParSHUM_pivot_list.h"

ParSHUM_pivot_list 
ParSHUM_pivot_list_create()
{
  ParSHUM_pivot_list self = calloc( (size_t) 1, (size_t) sizeof(*self));
   return self;
}


ParSHUM_pivot_set 
ParSHUM_pivot_set_create(ParSHUM_solver solver, int n, int m)
{
  ParSHUM_pivot_set self = calloc( (size_t) 1, (size_t) sizeof(*self));
  ParSHUM_solver_alloc_counters(solver, &self->cols_count, &self->rows_count);
  self->base = solver->done_pivots + 1;

   return self;
}

ParSHUM_pivot_cell
ParSHUM_pivot_cell_create(int row, int col, int marko)
{
  ParSHUM_pivot_cell self = calloc((size_t) 1, sizeof(*self));
  self->row    = row;
  self->col    = col;
  self->marko  = marko;
  self->next   = NULL;

  return self;
}


ParSHUM_pivot_list 
ParSHUM_pivot_list_insert_set(ParSHUM_pivot_list self, ParSHUM_pivot_set set)
{
  ParSHUM_pivot_set  tmp;

  tmp = self->sets;
  self->sets = set;
  set->next = tmp;
  self->nb_elem++;
  
  return self;
}


ParSHUM_pivot_set
get_next_merging_set(ParSHUM_pivot_list self)
{
  ParSHUM_pivot_set set;

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

void
print_pivot_list(ParSHUM_pivot_list self, char *mess)
{
  ParSHUM_pivot_set print_set = self->sets;

  printf("%s\n", mess);
  printf("List contains %d elements\n", self->nb_elem);
  while(print_set)
    {
      ParSHUM_pivot_cell cells = print_set->cells;
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
ParSHUM_pivot_cell_destroy(ParSHUM_pivot_cell self)
{
  free(self);
}

void
ParSHUM_pivot_set_destroy(ParSHUM_pivot_set self, ParSHUM_solver solver)
{
  ParSHUM_solver_dealloc_counters(solver, self->cols_count, self->rows_count);

  while(self->cells) {
    ParSHUM_pivot_cell tmp = self->cells;
    self->cells = self->cells->next;
    free(tmp);
  }
  free(self);
}

void
ParSHUM_pivot_list_destroy(ParSHUM_pivot_list self, ParSHUM_solver solver)
{
  ParSHUM_pivot_set set = self->sets;

  /* while the set is not empty */
  while(set)
    {
      ParSHUM_pivot_set  tmp = set->next;

      ParSHUM_pivot_set_destroy(set, solver);
      set = tmp;
    }

  //destroy list
  free(self);
}

