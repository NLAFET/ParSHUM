#include <stdlib.h>

#include "TP_verbose.h"


TP_verbose 
TP_verbose_create_V0(int verbosity, TP_matrix matrix)
{
  TP_verbose self = NULL;
  
  self = calloc((size_t) 1, sizeof(*self));
  
  self->verbosity = verbosity;
  self->nnz_input = self->nnz_final = matrix->nnz;
  self->n = matrix->n;
  self->m = matrix->m;

  return self;
}

void
TP_verbose_step_start_V0(TP_verbose self)
{
  per_step step = calloc((size_t) 1, sizeof(*step));
  
  if (!self->stats_first) {
    self->stats_first = self->stats_last = step;
  } else {
    self->stats_last->next = step;
    self->stats_last = step;
  }
  
  self->nb_steps++;
}

void 
TP_verbose_destroy_V0(TP_verbose self)
{
  per_step step = self->stats_first;
  
  while(step)
    {
      per_step tmp = step->next;
      free(step);
      step = tmp;
    }
  free(self);
}
