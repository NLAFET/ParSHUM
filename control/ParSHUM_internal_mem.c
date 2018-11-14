#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "ParSHUM_internal_mem.h"

typedef struct _free_space  *free_space;

struct _free_space {
  size_t size; 
  size_t offset; 
  free_space next;
  free_space previous;
};

struct _ParSHUM_internal_mem {
  void  **mem;
  
  free_space *unused;
  omp_lock_t lock;

  int nb_mems;
  size_t init_size;
  size_t alignment;
};


ParSHUM_internal_mem
ParSHUM_internal_mem_create(size_t size, size_t alignment)
{
  ParSHUM_internal_mem self = calloc(1, sizeof(*self));

  self->init_size = size; 

  self->mem = malloc(sizeof(*self->mem));
  posix_memalign(self->mem, alignment, size);
  self->unused = malloc(sizeof(*self->unused));
  *self->unused = calloc(1, sizeof(**self->unused));
  self->unused[0]->size = size;

  self->nb_mems = 1;
  self->alignment = alignment;
  omp_init_lock(&self->lock);

  return self;
}


void
ParSHUM_internal_mem_realloc(ParSHUM_internal_mem self)
{
  int i;
  size_t size = self->init_size;
  for( i = 1; i < self->nb_mems; i++)
    size *= 2;

  self->unused  = realloc((void *) self->unused, (size_t) (self->nb_mems + 1)   * sizeof(*self->unused));
  self->unused[self->nb_mems]  = calloc(1, sizeof(**self->unused));
  self->unused[self->nb_mems]->size = size;

  self->mem = realloc((void *) self->mem, (size_t) (self->nb_mems + 1) * sizeof(*self->mem));
  posix_memalign(&self->mem[self->nb_mems++], self->alignment, size);
}

size_t
ParSHUM_internal_mem_alloc(ParSHUM_internal_mem self, void **mem, size_t size)
{
  omp_set_lock(&self->lock);
  int i, nb = self->nb_mems; 

  /* first try to find memory in the current space  */
  for( i = 0; i < nb; i++)
    {
      free_space current_unused = self->unused[i];
      void *internal_mem = self->mem[i];
      while(current_unused) {
	if (current_unused->size >=size ) {
	  *mem = internal_mem + current_unused->offset;
	  current_unused->offset  += size;
	  current_unused->size    -= size;
	  omp_unset_lock(&self->lock);
	  return size;
	}
	current_unused = current_unused->next;
      }
    }

  /* if there is no spaace then we reallocate */
  while (1) {
    ParSHUM_internal_mem_realloc(self);
    free_space current_unused = self->unused[i];
    void *internal_mem = self->mem[i];
    if (current_unused->size >= size ) {
      *mem = internal_mem + current_unused->offset;
      current_unused->offset  += size;
      current_unused->size -= size;
      omp_unset_lock(&self->lock);
      return size;
    }
    i++;
  }
}


void
ParSHUM_internal_mem_destroy(ParSHUM_internal_mem self)
{
  int nb_mems = self->nb_mems, i ;
  
  for (i = 0; i < nb_mems; i++) 
    {
      free_space current_unused = self->unused[i];
      while (current_unused) {
	free_space tmp = current_unused->next;
	free(current_unused);
	current_unused = tmp;
      }
      free(self->mem[i]);
    }
  free(self->mem);
  free(self->unused);

  omp_destroy_lock(&self->lock);
  free(self);
}
