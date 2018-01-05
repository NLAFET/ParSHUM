#include <stdlib.h>
#include <string.h>
#include "TP_internal_mem.h"

typedef struct _free_space  *free_space;

struct _free_space {
  long nb_elem;
  long offset;
  free_space next;
  free_space previous;
};

struct _TP_internal_mem {
  double **val;
  int    **row;
  int    **col;
  
  free_space *unused_CSC;
  free_space *unused_CSR;
  omp_lock_t CSC_lock;
  omp_lock_t CSR_lock;

  int nb_CSC;
  int nb_CSR;
  long init_size;
};


TP_internal_mem
TP_internal_mem_create(long nnz)
{
  TP_internal_mem self = calloc(1, sizeof(*self));

  self->init_size   = nnz;
  self->nb_CSC      = 1;
  self->nb_CSR      = 1;

  self->val  = malloc((size_t) sizeof(*self->val));
  self->row  = malloc((size_t) sizeof(*self->row));
  self->col  = malloc((size_t) sizeof(*self->col));
  *self->val = malloc((size_t) nnz * sizeof(**self->val));
  *self->row = malloc((size_t) nnz * sizeof(**self->row));
  *self->col = malloc((size_t) nnz * sizeof(**self->col));
 
  self->unused_CSC  = malloc(sizeof(*self->unused_CSC));
  *self->unused_CSC = calloc(1, sizeof(**self->unused_CSC));
  self->unused_CSC[0]->nb_elem = nnz;
  
  self->unused_CSR  = malloc(sizeof(*self->unused_CSR));
  *self->unused_CSR = calloc(1, sizeof(**self->unused_CSR));
  self->unused_CSR[0]->nb_elem = nnz;

  omp_init_lock(&self->CSC_lock);
  omp_init_lock(&self->CSR_lock);
  
  return self;
}

/* TODO: nutex  and concurent execution */
void
TP_internal_mem_realloc_CSC(TP_internal_mem self)
{
  int i;
  long size = self->init_size;
  for( i = 1; i < self->nb_CSC; i++)
    size *= 2;

  self->nb_CSC++;
  
  self->unused_CSC  = realloc((void *) self->unused_CSC, (size_t) self->nb_CSC * sizeof(*self->unused_CSC));
  self->unused_CSC[self->nb_CSC - 1]  = calloc(1, sizeof(**self->unused_CSC));
  self->unused_CSC[self->nb_CSC - 1]->nb_elem = size;

  self->val = realloc((void *) self->val, (size_t) self->nb_CSC * sizeof(*self->val));
  self->row = realloc((void *) self->row, (size_t) self->nb_CSC * sizeof(*self->row));
  self->val[self->nb_CSC - 1] = malloc((size_t) size * sizeof(**self->val));
  self->row[self->nb_CSC - 1] = malloc((size_t) size * sizeof(**self->row));
}

/* TODO: mutex  and concurent execution */
void
TP_internal_mem_realloc_CSR(TP_internal_mem self)
{
  int i;
  long size = self->init_size;
  for( i = 1; i < self->nb_CSR; i++)
    size *= 2;

  self->nb_CSR++;
  
  self->unused_CSR  = realloc((void *) self->unused_CSR, (size_t) self->nb_CSR * sizeof(*self->unused_CSR));
  self->unused_CSR[self->nb_CSR - 1]  = calloc(1, sizeof(**self->unused_CSR));
  self->unused_CSR[self->nb_CSR - 1]->nb_elem = size;

  self->col = realloc((void *) self->col, (size_t) self->nb_CSR * sizeof(*self->col));
  self->col[self->nb_CSR - 1] = malloc((size_t) size * sizeof(**self->col));
}

void
TP_internal_mem_CSC_alloc(TP_internal_mem self, CSC_struct *CSC, long size)
{
  omp_set_lock(&self->CSC_lock);
  int i, nb = self->nb_CSC; 

  /* first try to find memory in the current space  */
  for( i = 0; i < nb; i++)
    {
      free_space current_unused = self->unused_CSC[i];
      while(current_unused) {
	if (current_unused->nb_elem >=size ) {
	  if (CSC->nb_elem)  {
	    memcpy(&self->val[i][current_unused->offset],CSC->val,  CSC->nb_elem * sizeof(*CSC->val));
	    memcpy(&self->row[i][current_unused->offset],CSC->row,  CSC->nb_elem * sizeof(*CSC->row));
	  }
	  CSC->nb_free = size - CSC->nb_elem;
	  CSC->val = &self->val[i][current_unused->offset];
	  CSC->row = &self->row[i][current_unused->offset];
	  current_unused->offset  += size;
	  current_unused->nb_elem -= size;
	  omp_unset_lock(&self->CSC_lock);
	  return;
	}
	current_unused = current_unused->next;
      }
    }

  /* if there is no spaace then we reallocate */
  while (1) {
    TP_internal_mem_realloc_CSC(self);
    free_space current_unused = self->unused_CSC[i];
    if (current_unused->nb_elem >=size ) {
      if (CSC->nb_elem)  {
	memcpy(&self->val[i][current_unused->offset], CSC->val,  CSC->nb_elem * sizeof(*CSC->val));
	memcpy(&self->row[i][current_unused->offset], CSC->row,  CSC->nb_elem * sizeof(*CSC->row));
      }
      CSC->nb_free = size - CSC->nb_elem;
      CSC->val = &self->val[i][current_unused->offset];
      CSC->row = &self->row[i][current_unused->offset];
      current_unused->offset  += size;
      current_unused->nb_elem -= size;
      omp_unset_lock(&self->CSC_lock);
      return;
    }
    i++;
  }
}

void
TP_internal_mem_CSR_alloc(TP_internal_mem self, CSR_struct *CSR, long size)
{
  omp_set_lock(&self->CSR_lock);
  int i, nb = self->nb_CSR; 

  /* first try to find memory in the current space  */
  for( i = 0; i < nb; i++)
    {
      free_space current_unused = self->unused_CSR[i];
      while(current_unused) {
	if (current_unused->nb_elem >=size ) {
	  if (CSR->nb_elem)  
	    memcpy(&self->col[i][current_unused->offset], CSR->col, CSR->nb_elem * sizeof(*CSR->col));
	  CSR->nb_free = size - CSR->nb_elem;
	  CSR->col = &self->col[i][current_unused->offset];
	  current_unused->offset  += size;
	  current_unused->nb_elem -= size;
	  omp_unset_lock(&self->CSR_lock);
	  return;
	}
	current_unused = current_unused->next;
      }
    }

  /* if there is no spaace then we reallocate */
  while (1) {
    TP_internal_mem_realloc_CSR(self);
    free_space current_unused = self->unused_CSR[self->nb_CSR - 1];
    if (current_unused->nb_elem >=size ) {
      if (CSR->nb_elem)  
	memcpy(&self->col[self->nb_CSR - 1][current_unused->offset], CSR->col, CSR->nb_elem * sizeof(*CSR->col));
      CSR->nb_free = size - CSR->nb_elem;
      CSR->col = &self->col[self->nb_CSR - 1][current_unused->offset];
      current_unused->offset  += size;
      current_unused->nb_elem -= size;
      omp_unset_lock(&self->CSR_lock);
      return;
    }
  }
}
#define TP_internal_mem_CSR_realloc(self, CSR) TP_internal_mem_CSR_alloc(self, CSR, 2*CSR->nb_elem)

void
TP_internal_mem_destroy(TP_internal_mem self)
{
  int i, nb_CSC = self->nb_CSC, nb_CSR = self->nb_CSR;

  for( i = 0; i < nb_CSC; i++)
    {
      free_space current_unused = self->unused_CSC[i];
      while( current_unused ) { 
	free_space tmp = current_unused->next;
	free(current_unused);
	current_unused = tmp;
      }
      free(self->val[i]);
      free(self->row[i]);
    }
  free(self->val);
  free(self->row);
  free(self->unused_CSC);

  for( i = 0; i < nb_CSR; i++)
    {
      free_space current_unused = self->unused_CSR[i];
      while( current_unused ) { 
	free_space tmp = current_unused->next;
	free(current_unused);
	current_unused = tmp;
      }
      free(self->col[i]);
    }
  free(self->col);
  free(self->unused_CSR);

  omp_destroy_lock(&self->CSC_lock);
  omp_destroy_lock(&self->CSR_lock);
  
  free(self);
}
