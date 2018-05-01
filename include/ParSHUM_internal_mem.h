#ifndef   _ParSHUM_INTERNAL_MEM_H 
#define   _ParSHUM_INTERNAL_MEM_H

#include <omp.h>

typedef struct _ParSHUM_internal_mem *ParSHUM_internal_mem;
typedef struct _CSC_struct CSC_struct;
typedef struct _CSR_struct CSR_struct;

struct _CSR_struct {
  int nb_elem;
  int nb_free;
  int *col;
};

struct _CSC_struct {
  double col_max;
  int    nb_elem;
  int    nb_free;
  int    nb_numerical_eligible;
  int    nb_eligible;

  double *val;
  int    *row;
};


ParSHUM_internal_mem ParSHUM_internal_mem_create(long nnz);
void                 ParSHUM_internal_mem_destroy(ParSHUM_internal_mem self);
void                 ParSHUM_internal_mem_CSC_alloc(ParSHUM_internal_mem self, CSC_struct *CSC, long size);
void                 ParSHUM_internal_mem_CSR_alloc(ParSHUM_internal_mem self, CSR_struct *CSR, long size);

#define ParSHUM_internal_mem_CSC_realloc(self, CSC) ParSHUM_internal_mem_CSC_alloc(self, CSC, 2*(CSC->nb_elem+CSC->nb_free))
#define ParSHUM_internal_mem_CSR_realloc(self, CSR) ParSHUM_internal_mem_CSR_alloc(self, CSR, 2*(CSR->nb_elem+CSR->nb_free))
 

#endif // _ParSHUM_INTERNAL_MEM_H
