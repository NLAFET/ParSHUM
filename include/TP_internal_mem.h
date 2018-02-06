#ifndef   _TP_INTERNAL_MEM_H 
#define   _TP_INTERNAL_MEM_H

#include <omp.h>

typedef struct _TP_internal_mem *TP_internal_mem;
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
  int    nb_eligible;

  double *val;
  int    *row;
};


TP_internal_mem TP_internal_mem_create(long nnz);
void            TP_internal_mem_destroy(TP_internal_mem self);
void            TP_internal_mem_CSC_alloc(TP_internal_mem self, CSC_struct *CSC, long size);
void            TP_internal_mem_CSR_alloc(TP_internal_mem self, CSR_struct *CSR, long size);

#define TP_internal_mem_CSC_realloc(self, CSC) TP_internal_mem_CSC_alloc(self, CSC, 2*(CSC->nb_elem+CSC->nb_free))
#define TP_internal_mem_CSR_realloc(self, CSR) TP_internal_mem_CSR_alloc(self, CSR, 2*(CSR->nb_elem+CSR->nb_free))
 

#endif // _TP_INTERNAL_MEM_H
