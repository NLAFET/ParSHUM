#ifndef _ParSHUM_ZOLTAN_H 
#define _ParSHUM_ZOLTAN_H 

#include <zoltan.h>
#include "ParSHUM_matrix.h" 

typedef struct _Zoltan_Hypergraph *Zoltan_Hypergraph; 

float ParSHUM_Zoltan_init(void);

Zoltan_Hypergraph  ParSHUM_Zoltan_create(MPI_Comm comm);

int  ParSHUM_Zoltan_init_distrubtion(Zoltan_Hypergraph self, ParSHUM_matrix matrix);
void ParSHUM_Zoltan_print_distribution(Zoltan_Hypergraph self);
void ParSHUM_Zoltan_parition(Zoltan_Hypergraph self, ParSHUM_matrix A);

void  ParSHUM_Zoltan_destroy(Zoltan_Hypergraph self);

#endif // _ParSHUM_ZOLTAN_H 
