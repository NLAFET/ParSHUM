#ifndef _ParSHUM_ZOLTAN_H 
#define _ParSHUM_ZOLTAN_H 

#include <zoltan.h>
#include "ParSHUM_SBBD_util.h"
#include "ParSHUM_solver.h"
#include "ParSHUM_schur_matrix.h" 

typedef struct _Zoltan_Hypergraph *Zoltan_Hypergraph; 

Zoltan_Hypergraph  ParSHUM_Zoltan_create(ParSHUM_MPI_info MPI_info);

int  ParSHUM_Zoltan_register_data(Zoltan_Hypergraph self, ParSHUM_schur_matrix matrix);
void ParSHUM_Zoltan_partition(Zoltan_Hypergraph self, ParSHUM_schur_matrix matrix);
void ParSHUM_Zoltan_get_row_blocks(Zoltan_Hypergraph hypergraph, row_block row_blocks);
ParSHUM_matrix ParSUM_Zoltan_distribute(ParSHUM_schur_matrix matrix, row_block row_blocks,
					col_block col_blocks, ParSHUM_solver solver, ParSHUM_MPI_info MPI_info);
void ParSHUM_Zoltan_print_distribution(Zoltan_Hypergraph self);
void  ParSHUM_Zoltan_destroy(Zoltan_Hypergraph self);

#endif // _ParSHUM_ZOLTAN_H 
