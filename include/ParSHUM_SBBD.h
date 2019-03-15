#ifndef _ParSHUM_SBBD_H 
#define _ParSHUM_SBBD_H 

#include "ParSHUM_Zoltan.h"
#include "ParSHUM_solver.h"
#include "ParSHUM_auxiliary.h"
#include "ParSHUM_SBBD_util.h"
#include "ParSHUM_SBBD_util.h"

struct _ParSHUM_SBBD {  
  ParSHUM_MPI_info MPI_info;
  Zoltan_Hypergraph hypergraph;
  row_block row_blocks;
  col_block col_blocks;
  ParSHUM_solver solver;
  ParSHUM_schur_matrix A;
  ParSHUM_matrix input_A;
  char *matrix_file;
  ParSHUM_dense_matrix Schur;
};

typedef struct _ParSHUM_SBBD *ParSHUM_SBBD; 

ParSHUM_SBBD ParSHUM_SBBD_create(MPI_Comm world);
  
void  ParSHUM_SBBD_parse_args(ParSHUM_SBBD self, int argc, char **argv);
void  ParSUHM_SBBD_read_matrix(ParSHUM_SBBD self);
void  ParSHUM_SBBD_partition(ParSHUM_SBBD self);
void  ParSHUM_SBBD_factorize(ParSHUM_SBBD self);
void  ParSHUM_SBBD_solve(ParSHUM_SBBD self, ParSHUM_vector RHS);
void  ParSHUM_SBBD_finalize(ParSHUM_SBBD self);
void  ParSHUM_SBBD_destroy(ParSHUM_SBBD self);

#endif // _ParSHUM_SBBD_H 
