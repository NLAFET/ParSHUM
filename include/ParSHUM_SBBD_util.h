#ifndef _ParSHUM_SBBD_UTIL_H 
#define _ParSHUM_SBBD_UTIL_H 

#include <mpi.h>
#include "ParSHUM_schur_matrix.h"

struct _row_block {
  int nb_blocks;
  int n;
  
  int *perms;
  int *invr_perms;
  int *sizes;
};

struct _col_block {
  int nb_blocks;
  int n;
  int nb_BB_cols;
  
  int *perms;
  int *invr_perms;
  int *sizes;
  int *nnz;
  int **BB_index;
  int *BB_size;
  int *local_schur_m;
};

struct _ParSHUM_MPI_info {
  MPI_Comm world;
  int MPI_size;
  int rank;
};

typedef struct _row_block  *row_block;
typedef struct _col_block  *col_block;
typedef struct _ParSHUM_MPI_info  *ParSHUM_MPI_info;

void  ParSHUM_get_col_blocks(ParSHUM_schur_matrix A, col_block col_blocks, row_block row_blocks);
ParSHUM_matrix  ParSHUM_get_block(ParSHUM_schur_matrix matrix, row_block row_blocks,
				  col_block col_blocks, int block);
void ParSHUM_collect_BB_block(double *local_schur, double *global_schur, col_block col_blocks,
			      row_block row_blocks, int m, int n, int BB_cols, ParSHUM_MPI_info MPI_info);

void ParSHUM_print_blocks(row_block row_blocks, col_block col_blocks);
void ParSHUM_check_blocks(ParSHUM_schur_matrix A, row_block row_blocks, col_block col_blocks);
void ParSHUM_blocks_print_stats(ParSHUM_schur_matrix A, row_block row_blocks, col_block col_blocks);


#endif // _ParSHUM_SBBD_UTIL_H
