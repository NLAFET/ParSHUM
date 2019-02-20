#include <limits.h>
#include "ParSHUM_Zoltan.h"
#include "ParSHUM_solver.h"
#include "ParSHUM_auxiliary.h"
#include "ParSHUM_SBBD_util.h"
#include "ParSHUM_SBBD.h"


struct _ParSHUM_SBBD {  
  ParSHUM_MPI_info MPI_info;
  Zoltan_Hypergraph hypergraph;
  row_block row_blocks;
  col_block col_blocks;
  ParSHUM_solver solver;
  ParSHUM_schur_matrix A;
  char *matrix_file;
  ParSHUM_dense_matrix Schur;
};

ParSHUM_SBBD
ParSHUM_SBBD_create(MPI_Comm world)
{
  ParSHUM_SBBD self = calloc(1, sizeof(*self));
  int is_init;

  self->MPI_info = calloc(1, sizeof(*self->MPI_info));
  MPI_Initialized(&is_init);
  if (!is_init) {
    MPI_Init(NULL, NULL);
    self->MPI_info->world = MPI_COMM_WORLD;
  } else {
    self->MPI_info->world = world;
  }
  MPI_Comm_rank(self->MPI_info->world, &self->MPI_info->rank);
  MPI_Comm_size(self->MPI_info->world, &self->MPI_info->MPI_size);
  
  self->row_blocks     = calloc(1, sizeof(*self->row_blocks));
  self->col_blocks     = calloc(1, sizeof(*self->col_blocks));
  self->hypergraph     = ParSHUM_Zoltan_create(self->MPI_info);
  self->solver         = ParSHUM_solver_create();

  return self;
}

void
ParSHUM_SBBD_parse_args(ParSHUM_SBBD self, int argc, char **argv)
{
  ParSHUM_solver_parse_args(self->solver, argc, argv, 1);
  self->matrix_file = self->solver->exe_parms->matrix_file;
}

void
ParSUHM_SBBD_read_matrix(ParSHUM_SBBD self)
{
  if (!self->MPI_info->rank) {
    ParSHUM_matrix input_A;
    char *file_ext = strrchr(self->matrix_file, '.');
    input_A = ParSHUM_matrix_create();

    if (!strcmp(file_ext, ".mtl"))
      ParSHUM_read_mtl_file(input_A, self->matrix_file);
#ifdef HAVE_SPRAL
    else if (!strcmp(file_ext, ".rb"))
      ParSHUM_read_rutherford_boeing(input_A, self->matrix_file);
#endif
    else
      ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"unsupported matrix file");

    self->A = ParSHUM_schur_matrix_create();

    ParSHUM_schur_matrix_allocate(self->A, input_A->n, input_A->m, input_A->nnz, 0, 0, 0, 0.0, 0.0);
    ParSHUM_schur_matrix_copy(input_A, self->A, 0.0);
    ParSHUM_matrix_destroy(input_A);
  }
}

void
ParSHUM_SBBD_partition(ParSHUM_SBBD self)
{
  char mess[2048];
  ParSHUM_Zoltan_register_data(self->hypergraph, self->A);

  ParSHUM_Zoltan_partition(self->hypergraph, self->A);

  ParSHUM_Zoltan_get_row_blocks(self->hypergraph, self->row_blocks);
  if (!self->MPI_info->rank)  {
    ParSHUM_get_col_blocks(self->A, self->col_blocks, self->row_blocks);
    
    ParSHUM_check_blocks(self->A, self->row_blocks, self->col_blocks);
    ParSHUM_blocks_print_stats(self->A, self->row_blocks, self->col_blocks);
    self->Schur = ParSHUM_dense_matrix_create(self->col_blocks->nb_BB_cols, self->col_blocks->nb_BB_cols);
    /* ParSHUM_print_blocks(self->row_blocks, self->col_blocks); */
  }
  self->solver->A = ParSUM_Zoltan_distribute(self->A, self->row_blocks, self->col_blocks, self->solver, self->MPI_info);
  /* for ( i = 0; i < self->MPI_info->MPI_size; i++)  */
  /*   if ( i == self->MPI_info->rank) { */
  /*     snprintf(mess, 2048, "matrix A on proc %d", self->MPI_info->rank); */
  /*     ParSHUM_matrix_print(self->solver->A, mess); */
  /*     MPI_Barrier(self->MPI_info->world); */
  /*   } */
  
  ParSHUM_solver_init(self->solver);
}

void 
ParSHUM_SBBD_factorize(ParSHUM_SBBD self)
{
  ParSHUM_vector X, B, SOL;
  ParSHUM_dense_matrix local_S;
  int nb_BB_cols = self->solver->BB_cols;
  int rank = self->MPI_info->rank;

  X    = ParSHUM_vector_create(self->solver->A->n);
  SOL  = ParSHUM_vector_create(self->solver->A->m);  
  B    = ParSHUM_vector_create(self->solver->A->m);  
  
  ParSHUM_vector_read_file(X, self->solver->exe_parms->RHS_file); 
  /* ParSHUM_vector_print(X, "the real X"); */
  fflush(stdout);
  ParSHUM_matrix_SpMV(self->solver->A, X, SOL);
  
  /* copy the vector SOL in B */
  ParSHUM_vector_copy(SOL, B);
  
  /* self->solver->debug |= ParSHUM_CHECK_PIVOTS; */
  /* self->solver->debug |= ParSHUM_CHECK_SCHUR_MEMORY; */
  /* self->solver->debug |= ParSHUM_CHECK_SCHUR_SYMETRY; */
  /* self->solver->debug |= ParSHUM_CHECK_COUNTERS; */
  /* self->solver->exe_parms->density_tolerance = 1.0; */
  /* self->solver->exe_parms->min_pivot_per_steps = 5; */

  ParSHUM_solver_factorize(self->solver);

  local_S = self->solver->S_dense;
  if (rank ) {
    ParSHUM_collect_BB_block(&local_S->val[(local_S->n - nb_BB_cols)  * local_S->m], NULL, NULL, NULL,
  			     local_S->m, local_S->n, nb_BB_cols, self->MPI_info);
  } else {
    ParSHUM_collect_BB_block(&local_S->val[(local_S->n - nb_BB_cols)  * local_S->m],
			     self->Schur->val, self->col_blocks, self->row_blocks,
			     local_S->m, local_S->n, nb_BB_cols, self->MPI_info);
    /* ParSHUM_dense_matrix_print(self->Schur, "dense global schur"); */

    int i, j, m = self->Schur->m, n = self->Schur->n; 
    printf("n = %d and m = %d\n", n, m);
    for(i = 0; i < n; i++) {
      double *tmp = &self->Schur->val[i*m];
      int found = 0;
      for(j = 0; j < m; j++) 
	if (tmp[j] != 0)
	  found++;
      if (!found) 
	printf("KO!\n");
    }
				     
    ParSHUM_dense_matrix_factorize(self->Schur, 0, self->solver->exe_parms->nb_threads);
  }



  /* Perform the solve operation */
  /* ParSHUM_solver_solve(self->solver, B); */
  /* ParSHUM_vector_print(B, "computed solution"); */
  /* fflush(stdout);       */
  /* compute the norms */
  /* ParSHUM_solver_compute_norms(self->solver, B, SOL); */
  
  /* ParSHUM_solver_finalize(self->solver); */
}

void
ParSHUM_SBBD_finalize(ParSHUM_SBBD self)
{
  /* ParSHUM_solver_finalize(self->solver); */
}

void
ParSHUM_SBBD_destroy(ParSHUM_SBBD self)
{
  if (!self->MPI_info->rank) 
    ParSHUM_schur_matrix_destroy(self->A);

  ParSHUM_Zoltan_destroy(self->hypergraph);

  ParSHUM_solver_destroy(self->solver);

  free(self->row_blocks);
  free(self->col_blocks);
  free(self->MPI_info);

  free(self);
}
