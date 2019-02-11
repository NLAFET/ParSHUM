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
  /* ParSHUM_solver_parse_args(self->solver, argc, argv, 0); */
  self->matrix_file = argv[1];
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
  ParSHUM_Zoltan_register_data(self->hypergraph, self->A);

  ParSHUM_Zoltan_partition(self->hypergraph, self->A);

  ParSHUM_Zoltan_get_row_blocks(self->hypergraph, self->row_blocks);
  if (!self->MPI_info->rank)  {
    ParSHUM_get_col_blocks(self->A, self->col_blocks, self->row_blocks);
    
    ParSHUM_check_blocks(self->A, self->row_blocks, self->col_blocks);
    ParSHUM_blocks_print_stats(self->A, self->row_blocks, self->col_blocks);
  }
  self->solver->A = ParSUM_Zoltan_distribute(self->A, self->row_blocks, self->col_blocks, self->MPI_info);
  ParSHUM_matrix_print(self->solver->A, "matrix A");
  /* ParSHUM_solver_init(self->solver); */
}

void 
ParSHUM_SBBD_factorize(ParSHUM_SBBD self)
{
  ParSHUM_solver_factorize(self->solver);
}

void
ParSHUM_SBBD_finalize(ParSHUM_SBBD self)
{
  ParSHUM_solver_finalize(self->solver);
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
