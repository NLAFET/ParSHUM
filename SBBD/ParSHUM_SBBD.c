#include <limits.h>
#include "ParSHUM_Zoltan.h"
#include "ParSHUM_solver.h"
#include "ParSHUM_auxiliary.h"
#include "ParSHUM_SBBD.h"

struct _ParSHUM_SBBD {
  MPI_Comm world;
  int MPI_size;
  int rank;
  
  Zoltan_Hypergraph hypergraph;
  ParSHUM_solver solver;
  ParSHUM_matrix A;
  char *matrix_file;

  float Zoltan_version;
};

ParSHUM_SBBD
ParSHUM_SBBD_create(MPI_Comm world)
{
  ParSHUM_SBBD self = calloc((size_t) 1, sizeof(*self));
  int is_init;

  MPI_Initialized(&is_init);
  if (!is_init) {
    MPI_Init(NULL, NULL);
    self->world = MPI_COMM_WORLD;
  } else {
    self->world = world;
  }
  MPI_Comm_rank(self->world, &self->rank);
  MPI_Comm_size(self->world, &self->MPI_size);
  
  self->Zoltan_version = ParSHUM_Zoltan_init();
  self->hypergraph     = ParSHUM_Zoltan_create(self->world);
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
  if (!self->rank) {
    char *file_ext = strrchr(self->matrix_file, '.');
    self->A        = ParSHUM_matrix_create();

    if (!strcmp(file_ext, ".rb"))
      ParSHUM_read_rutherford_boeing(self->A, self->matrix_file);
    else if (!strcmp(file_ext, ".mtl"))
      ParSHUM_read_mtl_file(self->A, self->matrix_file);
    else
      ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"unsupported matrix file");
  }
}

void
ParSHUM_SBBD_partition(ParSHUM_SBBD self)
{
  ParSHUM_Zoltan_init_distrubtion(self->hypergraph, self->A);
  ParSHUM_Zoltan_parition(self->hypergraph, self->A);
}

void
ParSHUM_SBBD_destroy(ParSHUM_SBBD self)
{
  if (!self->rank)
    ParSHUM_matrix_destroy(self->A);

  ParSHUM_Zoltan_destroy(self->hypergraph);
  
  free(self->solver);
  free(self);
}
