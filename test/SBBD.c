#include "ParSHUM_SBBD.h"
#include <mpi.h>

int
main(int argc, char **argv)
{
  ParSHUM_SBBD solver;
  
  solver = ParSHUM_SBBD_create(0);

  ParSHUM_SBBD_parse_args(solver, argc, argv);
  ParSUHM_SBBD_read_matrix(solver);
  ParSHUM_SBBD_partition(solver);

  ParSHUM_SBBD_factorize(solver);
  /* ParSHUM_SBBD_destroy(solver); */

  MPI_Finalize();

  return 0;
}
