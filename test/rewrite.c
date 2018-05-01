#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>
#include <math.h>

#include "ParSHUM_solver.h" 
#include "ParSHUM_matrix.h" 
#include "ParSHUM_dense.h"
#include "ParSHUM_enum.h"
#include "ParSHUM_pivot_list.h" 
#include "ParSHUM_auxiliary.h"


int
main(int argc, char **argv) 
{
  ParSHUM_matrix A = ParSHUM_matrix_create();
  int i;
  long j;
  FILE *file;
  char *file_ext = strrchr(argv[1], '.');

  if (!strcmp(file_ext, ".mtl"))
    ParSHUM_read_mtl_file(A, argv[1]);
  else 
    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"the input matrix should be in mtl format");

  file = fopen(argv[2], "r");
  if (!file)
    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"error while opening the  matrix file");
  
  fprintf(file,"%d\t%ld\n", A->n, A->nnz);
  for( i = 0; i < A->n; i++) 
    for( j = A->col_ptr[i]; j < A->col_ptr[i+1]; j++)
      fprintf(file,"%d\t%d\t%lf\n", i, A->row[j], A->val[j]);

  fclose(file);
  ParSHUM_matrix_destroy(A);

  return 0;
}
