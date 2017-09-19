#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>
#include <math.h>

#include "TP_solver.h" 
#include "TP_matrix.h" 
#include "TP_dense.h"
#include "TP_enum.h"
#include "TP_pivot_list.h" 
#include "TP_auxiliary.h"


int
main(int argc, char **argv) 
{
  TP_matrix A = TP_matrix_create();
  int i;
  long j;
  FILE *file;
  char *file_ext = strrchr(argv[1], '.');

  if (!strcmp(file_ext, ".mtl"))
    TP_read_mtl_file(A, argv[1]);
  else 
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"the input matrix should be in mtl format");

  file = fopen(argv[2], "r");
  if (!file)
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"error while opening the  matrix file");
  
  fprintf(file,"%d\t%ld\n", A->n, A->nnz);
  for( i = 0; i < A->n; i++) 
    for( j = A->col_ptr[i]; j < A->col_ptr[i+1]; j++)
      fprintf(file,"%d\t%d\t%lf\n", i, A->row[j], A->val[j]);

  fclose(file);
  TP_matrix_destroy(A);

  return 0;
}
