#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <plasma.h>

#include "TP_enum.h" 
#include "TP_matrix.h" 
#include "TP_dense.h" 
#include "TP_schur_matrix.h" 
#include "TP_pivot_list.h" 
#include "TP_auxiliary.h"
#include "TP_verbose.h"
#include "TP_solver.h" 

const char *usageStrign[] = {
  "usage test: [-h] [--matrix (matrix)] [--debug_mode] [--verbosity (level)] [--marko_tol (tol)] [--value_tol (tol)]",
  "            [--extra_space (factor)] [--nb_threads (#threads)] [--sparse_part (factor)] [--nb_init_blocks (#blocks)] ", 
  NULL,
};


void 
print_sthg(int i)
{
  printf("i = %d\n", i);
}

void 
print_sthg2(int i, int j)
{
  printf("i = %d and j = %d \n", i, j);
}

void 
print_string(char *mess)
{
  printf("mess = %s  \n", mess);
}

TP_solver 
TP_solver_create(int verbosity)
{
  TP_solver self = calloc(1, sizeof(*self));

  self->A    = TP_matrix_create();
  self->S    = TP_schur_matrix_create();
  self->L    = TP_matrix_create();
  self->U    = TP_matrix_create();
  self->D    = TP_matrix_create();

  self->verbosity             = verbosity;
  self->nb_threads            = 1;
  self->value_tol             = 0.1;
  self->marko_tol             = 4;
  self->extra_space           = 1.0;
  self->extra_space_inbetween = 1.0;
  self->sparse_part           = 0.5;
  self->nb_init_blocks        = 10;

  return self;
}

void 
parse_args(TP_solver self, 
	   int argc,
	   char **argv)
{
  int i;
  for (i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--verbosity")) {
      self->verbosity = atoi(argv[++i]);
      continue;
    } else if (!strcmp(argv[i], "--matrix")) {
      self->matrix_file = argv[++i];
      continue;
    } else if (!strcmp(argv[i], "--value_tol")) {
      double tmp = atof(argv[++i]);
      if ( tmp > 1.0 || tmp < 0.0) 
	TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "value tolerance should be between 0 and 1");
      self->value_tol = tmp;
      continue;
    } else if (!strcmp(argv[i], "--marko_tol")) {
      double tmp = atof(argv[++i]);
      if ( tmp < 1.0 ) 
	TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "marko tolerance should be larger then 1");
      self->marko_tol = tmp;
      continue;
    } else if (!strcmp(argv[i], "--nb_threads")) {
      int tmp = atoi(argv[++i]);
      self->nb_threads = tmp;
      continue;
    } else if (!strcmp(argv[i], "--extra_space")) {
      double tmp = atof(argv[++i]);
      self->extra_space = tmp;
      continue;
    } else if (!strcmp(argv[i], "--extra_space_inbetween")) {
      double tmp = atof(argv[++i]);
      self->extra_space_inbetween = tmp;
      continue;
    }else if (!strcmp(argv[i], "--sparse_part")) {
      double tmp = atof(argv[++i]);
      if ( tmp > 1.0 ) 
	TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "sparse part can not be larger then 1");
      if ( tmp < 0.0 ) 
	TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "sparse part can not be smaller then 1");
      self->sparse_part = tmp;
      continue;
    }else if (!strcmp(argv[i], "--nb_init_blocks")) {
      int tmp = atoi(argv[++i]);
      if (tmp < 1)
	TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "nb_init_blocks shouls be at least 1");
      self->nb_init_blocks = tmp;
      continue;
    } else if (!strcmp(argv[i], "--debug_mode")) {
      /* self->verbosity = TP_debug; */
      continue;
    } else if (!strcmp(argv[i], "--help")) {
      int i = 0;
      while( usageStrign[i] !=  NULL)
	printf("%s\n", usageStrign[i++]);
      exit(0);
    }// if (self->debug_mode) {
      //  //
    //}
  }
}

void
TP_solver_init(TP_solver self, int argc, char **argv)
{
  parse_args(self, argc, argv);
  TP_solver_read_matrix(self);
  plasma_init(self->nb_threads);

  self->verbose   = DD_verbose_create(self->verbosity, self->A)

  int_array_memset(self->row_perm, -1, self->A->m);
  int_array_memset(self->col_perm, -1, self->A->n);
  self->random_col = create_randomize(self->A->n);
}

void
TP_solver_read_matrix(TP_solver self)
{
  double total_extra_space = 1 + self->extra_space_inbetween + self->extra_space;

  /* TODO:  do this with verbosity  */
  /* if ( !self->debug_mode)  */
  read_rutherford_boeing(self->A, self->matrix_file);
  /* else  */
  /*   self->A = TP_matrix_create_random_matrix(3, 3); */
  
  self->A->type = TP_Rutherford_matrix;
  
  TP_schur_matrix_allocate(self->S, self->A->n, self->A->m, self->A->nnz, 
			   self->extra_space, self->extra_space_inbetween);
  TP_schur_matrix_copy(self->A, self->S);

  // For the sizes of L and U
  // we assume that all diagonal entries are non-zero and   
  TP_matrix_allocate(self->L, self->A->n, 0, (self->A->nnz - self->A->n) / 2, total_extra_space, TP_CSC_matrix);
  TP_matrix_allocate(self->U, 0, self->A->m, (self->A->nnz - self->A->m) / 2, total_extra_space, TP_CSR_matrix);
  TP_matrix_allocate(self->D, self->A->n, 0, 0, 1.0, TP_Diag_matrix);
  
  self->row_perm = malloc((size_t) self->A->m * sizeof(*self->row_perm));
  self->col_perm = malloc((size_t) self->A->n * sizeof(*self->col_perm));
  int_array_memset(self->row_perm, -1, self->A->m);
  int_array_memset(self->col_perm, -1, self->A->n);
}

void

TP_solver_find_pivot_set(TP_solver self)
{
  TP_pivot_list list;
  int new_pivots;

  list = get_possible_pivots(self->S, self->random_col, self->value_tol, self->marko_tol, self->nb_init_blocks);
  
  if( ! list->nb_elem ) 
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "no possible pivot were found");
  list = merge_pivot_sets(list, self->S);

  new_pivots = TP_pivot_get_pivots(list->first, &self->row_perm[self->done_pivots], &self->col_perm[self->done_pivots]);

  TP_pivot_list_destroy(list);
  self->found_pivots += new_pivots;
}

void
TP_solver_update_matrix(TP_solver self)
{
  TP_matrix L = self->L;
  TP_matrix D = self->D;
  TP_matrix U = self->U;
  TP_schur_matrix S = self->S;
  int nb_pivots = self->found_pivots - self->done_pivots;

  TP_schur_matrix_update_LD(S, L, D, &self->row_perm[self->done_pivots], &self->col_perm[self->done_pivots], nb_pivots);
  TP_schur_matrix_update_U (S, U,    &self->row_perm[self->done_pivots], &self->col_perm[self->done_pivots], nb_pivots);
  TP_schur_matrix_update_S (S, L, U, self->done_pivots, self->found_pivots);
  self->done_pivots = self->found_pivots;
  /* TP_print_matrix(L, "L"); */
  /* TP_print_matrix(U, "U"); */
  /* TP_print_matrix(D, "D"); */
  /* TP_schur_matrix_print(S, "S"); */
}


void
TP_solver_factorize(TP_solver self)
{
  int needed_pivots = self->A->n < self->A->m ? self->A->n : self->A->m;

  while (self->done_pivots < (int) needed_pivots * 0.8 )
    { 
      DD_verbose_step_start(self->verbose)
      TP_solver_find_pivot_set(self);
      TP_solver_update_matrix(self);
      /* TP_schur_matrix_check_perms(self->S, self->row_perm, self->col_perm, self->done_pivots); */

      /* printf("checking col_perms\n"); */
      /* check_vlaid_perms(self->col_perm, self->A->n); */
      /* printf("checking row_perms\n"); */
      /* check_vlaid_perms(self->row_perm, self->A->m); */
    }

  /* print_int(self->col_perm, self->A->n, "col_perms"); */
  /* print_int(self->row_perm, self->A->m, "row_perms"); */
  self->S_dense = TP_schur_matrix_convert(self->S, self->done_pivots);
  TP_dense_matrix_factorize(self->S_dense);
   
  memcpy(&self->row_perm[self->done_pivots], self->S_dense->original_rows,
  	 (self->A->m - self->done_pivots) * sizeof(*self->row_perm));
  memcpy(&self->col_perm[self->done_pivots], self->S_dense->original_cols,
  	 (self->A->n - self->done_pivots) * sizeof(*self->col_perm));
  fflush(stdout);
}


void
TP_solver_solve(TP_solver self, TP_vector X, TP_vector rhs)
{
  double *X_vals = X->vect;

  TP_vector_permute(X, self->row_perm);
  TP_matrix_solve_L(self->L, X, self->row_perm);
  
  plasma_dgetrs(self->S_dense->n, 1, self->S_dense->val, self->S_dense->n,
  		self->S_dense->pivots, &X_vals[self->done_pivots], self->S_dense->n);
  
  TP_matrix_solve_UD(self->U, self->D, X, self->col_perm);

  TP_vector_print(X, "x after the U solve");
}



void 
TP_solver_finalize(TP_solver self)
{
  plasma_finalize();
}


void
TP_solver_destroy(TP_solver self)
{
  TP_matrix_destroy(self->A);
  TP_schur_matrix_destroy(self->S);
  TP_matrix_destroy(self->L);
  TP_matrix_destroy(self->U);
  TP_matrix_destroy(self->D);

  if(self->S_dense)
    TP_dense_matrix_destroy(self->S_dense);

  DD_verbose_destroy(self->verbose)

  free(self->row_perm);
  free(self->col_perm);
  free(self->random_col);
  free(self);
}
