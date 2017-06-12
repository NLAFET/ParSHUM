#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h>
#include <plasma.h>

#include "TP_verbose.h"
#include "TP_enum.h" 
#include "TP_matrix.h" 
#include "TP_dense.h" 
#include "TP_schur_matrix.h" 
#include "TP_pivot_list.h" 
#include "TP_auxiliary.h"
#include "TP_solver.h" 
#include <limits.h>

const char *usageStrign[] = {
  "usage test: [-h] [--matrix matrix] [--debug_mode] [--verbosity level] [--marko_tol tol] [--value_tol tol]",
  "            [--extra_space factor] [--extra_space_inbetween factor] [--nb_threads #threads] [--nb_candidates_per_block #blocks] ", 
  "            [--ouput_dir dir] [--output_file file] [--nb_previous_pivots #pivtos] [--schur_density_tolerance tol]",
  "            [--min_pivot_per_steps #steps] [--output_dir dir] [--prog_name name ] [--debug_mode]  ",
  NULL,
};

int is_plasma_init;

TP_solver 
TP_solver_create()
{
  TP_solver self = calloc(1, sizeof(*self));
  self->exe_parms = calloc(1, sizeof(*self->exe_parms));

  self->verbosity               = 1;

  self->exe_parms->nb_threads              = 1;
  self->exe_parms->value_tol               = 0.1;
  self->exe_parms->marko_tol               = 4;
  self->exe_parms->extra_space             = 1.0;
  self->exe_parms->extra_space_inbetween   = 1.0;
  self->exe_parms->nb_candidates_per_block = 10;
  self->exe_parms->nb_previous_pivots      = 5;
  self->exe_parms->min_pivot_per_steps     = 20;
  self->exe_parms->density_tolerance       = 0.2;

  self->verbose = TP_verbose_create(self->exe_parms);
  return self;
}

void
TP_solver_dealloc(TP_solver self)
{
  TP_verbose_destroy(self->verbose);
  
  free(self->exe_parms);
  free(self);
}


int
check_TP_with_plasma_perm(int argc, char **argv)
{
  TP_solver plasma;
  TP_vector X, rhs, sol_plasma, sol_TP;
  int i;
  
  plasma = TP_solver_create();
  TP_solver_parse_args(plasma, argc, argv);
  plasma->exe_parms->density_tolerance   = 1.0;
  plasma->exe_parms->min_pivot_per_steps = 5;
  plasma->exe_parms->nb_previous_pivots  = 5;
  plasma->debug |= TP_CHECK_TP_W_PLASMA_PERM;
  TP_solver_read_matrix(plasma);
  TP_solver_init(plasma);
  
  X   = TP_vector_create(plasma->A->n);
  rhs = TP_vector_create(plasma->A->n);
  sol_plasma = TP_vector_create(plasma->A->n);
  sol_TP     = TP_vector_create(plasma->A->n);
  
  for(i=0; i < plasma->A->n; i++)
    X->vect[i] = i * 10;
  TP_vector_memset(X, 1.0);
  TP_vector_memset(rhs, 0.0);
  TP_matrix_SpMV(plasma->A, X, rhs);
  TP_vector_copy(rhs, sol_plasma);

  TP_solver_factorize(plasma);
  TP_solver_solve(plasma, sol_plasma);
  TP_solver_copmpute_norms(plasma, X, sol_plasma, rhs);
  TP_solver_iterative_refinement(plasma, X, rhs, 0.0);

  TP_solver_finalize(plasma);

  /* apply plasma row permutation to A */
  int *plasma_perms = TP_dense_get_row_perms(plasma->S_dense);
  TP_matrix debug_matrix = TP_matrix_permute(plasma->A, plasma->S_dense->original_cols, plasma_perms);
  TP_solver debug_solver = TP_solver_create();

  debug_solver->A = debug_matrix;
  debug_solver->debug |= TP_CHECK_TP_W_PLASMA_PERM;
  debug_solver->exe_parms->density_tolerance = 1.0;
  debug_solver->exe_parms->min_pivot_per_steps = 5;
  debug_solver->exe_parms->nb_previous_pivots  = 5;
  TP_solver_init(debug_solver);

  TP_vector_permute(rhs, plasma_perms);
  TP_vector_copy(rhs, sol_TP);
  TP_solver_factorize(debug_solver);
  TP_solver_solve(debug_solver, sol_TP);

  TP_solver_copmpute_norms(debug_solver, X, sol_TP, rhs);
  TP_solver_iterative_refinement(debug_solver, X, rhs, 0.0);
  
  TP_solver_finalize(debug_solver);

  free(plasma_perms);
  TP_vector_destroy(X);
  TP_vector_destroy(rhs);
  TP_vector_destroy(sol_plasma);
  TP_vector_destroy(sol_TP);
  TP_solver_destroy(debug_solver);
  TP_solver_destroy(plasma);

  return 0;
}


int
check_dense_with_TP_perm(int argc, char **argv)
{
  TP_solver self;
  TP_vector X, rhs, sol_TP, sol_dense;
  int i;
  
  self = TP_solver_create();
  TP_solver_parse_args(self, argc, argv);
  self->exe_parms->density_tolerance = 1.0;
  self->exe_parms->min_pivot_per_steps = 5;
  self->exe_parms->nb_previous_pivots  = 5;
  TP_solver_read_matrix(self);
  TP_solver_init(self);
  
  X   = TP_vector_create(self->A->n);
  rhs = TP_vector_create(self->A->n);
  sol_TP    = TP_vector_create(self->A->n);
  sol_dense = TP_vector_create(self->A->n);
  
  for(i=0; i < self->A->n; i++)
    X->vect[i] = i * 10;
  TP_vector_memset(rhs, 0.0);
  TP_matrix_SpMV(self->A, X, rhs);
  TP_vector_copy(rhs, sol_TP);

  TP_solver_factorize(self);
  TP_solver_solve(self, sol_TP);
  TP_solver_copmpute_norms(self, X, sol_TP, rhs);
  TP_solver_iterative_refinement(self, X, rhs, 0.0);

  TP_solver_finalize(self);
 
  TP_solver dense_solver = TP_solver_create();
  TP_solver_parse_args(dense_solver, argc, argv);

  dense_solver->A = self->A;

  dense_solver->debug |= TP_CHECK_DENSE_W_TP_PERM;
  dense_solver->row_perm = self->row_perm;
  dense_solver->col_perm = self->col_perm;
  dense_solver->invr_row_perm = self->invr_row_perm;
  dense_solver->invr_col_perm = self->invr_col_perm;
  dense_solver->verbose->parms->prog_name = "Dense solver";

  TP_solver_init(dense_solver);

  TP_vector_copy(rhs, sol_dense);
  TP_solver_factorize(dense_solver);

  TP_solver_solve(dense_solver, sol_dense);

  TP_solver_copmpute_norms(dense_solver, X, sol_dense, rhs);
  TP_solver_iterative_refinement(dense_solver, X, rhs, 0.0);

  TP_solver_finalize(dense_solver);
 
  TP_vector_destroy(X);
  TP_vector_destroy(rhs);
  TP_vector_destroy(sol_TP);
  TP_vector_destroy(sol_dense);
  TP_solver_destroy(dense_solver);
  TP_solver_destroy(self);

  return 0;
}

void 
TP_solver_parse_args(TP_solver self, int argc, char **argv)
{
  int i;
  for (i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--verbosity")) {
      self->verbose->parms->verbosity = atoi(argv[++i]);
      self->verbosity = atoi(argv[i]);
      continue;
    } else if (!strcmp(argv[i], "--matrix")) {
      self->exe_parms->matrix_file = argv[++i];
      continue;
    } else if (!strcmp(argv[i], "--value_tol")) {
      double tmp = atof(argv[++i]);
      if ( tmp > 1.0 || tmp < 0.0) 
	TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "value tolerance should be between 0 and 1");
      self->exe_parms->value_tol = tmp;
      continue;
    } else if (!strcmp(argv[i], "--marko_tol")) {
      double tmp = atof(argv[++i]);
      if ( tmp < 1.0 ) 
	TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "marko tolerance should be larger then 1");
      self->exe_parms->marko_tol = tmp;
      continue;
    } else if (!strcmp(argv[i], "--nb_threads")) {
      int tmp = atoi(argv[++i]);
      self->exe_parms->nb_threads = tmp;
      continue;
    } else if (!strcmp(argv[i], "--extra_space")) {
      double tmp = atof(argv[++i]);
      self->exe_parms->extra_space = tmp;
      continue;
    } else if (!strcmp(argv[i], "--extra_space_inbetween")) {
      double tmp = atof(argv[++i]);
      self->exe_parms->extra_space_inbetween = tmp;
      continue;
    } else if (!strcmp(argv[i], "--nb_candidates_per_block")) {
      int tmp = atoi(argv[++i]);
      if (tmp < 1)
	TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "nb_candidates_per_blocks should be at least 1");
      self->exe_parms->nb_candidates_per_block = tmp;
      continue;
    }  else if (!strcmp(argv[i], "--prog_name")) {
      self->verbose->parms->prog_name = argv[++i];
      continue;
    } else if (!strcmp(argv[i], "--output_dir")) {
      self->verbose->parms->output_dir      = argv[++i];
      self->verbose->parms->user_out_dir = 1;
      continue;
    } else if (!strcmp(argv[i], "--output_file")) {
      FILE *file  = fopen(argv[++i], "w+");
      if ( file ) { 
	self->verbose->parms->out_file = file;  
	self->verbose->parms->user_out_file = 1;
      } else {
	TP_warning(__FUNCTION__, __FILE__, __LINE__,"unable to open output file, the program will write on stdout instead");	  
	file = stdout;
      }
      continue;
    } else if (!strcmp(argv[i], "--nb_previous_pivots")) {
      int tmp = atoi(argv[++i]);
      if (tmp < 1)
	TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "nb_previous_pivots should be at least 1");
      self->exe_parms->nb_previous_pivots  = tmp;
      continue;
    } else if (!strcmp(argv[i], "--schur_density_tolerance")) {
      double tmp = atof(argv[++i]);
      if ( tmp > 1.0 ) 
	TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "schur density tolerance can not be larger then 1");
      if ( tmp < 0.0 ) 
	TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "schur density tolerance can not be smaller then 0");
      self->exe_parms->density_tolerance = tmp;
      continue;
    } else if (!strcmp(argv[i], "--min_pivot_per_steps")) {
      int tmp = atof(argv[++i]);
      if ( tmp < self->exe_parms->nb_previous_pivots ) 
	TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "min_pivot_per_steps should be at least nb_previous_pivots");
      self->exe_parms->min_pivot_per_steps = tmp;
      continue;
    } else if (!strcmp(argv[i], "--debug_mode")) {
      TP_warning(__FUNCTION__, __FILE__, __LINE__,"debug mode is not implemented");	  
      continue;
    }  else if (!strcmp(argv[i], "--check_pivots")) {
      self->debug |= TP_CHECK_PIVOTS;
      continue;
    }  else if (!strcmp(argv[i], "--check_schur_memory")) {
      self->debug |= TP_CHECK_SCHUR_MEMORY;
      continue;
    }  else if (!strcmp(argv[i], "--check_schur_symetry")) {
      self->debug |= TP_CHECK_SCHUR_SYMETRY;
      continue;
    }  else if (!strcmp(argv[i], "--check_TP_with_plasma_perm")) {
      TP_solver_dealloc(self);
      exit(check_TP_with_plasma_perm(argc, argv));
    }  else if (!strcmp(argv[i], "--check_dense_with_TP_perm")) {
      TP_solver_dealloc(self);
      exit(check_dense_with_TP_perm(argc, argv));
      continue;
    } else if (!strcmp(argv[i], "--help")) {
      int i = 0;
      while( usageStrign[i] !=  NULL)
	printf("%s\n", usageStrign[i++]);
      exit(0);
    } else {
      char mess[2048];
      snprintf(mess, 2048, "unrecognized option \"%s\" ", argv[i]);
      while( usageStrign[i] !=  NULL)
	printf("%s\n", usageStrign[i++]);
      TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, mess);
    }
  }
}


char *
get_outfile_prefix(TP_exe_parms exe_parms)
{
  char *self = malloc(2048*sizeof(*self));
  snprintf(self, 2048,"%s_%fValTol_%fMarkoTol_%dthreads_%dcandidates_%fdensityTol_%dminPivots",
	   basename(exe_parms->matrix_file), exe_parms->value_tol, exe_parms->marko_tol, 
	   exe_parms->nb_threads, exe_parms->nb_candidates_per_block, exe_parms->density_tolerance, 
	   exe_parms->min_pivot_per_steps);
  return self;
}


void
TP_solver_alloc_internal(TP_solver self) 
{
  double total_extra_space = 1 + self->exe_parms->extra_space_inbetween + self->exe_parms->extra_space;
  int needed_pivots = self->A->n < self->A->m ? self->A->n : self->A->m;

  self->S    = TP_schur_matrix_create();
  self->L    = TP_matrix_create();
  self->U    = TP_matrix_create();
  self->D    = TP_matrix_create();

  if (self->debug & (TP_DEBUG_VERBOSE_EACH_STEP | TP_DEBUG_GOSSIP_GIRL))
      TP_matrix_print(self->A, "A on input");
  TP_schur_matrix_allocate(self->S, self->A->n, self->A->m, self->A->nnz, self->debug,
			   self->exe_parms->extra_space, self->exe_parms->extra_space_inbetween);
  TP_schur_matrix_copy(self->A, self->S);

  if (self->debug & (TP_DEBUG_VERBOSE_EACH_STEP | TP_DEBUG_GOSSIP_GIRL))
      TP_schur_matrix_print(self->S, "S on input");
  TP_matrix_allocate(self->L, self->A->n, 0, (self->A->nnz - self->A->n) / 2, total_extra_space, TP_CSC_matrix);
  TP_matrix_allocate(self->U, 0, self->A->m, (self->A->nnz - self->A->m) / 2, total_extra_space, TP_CSR_matrix);
  TP_matrix_allocate(self->D, self->A->n, 0, 0, 1.0, TP_Diag_matrix);
  
  self->row_perm      = malloc((size_t) needed_pivots * sizeof(*self->invr_row_perm));
  self->col_perm      = malloc((size_t) needed_pivots * sizeof(*self->invr_col_perm));
  self->invr_row_perm = malloc((size_t) needed_pivots * sizeof(*self->invr_row_perm));
  self->invr_col_perm = malloc((size_t) needed_pivots * sizeof(*self->invr_col_perm));

  int_array_memset(self->col_perm,      TP_UNUSED_PIVOT, needed_pivots);
  int_array_memset(self->row_perm,      TP_UNUSED_PIVOT, needed_pivots);
  int_array_memset(self->invr_col_perm, TP_UNUSED_PIVOT, needed_pivots);
  int_array_memset(self->invr_row_perm, TP_UNUSED_PIVOT, needed_pivots);

  self->random_col = create_randomize(self->A->n);

  self->previous_pivots = malloc(self->exe_parms->nb_previous_pivots *  sizeof(*self->previous_pivots));
  int_array_memset(self->previous_pivots, INT_MAX / self->exe_parms->nb_previous_pivots, self->exe_parms->nb_previous_pivots);

  self->candidates = calloc(1, sizeof(*self->candidates));
  self->candidates->row        = malloc(self->A->n * sizeof(*self->candidates->row));
  self->candidates->marko      = malloc(self->A->n * sizeof(*self->candidates->marko));
  self->candidates->best_marko = malloc(self->exe_parms->nb_threads * sizeof(*self->candidates->best_marko));
}

void
TP_solver_init(TP_solver self)
{
  self->verbose->n = self->A->n;
  self->verbose->m = self->A->m;
  self->verbose->nnz_input = self->verbose->nnz_final = self->A->nnz;

  self->verbose->parms->outfiles_prefix = get_outfile_prefix(self->exe_parms);
  
  if (!is_plasma_init) {
    plasma_init(self->exe_parms->nb_threads);
    is_plasma_init = 1;
  }
  
  if ( ! (self->debug & TP_CHECK_DENSE_W_TP_PERM) )
       TP_solver_alloc_internal(self);
}

void
TP_solver_read_matrix(TP_solver self)
{
  char *file_ext = strrchr(self->exe_parms->matrix_file, '.');
  self->A        = TP_matrix_create();
  
  if (!strcmp(file_ext, ".rb")) 
    read_rutherford_boeing(self->A, self->exe_parms->matrix_file);
  else if (!strcmp(file_ext, ".mtl"))
    TP_read_mtl_file(self->A, self->exe_parms->matrix_file);
  else 
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"the input matrix is not sorted column-wise ");
}

void 
TP_solver_get_pivots(TP_solver self, TP_pivot_set set)
{
  TP_pivot_cell cells = set->cells;
  int new_pivots = 0, old_pivots = self->found_pivots;
  int *row_perms = self->row_perm;
  int *col_perms = self->col_perm;

  int *invr_row_perms = self->invr_row_perm;
  int *invr_col_perms = self->invr_col_perm;

  while (cells)
    {
      int current_pivot = new_pivots + old_pivots;
      row_perms[current_pivot] = cells->row;
      col_perms[current_pivot] = cells->col;

      invr_row_perms[cells->row] = current_pivot;
      invr_col_perms[cells->col] = current_pivot;

      new_pivots++;
      cells = cells->next;
    }

  self->found_pivots += new_pivots;
  TP_verbose_update_pivots(self->verbose, new_pivots);
}

void
TP_solver_find_pivot_set(TP_solver self)
{
  TP_pivot_list list;
  TP_exe_parms exe_parms = self->exe_parms;
  TP_verbose_per_step step = TP_verbose_get_step(self->verbose);
  int nb_threads = self->exe_parms->nb_threads;

  TP_verbose_start_timing(&step->timing_extracting_candidates);
  list = get_possible_pivots(self->S, self->random_col, self->candidates, nb_threads,
  			     exe_parms->value_tol, exe_parms->marko_tol,
  			     exe_parms->nb_candidates_per_block);
  TP_verbose_stop_timing(&step->timing_extracting_candidates);
  
  if( !list->nb_elem )
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "no possible pivot were found");

  TP_verbose_start_timing(&step->timing_merging_pivots);
  list = merge_pivot_sets(list, self->S);

  TP_solver_get_pivots(self, list->sets);

  if (self->debug & (TP_DEBUG_VERBOSE_EACH_STEP | TP_DEBUG_GOSSIP_GIRL))
    print_pivot_list(list, "found pivots");
  TP_pivot_list_destroy(list);

  TP_verbose_stop_timing(&step->timing_merging_pivots);
}

void
TP_solver_update_matrix(TP_solver self)
{
  TP_matrix L = self->L;
  TP_matrix D = self->D;
  TP_matrix U = self->U;
  TP_schur_matrix S = self->S;
  int nb_pivots = self->found_pivots - self->done_pivots;
  TP_verbose_per_step step = TP_verbose_get_step(self->verbose);

  TP_verbose_start_timing(&step->timing_update_LD);
  TP_schur_matrix_update_LD(S, L, D, &self->row_perm[self->done_pivots], &self->col_perm[self->done_pivots], nb_pivots);
  TP_verbose_stop_timing(&step->timing_update_LD);
  if (self->debug & TP_CHECK_SCHUR_MEMORY )
    TP_schur_matrix_memory_check(self->S);
  if (self->debug & TP_CHECK_SCHUR_SYMETRY )
    TP_schur_matrix_check_symetry(self->S);
  
  TP_verbose_start_timing(&step->timing_update_U);
  TP_schur_matrix_update_U (S, U,    &self->row_perm[self->done_pivots], &self->col_perm[self->done_pivots], nb_pivots);
  TP_verbose_stop_timing(&step->timing_update_U);
  if (self->debug & TP_CHECK_SCHUR_MEMORY )
    TP_schur_matrix_memory_check(self->S);
  if (self->debug & TP_CHECK_SCHUR_SYMETRY )
    TP_schur_matrix_check_symetry(self->S);
  
  TP_verbose_start_timing(&step->timing_update_S);
  TP_schur_matrix_update_S (S, L, U, self->done_pivots, self->found_pivots);
  TP_verbose_stop_timing(&step->timing_update_S);
  if (self->debug & TP_CHECK_SCHUR_MEMORY )
    TP_schur_matrix_memory_check(self->S);
  if (self->debug & TP_CHECK_SCHUR_SYMETRY )
    TP_schur_matrix_check_symetry(self->S);


  self->done_pivots = self->found_pivots;

  if (self->debug & TP_CHECK_PIVOTS )
    TP_schur_matrix_check_pivots(self->S,
				 self->row_perm,
				 self->col_perm,
				 self->invr_row_perm,
				 self->invr_col_perm,
				 self->done_pivots);

  if (self->debug & (TP_DEBUG_VERBOSE_EACH_STEP | TP_DEBUG_GOSSIP_GIRL)) {
    TP_schur_matrix_print(S, "S after update");
    TP_matrix_print(L, "L after update");
    TP_matrix_print(D, "D after update");
    TP_matrix_print(U, "U after update");
  }
}

int
TP_continue_pivot_search(TP_schur_matrix S,
			 int nb_done_pivots,
			 int nb_needed_pivots, 
			 int *previous_pivots,
			 int nb_previous_pivots, 
			 double density_tolerance,
			 int min_pivot_per_steps,
			 int debug,
			 TP_verbose verbose)
{
  long n_schur, m_schur, i, sum_pivots ; 

  if ( debug & TP_CHECK_TP_W_PLASMA_PERM ||
       debug & TP_CHECK_DENSE_W_TP_PERM ) {
    verbose->reason = TP_because;
    return 0;
  }

  n_schur = (long) S->n - nb_done_pivots;
  m_schur = (long) S->m - nb_done_pivots;
  if ( S->nnz > n_schur * m_schur * density_tolerance) {
    verbose->reason = TP_density;
    return 0;
  }

  for( i = 0, sum_pivots = 0; i < nb_previous_pivots; i++)
    sum_pivots += previous_pivots[i];
  if (sum_pivots < min_pivot_per_steps) {
    verbose->reason = TP_no_pivots;
    return 0;
  }

  if (nb_done_pivots >= nb_needed_pivots) {
    verbose->reason = TP_no_pivots;
    return 0;
  }

  return 1;
}

void
TP_solver_factorize(TP_solver self)
{
  TP_exe_parms exe_parms = self->exe_parms;
  TP_verbose verbose = self->verbose;
  int needed_pivots = self->A->n < self->A->m ? self->A->n : self->A->m;
  int *previous_pivots   = self->previous_pivots;
  int nb_previous_pivots = exe_parms->nb_previous_pivots;
  int nb_pivot_blocks = 0;
  
  TP_verbose_start_timing(&verbose->timing_facto);
  TP_verbose_start_timing(&verbose->timing_facto_sparse);
  while ( TP_continue_pivot_search(self->S, self->done_pivots, needed_pivots, 
				   previous_pivots, nb_previous_pivots,
				   exe_parms->density_tolerance,
				   exe_parms->min_pivot_per_steps,
				   self->debug, verbose) )
    { 
       TP_verbose_per_step  step = TP_verbose_step_start(verbose);
       TP_verbose_start_timing(&step->timing_step);
       TP_verbose_start_timing(&step->timing_pivot_search);
       TP_solver_find_pivot_set(self);
       previous_pivots[nb_pivot_blocks++ % nb_previous_pivots] = self->found_pivots - self->done_pivots;
       TP_verbose_stop_timing(&step->timing_pivot_search);
       
       TP_verbose_start_timing(&step->timing_apply_perms);
       TP_solver_update_matrix(self);
       TP_verbose_stop_timing(&step->timing_apply_perms);
       TP_verbose_stop_timing(&step->timing_step);
    }
  TP_verbose_stop_timing(&verbose->timing_facto_sparse);
  
  TP_verbose_start_timing(&verbose->timing_convert_schur);
  if ( self->debug & TP_CHECK_DENSE_W_TP_PERM )  {
    self->A_debug = TP_dense_2D_permute(TP_dense_2D_convert_sparse(self->A), 
					self->row_perm, self->col_perm);
    verbose->schur_density = 1.00;
  }  else {
    self->S_dense = TP_schur_matrix_convert(self->S, self->done_pivots);
    verbose->schur_density = 
      (double) self->S->nnz / ((self->S->n - self->done_pivots) * (self->S->m - self->done_pivots));

    if (self->debug & (TP_DEBUG_VERBOSE_EACH_STEP | TP_DEBUG_GOSSIP_GIRL)) 
      TP_dense_matrix_print(self->S_dense, "dense schur after conversion");
  }
  TP_verbose_stop_timing(&verbose->timing_convert_schur);

  TP_verbose_start_timing(&verbose->timing_facto_dense);
  if (self->debug & TP_CHECK_DENSE_W_TP_PERM)  {
    TP_dense_2D_facto(self->A_debug);
  } else {
    TP_dense_matrix_factorize(self->S_dense);

    // Handeling the pivots
    memcpy(&self->row_perm[self->done_pivots], self->S_dense->original_rows,
	   (self->A->m - self->done_pivots) * sizeof(*self->row_perm));
    memcpy(&self->col_perm[self->done_pivots], self->S_dense->original_cols,
	   (self->A->n - self->done_pivots) * sizeof(*self->col_perm));
    int i;
    for(i = self->done_pivots; i < needed_pivots; i++) {
      self->invr_row_perm[self->row_perm[i]] = i;
      self->invr_col_perm[self->col_perm[i]] = i;
    }

    TP_verbose_update_dense_pivots(verbose, needed_pivots - self->done_pivots);
    verbose->nnz_final   = self->L->nnz + self->U->nnz + self->D->n + self->S_dense->n * self->S_dense->m;
    verbose->nnz_L       = self->L->nnz  + ( ( self->S_dense->n * self->S_dense->m - self->S_dense->m) / 2 )  ;
    verbose->nnz_U       = self->U->nnz + self->D->n + ( ( self->S_dense->n * self->S_dense->m - self->S_dense->m) / 2    + self->S_dense->m ) ;
    verbose->nnz_S_dense = self->S_dense->n * self->S_dense->m;
  }
  TP_verbose_stop_timing(&verbose->timing_facto_dense);
  TP_verbose_stop_timing(&verbose->timing_facto);
}


void
TP_solver_solve(TP_solver self, TP_vector RHS)
{
  double *RHS_vals          = RHS->vect;
  TP_verbose verbose        = self->verbose;

  TP_verbose_start_timing(&verbose->timing_solve); 
  if (self->debug & TP_CHECK_DENSE_W_TP_PERM ) { 
    TP_vector_permute(RHS, self->row_perm);
    TP_dense_2D_solve(self->A_debug, RHS);
    TP_vector_permute(RHS, self->invr_col_perm);
  } else {
    TP_verbose_start_timing(&verbose->timing_solve_L);
    if (self->debug & (TP_DEBUG_VERBOSE_EACH_STEP | TP_DEBUG_GOSSIP_GIRL)) 
      TP_vector_print(RHS, "initial RHS");
    TP_vector_permute(RHS, self->row_perm);
    if (self->debug & TP_DEBUG_GOSSIP_GIRL) 
      TP_vector_print(RHS, "P RHS");
    TP_matrix_solve_L(self->L, RHS, self->invr_row_perm);
    if (self->debug & TP_DEBUG_GOSSIP_GIRL) 
      TP_vector_print(RHS, "after forward solve");
    TP_verbose_stop_timing(&verbose->timing_solve_L);
    
    TP_verbose_start_timing(&verbose->timing_solve_dense);
    if (self->S_dense->n)
      plasma_dgetrs(self->S_dense->n, 1, self->S_dense->val, self->S_dense->n,
    		    self->S_dense->pivots, &RHS_vals[self->done_pivots], self->S_dense->n);
    if (self->debug & TP_DEBUG_GOSSIP_GIRL) 
      TP_vector_print(RHS, "after dense solve");
    TP_verbose_stop_timing(&verbose->timing_solve_dense);
    
    TP_verbose_start_timing(&verbose->timing_solve_U);
    TP_matrix_solve_UD(self->U, self->D, RHS, self->invr_col_perm);
    if (self->debug & TP_DEBUG_GOSSIP_GIRL) 
      TP_vector_print(RHS, "after backward solve");
    TP_vector_permute(RHS, self->invr_col_perm);
    if (self->debug & (TP_DEBUG_VERBOSE_EACH_STEP | TP_DEBUG_GOSSIP_GIRL)) 
      TP_vector_print(RHS, "after solve operation");
    TP_verbose_stop_timing(&verbose->timing_solve_U);
  }

  TP_verbose_stop_timing(&verbose->timing_solve);
}

void
TP_solver_copmpute_norms(TP_solver self,       TP_vector X,
			 TP_vector X_computed, TP_vector rhs)
{
  double x_norm, A_norm, b_norm;
  TP_vector tmp = TP_vector_create(X->n);

  TP_vector_add(X, 1.00, X_computed, -1.00, tmp);
  self->verbose->forward_error = TP_vector_2norm(tmp) / TP_vector_2norm(X);
  
  /* || Ax - b || */
  TP_matrix_SpMV(self->A, X_computed, tmp);
  TP_vector_add(tmp, 1.00, rhs, -1.00, tmp);
  self->verbose->backward_error = TP_vector_2norm(tmp);

  A_norm  = TP_matrix_get_norm(self->A);
  x_norm  = TP_vector_2norm(X_computed);
  b_norm  = TP_vector_2norm(rhs);

  self->verbose->backward_error /= A_norm * x_norm + b_norm;

  TP_vector_destroy(tmp);
}


void
TP_solver_iterative_refinement(TP_solver self, 
			       TP_vector X, 
			       TP_vector rhs,
			       double wanted_precision)
{
  TP_vector sol = TP_vector_create(self->A->n);
  TP_vector tmp = TP_vector_create(self->A->n);
  int i = 0;

  TP_vector_copy(rhs, sol);
  TP_solver_solve(self, sol);
  TP_solver_copmpute_norms(self, X, sol, rhs);
  printf("iteration %d: backward error = %e and forward error = %e\n", 
	 i++,
	 self->verbose->backward_error,
	 self->verbose->forward_error);

  while ( i < 20 &&  self->verbose->backward_error > wanted_precision) {
    TP_matrix_SpMV(self->A, sol, tmp);
    TP_vector_add(rhs, 1.00, tmp, -1.00, tmp);
    TP_solver_solve(self, tmp);
    TP_vector_add(sol, 1.00, tmp, 1.00, sol);
	 
    TP_solver_copmpute_norms(self, X, sol, rhs);
    printf("iteration %d: backward error = %e and forward error = %e\n", 
	   i++,
	   self->verbose->backward_error,
	   self->verbose->forward_error);
  }

  TP_vector_destroy(sol);
  TP_vector_destroy(tmp);
}


void 
TP_solver_finalize(TP_solver self)
{
  if (is_plasma_init) {
    plasma_finalize();
    is_plasma_init = 0;
  }
  TP_verbose_print(self->verbose);
  TP_verbose_create_dirs(self->verbose->parms->output_dir);
  TP_verbose_draw_graph(self->verbose);
}

void
TP_solver_destroy(TP_solver self)
{
  if (self->debug & TP_CHECK_DENSE_W_TP_PERM) {
    TP_dense_2D_destroy(self->A_debug);
  } else   {
    TP_matrix_destroy(self->A);
    TP_schur_matrix_destroy(self->S);
    TP_matrix_destroy(self->L);
    TP_matrix_destroy(self->U);
    TP_matrix_destroy(self->D);
    
    if(self->S_dense)
      TP_dense_matrix_destroy(self->S_dense);

    free(self->candidates->row);
    free(self->candidates->marko);
    free(self->candidates->best_marko);
    free(self->candidates);
    
    free(self->row_perm);
    free(self->col_perm);
    free(self->invr_row_perm);
    free(self->invr_col_perm);
    free(self->random_col);
    free(self->previous_pivots);
  }
  
  TP_solver_dealloc(self);
}

