#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h>
#include <plasma.h>
#include <math.h>
#include <limits.h>

#include "TP_verbose.h"
#include "TP_enum.h" 
#include "TP_matrix.h" 
#include "TP_dense.h" 
#include "TP_schur_matrix.h" 
#include "TP_pivot_list.h" 
#include "TP_paje.h"
#include "TP_auxiliary.h"
#include "TP_solver.h" 

const char *usageStrign[] = {
  "usage test: [--help] [--matrix matrix] [--RHS_file file] [--debug_mode] [--verbosity level] [--marko_tol tol] [--value_tol tol]",
  "            [--extra_space factor] [--extra_space_inbetween factor] [--nb_threads #threads] [--nb_candidates_per_block #blocks] ", 
  "            [--output_dir dir] [--output_file file] [--nb_previous_pivots #pivtos] [--schur_density_tolerance tol]",
  "            [--min_pivot_per_steps #steps] [--prog_name name ] [--check_schur_symetry] [--check_schur_memory] [--check_pivots]",
  "            [--check_TP_with_plasma_perm] [--check_dense_with_TP_perm] [--print_each_step] [--check_GC]",
  "            [--group_run value_tol|marko_tol|schur_density|nb_candidates|min_pivots|nb_threads init inc nb_steps]",
  "            [--counters_size #double_counters] [--check_counters] [--check_schur_doubles] [--max_dense_schur size]",
  "            [--luby_algorithm] [--trace]", 
  NULL,
};

int is_plasma_init;

int  TP_solver_run_group(TP_solver solver, TP_parm_type type, 
			 void *init_val, int nb_steps, void *inc);

TP_solver 
TP_solver_create()
{
  TP_solver self = calloc(1, sizeof(*self));
  self->exe_parms = calloc(1, sizeof(*self->exe_parms));

  self->verbosity      = 1;
  self->size_counters = 100;

  self->exe_parms->nb_threads              = 1;
  self->exe_parms->value_tol               = 0.1;
  self->exe_parms->marko_tol               = 4;
  self->exe_parms->extra_space             = 1.0;
  self->exe_parms->extra_space_inbetween   = 1.0;
  self->exe_parms->nb_candidates_per_block = 10;
  self->exe_parms->nb_previous_pivots      = 5;
  self->exe_parms->min_pivot_per_steps     = 20;
  self->exe_parms->density_tolerance       = 0.2;
  self->exe_parms->max_dense_schur         = 10000;
  
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
  
  if (plasma->exe_parms->RHS_file)
    TP_vector_read_file(X, plasma->exe_parms->RHS_file);
  else 
    for(i=0; i < plasma->A->n; i++)
      X->vect[i] = (i+1) * 10;
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
  
  if (self->exe_parms->RHS_file)
    TP_vector_read_file(X, self->exe_parms->RHS_file);
  else 
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
  int i, run_args_start = 0;

  for (i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--verbosity")) {
      self->verbose->parms->verbosity = atoi(argv[++i]);
      self->verbosity = atoi(argv[i]);
      continue;
    } else if (!strcmp(argv[i], "--matrix")) {
      self->exe_parms->matrix_file = argv[++i];
      continue;
    }  else if (!strcmp(argv[i], "--RHS_file")) {
      self->exe_parms->RHS_file = argv[++i];
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
      int tmp = atoi(argv[++i]);
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
    }  else if (!strcmp(argv[i], "--check_counters")) {
      self->debug |= TP_CHECK_COUNTERS;
      continue;
    } else if (!strcmp(argv[i], "--print_each_step")) {
      self->debug |= TP_DEBUG_VERBOSE_EACH_STEP;
      continue;
    } else if (!strcmp(argv[i], "--verbose_gossip_girl")) {
      self->debug |= TP_DEBUG_GOSSIP_GIRL;
      continue;
    } else if (!strcmp(argv[i], "--check_schur_doubles")) {
      self->debug |= TP_CHECK_SCHUR_DOUBLES;
      continue;
    }  else if (!strcmp(argv[i], "--check_TP_with_plasma_perm")) {
      TP_solver_dealloc(self);
      exit(check_TP_with_plasma_perm(argc, argv));
    }  else if (!strcmp(argv[i], "--check_dense_with_TP_perm")) {
      TP_solver_dealloc(self);
      exit(check_dense_with_TP_perm(argc, argv));
      continue;
    }  else if (!strcmp(argv[i], "--group_run")) {
      run_args_start = ++i;
      i += 3;
      continue;
    } else if (!strcmp(argv[i], "--counters_size")) {
      self->size_counters = atoi(argv[++i]);
      continue;
    } else if (!strcmp(argv[i], "--max_dense_schur")) {
      int tmp = atoi(argv[++i]);
      if (tmp < 1)
	TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "max_dense_schur should be at least 1");
      self->exe_parms->max_dense_schur = tmp;
      continue;
    } else if (!strcmp(argv[i], "--luby_algorithm")) {
      self->exe_parms->luby_algo = 1;
    } else if (!strcmp(argv[i], "--trace")) {
      self->exe_parms->trace = 1;
    } else if (!strcmp(argv[i], "--help")) {
      int j = 0;
      while( usageStrign[j] !=  NULL)
	printf("%s\n", usageStrign[j++]);
      exit(0);
    } else {
      char mess[2048];
      snprintf(mess, 2048, "unrecognized option \"%s\" ", argv[i]);
      int j = 0;
      while( usageStrign[j] !=  NULL)
	printf("%s\n", usageStrign[j++]);
      TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, mess);
    }
  }
  
  if (run_args_start) {
    TP_parm_type type;
    if ( !strcmp(argv[run_args_start], "value_tol") ) {
      type = TP_value_tol;
      double init = atof(argv[++run_args_start]), inc = atof(argv[++run_args_start]);
      int nb_steps = atoi(argv[++run_args_start]);
      exit(TP_solver_run_group(self, type, (void *) &init, nb_steps, (void *) &inc)); 
    } else if ( !strcmp(argv[run_args_start], "marko_tol") ) {
      type = TP_marko_tol;
      double init = atof(argv[++run_args_start]), inc = atof(argv[++run_args_start]);
      int nb_steps = atoi(argv[++run_args_start]);
      exit(TP_solver_run_group(self, type, (void *) &init, nb_steps, (void *) &inc)); 
    } else if ( !strcmp(argv[run_args_start], "schur_density") ) {
      type = TP_schur_density;
      double init = atof(argv[++run_args_start]), inc = atof(argv[++run_args_start]);
      int nb_steps = atoi(argv[++run_args_start]);
      exit(TP_solver_run_group(self, type, (void *) &init, nb_steps, (void *) &inc)); 
    } else if ( !strcmp(argv[run_args_start], "nb_candidates") ) {
      type = TP_nb_candidates;
      int init = atoi(argv[++run_args_start]), inc = atoi(argv[++run_args_start]);
      int nb_steps = atoi(argv[++run_args_start]);
      exit(TP_solver_run_group(self, type, (void *) &init, nb_steps, (void *) &inc)); 
    } else if ( !strcmp(argv[run_args_start], "min_pivots") ) {
      type = TP_min_pivots;
      int init = atoi(argv[++run_args_start]), inc = atoi(argv[++run_args_start]);
      int nb_steps = atoi(argv[++run_args_start]);
      exit(TP_solver_run_group(self, type, (void *) &init, nb_steps, (void *) &inc)); 
    } else if ( !strcmp(argv[run_args_start], "nb_threads") ) {
      type = TP_nb_threads;
      int init = atoi(argv[++run_args_start]), inc = atoi(argv[++run_args_start]);
      int nb_steps = atoi(argv[++run_args_start]);
      exit(TP_solver_run_group(self, type, (void *) &init, nb_steps, (void *) &inc)); 
    } else  {
      int j = 0;
      while( usageStrign[j] !=  NULL)
	printf("%s\n", usageStrign[j++]);
      TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "for the group run, unrecognized type of argument is given" );
    }  
  }
}

void
update_exe_parms(TP_exe_parms parms, TP_parm_type type,
		 void *init_val, int step, void *val, void *inc)
{
  double *Dinit, *Dinc, *Dval;
  int *Iinit, *Iinc, *Ival; 

  switch (type) {
  case (TP_value_tol) :
    Dinit = (double *) init_val; Dinc = (double *) inc; Dval = (double *) val; 
    *Dval = *Dinit * pow(*Dinc, (double) step);
    parms->value_tol = *Dval;
    break;
  case (TP_marko_tol) :
    Dinit = (double *) init_val; Dinc = (double *) inc;  Dval = (double *) val; 
    *Dval = *Dinit + *Dinc * step;
    parms->marko_tol = *Dval;
    break;
  case (TP_schur_density) :
    Dinit = (double *) init_val; Dinc = (double *) inc;  Dval = (double *) val; 
    *Dval = *Dinit + *Dinc * step;
    parms->density_tolerance = *Dval;
    break;
  case (TP_nb_candidates) :
    Iinit = (int *) init_val, Iinc = (int *) inc;  Ival = (int *) val; 
    *Ival = *Iinit + *Iinc * step;
    parms->nb_candidates_per_block = *Ival;
    break;
  case (TP_min_pivots) :
    Iinit = (int *) init_val, Iinc = (int *) inc; Ival = (int *) val;
    *Ival = *Iinit + *Iinc * step;
    parms->min_pivot_per_steps =  *Ival;
    break;
  case (TP_nb_threads) :
    Iinit = (int *) init_val, Iinc = (int *) inc; Ival = (int *) val;
    *Ival = *Iinit + *Iinc * step;
    parms->nb_threads =  *Ival;
    break;
  default :
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "unrecognized type of exe_parms");
  }
}

char *
get_outfile_prefix(TP_exe_parms exe_parms, TP_parm_type type)
{
  char *self = calloc(PATH_LENGTH, sizeof(*self));
  size_t length = 0;
  *self = '\0';
  
  if (exe_parms->matrix_file)
    snprintf(self, 2058, "%s", basename(exe_parms->matrix_file));

  length = strnlen(self, PATH_LENGTH - length);
  if (type == TP_value_tol)
    snprintf(self + length, PATH_LENGTH - length,"_MULTIValTol");
  else
    snprintf(self + length, PATH_LENGTH - length,"_%fValTol", exe_parms->value_tol);

  length = strnlen(self, PATH_LENGTH - length);
  if (type == TP_marko_tol)
    snprintf(self + length, PATH_LENGTH - length,"_MULTIMarkoTol");
  else
    snprintf(self + length, PATH_LENGTH - length,"_%fMarkoTol", exe_parms->marko_tol);
  
  length = strnlen(self, PATH_LENGTH - length);
  if (type == TP_nb_threads)
    snprintf(self + length, PATH_LENGTH - length,"_MULTIthreads");
  else
    snprintf(self + length, PATH_LENGTH - length,"_%dthreads", exe_parms->nb_threads);
  
  length = strnlen(self, PATH_LENGTH - length);
  if (type == TP_nb_candidates)
    snprintf(self + length, PATH_LENGTH - length,"_MULTIcandidates");
  else
    snprintf(self + length, PATH_LENGTH - length,"_%dcandidates", exe_parms->nb_candidates_per_block);
  
  length = strnlen(self, PATH_LENGTH - length);
  if (type == TP_schur_density)
    snprintf(self + length, PATH_LENGTH - length,"_MULTIdensityTol");
  else
    snprintf(self + length, PATH_LENGTH - length,"_%fdensityTol", exe_parms->density_tolerance);
  
  length = strnlen(self, PATH_LENGTH - length);
  if (type == TP_min_pivots)
    snprintf(self + length, PATH_LENGTH - length,"_MULTIminPivots");
  else
    snprintf(self + length, PATH_LENGTH - length,"_%dminPivots", exe_parms->min_pivot_per_steps);
  length = strnlen(self, PATH_LENGTH - length);
  self[length] = '\0';
  
  return self;
}

int 
TP_solver_run_group(TP_solver solver, TP_parm_type type, 
		    void *init_val, int nb_steps, void *inc)
{
  FILE *file;
  int i;
  TP_matrix A = TP_matrix_create();
  TP_exe_parms exe_parms = solver->exe_parms;
  char *file_ext = strrchr(exe_parms->matrix_file, '.');
  TP_vector X, rhs, sol;
  char *output_runs_file = get_outfile_prefix(exe_parms, type);
  char filename[PATH_LENGTH];
  double current_val;

  TP_verbose_create_dirs(solver->verbose->parms->output_dir);
  snprintf(filename, PATH_LENGTH, "%s/data/MULTI_%s_raw.dat", solver->verbose->parms->output_dir, output_runs_file);
 
  file = fopen(filename, "w+");
       
  if (!strcmp(file_ext, ".rb"))
    TP_read_rutherford_boeing(A, exe_parms->matrix_file);
  else  if (!strcmp(file_ext, ".mtl"))
    TP_read_mtl_file(A, exe_parms->matrix_file);
  else 
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"unrecognized file type format");
  
  solver->A = A;
  X   = TP_vector_create(A->n);
  rhs = TP_vector_create(A->n);
  sol = TP_vector_create(A->n);

  if (solver->exe_parms->RHS_file)
    TP_vector_read_file(X, solver->exe_parms->RHS_file);
  else 
    for(i=0; i<X->n; i++)
      X->vect[i] = (1+i) * 200;
  TP_vector_memset(rhs, 0.0);
  TP_matrix_SpMV(A, X, rhs);
  
  TP_verbose_print_parms_raw(exe_parms, type, file);
  for( i = 0; i < nb_steps; i++)
    {
      TP_solver run = TP_solver_create();
      TP_matrix matrix = TP_matrix_create();
      TP_exe_parms run_exe_parms = malloc(sizeof(*run_exe_parms));

      *run_exe_parms = *exe_parms;
      TP_matrix_copy(A, matrix);
      
      free(run->exe_parms);
      run->A = matrix;
      run->exe_parms = run->verbose->exe_parms = run_exe_parms;
      free(run->verbose->parms->output_dir);
      run->verbose->parms->output_dir = solver->verbose->parms->output_dir;
      run->verbose->parms->user_out_dir = 1;
      update_exe_parms(run->exe_parms, type, init_val, i, (void *) &current_val, inc);
      
      TP_solver_init(run);
      TP_vector_copy(rhs, sol);
      TP_solver_factorize(run);
      TP_solver_solve(run, sol);
      
      TP_solver_copmpute_norms(run, X, sol, rhs);
      
      TP_solver_finalize(run);
      TP_verbose_print_group_run(run->verbose, type, (void *) &current_val, i, file);
      TP_solver_destroy(run);
    }
  fclose(file);
  
  TP_vector_destroy(X);
  TP_vector_destroy(rhs);
  TP_vector_destroy(sol);

  TP_matrix_destroy(A);
  free(output_runs_file);
  TP_solver_dealloc(solver);
  
  return 0;
}

/* TODO: FOR NOW WE ASUME THAT N == M */
void
TP_solver_alloc_counters(TP_solver solver, int **col_count, int **row_count)
{
  int i, total;
  TP_counters counter = NULL;
  pthread_mutex_lock(&solver->counters_lock);

  for( i = 0; i < solver->nb_counters; i++)
    if ( solver->counters[i]->nb_used_counters  < solver->size_counters) {
      counter = solver->counters[i];
      break;
    }

  if (!counter) {
    solver->nb_counters++;
    int new = solver->nb_counters - 1;
    
    solver->counters  = realloc(solver->counters, solver->nb_counters * sizeof(*solver->counters));
    solver->counters[new] = calloc(1, sizeof(**solver->counters));
    solver->counters[new]->array = calloc((size_t) solver->size_counters * 2 * solver->A->n,  sizeof(*solver->counters[new]->array));
    solver->counters[new]->used_counters = calloc((size_t) solver->size_counters, sizeof(*solver->counters[new]->used_counters));
    counter = solver->counters[new];
  }
  total = solver->size_counters;

  for(i = 0; i < total; i++)
    if(!counter->used_counters[i])  {
      *col_count = &counter->array[i * solver->A->n * 2 ];
      *row_count = &counter->array[i * solver->A->n * 2 + solver->A->n];
      counter->nb_used_counters++;
      counter->used_counters[i] = 1;
      pthread_mutex_unlock(&solver->counters_lock);
      return;
    }

  TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"this should not happened!");
}

void
TP_solver_dealloc_counters(TP_solver solver, int *col_count, int *row_count)
{
  pthread_mutex_lock(&solver->counters_lock);
  int i;
  TP_counters counter = NULL;
  long diff;
  for( i = 0; i < solver->nb_counters; i++) {
    diff = col_count - solver->counters[i]->array;
    if ( diff >= 0  && diff < solver->A->n * 2 * solver->size_counters){
      counter = solver->counters[i];
      break;
    }
  }
  if (!counter) 
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"puffff");

  counter->used_counters[ diff / ( solver->A->n * 2) ] = 0;

  counter->nb_used_counters--;
  pthread_mutex_unlock(&solver->counters_lock);
}

void
TP_solver_alloc_internal(TP_solver self) 
{
  double total_extra_space = 1 + self->exe_parms->extra_space_inbetween + self->exe_parms->extra_space;
  int needed_pivots = self->A->n < self->A->m ? self->A->n : self->A->m;
  int i;
 
  self->S    = TP_schur_matrix_create();
  self->D    = TP_matrix_create();

  if (self->debug & (TP_DEBUG_VERBOSE_EACH_STEP | TP_DEBUG_GOSSIP_GIRL))
      TP_matrix_print(self->A, "A on input");
  TP_schur_matrix_allocate(self->S, self->A->n, self->A->m, self->A->nnz, self->debug, self->verbose,
			   self->exe_parms->nb_threads, self->exe_parms->extra_space, self->exe_parms->extra_space_inbetween);
  TP_schur_matrix_copy(self->A, self->S, self->exe_parms->value_tol);

  if (self->debug & (TP_DEBUG_VERBOSE_EACH_STEP | TP_DEBUG_GOSSIP_GIRL))
      TP_schur_matrix_print(self->S, "S on input");
  TP_matrix_allocate(self->D, self->A->n, 0, 0, 1.0, TP_Diag_matrix);
  self->L = TP_L_matrix_create(self->A->n, ((self->A->nnz - self->A->n) / 2 ) * total_extra_space);
  self->U = TP_U_matrix_create(self->A, total_extra_space);
  
  self->row_perm      = malloc((size_t) needed_pivots * sizeof(*self->invr_row_perm));
  self->col_perm      = malloc((size_t) needed_pivots * sizeof(*self->invr_col_perm));
  self->invr_row_perm = malloc((size_t) needed_pivots * sizeof(*self->invr_row_perm));
  self->invr_col_perm = malloc((size_t) needed_pivots * sizeof(*self->invr_col_perm));

  int_array_memset(self->col_perm,      TP_UNUSED_PIVOT, needed_pivots);
  int_array_memset(self->row_perm,      TP_UNUSED_PIVOT, needed_pivots);
  int_array_memset(self->invr_col_perm, TP_UNUSED_PIVOT, needed_pivots);
  int_array_memset(self->invr_row_perm, TP_UNUSED_PIVOT, needed_pivots);

  self->previous_pivots = malloc((size_t) self->exe_parms->nb_previous_pivots *  sizeof(*self->previous_pivots));
  int_array_memset(self->previous_pivots, INT_MAX / self->exe_parms->nb_previous_pivots, self->exe_parms->nb_previous_pivots);

  self->candidates = calloc(1, sizeof(*self->candidates));
  self->candidates->row        = malloc((size_t) self->A->n * sizeof(*self->candidates->row));
  self->candidates->marko      = malloc((size_t) self->A->n * sizeof(*self->candidates->marko));
  self->candidates->best_marko = malloc((size_t) self->exe_parms->nb_threads * sizeof(*self->candidates->best_marko));
 
  self->counters              = calloc( 1, sizeof(*self->counters));
  *self->counters             = calloc( 1, sizeof(**self->counters));

  self->random_col = create_randomize(self->A->n);

  (*self->counters)->array    = calloc((size_t) self->size_counters * 2 * self->A->n,  sizeof(*self->counters[0]->array));
  (*self->counters)->used_counters = calloc((size_t) self->size_counters, sizeof(*self->counters[0]->used_counters));
  self->nb_counters = 1;

  self->allocated_U_struct = self->A->n / 100;  
  self->allocated_U_struct = self->allocated_U_struct ? self->allocated_U_struct : 1;
  self->U_struct = calloc((size_t) self->allocated_U_struct, sizeof(*self->U_struct));
  self->allocated_L_struct = self->A->m / 100;  
  self->allocated_L_struct = self->allocated_L_struct ? self->allocated_L_struct : 1;
  self->L_struct = calloc((size_t) self->allocated_L_struct, sizeof(*self->L_struct));
  
  self->Luby = TP_Luby_create(self->S);

  self->cols = malloc((size_t) self->A->n * sizeof(*self->cols));
  for( i = 0; i < self->A->n; i++)
    self->cols[i] = i;
  self->rows = malloc((size_t) self->A->m * sizeof(*self->rows));
  for( i = 0; i < self->A->m; i++)
    self->rows[i] = i;
  self->distributions = malloc((size_t) (self->exe_parms->nb_threads + 1)  * sizeof(*self->distributions));
  
  /* srand(time(NULL)); */
  self->seeds = malloc((size_t) self->exe_parms->nb_threads * sizeof(*self->seeds));
  for(i = 0; i < self->exe_parms->nb_threads; i++)
    /* TODO: Put an argument "determenistic" and call srand before if activated */
    self->seeds[i] = rand();

  self->workspace = malloc((size_t) self->exe_parms->nb_threads * sizeof(*self->workspace));
  for(i = 0; i < self->exe_parms->nb_threads; i++)
    self->workspace[i] = malloc((size_t) self->A->n * (sizeof(int) + sizeof(double)));
  
  self->logical_cols = calloc((size_t) self->A->n, sizeof(*self->logical_cols));
  self->logical_rows = calloc((size_t) self->A->n, sizeof(*self->logical_rows));
  
  pthread_mutex_init(&self->counters_lock, NULL); 
}

void
TP_solver_init(TP_solver self)
{
  self->verbose->n = self->A->n;
  self->verbose->m = self->A->m;
  self->verbose->nnz_input = self->verbose->nnz_final = self->A->nnz;
  self->verbose->Luby = self->exe_parms->luby_algo;
  
  self->verbose->parms->outfiles_prefix = get_outfile_prefix(self->exe_parms, TP_parm_none);

  if (self->exe_parms->trace) 
    self->verbose->paje = TP_paje_create(self->exe_parms->nb_threads);

  if (!is_plasma_init) {
    plasma_init(1);
    is_plasma_init = 1;
  }
  
  if ( !(self->debug & TP_CHECK_DENSE_W_TP_PERM) )
    TP_solver_alloc_internal(self);

  TP_verbose_create_dirs(self->verbose->parms->output_dir);
}

void
TP_solver_read_matrix(TP_solver self)
{
  char *file_ext = strrchr(self->exe_parms->matrix_file, '.');
  self->A        = TP_matrix_create();
  
  if (!strcmp(file_ext, ".rb"))
    TP_read_rutherford_boeing(self->A, self->exe_parms->matrix_file);
  else if (!strcmp(file_ext, ".mtl")) 
    TP_read_mtl_file(self->A, self->exe_parms->matrix_file);
  else 
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"unsupported matrix file");
}

void 
TP_solver_get_pivots(TP_solver self, TP_pivot_set set)
{
  TP_pivot_cell cells = set->cells;
  int n = self->A->n, i, k, nnz ;

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
                            
  for ( i = 0, k=0, nnz = 0; i < n; i++) {
    if (set->cols_count[i] < set->base)
      continue;
    if (self->invr_col_perm[i] != TP_UNUSED_PIVOT)  {
      if (set->cols_count[i] > set->base) {
	TP_warning(__FUNCTION__, __FILE__, __LINE__, "col is supposed to be a pivot bit is larger then base ");
      } else {
	continue;
      }
    }
    self->U_struct[k].col = i;
    self->U_struct[k].nb_elem = set->cols_count[i] - set->base + 1;
    nnz += self->U_struct[k].nb_elem;
    if( ++k >= self->allocated_U_struct)
      {
	self->allocated_U_struct *= 2;
	self->U_struct = realloc(self->U_struct, (size_t) self->allocated_U_struct * sizeof(*self->U_struct));
      }
  }
  self->n_U_structs = k;
  self->nnz_U_structs =  nnz;
  
  for ( i = 0, k=0, nnz = 0; i < n; i++) {
    if (set->rows_count[i] < set->base)
      continue;
    if (self->invr_row_perm[i] != TP_UNUSED_PIVOT)  {
      if (set->rows_count[i] > set->base) {
	TP_warning(__FUNCTION__, __FILE__, __LINE__, "row is supposed to be a pivot bit is larger then base ");
      } else {
	continue;
      }
    }
    self->L_struct[k].col = i;
    self->L_struct[k].nb_elem = set->rows_count[i] - set->base + 1;
    nnz += self->L_struct[k].nb_elem;
    if( ++k >= self->allocated_L_struct)
      {
	self->allocated_L_struct *= 2;
	self->L_struct = realloc(self->L_struct, (size_t) self->allocated_L_struct * sizeof(*self->L_struct));
      }
  }
  self->n_L_structs = k;
  self->nnz_L_structs = nnz;

  self->rows_count = set->rows_count;
  self->cols_count = set->cols_count;

  self->found_pivots += new_pivots;
  TP_verbose_update_pivots(self->verbose, new_pivots);
}



void 
TP_solver_get_Luby_pivots(TP_solver self, TP_Luby Luby, int new_Luby_pivots)
{
  TP_schur_matrix S = self->S;
  TP_verbose verbose = self->verbose;

  int n = S->n, nb_cols = S->n - self->done_pivots;
  int done_pivots = self->done_pivots;
  int new_pivots =  new_Luby_pivots + self->nb_col_singletons + self->nb_row_singletons;
  int all_pivots = new_pivots + done_pivots;
  int i, base = self->step + 1;
  int nb_threads = self->exe_parms->nb_threads;
  int distribution_perms[nb_threads+1], distribution_n[nb_threads+1];
  int row_sizes[nb_threads], col_sizes[nb_threads];

  int *logical_cols = self->logical_cols;
  int *logical_rows = self->logical_rows;

  int *row_perms = self->row_perm;
  int *col_perms = self->col_perm;
  int *invr_row_perms = self->invr_row_perm;
  int *invr_col_perms = self->invr_col_perm;

  int *cols = self->cols;
  int *rows = self->rows;
  
  self->previous_step_pivots = new_pivots;

  self->n_U_structs = 0;
  self->nnz_U_structs = 0;
  self->n_L_structs = 0;
  self->nnz_L_structs = 0;

  *distribution_perms = done_pivots + self->nb_col_singletons + self->nb_row_singletons;
  for(i = 1; i < nb_threads; i++)
    distribution_perms[i] = *distribution_perms + (new_Luby_pivots / nb_threads) * i;
  distribution_perms[nb_threads] = all_pivots;

  *distribution_n = 0;
  for(i = 1; i < nb_threads; i++)
    distribution_n[i] = (nb_cols / nb_threads) * i;
  distribution_n[nb_threads] = nb_cols;

#pragma omp parallel num_threads(nb_threads) shared(distribution_perms, distribution_n, base, nb_cols, col_sizes, row_sizes, n, verbose)
  {
  TP_verbose_trace_start_event(verbose, TP_GETTING_PIVOTS);
  int i, j;
  int me = omp_get_thread_num();
  int start = distribution_perms[me] ;
  int end   = distribution_perms[me+1];
  int *my_col_perms = (int *) self->workspace[me];
  int *my_row_perms = &my_col_perms[nb_cols];
  int my_nb_cols = 0, my_nb_rows = 0;
  
  for( i = start; i < end; i++)
    {
      invr_col_perms[col_perms[i]] = i;
      invr_row_perms[row_perms[i]] = i;
    }

#pragma omp barrier

  for ( i = start; i < end; i++)
    {
      CSC_struct *CSC = &S->CSC[col_perms[i]];
      int *rows = CSC->row; 
      int nb_elem = CSC->nb_elem;
      for ( j = 0; j < nb_elem; j++) {
	logical_rows[rows[j]] = base;
      }
      CSR_struct *CSR = &S->CSR[row_perms[i]];
      int *cols = CSR->col;
      nb_elem = CSR->nb_elem;
      for ( j = 0; j < nb_elem; j++) {
	logical_cols[cols[j]] = base;
      }
    }
#pragma omp barrier
#pragma omp single
  {
  start = self->done_pivots + self->nb_row_singletons;
  end   = self->done_pivots + self->nb_row_singletons + self->nb_col_singletons;
  for ( i = start; i < end; i++)
    {
      CSR_struct *CSR = &S->CSR[row_perms[i]];
      int *cols = CSR->col;
      int nb_elem = CSR->nb_elem;
      for ( j = 0; j < nb_elem; j++) {
	int col = cols[j];
	int perm = invr_col_perms[col];

	if (perm != TP_UNUSED_PIVOT && perm > end) 
	  if (col_perms[perm] < n)
	    col_perms[perm] += n;
	
	logical_cols[col] = base;
      }
    }      
  
  *distribution_perms = self->done_pivots;
  for(i = 1; i < nb_threads; i++)
    distribution_perms[i] = *distribution_perms + (self->nb_row_singletons / nb_threads) * i;
  distribution_perms[nb_threads] = *distribution_perms + self->nb_row_singletons;
  }
#pragma omp barrier

  start = distribution_perms[me] ;
  end   = distribution_perms[me+1];

  for ( i = start; i < end; i++)
    {
      CSC_struct *CSC = &S->CSC[col_perms[i]];
      int *rows = CSC->row; 
      int nb_elem = CSC->nb_elem;
      for ( j = 0; j < nb_elem; j++) {
	logical_rows[rows[j]] = base;
      }
    }
#pragma omp barrier

  start = distribution_n[me];
  end = distribution_n[me+1];

  for( i = start; i < end; i++)  {
    int col = cols[i];
    if(logical_cols[col] == base && invr_col_perms[col] == TP_UNUSED_PIVOT) 
      my_col_perms[my_nb_cols++] = col;
  }

  for( i = start; i < end; i++)  {
    int row = rows[i];
    if(logical_rows[row] == base && invr_row_perms[row] == TP_UNUSED_PIVOT) 
      my_row_perms[my_nb_rows++] = row;
  }

  col_sizes[me] = my_nb_cols;
  row_sizes[me] = my_nb_rows;

#pragma omp barrier
#pragma omp single 
  {
    for (i = 1; i < nb_threads; i++) {
      col_sizes[i] += col_sizes[i-1];
      row_sizes[i] += row_sizes[i-1];
    }
    self->n_U_structs = col_sizes[nb_threads-1];
    self->n_L_structs = row_sizes[nb_threads-1];
  }
#pragma omp barrier
  if (me) {
    memcpy(&col_perms[all_pivots + col_sizes[me-1]], my_col_perms,
	   (size_t) (col_sizes[me] - col_sizes[me-1]) * sizeof(*col_perms));
    
    memcpy(&row_perms[all_pivots + row_sizes[me-1]], my_row_perms,
	   (size_t) (row_sizes[me] - row_sizes[me-1]) * sizeof(*row_perms));
  } else { 
    memcpy(&col_perms[all_pivots], my_col_perms, (size_t) col_sizes[me] * sizeof(*col_perms));
    memcpy(&row_perms[all_pivots], my_row_perms, (size_t) row_sizes[me] * sizeof(*row_perms));
  }
  TP_verbose_trace_stop_event(verbose);
  }

  self->found_pivots = all_pivots;
  TP_verbose_update_pivots(self->verbose, new_pivots);
}

void
TP_solver_check_counters(TP_solver self)
{
  int i; 
  
  for (i = 0; i < self->nb_counters; i++) 
    TP_check_counters(i, self->counters[i]->array, self->counters[i]->used_counters,
		      self->done_pivots + 1, self->size_counters, self->A->n);
}

void
TP_solver_find_pivot_set(TP_solver self)
{
  TP_pivot_list list;
  TP_exe_parms exe_parms = self->exe_parms;
  TP_verbose_per_step step = TP_verbose_get_step(self->verbose);
  TP_verbose verbose = self->verbose;
  int nb_threads = self->exe_parms->nb_threads;
    
  TP_verbose_trace_start_event(verbose, TP_GET_SINGELTONS);
  TP_schur_get_singletons(self->S, self->done_pivots, self->previous_step_pivots,
  			  &self->nb_col_singletons, &self->nb_row_singletons,
			  self->cols, self->rows, self->distributions,
  			  self->done_pivots, self->col_perm, self->row_perm,
  			  self->invr_col_perm, self->invr_row_perm);
  TP_verbose_trace_stop_event(verbose);

  if (exe_parms->luby_algo) {
    int new_pivots = 0, best_marko, nb_cols = self->S->n - self->done_pivots;
    int *distributions = self->distributions, best_markos[nb_threads];
    int candidates[nb_threads];
    int i;

    *distributions = 0;
    for (i = 1; i < nb_threads; i++) 
      distributions[i] = (nb_cols / nb_threads) * i;
    distributions[nb_threads] = nb_cols;

#pragma omp parallel num_threads(nb_threads) shared(nb_cols, best_marko, distributions)
    {
    TP_verbose_trace_start_event(verbose, TP_GET_ELIGEBLE);
    int me = omp_get_thread_num();
    int *my_col_perms = (int *) self->workspace[me];
    int *my_row_perms = &my_col_perms[nb_cols];
    int *global_col_perms = (int *) self->workspace[0];
    int *global_row_perms = &global_col_perms[nb_cols];

    TP_verbose_start_timing(&step->timing_extracting_candidates);
    
    best_markos[me] = TP_Luby_get_eligible(self->S, self->Luby, exe_parms->value_tol, self->invr_col_perm, self->invr_row_perm, self->cols, distributions[me], distributions[me + 1], 0);
    TP_verbose_trace_stop_event(verbose);
#pragma omp barrier    
#pragma omp single 
    {
    TP_verbose_trace_start_event(verbose, TP_AUXILIARY);
    best_marko = *best_markos;
    for(i = 1; i < nb_threads; i++)
      best_marko = best_marko > best_markos[i] ? best_markos[i] : best_marko;
    best_marko = !best_marko ? 1 : best_marko;
    best_marko *= exe_parms->marko_tol;
    TP_verbose_trace_stop_event(verbose);
    }
#pragma omp barrier    
    TP_verbose_trace_start_event(verbose, TP_ASSIGN_SCORE);
    candidates[me] = TP_Luby_assign_score(self->Luby, self->S, best_marko, &self->seeds[me], my_col_perms, my_row_perms, self->cols, distributions[me], distributions[me + 1]);
    TP_verbose_trace_stop_event(verbose);
#pragma omp barrier
    TP_verbose_trace_start_event(verbose, TP_AUXILIARY);
#pragma omp single
    {
      for( i = 1; i < nb_threads; i++)
	candidates[i] += candidates[i-1];
      step->nb_candidates = candidates[nb_threads-1];
    } 
#pragma omp barrier    

    if (me) {
      memcpy(&global_col_perms[candidates[me - 1]],  my_col_perms, (size_t) (candidates[me] - candidates[me-1]) * sizeof(*self->col_perm));
      memcpy(&global_row_perms[candidates[me - 1]], my_row_perms, (size_t) (candidates[me] - candidates[me-1]) * sizeof(*self->row_perm));
    }
    TP_verbose_trace_stop_event(verbose);
    }
    
    TP_verbose_stop_timing(&step->timing_extracting_candidates);

    TP_verbose_start_timing(&step->timing_merging_pivots);
    *distributions = 0;
    for (i = 1; i < nb_threads; i++) 
      distributions[i] = (step->nb_candidates / nb_threads) * i;
    distributions[nb_threads] = step->nb_candidates;
    
#pragma omp parallel num_threads(nb_threads) shared(nb_cols)
    {
    TP_verbose_trace_start_event(verbose, TP_FIRST_PASS);
    int me = omp_get_thread_num();
    int my_size = distributions[me + 1] - distributions[me];
    int *tmp = (int *) self->workspace[0];
    int *my_col_perms = &tmp[distributions[me]];
    int *my_row_perms = &tmp[nb_cols + distributions[me]];
    int *global_col_perms = &self->col_perm[self->done_pivots + self->nb_col_singletons + self->nb_row_singletons];
    int *global_row_perms = &self->row_perm[self->done_pivots + self->nb_col_singletons + self->nb_row_singletons];

    TP_Luby_first_pass(self->Luby, self->S, my_col_perms, my_row_perms, my_size);
    TP_verbose_trace_stop_event(verbose);

#pragma omp barrier
    TP_verbose_trace_start_event(verbose, TP_DISCARDING_PIVOTS);
    candidates[me] = TP_Luby_second_pass(self->S, self->Luby, my_col_perms, my_row_perms, my_size);
    TP_verbose_trace_stop_event(verbose);
#pragma omp barrier
    TP_verbose_trace_start_event(verbose, TP_AUXILIARY);
#pragma omp single
    {
      self->Luby->chosen_base += 3;
      for( i = 1; i < nb_threads; i++)
	candidates[i] += candidates[i-1];
    }
#pragma omp barrier
    if (me) {
      memcpy(&global_col_perms[candidates[me-1]], my_col_perms, (size_t) (candidates[me] - candidates[me-1]) * sizeof(*my_col_perms));
      
      memcpy(&global_row_perms[candidates[me-1]], my_row_perms,(size_t) (candidates[me] - candidates[me-1]) * sizeof(*self->row_perm));
    } else {
      memcpy(global_col_perms, my_col_perms, (size_t) candidates[me] * sizeof(*my_col_perms));
      memcpy(global_row_perms, my_row_perms, (size_t) candidates[me] * sizeof(*my_row_perms));
    }
    TP_verbose_trace_stop_event(verbose);
    }
    new_pivots = candidates[nb_threads-1];
    
    if (self->debug & (TP_DEBUG_VERBOSE_EACH_STEP | TP_DEBUG_GOSSIP_GIRL)) {
      print_int_array(self->col_perm, self->A->n, "col_perms");
      print_int_array(self->row_perm, self->A->n, "row_perms");
    }
    
    TP_solver_get_Luby_pivots(self, self->Luby, new_pivots);
    TP_verbose_stop_timing(&step->timing_merging_pivots);

  } else {
    
    if (self->debug & TP_CHECK_COUNTERS )
      TP_solver_check_counters(self);
    
    TP_verbose_start_timing(&step->timing_extracting_candidates);
    list = get_possible_pivots(self, self->S, self->random_col, self->candidates, nb_threads,
			       exe_parms->value_tol, exe_parms->marko_tol,
			       exe_parms->nb_candidates_per_block);
    TP_verbose_stop_timing(&step->timing_extracting_candidates);
    
    if( !list->nb_elem )
      TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "no possible pivot were found");
    
    TP_verbose_start_timing(&step->timing_merging_pivots);
    list = merge_pivot_sets(list, self->S, self);
    
    if (list->nb_elem > 1 )
      TP_warning(__FUNCTION__, __FILE__, __LINE__, "not unique set");
    
    TP_solver_get_pivots(self, list->sets);
    
    if (self->debug & TP_CHECK_COUNTERS )
      TP_check_current_counters(self->S, &self->col_perm[self->done_pivots], &self->row_perm[self->done_pivots],
				self->found_pivots - self->done_pivots, 
				self->cols_count, self->rows_count, self->done_pivots + 1);
    
    if (self->debug & (TP_DEBUG_VERBOSE_EACH_STEP | TP_DEBUG_GOSSIP_GIRL))
      print_pivot_list(list, "found pivots");
    TP_pivot_list_destroy(list, self);
    
    TP_verbose_stop_timing(&step->timing_merging_pivots);
  }
}

void
TP_solver_update_matrix(TP_solver self)
{
  TP_L_matrix L = self->L;
  TP_matrix D = self->D;
  TP_U_matrix U = self->U;
  TP_schur_matrix S = self->S;
  int nb_pivots = self->found_pivots - self->done_pivots;
  TP_verbose_per_step step = TP_verbose_get_step(self->verbose);

  TP_verbose_start_timing(&step->timing_update_LD);
  TP_schur_matrix_update_LD(S, L, U, D, &self->row_perm[self->done_pivots],
			    &self->col_perm[self->done_pivots], nb_pivots,
			    self->invr_row_perm, self->nb_row_singletons,
			    self->nb_col_singletons, self->workspace);
  TP_verbose_stop_timing(&step->timing_update_LD);
  
  TP_verbose_start_timing(&step->timing_update_S);
  TP_schur_matrix_update_S(S, L, U, &self->col_perm[self->found_pivots],
			   self->n_U_structs, self->nnz_U_structs, self->invr_row_perm,
			   nb_pivots, &self->row_perm[self->done_pivots], self->workspace,
			   self->exe_parms->value_tol);
  TP_verbose_stop_timing(&step->timing_update_S);
  
  TP_verbose_start_timing(&step->timing_update_U);
  TP_schur_matrix_update_S_rows(S, &self->row_perm[self->found_pivots],
				self->n_L_structs, self->nnz_L_structs,
				self->invr_col_perm, nb_pivots, self->row_perm, self->done_pivots);
  TP_verbose_stop_timing(&step->timing_update_U);

  self->done_pivots = self->found_pivots;
  self->step++;
    
  if (self->debug & TP_CHECK_SCHUR_SYMETRY )
    TP_schur_matrix_check_symetry(self->S);

  if (self->debug &  TP_CHECK_PIVOTS)
    TP_schur_matrix_check_pivots(self->S, self->row_perm, self->col_perm, 
				 self->invr_row_perm, self->invr_col_perm,
				 self->done_pivots);
  
  if (self->debug & (TP_DEBUG_VERBOSE_EACH_STEP | TP_DEBUG_GOSSIP_GIRL)) {
    TP_schur_matrix_print(S, "S after update");
    TP_L_matrix_print(L, "L after update");
    TP_matrix_print(D, "D after update");
    TP_U_matrix_print(U, "U after update");
  }
  
  if (self->debug & TP_CHECK_SCHUR_DOUBLES) { 
    TP_schur_check_doubles(S);
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
			 int max_dense_schur,
			 int debug,
			 TP_verbose verbose)
{
  long n_schur, m_schur, i, sum_pivots;
  int retval = 1;

  if ( debug & TP_CHECK_TP_W_PLASMA_PERM ||
       debug & TP_CHECK_DENSE_W_TP_PERM ) {
    verbose->reason = TP_reason_because;
    return 0;
  }

  n_schur = (long) S->n - nb_done_pivots;
  m_schur = (long) S->m - nb_done_pivots;
  if ( S->nnz > (long) n_schur * m_schur * density_tolerance) {
    verbose->reason = TP_reason_density;
    retval = 0;
  }

  if (retval)
    {
      for( i = 0, sum_pivots = 0; i < nb_previous_pivots; i++)
	sum_pivots += previous_pivots[i];
      if (sum_pivots < min_pivot_per_steps) {
	verbose->reason = TP_reason_no_pivots;
	retval = 0;
      }
    }

  if (nb_done_pivots >= nb_needed_pivots && retval) {
    verbose->reason = TP_reason_no_pivots;
    retval = 0;
  }
  
  if ( !retval && max_dense_schur < (nb_needed_pivots - nb_done_pivots) ) {
    verbose->reason |= TP_reason_dense_too_large;
  }
  return retval;
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

  if (self->debug &  TP_CHECK_SCHUR_DOUBLES)
    TP_schur_check_doubles(self->S);  

  TP_verbose_start_timing(&verbose->timing_facto);
  TP_verbose_start_timing(&verbose->timing_facto_sparse);
  while ( TP_continue_pivot_search(self->S, self->done_pivots, needed_pivots, 
				   previous_pivots, nb_previous_pivots,
				   exe_parms->density_tolerance,
				   exe_parms->min_pivot_per_steps,
				   exe_parms->max_dense_schur,
				   self->debug, verbose) )
    { 
       TP_verbose_per_step  step = TP_verbose_step_start(verbose);
       TP_verbose_start_timing(&step->timing_step);
       TP_verbose_start_timing(&step->timing_pivot_search);
       TP_solver_find_pivot_set(self);
       /* if ( !self->nb_row_singletons && !self->nb_col_singletons) */
       previous_pivots[nb_pivot_blocks++ % nb_previous_pivots] = self->found_pivots - self->done_pivots;
       TP_verbose_stop_timing(&step->timing_pivot_search);
       
       TP_verbose_start_timing(&step->timing_apply_perms);
       TP_solver_update_matrix(self);
       TP_verbose_stop_timing(&step->timing_apply_perms);
       TP_verbose_stop_timing(&step->timing_step);
    }
  TP_verbose_stop_timing(&verbose->timing_facto_sparse);

  if (verbose->reason & TP_reason_dense_too_large)
    return;
  
  TP_verbose_start_timing(&verbose->timing_convert_schur);
  TP_verbose_trace_start_event(verbose, TP_CONVERT_MATRIX);
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
  TP_verbose_trace_stop_event(verbose);
  TP_verbose_stop_timing(&verbose->timing_convert_schur);

  TP_verbose_start_timing(&verbose->timing_facto_dense);
  TP_verbose_trace_start_event(verbose, TP_DENSE_FACTORIZATION);
  if (self->debug & TP_CHECK_DENSE_W_TP_PERM)  {
    TP_dense_2D_facto(self->A_debug);
  } else {
    TP_dense_matrix_factorize(self->S_dense, exe_parms->nb_threads);
    // Handeling the pivots
    memcpy(&self->row_perm[self->done_pivots], self->S_dense->original_rows,
	   (size_t) (self->A->m - self->done_pivots) * sizeof(*self->row_perm));
    memcpy(&self->col_perm[self->done_pivots], self->S_dense->original_cols,
	   (size_t) (self->A->n - self->done_pivots) * sizeof(*self->col_perm));
    int i;
    for(i = self->done_pivots; i < needed_pivots; i++) {
      self->invr_row_perm[self->row_perm[i]] = i;
      self->invr_col_perm[self->col_perm[i]] = i;
    }
    self->dense_pivots = needed_pivots - self->done_pivots;
    TP_verbose_update_dense_pivots(verbose, self->dense_pivots);
    verbose->nnz_final   = self->L->nnz + self->U->nnz + self->D->n + self->S_dense->n * self->S_dense->m;
    verbose->nnz_L       = self->L->nnz  + ( ( self->S_dense->n * self->S_dense->m - self->S_dense->m) / 2 )  ;
    verbose->nnz_U       = self->U->nnz + self->D->n + (( self->S_dense->n * self->S_dense->m - self->S_dense->m) / 2 + self->S_dense->m ) ;
    verbose->nnz_S_dense = self->S_dense->n * self->S_dense->m;
  }

  TP_verbose_trace_stop_event(verbose);
  TP_verbose_stop_timing(&verbose->timing_facto_dense);
  TP_verbose_stop_timing(&verbose->timing_facto);
}

void
TP_solver_solve(TP_solver self, TP_vector RHS)
{
  if (self->verbose->reason & TP_reason_dense_too_large)
    return;
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
    TP_L_matrix_solve(self->L, RHS, self->invr_row_perm);
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
    TP_U_matrix_solve(self->U, self->D, RHS, self->col_perm, self->invr_row_perm, self->dense_pivots);
    /* TP_matrix_solve_UD(self->U, self->D, RHS, self->invr_col_perm); */
    if (self->debug & TP_DEBUG_GOSSIP_GIRL) 
      TP_vector_print(RHS, "after backward solve");
    TP_vector_permute(RHS, self->invr_col_perm);
    if (self->debug & (TP_DEBUG_VERBOSE_EACH_STEP | TP_DEBUG_GOSSIP_GIRL)) 
      TP_vector_print(RHS, "after solve operation");
    TP_verbose_stop_timing(&verbose->timing_solve_U);
  }
  TP_verbose_computed_norms(verbose);
  TP_verbose_stop_timing(&verbose->timing_solve);
}

void
TP_solver_copmpute_norms(TP_solver self,       TP_vector X,
			 TP_vector X_computed, TP_vector rhs)
{
  if (self->verbose->reason & TP_reason_dense_too_large)
    return;
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
  if (self->verbose->reason & TP_reason_dense_too_large)
    return;
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
  TP_verbose_draw_graph(self->verbose);
}

void
TP_solver_destroy(TP_solver self)
{
  int i;

  if (self->debug & TP_CHECK_DENSE_W_TP_PERM) {
    TP_dense_2D_destroy(self->A_debug);
  } else {
    TP_matrix_destroy(self->A);
    TP_schur_matrix_destroy(self->S);
    TP_L_matrix_destroy(self->L);
    TP_matrix_destroy(self->D);
    TP_U_matrix_destroy(self->U);
    if(self->S_dense)
      TP_dense_matrix_destroy(self->S_dense);

    free(self->candidates->row);
    free(self->candidates->marko);
    free(self->candidates->best_marko);
    free(self->candidates);
    
    for(i = 0; i < self->nb_counters; i++) {
      free(self->counters[i]->array);
      free(self->counters[i]->used_counters);
      free(self->counters[i]);
    }

    free(self->counters);
    free(self->U_struct);
    free(self->L_struct);
    
    pthread_mutex_destroy(&self->counters_lock);
    
    free(self->row_perm);
    free(self->col_perm);
    free(self->invr_row_perm);
    free(self->invr_col_perm);
    free(self->random_col);
    free(self->previous_pivots);

    TP_Luby_destroy(self->Luby);
  }
  
  TP_solver_dealloc(self);
}
