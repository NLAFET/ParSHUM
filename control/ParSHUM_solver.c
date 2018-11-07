#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h>
#include <plasma.h>
#include <math.h>
#include <limits.h>

#include "ParSHUM_verbose.h"
#include "ParSHUM_enum.h" 
#include "ParSHUM_matrix.h" 
#include "ParSHUM_dense.h" 
#include "ParSHUM_schur_matrix.h" 
#include "ParSHUM_pivot_list.h" 
#include "ParSHUM_paje.h"
#include "ParSHUM_auxiliary.h"
#include "ParSHUM_solver.h" 

const char *usageStrign[] = {
  "usage test: [--help] [--matrix matrix] [--RHS_file file] [--debug_mode] [--verbosity level] [--marko_tol tol] [--value_tol tol]",
  "            [--extra_space factor] [--extra_space_inbetween factor] [--nb_threads #threads] [--nb_candidates_per_block #blocks] ", 
  "            [--output_dir dir] [--output_file file] [--nb_previous_pivots #pivtos] [--schur_density_tolerance tol]",
  "            [--min_pivot_per_steps #steps] [--prog_name name ] [--check_schur_symetry] [--check_schur_memory] [--check_pivots]",
  "            [--check_ParSHUM_with_plasma_perm] [--check_dense_with_ParSHUM_perm] [--print_each_step] [--check_GC]",
  "            [--group_run value_tol|marko_tol|schur_density|nb_candidates|min_pivots|nb_threads init inc nb_steps]",
  "            [--counters_size #double_counters] [--check_counters] [--check_schur_doubles] [--max_dense_schur size]",
  "            [--luby_algorithm] [--singeltons_relaxation tol] [--trace]", 
  NULL,
};

int is_plasma_init;

int  ParSHUM_solver_run_group(ParSHUM_solver solver, ParSHUM_parm_type type, 
			 void *init_val, int nb_steps, void *inc);

ParSHUM_solver 
ParSHUM_solver_create()
{
  ParSHUM_solver self = calloc(1, sizeof(*self));
  self->exe_parms = calloc(1, sizeof(*self->exe_parms));

  self->verbosity      = 1;
  self->size_counters = 100;

  self->exe_parms->nb_threads              = 1;
  self->exe_parms->value_tol               = 0.1;
  self->exe_parms->singeltons_relaxation   = 0.01;
  self->exe_parms->marko_tol               = 4;
  self->exe_parms->extra_space             = 1.0;
  self->exe_parms->extra_space_inbetween   = 1.0;
  self->exe_parms->nb_candidates_per_block = 10;
  self->exe_parms->nb_previous_pivots      = 5;
  self->exe_parms->min_pivot_per_steps     = 20;
  self->exe_parms->density_tolerance       = 0.2;
  self->exe_parms->max_dense_schur         = 10000;
  
  self->verbose = ParSHUM_verbose_create(self->exe_parms);
  return self;
}

void
ParSHUM_solver_dealloc(ParSHUM_solver self)
{
  ParSHUM_verbose_destroy(self->verbose);
  
  free(self->exe_parms);
  free(self);
}

int
check_ParSHUM_with_plasma_perm(int argc, char **argv)
{
  ParSHUM_solver plasma;
  ParSHUM_vector X, sol_plasma, sol_ParSHUM;
  
  plasma = ParSHUM_solver_create();
  ParSHUM_solver_parse_args(plasma, argc, argv, 1);
  plasma->exe_parms->density_tolerance   = 1.0;
  plasma->exe_parms->min_pivot_per_steps = 5;
  plasma->exe_parms->nb_previous_pivots  = 5;
  plasma->debug |= ParSHUM_CHECK_ParSHUM_W_PLASMA_PERM;
  ParSHUM_solver_read_matrix(plasma);
  ParSHUM_solver_init(plasma);
  
  X   = ParSHUM_vector_create(plasma->A->n);
  sol_plasma  = ParSHUM_vector_create(plasma->A->n);
  sol_ParSHUM = ParSHUM_vector_create(plasma->A->n);
  
  ParSHUM_vector_read_file(X, plasma->exe_parms->RHS_file);
  ParSHUM_vector_copy(X, sol_plasma);

  ParSHUM_solver_factorize(plasma);
  ParSHUM_solver_solve(plasma, sol_plasma);

  ParSHUM_solver_compute_norms(plasma, X, sol_plasma);

  ParSHUM_solver_finalize(plasma);

  /* apply plasma row permutation to A */
  int *plasma_perms = ParSHUM_dense_get_row_perms(plasma->S_dense, plasma->row_perm);
  ParSHUM_matrix debug_matrix = ParSHUM_matrix_permute(plasma->A, plasma->col_perm, plasma_perms);
  ParSHUM_solver debug_solver = ParSHUM_solver_create();

  debug_solver->A = debug_matrix;
  debug_solver->debug |= ParSHUM_CHECK_ParSHUM_W_PLASMA_PERM;
  debug_solver->exe_parms->density_tolerance = 1.0;
  debug_solver->exe_parms->min_pivot_per_steps = 5;
  debug_solver->exe_parms->nb_previous_pivots  = 5;
  ParSHUM_solver_init(debug_solver);

  ParSHUM_vector_permute(X, plasma_perms, X->n);
  ParSHUM_vector_copy(X, sol_ParSHUM);

  ParSHUM_solver_factorize(debug_solver);
  ParSHUM_solver_solve(debug_solver, sol_ParSHUM);

  ParSHUM_solver_compute_norms(debug_solver, X, sol_ParSHUM);

  ParSHUM_solver_finalize(debug_solver);

  free(plasma_perms);
  ParSHUM_vector_destroy(X);
  ParSHUM_vector_destroy(sol_plasma);
  ParSHUM_vector_destroy(sol_ParSHUM);
  ParSHUM_solver_destroy(debug_solver);
  ParSHUM_solver_destroy(plasma);

  return 0;
}

int
check_dense_with_ParSHUM_perm(int argc, char **argv)
{
  ParSHUM_solver self;
  ParSHUM_vector X, sol_ParSHUM, sol_dense;
  
  self = ParSHUM_solver_create();
  ParSHUM_solver_parse_args(self, argc, argv, 1);
  self->exe_parms->density_tolerance = 1.0;
  self->exe_parms->min_pivot_per_steps = 5;
  self->exe_parms->nb_previous_pivots  = 5;
  ParSHUM_solver_read_matrix(self);
  ParSHUM_solver_init(self);
  
  X   = ParSHUM_vector_create(self->A->n);
  sol_ParSHUM = ParSHUM_vector_create(self->A->n);
  sol_dense   = ParSHUM_vector_create(self->A->n);
  
  ParSHUM_vector_read_file(X, self->exe_parms->RHS_file);
  ParSHUM_vector_copy(X, sol_ParSHUM);
  ParSHUM_vector_copy(X, sol_dense);

  ParSHUM_solver_factorize(self);
  ParSHUM_solver_solve(self, sol_ParSHUM);
  ParSHUM_solver_compute_norms(self, X, sol_ParSHUM);

  ParSHUM_solver_finalize(self);
 
  ParSHUM_solver dense_solver = ParSHUM_solver_create();
  ParSHUM_solver_parse_args(dense_solver, argc, argv, 1);

  dense_solver->A = self->A;

  dense_solver->debug |= ParSHUM_CHECK_DENSE_W_ParSHUM_PERM;
  dense_solver->row_perm = self->row_perm;
  dense_solver->col_perm = self->col_perm;
  dense_solver->invr_row_perm = self->invr_row_perm;
  dense_solver->invr_col_perm = self->invr_col_perm;
  dense_solver->verbose->parms->prog_name = "Dense solver";

  ParSHUM_solver_init(dense_solver);

  ParSHUM_solver_factorize(dense_solver);

  ParSHUM_solver_solve(dense_solver, sol_dense);

  ParSHUM_solver_compute_norms(dense_solver, X, sol_dense);

  ParSHUM_solver_finalize(dense_solver);
 
  ParSHUM_vector_destroy(X);
  ParSHUM_vector_destroy(sol_ParSHUM);
  ParSHUM_vector_destroy(sol_dense);
  ParSHUM_solver_destroy(dense_solver);
  ParSHUM_solver_destroy(self);
  
  return 0;
}

void 
ParSHUM_solver_parse_args(ParSHUM_solver self, int argc, char **argv, int exit_on_notFound)
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
	ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__, "value tolerance should be between 0 and 1");
      self->exe_parms->value_tol = tmp;
      continue;
    } else if (!strcmp(argv[i], "--marko_tol")) {
      double tmp = atof(argv[++i]);
      if ( tmp < 1.0 ) 
	ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__, "marko tolerance should be larger then 1");
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
	ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__, "nb_candidates_per_blocks should be at least 1");
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
	ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__,"unable to open output file, the program will write on stdout instead");	  
	file = stdout;
      }
      continue;
    } else if (!strcmp(argv[i], "--nb_previous_pivots")) {
      int tmp = atoi(argv[++i]);
      if (tmp < 1)
	ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__, "nb_previous_pivots should be at least 1");
      self->exe_parms->nb_previous_pivots  = tmp;
      continue;
    } else if (!strcmp(argv[i], "--schur_density_tolerance")) {
      double tmp = atof(argv[++i]);
      if ( tmp > 1.0 ) 
	ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__, "schur density tolerance can not be larger then 1");
      if ( tmp < 0.0 ) 
	ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__, "schur density tolerance can not be smaller then 0");
      self->exe_parms->density_tolerance = tmp;
      continue;
    } else if (!strcmp(argv[i], "--min_pivot_per_steps")) {
      int tmp = atoi(argv[++i]);
      if ( tmp < self->exe_parms->nb_previous_pivots ) 
	ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__, "min_pivot_per_steps should be at least nb_previous_pivots");
      self->exe_parms->min_pivot_per_steps = tmp;
      continue; 
    } else if (!strcmp(argv[i], "--debug_mode")) {
      ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__,"debug mode is not implemented");	  
      continue;
    }  else if (!strcmp(argv[i], "--check_pivots")) {
      self->debug |= ParSHUM_CHECK_PIVOTS;
      continue;
    }  else if (!strcmp(argv[i], "--check_schur_memory")) {
      self->debug |= ParSHUM_CHECK_SCHUR_MEMORY;
      continue;
    }  else if (!strcmp(argv[i], "--check_schur_symetry")) {
      self->debug |= ParSHUM_CHECK_SCHUR_SYMETRY;
      continue;
    }  else if (!strcmp(argv[i], "--check_counters")) {
      self->debug |= ParSHUM_CHECK_COUNTERS;
      continue;
    } else if (!strcmp(argv[i], "--print_each_step")) {
      self->debug |= ParSHUM_DEBUG_VERBOSE_EACH_STEP;
      continue;
    } else if (!strcmp(argv[i], "--verbose_gossip_girl")) {
      self->debug |= ParSHUM_DEBUG_GOSSIP_GIRL;
      continue;
    } else if (!strcmp(argv[i], "--check_schur_doubles")) {
      self->debug |= ParSHUM_CHECK_SCHUR_DOUBLES;
      continue;
    }  else if (!strcmp(argv[i], "--check_ParSHUM_with_plasma_perm")) {
      ParSHUM_solver_dealloc(self);
      exit(check_ParSHUM_with_plasma_perm(argc, argv));
    }  else if (!strcmp(argv[i], "--check_dense_with_ParSHUM_perm")) {
      ParSHUM_solver_dealloc(self);
      exit(check_dense_with_ParSHUM_perm(argc, argv));
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
	ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__, "max_dense_schur should be at least 1");
      self->exe_parms->max_dense_schur = tmp;
      continue;
    } else if (!strcmp(argv[i], "--singeltons_relaxation")) {
      self->exe_parms->singeltons_relaxation = atof(argv[++i]);
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
      if (exit_on_notFound) {
	char mess[2048];
	snprintf(mess, 2048, "unrecognized option \"%s\" ", argv[i]);
	int j = 0;
	while( usageStrign[j] !=  NULL)
	  printf("%s\n", usageStrign[j++]);
	ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__, mess);
      }
    }
  }
  
  if (run_args_start) {
    ParSHUM_parm_type type;
    if ( !strcmp(argv[run_args_start], "value_tol") ) {
      type = ParSHUM_value_tol;
      double init = atof(argv[++run_args_start]), inc = atof(argv[++run_args_start]);
      int nb_steps = atoi(argv[++run_args_start]);
      exit(ParSHUM_solver_run_group(self, type, (void *) &init, nb_steps, (void *) &inc)); 
    } else if ( !strcmp(argv[run_args_start], "marko_tol") ) {
      type = ParSHUM_marko_tol;
      double init = atof(argv[++run_args_start]), inc = atof(argv[++run_args_start]);
      int nb_steps = atoi(argv[++run_args_start]);
      exit(ParSHUM_solver_run_group(self, type, (void *) &init, nb_steps, (void *) &inc)); 
    } else if ( !strcmp(argv[run_args_start], "schur_density") ) {
      type = ParSHUM_schur_density;
      double init = atof(argv[++run_args_start]), inc = atof(argv[++run_args_start]);
      int nb_steps = atoi(argv[++run_args_start]);
      exit(ParSHUM_solver_run_group(self, type, (void *) &init, nb_steps, (void *) &inc)); 
    } else if ( !strcmp(argv[run_args_start], "nb_candidates") ) {
      type = ParSHUM_nb_candidates;
      int init = atoi(argv[++run_args_start]), inc = atoi(argv[++run_args_start]);
      int nb_steps = atoi(argv[++run_args_start]);
      exit(ParSHUM_solver_run_group(self, type, (void *) &init, nb_steps, (void *) &inc)); 
    } else if ( !strcmp(argv[run_args_start], "min_pivots") ) {
      type = ParSHUM_min_pivots;
      int init = atoi(argv[++run_args_start]), inc = atoi(argv[++run_args_start]);
      int nb_steps = atoi(argv[++run_args_start]);
      exit(ParSHUM_solver_run_group(self, type, (void *) &init, nb_steps, (void *) &inc)); 
    } else if ( !strcmp(argv[run_args_start], "nb_threads") ) {
      type = ParSHUM_nb_threads;
      int init = atoi(argv[++run_args_start]), inc = atoi(argv[++run_args_start]);
      int nb_steps = atoi(argv[++run_args_start]);
      exit(ParSHUM_solver_run_group(self, type, (void *) &init, nb_steps, (void *) &inc)); 
    } else  {
      int j = 0;
      while( usageStrign[j] !=  NULL)
	printf("%s\n", usageStrign[j++]);
      ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__, "for the group run, unrecognized type of argument is given" );
    }  
  }
}

void
update_exe_parms(ParSHUM_exe_parms parms, ParSHUM_parm_type type,
		 void *init_val, int step, void *val, void *inc)
{
  double *Dinit, *Dinc, *Dval;
  int *Iinit, *Iinc, *Ival; 

  switch (type) {
  case (ParSHUM_value_tol) :
    Dinit = (double *) init_val; Dinc = (double *) inc; Dval = (double *) val; 
    *Dval = *Dinit * pow(*Dinc, (double) step);
    parms->value_tol = *Dval;
    break;
  case (ParSHUM_marko_tol) :
    Dinit = (double *) init_val; Dinc = (double *) inc;  Dval = (double *) val; 
    *Dval = *Dinit + *Dinc * step;
    parms->marko_tol = *Dval;
    break;
  case (ParSHUM_schur_density) :
    Dinit = (double *) init_val; Dinc = (double *) inc;  Dval = (double *) val; 
    *Dval = *Dinit + *Dinc * step;
    parms->density_tolerance = *Dval;
    break;
  case (ParSHUM_nb_candidates) :
    Iinit = (int *) init_val, Iinc = (int *) inc;  Ival = (int *) val; 
    *Ival = *Iinit + *Iinc * step;
    parms->nb_candidates_per_block = *Ival;
    break;
  case (ParSHUM_min_pivots) :
    Iinit = (int *) init_val, Iinc = (int *) inc; Ival = (int *) val;
    *Ival = *Iinit + *Iinc * step;
    parms->min_pivot_per_steps =  *Ival;
    break;
  case (ParSHUM_nb_threads) :
    Iinit = (int *) init_val, Iinc = (int *) inc; Ival = (int *) val;
    *Ival = *Iinit + *Iinc * step;
    parms->nb_threads =  *Ival;
    break;
  default :
    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__, "unrecognized type of exe_parms");
  }
}

char *
get_outfile_prefix(ParSHUM_exe_parms exe_parms, ParSHUM_parm_type type)
{
  char *self = calloc(PATH_LENGTH, sizeof(*self));
  size_t length = 0;
  *self = '\0';
  
  if (exe_parms->matrix_file)
    snprintf(self, 2058, "%s", basename(exe_parms->matrix_file));

  length = strnlen(self, PATH_LENGTH - length);
  if (type == ParSHUM_value_tol)
    snprintf(self + length, PATH_LENGTH - length,"_MULTIValTol");
  else
    snprintf(self + length, PATH_LENGTH - length,"_%fValTol", exe_parms->value_tol);

  length = strnlen(self, PATH_LENGTH - length);
  if (type == ParSHUM_marko_tol)
    snprintf(self + length, PATH_LENGTH - length,"_MULTIMarkoTol");
  else
    snprintf(self + length, PATH_LENGTH - length,"_%fMarkoTol", exe_parms->marko_tol);
  
  length = strnlen(self, PATH_LENGTH - length);
  if (type == ParSHUM_nb_threads)
    snprintf(self + length, PATH_LENGTH - length,"_MULTIthreads");
  else
    snprintf(self + length, PATH_LENGTH - length,"_%dthreads", exe_parms->nb_threads);
  
  length = strnlen(self, PATH_LENGTH - length);
  if (type == ParSHUM_nb_candidates)
    snprintf(self + length, PATH_LENGTH - length,"_MULTIcandidates");
  else
    snprintf(self + length, PATH_LENGTH - length,"_%dcandidates", exe_parms->nb_candidates_per_block);
  
  length = strnlen(self, PATH_LENGTH - length);
  if (type == ParSHUM_schur_density)
    snprintf(self + length, PATH_LENGTH - length,"_MULTIdensityTol");
  else
    snprintf(self + length, PATH_LENGTH - length,"_%fdensityTol", exe_parms->density_tolerance);
  
  length = strnlen(self, PATH_LENGTH - length);
  if (type == ParSHUM_min_pivots)
    snprintf(self + length, PATH_LENGTH - length,"_MULTIminPivots");
  else
    snprintf(self + length, PATH_LENGTH - length,"_%dminPivots", exe_parms->min_pivot_per_steps);
  length = strnlen(self, PATH_LENGTH - length);
  self[length] = '\0';
  
  return self;
}

int 
ParSHUM_solver_run_group(ParSHUM_solver solver, ParSHUM_parm_type type, 
			 void *init_val, int nb_steps, void *inc)
{
  FILE *file;
  int i;
  ParSHUM_matrix A = ParSHUM_matrix_create();
  ParSHUM_exe_parms exe_parms = solver->exe_parms;
  char *file_ext = strrchr(exe_parms->matrix_file, '.');
  ParSHUM_vector X, rhs;
  char *output_runs_file = get_outfile_prefix(exe_parms, type);
  char filename[PATH_LENGTH];
  double current_val;

  ParSHUM_verbose_create_dirs(solver->verbose->parms->output_dir);
  snprintf(filename, PATH_LENGTH, "%s/data/MULTI_%s_raw.dat", solver->verbose->parms->output_dir, output_runs_file);
 
  file = fopen(filename, "w+");
       
  if (!strcmp(file_ext, ".mtl"))
    ParSHUM_read_mtl_file(A, exe_parms->matrix_file);
#ifdef HAVE_SPRAL
  else  if (!strcmp(file_ext, ".rb"))
      ParSHUM_read_rutherford_boeing(A, exe_parms->matrix_file);
#endif
  else 
    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"unrecognized file type format");
  
  solver->A = A;
  X   = ParSHUM_vector_create(A->n);
  rhs = ParSHUM_vector_create(A->n);

  ParSHUM_vector_read_file(X, solver->exe_parms->RHS_file);
  
  ParSHUM_verbose_print_parms_raw(exe_parms, type, file);
  for( i = 0; i < nb_steps; i++)
    {
      ParSHUM_solver run = ParSHUM_solver_create();
      ParSHUM_matrix matrix = ParSHUM_matrix_create();
      ParSHUM_exe_parms run_exe_parms = malloc(sizeof(*run_exe_parms));

      *run_exe_parms = *exe_parms;
      ParSHUM_matrix_copy(A, matrix);
      
      free(run->exe_parms);
      run->A = matrix;
      run->exe_parms = run->verbose->exe_parms = run_exe_parms;
      free(run->verbose->parms->output_dir);
      run->verbose->parms->output_dir = solver->verbose->parms->output_dir;
      run->verbose->parms->user_out_dir = 1;
      update_exe_parms(run->exe_parms, type, init_val, i, (void *) &current_val, inc);
      
      ParSHUM_solver_init(run);
      ParSHUM_vector_copy(X, rhs);
      ParSHUM_solver_factorize(run);
      ParSHUM_solver_solve(run, rhs);
      
      ParSHUM_solver_compute_norms(run, X, rhs);
      
      ParSHUM_solver_finalize(run);
      ParSHUM_verbose_print_group_run(run->verbose, type, (void *) &current_val, i, file);
      ParSHUM_solver_destroy(run);
    }
  fclose(file);
  
  ParSHUM_vector_destroy(X);
  ParSHUM_vector_destroy(rhs);

  ParSHUM_matrix_destroy(A);
  free(output_runs_file);
  ParSHUM_solver_dealloc(solver);
  
  return 0;
}

/* TODO: FOR NOW WE ASUME THAT N == M */
void
ParSHUM_solver_alloc_counters(ParSHUM_solver solver, int **col_count, int **row_count)
{
  int i, total;
  ParSHUM_counters counter = NULL;
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

  ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"this should not happened!");
}

void
ParSHUM_solver_dealloc_counters(ParSHUM_solver solver, int *col_count, int *row_count)
{
  pthread_mutex_lock(&solver->counters_lock);
  int i;
  ParSHUM_counters counter = NULL;
  long diff;
  for( i = 0; i < solver->nb_counters; i++) {
    diff = col_count - solver->counters[i]->array;
    if ( diff >= 0  && diff < solver->A->n * 2 * solver->size_counters){
      counter = solver->counters[i];
      break;
    }
  }
  if (!counter) 
    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"puffff");

  counter->used_counters[ diff / ( solver->A->n * 2) ] = 0;

  counter->nb_used_counters--;
  pthread_mutex_unlock(&solver->counters_lock);
}

void
ParSHUM_solver_alloc_internal(ParSHUM_solver self) 
{
  double total_extra_space = 1 + self->exe_parms->extra_space_inbetween + self->exe_parms->extra_space;
  int n = self->A->n, m = self->A->m, i;
  int needed_pivots = n < m ? n : m;
  int larger_size = n > m ? n : m;
 
  self->S    = ParSHUM_schur_matrix_create();
  self->D    = ParSHUM_matrix_create();

  if (self->debug & (ParSHUM_DEBUG_VERBOSE_EACH_STEP | ParSHUM_DEBUG_GOSSIP_GIRL))
      ParSHUM_matrix_print(self->A, "A on input");
  ParSHUM_schur_matrix_allocate(self->S, n, m, self->A->nnz,
				self->debug, self->verbose, self->exe_parms->nb_threads,
				self->exe_parms->extra_space, self->exe_parms->extra_space_inbetween);
  ParSHUM_schur_matrix_copy(self->A, self->S, self->exe_parms->value_tol);

  if (self->debug & (ParSHUM_DEBUG_VERBOSE_EACH_STEP | ParSHUM_DEBUG_GOSSIP_GIRL))
      ParSHUM_schur_matrix_print(self->S, "S on input");
  ParSHUM_matrix_allocate(self->D, needed_pivots, 0, 0, 1.0, ParSHUM_Diag_matrix);
  self->L = ParSHUM_L_matrix_create(needed_pivots);
  self->U = ParSHUM_U_matrix_create(self->A, total_extra_space);
  
  self->row_perm      = malloc((size_t) m * sizeof(*self->row_perm));
  self->invr_row_perm = malloc((size_t) m * sizeof(*self->invr_row_perm));
  self->col_perm      = malloc((size_t) n * sizeof(*self->col_perm));
  self->invr_col_perm = malloc((size_t) n * sizeof(*self->invr_col_perm));
  
  int_array_memset(self->row_perm,      ParSHUM_UNUSED_PIVOT, m);
  int_array_memset(self->invr_row_perm, ParSHUM_UNUSED_PIVOT, m);
  int_array_memset(self->col_perm,      ParSHUM_UNUSED_PIVOT, n);
  int_array_memset(self->invr_col_perm, ParSHUM_UNUSED_PIVOT, n);

  self->previous_pivots = malloc((size_t) self->exe_parms->nb_previous_pivots *  sizeof(*self->previous_pivots));
  int_array_memset(self->previous_pivots, INT_MAX / self->exe_parms->nb_previous_pivots, self->exe_parms->nb_previous_pivots);
  /* TODO: check if this is correct for the TP algo */
  self->candidates = calloc(1, sizeof(*self->candidates));
  self->candidates->row        = malloc((size_t) n * sizeof(*self->candidates->row));
  self->candidates->marko      = malloc((size_t) n * sizeof(*self->candidates->marko));
  self->candidates->best_marko = malloc((size_t) self->exe_parms->nb_threads * sizeof(*self->candidates->best_marko));
 
  self->counters              = calloc( 1, sizeof(*self->counters));
  *self->counters             = calloc( 1, sizeof(**self->counters));

  self->random_col = create_randomize(n);

  /* TODO: check if this is correct for the TP algo */
  (*self->counters)->array    = calloc((size_t) self->size_counters * 2 * n,  sizeof(*self->counters[0]->array));
  (*self->counters)->used_counters = calloc((size_t) self->size_counters, sizeof(*self->counters[0]->used_counters));
  self->nb_counters = 1;

  self->allocated_U_struct = n / 100;  
  self->allocated_U_struct = self->allocated_U_struct ? self->allocated_U_struct : 1;
  self->U_struct = calloc((size_t) self->allocated_U_struct, sizeof(*self->U_struct));
  self->allocated_L_struct = m / 100;  
  self->allocated_L_struct = self->allocated_L_struct ? self->allocated_L_struct : 1;
  self->L_struct = calloc((size_t) self->allocated_L_struct, sizeof(*self->L_struct));

  self->Luby = ParSHUM_Luby_create(self->S);
  
  self->cols = malloc((size_t) n * sizeof(*self->cols));
  for( i = 0; i < n; i++)
    self->cols[i] = i;
  self->rows = malloc((size_t) m * sizeof(*self->rows));
  for( i = 0; i < m; i++)
    self->rows[i] = i;
  self->distributions = malloc((size_t) (self->exe_parms->nb_threads + 1)  * sizeof(*self->distributions));
  
  /* srand(time(NULL)); */
  self->seeds = malloc((size_t) self->exe_parms->nb_threads * sizeof(*self->seeds));
  for(i = 0; i < self->exe_parms->nb_threads; i++)
    /* TODO: Put an argument "determenistic" and call srand before if activated */
    self->seeds[i] = rand();

  self->workspace = malloc((size_t) self->exe_parms->nb_threads * sizeof(*self->workspace));
  for(i = 0; i < self->exe_parms->nb_threads; i++)
    self->workspace[i] = malloc((size_t) larger_size * (sizeof(int) + sizeof(double)));
  
  self->logical_rows = calloc((size_t) m, sizeof(*self->logical_rows));
  self->logical_cols = calloc((size_t) n, sizeof(*self->logical_cols));
  
  pthread_mutex_init(&self->counters_lock, NULL); 
}

void
ParSHUM_solver_init(ParSHUM_solver self)
{
  self->verbose->n = self->A->n;
  self->verbose->m = self->A->m;
  self->verbose->nnz_input = self->verbose->nnz_final = self->A->nnz;
  self->verbose->Luby = self->exe_parms->luby_algo;
  
  self->verbose->parms->outfiles_prefix = get_outfile_prefix(self->exe_parms, ParSHUM_parm_none);

  if (self->exe_parms->trace) 
    self->verbose->paje = ParSHUM_paje_create(self->exe_parms->nb_threads);

  if (!is_plasma_init) {
    plasma_init(1);
    is_plasma_init = 1;
  }
  
  if ( !(self->debug & ParSHUM_CHECK_DENSE_W_ParSHUM_PERM) )
    ParSHUM_solver_alloc_internal(self);

  ParSHUM_verbose_create_dirs(self->verbose->parms->output_dir);
}

void
ParSHUM_solver_read_matrix(ParSHUM_solver self)
{
  char *file_ext = strrchr(self->exe_parms->matrix_file, '.');
  self->A        = ParSHUM_matrix_create();
  
  if (!strcmp(file_ext, ".mtl")) 
    ParSHUM_read_mtl_file(self->A, self->exe_parms->matrix_file);
#ifdef HAVE_SPRAL
  else if (!strcmp(file_ext, ".rb"))
    ParSHUM_read_rutherford_boeing(self->A, self->exe_parms->matrix_file);
#endif
  else 
    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"unsupported matrix file");
}

void 
ParSHUM_solver_get_pivots(ParSHUM_solver self, ParSHUM_pivot_set set)
{
  ParSHUM_pivot_cell cells = set->cells;
  int n = self->A->n, i, k, nnz ;

  int new_pivots = 0, old_pivots = self->found_pivots + self->nb_col_singletons + self->nb_row_singletons;
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
    if (self->invr_col_perm[i] != ParSHUM_UNUSED_PIVOT)  {
      if (set->cols_count[i] > set->base) {
	ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, "col is supposed to be a pivot but is larger then base ");
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
    if (self->invr_row_perm[i] != ParSHUM_UNUSED_PIVOT)  {
      if (set->rows_count[i] > set->base) {
	ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, "row is supposed to be a pivot but is larger then base ");
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
  ParSHUM_verbose_update_pivots(self->verbose, new_pivots);
}

void 
ParSHUM_solver_get_Luby_pivots(ParSHUM_solver self, ParSHUM_Luby Luby, int new_Luby_pivots)
{
  ParSHUM_schur_matrix S = self->S;
  ParSHUM_verbose verbose = self->verbose;

  int n = S->n, nb_cols = S->n - self->done_pivots, nb_rows = S->m - self->done_pivots;
  int done_pivots = self->done_pivots;
  int new_pivots =  new_Luby_pivots + self->nb_col_singletons + self->nb_row_singletons;
  int all_pivots = new_pivots + done_pivots;
  int i, base = self->step + 1;
  int nb_threads = self->exe_parms->nb_threads;
  int distribution_perms[nb_threads+1], distribution_n[nb_threads+1], distribution_m[nb_threads+1];
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

  *distribution_m = 0;
  for(i = 1; i < nb_threads; i++)
    distribution_m[i] = (nb_rows / nb_threads) * i;
  distribution_m[nb_threads] = nb_rows;

#pragma omp parallel num_threads(nb_threads) shared(distribution_perms, distribution_n, base, nb_cols, col_sizes, row_sizes, n, verbose)
  {
  ParSHUM_verbose_trace_start_event(verbose, ParSHUM_GETTING_PIVOTS);
  int i, j;
  int me = omp_get_thread_num();
  int start = distribution_perms[me] ;
  int end   = distribution_perms[me+1];
  int *my_col_perms = (int *) self->workspace[me];
  int *my_row_perms = &my_col_perms[nb_cols];
  int my_nb_cols = 0, my_nb_rows = 0;
  
  /* Updte the invr_*_perms arrays */
  for( i = start; i < end; i++)
    {
      invr_col_perms[col_perms[i]] = i;
      invr_row_perms[row_perms[i]] = i;
    }

#pragma omp barrier

  /* Updte the logical arrays with the  non-singelton pivots */
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
  /* 
     Sequential updte the logical arrays with the col_singeltons. 
     This is done seq  beacuse the update in parallel will need to use atomics. 
     We need to treat the col singeltons separately, because we need to indetify 
     entries that are in the pivotal square matrix in order to put them in U. 
  */
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

	if (perm != ParSHUM_UNUSED_PIVOT && perm > end) 
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
  end   = distribution_n[me+1];
  for( i = start; i < end; i++)  {
    int col = cols[i];
    if(logical_cols[col] == base && invr_col_perms[col] == ParSHUM_UNUSED_PIVOT) 
      my_col_perms[my_nb_cols++] = col;
  }

  start = distribution_m[me];
  end   = distribution_m[me+1];
  for( i = start; i < end; i++)  {
    int row = rows[i];
    if(logical_rows[row] == base && invr_row_perms[row] == ParSHUM_UNUSED_PIVOT) 
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
  ParSHUM_verbose_trace_stop_event(verbose);
  }

  self->found_pivots = all_pivots;
  ParSHUM_verbose_update_pivots(self->verbose, new_pivots);
}

void
ParSHUM_solver_check_counters(ParSHUM_solver self)
{
  int i; 
   
  for (i = 0; i < self->nb_counters; i++) 
    ParSHUM_check_counters(i, self->counters[i]->array, self->counters[i]->used_counters,
		      self->done_pivots + 1, self->size_counters, self->A->n);
}

void
ParSHUM_solver_find_pivot_set(ParSHUM_solver self)
{
  ParSHUM_pivot_list list;
  ParSHUM_exe_parms exe_parms = self->exe_parms;
  ParSHUM_verbose_per_step step = ParSHUM_verbose_get_step(self->verbose);
  ParSHUM_verbose verbose = self->verbose;
  int nb_threads = self->exe_parms->nb_threads;
  int needed_pivots = self->A->n < self->A->m ? self->A->n : self->A->m;
  needed_pivots -= self->done_pivots;

  int new_pivots = 0;
  ParSHUM_verbose_trace_start_event(verbose, ParSHUM_GET_SINGELTONS);
  ParSHUM_schur_get_singletons(self->S, self->done_pivots, self->previous_step_pivots,
			       self->exe_parms->value_tol * self->exe_parms->singeltons_relaxation,
			       &self->nb_col_singletons, &self->nb_row_singletons, self->cols,
			       self->rows, self->distributions, self->done_pivots, self->col_perm,
			       self->row_perm, self->invr_col_perm, self->invr_row_perm, self->workspace);
  needed_pivots -= self->nb_col_singletons + self->nb_row_singletons;
  ParSHUM_verbose_trace_stop_event(verbose);

  if (exe_parms->luby_algo) {
  if (needed_pivots > 0) {
    int  best_marko, nb_cols = self->S->n - self->done_pivots;
    int *distributions = self->distributions, best_markos[nb_threads];
    int candidates[nb_threads];
    int i;
    
    *distributions = 0;
    for (i = 1; i < nb_threads; i++) 
      distributions[i] = (nb_cols / nb_threads) * i;
    distributions[nb_threads] = nb_cols;

#pragma omp parallel num_threads(nb_threads) shared(nb_cols, best_marko, distributions)
    {
    ParSHUM_verbose_trace_start_event(verbose, ParSHUM_GET_ELIGEBLE);
    int me = omp_get_thread_num();
    int *my_col_perms = (int *) self->workspace[me];
    int *my_row_perms = &my_col_perms[nb_cols];
    int *global_col_perms = (int *) self->workspace[0];
    int *global_row_perms = &global_col_perms[nb_cols];
    int max_col_length = self->S->nnz /nb_cols ;

    ParSHUM_verbose_start_timing(&step->timing_extracting_candidates);
    
    best_markos[me] = ParSHUM_Luby_get_eligible(self->S, self->Luby, exe_parms->value_tol, self->invr_col_perm, self->invr_row_perm, self->cols, distributions[me], distributions[me + 1], max_col_length );
    ParSHUM_verbose_trace_stop_event(verbose);
#pragma omp barrier    
#pragma omp single 
    {
    ParSHUM_verbose_trace_start_event(verbose, ParSHUM_AUXILIARY);
    best_marko = *best_markos;
    for(i = 1; i < nb_threads; i++)
      best_marko = best_marko > best_markos[i] ? best_markos[i] : best_marko;
    best_marko = !best_marko ? 1 : best_marko;
    best_marko *= exe_parms->marko_tol;
    ParSHUM_verbose_trace_stop_event(verbose);
    }
#pragma omp barrier    
    ParSHUM_verbose_trace_start_event(verbose, ParSHUM_ASSIGN_SCORE);
    candidates[me] = ParSHUM_Luby_assign_score(self->Luby, self->S, best_marko, &self->seeds[me], my_col_perms, my_row_perms, self->cols, distributions[me], distributions[me + 1]);
    ParSHUM_verbose_trace_stop_event(verbose);
#pragma omp barrier
    ParSHUM_verbose_trace_start_event(verbose, ParSHUM_AUXILIARY);
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
    ParSHUM_verbose_trace_stop_event(verbose);
    }
    
    ParSHUM_verbose_stop_timing(&step->timing_extracting_candidates);

    ParSHUM_verbose_start_timing(&step->timing_merging_pivots);
    /*
      This restrictoin most be done for the recantular (overderminated)
      case, so that we do not find too many pivots
    */
    if (step->nb_candidates > needed_pivots)  
      step->nb_candidates = needed_pivots;
    *distributions = 0;
    for (i = 1; i < nb_threads; i++) 
      distributions[i] = (step->nb_candidates / nb_threads) * i;
    distributions[nb_threads] = step->nb_candidates;
    
#pragma omp parallel num_threads(nb_threads) shared(nb_cols)
    {
    ParSHUM_verbose_trace_start_event(verbose, ParSHUM_FIRST_PASS);
    int me = omp_get_thread_num();
    int my_size = distributions[me + 1] - distributions[me];
    int *tmp = (int *) self->workspace[0];
    int *my_col_perms = &tmp[distributions[me]];
    int *my_row_perms = &tmp[nb_cols + distributions[me]];
    int *global_col_perms = &self->col_perm[self->done_pivots + self->nb_col_singletons + self->nb_row_singletons];
    int *global_row_perms = &self->row_perm[self->done_pivots + self->nb_col_singletons + self->nb_row_singletons];

    ParSHUM_Luby_first_pass(self->Luby, self->S, my_col_perms, my_row_perms, my_size);
    ParSHUM_verbose_trace_stop_event(verbose);

#pragma omp barrier
    ParSHUM_verbose_trace_start_event(verbose, ParSHUM_DISCARDING_PIVOTS);
    candidates[me] = ParSHUM_Luby_second_pass(self->S, self->Luby, my_col_perms, my_row_perms, my_size);
    ParSHUM_verbose_trace_stop_event(verbose);
#pragma omp barrier
    ParSHUM_verbose_trace_start_event(verbose, ParSHUM_AUXILIARY);
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
    ParSHUM_verbose_trace_stop_event(verbose);
    }
    new_pivots = candidates[nb_threads-1];
    
  } else { // needed_pivots
    ParSHUM_verbose_start_timing(&step->timing_merging_pivots);
  }
  
  ParSHUM_solver_get_Luby_pivots(self, self->Luby, new_pivots);
  ParSHUM_verbose_stop_timing(&step->timing_merging_pivots);
  if (self->debug & (ParSHUM_DEBUG_VERBOSE_EACH_STEP | ParSHUM_DEBUG_GOSSIP_GIRL)) {
    char mess[2048];
    snprintf(mess, 2048,"%d row singeltons, %d col_singeltons and %d luby pivots were found\ncol pemrs",
    	     self->nb_row_singletons, self->nb_col_singletons, new_pivots);
    print_int_array(self->col_perm, self->A->n, mess);
    print_int_array(self->row_perm, self->A->m, "row_perms");
    print_int_array(self->invr_col_perm, self->A->m, "invr_col_perms");
    print_int_array(self->invr_row_perm, self->A->m, "invr_row_perms");
  }

  } else {
    
    if (self->debug & ParSHUM_CHECK_COUNTERS )
      ParSHUM_solver_check_counters(self);
    
    ParSHUM_verbose_start_timing(&step->timing_extracting_candidates);
    list = get_possible_pivots(self, self->S, self->random_col, self->candidates, nb_threads,
			       exe_parms->value_tol, exe_parms->marko_tol,
			       exe_parms->nb_candidates_per_block);
    ParSHUM_verbose_stop_timing(&step->timing_extracting_candidates);
    
    if( !list->nb_elem )
      ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__, "no possible pivot were found");
    
    ParSHUM_verbose_start_timing(&step->timing_merging_pivots);
    list = merge_pivot_sets(list, self->S, self);
    
    if (list->nb_elem > 1 )
      ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, "not unique set");
    
    ParSHUM_solver_get_pivots(self, list->sets);
    
    if (self->debug & ParSHUM_CHECK_COUNTERS )
      ParSHUM_check_current_counters(self->S, &self->col_perm[self->done_pivots], &self->row_perm[self->done_pivots],
				self->found_pivots - self->done_pivots, 
				self->cols_count, self->rows_count, self->done_pivots + 1);
    
    if (self->debug & (ParSHUM_DEBUG_VERBOSE_EACH_STEP | ParSHUM_DEBUG_GOSSIP_GIRL))
      print_pivot_list(list, "found pivots");
    ParSHUM_pivot_list_destroy(list, self);
    
    ParSHUM_verbose_stop_timing(&step->timing_merging_pivots);
  }
}

void
ParSHUM_solver_update_matrix(ParSHUM_solver self)
{
  ParSHUM_L_matrix L = self->L;
  ParSHUM_matrix D = self->D;
  ParSHUM_U_matrix U = self->U;
  ParSHUM_schur_matrix S = self->S;
  int nb_pivots = self->found_pivots - self->done_pivots;
  ParSHUM_verbose_per_step step = ParSHUM_verbose_get_step(self->verbose);

  ParSHUM_verbose_start_timing(&step->timing_update_LD);
  ParSHUM_schur_matrix_update_LD(S, L, U, D, &self->row_perm[self->done_pivots],
			    &self->col_perm[self->done_pivots], nb_pivots,
			    self->invr_row_perm, self->nb_row_singletons,
			    self->nb_col_singletons, self->workspace);
  ParSHUM_verbose_stop_timing(&step->timing_update_LD);
  
  ParSHUM_verbose_start_timing(&step->timing_update_S);
  ParSHUM_schur_matrix_update_S(S, L, U, &self->col_perm[self->found_pivots],
			   self->n_U_structs, self->invr_row_perm, nb_pivots,
			   &self->row_perm[self->done_pivots], self->workspace,
			   self->exe_parms->value_tol);
  ParSHUM_verbose_stop_timing(&step->timing_update_S);
  
  ParSHUM_verbose_start_timing(&step->timing_update_U);
  ParSHUM_schur_matrix_update_S_rows(S, &self->row_perm[self->found_pivots],
				self->n_L_structs, self->nnz_L_structs,
				self->invr_col_perm, nb_pivots, self->row_perm,
				self->done_pivots, self->workspace);
  ParSHUM_verbose_stop_timing(&step->timing_update_U);

  self->done_pivots = self->found_pivots;
  self->step++;
    
  if (self->debug & ParSHUM_CHECK_SCHUR_SYMETRY )
    ParSHUM_schur_matrix_check_symetry(self->S);

  if (self->debug &  ParSHUM_CHECK_PIVOTS)
    ParSHUM_schur_matrix_check_pivots(self->S, self->row_perm, self->col_perm, 
				      self->invr_row_perm, self->invr_col_perm,
				      self->done_pivots);
  
  if (self->debug & (ParSHUM_DEBUG_VERBOSE_EACH_STEP | ParSHUM_DEBUG_GOSSIP_GIRL)) {
    ParSHUM_schur_matrix_print(S, "S after update");
    ParSHUM_L_matrix_print(L, "L after update");
    ParSHUM_matrix_print(D, "D after update");
    ParSHUM_U_matrix_print(U, "U after update");
  }
  
  if (self->debug & ParSHUM_CHECK_SCHUR_DOUBLES) { 
    ParSHUM_schur_check_doubles(S);
  }
}

int
ParSHUM_continue_pivot_search(ParSHUM_schur_matrix S,
			      int nb_done_pivots,
			      int nb_needed_pivots,
			      int *previous_pivots,
			      int nb_previous_pivots,
			      double density_tolerance,
			      int min_pivot_per_steps,
			      int max_dense_schur,
			      int debug,
			      ParSHUM_verbose verbose)
{
  long n_schur, m_schur, i, sum_pivots;
  int retval = 1;

  if ( debug & ParSHUM_CHECK_ParSHUM_W_PLASMA_PERM ||
       debug & ParSHUM_CHECK_DENSE_W_ParSHUM_PERM ) {
    verbose->reason = ParSHUM_reason_because;
    return 0;
  }

  n_schur = (long) S->n - nb_done_pivots;
  m_schur = (long) S->m - nb_done_pivots;
  if ( S->nnz > (long) n_schur * m_schur * density_tolerance) {
    verbose->reason = ParSHUM_reason_density;
    retval = 0;
  }

  if (retval)
    {
      for( i = 0, sum_pivots = 0; i < nb_previous_pivots; i++)
	sum_pivots += previous_pivots[i];
      if (sum_pivots < min_pivot_per_steps) {
	verbose->reason = ParSHUM_reason_no_pivots;
	retval = 0;
      }
    }

  if (nb_done_pivots >= nb_needed_pivots && retval) {
    verbose->reason = ParSHUM_reason_no_pivots;
    retval = 0;
  }
  
  if ( !retval && max_dense_schur < (nb_needed_pivots - nb_done_pivots) ) {
    verbose->reason |= ParSHUM_reason_dense_too_large;
  }
  return retval;
}

void
ParSHUM_solver_factorize(ParSHUM_solver self)
{
  ParSHUM_exe_parms exe_parms = self->exe_parms;
  ParSHUM_verbose verbose = self->verbose;
  int n = self->S->n, m = self->S->m, i;
  int needed_pivots = n < m ? n : m;
  int *previous_pivots   = self->previous_pivots;
  int nb_previous_pivots = exe_parms->nb_previous_pivots;
  int nb_pivot_blocks = 0;

  if (self->debug &  ParSHUM_CHECK_SCHUR_DOUBLES)
    ParSHUM_schur_check_doubles(self->S);  

  ParSHUM_verbose_start_timing(&verbose->timing_facto);
  ParSHUM_verbose_start_timing(&verbose->timing_facto_sparse);
  while ( ParSHUM_continue_pivot_search(self->S, self->done_pivots, needed_pivots, 
					previous_pivots, nb_previous_pivots,
					exe_parms->density_tolerance,
					exe_parms->min_pivot_per_steps,
					exe_parms->max_dense_schur,
					self->debug, verbose) )
    { 
       ParSHUM_verbose_per_step  step = ParSHUM_verbose_step_start(verbose);
       ParSHUM_verbose_start_timing(&step->timing_step);
       ParSHUM_verbose_start_timing(&step->timing_pivot_search);
       ParSHUM_solver_find_pivot_set(self);
       previous_pivots[nb_pivot_blocks++ % nb_previous_pivots] = self->found_pivots - self->done_pivots;
       ParSHUM_verbose_stop_timing(&step->timing_pivot_search);

       ParSHUM_verbose_start_timing(&step->timing_apply_perms);
       ParSHUM_solver_update_matrix(self);
       ParSHUM_verbose_stop_timing(&step->timing_apply_perms);
       ParSHUM_verbose_stop_timing(&step->timing_step);
    }
  ParSHUM_verbose_stop_timing(&verbose->timing_facto_sparse);

  if (verbose->reason & ParSHUM_reason_dense_too_large)
    return;

  ParSHUM_verbose_start_timing(&verbose->timing_convert_schur);
  ParSHUM_verbose_trace_start_event(verbose, ParSHUM_CONVERT_MATRIX);
  if ( self->debug & ParSHUM_CHECK_DENSE_W_ParSHUM_PERM )  {
    self->A_debug = ParSHUM_dense_2D_permute(ParSHUM_dense_2D_convert_sparse(self->A), 
					self->row_perm, self->col_perm);
    verbose->schur_density = 1.00;
  }  else {
    self->S_dense = ParSHUM_schur_matrix_convert(self->S, self->done_pivots, 
						 self->col_perm, self->invr_col_perm,
						 self->row_perm, self->invr_row_perm);
    verbose->schur_density = 
      (double) self->S->nnz / ((n - self->done_pivots) * (m - self->done_pivots));

    if (self->debug & (ParSHUM_DEBUG_VERBOSE_EACH_STEP | ParSHUM_DEBUG_GOSSIP_GIRL)) 
      ParSHUM_dense_matrix_print(self->S_dense, "dense schur after conversion");
  }
  ParSHUM_verbose_trace_stop_event(verbose);
  ParSHUM_verbose_stop_timing(&verbose->timing_convert_schur);

  ParSHUM_verbose_start_timing(&verbose->timing_facto_dense);
  ParSHUM_verbose_trace_start_event(verbose, ParSHUM_DENSE_FACTORIZATION);
  if (self->debug & ParSHUM_CHECK_DENSE_W_ParSHUM_PERM)  {
    ParSHUM_dense_2D_facto(self->A_debug);
  } else {
    ParSHUM_dense_matrix_factorize(self->S_dense, exe_parms->nb_threads);

    self->dense_pivots = needed_pivots - self->done_pivots;
    ParSHUM_verbose_update_dense_pivots(verbose, self->dense_pivots);
    verbose->nnz_final   = self->L->nnz + self->U->nnz + self->D->n + self->S_dense->n * self->S_dense->m;
    verbose->nnz_L       = self->L->nnz  + ( ( self->S_dense->n * self->S_dense->m - self->S_dense->m) / 2 )  ;
    verbose->nnz_U       = self->U->nnz + self->D->n + (( self->S_dense->n * self->S_dense->m - self->S_dense->m) / 2 + self->S_dense->m ) ;
    verbose->nnz_S_dense = self->S_dense->n * self->S_dense->m;
  }

  ParSHUM_verbose_trace_stop_event(verbose);
  ParSHUM_verbose_stop_timing(&verbose->timing_facto_dense);
  ParSHUM_verbose_stop_timing(&verbose->timing_facto);
}

void
ParSHUM_solver_solve(ParSHUM_solver self, ParSHUM_vector RHS)
{
  if (self->verbose->reason & ParSHUM_reason_dense_too_large)
    return;
  double *RHS_vals          = RHS->vect;
  ParSHUM_verbose verbose        = self->verbose;


  ParSHUM_verbose_start_timing(&verbose->timing_solve);
  if (self->debug & ParSHUM_CHECK_DENSE_W_ParSHUM_PERM ) {
    ParSHUM_vector_permute(RHS, self->row_perm, self->A_debug->m);
    ParSHUM_dense_2D_solve(self->A_debug, RHS);
    ParSHUM_vector_permute(RHS, self->invr_col_perm, self->A_debug->m);
  } else {
    ParSHUM_verbose_start_timing(&verbose->timing_solve_L);

    if (self->debug & ParSHUM_DEBUG_GOSSIP_GIRL)
      ParSHUM_vector_print(RHS, "Init RHS");
    ParSHUM_L_matrix_solve(self->L, RHS, self->row_perm);

    if (self->debug & ParSHUM_DEBUG_GOSSIP_GIRL)
      ParSHUM_vector_print(RHS, "after forward solve");
    ParSHUM_verbose_stop_timing(&verbose->timing_solve_L);
    
    ParSHUM_verbose_start_timing(&verbose->timing_solve_dense);
    if (self->S_dense->n && self->S_dense->m) {
      if (self->S_dense->n == self->S_dense->m) {
  	double *dense_RHS = (double *) *self->workspace;
  	ParSHUM_dense_matrix_get_RHS(self->S_dense, dense_RHS, &self->row_perm[self->done_pivots], RHS_vals, ParSHUM_perm_global);

	plasma_dgetrs(self->S_dense->n, 1, self->S_dense->val, self->S_dense->n,
		      self->S_dense->pivots, dense_RHS, self->S_dense->n);
  	ParSHUM_dense_matrix_update_RHS(self->S_dense, dense_RHS, &self->row_perm[self->done_pivots], RHS_vals);
      } else if ( self->S_dense->n < self->S_dense->m )  {
  	int diff_size = self->S_dense->m - self->S_dense->n;
  	double *dense_RHS = (double *) *self->workspace;
	
  	ParSHUM_dense_matrix_get_RHS(self->S_dense, dense_RHS, &self->row_perm[self->done_pivots], RHS_vals, ParSHUM_perm_both);

  	/* print_int_array(self->S_dense->pivots, self->S_dense->m, "plasma pivots"); */
  	plasma_dtrsm(PlasmaLeft, PlasmaLower, PlasmaNoTrans, PlasmaUnit,  self->S_dense->n, 1, 1.0, self->S_dense->val,self->S_dense->m, dense_RHS, self->S_dense->m);

  	plasma_dgemm(PlasmaNoTrans, PlasmaNoTrans, diff_size, 1, self->S_dense->n, -1.0, &self->S_dense->val[self->S_dense->m - diff_size], self->S_dense->m,  dense_RHS, self->S->m, 1.0,  &dense_RHS[self->S_dense->n], self->S->m);
  	plasma_dtrsm(PlasmaLeft, PlasmaUpper, PlasmaNoTrans, PlasmaNonUnit,  self->S_dense->n, 1, 1.0, self->S_dense->val,  self->S_dense->m, dense_RHS, self->S_dense->m);

  	ParSHUM_dense_matrix_update_RHS(self->S_dense, dense_RHS, &self->row_perm[self->done_pivots], RHS_vals);
      } else {
  	ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__, "not implemeted");
      }
    }


    if (self->debug & ParSHUM_DEBUG_GOSSIP_GIRL)
      ParSHUM_vector_print(RHS, "after dense solve");
    ParSHUM_verbose_stop_timing(&verbose->timing_solve_dense);
    
    ParSHUM_verbose_start_timing(&verbose->timing_solve_U);

    ParSHUM_U_matrix_solve(self->U, self->D, RHS, self->col_perm, self->row_perm, self->dense_pivots);
    if (self->debug & (ParSHUM_DEBUG_VERBOSE_EACH_STEP | ParSHUM_DEBUG_GOSSIP_GIRL))
      ParSHUM_vector_print(RHS, "after backward solve");

    ParSHUM_vector_permute(RHS, self->row_perm, self->S->m);
    if (self->debug & ParSHUM_DEBUG_GOSSIP_GIRL)
      ParSHUM_vector_print(RHS, "P RHS");

    ParSHUM_vector_permute(RHS, self->invr_col_perm, self->S->n);
    if (self->debug & (ParSHUM_DEBUG_VERBOSE_EACH_STEP | ParSHUM_DEBUG_GOSSIP_GIRL))
      ParSHUM_vector_print(RHS, "after solve operation");
    ParSHUM_verbose_stop_timing(&verbose->timing_solve_U);
  }
  ParSHUM_verbose_computed_norms(verbose);
  ParSHUM_verbose_stop_timing(&verbose->timing_solve);
}

void
ParSHUM_solver_compute_norms(ParSHUM_solver self,
			     ParSHUM_vector X,
			     ParSHUM_vector rhs)
{
  if (self->verbose->reason & ParSHUM_reason_dense_too_large)
    return;
  double x_norm, A_norm, b_norm;
  ParSHUM_vector r = ParSHUM_vector_create(X->n);

  /* || r = Ax - b || */
  ParSHUM_matrix_SpMV(self->A, X, r);
  ParSHUM_vector_add(r, 1.00, rhs, -1.00, r);
  self->verbose->backward_error = ParSHUM_vector_2norm(r);

  A_norm  = ParSHUM_matrix_get_norm(self->A);
  x_norm  = ParSHUM_vector_2norm(X);
  b_norm  = ParSHUM_vector_2norm(rhs);

  self->verbose->backward_error /= A_norm * x_norm + b_norm;

  ParSHUM_vector_destroy(r);
}

/* void */
/* ParSHUM_solver_iterative_refinement(ParSHUM_solver self,  */
/* 				    ParSHUM_vector X,  */
/* 				    ParSHUM_vector rhs, */
/* 				    double wanted_precision) */
/* { */
/*   if (self->verbose->reason & ParSHUM_reason_dense_too_large) */
/*     return; */
/*   ParSHUM_vector sol = ParSHUM_vector_create(self->A->n); */
/*   ParSHUM_vector tmp = ParSHUM_vector_create(self->A->n); */
/*   int i = 0; */

/*   ParSHUM_vector_copy(rhs, sol); */
/*   ParSHUM_solver_solve(self, sol); */
/*   ParSHUM_solver_copmpute_norms(self, X, sol, rhs); */
/*   printf("iteration %d: backward error = %e and forward error = %e\n",  */
/* 	 i++, */
/* 	 self->verbose->backward_error, */
/* 	 self->verbose->forward_error); */

/*   while ( i < 20 &&  self->verbose->backward_error > wanted_precision) { */
/*     ParSHUM_matrix_SpMV(self->A, sol, tmp); */
/*     ParSHUM_vector_add(rhs, 1.00, tmp, -1.00, tmp); */
/*     ParSHUM_solver_solve(self, tmp); */
/*     ParSHUM_vector_add(sol, 1.00, tmp, 1.00, sol); */
    
/*     ParSHUM_solver_copmpute_norms(self, X, sol, rhs); */
/*     printf("iteration %d: backward error = %e and forward error = %e\n",  */
/* 	   i++, */
/* 	   self->verbose->backward_error, */
/* 	   self->verbose->forward_error); */
/*   } */

/*   ParSHUM_vector_destroy(sol); */
/*   ParSHUM_vector_destroy(tmp); */
/* } */

void 
ParSHUM_solver_finalize(ParSHUM_solver self)
{
  if (is_plasma_init) {
    plasma_finalize();
    is_plasma_init = 0;
  }
  
  ParSHUM_verbose_print(self->verbose);
  ParSHUM_verbose_draw_graph(self->verbose);
}

void
ParSHUM_solver_destroy(ParSHUM_solver self)
{
  int i;

  if (self->debug & ParSHUM_CHECK_DENSE_W_ParSHUM_PERM) {
    ParSHUM_dense_2D_destroy(self->A_debug);
  } else {
    ParSHUM_matrix_destroy(self->A);
    ParSHUM_schur_matrix_destroy(self->S);
    ParSHUM_L_matrix_destroy(self->L);
    ParSHUM_matrix_destroy(self->D);
    ParSHUM_U_matrix_destroy(self->U);
    if(self->S_dense)
      ParSHUM_dense_matrix_destroy(self->S_dense);

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

    ParSHUM_Luby_destroy(self->Luby);
  }
  
  ParSHUM_solver_dealloc(self);
}
