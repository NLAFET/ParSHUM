#ifndef   _ParSHUM_VERBOSE_H
#define   _ParSHUM_VERBOSE_H

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "ParSHUM_matrix.h"
#include "ParSHUM_paje.h"

typedef struct _ParSHUM_verbose *ParSHUM_verbose;
typedef struct _ParSHUM_verbose_per_step *ParSHUM_verbose_per_step;
typedef struct _ParSHUM_verbose_parms *ParSHUM_verbose_parms;
typedef struct _ParSHUM_exe_parms *ParSHUM_exe_parms;

struct _ParSHUM_verbose_parms {
  char *prog_name;
  char *output_dir;
  int user_out_dir;
  char *outfiles_prefix;
  
  FILE *out_file;
  int user_out_file;
  int verbosity;
};

struct _ParSHUM_exe_parms {
  double value_tol;
  double singeltons_relaxation;
  double marko_tol;
  double extra_space;
  double extra_space_inbetween;
  double density_tolerance;

  char *matrix_file;
  char *RHS_file;

  int min_pivot_per_steps;
  int nb_threads;
  int nb_candidates_per_block;
  int nb_previous_pivots;
  int max_dense_schur;
  int luby_algo;
  int trace;
};

struct _ParSHUM_verbose_per_step {
  double timing_step;

  double timing_pivot_search;
  double timing_extracting_candidates;
  double timing_merging_pivots;

  /* total time for updating the Schur with the new pivot set */
  double timing_apply_perms;
  double timing_update_LD;
  double timing_update_U;
  double timing_update_S;

  long new_nnz;
  long nb_flops;
  int nb_pivots;
  int nb_candidates;

  ParSHUM_verbose_per_step next;
};

struct _ParSHUM_verbose {
  double timing_facto;
  double timing_facto_sparse;
  double timing_facto_dense;

  /* time needed to convert the Schur from schur matrix to a dense matrix */
  double timing_convert_schur;
  double timing_total_pivot_search;
  double timing_total_extracting_candidates;
  double timing_total_merging_pivots;

  /* LUBY's algorithm */
  double timing_total_Luby_assign_scores;
  double timing_total_Luby_first_pass;
  double timing_total_Luby_second_pass;

  double timing_solve;
  double timing_solve_L;
  double timing_solve_dense;
  double timing_solve_U;
  double timing_total_update;
  double timing_total_update_LD;
  double timing_total_update_U;
  double timing_total_update_S;

  double forward_error;
  double backward_error;
  
  double schur_density;

  int n;
  int m;
  /* #pivots eliminated by our algorithm */
  int sparse_pivots;
  /* #pivtos eliminated by the dense solver */
  int dense_pivots;
  long nnz_input;
  long nnz_final;
  long nnz_L;
  long nnz_U;
  long nnz_S_dense;
  int computed_norms;

  int nb_steps;

  int reason;
  
  int Luby;

  ParSHUM_verbose_per_step stats_first;
  ParSHUM_verbose_per_step stats_last;

  ParSHUM_paje paje;

  ParSHUM_verbose_parms parms;
  ParSHUM_exe_parms exe_parms;
};

/* TODO: this should not be here! make macros for all of them*/
ParSHUM_verbose          ParSHUM_verbose_create_V0(ParSHUM_exe_parms exe_parms);
ParSHUM_verbose_per_step ParSHUM_verbose_step_start_V0(ParSHUM_verbose self);
void                ParSHUM_verbose_print_V0(ParSHUM_verbose self);
void                ParSHUM_verbose_create_dirs_V0(char *dir);
void                ParSHUM_verbose_draw_graph_V0(ParSHUM_verbose verbose);
void                ParSHUM_verbose_destroy_V0(ParSHUM_verbose self);
void                ParSHUM_verbose_print_parms_raw(ParSHUM_exe_parms exe_parms, ParSHUM_parm_type type, FILE *file);
void                ParSHUM_verbose_print_group_run(ParSHUM_verbose verbose, ParSHUM_parm_type type, void *val,
					       int current_run, FILE *file);
void                ParSHUM_verbose_trace_start_event_V0(ParSHUM_verbose self, int id);
void                ParSHUM_verbose_trace_stop_event_V0(ParSHUM_verbose self);

static inline void
ParSHUM_verbose_timing_start_V0(double *timing)
{
  struct timeval time; 
  gettimeofday(&time, NULL);
  *timing =  (double) (time.tv_sec * 1e6 + time.tv_usec);
}

static inline void
ParSHUM_verbose_timing_stop_V0(double *timing)
{
  struct timeval time; 
  gettimeofday(&time, NULL);
  *timing =  (double) (time.tv_sec * 1e6 + time.tv_usec) - *timing;
}

#ifdef ParSHUM_VERBOSITY 
#if ParSHUM_VERBOSITY == 1

#define ParSHUM_verbose_create(E)                ParSHUM_verbose_create_V0(E)
#define ParSHUM_verbose_step_start(self)         ParSHUM_verbose_step_start_V0(self)
#define ParSHUM_verbose_destroy(self)            ParSHUM_verbose_destroy_V0(self)
#define ParSHUM_verbose_start_timing(time)       ParSHUM_verbose_timing_start_V0(time)
#define ParSHUM_verbose_stop_timing(time)        ParSHUM_verbose_timing_stop_V0(time)
#define ParSHUM_verbose_get_step(V)              V->stats_last
#define ParSHUM_verbose_update_pivots(V,n)       V->sparse_pivots += n; V->stats_last->nb_pivots = n
#define ParSHUM_verbose_computed_norms(V)        V->computed_norms = 1
#define ParSHUM_verbose_print(V)                 ParSHUM_verbose_print_V0(V)
#define ParSHUM_verbose_update_dense_pivots(V,n) V->dense_pivots = n
#define ParSHUM_verbose_create_dirs(S)           ParSHUM_verbose_create_dirs_V0(S)
#define ParSHUM_verbose_draw_graph(V)            ParSHUM_verbose_draw_graph_V0(V)
#define ParSHUM_verbose_trace_start_event(V,id)  ParSHUM_verbose_trace_start_event_V0(V,id)
#define ParSHUM_verbose_trace_stop_event(V)      ParSHUM_verbose_trace_stop_event_V0(V)
/* #elif ParSHUM_VERBOSITY == 1  */

#endif
#else 

#define ParSHUM_verbose_create(E)                NULL
#define ParSHUM_verbose_step_start(self)         NULL
#define ParSHUM_verbose_destroy(self)            (void) 0
#define ParSHUM_verbose_start_timing(time)       (void) 0
#define ParSHUM_verbose_stop_timing(time)        (void) 0
#define ParSHUM_verbose_get_step(V)              NULL
#define ParSHUM_verbose_update_pivots(V,n)       (void) 0
#define ParSHUM_verbose_computed_norms(V)        (void) 0
#define ParSHUM_verbose_print(V)                 (void) 0
#define ParSHUM_verbose_update_dense_pivots(V,n) (void) 0
#define ParSHUM_verbose_create_dirs(S)           (void) 0
#define ParSHUM_verbose_draw_graph(S)            (void) 0
#define ParSHUM_verbose_trace_start_event(V,id)  (void) 0
#define ParSHUM_verbose_trace_stop_event(V)      (void) 0
#endif 

#endif // _ParSHUM_VERBOSE_H
