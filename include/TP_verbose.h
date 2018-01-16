#ifndef   _TP_VERBOSE_H
#define   _TP_VERBOSE_H

#include <sys/time.h>
#include "TP_matrix.h"

typedef struct _TP_verbose *TP_verbose;
typedef struct _TP_verbose_per_step *TP_verbose_per_step;
typedef struct _TP_verbose_parms *TP_verbose_parms;
typedef struct _TP_exe_parms *TP_exe_parms;
typedef struct _TP_Luby_step *TP_Luby_step; 

struct _TP_verbose_parms {
  char *prog_name;
  char *output_dir;
  int user_out_dir;
  char *outfiles_prefix;
  
  FILE *out_file;
  int user_out_file;
  int verbosity;
};

struct _TP_exe_parms {
  double value_tol;
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
};

struct _TP_verbose_per_step {
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
  int nb_Luby_steps;

  TP_Luby_step Luby_init_phase;

  TP_Luby_step Luby_step_first;
  TP_Luby_step Luby_step_last;

  TP_verbose_per_step next;
};

struct _TP_Luby_step {
  double timing_max;
  double timing_first_pass;
  double timing_second_pass;
  double timing_discarding;

  int nb_candidates;
  int nb_pivots;

  TP_Luby_step next;
};

struct _TP_verbose {
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

  int nb_steps;

  int reason;
  
  int Luby;

  TP_verbose_per_step stats_first;
  TP_verbose_per_step stats_last;

  TP_verbose_parms parms;
  TP_exe_parms exe_parms;
};

/* TODO: this should not be here! make macros for all of them*/
TP_verbose          TP_verbose_create_V0(TP_exe_parms exe_parms);
TP_verbose_per_step TP_verbose_step_start_V0(TP_verbose self);
void                TP_verbose_print_V0(TP_verbose self);
void                TP_verbose_create_dirs_V0(char *dir);
void                TP_verbose_draw_graph_V0(TP_verbose verbose);
void                TP_verbose_destroy_V0(TP_verbose self);
void                TP_verbose_print_parms_raw(TP_exe_parms exe_parms, TP_parm_type type, FILE *file);
void                TP_verbose_print_group_run(TP_verbose verbose, TP_parm_type type, void *val,
					       int current_run, FILE *file);
void                TP_verbose_update_Luby_step_V0(TP_verbose_per_step step, TP_Luby_step Luby_step);

static inline void
TP_verbose_timing_start_V0(double *timing)
{
  struct timeval time; 
  gettimeofday(&time, NULL);
  *timing =  (double) (time.tv_sec * 1e6 + time.tv_usec);
}

static inline void
TP_verbose_timing_stop_V0(double *timing)
{
  struct timeval time; 
  gettimeofday(&time, NULL);
  *timing =  (double) (time.tv_sec * 1e6 + time.tv_usec) - *timing;
}

#ifdef TP_VERBOSITY 
#if TP_VERBOSITY == 1

#define TP_verbose_create(E)                TP_verbose_create_V0(E)
#define TP_verbose_step_start(self)         TP_verbose_step_start_V0(self)
#define TP_verbose_destroy(self)            TP_verbose_destroy_V0(self)
#define TP_verbose_start_timing(time)       TP_verbose_timing_start_V0(time)
#define TP_verbose_stop_timing(time)        TP_verbose_timing_stop_V0(time)
#define TP_verbose_get_step(V)              V->stats_last
#define TP_verbose_update_pivots(V,n)       V->sparse_pivots += n; V->stats_last->nb_pivots = n
#define TP_verbose_print(V)                 TP_verbose_print_V0(V)
#define TP_verbose_update_dense_pivots(V,n) V->dense_pivots = n
#define TP_verbose_create_dirs(S)           TP_verbose_create_dirs_V0(S)
#define TP_verbose_draw_graph(V)            TP_verbose_draw_graph_V0(V)
#define TP_verbose_update_Luby_step(s,l)    TP_verbose_update_Luby_step_V0(s, l)
/* #elif TP_VERBOSITY == 1  */

#endif
#else 

#define TP_verbose_create(E)                NULL
#define TP_verbose_step_start(self)         NULL
#define TP_verbose_destroy(self)            (void) 0
#define TP_verbose_start_timing(time)       (void) 0
#define TP_verbose_stop_timing(time)        (void) 0
#define TP_verbose_get_step(V)              NULL
#define TP_verbose_update_pivots(V,n)       (void) 0
#define TP_verbose_print(V)                 (void) 0
#define TP_verbose_update_dense_pivots(V,n) (void) 0
#define TP_verbose_create_dirs(S)           (void) 0
#define TP_verbose_draw_graph(S)            (void) 0
#define TP_verbose_update_Luby_step(s,l)    (void) 0
#endif 

#endif // _TP_VERBOSE_H
