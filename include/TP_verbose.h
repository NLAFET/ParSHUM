#ifndef   _TP_VERBOSE_H
#define   _TP_VERBOSE_H

#include <sys/time.h>
#include "TP_matrix.h"

typedef struct _TP_verbose *TP_verbose;
typedef struct _TP_verbose_per_step *TP_verbose_per_step;
typedef struct _TP_verbose_parms *TP_verbose_parms;
typedef struct _TP_exe_parms *TP_exe_parms;

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
  int  read_matrix;

  int min_pivot_per_steps;
  int nb_threads;
  int nb_candidates_per_block;
  int nb_previous_pivots;
};


struct _TP_verbose_per_step {
  double timing_step;

  /* Total time for finding a set of independent pivots  */
  double timing_pivot_search;
  /* Time  to extract initial  pivot sets */
  double timing_extracting_candidates;
  /* Time for merging tbe pivot sets */
  double timing_merging_pivots;

  /* total time for updating the Schur with the new pivot set */
  double timing_apply_perms;
  double timing_update_LD;
  double timing_update_U;
  double timing_update_S;

  long new_nnz;
  long nb_flops;
  int nb_pivots;

  TP_verbose_per_step next;
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


  double timing_solve;
  double timing_solve_L;
  double timing_solve_dense;
  double timing_solve_U;
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

  TP_why_KO reason;

  TP_verbose_per_step stats_first;
  TP_verbose_per_step stats_last;

  TP_verbose_parms parms;
  TP_exe_parms exe_parms;
};


TP_verbose          TP_verbose_create_V0(TP_exe_parms exe_parms);
TP_verbose_per_step TP_verbose_step_start_V0(TP_verbose self);
void                TP_verbose_print_V0(TP_verbose self);
void                TP_verbose_create_dirs_V0(char *dir);
void                TP_verbose_draw_graph_V0(TP_verbose verbose);
void                TP_verbose_destroy_V0(TP_verbose self);

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

#endif 

#endif // _TP_VERBOSE_H
