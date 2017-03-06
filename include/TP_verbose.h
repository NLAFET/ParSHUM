#ifndef   _TP_VERBOSE_H
#define   _TP_VERBOSE_H

#include <sys/time.h>

#include "TP_matrix.h"

typedef struct _TP_verbose *TP_verbose;
typedef struct _per_step *per_step;

struct _per_step {
  double timing_pivots;
  double extracting_candidates;
  double merging_pivots;
  double timing_update;
  double update_LD;
  double update_U;

  long new_nnz;
  long nb_flops;
  int nb_pivots;

  per_step next;
};

struct _TP_verbose {
  double timing_search_pivots;
  double timing_update;

  int n;
  int m;
  long nnz_input;
  long nnz_final;

  per_step stats_first;
  per_step stats_last;

  int nb_steps;
  int verbosity;
};

#ifdef TP_VERBOSITY 
#if TP_VERBOSITY == 0 

#define _verbose_create(V, M)         TP_verbose_create_V0(V, M);
#define _verbose_step_start(self)     TP_verbose_step_start_V0(self);
#define _verbose_destroy(self)        TP_verbose_destroy_V0(self);
#define _verbose_start_timing(time)   TP_verbose_timing_start(time);
#define _verbose_stop_timing(time)    TP_verbose_timing_stop(time);

/* #elif TP_VERBOSITY == 1  */

/* #define  */
/* #define  */
/* #define  */

#endif
#else 

#define DD_verbose_create(V, M)  NULL;
#define DD_verbose_step_start(self)
#define DD_verbose_destroy(self)
#define _verbose_start_timing(time)
#define _verbose_stop_timing(time)

#endif 


TP_verbose TP_verbose_create_V0(int verbosity, TP_matrix matrix);
void       TP_verbose_step_start_V0(TP_verbose self);
void       TP_verbose_destroy_V0(TP_verbose self);

static inline void
TP_verbose_timing_start(double *timing)
{
  struct timeval time; 
  gettimeofday(&time, NULL);
  *timing =  (double) (time.tv_sec * 1e6 + time.tv_usec);
}

static inline void
TP_verbose_timing_stop(double *timing)
{
  struct timeval time; 
  gettimeofday(&time, NULL);
  *timing -=  (double) (time.tv_sec * 1e6 + time.tv_usec);
}

#endif // _TP_VERBOSE_H
