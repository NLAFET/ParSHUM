#ifndef _ParSHUM_ENUM_H 
#define  _ParSHUM_ENUM_H 

enum UM_verbosity_level
  {
    ParSHUM_silent, 
    ParSHUM_default, 
    ParSHUM_debug
  };

enum ParSHUM_return_code
  {
    ParSHUM_sucess, 
    ParSHUM_uninitialized_data
  };

typedef  enum _ParSHUM_matrix_type 
  {
    ParSHUM_CSC_matrix,    
    ParSHUM_CSR_matrix,
    ParSHUM_Diag_matrix,
    ParSHUM_Rutherford_matrix
  } ParSHUM_matrix_type;

#define  ParSHUM_reason_unknown              0 
#define  ParSHUM_reason_density          (1<<1)
#define  ParSHUM_reason_no_pivots        (1<<2)
#define  ParSHUM_reason_dense_too_large  (1<<3)
#define  ParSHUM_reason_no_switch        (1<<4)
#define  ParSHUM_reason_because          (1<<5)

typedef enum _ParSHUM_overlaps
  {
    ParSHUM_overlap_none,
    ParSHUM_overlap_begin,
    ParSHUM_overlap_end, 
    ParSHUM_overlap_total
  } ParSHUM_overlaps;

typedef enum _ParSHUM_parm_type
  {
    ParSHUM_parm_none,
    ParSHUM_value_tol,
    ParSHUM_marko_tol,
    ParSHUM_schur_density,
    ParSHUM_nb_candidates,
    ParSHUM_min_pivots,
    ParSHUM_nb_threads
  } ParSHUM_parm_type;

typedef enum _ParSHUM_perm_type
  {
    ParSHUM_perm_none,
    ParSHUM_perm_global,
    ParSHUM_perm_both
  } ParSHUM_perm_type;

#define  ParSHUM_DEBUG_NONE              0
#define  ParSHUM_CHECK_PIVOTS            (1<<0)
#define  ParSHUM_CHECK_SCHUR_MEMORY      (1<<1)
#define  ParSHUM_CHECK_SCHUR_SYMETRY     (1<<2)
#define  ParSHUM_CHECK_COUNTERS          (1<<3)
#define  ParSHUM_DEBUG_VERBOSE_EACH_STEP (1<<4)
#define  ParSHUM_DEBUG_GOSSIP_GIRL       (1<<5)
#define  ParSHUM_DEBUG_GARBAGE_COLLECTOR (1<<6)
#define  ParSHUM_CHECK_SCHUR_DOUBLES     (1<<7)
#define  ParSHUM_CHECK_ParSHUM_W_PLASMA_PERM  (1<<10)
#define  ParSHUM_CHECK_DENSE_W_ParSHUM_PERM   (1<<11)

#define ParSHUM_UNUSED_PIVOT  -100

#define PATH_LENGTH 2048

#endif  // _ParSHUM_ENUM_H 
