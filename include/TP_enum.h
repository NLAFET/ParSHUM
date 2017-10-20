#ifndef _TP_ENUM_H 
#define  _TP_ENUM_H 

enum UM_verbosity_level
  {
    TP_silent, 
    TP_default, 
    TP_debug
  };

enum TP_return_code
  {
    TP_sucess, 
    TP_uninitialized_data
  };

typedef  enum _TP_matrix_type 
  {
    TP_CSC_matrix,    
    TP_CSR_matrix,
    TP_Diag_matrix,
    TP_Rutherford_matrix
  } TP_matrix_type;

#define  TP_reason_unknown              0 
#define  TP_reason_density          (1<<1)
#define  TP_reason_no_pivots        (1<<2)
#define  TP_reason_dense_too_large  (1<<3)
#define  TP_reason_because          (1<<4)

typedef enum _TP_overlaps
  {
    TP_overlap_none,
    TP_overlap_begin,
    TP_overlap_end, 
    TP_overlap_total
  } TP_overlaps;

typedef enum _TP_parm_type
  {
    TP_parm_none,
    TP_value_tol,
    TP_marko_tol,
    TP_schur_density,
    TP_nb_candidates,
    TP_min_pivots,
    TP_nb_threads
  } TP_parm_type;

#define  TP_DEBUG_NONE              0
#define  TP_CHECK_PIVOTS            (1<<0)
#define  TP_CHECK_SCHUR_MEMORY      (1<<1)
#define  TP_CHECK_SCHUR_SYMETRY     (1<<2)
#define  TP_CHECK_COUNTERS          (1<<3)
#define  TP_DEBUG_VERBOSE_EACH_STEP (1<<4)
#define  TP_DEBUG_GOSSIP_GIRL       (1<<5)
#define  TP_DEBUG_GARBAGE_COLLECTOR (1<<6)
#define  TP_CHECK_SCHUR_DOUBLES     (1<<7)
#define  TP_CHECK_TP_W_PLASMA_PERM  (1<<10)
#define  TP_CHECK_DENSE_W_TP_PERM   (1<<11)

#define TP_UNUSED_PIVOT  -100

#define PATH_LENGTH 2048

#endif  // _TP_ENUM_H 
