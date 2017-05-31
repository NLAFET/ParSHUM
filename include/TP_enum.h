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

enum TP_matrix_type 
  {
    TP_CSC_matrix,    
    TP_CSR_matrix,
    TP_Diag_matrix,
    TP_Rutherford_matrix
  };

typedef enum _TP_why_KO 
  {
    TP_reason_unknown, 
    TP_density, 
    TP_no_pivots, 
    TP_rectangle,
    TP_because
  } TP_why_KO;


typedef enum _TP_overlaps
  {
    TP_overlap_none,
    TP_overlap_begin, 
    TP_overlap_end, 
    TP_overlap_total
  } TP_overlaps;


#define  TP_DEBUG_NONE              0
#define  TP_CHECK_PIVOTS            (1<<0)
#define  TP_CHECK_SCHUR_MEMORY      (1<<1)
#define  TP_CHECK_SCHUR_SYMETRY     (1<<2)
#define  TP_DEBUG_VERBOSE_EACH_STEP (1<<3)
#define  TP_DEBUG_GOSSIP_GIRL       (1<<4)
#define  TP_DEBUG_GARBAGE_COLLECTOR (1<<4)
#define  TP_CHECK_TP_W_PLASMA_PERM  (1<<10)
#define  TP_CHECK_DENSE_W_TP_PERM   (1<<11)

#define TP_UNUSED_PIVOT     -100

#define PATH_LENGTH 2048  

#endif  // _TP_ENUM_H 
