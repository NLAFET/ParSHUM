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

#define _TP_ENUM_H

#endif //_TP_ENUM_H 
