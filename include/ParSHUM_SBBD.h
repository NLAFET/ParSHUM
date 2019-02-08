#ifndef _ParSHUM_SBBD_H 
#define _ParSHUM_SBBD_H 

#include "ParSHUM_SBBD_util.h"

typedef struct _ParSHUM_SBBD *ParSHUM_SBBD; 

ParSHUM_SBBD ParSHUM_SBBD_create(MPI_Comm world);
  
void  ParSHUM_SBBD_parse_args(ParSHUM_SBBD self, int argc, char **argv);
void  ParSUHM_SBBD_read_matrix(ParSHUM_SBBD self);
void  ParSHUM_SBBD_partition(ParSHUM_SBBD self);
void  ParSHUM_SBBD_factorize(ParSHUM_SBBD self);
void  ParSHUM_SBBD_destroy(ParSHUM_SBBD self);

#endif // _ParSHUM_SBBD_H 
