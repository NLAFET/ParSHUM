#ifndef _ParSHUM_PAJE_H 
#define _ParSHUM_PAJE_H 

typedef struct _ParSHUM_paje *ParSHUM_paje;

#define ParSHUM_GET_SINGELTONS       0 
#define ParSHUM_GET_ELIGEBLE         1 
#define ParSHUM_AUXILIARY            2 
#define ParSHUM_ASSIGN_SCORE         3 
#define ParSHUM_FIRST_PASS           4 
#define ParSHUM_DISCARDING_PIVOTS    5 
#define ParSHUM_GETTING_PIVOTS       6 
#define ParSHUM_UPDATE_L             7 
#define ParSHUM_UPDATE_S_COLS        8 
#define ParSHUM_UPDATE_S_ROWS        9 
#define ParSHUM_CONVERT_MATRIX      10 
#define ParSHUM_DENSE_FACTORIZATION 11 

ParSHUM_paje ParSHUM_paje_create(int nb_threads);
void         ParSHUM_paje_start_event(ParSHUM_paje self, int id);
void         ParSHUM_paje_stop_event(ParSHUM_paje self);
void         ParSHUM_paje_create_file(ParSHUM_paje self, char *dir, char *prefix);
void         ParSHUM_paje_destroy(ParSHUM_paje self);

#endif // _ParSHUM_PAJE_H 
