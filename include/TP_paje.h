#ifndef _TP_PAJE_H 
#define _TP_PAJE_H 

typedef struct _TP_paje *TP_paje;

#define TP_GET_SINGELTONS       0 
#define TP_GET_ELIGEBLE         1 
#define TP_AUXILIARY            2 
#define TP_ASSIGN_SCORE         3 
#define TP_FIRST_PASS           4 
#define TP_DISCARDING_PIVOTS    5 
#define TP_GETTING_PIVOTS       6 
#define TP_UPDATE_L             7 
#define TP_UPDATE_S_COLS        8 
#define TP_UPDATE_S_ROWS        9 
#define TP_CONVERT_MATRIX      10 
#define TP_DENSE_FACTORIZATION 11 

TP_paje TP_paje_create(int nb_threads);
void    TP_paje_start_event(TP_paje self, int id);
void    TP_paje_stop_event(TP_paje self);
void    TP_paje_create_file(TP_paje self, char *dir, char *prefix);
void    TP_paje_destroy(TP_paje self);

#endif // _TP_PAJE_H 
