#ifndef _TP_AUXILIARY_H 
#define _TP_AUXILIARY_H 


void  TP_fatal_error(const char *func, char *filename, const int line, const char *msg);
void  TP_warning    (const char *func, char *filename, const int line, const char *msg);

double get_max_double(double *array, int nb_elem);
void   print_int(int *array, int n, char *mess);
void   update_counter(int *counter, int *index, int n);
void   check_vlaid_perms(int *perms, int n);
void   int_array_memset(int *array, int val, int n);
int    *create_randomize(int n);


#endif // _TP_AUXILIARY_H 
