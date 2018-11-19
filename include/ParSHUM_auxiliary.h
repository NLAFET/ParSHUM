#ifndef _ParSHUM_AUXILIARY_H 
#define _ParSHUM_AUXILIARY_H 

#include <ParSHUM_enum.h>
#include <signal.h>

#define GDB_BREAK   /* raise(SIGINT) */

void  ParSHUM_fatal_error(const char *func, char *filename, const int line, const char *msg);
void  ParSHUM_warning    (const char *func, char *filename, const int line, const char *msg);

double get_max_double(double *array, int nb_elem);
void   print_int_array(int *array, int n, char *mess);
void   update_counter(int *counter, int *index, int n, int base);
void   check_vlaid_perms(int *perms, int needed_pivots, int nb_pivots);
void   check_perms_and_invr_perms(int *perms, int *invr_perms, int nb_pivots, char *name);

void   int_array_memset(int *array, int val, int n);
int    *create_randomize(int n);
ParSHUM_overlaps check_overalping_regions(long r1_start, long r1_end, 
				     long r2_start, long r2_end);
void    ParSHUM_check_counters(int counter, int *array, int *used,
			  int base, int size, int n);
int    ParSHUM_rand_int(int *seed, int size);
double ParSHUM_rand_double(int *seed);

#endif // _ParSHUM_AUXILIARY_H 
