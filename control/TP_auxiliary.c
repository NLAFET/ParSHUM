#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <libgen.h>
#include <string.h>

#include "TP_auxiliary.h"

#define TP_RAND_MAX 32767

void 
TP_fatal_error(const char *func, char *filename, const int line, const char *msg)
{
  printf("TP_FATAL ERROR in %s(%s:%d) \n%s\n", func, basename(filename), line, msg);
  GDB_BREAK;
  abort();
}

void 
TP_warning(const char *func, char *filename, const int line, const char *msg)
{
  printf("TP_warning in %s(%s:%d) \n%s\n", func, basename(filename), line, msg);
  GDB_BREAK;
}

double 
get_max_double(double *array, int nb_elem)
{
  if (nb_elem == 0) 
    return 0.0;
  double max = fabs(array[0]);
  int i;
  
  for(i = 1; i < nb_elem; i++)
    if(fabs(array[i]) > max)
      max  = fabs(array[i]);

  return max;
}

void 
print_int_array(int *array, int n, char *mess)
{
  int i;
  printf("%s\n", mess);
  for(i = 0; i < n; i++)
    printf("%d   ", array[i]);
  printf("\n");
}

void
update_counter(int *counter, int *index, int n, int base)
{
  int i; 
  for(i = 0; i < n; i++) {
    int tmp = index[i];
    if ( counter[tmp] < base )
      counter[tmp] = base;
    else
      counter[tmp]++;
  }
}

void 
check_vlaid_perms(int *perms, int needed_pivots, int nb_pivots)
{
  int i;
  
  for(i = 0; i < needed_pivots; i++) {
    if (perms[i] == TP_UNUSED_PIVOT) {
      continue;
    } else if (perms[i] < 0) {
      TP_warning(__FUNCTION__, __FILE__, __LINE__, "one of the pivots is negative ");
      continue;
    } else if (perms[i] >= needed_pivots) { 
      TP_warning(__FUNCTION__, __FILE__, __LINE__, "one of the pivots is larger then n ");
      continue;
    }      
  }
}

void 
check_perms_and_invr_perms(int *perms, int *invr_perms, int nb_pivots, char *name)
{
  int i;
  char mess[2048];
  
  for(i = 0; i < nb_pivots; i++) 
    if (invr_perms[perms[i]] != i ) {
      snprintf(mess, 2048, "%s_perm[%d] = %d but invr_%s_perm[%d] = %d",
	       name, i, perms[i], name, perms[i], invr_perms[perms[i]]);
      TP_warning(__FUNCTION__, __FILE__, __LINE__, mess);
    }
}

void
int_array_memset(int *array, int val, int n)
{
  int i;
  for (i = 0; i < n; i++)
    array[i] = val;
}

int *
create_randomize(int n)
{
  int *vect = malloc((size_t) n * sizeof(*vect));
  int i;

  for (i = 0; i < n; i++) 
    vect[i] = i;

  for (i = 0; i < n; i++) {
    int swap = rand() % n, tmp;
    tmp  = vect[swap];
    vect[swap] = vect[i];
    vect[i] = tmp;
  }

  return vect;
}

TP_overlaps
check_overalping_regions(long r1_start, long r1_end,
                         long r2_start, long r2_end)
{
  if(r1_start <= r2_start && r1_end > r2_start) {
    if(r1_end >= r2_end) 
      return TP_overlap_total;
    else 
      return TP_overlap_begin;
  }
  
  if(r1_start >= r2_start && r1_end <= r2_end)
    return TP_overlap_total;
  
  if(r1_start < r2_end && r1_end >= r2_end) {
    if(r1_start <= r2_start) 
      return TP_overlap_total;
    else 
      return TP_overlap_end;
  }

  return TP_overlap_none;
}


void
TP_check_counters(int counter, int *array, int *used,
		  int base, int size, int n)
{
  int i, nb_elem = n * 2 * size;
  char mess[2048];
  
  for( i = 0; i < size; i++) 
    if (used[i])  {
      snprintf(mess, 2048, "in counter %d : used %d is still active", counter, i);
      TP_warning(__FUNCTION__, __FILE__, __LINE__, mess);
    }
  
  for( i = 0; i < nb_elem; i++) 
    if ( array[i] >= base) {
      snprintf(mess, 2048, "in counter %d : array[%d] = (%d) larger then base (%d)",
	       counter, i, array[i], base);
      TP_warning(__FUNCTION__, __FILE__, __LINE__, mess);
    }  
}

int 
TP_rand_int(int *seed, int size)
{
  *seed = *seed * 1103515245 + 12345;
  return (*seed /65536) % size;
}

double
TP_rand_double(int *seed) 
{
  return fabs((double) TP_rand_int(seed, TP_RAND_MAX) / (double) TP_RAND_MAX);
}
