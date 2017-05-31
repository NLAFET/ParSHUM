#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <libgen.h>
#include <string.h>

#include "TP_auxiliary.h"

void 
TP_fatal_error(const char *func, char *filename, const int line, const char *msg)
{
  printf("TP_FATAL ERROR in %s(%s:%d) \n%s\n", func, basename(filename), line, msg);
  abort();
}

void 
TP_warning(const char *func, char *filename, const int line, const char *msg)
{
  printf("TP_warning in %s(%s:%d) \n%s\n", func, basename(filename), line, msg);
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
update_counter(int *counter, int *index, int n)
{
  int i; 
  for(i = 0; i < n; i++)
    counter[index[i]]++;
}

void 
check_vlaid_perms(int *perms, int needed_pivots, int nb_pivots)
{
  int i, _nb_pivots = 0;
  int tmp[needed_pivots];
  
  memset((void *) tmp, 0, (size_t) needed_pivots*sizeof(*tmp));
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
    tmp[perms[i]]++;
    _nb_pivots++;
  }
  
  if (nb_pivots != _nb_pivots) 
      TP_warning(__FUNCTION__, __FILE__, __LINE__, "found pivots is not equal to nb_pivots");
     
  for(i = 0; i < needed_pivots; i++) 
    if ( tmp[i] > 1) {
      TP_warning(__FUNCTION__, __FILE__, __LINE__, "found a double pivot");
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
