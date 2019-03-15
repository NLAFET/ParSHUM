#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <libgen.h>
#include <string.h>

#include "ParSHUM_auxiliary.h"

#define ParSHUM_RAND_MAX 32767

void 
ParSHUM_fatal_error(const char *func, char *filename, const int line, const char *msg)
{
  printf("ParSHUM_FATAL ERROR in %s(%s:%d) \n%s\n", func, basename(filename), line, msg);
  GDB_BREAK;
  abort();
}

void 
ParSHUM_warning(const char *func, char *filename, const int line, const char *msg)
{
  printf("ParSHUM_warning in %s(%s:%d) \n%s\n", func, basename(filename), line, msg);
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
print_double_array(double *array, int n, char *mess)
{
  int i;
  printf("%s\n", mess);
  for(i = 0; i < n; i++)
    printf("%f   ", array[i]);
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
check_vlaid_perms(int *perms, int *invr_perms, int n, int nb_pivots, char *name)
{
  int i, j, found_pivots = 0;
  char mess[2048];

  for(i = 0; i < nb_pivots; i++) {
    int pivot = perms[i];
    if (pivot < 0) {
      snprintf(mess, 2048, "%s_perm[%d] = %d", name, i, pivot);
      ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
      continue;
    } else if (pivot >= n) { 
      snprintf(mess, 2048, "%s_perm[%d] = %d", name, i, pivot);
      ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
      continue;
    } else { 
      for( j = i+1; j < nb_pivots; j++) 
	if (perms[j] == pivot) 
	  snprintf(mess, 2048, "%s_perm[%d] = %d and %s_perm[%d] = %d",
		   name, i, pivot, name, j, perms[j]);
    }

  }

  for(i = 0; i < n; i++) {
    int pivot = invr_perms[i];
    if (pivot == ParSHUM_UNUSED_PIVOT) {
      continue;
    } else {
      found_pivots++;
      if (pivot < 0) {
	snprintf(mess, 2048, "invr_%s_perm[%d] = %d", name, i, pivot);
	ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
	continue;
      } else if (pivot >= n) { 
	snprintf(mess, 2048, "invr_%s_perm[%d] = %d", name, i, pivot);
	ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
	continue;
      }
      for( j = i+1; j < n; j++) 
	if (invr_perms[j] == pivot) 
	  snprintf(mess, 2048, "invr_%s_perm[%d] = %d and invr_%s_perm[%d] = %d",
		   name, i, pivot, name, j, invr_perms[j]);
    }
  }

  if (found_pivots != nb_pivots) 
    ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, "did not found all the pivots in the invr_pivots");

  for(i = 0; i < nb_pivots; i++) 
    if (invr_perms[perms[i]] != i ) {
      snprintf(mess, 2048, "%s_perm[%d] = %d but invr_%s_perm[%d] = %d",
	       name, i, perms[i], name, perms[i], invr_perms[perms[i]]);
      ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
    }
}

void 
check_vect_doubles(int *vect, int n, char *name)
{
  int i, j;
  char mess[2048];

  for(i = 0; i < n; i++) {
    int tmp = vect[i];
    for(j = i+1; j < n; j++) 
      if (vect[j] == tmp) {
	snprintf(mess, 2048, "in %s the entrie %d is doubled (vect[%d] = %d and vect[%d] = %d)", 
		 name, tmp, i, tmp, j, vect[j]);
	ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
      }
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

ParSHUM_overlaps
check_overalping_regions(long r1_start, long r1_end,
                         long r2_start, long r2_end)
{
  if(r1_start <= r2_start && r1_end > r2_start) {
    if(r1_end >= r2_end) 
      return ParSHUM_overlap_total;
    else 
      return ParSHUM_overlap_begin;
  }
  
  if(r1_start >= r2_start && r1_end <= r2_end)
    return ParSHUM_overlap_total;
  
  if(r1_start < r2_end && r1_end >= r2_end) {
    if(r1_start <= r2_start) 
      return ParSHUM_overlap_total;
    else 
      return ParSHUM_overlap_end;
  }

  return ParSHUM_overlap_none;
}


void
ParSHUM_check_counters(int counter, int *array, int *used,
		  int base, int size, int n)
{
  int i, nb_elem = n * 2 * size;
  char mess[2048];
  
  for( i = 0; i < size; i++) 
    if (used[i])  {
      snprintf(mess, 2048, "in counter %d : used %d is still active", counter, i);
      ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
    }
  
  for( i = 0; i < nb_elem; i++) 
    if ( array[i] >= base) {
      snprintf(mess, 2048, "in counter %d : array[%d] = (%d) larger then base (%d)",
	       counter, i, array[i], base);
      ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
    }  
}

int 
ParSHUM_rand_int(int *seed, int size)
{
  *seed = *seed * 1103515245 + 12345;
  return (*seed /65536) % size;
}

double
ParSHUM_rand_double(int *seed) 
{
  return fabs((double) ParSHUM_rand_int(seed, ParSHUM_RAND_MAX) / (double) ParSHUM_RAND_MAX);
}
