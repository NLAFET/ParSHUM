#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <omp.h>
#include "TP_auxiliary.h"
#include "TP_schur_matrix.h"

struct _free_space {
  long nb_elem;
  long offset;
  free_space next;
  free_space previous;
};

struct _schur_mem {
  double **val;
  int    **row;
  int    **col;
  
  free_space *unused_CSC;
  free_space *unused_CSR;
  omp_lock_t CSC_lock;
  omp_lock_t CSR_lock;

  int nb_CSC;
  int nb_CSR;
  long init_size;
};

TP_schur_matrix
TP_schur_matrix_create()
{
  TP_schur_matrix self = calloc((size_t) 1, sizeof(*self));

  return self;
}

schur_mem
TP_schur_mem_create(long nnz)
{
  schur_mem self = calloc(1, sizeof(*self));

  self->init_size   = nnz;
  self->nb_CSC      = 1;
  self->nb_CSR      = 1;

  self->val  = malloc((size_t) sizeof(*self->val));
  self->row  = malloc((size_t) sizeof(*self->row));
  self->col  = malloc((size_t) sizeof(*self->col));
  *self->val = malloc((size_t) nnz * sizeof(**self->val));
  *self->row = malloc((size_t) nnz * sizeof(**self->row));
  *self->col = malloc((size_t) nnz * sizeof(**self->col));
 
  self->unused_CSC  = malloc(sizeof(*self->unused_CSC));
  *self->unused_CSC = calloc(1, sizeof(**self->unused_CSC));
  self->unused_CSC[0]->nb_elem = nnz;
  
  self->unused_CSR  = malloc(sizeof(*self->unused_CSR));
  *self->unused_CSR = calloc(1, sizeof(**self->unused_CSR));
  self->unused_CSR[0]->nb_elem = nnz;

  omp_init_lock(&self->CSC_lock);
  omp_init_lock(&self->CSR_lock);
  
  return self;
}

void
TP_schur_mem_destroy(schur_mem self)
{
  int i, nb_CSC = self->nb_CSC, nb_CSR = self->nb_CSR;

  for( i = 0; i < nb_CSC; i++)
    {
      free_space current_unused = self->unused_CSC[i];
      while( current_unused ) { 
	free_space tmp = current_unused->next;
	free(current_unused);
	current_unused = tmp;
      }
      free(self->val[i]);
      free(self->row[i]);
    }
  free(self->val);
  free(self->row);
  free(self->unused_CSC);

  for( i = 0; i < nb_CSR; i++)
    {
      free_space current_unused = self->unused_CSR[i];
      while( current_unused ) { 
	free_space tmp = current_unused->next;
	free(current_unused);
	current_unused = tmp;
      }
      free(self->col[i]);
    }
  free(self->col);
  free(self->unused_CSR);

  omp_destroy_lock(&self->CSC_lock);
  omp_destroy_lock(&self->CSR_lock);
  
  free(self);
}

/* TODO: nutex  and concurent execution */
void
TP_schur_mem_realloc_CSC(schur_mem self)
{
  int i;
  long size = self->init_size;
  for( i = 1; i < self->nb_CSC; i++)
    size *= 2;

  self->nb_CSC++;
  
  self->unused_CSC  = realloc((void *) self->unused_CSC, (size_t) self->nb_CSC * sizeof(*self->unused_CSC));
  self->unused_CSC[self->nb_CSC - 1]  = calloc(1, sizeof(**self->unused_CSC));
  self->unused_CSC[self->nb_CSC - 1]->nb_elem = size;

  self->val = realloc((void *) self->val, (size_t) self->nb_CSC * sizeof(*self->val));
  self->row = realloc((void *) self->row, (size_t) self->nb_CSC * sizeof(*self->row));
  self->val[self->nb_CSC - 1] = malloc((size_t) size * sizeof(**self->val));
  self->row[self->nb_CSC - 1] = malloc((size_t) size * sizeof(**self->row));
}


/* TODO: nutex  and concurent execution */
void
TP_schur_mem_realloc_CSR(schur_mem self)
{
  int i;
  long size = self->init_size;
  for( i = 1; i < self->nb_CSR; i++)
    size *= 2;

  self->nb_CSR++;
  
  self->unused_CSR  = realloc((void *) self->unused_CSR, (size_t) self->nb_CSR * sizeof(*self->unused_CSR));
  self->unused_CSR[self->nb_CSR - 1]  = calloc(1, sizeof(**self->unused_CSR));
  self->unused_CSR[self->nb_CSR - 1]->nb_elem = size;

  self->col = realloc((void *) self->col, (size_t) self->nb_CSR * sizeof(*self->col));
  self->col[self->nb_CSR - 1] = malloc((size_t) size * sizeof(**self->col));
}

void
TP_schur_mem_CSC_alloc(schur_mem self, CSC_struct *CSC, long size)
{
  omp_set_lock(&self->CSC_lock);
  int i, nb = self->nb_CSC; 

  /* first try to find memory in the current space  */
  for( i = 0; i < nb; i++)
    {
      free_space current_unused = self->unused_CSC[i];
      while(current_unused) {
	if (current_unused->nb_elem >=size ) {
	  if (CSC->nb_elem)  {
	    memcpy(&self->val[i][current_unused->offset],CSC->val,  CSC->nb_elem * sizeof(*CSC->val));
	    memcpy(&self->row[i][current_unused->offset],CSC->row,  CSC->nb_elem * sizeof(*CSC->row));
	  }
	  CSC->nb_free = size - CSC->nb_elem;
	  CSC->val = &self->val[i][current_unused->offset];
	  CSC->row = &self->row[i][current_unused->offset];
	  current_unused->offset  += size;
	  current_unused->nb_elem -= size;
	  omp_unset_lock(&self->CSC_lock);
	  return;
	}
	current_unused = current_unused->next;
      }
    }

  /* if there is no spaace then we reallocate */
  while (1) {
    TP_schur_mem_realloc_CSC(self);
    free_space current_unused = self->unused_CSC[i];
    if (current_unused->nb_elem >=size ) {
      if (CSC->nb_elem)  {
	memcpy(&self->val[i][current_unused->offset], CSC->val,  CSC->nb_elem * sizeof(*CSC->val));
	memcpy(&self->row[i][current_unused->offset], CSC->row,  CSC->nb_elem * sizeof(*CSC->row));
      }
      CSC->nb_free = size - CSC->nb_elem;
      CSC->val = &self->val[i][current_unused->offset];
      CSC->row = &self->row[i][current_unused->offset];
      current_unused->offset  += size;
      current_unused->nb_elem -= size;
      omp_unset_lock(&self->CSC_lock);
      return;
    }
    i++;
  }
}
#define TP_schur_mem_CSC_realloc(self, CSC) TP_schur_mem_CSC_alloc(self, CSC, 2*CSC->nb_elem)


void
TP_schur_mem_CSR_alloc(schur_mem self, CSR_struct *CSR, long size)
{
  omp_set_lock(&self->CSR_lock);
  int i, nb = self->nb_CSR; 

  /* first try to find memory in the current space  */
  for( i = 0; i < nb; i++)
    {
      free_space current_unused = self->unused_CSR[i];
      while(current_unused) {
	if (current_unused->nb_elem >=size ) {
	  if (CSR->nb_elem)  
	    memcpy(&self->col[i][current_unused->offset], CSR->col, CSR->nb_elem * sizeof(*CSR->col));
	  CSR->nb_free = size - CSR->nb_elem;
	  CSR->col = &self->col[i][current_unused->offset];
	  current_unused->offset  += size;
	  current_unused->nb_elem -= size;
	  omp_unset_lock(&self->CSR_lock);
	  return;
	}
	current_unused = current_unused->next;
      }
    }

  /* if there is no spaace then we reallocate */
  while (1) {
    TP_schur_mem_realloc_CSR(self);
    free_space current_unused = self->unused_CSR[self->nb_CSR - 1];
    if (current_unused->nb_elem >=size ) {
      if (CSR->nb_elem)  
	memcpy(&self->col[self->nb_CSR - 1][current_unused->offset], CSR->col, CSR->nb_elem * sizeof(*CSR->col));
      CSR->nb_free = size - CSR->nb_elem;
      CSR->col = &self->col[self->nb_CSR - 1][current_unused->offset];
      current_unused->offset  += size;
      current_unused->nb_elem -= size;
      omp_unset_lock(&self->CSR_lock);
      return;
    }
  }
}
#define TP_schur_mem_CSR_realloc(self, CSR) TP_schur_mem_CSR_alloc(self, CSR, 2*CSR->nb_elem)

void
TP_schur_matrix_allocate(TP_schur_matrix self, int n, int m, long nnz, int debug, 
			 int nb_threads, double extra_space, double extra_space_inbetween)
{
  int i;

  self->nb_threads = nb_threads;
  self->n = n;
  self->m = m;

  self->CSC = calloc( (size_t) n, sizeof(*self->CSC));
  self->CSR = calloc( (size_t) m, sizeof(*self->CSR));

  self->extra_space = extra_space_inbetween;

  self->internal_mem = TP_schur_mem_create((long) nnz * (1 + extra_space + extra_space_inbetween));

  self->nnz   = 0;
  self->debug = debug;

  self->row_struct = malloc((size_t) nb_threads * sizeof(*self->row_struct));
  for( i = 0; i < nb_threads; i++)
    self->row_struct[i] = malloc((size_t) n * sizeof(**self->row_struct));
  
  self->row_locks = malloc(m * sizeof(*self->row_locks));
  for( i = 0; i < m; i++)
    omp_init_lock(&self->row_locks[i]);

  self->col_locks = malloc(m * sizeof(*self->col_locks));
  for( i = 0; i < m; i++)
    omp_init_lock(&self->col_locks[i]);
  
  if ( self->debug & (TP_DEBUG_GOSSIP_GIRL | TP_DEBUG_GARBAGE_COLLECTOR))
    TP_print_GB(self, "GB: after allocating");
}

void
TP_schur_get_singletons(TP_schur_matrix self, int done_pivots, 
			int *nb_col_singletons, int *nb_row_singletons,
			int *col_perm, int *row_perm, 
			int *invr_col_perm, int *invr_row_perm)
{
  int n = self->n, m = self->m, i;
  int _nb_col_singletons = 0, _nb_row_singletons = 0;

  for(i = 0; i < m; i++) 
    if ( self->CSR[i].nb_elem == 1 &&
	 self->CSC[self->CSR[i].col[0]].nb_elem > 0)
      {
	int col = self->CSR[i].col[0];
	row_perm[done_pivots + _nb_row_singletons] = i;
	invr_row_perm[i]  = _nb_row_singletons + done_pivots;
	col_perm[done_pivots +_nb_row_singletons] = col;
	invr_col_perm[col] = _nb_row_singletons + done_pivots;
	_nb_row_singletons++;
      }
  done_pivots += _nb_row_singletons;

  for(i = 0; i < n; i++)
    if ( self->CSC[i].nb_elem == 1 && self->CSR[self->CSC[i].row[0]].nb_elem > 1)
      {
	int row = self->CSC[i].row[0];
	col_perm[done_pivots + _nb_col_singletons] = i;
	invr_col_perm[i]  = _nb_col_singletons + done_pivots;
	row_perm[done_pivots + _nb_col_singletons] = row;
	invr_row_perm[row] = _nb_col_singletons + done_pivots;
	_nb_col_singletons++;
      }

  *nb_col_singletons = _nb_col_singletons;
  *nb_row_singletons = _nb_row_singletons;
}

void 
TP_schur_matrix_init_ptr(TP_schur_matrix self, long *col_ptr, int *row_sizes)
{
  int i;
  int n = self->n;
  int m = self->m;
  schur_mem memory = self->internal_mem;

  // handeling the CSC part
  for(i = 0; i < n; i++)
    TP_schur_mem_CSC_alloc(memory, &self->CSC[i], (long) (1 + self->extra_space) * (col_ptr[i+1] - col_ptr[i])); 
  
  // handeling the CSR part
  for(i = 0; i < m; i++)
    TP_schur_mem_CSR_alloc(memory, &self->CSR[i], (long) (1 + self->extra_space) * row_sizes[i]);

  /* TODO: PRINT GB */
  /* if ( self->debug & (TP_DEBUG_GOSSIP_GIRL | TP_DEBUG_GARBAGE_COLLECTOR)) */
  /*   TP_print_GB(self, "GB: after initializing pointers"); */
}

void
TP_schur_matrix_copy(TP_matrix A, TP_schur_matrix self)
{
  int i, j;
  int *row_sizes;
  int n = self->n;
  
  row_sizes = TP_matrix_rows_sizes(A);
  TP_schur_matrix_init_ptr(self, A->col_ptr, row_sizes);
  free(row_sizes);
  self->nnz = A->nnz;

  for( i = 0; i < n; i++) 
    {
      long    A_col_start = A->col_ptr[i];
      long    A_col_end   = A->col_ptr[i+1];
      long    col_length  = A_col_end - A_col_start;
      double *CSC_vals    = self->CSC[i].val;
      int    *CSC_rows    = self->CSC[i].row;
      
      // handle the copy  of the column into the CSC structure
      memcpy((void *) CSC_rows,
	     (void *) &A->row[A_col_start],
	     col_length * sizeof(*A->row));

      memcpy((void *) CSC_vals,
	     (void *) &A->val[A_col_start],
	     col_length * sizeof(*A->val));

      self->CSC[i].col_max = get_max_double(CSC_vals, col_length);
      
      self->CSC[i].nb_elem += col_length; 
      self->CSC[i].nb_free -= col_length; 

      // handle the copy  of the column into the CSR structure
      for(j = A_col_start; j < A_col_end; j++)
	{       
	  int row = A->row[j];
	  CSR_struct *CSR = &self->CSR[row];
	  int *CSR_cols = CSR->col;
	  
	  CSR_cols[CSR->nb_elem++] = i;
	  CSR->nb_free--;
	}
    }
}

void
delete_entry_from_CSR(TP_schur_matrix self, int col, int row)
{
  CSR_struct *CSR;
  int i, nb_elem, found, *cols;

  CSR     = &self->CSR[row];
  cols    = CSR->col;
  nb_elem = CSR->nb_elem;
  found = 0;

  if(nb_elem < 1) {
    TP_warning(__FUNCTION__, __FILE__, __LINE__,"tring to delete an entry in CSR with zero elems");
    if(cols[0] == col) 
    TP_warning(__FUNCTION__, __FILE__, __LINE__,"even better, the entry is there!");
  }

  /* TODO: this is stupid, we should somehow mark the unused enteris, not move all the memory */
  /* TODO: instead od moving the rest, just move the last one to the free position */
  for(i = 0; i < nb_elem; i++) 
    if (cols[i] == col) {
	found = 1;
	break;
    }

  if ( !found )
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"tring to delete an non existing entry in CSR");

  CSR->nb_elem--;
  CSR->nb_free++;
  cols[i] = cols[CSR->nb_elem];
}

double
delete_entry_from_CSC(TP_schur_matrix self, int col, int row)
{
  CSC_struct *CSC;
  int i, nb_elem, found, *rows;
  double *vals, return_val;

  CSC     = &self->CSC[col];
  rows    = CSC->row;
  vals    = CSC->val;
  nb_elem = CSC->nb_elem;
  found = 0;

  for(i = 0; i < nb_elem; i++) 
    if (rows[i] == row) {
      found = 1;
      return_val = vals[i];
      break;
    }

  if ( !found ) 
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"tring to delete an non existing entry in CSC");
  
  CSC->nb_elem--;
  CSC->nb_free++;
  rows[i] = rows[CSC->nb_elem];
  vals[i] = vals[CSC->nb_elem];

  return return_val;
}

inline double
delete_entry(TP_schur_matrix self, int col, int row)
{
  delete_entry_from_CSR(self, col, row);
  return delete_entry_from_CSC(self, col, row);
}

void
TP_schur_matrix_update_LD(TP_schur_matrix self, TP_matrix L, TP_matrix D,
			  int *row_perm, int *col_perm, int nb_pivots)
{
  int pivot, nb_threads = self->nb_threads;
  int nb_steps = (nb_pivots +  nb_threads - 1 ) / nb_threads, step; 

  int L_input_size = L->n;
  int D_input_size = D->n;
  
  if ( D->n + nb_pivots > D->allocated ) 
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "not enought memory in D matrix. this should never happen, so something went wrong");

  for(pivot = 0; pivot < nb_pivots; pivot++)  {
    int nb_elem = self->CSC[col_perm[pivot]].nb_elem - 1;
    int new_col_ptr = L->col_ptr[L->n] + nb_elem;
    if ( new_col_ptr > L->allocated)
      TP_matrix_realloc(L);

    L->n++;
    L->col_ptr[L->n] = new_col_ptr;
    L->nnz += nb_elem;
  }
  D->n += nb_pivots;

#pragma omp parallel num_threads(nb_threads) private(step)
  {
    long S_nnz = 0;
    for(step = 0; step < nb_steps; step++) 
      {
	int me =  omp_get_thread_num();
	int current_pivot = step * nb_threads + me;
	if ( current_pivot < nb_pivots)  {
	  CSC_struct *CSC;
	  int i, nb_elem, L_current_col, row, col;
	  int *rows, *L_rows;
	  double *vals, *L_vals, pivot_val;
	  
	  col     = col_perm[current_pivot];
	  row     = row_perm[current_pivot];
	  CSC     = &self->CSC[col];
	  nb_elem = CSC->nb_elem;
	  vals    = CSC->val;
	  rows    = CSC->row;
	  L_current_col = L->col_ptr[L_input_size + current_pivot];
	  L_rows        = L->row;
	  L_vals        = L->val;
	  
	  if ( L->nnz + nb_elem > L->allocated ) {
	    TP_matrix_realloc(L);
	    L_rows        = L->row;
	    L_vals        = L->val;
	  }
	  
	  for(i = 0; i < nb_elem; i++)
	    {
	      if ( rows[i] != row) {
		L_rows[L_current_col] = rows[i];
		L_vals[L_current_col] = vals[i];
		L_current_col++;
	      } else {
		D->val[D_input_size + current_pivot] = pivot_val = vals[i];
	      }
	      omp_set_lock(&self->row_locks[rows[i]]);
	      delete_entry_from_CSR(self, col, rows[i]);
	      omp_unset_lock(&self->row_locks[rows[i]]);
	    }
	  
	  /* TODO: we could split the previopud for in two fors: one before we found the pivot, update the begining, and then do the rest */
	  /* TODO do the delete_entru_from_CSR  in a seperate loop maybe better???  try that option */
	  for( i = L->col_ptr[L_input_size + current_pivot]; i < L->col_ptr[L_input_size + current_pivot + 1]; i++)
	    L_vals[i] /= pivot_val;
	  
	  /* TODO: recycle the col's memory */
	  CSC->nb_elem = 0;
	  CSC->nb_free = 0;
	  S_nnz += nb_elem;
	}
      }
#pragma omp atomic
    self->nnz   -= S_nnz;
  }
}

void
TP_schur_matrix_update_U_singletons(TP_schur_matrix S, TP_U_matrix U, 
				    TP_matrix D, TP_matrix L, int nb_pivots,
				    int *col_perm, int *row_perm)
{
  int nb_threads = S->nb_threads, d, sthg = L->col_ptr[L->n];
  int nb_steps = ( nb_pivots + nb_threads -1 ) / nb_threads, step;
  if ( D->n + nb_pivots > D->allocated ) 
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "not enought memory in D matrix. this should never happen, so something went wrong");

  for ( d = 0; d < nb_pivots; d++) {
    L->n++;
    L->col_ptr[L->n] = sthg;
  }

#pragma omp parallel num_threads(nb_threads) private(step)
  {
    int me =  omp_get_thread_num();
    for ( step = 0; step < nb_steps; step++) 
      {
	int current_pivot = step * nb_threads + me;
	if ( current_pivot < nb_pivots)  {
	  CSC_struct *CSC;
	  CSR_struct *CSR;
	  int col, row, D_indice, row_n, i;
	  int *row_cols;

	  col = col_perm[current_pivot];
	  row = row_perm[current_pivot];

	  CSC = &S->CSC[col]; 
	  CSR = &S->CSR[row];
	  row_cols = CSR->col;

	  if (CSC->nb_elem != 1) { 
	    printf("nb_elem = %d\n", CSC->nb_elem);
	    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "The pivot is not row singelton");
	  }
	  if (CSC->row[0] != row) 
	    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "The pivot is not the same as before");
	  D_indice = __atomic_fetch_add(&D->n, 1, __ATOMIC_SEQ_CST);
	  D->val[D_indice] = CSC->val[0];

	  delete_entry_from_CSR(S, col, row);	  
	  row_n = CSR->nb_elem;

	  for ( i = 0; i < row_n; i++)
	    {
	      int col1 = row_cols[i]; 
	      U_col *u_col = &U->col[col1];
	      double val;

	      omp_set_lock(&S->col_locks[col1]);
	      val = delete_entry_from_CSC(S, col1, row);
	      omp_unset_lock(&S->col_locks[col1]);
	      
	      omp_set_lock(&u_col->lock);
	      if (u_col->nb_elem == u_col->allocated) 
		TP_U_col_realloc(u_col);
	      u_col->row[u_col->nb_elem] = row;
	      u_col->val[u_col->nb_elem] = val;
	      u_col->nb_elem++;
	      omp_unset_lock(&u_col->lock);
	    }

	  CSR->nb_elem = 0;
	  CSR->nb_free = 0;
	  CSC->nb_elem = 0;
	  CSC->nb_free = 0;
#pragma omp atomic
	  S->nnz -= row_n + 1;
	}
      }
  }
}

void
TP_schur_matrix_update_U(TP_schur_matrix S, TP_U_matrix U, TP_matrix L, 
			 int nb_pivots, int *row_perm,
			 TP_U_struct *U_struct, int U_new_n, int U_new_nnz)
{
  int pivot, i, j;
  int nb_threads = S->nb_threads;
  int indices[nb_threads+1];
  int nnz_part = U_new_nnz / nb_threads;
  
  if ( nb_pivots <= nb_threads) {
    nb_threads = nb_pivots;
    for( i = 0; i <= nb_threads; i++)
      indices[i] = i;
  } else { 
    *indices = 0;
    for ( i = 1, j = 0; i < nb_threads; i++) 
      {
	int part = 0;
	while ( part < nnz_part &&  j < nb_pivots )
	  part += S->CSR[row_perm[j++]].nb_elem;
	indices[i] = j;
      }
    indices[nb_threads] = nb_pivots;
  }
  
  for( i = 0; i < U_new_n; i++) {
    int col = U_struct[i].col;
    int nb_elem = U_struct[i].nb_elem;
    U_col *u_col = &U->col[col];
    
    u_col->cost = 0;
    while(u_col->allocated - u_col->nb_elem  < nb_elem)
      TP_U_col_realloc(u_col);
  }
  
#pragma omp parallel num_threads(nb_threads) private(pivot, i)
  {
    int me =  omp_get_thread_num();
    long S_nnz = 0;
    for( pivot = indices[me]; pivot < indices[me+1]; pivot++)
      {
	int row = row_perm[pivot];
	
	CSR_struct *CSR = &S->CSR[row];
	int row_nb_elem = CSR->nb_elem;
	int *cols = CSR->col;
	
	for( i = 0; i < row_nb_elem; i++) {
	  int current_col = cols[i];
	  U_col  *u_col = &U->col[current_col];
	  int indice = __atomic_fetch_add(&u_col->nb_elem, 1, __ATOMIC_SEQ_CST);
	  
	  u_col->row[indice] = row;
	}
	
	CSR->nb_elem = 0;
	CSR->nb_free = 0;
	S_nnz    += row_nb_elem;
      }
#pragma omp atomic
    S->nnz -= S_nnz;
  }
}

/* NAREDNO: VIDI SO KE BIDE NAREDNATA STRUKTURA I NAPRAVI CHECK SO NEA.
 args: Schur i structurata i rows struct  */

/* void  */
/* TP_check_U_update(TP_U_matrix U, int n, int base, int *rows_struct) */
/* { */
/*   for( i = 0; i < n; i++) */
/*     int row_new_elems = ( rows_struct[i] - base_value + 1; */
/*     if (expected_new[i] >= 0) { */
/*       if (expected_new[i] != U->col[i].nb_new) */
/*         if (expected_new[i] == 1 && U->col[i].nb_new ==0) */
/*           pivots_in_U++; */
/*         else */
/*           printf("for the %dth col: expected = %d  and new = %d with basevaue of %d \n",
	     i, expected_new[i], U->col[i].nb_new, base_value); */
/*     } else { */
/*       if ( U->col[i].nb_new > 0) */
/*         printf("for the %dth col: non are expected  and there is %d new\n", i, U->col[i].nb_new); */
/*     } */
/*   if ( pivots_in_U != nb_pivots)  */
/*     printf("expected %d pivots but we have found %d\n", nb_pivots, pivots_in_U); */
    
/*   free(expected_new); */
/* } */

void
TP_schur_matrix_add_to_entry(TP_schur_matrix self, int row, int col, double val)
{
  CSC_struct *CSC = &self->CSC[col];
  int i;
  int nb_elem  = CSC->nb_elem;
  int *rows    = CSC->row;
  double *vals = CSC->val;

  for( i = 0; i < nb_elem; i++)
    if ( rows[i] == row) { 
#pragma omp atomic
      vals[i] += val;
      return;
    }
  
  TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "entry not found");
}

void
TP_schur_matrix_insert_entry(TP_schur_matrix self, int row, int col, double val)
{
  CSC_struct *CSC = &self->CSC[col];
  CSR_struct *CSR = &self->CSR[row];

  if (CSC->nb_free <= 0) 
    TP_schur_mem_CSC_realloc(self->internal_mem, CSC);
  CSC->val[CSC->nb_elem  ] = val;
  CSC->row[CSC->nb_elem++] = row;
  CSC->nb_free--;
  
  omp_set_lock(&self->row_locks[row]);
  if (CSR->nb_free <= 0) 
    TP_schur_mem_CSR_realloc(self->internal_mem, CSR);
  CSR->col[CSR->nb_elem++] = col;
  CSR->nb_free--;
  omp_unset_lock(&self->row_locks[row]);
  
#pragma omp atomic 
  self->nnz++;
}

void
TP_schur_matrix_update_S(TP_schur_matrix S, TP_matrix L, TP_U_matrix U,
			 TP_U_struct *U_struct, int U_new_n, int *invr_row_perm)
{
  int pivot, nb_threads = S->nb_threads;
  int nb_steps = (U_new_n + nb_threads - 1) / nb_threads, step; 

  long *L_col_ptr = L->col_ptr;
  int  *L_rows    = L->row;
  double *L_vals  = L->val;
  

#pragma omp parallel num_threads(nb_threads) private(step) 
  {
  int me =  omp_get_thread_num();
  int n = S->n;
  int  i, j, k, l;
  long *schur_row_struct = S->row_struct[me];

  for ( i = 0; i < n; i++)
    schur_row_struct[i] = -1;
  
  for(step = 0; step < nb_steps; step++)
    {
      i = step * nb_threads + me;
      if ( i >= U_new_n) 
	continue;
      TP_U_struct U_col_struct = U_struct[i];
      int col = U_col_struct.col;
      int U_col_new = U_col_struct.nb_elem;
      int U_col_new_unchanged = U_col_new ;

      CSC_struct *CSC = &S->CSC[col];
      double *U_vals  = U->col[col].val + U->col[col].nb_elem - U_col_new;
      int    *U_rows  = U->col[col].row + U->col[col].nb_elem - U_col_new;

      int S_col_nb_elem      = CSC->nb_elem;
      int *S_rows            = CSC->row;
      double *S_vals         = CSC->val;
      
      /* constructing schur_row_struct and discovering the vals of U  */
      for ( k = 0; k < S_col_nb_elem; ) {
	int S_row = S_rows[k];
	int found = 0;

	for ( j = 0; j < U_col_new; j++)  {
	  /* If we find an entry that needs to go in U, put it in the
	     ebginig  of U and increment the pointers of U. additionaly  
	     decrease U_col_new (we got the original one in backup in U_col_new_unchanged).
	     For S, just put the last elemnt in the k^th place and redo 
	     k again. Decrease the number of elements of course. */
	  if (U_rows[j] == S_row ) { 

	    U_rows[j] = U_rows[0];

	    U_rows[0] = S_row;
	    U_vals[0] = S_vals[k];

	    S_rows[k] = S_rows[S_col_nb_elem - 1];
	    S_vals[k] = S_vals[S_col_nb_elem - 1];

	    S_col_nb_elem--;
	    U_col_new--;
	    found = 1;
	    U_rows++;
	    U_vals++;

	    break;
	  }
	}

 	if ( !found )
	  schur_row_struct[S_row] = (long) col*n + k++;
      }

      if (U_col_new) 
	TP_warning(__FUNCTION__, __FILE__, __LINE__,"All vals for U were not found");
      
      if (S_col_nb_elem + U_col_new_unchanged != CSC->nb_elem) 
	TP_warning(__FUNCTION__, __FILE__, __LINE__,"Somethign went wrong");

      /* UPDATE ALL THE POINTERS AND STRUCTURES  */
      U_col_new = U_col_new_unchanged;
      CSC->nb_elem -= U_col_new;
      CSC->nb_free += U_col_new;
      U_vals  = U->col[col].val + U->col[col].nb_elem - U_col_new;
      U_rows  = U->col[col].row + U->col[col].nb_elem - U_col_new;
	
      for ( l = 0; l < U_col_new; l++) {
        int L_col  = invr_row_perm[U_rows[l]];
        double U_val = U_vals[l];

        for ( k = L_col_ptr[L_col]; k < L_col_ptr[L_col + 1]; k++) {
          int row = L_rows[k];
          double val = -U_val * L_vals[k];
          if (schur_row_struct[row] >= (long) col*n &&
              schur_row_struct[row] < (long) (col+1)*n)  {
            S_vals[schur_row_struct[row] % n] += val;
          } else {
	    /* MAKE THIS FNCTION TELL U IF STHG HAS BEEN UPDATED */
            TP_schur_matrix_insert_entry(S, row, col, val);
	    S_vals         = CSC->val;
	    S_rows         = CSC->row;
	    /* maybe do not use S_col_nb_elem, cause is just updates it. 
	       it is maybe better just to use CSC->nb_elem*/
            schur_row_struct[row] = (long) col*n + CSC->nb_elem - 1;
          }
        }
      }
      CSC->col_max = get_max_double(CSC->val, CSC->nb_elem);
    }
  }
}


TP_dense_matrix 
TP_schur_matrix_convert(TP_schur_matrix S, int done_pivots)
{
  TP_dense_matrix self; 
  int col, k, i;
  int n = S->n;
  int m = S->m;
  int n_schur = n - done_pivots;
  int m_schur = m - done_pivots;
  int invr_row[m];

  self = TP_dense_matrix_create(n_schur, m_schur);
  
  for(i = 0, k=0; i < m; i++)
    {
      CSR_struct *CSR = &S->CSR[i];
      if ( !CSR->nb_elem )
	continue;
      invr_row[i] = k;
      self->original_rows[k++] = i;
    }
  
  for(col = 0, k=0;  col < n; col++)
    {
      CSC_struct *CSC = &S->CSC[col];
      int nb_elem = CSC->nb_elem;
      if ( !nb_elem ) 
	continue;
      self->original_cols[k] =  col;
      double *CSC_vals = CSC->val;
      int    *CSC_rows = CSC->row;
      for( i=0; i < nb_elem; i++) 
	{ 
	  int schur_row = invr_row[CSC_rows[i]];
	  self->val[k*m_schur + schur_row] =  CSC_vals[i];
	}
      k++;
    }

  return self;
}


void
TP_schur_matrix_print(TP_schur_matrix self, char *mess)
{
  int n = self->n;
  int m = self->m;
  
  printf("%s\n", mess);
  printf("PRINTING THE CSC PART\n");
  for(int i = 0; i < n; i++)
    {
      CSC_struct *CSC = &self->CSC[i];
      int *rows    = CSC->row;
      double *vals = CSC->val;
      int nb_elem = CSC->nb_elem;
      printf("================%d======================\n", i);
      printf("Colum's max is %f\n", CSC->col_max);
      for(int j = 0; j < nb_elem; j++)
	printf("%d:(%e)  ", rows[j], vals[j]);
      printf("\n");
    }
  
  printf("\n\nPRINTING THE CSR PART\n");
  for(int i = 0; i < m; i++)
    {
      CSR_struct *CSR = &self->CSR[i];
      int *cols = CSR->col;
       int nb_elem = CSR->nb_elem;
       printf("================%d======================\n", i);
       for(int j = 0; j < nb_elem; j++)
	 printf("(%d)  ", cols[j]);
       printf("\n");
    }
  printf("\n");
}

void
TP_schur_matrix_destroy(TP_schur_matrix self)
{
  int i;
  
  free(self->CSC);
  free(self->CSR);

  TP_schur_mem_destroy(self->internal_mem);

  for( i = 0; i < self->nb_threads; i++)
    free(self->row_struct[i]);
  free(self->row_struct);
  
  for( i = 0; i < self->m; i++)
    omp_destroy_lock(&self->row_locks[i]);
  free(self->row_locks);
  for( i = 0; i < self->n; i++)
    omp_destroy_lock(&self->col_locks[i]);
  free(self->col_locks);

  free(self);
}


/* ********************************************************************************************* */
/* ********************************************************************************************* */
/* DEBBUG */
/* DEBBUG */
/* ********************************************************************************************* */
/* ********************************************************************************************* */
/* ********************************************************************************************* */
void
TP_schur_check_doubles(TP_schur_matrix self)
{
  int col, row, i, j, n = self->n, m = self->m;
  char mess[2048];

  for(col = 0; col < n; col++)
    {
      CSC_struct *CSC = &self->CSC[col];
      int *rows = CSC->row;
      int nb_elem = CSC->nb_elem;

      for (i = 0; i < nb_elem; i++) {
	row = rows[i];
	for(j = i + 1; j <  nb_elem; j++) 
	  if (rows[j] == row) {
	    snprintf(mess, 2048, "in column %d, row %d is ducplicated on positions %d and %d",
		     col, row, i, j);
	    TP_warning(__FUNCTION__, __FILE__, __LINE__, mess);
	  }
      }
    }
  
  for(row = 0; row < m; row++)
    {
      CSR_struct *CSR = &self->CSR[row];
      int *cols   = CSR->col;
      int nb_elem = CSR->nb_elem;

      for (i = 0; i < nb_elem; i++) {
	int col = cols[i];
	for(j = i + 1; j <  nb_elem; j++) 
	  if (cols[j] == col) {
	    snprintf(mess, 2048, "in row %d, col %d is ducplicated on positions %d and %d",
		     row, col, i, j);
	    TP_warning(__FUNCTION__, __FILE__, __LINE__, mess);
	  }
      }
    }

}

void
TP_schur_matrix_check_pivots(TP_schur_matrix self,
			     int *row_perms, int *col_perms,
			     int *invr_row_perms, int *invr_col_perms,
			     int nb_pivots)
{
  int  i, j, k, n = self->n, m = self->m;
  int needed_pivots = self->n < self->m ? self->n : self->m;
  char mess[2048];

  check_vlaid_perms(row_perms,      needed_pivots, nb_pivots);
  check_vlaid_perms(col_perms,      needed_pivots, nb_pivots);
  check_vlaid_perms(invr_row_perms, needed_pivots, nb_pivots);
  check_vlaid_perms(invr_row_perms, needed_pivots, nb_pivots);

  for( i = 0; i < nb_pivots; i++)
    {
      int row = row_perms[i], col = col_perms[i];
      
      if ( self->CSC[col].nb_elem ) {
	snprintf(mess, 2048, "column %d is a pivot, but not empty in S with nb_elem %d and nb_free %d",
		 col, self->CSC[col].nb_elem, self->CSC[col].nb_free);
	TP_warning(__FUNCTION__, __FILE__, __LINE__, mess);
      }
      
      if ( self->CSR[row].nb_elem ) {
	snprintf(mess, 2048, "row %d is a pivot, but not empty in S with nb_elem %d and nb_free %d",
		 row, self->CSR[row].nb_elem, self->CSR[row].nb_free);
	TP_warning(__FUNCTION__, __FILE__, __LINE__, mess);
      }
    }

  for ( i = 0; i < n; i++) 
    {
      CSC_struct *CSC = &self->CSC[i];
      int nb_elem = CSC->nb_elem;
      int *rows = CSC->row;
      for ( j = 0; j < nb_elem; j++) 
	{
	  int row = rows[j];
	  for(k = 0; k < nb_pivots; k++)
	    if ( row_perms[k] == row ) {
	      snprintf(mess, 2048, "in col %d, %d is present, but %d is a row pivot\n", i, row, row);
	      TP_warning(__FUNCTION__, __FILE__, __LINE__, mess);
	    }
	}
    }

  for ( i = 0; i < m; i++) 
    {
      CSR_struct *CSR = &self->CSR[i];
      int nb_elem = CSR->nb_elem;
      int *cols = CSR->col;
      for ( j = 0; j < nb_elem; j++) 
	{
	  int col = cols[j];
	  for(k = 0; k < nb_pivots; k++)
	    if ( col_perms[k] == col ) {
	      snprintf(mess, 2048, "in row %d, %d is present, but %d is a row pivot\n", i, col, col);
	      TP_warning(__FUNCTION__, __FILE__, __LINE__, mess);
	    }
	}
    }
}

void
TP_schur_matrix_memory_check(TP_schur_matrix self)
{
  int i, n = self->n, m = self->m;
  char mess[2048];
  
  for ( i = 0; i < n; i++)
    if (self->CSC[i].nb_free < 0 )  {
      snprintf(mess, 2048, "error on the column %d  with nb_free %d\n", i, self->CSC[i].nb_free);
      TP_warning(__FUNCTION__, __FILE__, __LINE__, mess);
    }
  
  for ( i = 0; i < m; i++)
    if (self->CSR[i].nb_free < 0 ) {
      snprintf(mess, 2048, "error on the row %d with nb_free %d\n", i, self->CSR[i].nb_free);
      TP_warning(__FUNCTION__, __FILE__, __LINE__, mess);
    }

  /* TODO: rewrtie this */
  /* for(i = 0; i < n; i++) */
  /*   { */
  /*     CSC_struct *CSC = &self->CSC[i]; */
  /*     free_space unused_CSC  = self->unused_CSC; */
  /*     long CSC_begin = CSC->offset, CSC_end = CSC->offset + CSC->nb_elem + CSC->nb_free; */
  /*     int current_unused = 0, j; */
  /*     if (CSC_begin == CSC_end) */
  /* 	continue; */
      
  /*     while(unused_CSC) */
  /* 	{ */
  /* 	  long free_begin = unused_CSC->offset, free_end = unused_CSC->offset + unused_CSC->nb_elem; */
  /* 	  int print = 0; */
	  
  /* 	  TP_overlaps overlaped = check_overalping_regions(free_begin, free_end, CSC_begin, CSC_end); */
  /* 	  switch (overlaped) { */
  /* 	  case (TP_overlap_none) : */
  /* 	    break; */
  /* 	  case (TP_overlap_begin) : */
  /* 	    snprintf(mess, 2048, "The %d^th free space and the %d^th column are overlaping in the begining of the col (col start %ld end %ld, free starts on %ld ends on %ld).", */
  /* 		     current_unused, i, CSC_begin, CSC_end, free_begin, free_end); */
  /* 	    print = 1; */
  /* 	    break; */
  /* 	  case (TP_overlap_end) : */
  /* 	    snprintf(mess, 2048, "The %d^th free space and the %d^th column are overlaping in the end of the col (col start %ld end %ld, free starts on %ld ends on %ld).", */
  /* 		     current_unused, i, CSC_begin, CSC_end, free_begin, free_end); */
  /* 	    print = 1; */
  /* 	    break; */
  /* 	  case (TP_overlap_total) : */
  /* 	    snprintf(mess, 2048, "The %d^th free space and the %d^th column are overlapping (col starts at %ld and ends on %ld; free starts on %ld and ends on %ld).", */
  /* 		     current_unused, i, CSC_begin, CSC_end, free_begin, free_end); */
  /* 	    print = 1; */
  /* 	    break; */
  /* 	   default: */
  /* 	     break; */
  /* 	  } */
  /* 	  if (print) */
  /* 	    TP_warning(__FUNCTION__, __FILE__, __LINE__, mess); */
  /* 	  current_unused++; */
  /* 	  unused_CSC = unused_CSC->next; */
  /* 	} */
      
  /*     for(j = 0; j < n; j++) */
  /* 	{ */
  /* 	  if ( i == j) */
  /* 	    continue; */
  /* 	  CSC_struct *current_CSC = &self->CSC[j]; */
  /* 	  long current_CSC_begin = current_CSC->offset, current_CSC_end = current_CSC->offset + current_CSC->nb_elem + current_CSC->nb_free; */
  /* 	  int print = 0; */
  /* 	  if( current_CSC_begin == current_CSC_end) */
  /* 	    continue; */
	  
  /* 	  TP_overlaps overlaped = check_overalping_regions(current_CSC_begin, current_CSC_end, CSC_begin, CSC_end); */
  /* 	  switch (overlaped) { */
  /* 	  case (TP_overlap_none) : */
  /* 	    break; */
  /* 	  case (TP_overlap_begin) : */
  /* 	    snprintf(mess, 2048, "The %d^th column and the %d^th column are overlaping in the begining of the col (col start %ld end %ld; col starts %ld ends on %ld).", */
  /* 		     j, i, current_CSC_begin , current_CSC_end, CSC_begin, CSC_end); */
  /* 	    print = 1; */
  /* 	    break; */
  /* 	  case (TP_overlap_end) : */
  /* 	    snprintf(mess, 2048, "The %d^th column and the %d^th column are overlaping in the end of the col (col start %ld end %ld; col starts %ld ends on %ld).", */
  /* 		     j, i, current_CSC_end, current_CSC_end, CSC_begin, CSC_end); */
  /* 	    print = 1; */
  /* 	    break; */
  /* 	  case (TP_overlap_total) : */
  /* 	    snprintf(mess, 2048, "The %d^th column and the %d^th column are overlapping (col starts at %ld and ends on %ld; col starts %ld ends on %ld).", */
  /* 		     j, i, current_CSC_end, current_CSC_end, CSC_begin, CSC_end); */
  /* 	    print = 1; */
  /* 	    break; */
  /* 	  default: */
  /* 	    break; */
  /* 	  } */
  /* 	  if (print) */
  /* 	    TP_warning(__FUNCTION__, __FILE__, __LINE__, mess); */
  /* 	} */
  /*   } */
  
  /* for(i = 0; i < m; i++) */
  /*   { */
  /*     CSR_struct *CSR = &self->CSR[i]; */
  /*     free_space unused_CSR  = self->unused_CSR; */
  /*     long CSR_begin = CSR->offset, CSR_end = CSR->offset + CSR->nb_elem + CSR->nb_free; */
  /*     int current_unused = 0, j; */
      
  /*     if (CSR_begin == CSR_end) */
  /* 	continue; */
      
  /*     while(unused_CSR) */
  /* 	{ */
  /* 	  long free_begin = unused_CSR->offset, free_end = unused_CSR->offset + unused_CSR->nb_elem; */
  /* 	  int print = 0; */
	  
  /* 	  TP_overlaps overlaped = check_overalping_regions(free_begin, free_end, CSR_begin, CSR_end); */
  /* 	  switch (overlaped) { */
  /* 	  case (TP_overlap_none)  : */
  /* 	    break; */
  /* 	  case (TP_overlap_begin) : */
  /* 	    snprintf(mess, 2048, "The %d^th free space and the %d^th row are overlaping in the begining of the col (row start %ld end %ld, free starts on  %ld ends on %ld).", */
  /* 		     current_unused, i, CSR_begin, CSR_end, free_begin, free_end); */
  /* 	    print = 1; */
  /* 	    break; */
  /* 	  case (TP_overlap_end)   : */
  /* 	    snprintf(mess, 2048, "The %d^th free space and the %d^th row are overlaping in the end of the col (row start %ld end %ld, free starts on %ld ends on %ld).", */
  /* 		      current_unused, i, CSR_begin, CSR_end, free_begin, free_end); */
  /* 	    print = 1; */
  /* 	    break; */
  /* 	  case (TP_overlap_total) : */
  /* 	    snprintf(mess, 2048, "The %d^th free space and the %d^th column are overlapping (col starts at %ld and ends on %ld; free starts on %ld and ends on %ld).", */
  /* 		     current_unused, i, CSR_begin, CSR_end, free_begin, free_end); */
  /* 	    print = 1; */
  /* 	    break; */
  /* 	  default: */
  /* 	    break; */
  /* 	  } */
  /* 	  if (print) */
  /* 	    TP_warning(__FUNCTION__, __FILE__, __LINE__, mess); */
  /* 	  current_unused++; */
  /* 	  unused_CSR = unused_CSR->next; */
  /* 	} */
      
  /*     for(j = 0; j < m; j++) */
  /* 	{ */
  /* 	   if ( i == j) */
  /* 	     continue; */
  /* 	   CSR_struct *current_CSR = &self->CSR[j]; */
  /* 	   long current_CSR_begin = current_CSR->offset, current_CSR_end = current_CSR->offset + current_CSR->nb_elem + current_CSR->nb_free; */
  /* 	   int print = 0; */
  /* 	   if( current_CSR_begin == current_CSR_end) */
  /* 	     continue; */
	   
  /* 	   TP_overlaps overlaped = check_overalping_regions(current_CSR_begin, current_CSR_end, CSR_begin, CSR_end); */
  /* 	   switch (overlaped) { */
  /* 	   case (TP_overlap_none) : */
  /* 	     break; */
  /* 	   case (TP_overlap_begin) : */
  /* 	     snprintf(mess, 2048, "The %d^th row and the %d^th column are overlaping in the begining of the row (row start %ld end %ld; row starts %ld ends on %ld).", */
  /* 		      j, i, current_CSR_begin , current_CSR_end, CSR_begin, CSR_end); */
  /* 	     print = 1; */
  /* 	     break; */
  /* 	   case (TP_overlap_end) : */
  /* 	     snprintf(mess, 2048, "The %d^th row and the %d^th row are overlaping in the end of the row (row start %ld end %ld; row starts %ld ends on %ld).", */
  /* 		      j, i, current_CSR_end, current_CSR_end, CSR_begin, CSR_end); */
  /* 	     print = 1; */
  /* 	     break; */
  /* 	   case (TP_overlap_total) : */
  /* 	     snprintf(mess, 2048, "The %d^th column and the %d^th column are overlapping (col starts at %ld and ends on %ld; col starts %ld ends on %ld).", */
  /* 		      j, i, current_CSR_end, current_CSR_end, CSR_begin, CSR_end); */
  /* 	     print = 1; */
  /* 	     break; */
  /* 	   default: */
  /* 	     break; */
  /* 	   } */
  /* 	   if (print) */
  /* 	     TP_warning(__FUNCTION__, __FILE__, __LINE__, mess); */
  /* 	} */
  /*   } */
}

void
TP_schur_matrix_check_symetry(TP_schur_matrix self)
{
  long CSC_nnz = 0, CSR_nnz = 0;
  int i, j, n = self->n, m = self->m, row, col;
  char mess[2048];
  
  for(i = 0; i < n; i++)
    CSC_nnz += self->CSC[i].nb_elem;
  
  for(i = 0; i < m; i++)
    CSR_nnz += self->CSR[i].nb_elem;
  
  if (CSC_nnz != CSR_nnz) {
    snprintf(mess, 2048, "CSR and CSC nnz are not the same, CSC = %ld and CSR = %ld", CSC_nnz, CSR_nnz);
    TP_warning(__FUNCTION__, __FILE__, __LINE__, mess);
  }
  
  if (CSC_nnz != self->nnz) {
    snprintf(mess, 2048, "CSC and S nnz are not the same, CSC = %ld and S_nnz = %ld", CSC_nnz, self->nnz);
    TP_warning(__FUNCTION__, __FILE__, __LINE__, mess);
  }
  
  if (CSR_nnz != self->nnz) {
    snprintf(mess, 2048, "CSR and S nnz are not the same, CSR = %ld and S_nnz = %ld", CSR_nnz, self->nnz);
    TP_warning(__FUNCTION__, __FILE__, __LINE__, mess);
  }
  
  for(col = 0; col < n; col++)
    {
      CSC_struct *CSC = &self->CSC[col];
      int *rows   = CSC->row;
      int col_nb_elem = CSC->nb_elem;
      
      for(i = 0; i < col_nb_elem; i++)
	{
	  row = rows[i];
	  CSR_struct *CSR = &self->CSR[row];
	  int row_nb_elem = CSR->nb_elem;
	  int *cols = CSR->col;
	  int found = 0;
	  for(j = 0; j < row_nb_elem; j++) {
	    if(cols[j] == col) {
	      found = 1;
	      break;
	    }
	  }
	  
	  if (!found) {
	    snprintf(mess, 2048, "In CSC in col %d row %d exists, but in CSR in row %d,  col %d does not exist",
		     col, row, row, col);
	    TP_warning(__FUNCTION__, __FILE__, __LINE__, mess);
	    /* TP_schur_matrix_print(self, "before exiting"); */
	  }
	} /* for I */
    }
}

/* TODO: addapt this to the new thing */
void
TP_print_GB(TP_schur_matrix self, char *mess)
{
  /* free_space CSC = self->unused_CSC; */
  /* free_space CSR = self->unused_CSR; */
  
  /* fprintf(stdout,"%s\n", mess); */
  /* if (CSC->previous) */
  /*   TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "the first CSC free memory has a predecessor"); */
  /* if (CSR->previous) */
  /*   TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "the first CSR free memory has a predecessor"); */
  
  /* fprintf(stdout, "|================||=================|\n"); */
  /* fprintf(stdout, "|      CSC       ||      CSR        |\n"); */
  /* fprintf(stdout, "|================||=================|\n"); */
  /* fprintf(stdout, "|    address     ||    address      |\n"); */
  /* fprintf(stdout, "|    nb_elem     ||    nb_elem      |\n"); */
  /* fprintf(stdout, "|    offset      ||    offset       |\n"); */
  /* fprintf(stdout, "|================||=================|\n\n"); */
  /* fprintf(stdout, "|===================================|\n"); */
  
  /* while( CSC || CSR) { */
  /*   if (CSC) */
  /*     fprintf(stdout, "| %14p |", &CSC); */
  /*   else */
  /*     fprintf(stdout, "                  "); */
    
  /*   if (CSR) */
  /*     fprintf(stdout, "| %14p  |\n", &CSR); */
  /*   else */
  /*     fprintf(stdout, "\n"); */
    
  /*   if (CSC) */
  /*     fprintf(stdout, "|  %11ld   |", CSC->nb_elem); */
  /*   else */
  /*     fprintf(stdout, "                  "); */
    
  /*   if (CSR) */
  /*     fprintf(stdout, "|   %11ld   |\n", CSR->nb_elem); */
  /*    else */
  /*      fprintf(stdout, "\n"); */
    
  /*   if (CSC) */
  /*     fprintf(stdout, "|  %11ld   |", CSC->offset); */
  /*   else */
  /*     fprintf(stdout, "                  "); */
    
  /*   if (CSR) */
  /*     fprintf(stdout, "|   %11ld   |\n", CSR->offset); */
  /*   else */
  /*     fprintf(stdout, "\n"); */
  /*   fprintf(stdout, "|===================================|\n"); */
    
  /*   if (CSC) */
  /*     if (CSC->next) */
  /* 	if (CSC->next->previous != CSC) */
  /* 	  TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "the CSC free memory's next cell, has a predecessor different from CSC"); */

  /*   if (CSR) */
  /*     if (CSR->next) */
  /* 	if (CSR->next->previous != CSR) */
  /* 	  TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "the CSR free memory's next cell, has a predecessor different from CSR"); */
    
  /*   if (CSC) */
  /*     CSC=CSC->next; */
  /*   if (CSR) */
  /*     CSR=CSR->next; */
  /* } */
}

/* TODO: addapt this to the new thing */
void
TP_print_single_GB(free_space self, char *mess)
{
  /* fprintf(stdout,"%s\n", mess); */
  
  /* if (self->previous) */
  /*   TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "the first CSR free memory has a predecessor"); */
    
  /* fprintf(stdout, "|================|\n"); */
  /* fprintf(stdout, "|    address     |\n"); */
  /* fprintf(stdout, "|    nb_elem     |\n"); */
  /* fprintf(stdout, "|    offset      |\n"); */
  /* fprintf(stdout, "|================|\n\n"); */
  /* fprintf(stdout, "|================|\n"); */
  
  /* while(self) */
  /*   { */
  /*     fprintf(stdout, "| %14p |\n", &self); */
  /*     fprintf(stdout, "|  %11ld   |\n", self->nb_elem); */
  /*     fprintf(stdout, "|  %11ld   |\n", self->offset); */
  /*     fprintf(stdout, "|================|\n"); */
      
  /*     if (self->next) */
  /* 	if (self->next->previous != self) */
  /* 	  TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "the free memory's next cell, has a predecessor different from CSC"); */
      
  /*     self = self->next; */
  /*   } */
}


void
TP_check_current_counters(TP_schur_matrix self,
			  int *col_perm, int *row_perm, int nb_perms, 
			  int *col_count, int *row_count, int base)
{
  int pivot, i, n = self->n, m = self->m; 
  int *_col_count = calloc( n, sizeof(*_col_count));
  int *_row_count = calloc( m, sizeof(*_row_count));
  char mess[2048];

  for( pivot = 0; pivot < nb_perms; pivot++) 
    {
      int col = col_perm[pivot];
      CSC_struct *CSC = &self->CSC[col];
      int *rows              = CSC->row;
      int nb_elem            = CSC->nb_elem;
      for( i = 0; i < nb_elem; i++) 
	_row_count[rows[i]]++;

      int row = row_perm[pivot];
      CSR_struct *CSR = &self->CSR[row];
      int *cols              = CSR->col;
      nb_elem                = CSR->nb_elem;
      for( i = 0; i < nb_elem; i++) 
	_col_count[cols[i]]++;
    }

  for( pivot = 0; pivot < nb_perms; pivot++) 
    {
      int col = col_perm[pivot];
      int row = row_perm[pivot];
      
      if ( _row_count[row] != 1) {
	snprintf(mess, 2048, "calculated row_count[%d] = %d, but %d is a pivot",
		 row, _row_count[i], row);
	TP_warning(__FUNCTION__, __FILE__, __LINE__, mess);
      }

      if ( row_count[row] != base) {
	snprintf(mess, 2048, "row_count[%d] = %d, base = %d, but %d is a pivot",
		 row, row_count[i], base, row);
	TP_warning(__FUNCTION__, __FILE__, __LINE__, mess);
      }

      if ( _col_count[col] != 1) {
	snprintf(mess, 2048, "calculated col_count[%d] = %d, but %d is a pivot",
		 col, _col_count[i], row);
	TP_warning(__FUNCTION__, __FILE__, __LINE__, mess);
      }

      if ( row_count[row] != base) {
	snprintf(mess, 2048, "col_count[%d] = %d, base = %d, but %d is a pivot",
		 col, col_count[i], base, col);
	TP_warning(__FUNCTION__, __FILE__, __LINE__, mess);
      }
    }  

  for( i = 0; i < n; i++) 
    {
      if (col_count[i] >=  base) 
	if ( ( col_count[i] - base + 1 ) !=  _col_count[i]) {
	  snprintf(mess, 2048, "error on col counter %d : col_count(%d) base(%d) and calculated col_count(%d)",
		   i, col_count[i], base, _col_count[i]);
	  TP_warning(__FUNCTION__, __FILE__, __LINE__, mess);
	}
      
      if (row_count[i] >=  base) 
	if ( ( row_count[i] - base + 1 ) !=  _row_count[i]) {
	  snprintf(mess, 2048, "error on row counter %d : row_count(%d) base(%d) and calculated row_count(%d)",
		   i, row_count[i], base, _row_count[i]);
	  TP_warning(__FUNCTION__, __FILE__, __LINE__, mess);
	}
    }	 
  
  free(_col_count);
  free(_row_count);
}
