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

TP_schur_matrix
TP_schur_matrix_create()
{
  TP_schur_matrix self = calloc((size_t) 1, sizeof(*self));

  return self;
}

void
TP_schur_matrix_allocate(TP_schur_matrix self, int n, int m, long nnz, int debug, 
			 int nb_threads, double extra_space, double extra_space_inbetween)
{
  long allocating;
  int i;

  self->nb_threads = nb_threads;
  self->n = n;
  self->m = m;

  self->CSC = calloc( (size_t) n, sizeof(*self->CSC));
  self->CSR = calloc( (size_t) m, sizeof(*self->CSR));

  self->extra_space           = extra_space;
  self->extra_space_inbetween = extra_space_inbetween;
  allocating = nnz * (1 + extra_space + extra_space_inbetween);
  
  self->nnz   = 0;
  self->debug = debug;
  self->allocated_CSC = allocating;
  self->allocated_CSR = allocating;

  self->val = malloc((size_t) allocating * sizeof(*self->val));
  self->row = malloc((size_t) allocating * sizeof(*self->row));
  self->col = malloc((size_t) allocating * sizeof(*self->col));

  self->unused_CSC = calloc(1, sizeof(*self->unused_CSC));
  self->unused_CSC->nb_elem = allocating;

  self->unused_CSR = calloc(1, sizeof(*self->unused_CSR));
  self->unused_CSR->nb_elem = allocating;

  self->row_locks = malloc(m * sizeof(*self->row_locks));
  for( i = 0; i < m; i++)
    pthread_mutex_init(&self->row_locks[i], NULL);
  self->col_locks = malloc(n * sizeof(*self->col_locks));
  for( i = 0; i < n; i++)
    pthread_mutex_init(&self->col_locks[i], NULL);
  
  if ( self->debug & (TP_DEBUG_GOSSIP_GIRL | TP_DEBUG_GARBAGE_COLLECTOR))
    TP_print_GB(self, "GB: after allocating");
}

free_space
add_unused_CSC(free_space unused_space, long nb_elem, long offset)
{
  free_space self = malloc(sizeof(*self));
  assert(nb_elem > 0);

  self->nb_elem  = nb_elem;
  self->offset   = offset;
  self->next     = unused_space;
  if(unused_space)
    unused_space->previous = self;
  self->previous = NULL;

  return self;
}

free_space
add_unused_CSR(free_space unused_space, long nb_elem, long offset)
{
  free_space self = malloc(sizeof(*self));
  assert(nb_elem > 0);

  self->nb_elem  = nb_elem;
  self->offset   = offset;
  self->next     = unused_space;
  if(unused_space)
    unused_space->previous = self;
  self->previous = NULL;

  return self;
}

free_space
TP_schur_matrix_realloc_CSC(free_space unused, double **val,
			    int **row, long *allocated)
{
  long allocating = *allocated * 2;

  *val = realloc((void *) *val, allocating * sizeof(**val));
  *row = realloc((void *) *row, allocating * sizeof(**row));

  unused = add_unused_CSC (unused, *allocated, *allocated);

  *allocated  = allocating;

  return unused;
}

free_space
TP_schur_matrix_realloc_CSR(free_space unused, int **col,
			    long *allocated)
{
  long allocating = *allocated * 2;

  *col = realloc((void *) *col, allocating * sizeof(**col));

  unused = add_unused_CSR(unused, *allocated, *allocated);
  *allocated  = allocating;

  return unused;
}


free_space
CSC_find_free_memory(free_space unused_CSC,  struct CSC_struct *CSC,
		     long nb_elem, double **vals, int **rows, long *allocated_CSC)
{
  free_space tmp = unused_CSC;

  if (CSC->nb_free > 0 )
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"acquiring new CSR memory to non-full memory ");

  if (CSC->nb_elem >= nb_elem )
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"the acquired CSR memory is smaller or equal then the current one");

  while(tmp) {
    if (tmp->nb_elem  < nb_elem) 
      tmp = tmp->next;
    else 
      break;
  }

  if (!tmp) { 
    tmp = 
    /* UNUSED_CSC IS THE RETURN VLUE (THE START OF THE FREE_SPACE CHAIN) */
      unused_CSC = 
      TP_schur_matrix_realloc_CSC(unused_CSC, vals, rows, allocated_CSC);
  }

  if (CSC->nb_elem) {
    double *CSC_vals  = *vals + CSC->offset;
    double *free_vals = *vals + tmp->offset;
    int *CSC_rows     = *rows + CSC->offset;
    int *free_rows    = *rows + tmp->offset;
    
    memcpy(free_vals, CSC_vals, (size_t) CSC->nb_elem * sizeof(**vals));
    memcpy(free_rows, CSC_rows, (size_t) CSC->nb_elem * sizeof(**rows));
  }

  CSC->offset   = tmp->offset;
  CSC->nb_free  = nb_elem - CSC->nb_elem;
  tmp->offset  += nb_elem;
  tmp->nb_elem -= nb_elem;

  /* is we need to destroy the tmp */
  if (!tmp->nb_elem) 
    {
      /* tmp == unused_CSC -> need to update  unused CSC */
      if (tmp == unused_CSC) 
	{
	  if (tmp->next) {
	    /* if it is not the only one, easy */
	    unused_CSC = tmp->next;
	    unused_CSC->previous = NULL;
	  } else  {
	    /* if not, there is no more memory, so realocate , allocate more memory  */
	    unused_CSC = TP_schur_matrix_realloc_CSC(NULL, vals, rows, allocated_CSC);
	  }
	} else 
	{
	  if ( tmp->previous ) 
	    tmp->previous->next = tmp->next;
	  if ( tmp->next ) 
	    tmp->next->previous = tmp->previous;
	}      
      free(tmp);
    } 

  if (!unused_CSC) 
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"something strange is happening in S (no more memory  apperently)");

  return unused_CSC;
}

free_space
CSR_find_free_memory(free_space unused_CSR,  struct CSR_struct *CSR,
		     long nb_elem, int **cols, long *allocated_CSR)
{
  free_space tmp = unused_CSR;

  if (CSR->nb_free > 0 )
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"acquiring new CSR memory to non-full memory ");

  if (CSR->nb_elem > nb_elem )
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"the acquired CSR memory is smaller than the current one");

  while(tmp) 
    if (tmp->nb_elem  < nb_elem) 
      tmp = tmp->next;
    else 
      break;

  if (!tmp)
    tmp = 
    /* UNUSED_CSR IS THE RETURN VLUE (THE START OF THE FREE_SPACE CHAIN)  */
      unused_CSR = 
      TP_schur_matrix_realloc_CSR(unused_CSR, cols, allocated_CSR);

  if (nb_elem) {
    int *CSR_cols = *cols + CSR->offset;
    int *tmp_cols = *cols + tmp->offset;

    memcpy(tmp_cols, CSR_cols, (size_t) CSR->nb_elem * sizeof(**cols));
  }

  CSR->offset   = tmp->offset;
  CSR->nb_free  = nb_elem - CSR->nb_elem;
  tmp->offset  += nb_elem;
  tmp->nb_elem -= nb_elem;
  
  /* is we need to destroy the tmp */
  if (!tmp->nb_elem) 
    {
      /* tmp == unused_CSR -> need just an update  on unused CSR */
      if (tmp == unused_CSR) 
	{
	  if (tmp->next) {
	    /* if it is not the only one, easy */
	    unused_CSR = tmp->next;
	    unused_CSR->previous = NULL;
	  } else {
	    /* if not, there is no more memory, so realocate  */
	    unused_CSR = TP_schur_matrix_realloc_CSR(NULL, cols, allocated_CSR);
	  }
	} else 
	{
	  /* just adapt around tmp */
	  if ( tmp->previous ) 
	    tmp->previous->next = tmp->next;
	  if ( tmp->next ) 
	    tmp->next->previous = tmp->previous;
	}
      
      free(tmp);
    } 

  if (!unused_CSR) 
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"something strange is happening in S (no more memory  apperently)");

  return unused_CSR;
}



void 
TP_schur_matrix_init_ptr(TP_schur_matrix self, long *col_ptr, int *row_sizes)
{
  int i;
  int n = self->n;
  int m = self->m;
  
  // handeling the CSC part
  for(i = 0; i < n; i++)
    self->unused_CSC = CSC_find_free_memory(self->unused_CSC, &self->CSC[i], 
					    (long) (1 + self->extra_space) * (col_ptr[i+1] - col_ptr[i]), 
					    &self->val, &self->row, &self->allocated_CSC);
  
  // handeling the CSR part
  for(i = 0; i < m; i++)
    self->unused_CSR = CSR_find_free_memory(self->unused_CSR, &self->CSR[i],
					    (long) (1 + self->extra_space) * row_sizes[i],
					    &self->col, &self->allocated_CSR);

  if ( self->debug & (TP_DEBUG_GOSSIP_GIRL | TP_DEBUG_GARBAGE_COLLECTOR))
    TP_print_GB(self, "GB: after initializing pointers");
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
      double *CSC_vals    = self->val + self->CSC[i].offset;
      int    *CSC_rows    = self->row + self->CSC[i].offset;
      
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
	  struct CSR_struct *CSR = &self->CSR[row];
	  int *CSR_cols = self->col + CSR->offset;
	  
	  CSR_cols[CSR->nb_elem++] = i;
	  CSR->nb_free--;
	}
    }
}

void
delete_entry_from_CSR(TP_schur_matrix self, int col, int row)
{
  struct CSR_struct *CSR;
  int i, nb_elem, found, *cols;

  CSR     = &self->CSR[row];
  cols    = self->col + CSR->offset;
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
  for(step = 0; step < nb_steps; step++) 
    {
      int me =  omp_get_thread_num();
      int current_pivot = step * nb_threads + me;
      if ( current_pivot < nb_pivots)  {
	
	struct CSC_struct *CSC;
	int i, nb_elem, L_current_col, row, col;
	int *rows, *L_rows;
	double *vals, *L_vals, pivot_val;
	
	col     = col_perm[current_pivot];
	row     = row_perm[current_pivot];
	CSC     = &self->CSC[col];
	nb_elem = CSC->nb_elem;
	vals    = self->val + CSC->offset;
	rows    = self->row + CSC->offset;
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
	    pthread_mutex_lock(&self->row_locks[rows[i]]);
	    delete_entry_from_CSR(self, col, rows[i]);
	    pthread_mutex_unlock(&self->row_locks[rows[i]]);
	  }
	
	/* TODO: we could split the previopud for in two fors: one before we found the pivot, update the begining, and then do the rest */
	/* TODO do teh delete_entru_from_CSR  in a seperate loop maybe better???  try that option */
	for( i = L->col_ptr[L_input_size + current_pivot]; i < L->col_ptr[L_input_size + current_pivot + 1]; i++)
	  L_vals[i] /= pivot_val;
	
	/* TODO: recycle the col's memory */
	CSC->nb_elem = 0;
	CSC->nb_free = 0;
	self->nnz   -= nb_elem;
      }
    }
}



double
delete_entry_from_CSC(TP_schur_matrix self, int col, int row)
{
  struct CSC_struct *CSC;
  int i, nb_elem, found, *rows;
  double *vals, return_val;

  CSC     = &self->CSC[col];
  rows    = self->row + CSC->offset;
  vals    = self->val + CSC->offset;
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

void
TP_schur_matrix_update_U(TP_schur_matrix self, TP_matrix U,
			 int *row_perm, int *col_perm, int nb_pivots)
{
  int pivot, nb_threads = self->nb_threads;
  int nb_steps = (nb_pivots + nb_threads - 1 ) / nb_threads, step; 
  int U_input_size = U->m;

  for(pivot = 0; pivot < nb_pivots; pivot++)  {
    int nb_elem = self->CSR[row_perm[pivot]].nb_elem;
    int new_row_ptr = U->row_ptr[U->m] + nb_elem;
    if ( new_row_ptr > U->allocated)
      TP_matrix_realloc(U);

    U->m++;
    U->row_ptr[U->m] = new_row_ptr;
    U->nnz += nb_elem;
  }


#pragma omp parallel num_threads(nb_threads) private(step)
  for(step = 0; step < nb_steps; step++) 
    {
      int me =  omp_get_thread_num();
      struct CSR_struct *CSR;
      int i, nb_elem, U_current_row, row, col;
      int *cols, current_pivot = step * nb_threads + me ;

      if (current_pivot < nb_pivots ) { 
	col     = col_perm[current_pivot];
	row     = row_perm[current_pivot];
	CSR     = &self->CSR[row];
	nb_elem = CSR->nb_elem;
	cols    = self->col + CSR->offset;
	U_current_row = U->row_ptr[U_input_size + current_pivot];
  
	for(i = 0; i < nb_elem; i++)
	  {
	    int current_col = cols[i];
	    double val;
	    /* PROTEC THIS TEST WITH SOME KIND OF DEBUG */
	    if (current_col == col )
	      TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "this entry should of been already copied to the D");

	    pthread_mutex_lock(&self->col_locks[current_col]);
	    val = delete_entry_from_CSC(self, current_col, row);
	    pthread_mutex_unlock(&self->col_locks[current_col]);

	    U->col[U_current_row  ] = current_col;
	    U->val[U_current_row++] = val;
	  }
	
	/* TODO: recycle the row's memory */
	CSR->nb_elem = 0;
	CSR->nb_free = 0;
	self->nnz   -= nb_elem;
      }
    }
}

void
TP_schur_matrix_add_to_entry(TP_schur_matrix self, int row, int col, double val)
{
  struct CSC_struct *CSC = &self->CSC[col];
  int i;
  int nb_elem  = CSC->nb_elem;
  int *rows    = self->row + CSC->offset;
  double *vals = self->val + CSC->offset;

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
  struct CSC_struct *CSC = &self->CSC[col];
  struct CSR_struct *CSR = &self->CSR[row];
  double       *CSC_vals = self->val + CSC->offset;
  int          *CSC_rows = self->row + CSC->offset;
  int          *CSR_cols = self->col + CSR->offset;

  if (CSC->nb_free <= 0) {
    self->unused_CSC = CSC_find_free_memory(self->unused_CSC, CSC, 
					    CSC->nb_elem * 2, 
					    &self->val, &self->row, &self->allocated_CSC);
    CSC_vals = self->val + CSC->offset;
    CSC_rows = self->row + CSC->offset;
    
    if ( self->debug & (TP_DEBUG_GOSSIP_GIRL | TP_DEBUG_GARBAGE_COLLECTOR))
      TP_print_single_GB(self->unused_CSC, "CSC GB: after reallocing col");
  }
  
  if (CSR->nb_free <= 0) {
    self->unused_CSR = CSR_find_free_memory(self->unused_CSR, CSR,
					    CSR->nb_elem * 2, 
					    &self->col, &self->allocated_CSR);
    CSR_cols = self->col + CSR->offset;
    
    if ( self->debug & (TP_DEBUG_GOSSIP_GIRL | TP_DEBUG_GARBAGE_COLLECTOR))
      TP_print_single_GB(self->unused_CSR, "CSR GB: after reallocing col");
  }
  
  CSC_vals[CSC->nb_elem  ] = val;
  CSC_rows[CSC->nb_elem++] = row;
  
  CSR_cols[CSR->nb_elem++] = col;
  
  CSC->nb_free--;
  CSR->nb_free--;
  self->nnz++;
}


void
TP_schur_matrix_update_colmax(TP_schur_matrix self)
{
  int n = self->n, nb_threads = self->nb_threads;
  int nb_steps = ( n + nb_threads - 1 ) / nb_threads, step;

#pragma omp parallel num_threads(nb_threads) private(step)
  {
    int me =  omp_get_thread_num();
    for( step = 0; step < nb_steps; step++) { 
      int my_col = step * nb_threads + me; 
      if (my_col < n)
	self->CSC[my_col].col_max = get_max_double(self->val + self->CSC[my_col].offset, self->CSC[my_col].nb_elem);
    }
  }
}


void
TP_schur_matrix_update_S(TP_schur_matrix S, TP_matrix L, TP_matrix U, int start, int end)
{
  int n = S->n;
  int  i, j, k, count;
  long schur_row_struct[n];
  long *L_col_ptr = L->col_ptr;
  int  *L_rows    = L->row;
  double *L_vals  = L->val;
  long *U_row_ptr = U->row_ptr;
  int  *U_cols    = U->col;
  double *U_vals  = U->val;
  
  for ( i = 0; i < n; i++)
    schur_row_struct[i] = -1;

/* THIS IS A BAD IMPLEMENTATION OF UPDATE S */
/* #pragma omp parallel num_threads(nb_threads) */
/*   { */
/* #pragma omp single */
/*     { */
      for( count = 0, i = start; i < end; i++, count++) 
	{
	  for(j = U_row_ptr[i]; j < U_row_ptr[i+1]; j++)
	    {
	      int col        = U_cols[j];
	      double col_val = U_vals[j];
	      
	      struct CSC_struct *CSC = &S->CSC[col];
	      int *S_rows            = S->row + CSC->offset;
	      int S_col_nb_elem      = CSC->nb_elem;
	      
	      for ( k = 0; k < S_col_nb_elem; k++)
		schur_row_struct[S_rows[k]] = count*n + col;	  
	      for ( k = L_col_ptr[i]; k < L_col_ptr[i+1]; k++) { 
		if (schur_row_struct[L_rows[k]] == count*n + col) {
/* #pragma omp task */
/* 		  { */
		    int task_row = L_rows[k], task_col = col; 
		    double val = -col_val * L_vals[k];
		    TP_schur_matrix_add_to_entry(S, task_row, task_col, val);
		  /* } */
		} else {
		  TP_schur_matrix_insert_entry(S, L_rows[k], col, -col_val * L_vals[k]);
		}
	      }
	    } // for j 
	}
  /*   } // single */
  /* } // parallel */
  
  TP_schur_matrix_update_colmax(S);
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
      struct CSR_struct *CSR = &S->CSR[i];
      if ( !CSR->nb_elem )
	continue;
      invr_row[i] = k;
      self->original_rows[k++] = i;
    }
  
  for(col = 0, k=0;  col < n; col++)
    {
      struct CSC_struct *CSC = &S->CSC[col];
      int nb_elem = CSC->nb_elem;
      if ( !nb_elem ) 
	continue;
      self->original_cols[k] =  col;
      double *CSC_vals = S->val + CSC->offset;
      int    *CSC_rows = S->row + CSC->offset;
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
      struct CSC_struct *CSC = &self->CSC[i];
      int *rows    = self->row + CSC->offset;
      double *vals = self->val + CSC->offset;
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
      struct CSR_struct *CSR = &self->CSR[i];
      int *cols = self->col + CSR->offset;
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
  free_space unused_CSC = self->unused_CSC;
  free_space unused_CSR = self->unused_CSR;
  
  free(self->CSC);
  free(self->CSR);

  while (unused_CSC) {
    free_space tmp = unused_CSC->next;
    free(unused_CSC);
    unused_CSC = tmp;
  }
  while (unused_CSR) {
    free_space tmp = unused_CSR->next;
    free(unused_CSR);
    unused_CSR = tmp;
  }
  
  for( i = 0; i < self->m; i++)
    pthread_mutex_destroy(&self->row_locks[i]);
  for( i = 0; i <  self->n; i++)
    pthread_mutex_destroy(&self->col_locks[i]);
  free(self->row_locks);
  free(self->col_locks);

  free(self->val);
  free(self->row);
  free(self->col);

  free(self);
}


void
TP_schur_matrix_check_pivots(TP_schur_matrix self,
			     int *row_perms, int *col_perms,
			     int *invr_row_perms, int *invr_col_perms,
			     int nb_pivots)
{
  int  i;
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
  
  for(i = 0; i < n; i++)
    {
      struct CSC_struct *CSC = &self->CSC[i];
      free_space unused_CSC  = self->unused_CSC;
      long CSC_begin = CSC->offset, CSC_end = CSC->offset + CSC->nb_elem + CSC->nb_free;
      int current_unused = 0, j;
      if (CSC_begin == CSC_end)
	continue;
      
      while(unused_CSC)
	{
	  long free_begin = unused_CSC->offset, free_end = unused_CSC->offset + unused_CSC->nb_elem;
	  int print = 0;
	  
	  TP_overlaps overlaped = check_overalping_regions(free_begin, free_end, CSC_begin, CSC_end);
	  switch (overlaped) {
	  case (TP_overlap_none) :
	    break;
	  case (TP_overlap_begin) :
	    snprintf(mess, 2048, "The %d^th free space and the %d^th column are overlaping in the begining of the col (col start %ld end %ld, free starts on %ld ends on %ld).",
		     current_unused, i, CSC_begin, CSC_end, free_begin, free_end);
	    print = 1;
	    break;
	  case (TP_overlap_end) :
	    snprintf(mess, 2048, "The %d^th free space and the %d^th column are overlaping in the end of the col (col start %ld end %ld, free starts on %ld ends on %ld).",
		     current_unused, i, CSC_begin, CSC_end, free_begin, free_end);
	    print = 1;
	    break;
	  case (TP_overlap_total) :
	    snprintf(mess, 2048, "The %d^th free space and the %d^th column are overlapping (col starts at %ld and ends on %ld; free starts on %ld and ends on %ld).",
		     current_unused, i, CSC_begin, CSC_end, free_begin, free_end);
	    print = 1;
	    break;
	   default:
	     break;
	  }
	  if (print)
	    TP_warning(__FUNCTION__, __FILE__, __LINE__, mess);
	  current_unused++;
	  unused_CSC = unused_CSC->next;
	}
      
      for(j = 0; j < n; j++)
	{
	  if ( i == j)
	    continue;
	  struct CSC_struct *current_CSC = &self->CSC[j];
	  long current_CSC_begin = current_CSC->offset, current_CSC_end = current_CSC->offset + current_CSC->nb_elem + current_CSC->nb_free;
	  int print = 0;
	  if( current_CSC_begin == current_CSC_end)
	    continue;
	  
	  TP_overlaps overlaped = check_overalping_regions(current_CSC_begin, current_CSC_end, CSC_begin, CSC_end);
	  switch (overlaped) {
	  case (TP_overlap_none) :
	    break;
	  case (TP_overlap_begin) :
	    snprintf(mess, 2048, "The %d^th column and the %d^th column are overlaping in the begining of the col (col start %ld end %ld; col starts %ld ends on %ld).",
		     j, i, current_CSC_begin , current_CSC_end, CSC_begin, CSC_end);
	    print = 1;
	    break;
	  case (TP_overlap_end) :
	    snprintf(mess, 2048, "The %d^th column and the %d^th column are overlaping in the end of the col (col start %ld end %ld; col starts %ld ends on %ld).",
		     j, i, current_CSC_end, current_CSC_end, CSC_begin, CSC_end);
	    print = 1;
	    break;
	  case (TP_overlap_total) :
	    snprintf(mess, 2048, "The %d^th column and the %d^th column are overlapping (col starts at %ld and ends on %ld; col starts %ld ends on %ld).",
		     j, i, current_CSC_end, current_CSC_end, CSC_begin, CSC_end);
	    print = 1;
	    break;
	  default:
	    break;
	  }
	  if (print)
	    TP_warning(__FUNCTION__, __FILE__, __LINE__, mess);
	}
    }
  
  for(i = 0; i < m; i++)
    {
      struct CSR_struct *CSR = &self->CSR[i];
      free_space unused_CSR  = self->unused_CSR;
      long CSR_begin = CSR->offset, CSR_end = CSR->offset + CSR->nb_elem + CSR->nb_free;
      int current_unused = 0, j;
      
      if (CSR_begin == CSR_end)
	continue;
      
      while(unused_CSR)
	{
	  long free_begin = unused_CSR->offset, free_end = unused_CSR->offset + unused_CSR->nb_elem;
	  int print = 0;
	  
	  TP_overlaps overlaped = check_overalping_regions(free_begin, free_end, CSR_begin, CSR_end);
	  switch (overlaped) {
	  case (TP_overlap_none)  :
	    break;
	  case (TP_overlap_begin) :
	    snprintf(mess, 2048, "The %d^th free space and the %d^th row are overlaping in the begining of the col (row start %ld end %ld, free starts on  %ld ends on %ld).",
		     current_unused, i, CSR_begin, CSR_end, free_begin, free_end);
	    print = 1;
	    break;
	  case (TP_overlap_end)   :
	    snprintf(mess, 2048, "The %d^th free space and the %d^th row are overlaping in the end of the col (row start %ld end %ld, free starts on %ld ends on %ld).",
		      current_unused, i, CSR_begin, CSR_end, free_begin, free_end);
	    print = 1;
	    break;
	  case (TP_overlap_total) :
	    snprintf(mess, 2048, "The %d^th free space and the %d^th column are overlapping (col starts at %ld and ends on %ld; free starts on %ld and ends on %ld).",
		     current_unused, i, CSR_begin, CSR_end, free_begin, free_end);
	    print = 1;
	    break;
	  default:
	    break;
	  }
	  if (print)
	    TP_warning(__FUNCTION__, __FILE__, __LINE__, mess);
	  current_unused++;
	  unused_CSR = unused_CSR->next;
	}
      
      for(j = 0; j < m; j++)
	{
	   if ( i == j)
	     continue;
	   struct CSR_struct *current_CSR = &self->CSR[j];
	   long current_CSR_begin = current_CSR->offset, current_CSR_end = current_CSR->offset + current_CSR->nb_elem + current_CSR->nb_free;
	   int print = 0;
	   if( current_CSR_begin == current_CSR_end)
	     continue;
	   
	   TP_overlaps overlaped = check_overalping_regions(current_CSR_begin, current_CSR_end, CSR_begin, CSR_end);
	   switch (overlaped) {
	   case (TP_overlap_none) :
	     break;
	   case (TP_overlap_begin) :
	     snprintf(mess, 2048, "The %d^th row and the %d^th column are overlaping in the begining of the row (row start %ld end %ld; row starts %ld ends on %ld).",
		      j, i, current_CSR_begin , current_CSR_end, CSR_begin, CSR_end);
	     print = 1;
	     break;
	   case (TP_overlap_end) :
	     snprintf(mess, 2048, "The %d^th row and the %d^th row are overlaping in the end of the row (row start %ld end %ld; row starts %ld ends on %ld).",
		      j, i, current_CSR_end, current_CSR_end, CSR_begin, CSR_end);
	     print = 1;
	     break;
	   case (TP_overlap_total) :
	     snprintf(mess, 2048, "The %d^th column and the %d^th column are overlapping (col starts at %ld and ends on %ld; col starts %ld ends on %ld).",
		      j, i, current_CSR_end, current_CSR_end, CSR_begin, CSR_end);
	     print = 1;
	     break;
	   default:
	     break;
	   }
	   if (print)
	     TP_warning(__FUNCTION__, __FILE__, __LINE__, mess);
	}
    }
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
      struct CSC_struct *CSC = &self->CSC[col];
      int *rows   = self->row + CSC->offset;
      int col_nb_elem = CSC->nb_elem;
      
      for(i = 0; i < col_nb_elem; i++)
	{
	  row = rows[i];
	  struct CSR_struct *CSR = &self->CSR[row];
	  int row_nb_elem = CSR->nb_elem;
	  int *cols = self->col + CSR->offset;
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

void
TP_print_GB(TP_schur_matrix self, char *mess)
{
  free_space CSC = self->unused_CSC;
  free_space CSR = self->unused_CSR;
  
  fprintf(stdout,"%s\n", mess);
  if (CSC->previous)
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "the first CSC free memory has a predecessor");
  if (CSR->previous)
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "the first CSR free memory has a predecessor");
  
  fprintf(stdout, "|================||=================|\n");
  fprintf(stdout, "|      CSC       ||      CSR        |\n");
  fprintf(stdout, "|================||=================|\n");
  fprintf(stdout, "|    address     ||    address      |\n");
  fprintf(stdout, "|    nb_elem     ||    nb_elem      |\n");
  fprintf(stdout, "|    offset      ||    offset       |\n");
  fprintf(stdout, "|================||=================|\n\n");
  fprintf(stdout, "|===================================|\n");
  
  while( CSC || CSR) {
    if (CSC)
      fprintf(stdout, "| %14p |", &CSC);
    else
      fprintf(stdout, "                  ");
    
    if (CSR)
      fprintf(stdout, "| %14p  |\n", &CSR);
    else
      fprintf(stdout, "\n");
    
    if (CSC)
      fprintf(stdout, "|  %11ld   |", CSC->nb_elem);
    else
      fprintf(stdout, "                  ");
    
    if (CSR)
      fprintf(stdout, "|   %11ld   |\n", CSR->nb_elem);
     else
       fprintf(stdout, "\n");
    
    if (CSC)
      fprintf(stdout, "|  %11ld   |", CSC->offset);
    else
      fprintf(stdout, "                  ");
    
    if (CSR)
      fprintf(stdout, "|   %11ld   |\n", CSR->offset);
    else
      fprintf(stdout, "\n");
    fprintf(stdout, "|===================================|\n");
    
    if (CSC)
      if (CSC->next)
	if (CSC->next->previous != CSC)
	  TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "the CSC free memory's next cell, has a predecessor different from CSC");

    if (CSR)
      if (CSR->next)
	if (CSR->next->previous != CSR)
	  TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "the CSR free memory's next cell, has a predecessor different from CSR");
    
    if (CSC)
      CSC=CSC->next;
    if (CSR)
      CSR=CSR->next;
  }
}

void
TP_print_single_GB(free_space self, char *mess)
{
  fprintf(stdout,"%s\n", mess);
  
  if (self->previous)
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "the first CSR free memory has a predecessor");
    
  
  
  fprintf(stdout, "|================|\n");
  fprintf(stdout, "|    address     |\n");
  fprintf(stdout, "|    nb_elem     |\n");
  fprintf(stdout, "|    offset      |\n");
  fprintf(stdout, "|================|\n\n");
  fprintf(stdout, "|================|\n");
  
  while(self)
    {
      fprintf(stdout, "| %14p |\n", &self);
      fprintf(stdout, "|  %11ld   |\n", self->nb_elem);
      fprintf(stdout, "|  %11ld   |\n", self->offset);
      fprintf(stdout, "|================|\n");
      
      if (self->next)
	if (self->next->previous != self)
	  TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "the free memory's next cell, has a predecessor different from CSC");
      
      self = self->next;
    }
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
      struct CSC_struct *CSC = &self->CSC[col];
      int *rows              = self->row + CSC->offset;
      int nb_elem            = CSC->nb_elem;
      for( i = 0; i < nb_elem; i++) 
	_row_count[rows[i]]++;

      int row = row_perm[pivot];
      struct CSR_struct *CSR = &self->CSR[row];
      int *cols              = self->col + CSR->offset;
      nb_elem                = CSR->nb_elem;
      for( i = 0; i < nb_elem; i++) 
	_col_count[cols[i]]++;
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
