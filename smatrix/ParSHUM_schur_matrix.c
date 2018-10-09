#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <strings.h>
#include <math.h>
#include <omp.h>
#include "ParSHUM_auxiliary.h"
#include "ParSHUM_schur_matrix.h"


ParSHUM_schur_matrix
ParSHUM_schur_matrix_create()
{
  ParSHUM_schur_matrix self = calloc((size_t) 1, sizeof(*self));

  return self;
}

void
ParSHUM_schur_matrix_allocate(ParSHUM_schur_matrix self, int n, int m, long nnz, int debug,
			      ParSHUM_verbose verbose, int nb_threads, double extra_space,
			      double extra_space_inbetween)
{
  int i, larger_size;

  self->nb_threads = nb_threads;
  self->n = n;
  self->m = m;
  self->verbose = verbose;
  larger_size = n > m ? n : m;

  self->CSC = calloc( (size_t) n, sizeof(*self->CSC));
  self->CSR = calloc( (size_t) m, sizeof(*self->CSR));

  self->extra_space = extra_space_inbetween;

  self->internal_mem = ParSHUM_internal_mem_create((long) nnz * (1 + extra_space + extra_space_inbetween));

  self->nnz   = 0;
  self->debug = debug;

  self->data_struct = malloc((size_t) nb_threads * sizeof(*self->data_struct));
  for( i = 0; i < nb_threads; i++) 
    self->data_struct[i] = calloc((size_t) larger_size, sizeof(**self->data_struct));
  self->base = malloc((size_t) nb_threads * sizeof(*self->base));
  int_array_memset(self->base, 1, nb_threads);

  self->row_locks = malloc((size_t) m * sizeof(*self->row_locks));
  for( i = 0; i < m; i++)
    omp_init_lock(&self->row_locks[i]);

  self->col_locks = malloc((size_t) m * sizeof(*self->col_locks));
  for( i = 0; i < m; i++)
    omp_init_lock(&self->col_locks[i]);
}

void
ParSHUM_schur_get_singletons(ParSHUM_schur_matrix self, int done_pivots, int previous_step_pivots,
			     int *nb_col_singletons, int *nb_row_singletons, int *cols, int *rows,
			     int *distributions, int nb_doane_pivots, int *col_perm, int *row_perm, 
			     int *invr_col_perm, int *invr_row_perm)
{
  int n = self->n - done_pivots + previous_step_pivots;
  int m = self->m - done_pivots + previous_step_pivots;
  int i, _done_pivots = done_pivots;
  int needed_pivots = self->n < self->m ? self->n : self->m;
  needed_pivots -= done_pivots;
  int _nb_col_singletons = 0, _nb_row_singletons = 0;

  for(i = 0; i < m; )  {
    int row = rows[i];

    if (invr_row_perm[row] != ParSHUM_UNUSED_PIVOT) { 
      rows[i] = rows[--m];
      continue;
    }
    
    /* We need the first check only for rectangular matrices  */
    /* if ( _nb_row_singletons < needed_pivots   && */
    /* 	 self->CSR[row].nb_elem == 1 ) */
    /*   { */
    /* 	int col = self->CSR[row].col[0]; */
    /* 	if (invr_col_perm[col] == ParSHUM_UNUSED_PIVOT) { */
    /* 	  row_perm[done_pivots + _nb_row_singletons] = row; */
    /* 	  invr_col_perm[col]  = _nb_row_singletons + done_pivots; */
    /* 	  col_perm[done_pivots +_nb_row_singletons] = col; */
    /* 	  invr_row_perm[row] = _nb_row_singletons++ + done_pivots; */
    /* 	} */
    /*   } */
    i++;
  }
  done_pivots += _nb_row_singletons;
  if (m != self->m - _done_pivots) 
    ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, "not all the rows are taken out from rows");

  for(i = 0; i < n; ) {
    int col = cols[i];
    
    if (invr_col_perm[col] != ParSHUM_UNUSED_PIVOT )  {
      if (invr_col_perm[col] < _done_pivots) 
	cols[i] = cols[--n];
      else 
	i++;
      continue;
    }


    /* /\* We need the first check only for rectangular matrices  *\/ */
    /* if ( _nb_col_singletons + done_pivots < needed_pivots && */
    /* 	 self->CSC[col].nb_elem == 1 ) */
    /*   { */
    /* 	int row = self->CSC[col].row[0]; */
    /* 	if (invr_row_perm[row] == ParSHUM_UNUSED_PIVOT) { */
    /* 	  col_perm[done_pivots + _nb_col_singletons] = col; */
    /* 	  invr_col_perm[col]  = _nb_col_singletons + done_pivots; */
    /* 	  row_perm[done_pivots + _nb_col_singletons] = row; */
    /* 	  invr_row_perm[row] = _nb_col_singletons++ + done_pivots; */
    /* 	} */
    /*   } */
    i++;
  }

  if (n != self->n - _done_pivots) 
    ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, "not all the cols are taken out from cols");
  
  *nb_col_singletons = _nb_col_singletons;
  *nb_row_singletons = _nb_row_singletons;
}

void 
ParSHUM_schur_matrix_init_ptr(ParSHUM_schur_matrix self, long *col_ptr, int *row_sizes)
{
  int i;
  int n = self->n;
  int m = self->m;
  ParSHUM_internal_mem memory = self->internal_mem;

  // handeling the CSC part
  for(i = 0; i < n; i++)
    ParSHUM_internal_mem_CSC_alloc(memory, &self->CSC[i],
				   (long) (1 + self->extra_space) * (col_ptr[i+1] - col_ptr[i])); 
  
  // handeling the CSR part
  for(i = 0; i < m; i++)
    ParSHUM_internal_mem_CSR_alloc(memory, &self->CSR[i],
				   (long) (1 + self->extra_space) * row_sizes[i]);
}

void
ParSHUM_CSC_update_col_max(CSC_struct *CSC, double value_tol)
{
  int i;
  int nb_elem = CSC->nb_elem;
  double max = 0.0;
  double *vals = CSC->val;
  int *rows = CSC->row;

  for(i = 0; i < nb_elem; i++)  {
    double tmp = fabs(vals[i]);
    if (tmp > max) 
      max = tmp;
  }
  CSC->col_max = max;
  max *= value_tol;

  for(i = 0; i < nb_elem; ) { 
    double val = vals[i];
    if ( fabs(vals[i]) < max )  {
      int row = rows[i];
      vals[i] = vals[--nb_elem];
      vals[nb_elem] = val;
      rows[i] = rows[nb_elem];
      rows[nb_elem] = row;
    } else {
      i++;
    }
  }

  CSC->nb_numerical_eligible = i;
}


void
ParSHUM_schur_matrix_copy(ParSHUM_matrix A, ParSHUM_schur_matrix self, double value_tol)
{
  int i, j;
  int *row_sizes;
  int n = self->n;
  
  row_sizes = ParSHUM_matrix_rows_sizes(A);
  ParSHUM_schur_matrix_init_ptr(self, A->col_ptr, row_sizes);
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

      self->CSC[i].nb_elem += col_length; 
      self->CSC[i].nb_free -= col_length; 

      ParSHUM_CSC_update_col_max(&self->CSC[i], value_tol);

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
delete_entry_from_CSR(ParSHUM_schur_matrix self, int col, int row)
{
  CSR_struct *CSR;
  int i, nb_elem, found, *cols;

  CSR     = &self->CSR[row];
  cols    = CSR->col;
  nb_elem = CSR->nb_elem;
  found = 0;

  if(nb_elem < 1) {
    ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__,"tring to delete an entry in CSR with zero elems");
    if(cols[0] == col) 
    ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__,"even better, the entry is there!");
  }

  /* TODO: this is stupid, we should somehow mark the unused enteris, not move all the memory */
  /* TODO: instead od moving the rest, just move the last one to the free position */
  for(i = 0; i < nb_elem; i++) 
    if (cols[i] == col) {
	found = 1;
	break;
    }

  if ( !found )
    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"tring to delete an non existing entry in CSR");

  CSR->nb_elem--;
  CSR->nb_free++;
  cols[i] = cols[CSR->nb_elem];
}

double
delete_entry_from_CSC(ParSHUM_schur_matrix self, int col, int row)
{
  CSC_struct *CSC;
  int i, nb_elem, found, *rows;
  double *vals, return_val = NAN;

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
    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"tring to delete an non existing entry in CSC");
  
  CSC->nb_elem--;
  CSC->nb_free++;
  rows[i] = rows[CSC->nb_elem];
  vals[i] = vals[CSC->nb_elem];

  return return_val;
}

inline double
delete_entry(ParSHUM_schur_matrix self, int col, int row)
{
  delete_entry_from_CSR(self, col, row);
  return delete_entry_from_CSC(self, col, row);
}

void
ParSHUM_schur_matrix_update_LD_singeltons(ParSHUM_schur_matrix self, ParSHUM_matrix L, ParSHUM_matrix D,
					  int *row_perm, int *col_perm, int *invr_col_perm, int nb_pivots)
{
  int pivot, nb_threads = self->nb_threads;
  int nb_steps = (nb_pivots +  nb_threads - 1 ) / nb_threads, step; 

  int L_input_size = L->n;
  int D_input_size = D->n;
  
  if ( D->n + nb_pivots > D->allocated ) 
    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__, "not enought memory in D matrix. this should never happen, so something went wrong");

  for(pivot = 0; pivot < nb_pivots; pivot++)  {
    int nb_elem = self->CSC[col_perm[pivot]].nb_elem - 1;
    int new_col_ptr = L->col_ptr[L->n] + nb_elem;
    if ( new_col_ptr > L->allocated)
      ParSHUM_matrix_realloc(L);

    L->n++;
    L->col_ptr[L->n] = new_col_ptr;
    L->nnz += nb_elem;
  }
  D->n += nb_pivots;

#pragma omp parallel num_threads(nb_threads) private(step)
  {
    long S_nnz = 0;
    int me =  omp_get_thread_num();
    for(step = 0; step < nb_steps; step++) 
      {
	int current_pivot = step * nb_threads + me;
	if ( current_pivot < nb_pivots)  {
	  CSC_struct *CSC;
	  int i, nb_elem, L_current_col, row, col;
	  int *rows, *L_rows;
	  double *vals, *L_vals, pivot_val = NAN; 
     
	  col     = col_perm[current_pivot];
	  row     = row_perm[current_pivot];
	  CSC     = &self->CSC[col];
	  nb_elem = CSC->nb_elem;
	  vals    = CSC->val;
	  rows    = CSC->row;
	  L_current_col = L->col_ptr[L_input_size + current_pivot];
	  L_rows        = L->row;
	  L_vals        = L->val;
	  
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
ParSHUM_schur_matrix_update_LD(ParSHUM_schur_matrix self, ParSHUM_L_matrix L, ParSHUM_U_matrix U, ParSHUM_matrix D,
			       int *row_perm, int *col_perm, int nb_pivots, int *invr_row_perm,
			       int nb_row_singeltons, int nb_col_singeltons, void **workspace)
{
  int nb_threads = self->nb_threads;
  long S_nnz = 0;
  long L_nnz = 0;

  int L_input_size = L->n;
  int D_input_size = D->n;
  
  if ( D->n + nb_pivots > D->allocated ) 
    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__, "not enought memory in D matrix. this should never happen, so something went wrong");

  D->n += nb_pivots;
  L->n += nb_pivots;

#pragma omp parallel for num_threads(nb_threads) shared(nb_pivots,L_input_size,D_input_size,nb_row_singeltons) reduction(+:S_nnz) reduction(+:L_nnz)
    for( int current_pivot = 0; current_pivot < nb_pivots; current_pivot++) {
      /* TODO: ParSHUM_verbose_trace_start_event(verbose, ParSHUM_UPDATE_L); */
      int me =  omp_get_thread_num();
      int n = self->n;
      int m = self->m;
      int *U_rows = (int *) workspace[me];
      int *tmp = &U_rows[m];
      double *U_vals = (double *) tmp;
      CSC_struct *CSC;
      int i, nb_elem, col, row;
      int L_indice = L_input_size + current_pivot;
      int *rows;
      double *vals, *L_vals, pivot_val = NAN;
      
      col = col_perm[current_pivot];
      row = row_perm[current_pivot];
      
      if (current_pivot < nb_row_singeltons ) {
	CSC     = &self->CSC[col];
	L->col[L_indice] = *CSC;
	vals    = CSC->val;
	rows    = CSC->row;
	nb_elem = CSC->nb_elem;
	S_nnz += CSC->nb_elem;
	for ( i = 0; i <  nb_elem; i++) 
	  if (rows[i] == row)  {
	    D->val[D_input_size + current_pivot] = pivot_val = vals[i];
	    vals[i] = vals[--L->col[L_indice].nb_elem];
	    rows[i] = rows[  L->col[L_indice].nb_elem];
	    break;
	  }
      }  else if ( col < n)  { 
	CSC     = &self->CSC[col];
	L->col[L_indice] = *CSC;
	vals    = CSC->val;
	S_nnz += CSC->nb_elem;
	
	D->val[D_input_size + current_pivot] = pivot_val = vals[--L->col[L_indice].nb_elem];
      } else { 
	col_perm[current_pivot] %= n;
	col = col_perm[current_pivot];
	CSC     = &self->CSC[col];
	L->col[L_indice] = *CSC;
	nb_elem = CSC->nb_elem - 1; 
	vals    = CSC->val;
	rows    = CSC->row;
	S_nnz += CSC->nb_elem;
	
	D->val[D_input_size + current_pivot] = pivot_val = vals[--L->col[L_indice].nb_elem];
	
	U_col *U_col = &U->col[col];
	int U_col_new = 0;
	for(i = 0; i < nb_elem; )
	  {
	    int tmp_row = rows[i];
	    if ( invr_row_perm[tmp_row] != ParSHUM_UNUSED_PIVOT ) {
	      U_rows[U_col_new] = tmp_row;
	      U_vals[U_col_new++] = vals[i];
	      vals[i] = vals[--nb_elem];
	      rows[i] = rows[  nb_elem];
	      } else {
	      i++;
	    }
	  }
	L->col[L_indice].nb_elem = nb_elem;
	
	while(U_col->allocated - U_col->nb_elem  < U_col_new)
	  ParSHUM_U_col_realloc(U_col);
	memcpy(&U_col->val[U_col->nb_elem], U_vals, U_col_new * sizeof(*U_vals));
	  memcpy(&U_col->row[U_col->nb_elem], U_rows, U_col_new * sizeof(*U_rows));
	  U_col->nb_elem += U_col_new;
      }
      
      L_vals = L->col[L_indice].val;
      L_nnz += nb_elem;
      nb_elem = L->col[L_indice].nb_elem;
      
      for( i = 0; i < nb_elem; i++)
	L_vals[i] /= pivot_val;
      
      /* TODO: recycle the col's memory */
      CSC->nb_elem = 0;
      CSC->nb_free = 0;
    }
    /* TODO: ParSHUM_verbose_trace_stop_event(verbose); */
    self->nnz -= S_nnz;
    L->nnz += S_nnz;
}

void
ParSHUM_schur_matrix_update_U_singletons(ParSHUM_schur_matrix S, ParSHUM_U_matrix U, 
					 ParSHUM_matrix D, ParSHUM_matrix L, int nb_pivots,
					 int *col_perm, int *row_perm)
{
  int nb_threads = S->nb_threads, d, sthg = L->col_ptr[L->n];
  int nb_steps = ( nb_pivots + nb_threads -1 ) / nb_threads, step;
  if ( D->n + nb_pivots > D->allocated ) 
    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__, "not enought memory in D matrix. this should never happen, so something went wrong");
 
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
	    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__, "The pivot is not row singelton");
	  }
	  if (CSC->row[0] != row) 
	    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__, "The pivot is not the same as before");
	  /* D_indice = __atomic_fetch_add(&D->n, 1, __ATOMIC_SEQ_CST); */
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
		ParSHUM_U_col_realloc(u_col);
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
#pragma omp atomic
	  U->nnz += row_n;
	}
      }
  }
}

void
ParSHUM_schur_matrix_update_U(ParSHUM_schur_matrix S, ParSHUM_U_matrix U,
			      ParSHUM_matrix L, int nb_pivots, int *row_perm,
			      ParSHUM_U_struct *U_struct, int U_new_n, int U_new_nnz)
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
      ParSHUM_U_col_realloc(u_col);
  }
  U->nnz += U_new_nnz;
  
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
	  /* int indice = __atomic_fetch_add(&u_col->nb_elem, 1, __ATOMIC_SEQ_CST); */
	  
	  /* u_col->row[indice] = row; */
	}
	
	CSR->nb_elem = 0;
	CSR->nb_free = 0;
	S_nnz    += row_nb_elem;
      }
#pragma omp atomic
    S->nnz -= S_nnz;
  }
}

void
ParSHUM_schur_matrix_add_to_entry(ParSHUM_schur_matrix self, int row, int col, double val)
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

  ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__, "entry not found");
}

void
ParSHUM_CSC_update_col_max_array(CSC_struct *CSC, double *vals, int *rows, 
				 int nb_elem, double value_tol)
{
  int i, nb_eligeble = 0, nb_ineligible = nb_elem; 
  int *CSC_rows = CSC->row;
  double *CSC_vals = CSC->val;
  double max = 0.0;
  
  for( i = 0; i < nb_elem; i++) {
    double tmp = fabs(vals[i]);
    if (tmp > max) 
      max = tmp;
  }
  
  CSC->col_max = max;
  max *= value_tol;

  for(i = 0; i < nb_elem; i++) { 
    double val = vals[i];
    if ( fabs(val) >= max )  {
      CSC_rows[nb_eligeble  ] = rows[i];
      CSC_vals[nb_eligeble++] = val;
    } else {
      CSC_rows[--nb_ineligible] = rows[i];
      CSC_vals[  nb_ineligible] = val;
    }
  }

  CSC->nb_numerical_eligible = nb_eligeble;
  CSC->nb_elem  = nb_elem;
  CSC->nb_free -= nb_elem;
};

void
ParSHUM_schur_matrix_update_S(ParSHUM_schur_matrix S, ParSHUM_L_matrix L, ParSHUM_U_matrix U,
			      int *U_struct, int U_new_n, int *invr_row_perm, int nb_pivots,
			      int *row_perms, void **workspace, double value_tol)
{
  int nb_threads = S->nb_threads;
  long S_new_nnz = 0;
  long U_new_nnz = 0;
  
#pragma omp  parallel for num_threads(nb_threads) reduction(+:S_new_nnz) reduction(+:U_new_nnz)
  for ( int i = 0; i < U_new_n; i++) {
  int me =  omp_get_thread_num();
  int m = S->m;
  int k, l;
  int *schur_row_struct = S->data_struct[me];
  int base = S->base[me];
  int *tmp_rows = (int *) workspace[me];
  int *tmp = &tmp_rows[m];
  double *tmp_vals = (double *) tmp;
  
  int col = U_struct[i];
  U_col *U_col = &U->col[col];
  
  CSC_struct *CSC = &S->CSC[col];
  int S_col_nb_elem = CSC->nb_elem;
  int U_nb_elem     = m;
  int S_nb_elem     = 0;
  int *S_rows       = CSC->row;
  double *S_vals    = CSC->val;
  int *U_rows;
  double *U_vals;
  long needed_size;
  
  for ( k = 0; k < S_col_nb_elem; k++) {
    int S_row = S_rows[k];
    if (invr_row_perm[S_row] != ParSHUM_UNUSED_PIVOT ) {
      tmp_rows[--U_nb_elem] = S_row;
      tmp_vals[  U_nb_elem] = S_vals[k];
    } else { 
      tmp_rows[S_nb_elem] = S_row;
      tmp_vals[S_nb_elem] = S_vals[k];
      schur_row_struct[S_row] = base + S_nb_elem++;
    }
  }

  U_vals = &tmp_vals[U_nb_elem];
  U_rows = &tmp_rows[U_nb_elem];
  S_vals = tmp_vals;
  S_rows = tmp_rows;
  U_nb_elem = m - U_nb_elem; 
  
  if (S_col_nb_elem != ( U_nb_elem + S_nb_elem ))
    ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__,"Something went wrong in the schur update");
  
  U_new_nnz += (long) U_nb_elem;
	
  for ( l = 0; l < U_nb_elem; l++) {
    int L_col  = invr_row_perm[U_rows[l]];
    CSC_struct *L_CSC = &L->col[L_col];
    double U_val = U_vals[l];
    
    int    *L_rows = L_CSC->row;
    double *L_vals = L_CSC->val;
    int  L_nb_elem = L_CSC->nb_elem; 
    
    for ( k = 0; k < L_nb_elem; k++) {
      int row = L_rows[k];
      double val = -U_val * L_vals[k];
      int indice = schur_row_struct[row] - base;
      if ( indice  >= 0 )  {
	S_vals[indice] += val;
      } else {
	S_vals[S_nb_elem] = val;
	S_rows[S_nb_elem] = row;
	schur_row_struct[row] = base + S_nb_elem++;
      }
    }
  }
  
  while(U_col->allocated - U_col->nb_elem  < U_nb_elem)
    ParSHUM_U_col_realloc(U_col);
  memcpy(&U_col->val[U_col->nb_elem], U_vals, U_nb_elem * sizeof(*U_vals));
  memcpy(&U_col->row[U_col->nb_elem], U_rows, U_nb_elem * sizeof(*U_rows));
  
  U_col->nb_elem += U_nb_elem;
  S_new_nnz += (long) S_nb_elem - (long) CSC->nb_elem;

  CSC->nb_free += CSC->nb_elem;
  CSC->nb_elem = 0;
  needed_size = (long) CSC->nb_free;

  if (needed_size < (long) S_nb_elem)  {
    while( needed_size < (long) S_nb_elem  )
      needed_size *= 2;
    ParSHUM_internal_mem_CSC_alloc(S->internal_mem, CSC, needed_size);
  }
  
  ParSHUM_CSC_update_col_max_array(CSC, S_vals, S_rows, S_nb_elem, value_tol);
  
  int new = base + S_nb_elem;
  if (new < base ) {
    base = 1;
    bzero((void *) schur_row_struct, (size_t) S->n * sizeof(*schur_row_struct));
  } else {
    base = new;
  }
  S->base[me] = base;
  } // for

  S->nnz += S_new_nnz;
  U->nnz += U_new_nnz;
}


void
ParSHUM_schur_matrix_update_S_rows(ParSHUM_schur_matrix S, int *L_struct, int L_new_n,
				   int L_new_nnz, int *invr_col_perm, int nb_pivots,
				   int *row_perms, int done_pivots, void **workspace)
{
  int nb_threads = S->nb_threads;

#pragma omp parallel for num_threads(nb_threads) 
  for ( int i = 0; i < L_new_n; i++) {
  int me = omp_get_thread_num();
  int k, l;
  int *schur_row_struct = S->data_struct[me];
  int base = S->base[me];
  
  int row = L_struct[i];
  
  CSR_struct *CSR = &S->CSR[row];
  int n = S->n;
  int *tmp = workspace[me];
  int S_row_nb_elem = CSR->nb_elem;
  int *S_cols       = CSR->col;
  int *L_cols;
  int S_nb_elem = 0 , L_nb_elem = n ;
  long needed_size;
  
  /* constructing schur_row_struct and discovering the vals of U  */
  for (k = 0; k < S_row_nb_elem; k++) {
    int S_col = S_cols[k];
    if (invr_col_perm[S_col] != ParSHUM_UNUSED_PIVOT ) {
      tmp[--L_nb_elem] = S_col;
    } else {
      tmp[S_nb_elem] = S_col;
      schur_row_struct[S_col] = base + S_nb_elem++;
    }
  }
  L_cols = &tmp[L_nb_elem];
  S_cols = tmp;
  L_nb_elem = n - L_nb_elem;
  
  for ( l = 0; l < L_nb_elem; l++) {
    int U_row = row_perms[invr_col_perm[L_cols[l]]];
    int *U_cols   = S->CSR[U_row].col;
    int U_nb_elem = S->CSR[U_row].nb_elem;
    
    for ( k = 0; k < U_nb_elem; k++) {
      int col = U_cols[k];
      int indice = schur_row_struct[col] - base;
      
      if (invr_col_perm[col] != ParSHUM_UNUSED_PIVOT) {
	continue;
      }
      
      if ( indice  < 0 )  {
	S_cols[S_nb_elem] = col;
	schur_row_struct[col] = base + S_nb_elem++;
      }
    }
  }
  CSR->nb_free += CSR->nb_elem;
  CSR->nb_elem = 0;
  needed_size = (long) CSR->nb_free;
  
  if (needed_size < (long) S_nb_elem)  {
    while( needed_size < (long) S_nb_elem  )
      needed_size *= 2;
    ParSHUM_internal_mem_CSR_alloc(S->internal_mem, CSR, needed_size);
  }
  memcpy(CSR->col, S_cols, S_nb_elem * sizeof(*S_cols));
  CSR->nb_free -= S_nb_elem;
  CSR->nb_elem += S_nb_elem;
  
  int new = base + CSR->nb_elem;
  if (new < base ) {
    base = 1;
    bzero((void *) schur_row_struct, (size_t) S->n * sizeof(*schur_row_struct));
  } else {
    base = new;
  }
  S->base[me] = base;
  }

  int start = done_pivots, end = done_pivots + nb_pivots;

#pragma omp parallel for num_threads(nb_threads) shared(start, end)
  for(int i = start; i < end; i++)  {
    S->CSR[row_perms[i]].nb_free += S->CSR[row_perms[i]].nb_elem;
    S->CSR[row_perms[i]].nb_elem = 0;
  }
}

ParSHUM_dense_matrix 
ParSHUM_schur_matrix_convert(ParSHUM_schur_matrix S, int *invr_col_perm,
			     int *invr_row_perm, int done_pivots)
{
  ParSHUM_dense_matrix self; 
  int col, k, i;
  int n = S->n;
  int m = S->m;
  int n_schur = n - done_pivots;
  int m_schur = m - done_pivots;
  int invr_row[m];
  int bb = 0 ;
  self = ParSHUM_dense_matrix_create(n_schur, m_schur);

  for(i = 0, k=0; i < m; i++)
    if (invr_row_perm[i] == ParSHUM_UNUSED_PIVOT ) { 
      invr_row[i] = k;
      self->original_rows[k++] = i;
    }
  
  for(col = 0, k=0;  col < n; col++)
    if (invr_col_perm[col] == ParSHUM_UNUSED_PIVOT ) { 
      CSC_struct *CSC = &S->CSC[col];

      int nb_elem = CSC->nb_elem;
      self->original_cols[k] = col;
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
ParSHUM_schur_matrix_print(ParSHUM_schur_matrix self, char *mess)
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
ParSHUM_schur_matrix_destroy(ParSHUM_schur_matrix self)
{
  int i;
  
  free(self->CSC);
  free(self->CSR);

  ParSHUM_internal_mem_destroy(self->internal_mem);

  for( i = 0; i < self->nb_threads; i++)
    free(self->data_struct[i]);
  free(self->data_struct);
  free(self->base);

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
/* ********************************************************************************************* */
/* ********************************************************************************************* */
/* ********************************************************************************************* */
void
ParSHUM_schur_check_doubles(ParSHUM_schur_matrix self)
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
	    ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
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
	    ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
	  }
      }
    }
}

void
ParSHUM_schur_matrix_check_pivots(ParSHUM_schur_matrix self,
				  int *row_perms, int *col_perms,
				  int *invr_row_perms, int *invr_col_perms,
				  int nb_pivots)
{
  int  i, j, n = self->n, m = self->m;
  char mess[2048];
  
  check_vlaid_perms(row_perms,      m, nb_pivots);
  check_vlaid_perms(invr_row_perms, m, nb_pivots);
  check_vlaid_perms(col_perms,      n, nb_pivots);
  check_vlaid_perms(invr_row_perms, n, nb_pivots);
  check_perms_and_invr_perms(col_perms, invr_col_perms, nb_pivots, "col");
  check_perms_and_invr_perms(row_perms, invr_row_perms, nb_pivots, "row");
  
  for( i = 0; i < nb_pivots; i++)
    {
      int row = row_perms[i], col = col_perms[i];
      if ( self->CSC[col].nb_elem ) {
	snprintf(mess, 2048, "column %d is a pivot, but not empty in S with nb_elem %d and nb_free %d",
		 col, self->CSC[col].nb_elem, self->CSC[col].nb_free);
	ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
      }
      
      if ( self->CSR[row].nb_elem ) {
	snprintf(mess, 2048, "row %d is a pivot, but not empty in S with nb_elem %d and nb_free %d",
		 row, self->CSR[row].nb_elem, self->CSR[row].nb_free);
	ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
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
	    if ( invr_row_perms[row] != ParSHUM_UNUSED_PIVOT ) {
	      snprintf(mess, 2048, "in col %d, %d is present, but %d is a row pivot\n", i, row, row);
	      ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
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
	  if ( invr_col_perms[col] != ParSHUM_UNUSED_PIVOT ) {
	    snprintf(mess, 2048, "in row %d, %d is present, but %d is a col pivot\n", i, col, col);
	    ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
	  }
	}
    }
}

void
ParSHUM_schur_matrix_memory_check(ParSHUM_schur_matrix self)
{
  int i, n = self->n, m = self->m;
  char mess[2048];
  
  for ( i = 0; i < n; i++)
    if (self->CSC[i].nb_free < 0 )  {
      snprintf(mess, 2048, "error on the column %d  with nb_free %d\n", i, self->CSC[i].nb_free);
      ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
    }
  
  for ( i = 0; i < m; i++)
    if (self->CSR[i].nb_free < 0 ) {
      snprintf(mess, 2048, "error on the row %d with nb_free %d\n", i, self->CSR[i].nb_free);
      ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
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
	  
  /* 	  ParSHUM_overlaps overlaped = check_overalping_regions(free_begin, free_end, CSC_begin, CSC_end); */
  /* 	  switch (overlaped) { */
  /* 	  case (ParSHUM_overlap_none) : */
  /* 	    break; */
  /* 	  case (ParSHUM_overlap_begin) : */
  /* 	    snprintf(mess, 2048, "The %d^th free space and the %d^th column are overlaping in the begining of the col (col start %ld end %ld, free starts on %ld ends on %ld).", */
  /* 		     current_unused, i, CSC_begin, CSC_end, free_begin, free_end); */
  /* 	    print = 1; */
  /* 	    break; */
  /* 	  case (ParSHUM_overlap_end) : */
  /* 	    snprintf(mess, 2048, "The %d^th free space and the %d^th column are overlaping in the end of the col (col start %ld end %ld, free starts on %ld ends on %ld).", */
  /* 		     current_unused, i, CSC_begin, CSC_end, free_begin, free_end); */
  /* 	    print = 1; */
  /* 	    break; */
  /* 	  case (ParSHUM_overlap_total) : */
  /* 	    snprintf(mess, 2048, "The %d^th free space and the %d^th column are overlapping (col starts at %ld and ends on %ld; free starts on %ld and ends on %ld).", */
  /* 		     current_unused, i, CSC_begin, CSC_end, free_begin, free_end); */
  /* 	    print = 1; */
  /* 	    break; */
  /* 	   default: */
  /* 	     break; */
  /* 	  } */
  /* 	  if (print) */
  /* 	    ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess); */
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
	  
  /* 	  ParSHUM_overlaps overlaped = check_overalping_regions(current_CSC_begin, current_CSC_end, CSC_begin, CSC_end); */
  /* 	  switch (overlaped) { */
  /* 	  case (ParSHUM_overlap_none) : */
  /* 	    break; */
  /* 	  case (ParSHUM_overlap_begin) : */
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
  /* 	    ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess); */
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
	  
  /* 	  ParSHUM_overlaps overlaped = check_overalping_regions(free_begin, free_end, CSR_begin, CSR_end); */
  /* 	  switch (overlaped) { */
  /* 	  case (ParSHUM_overlap_none)  : */
  /* 	    break; */
  /* 	  case (ParSHUM_overlap_begin) : */
  /* 	    snprintf(mess, 2048, "The %d^th free space and the %d^th row are overlaping in the begining of the col (row start %ld end %ld, free starts on  %ld ends on %ld).", */
  /* 		     current_unused, i, CSR_begin, CSR_end, free_begin, free_end); */
  /* 	    print = 1; */
  /* 	    break; */
  /* 	  case (ParSHUM_overlap_end)   : */
  /* 	    snprintf(mess, 2048, "The %d^th free space and the %d^th row are overlaping in the end of the col (row start %ld end %ld, free starts on %ld ends on %ld).", */
  /* 		      current_unused, i, CSR_begin, CSR_end, free_begin, free_end); */
  /* 	    print = 1; */
  /* 	    break; */
  /* 	  case (ParSHUM_overlap_total) : */
  /* 	    snprintf(mess, 2048, "The %d^th free space and the %d^th column are overlapping (col starts at %ld and ends on %ld; free starts on %ld and ends on %ld).", */
  /* 		     current_unused, i, CSR_begin, CSR_end, free_begin, free_end); */
  /* 	    print = 1; */
  /* 	    break; */
  /* 	  default: */
  /* 	    break; */
  /* 	  } */
  /* 	  if (print) */
  /* 	    ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess); */
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
	   
  /* 	   ParSHUM_overlaps overlaped = check_overalping_regions(current_CSR_begin, current_CSR_end, CSR_begin, CSR_end); */
  /* 	   switch (overlaped) { */
  /* 	   case (ParSHUM_overlap_none) : */
  /* 	     break; */
  /* 	   case (ParSHUM_overlap_begin) : */
  /* 	     snprintf(mess, 2048, "The %d^th row and the %d^th column are overlaping in the begining of the row (row start %ld end %ld; row starts %ld ends on %ld).", */
  /* 		      j, i, current_CSR_begin , current_CSR_end, CSR_begin, CSR_end); */
  /* 	     print = 1; */
  /* 	     break; */
  /* 	   case (ParSHUM_overlap_end) : */
  /* 	     snprintf(mess, 2048, "The %d^th row and the %d^th row are overlaping in the end of the row (row start %ld end %ld; row starts %ld ends on %ld).", */
  /* 		      j, i, current_CSR_end, current_CSR_end, CSR_begin, CSR_end); */
  /* 	     print = 1; */
  /* 	     break; */
  /* 	   case (ParSHUM_overlap_total) : */
  /* 	     snprintf(mess, 2048, "The %d^th column and the %d^th column are overlapping (col starts at %ld and ends on %ld; col starts %ld ends on %ld).", */
  /* 		      j, i, current_CSR_end, current_CSR_end, CSR_begin, CSR_end); */
  /* 	     print = 1; */
  /* 	     break; */
  /* 	   default: */
  /* 	     break; */
  /* 	   } */
  /* 	   if (print) */
  /* 	     ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess); */
  /* 	} */
  /*   } */
}

void
ParSHUM_schur_matrix_check_symetry(ParSHUM_schur_matrix self)
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
    ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
  }
  
  if (CSC_nnz != self->nnz) {
    snprintf(mess, 2048, "CSC and S nnz are not the same, CSC = %ld and S_nnz = %ld", CSC_nnz, self->nnz);
    ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
  }
  
  if (CSR_nnz != self->nnz) {
    snprintf(mess, 2048, "CSR and S nnz are not the same, CSR = %ld and S_nnz = %ld", CSR_nnz, self->nnz);
    ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
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
	    snprintf(mess, 2048, "In CSC in col %d row %d exists, but in CSR in row %d, col %d does not exist",
		     col, row, row, col);
	    ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
	  }
	} /* for I */
    }

  for(row = 0; row < n; row++)
    {
      CSR_struct *CSR = &self->CSR[row];
      int *cols   = CSR->col;
      int row_nb_elem = CSR->nb_elem;
      
      for(i = 0; i < row_nb_elem; i++)
	{
	  col = cols[i];
	  CSC_struct *CSC = &self->CSC[col];
	  int col_nb_elem = CSC->nb_elem;
	  int *rows = CSC->row;
	  int found = 0;
	  for(j = 0; j < col_nb_elem; j++) {
	    if(rows[j] == row) {
	      found = 1;
	      break;
	    }
	  }
	  
	  if (!found) {
	    snprintf(mess, 2048, "In CSR in row %d col %d exists, but in CSC in col %d, row %d does not exist",
		     row, col, col, row);
	    ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
	  }
	} /* for I */
    }
}

/* TODO: addapt this to the new thing */
void
ParSHUM_print_GB(ParSHUM_schur_matrix self, char *mess)
{
  /* free_space CSC = self->unused_CSC; */
  /* free_space CSR = self->unused_CSR; */
  
  /* fprintf(stdout,"%s\n", mess); */
  /* if (CSC->previous) */
  /*   ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__, "the first CSC free memory has a predecessor"); */
  /* if (CSR->previous) */
  /*   ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__, "the first CSR free memory has a predecessor"); */
  
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
  /* 	  ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__, "the CSC free memory's next cell, has a predecessor different from CSC"); */

  /*   if (CSR) */
  /*     if (CSR->next) */
  /* 	if (CSR->next->previous != CSR) */
  /* 	  ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__, "the CSR free memory's next cell, has a predecessor different from CSR"); */
    
  /*   if (CSC) */
  /*     CSC=CSC->next; */
  /*   if (CSR) */
  /*     CSR=CSR->next; */
  /* } */
}

/* TODO: addapt this to the new thing */
/* void */
/* ParSHUM_print_single_GB(free_space self, char *mess) */
/* { */
  /* fprintf(stdout,"%s\n", mess); */
  
  /* if (self->previous) */
  /*   ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__, "the first CSR free memory has a predecessor"); */
    
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
  /* 	  ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__, "the free memory's next cell, has a predecessor different from CSC"); */
      
  /*     self = self->next; */
  /*   } */
/* } */


void
ParSHUM_check_current_counters(ParSHUM_schur_matrix self,
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
	ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
      }

      if ( row_count[row] != base) {
	snprintf(mess, 2048, "row_count[%d] = %d, base = %d, but %d is a pivot",
		 row, row_count[i], base, row);
	ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
      }

      if ( _col_count[col] != 1) {
	snprintf(mess, 2048, "calculated col_count[%d] = %d, but %d is a pivot",
		 col, _col_count[i], row);
	ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
      }

      if ( row_count[row] != base) {
	snprintf(mess, 2048, "col_count[%d] = %d, base = %d, but %d is a pivot",
		 col, col_count[i], base, col);
	ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
      }
    }  

  for( i = 0; i < n; i++) 
    {
      if (col_count[i] >=  base) 
	if ( ( col_count[i] - base + 1 ) !=  _col_count[i]) {
	  snprintf(mess, 2048, "error on col counter %d : col_count(%d) base(%d) and calculated col_count(%d)",
		   i, col_count[i], base, _col_count[i]);
	  ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
	}
      
      if (row_count[i] >=  base) 
	if ( ( row_count[i] - base + 1 ) !=  _row_count[i]) {
	  snprintf(mess, 2048, "error on row counter %d : row_count(%d) base(%d) and calculated row_count(%d)",
		   i, row_count[i], base, _row_count[i]);
	  ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
	}
    }	 
  
  free(_col_count);
  free(_row_count);
}
