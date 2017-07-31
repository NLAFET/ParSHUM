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

  self->row_struct = malloc((size_t) nb_threads * sizeof(*self->row_struct));
  for( i = 0; i < nb_threads; i++)
    self->row_struct[i] = malloc((size_t) n * sizeof(**self->row_struct));

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

inline double
delete_entry(TP_schur_matrix self, int col, int row)
{
  delete_entry_from_CSR(self, col, row);
  return delete_entry_from_CSC(self, col, row);
}

/* void */
/* TP_schur_matrix_update_LD_V2(TP_schur_matrix self, TP_matrix L, TP_matrix D, */
/* 			     int *row_perm, int *col_perm, int nb_pivots) */
/* { */
/*   int pivot; */
/*   int L_current = L->n; */
/*   int D_current = D->n; */
 

/*   for(pivot = 0; pivot < nb_pivots; pivot++)  { */
/*     int col     = col_perm[current_pivot]; */
/*     int row     = row_perm[current_pivot]; */
/*     struct CSC_struct *S_col = &self->CSC[col]; */
/*     struct CSC_struct *L_col = &L->cols[L_current + pivot]; */

/*     double  pivot_val = delete_entry(self, col, row); */
/*     /\* podeli a S_col so pivot_val i posle stavi a u L_col  *\/ */
/*     /\*  AAA da, i stavi go pivot_val u D *\/ */
/*   } */
/* } */

/* PETOK: vaa ke treba da prais poinaku i da go izrecikliras  rows_counters. 
   PETOK: ideata e deka u CSR, voa entry pa ke se vrate samo so dr col. razmisli ubavo za voa, 
   PETOK: sea veke ne me biva za mislenje :D :D :D :D: D: (.Y.) (.Y.)(.Y.)(.Y.)(.Y.)(.Y.)(.Y.)(.Y.)(.Y.)(.Y.)(.Y.)
   PETOK: (.Y.)(.Y.)(.Y.)(.Y.)(.Y.)(.Y.)(.Y.)(.Y.)(.Y.)  
   PETOK: mozis da go predvidis duri koku  maximum sekoj row ke se sire  u S (ako a imas strukturata na U) razmisli 
   PETOK: ubavo ubavo!!!!!! (.Y.)
*/

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
	/* TODO do the delete_entru_from_CSR  in a seperate loop maybe better???  try that option */
	for( i = L->col_ptr[L_input_size + current_pivot]; i < L->col_ptr[L_input_size + current_pivot + 1]; i++)
	  L_vals[i] /= pivot_val;
	
	/* TODO: recycle the col's memory */
	CSC->nb_elem = 0;
	CSC->nb_free = 0;
	self->nnz   -= nb_elem;
      }
    }
}

/* NAREDNO: instead of going row by row, go col by col using rows struct, 
   construct the new part of U, save the coresponding colums of L 
   and the cost of each col for the spmm */
void
TP_schur_matrix_update_U(TP_schur_matrix S, TP_U_matrix U,
			 int nb_pivots, int *row_perm,
			 TP_U_struct *U_struct, int U_new_n, int U_new_nnz)
{
  int pivot, i, j;
  int nb_threads = S->nb_threads;
  int indices[nb_threads+1];
  int nnz_part = U_new_nnz / nb_threads;
  
  if ( U_new_n <= nb_threads) {
    nb_threads = U_new_n;
    for( i = 0; i < nb_threads; i++)
      indices[i] = i;
  } else { 
    *indices = 0;
    for ( i = 1, j = 0; i < nb_threads; i++) 
      {
	int part = 0;
	while ( part < nnz_part )
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
  
#pragma omp parallel num_threads(nb_threads) private(pivot)
  {
    int me =  omp_get_thread_num();
    for( pivot = indices[me]; pivot < indices[me + 1]; pivot++)
      {
	int row = row_perm[pivot];
	
	struct CSR_struct *CSR = &S->CSR[row];
	int row_nb_elem = CSR->nb_elem;
	int *cols = S->col + CSR->offset;
	
	for( i = 0; i < row_nb_elem; i++) {
	  int current_col = cols[i];
	  U_col  *u_col = &U->col[current_col];
	  
	  

	  /* PETOK: stavi atomic_fetch_and add  za indiceot. 
	     PETOK: posle mozis duri i od L da go najdis costot, tvoe e :D :D :D */
	  u_col->row[u_col->nb_elem++] = row;
	  u_col->cost += row;
	}
	
	/* TODO: recycle the row's memory */
	CSR->nb_elem = 0;
	CSR->nb_free = 0;
	S->nnz      -= row_nb_elem;
      }
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
TP_schur_matrix_update_S(TP_schur_matrix S, TP_matrix L, TP_U_matrix U,
			 TP_U_struct *U_struct, int U_new_n, int *invr_row_perm)
{
  int n = S->n;
  int  i, j, k, l;
  long *schur_row_struct = S->row_struct[0];
  long *L_col_ptr = L->col_ptr;
  int  *L_rows    = L->row;
  double *L_vals  = L->val;
  
  for ( i = 0; i < n; i++)
    schur_row_struct[i] = -1;

  for(i  = 0; i < U_new_n; i++)
    {
      TP_U_struct U_col_struct = U_struct[i];
      int col = U_col_struct.col;
      int U_col_new = U_col_struct.nb_elem;
      int U_col_new_unchanged = U_col_new ;

      struct CSC_struct *CSC = &S->CSC[col];
      double *U_vals  = U->col[col].val + U->col[col].nb_elem - U_col_new;
      int    *U_rows  = U->col[col].row + U->col[col].nb_elem - U_col_new;

      int S_col_nb_elem      = CSC->nb_elem;
      int *S_rows            = S->row + CSC->offset;
      double *S_vals         = S->val + CSC->offset;
      
      /* updating schur_row_struct and updating the vals of U  */
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
	    /* MAKE THIS FNCTION TELL U IF STHG NEEDS AN UPDATE */
            TP_schur_matrix_insert_entry(S, row, col, val);
	    S_vals         = S->val + CSC->offset;
	    S_rows         = S->row + CSC->offset;
	    /* maybe do not use S_col_nb_elem, cause is just updates it. 
	       it is maybe better just to use CSC->nb_elem*/
	    S_col_nb_elem   = CSC->nb_elem;
            schur_row_struct[row] = (long) col*n +  S_col_nb_elem - 1;
          }
        }
      }
    
    }

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
  
  for( i = 0; i < self->nb_threads; i++)
    free(self->row_struct[i]);
  free(self->row_struct);
  
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
TP_schur_check_doubles(TP_schur_matrix self)
{
  int col, row, i, j, n = self->n, m = self->m;
  char mess[2048];

  for(col = 0; col < n; col++)
    {
      struct CSC_struct *CSC = &self->CSC[col];
      int *rows = self->row + CSC->offset;
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
      struct CSR_struct *CSR = &self->CSR[row];
      int *cols = self->col + CSR->offset;
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
