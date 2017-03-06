#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "TP_auxiliary.h"
#include "TP_schur_matrix.h"


struct _free_space_CSC {
  long   nb_elem;
  int    *row;
  double *val;
  free_space_CSC next;
  free_space_CSC previous;
};

struct _free_space_CSR {
  long   nb_elem;
  int    *col;
  free_space_CSR next;
  free_space_CSR previous;
};

TP_schur_matrix
TP_schur_matrix_create()
{
  TP_schur_matrix self = calloc((size_t) 1, sizeof(*self));

  return self;
}

void
TP_schur_matrix_allocate(TP_schur_matrix self, int n, int m, long nnz, 
			 double extra_space, double extra_space_inbetween)
{
  long allocating;

  self->n = n;
  self->m = m;

  self->CSC = calloc( (size_t) n, sizeof(*self->CSC));
  self->CSR = calloc( (size_t) m, sizeof(*self->CSR));

  self->extra_space           = extra_space;
  self->extra_space_inbetween = extra_space_inbetween;
  allocating = nnz * (1 + extra_space + extra_space_inbetween);

  self->val    = malloc(sizeof (*self->val));
  self->val[0] = malloc((size_t) allocating * sizeof(**self->val));
  self->row    = malloc(sizeof(*self->row));
  self->row[0] = malloc((size_t) allocating * sizeof(**self->row));
  self->col    = malloc(sizeof(*self->col));
  self->col[0] = malloc((size_t) allocating * sizeof(**self->col));
  self->nb_CSC_memories = 1;
  self->nb_CSR_memories = 1;

  self->unused_CSC = calloc(1, sizeof(*self->unused_CSC));
  self->unused_CSR = calloc(1, sizeof(*self->unused_CSR));
  self->unused_CSC->nb_elem = allocating;
  self->unused_CSC->row     = self->row[0];
  self->unused_CSC->val     = self->val[0];
  self->unused_CSR->nb_elem = allocating;
  self->unused_CSR->col     = self->col[0];
}


void
TP_schur_matrix_destroy(TP_schur_matrix self)
{
  int i;
  free_space_CSC unused_CSC = self->unused_CSC;
  free_space_CSR unused_CSR = self->unused_CSR;
  
  free(self->CSC);
  free(self->CSR);

  while (unused_CSC) {
    free_space_CSC tmp = unused_CSC->next;
    free(unused_CSC);
    unused_CSC = tmp;
  }
  while (unused_CSR) {
    free_space_CSR tmp = unused_CSR->next;
    free(unused_CSR);
    unused_CSR = tmp;
  }
  
  for( i = 0; i < self->nb_CSC_memories; i++)
    {
      free(self->val[i]);
      free(self->row[i]);
    }
  free(self->val);
  free(self->row);

  for( i = 0; i < self->nb_CSR_memories; i++)
    free(self->col[i]);
  free(self->col);

  free(self);
}


free_space_CSC
CSC_find_free_memory(free_space_CSC free_space,  struct CSC_struct *CSC, long nb_elem)
{
  free_space_CSC tmp = free_space;

  if (CSC->nb_free > 0 )
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"acquiring new CSR memory to non-full memory ");

  if (CSC->nb_elem > nb_elem )
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"the acquired CSR memory is smaller than the current one");
      
  while(tmp) {
    if (tmp->nb_elem  <= nb_elem) 
      tmp  = tmp->next;
    else 
      break;
  }

  if (!tmp)
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"not enought CSC memory in the Schur matrix!  TODO: implement a new allocation");

  memcpy(tmp->val, CSC->val, (size_t) CSC->nb_elem * sizeof(*CSC->val));
  memcpy(tmp->row, CSC->row, (size_t) CSC->nb_elem * sizeof(*CSC->row));

  CSC->val      = tmp->val;
  CSC->row      = tmp->row;
  CSC->nb_free  = nb_elem - CSC->nb_elem;

  if (tmp->nb_elem - nb_elem == 0) {
    if (tmp == free_space) {
      if (!tmp->next)
	TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"not enought CSC memory in the Schur matrix!  TODO: implement a new allocation");
      else
	free_space = tmp->next;
    }
    if ( !tmp->previous ) 
      tmp->previous->next = tmp->next;
    if ( !tmp->previous ) 
      tmp->next->previous = tmp->previous;
    free(tmp);
  }else {
    tmp->val     += nb_elem;
    tmp->row     += nb_elem;
    tmp->nb_elem -= nb_elem;
  }

  return free_space;
}

free_space_CSR
CSR_find_free_memory(free_space_CSR free_space,  struct CSR_struct *CSR, long nb_elem)
{
  free_space_CSR tmp = free_space;

  if (CSR->nb_free > 0 )
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"acquiring new CSR memory to non-full memory ");

  if (CSR->nb_elem > nb_elem )
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"the acquired CSR memory is smaller than the current one");
      
  while(tmp) {
    if (tmp->nb_elem  <= nb_elem) 
      tmp = tmp->next;
    else 
      break;
  }

  if (!tmp)
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"not enought CSC memory in the Schur matrix!  TODO: implement a new allocation");

  memcpy(tmp->col, CSR->col, (size_t) CSR->nb_elem * sizeof(*CSR->col));
 
  CSR->col      = tmp->col;
  CSR->nb_free  = nb_elem - CSR->nb_elem;

  if (tmp->nb_elem - nb_elem == 0) {
    if (tmp == free_space) { 
      if (!tmp->next)
	TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"not enought CSC memory in the Schur matrix!  TODO: implement a new allocation");
      else
	free_space = tmp->next;
    }
    if ( !tmp->previous ) 
      tmp->previous->next = tmp->next;
    if ( !tmp->previous ) 
      tmp->next->previous = tmp->previous;
    free(tmp);
  }else {
    tmp->col     += nb_elem;
    tmp->nb_elem -= nb_elem;
  }

  return free_space;
}


void 
TP_schur_matrix_init_ptr(TP_schur_matrix self, long *col_ptr, int *row_sizes)
{
  int i;
  int n = self->n;
  int m = self->m;
  
  // handeling the CSC part
  for(i = 0; i < n; i++)
    self->unused_CSC = CSC_find_free_memory(self->unused_CSC, &self->CSC[i], col_ptr[i+1] - col_ptr[i]);
  
  // handeling the CSR part
  for(i = 0; i < m; i++)
    self->unused_CSR = CSR_find_free_memory(self->unused_CSR, &self->CSR[i], (long) row_sizes[i]);
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

  for( i = 0; i < n; i++) 
    {
      long   A_col_start = A->col_ptr[i];
      long   A_col_end   = A->col_ptr[i+1];
      long   col_length  = A_col_end - A_col_start;
      
      // handle the copy  of the column into the CSC structure
      memcpy((void *) self->CSC[i].row,
	     (void *) &A->row[A_col_start],
	     col_length * sizeof(*A->row));

      memcpy((void *) self->CSC[i].val,
	     (void *) &A->val[A_col_start],
	     col_length * sizeof(*A->val));

      self->CSC[i].col_max = get_max_double(self->CSC[i].val, col_length);
      
      self->CSC[i].nb_elem += col_length; 
      self->CSC[i].nb_free -= col_length; 

      // handle the copy  of the column into the CSR structure
      for(j = A_col_start; j < A_col_end; j++)
	{       
	  int row = A->row[j];
	  struct CSR_struct *row_CSR = &self->CSR[row];
	  
	  row_CSR->col[row_CSR->nb_elem++] = i;
	  row_CSR->nb_free--;
	}
    }
}

void
delete_entry_from_CSR(TP_schur_matrix self, int col, int row)
{
  struct CSR_struct *CSR;
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

  for(i = 0; i < nb_elem; i++) {
    if (!found)  {
      if (cols[i] == col) 
	found = 1;
    } else {
      cols[i-1] = cols[i];
    }
  }

  if ( !found )
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"tring to delete an non existing entry in CSR");

  CSR->nb_elem--;
  CSR->nb_free++;
}

free_space_CSC
add_unused_CSC(free_space_CSC unused_space, long nb_elem, int *row, double *val)
{
  free_space_CSC self = malloc(sizeof(*self));

  self->nb_elem  = nb_elem;
  self->row      = row;
  self->val      = val;
  self->next     = unused_space;
  self->previous = NULL;

  return self;
}

void
TP_schur_matrix_update_LD(TP_schur_matrix self, TP_matrix L, TP_matrix D,
			  int *row_perm, int *col_perm, int nb_pivots)
{
  int pivot; 
  for(pivot = 0; pivot < nb_pivots; pivot++) 
    {
      struct CSC_struct *CSC;
      int i, nb_elem, L_current_col, row, col;
      int *rows, *L_rows;
      double *vals, *L_vals, pivot_val;
      
      col     = col_perm[pivot];
      row     = row_perm[pivot];
      CSC     = &self->CSC[col];
      nb_elem = CSC->nb_elem;
      vals    = CSC->val;
      rows    = CSC->row;
      L_current_col = L->col_ptr[L->n];
      L_rows        = L->row;
      L_vals        = L->val;

      if ( D->n + 1 > D->allocated ) 
	TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "not enought memory in D matrix. this should never happen, so something went wrong");
      if ( L->nnz + nb_elem > L->allocated ) 
	TP_matrix_realloc(L, self->extra_space + self->extra_space_inbetween + 1);
  
      for(i = 0; i < nb_elem; i++)
	{
	  if ( rows[i] != row) {
	    L_rows[L_current_col  ] = rows[i];
	    L_vals[L_current_col++] = vals[i];
	  } else {
	    D->val[D->n++] = pivot_val = vals[i];
	  }
	  delete_entry_from_CSR(self, col, rows[i]);
	}
      
      L->n++;
      L->col_ptr[L->n] = L_current_col;
      L->nnz += nb_elem;
      
      // dividing the L part with the chosen pivot
      //TODO: do this in a more clever way
      for( i = L->col_ptr[L->n-1]; i < L->col_ptr[L->n]; i++)
	L_vals[i] /= pivot_val;
      
      self->unused_CSC = add_unused_CSC(self->unused_CSC, (long ) nb_elem + CSC->nb_free, rows, vals);
      CSC->nb_elem = 0;
      CSC->nb_free = 0;
      CSC->row     = NULL;
      CSC->val     = NULL;
    }
}
 

free_space_CSR
add_unused_CSR(free_space_CSR unused_space, long nb_elem, int *col)
{
  free_space_CSR self = malloc(sizeof(*self));

  self->nb_elem  = nb_elem;
  self->col      = col;
  self->next     = unused_space;
  self->previous = NULL;

  return self;
}

double
delete_entry_from_CSC(TP_schur_matrix self, int col, int row)
{
  struct CSC_struct *CSC;
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
      return_val = vals[i++];
      break;
    }

  if ( !found ) 
    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__,"tring to delete an non existing entry in CSC");

  // move one elem back the rest of the arrays (couldn't use memcpy because of overlaping areas and
  // donr want to move memmove  beacuse it uses and additional vector)
  for( ; i < nb_elem; i++) {
    rows[i-1] = rows[i];
    vals[i-1] = vals[i];
  }
  
  CSC->nb_elem--;
  CSC->nb_free++;
  
  return return_val;
}

void
TP_schur_matrix_update_U(TP_schur_matrix self, TP_matrix U,
			 int *row_perm, int *col_perm, int nb_pivots)
{
  int pivot; 
  for(pivot = 0; pivot < nb_pivots; pivot++) 
    {
      struct CSR_struct *CSR;
      int i, nb_elem, U_current_row, row, col;
      int *cols;
      
      col     = col_perm[pivot];
      row     = row_perm[pivot];
      CSR     = &self->CSR[row];
      nb_elem = CSR->nb_elem;
      cols    = CSR->col;
      U_current_row = U->row_ptr[U->m];

      if ( U->nnz + nb_elem > U->allocated ) 
	TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "not enought memory in U matrix. TODO: implement realloc");
  
      for(i = 0; i < nb_elem; i++)
	{
	  int current_col = cols[i];
	  double val;
	  if (current_col == col )
	    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "this entry should of been already copied to the D");
	  
	  val = delete_entry_from_CSC(self, current_col, row);
	  U->col[U_current_row  ] = current_col;
          U->val[U_current_row++] = val;
	}

      U->m++;
      U->row_ptr[U->m] = U_current_row;
      U->nnz += nb_elem;

      self->unused_CSR = add_unused_CSR(self->unused_CSR, (long ) nb_elem + CSR->nb_free, cols);
      CSR->nb_elem = 0;
      CSR->nb_free = 0;
      CSR->col     = NULL;
    }
}

void
TP_schur_matrix_add_to_entry(TP_schur_matrix self, int row, int col, double val)
{
  struct CSC_struct *CSC = &self->CSC[col];
  int i;
  int *rows    = CSC->row;
  int nb_elem  = CSC->nb_elem;
  double *vals = CSC->val;

  for( i = 0; i < nb_elem; i++)
    if ( rows[i] == row) { 
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
  
  if (CSC->nb_free <= 0)
    self->unused_CSC = CSC_find_free_memory(self->unused_CSC, CSC, CSC->nb_elem * (1 + self->extra_space_inbetween));
  if (CSR->nb_free <= 0)
    self->unused_CSR = CSR_find_free_memory(self->unused_CSR, CSR, CSR->nb_elem * (1 + self->extra_space_inbetween));

  
  CSC->val[CSC->nb_elem  ] = val;
  CSC->row[CSC->nb_elem++] = row;

  CSR->col[CSR->nb_elem++] = col;
  
  CSC->nb_free--;
  CSR->nb_free--;
}


void
TP_schur_matrix_update_colmax(TP_schur_matrix self)
{
  int n = self->n;

  for( int i = 0; i < n; i++) 
    {
      struct CSC_struct *CSC = &self->CSC[i];
      if ( CSC->nb_elem )
	CSC->col_max = get_max_double(CSC->val, CSC->nb_elem);
      else
	CSC->col_max = 0.0;
    }
}


void
TP_schur_matrix_update_S(TP_schur_matrix S, TP_matrix L, TP_matrix U, int start, int end)
{
  int m = S->m;
  int  i, j, k, count;
  long *schur_row_struct;
  long *L_col_ptr = L->col_ptr;
  int  *L_rows    = L->row;
  double *L_vals  = L->val;
  long *U_row_ptr = U->row_ptr;
  int  *U_cols    = U->col;
  double *U_vals  = U->val;
  
  schur_row_struct = malloc((size_t) m * sizeof(*schur_row_struct));
  for ( i = 0; i < m; i++)
    schur_row_struct[i] = -1;

  for( count = 0, i = start; i < end; i++, count++) 
    {
      for(j = L_col_ptr[i]; j < L_col_ptr[i+1]; j++)
	{
	  int    row   = L_rows[j];
	  double L_val = L_vals[j];

	  struct CSR_struct *CSR = &S->CSR[row];
	  int *S_col        = CSR->col;
	  int S_row_nb_elem = CSR->nb_elem;
	  
	  for ( k = 0; k < S_row_nb_elem; k++)
	    schur_row_struct[S_col[k]] = count*m + row ;
	  
	  for ( k = U_row_ptr[i]; k < U_row_ptr[i+1]; k++)
	    if (schur_row_struct[U_cols[k]] == count*m + row) 
	      TP_schur_matrix_add_to_entry(S, row, U_cols[k], -U_vals[k] * L_val);
	    else
	      TP_schur_matrix_insert_entry(S, row, U_cols[k], -U_vals[k] * L_val);
	}
    }
  free(schur_row_struct);

  TP_schur_matrix_update_colmax(S);
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
      printf("================%d======================\n", i);
      int nb_elem = self->CSC[i].nb_elem;
      printf("Colum's max is %f\n", self->CSC[i].col_max);
      for(int j = 0; j < nb_elem; j++)
	{	
	  printf("%d:(%f)  ", self->CSC[i].row[j], self->CSC[i].val[j]);
	}
      printf("\n");
    }

  printf("\n\nPRINTING THE CSR PART\n");
  for(int i = 0; i < m; i++)
    {
      printf("================%d======================\n", i);
      int nb_elem = self->CSR[i].nb_elem;
      for(int j = 0; j < nb_elem; j++)
	{	
	  printf("(%d)  ", self->CSR[i].col[j]);
	}
      printf("\n");
    }
  printf("\n");
}


TP_dense_matrix 
TP_schur_matrix_convert(TP_schur_matrix S, int done_pivots)
{
  TP_dense_matrix self; 
  int col, k, i, j;
  int n = S->n;
  int m = S->m;
  int n_schur = n - done_pivots;
  int m_schur = m - done_pivots;
  self = TP_dense_matrix_create(n_schur, m_schur);

  for(i = 0, k=0; i < m; i++)
    {
      struct CSR_struct *CSR = &S->CSR[i];
      int nb_elem = CSR->nb_elem;
      if ( !nb_elem )
        continue;
      self->original_rows[k++] = i;
    }
  
  for(col = 0, k=0;  col < n; col++)
    {
      struct CSC_struct *CSC = &S->CSC[col];
      int nb_elem = CSC->nb_elem;
      if ( !nb_elem ) 
	continue;
      self->original_cols[k] =  col;
      double *CSC_vals = CSC->val;
      int    *CSC_rows = CSC->row;
      int     found = 0; 
      for( i=0; i < nb_elem; i++) 
	{ 
	  int row = CSC_rows[i];
	  for(j = 0; j < m_schur; j++)
	    if(self->original_rows[j] == row) {
	      found = 1;
	      break;
	    }
	  if ( !found) 
	    TP_fatal_error(__FUNCTION__, __FILE__, __LINE__, "cant find a row index in the ");
	  self->val[k*m_schur + j] =  CSC_vals[i];
	}
      k++;
    }

  return self;
}

void 
TP_schur_matrix_check_perms(TP_schur_matrix self, int *row_perms, 
			    int *col_perms, int nb_pivots)
{
  int n = self->n, m = self->m, i;
  
  for ( i = 0; i < n; i++) 
    if (self->CSC[i].nb_free < 0 ) 
      printf("error on the column %d  with nb_free %d\n", i, self->CSC[i].nb_free);

  for ( i = 0; i < m; i++) 
    if (self->CSR[i].nb_free < 0 ) 
      printf("error on the row %d with nb_free %d\n", i, self->CSR[i].nb_free);

  for( i = 0; i < nb_pivots; i++)
    {
      char *mess;
      int row = row_perms[i], col = col_perms[i];

      if ( self->CSC[col].nb_elem ) {
	asprintf(&mess, "column %d is a pivot, but not empty in S with nb_elem %d and nb_free %d", col, self->CSC[col].nb_elem, self->CSC[col].nb_free);
	TP_warning(__FUNCTION__, __FILE__, __LINE__, mess);
      }      

      if ( self->CSR[row].nb_elem ) {
	asprintf(&mess, "row %d is a pivot, but not empty in S with nb_elem %d and nb_free %d", row, self->CSR[row].nb_elem, self->CSR[row].nb_free);
	TP_warning(__FUNCTION__, __FILE__, __LINE__, mess);
      }      
    }
}
