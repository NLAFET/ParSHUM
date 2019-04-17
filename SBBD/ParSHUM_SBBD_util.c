#include <math.h>
#include <string.h>
#include <limits.h>
#include "ParSHUM_enum.h" 
#include "ParSHUM_auxiliary.h" 
#include "ParSHUM_SBBD_util.h" 


/* TODO:  instead od using proc 0 as master, define a master
          and use it in the rest of the code. This makes it possible to 
	  call ParSHUM SBBD in a subset of procs in a simulation for exmaple.
	  This should be done by adding a root node in the MPI_info struct.  */


void 
ParSHUM_get_col_blocks(ParSHUM_schur_matrix A, col_block col_blocks, row_block row_blocks)
{
  int n = A->n, nb_blocks = row_blocks->nb_blocks, block, i, j, k;
  int *tmp, *BB_index; 
  
  col_blocks->nb_blocks = nb_blocks;
  col_blocks->nb_BB_cols = n;
  col_blocks->n = n;
  col_blocks->perms      = malloc((size_t) n * sizeof(*col_blocks->perms));
  col_blocks->invr_perms = malloc((size_t) n * sizeof(*col_blocks->invr_perms));
  col_blocks->sizes      = calloc((size_t) (nb_blocks + 2 ), sizeof(*col_blocks->sizes));
  col_blocks->nnz        = calloc((size_t) nb_blocks, sizeof(*col_blocks->nnz));
  col_blocks->BB_index   = calloc((size_t) nb_blocks, sizeof(*col_blocks->BB_index));
  col_blocks->BB_size    = calloc((size_t) nb_blocks, sizeof(*col_blocks->BB_size));
  col_blocks->local_schur_m = calloc((size_t) nb_blocks + 1, sizeof(*col_blocks->local_schur_m));
  tmp = malloc((size_t) col_blocks->n * sizeof(*tmp));
  BB_index = malloc((size_t) col_blocks->n * sizeof(*BB_index));

  int_array_memset(col_blocks->invr_perms, ParSHUM_UNUSED_PIVOT, n);
  int_array_memset(tmp, ParSHUM_UNUSED_PIVOT, n);

  for( block = 0; block < nb_blocks; block++)
    {
      int col_block_size = col_blocks->sizes[block];
      int start_block    = row_blocks->sizes[block];
      int end_block      = row_blocks->sizes[block+1];
      
      for( i = start_block; i < end_block; i++)
	{
	  int row = row_blocks->perms[i];
	  CSR_struct *CSR = &A->CSR[row];
	  int *cols       = CSR->col;
	  int row_nb_elem = CSR->nb_elem;
	  col_blocks->nnz[block] += row_nb_elem;

	  for ( j = 0; j < row_nb_elem; j++)
	    {
	      int col = cols[j];
	      if (col_blocks->invr_perms[col] != ParSHUM_UNUSED_PIVOT) {
		if (col_blocks->invr_perms[col] >= col_blocks->nb_BB_cols) {
		  if (tmp[col] > block || tmp[col] < 0) 
		    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"something is wrong");
		  if ( tmp[col] < block) { 
		    BB_index[col_blocks->BB_size[block]++] = col;
		    tmp[col] = block;
		  }
		}
		continue;
	      }
	      CSC_struct *CSC = &A->CSC[col];
	      int *rows       = CSC->row;
	      int col_nb_elem = CSC->nb_elem;
	      int BB_col      = 0;
	      
	      for ( k = 0; k < col_nb_elem; k++) 
		if (row_blocks->invr_perms[rows[k]] < start_block ||
		    row_blocks->invr_perms[rows[k]] >= end_block) {
		  BB_col = 1;
		  break;
		}

	      if (BB_col) {
		col_blocks->perms[--col_blocks->nb_BB_cols] = col;
		col_blocks->invr_perms[col] = col_blocks->nb_BB_cols;
		BB_index[col_blocks->BB_size[block]++] = col;
		tmp[col] = block;
	      } else {
		col_blocks->perms[col_block_size] = col;
		col_blocks->invr_perms[col] = col_block_size++;
	      }
	    }
	}
      col_blocks->sizes[block+1] = col_block_size;
      col_blocks->BB_index[block] = malloc((size_t) col_blocks->BB_size[block] * sizeof(BB_index));
      memcpy((void *) col_blocks->BB_index[block], BB_index, (size_t) col_blocks->BB_size[block] * sizeof(BB_index));
      col_blocks->local_schur_m[block+1] = col_blocks->local_schur_m[block] + 
	(row_blocks->sizes[block+1] - row_blocks->sizes[block]) - 
	(col_blocks->sizes[block+1] - col_blocks->sizes[block]);
    }
  col_blocks->nb_BB_cols = n - col_blocks->nb_BB_cols;
  col_blocks->sizes[nb_blocks+1] = n;

  /* make BB_index local to BB */
  for(block = 0; block < nb_blocks; block++) {   
    int base = n - col_blocks->nb_BB_cols;
    int *_BB_index = col_blocks->BB_index[block];
    for (i = 0; i < col_blocks->BB_size[block]; i++)
      _BB_index[i] = col_blocks->invr_perms[_BB_index[i]] - base;
  }

  for(block = 0; block < nb_blocks; block++)  {
    int *_BB_index = col_blocks->BB_index[block];
    for (i = 0; i < col_blocks->BB_size[block]; i++) {
      int col = _BB_index[i];
      if (col >= col_blocks->nb_BB_cols || col < 0) 
	ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"BB_index is not correct");
    }
  }

  if ( col_blocks->local_schur_m[nb_blocks] != col_blocks->nb_BB_cols) 
    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"local_schur_m is not correct");

  if (col_blocks->sizes[nb_blocks] + col_blocks->nb_BB_cols != n)
    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"col_block sizes are not correct");

  for (i = 0; i < n - col_blocks->nb_BB_cols; i++) 
    if(tmp[col_blocks->perms[i]] != ParSHUM_UNUSED_PIVOT)
      ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"tmp is wrong before");
  
  for (i = n - col_blocks->nb_BB_cols; i < n; i++) 
    if(tmp[col_blocks->perms[i]] == ParSHUM_UNUSED_PIVOT || tmp[col_blocks->perms[i]] >= nb_blocks)
      ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"tmp is wrong before");

  free(tmp);
  free(BB_index);
}


ParSHUM_matrix 
ParSHUM_get_block(ParSHUM_schur_matrix matrix, row_block row_blocks,
		  col_block col_blocks, int block)
{
  int i, j, k, l, start, end, n;

  int nnz = col_blocks->nnz[block], *rows;
  int global_BB = col_blocks->nb_BB_cols, local_BB = col_blocks->BB_size[block];
  int *col_perms = col_blocks->perms;
  int *invr_row_perms = row_blocks->invr_perms;
  int *BB_indices = col_blocks->BB_index[block];
  int start_block = row_blocks->sizes[block];
  int end_block   = row_blocks->sizes[block+1];
  double *vals;
  long *col_ptr; 
  ParSHUM_matrix self = ParSHUM_matrix_create();  

  start = col_blocks->sizes[block];
  end   = col_blocks->sizes[block+1];
  n     = end - start;

  self->n = n + local_BB;
  self->m = row_blocks->sizes[block+1] - row_blocks->sizes[block];
  self->nnz = nnz;
  ParSHUM_matrix_allocate(self, self->n, self->m, self->nnz, 1.0, ParSHUM_CSC_matrix);

  rows    = self->row;
  vals    = self->val;
  col_ptr = self->col_ptr; 
  *col_ptr = 0;
  
  for( i = start, l = 0; i < end; i++) {
    CSC_struct *CSC  = &matrix->CSC[col_perms[i]];
    int *CSC_rows    = CSC->row;
    double *CSC_vals = CSC->val;
    int nb_elem      = CSC->nb_elem;
    int  start_ = (int) col_ptr[l], end_ = (int) col_ptr[l] + nb_elem;

    for( j = start_, k = 0; j < end_; j++, k++) {
      int local_row = invr_row_perms[CSC_rows[k]];
      if (local_row < start_block || local_row >= end_block) 
	ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"independent column has row entry outside the block");
      rows[j] = local_row - start_block;
    }
    memcpy((void *) &vals[start_], (void *) CSC_vals, (size_t) nb_elem*sizeof(*CSC_vals));
    col_ptr[++l] = (long) end_;
  }

  for(j = 0; j < local_BB; j++){
    CSC_struct *CSC  = &matrix->CSC[col_perms[BB_indices[j] + col_blocks->sizes[col_blocks->nb_blocks]]];
    int *CSC_rows    = CSC->row;
    double *CSC_vals = CSC->val;
    int nb_elem      = CSC->nb_elem;
    long col_index = col_ptr[l];

    for(k = 0 ; k < nb_elem; k++) {
      int local_row = invr_row_perms[CSC_rows[k]];
      if (local_row >= start_block && local_row < end_block) {
	if (col_index + 1 > (long) nnz) 
	  ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"block has more entries then expected");
	rows[col_index  ] = local_row - start_block;
	vals[col_index++] = CSC_vals[k];
      }
    }

    col_ptr[++l] = col_index;
  }
  
  /* for( j = start; j < end; j++)  { */
  /*   CSC_struct *CSC  = &matrix->CSC[col_perms[j]]; */
  /*   int *CSC_rows    = CSC->row; */
  /*   double *CSC_vals = CSC->val; */
  /*   int nb_elem      = CSC->nb_elem; */
  /*   int BB_col       = 0; */

  /*   for( k = 0; k < nb_elem; k++) { */
  /*     int local_row = invr_row_perms[CSC_rows[k]]; */
  /*     if (local_row >= start_block && local_row < end_block) { */
  /* 	BB_col = 1; */
  /* 	break; */
  /*     } */
  /*   } */
  /*   if (BB_col) { */
  /*     long col_index = col_ptr[l]; */
  /*     for( ; k < nb_elem; k++) { */
  /* 	int local_row = invr_row_perms[CSC_rows[k]]; */
  /* 	if (local_row >= start_block && local_row < end_block) { */
  /* 	  if (col_index + 1 > (long) nnz) */
  /* 	    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"block has more entries then expected"); */
  /* 	  rows[col_index  ] = local_row - start_block; */
  /* 	  vals[col_index++] = CSC_vals[k]; */
  /* 	} */
  /*     } */
  /*     if (++l > local_BB + n) */
  /* 	ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"block has more columns then expected"); */
  /*     col_ptr[l] = col_index; */
  /*   } */
  /* } */

  if (l != n + local_BB) 
    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"did not found all the columns in the block matrix");
  if(nnz != col_ptr[l]) 
    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"did not found all the enteries in the block matrix");

  return self;
}

void
ParSHUM_print_blocks(row_block row_blocks, col_block col_blocks)
{
  printf("row_block:\nn =  %d \nnb_block %d \n",
	 row_blocks->n, row_blocks->nb_blocks);
  print_int_array(row_blocks->perms,      row_blocks->n, "row permutation");
  print_int_array(row_blocks->invr_perms, row_blocks->n, "row inverse permutation");
  print_int_array(row_blocks->sizes,      row_blocks->nb_blocks + 1, "row block sizes");

  printf("\n\ncol_block:\nn= %d \nnb_blocks= %d\nnb_BB_cols = %d\n",
	 col_blocks->n, col_blocks->nb_blocks, col_blocks->nb_BB_cols);
  print_int_array(col_blocks->perms,      col_blocks->n, "col permutation");
  print_int_array(col_blocks->invr_perms, col_blocks->n, "col inverse permutation");
  print_int_array(col_blocks->nnz,        col_blocks->nb_blocks, "col block non-zero elements");
  print_int_array(col_blocks->sizes,      col_blocks->nb_blocks+2, "col block sizes");
  print_int_array(col_blocks->BB_size,    col_blocks->nb_blocks, "col local border-block sizes");
}

void
ParSHUM_blocks_print_stats(ParSHUM_schur_matrix A, row_block row_blocks, col_block col_blocks)
{
  int i, j, m = A->m, n = A->n ;
  int nb_blocks = col_blocks->nb_blocks;
  int min_n = INT_MAX, max_n = 0, min_nnz = INT_MAX, max_nnz = 0;
  int min_m = INT_MAX, max_m = 0;
  int BB_n, BB_nnz = 0;
  double avg_n, std_n = 0, avg_nnz = 0, std_nnz = 0, total_nnz = 0;
  double avg_m, std_m = 0;
  double *block_nnz;
  block_nnz = malloc(n*sizeof(*block_nnz));
  
  
  avg_m = m / nb_blocks;
  for( i = 0; i < nb_blocks; i++) {
    int block_size = row_blocks->sizes[i+1] - row_blocks->sizes[i];
    std_m += (avg_m - block_size) * (avg_m - block_size);
    max_m = (block_size > max_m) ? block_size : max_m;
    min_m = (block_size < min_m) ? block_size : min_m;
  }
  std_m /= nb_blocks;
  std_m = sqrt(std_m);

  avg_n =  col_blocks->sizes[nb_blocks] / nb_blocks;
  for( i = 0; i < nb_blocks; i++) {
    int block_n = col_blocks->sizes[i+1] - col_blocks->sizes[i];
    block_nnz[i] = 0;
    for( j  = col_blocks->sizes[i]; j < col_blocks->sizes[i+1]; j++)
      block_nnz[i] += A->CSC[col_blocks->perms[j]].nb_elem;
    
    std_n += (avg_n - block_n) * (avg_n - block_n);
    max_n = (block_n > max_n) ? block_n : max_n;
    min_n = (block_n < min_n) ? block_n : min_n;
  }
 
  for( i = 0; i < nb_blocks; i++)  {
    total_nnz += block_nnz[i];
    max_nnz = (block_nnz[i] > max_nnz) ? block_nnz[i] : max_nnz;
    min_nnz = (block_nnz[i] < min_nnz) ? block_nnz[i] : min_nnz;
  }
  avg_nnz = total_nnz / nb_blocks;
  for( i = 0; i < nb_blocks; i++)
    std_nnz += (avg_nnz - block_nnz[i]) * (avg_nnz - block_nnz[i]);

  BB_n = col_blocks->nb_BB_cols;
  for( i = col_blocks->sizes[nb_blocks]; i < col_blocks->sizes[nb_blocks+1]; i++)
    BB_nnz += A->CSC[col_blocks->perms[i]].nb_elem;

  std_nnz /= nb_blocks;
  std_nnz  = sqrt(std_nnz);
  std_n   /= nb_blocks;
  std_n    = sqrt(std_n);
  printf("#blocks\tavg_m\t\tstd_m\t\tmax_m\tmin_m\tavg_n\t\tstd_n\t\tmax_n\tmin_n\tavg_nnz\t\tstd_nnz\t\tmax_nnz\tmin_nnz\tBB_n\tBB_nnz\n");
  printf("%d\t%e\t%e\t%d\t%d\t%e\t%e\t%d\t%d\t%e\t%e\t%d\t%d\t%d\t%d\n", nb_blocks, avg_m,   std_m,   max_m,   min_m,
	                                                                           avg_n,   std_n,   max_n,   min_n,
	                                                                           avg_nnz, std_nnz, max_nnz,min_nnz,
	                                                                           BB_n,    BB_nnz);
}

void
ParSHUM_check_blocks(ParSHUM_schur_matrix A, row_block row_blocks, col_block col_blocks)
{
  int  nb_blocks = row_blocks->nb_blocks, n = col_blocks->n;
  int  m = row_blocks->n, block, i, j;
  char mess[2048];

  /* check_vlaid_perms(col_blocks->perms, col_blocks->invr_perms, n, n, "col"); */
  /* check_vlaid_perms(row_blocks->perms, row_blocks->invr_perms, m, m, "row"); */
  
  if (*row_blocks->sizes)
    ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, "the first row_block size is not zero");
  for( block = 1; block <= nb_blocks; block++)
    if (row_blocks->sizes[block] < row_blocks->sizes[block-1]) {
      snprintf(mess, 2048, "row_block_sizes[%d] = %d is larger then  row_block_sizes[%d] = %d",
  	       block, row_blocks->sizes[block], block - 1, row_blocks->sizes[block - 1]);
      ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
    }
  if (row_blocks->sizes[nb_blocks] != n)
    ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, "the last row_block size is not same as the matrix size");
  
  if (*col_blocks->sizes)
    ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, "the first col_block size is not zero");
  for( block = 1; block <= nb_blocks; block++)
    if (col_blocks->sizes[block] < col_blocks->sizes[block-1]) {
      snprintf(mess, 2048, "col_block_sizes[%d] = %d is larger then  row_col_sizes[%d] = %d",
  	       block, col_blocks->sizes[block], block - 1, col_blocks->sizes[block - 1]);
      ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
    }
  if (col_blocks->sizes[nb_blocks+1] != n)
    ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, "the last col_block size is not same as the matrix size");
  if ((col_blocks->sizes[nb_blocks + 1] - col_blocks->sizes[nb_blocks]) != col_blocks->nb_BB_cols)
    ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__,
  		    "the difference the last two col_block sizes are not same as the nb_BB_cols");

  for (block = 0; block < nb_blocks; block++) {
    int start_block = row_blocks->sizes[block];
    int end_block   = row_blocks->sizes[block+1];

    for (i = col_blocks->sizes[block];  i < col_blocks->sizes[block+1]; i++)
      {
  	int col = col_blocks->perms[i];
  	CSC_struct *CSC = &A->CSC[col];
  	int *rows       = CSC->row;
  	int col_nb_elem = CSC->nb_elem;
  	for ( j = 0; j < col_nb_elem; j++)
  	  if (row_blocks->invr_perms[rows[j]] <  start_block ||
	      row_blocks->invr_perms[rows[j]] >= end_block) {
  	    snprintf(mess, 2048, "col %d belongs to col_block %d, but row %d is not in the same block ",
  		     col, block, rows[j]);
  	    ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
  	  }
      }
  }

  /* Checking the BB cols */
  for (i = col_blocks->sizes[nb_blocks];  i < col_blocks->sizes[nb_blocks + 1]; i++)
    {
      int col = col_blocks->perms[i];
      CSC_struct *CSC    = &A->CSC[col];
      int *rows          = CSC->row;
      int col_nb_elem    = CSC->nb_elem;
      int first_block    = -1;
      int first_invr_row = row_blocks->invr_perms[*rows];
      int is_BB_col = 0, start_block, end_block;

      for ( block = 1; block <= nb_blocks; block++)
	if ( first_invr_row < row_blocks->sizes[block]) {
	  first_block = block-1;
	  start_block = row_blocks->sizes[first_block];
	  end_block   = row_blocks->sizes[block];
	  break;
	}
      
      for ( j = 1; j < col_nb_elem; j++)
	if (row_blocks->invr_perms[rows[j]] <  start_block ||
	    row_blocks->invr_perms[rows[j]] >= end_block) {
  	  is_BB_col = 1;
  	  break;
  	}
      if (!is_BB_col) {
  	snprintf(mess, 2048, "col %d is in the BB block, but all of its entries belongs to %d",
  		 col, first_block);
  	ParSHUM_warning(__FUNCTION__, __FILE__, __LINE__, mess);
      }
    }
}

void
ParSHUM_collect_BB_block(double *local_schur, double *global_schur, col_block col_blocks,
			 row_block row_blocks,int m, int n, int BB_cols, ParSHUM_MPI_info MPI_info) 
{
  double *buff = NULL;
  int i;
  int rank = MPI_info->rank;
  MPI_Comm comm = MPI_info->world;

  if (rank) {
    int schur_m = m - n + BB_cols, j;
    buff         = malloc((size_t) BB_cols * schur_m * sizeof(*buff));

    for( i = 0; i < BB_cols; i++)
      memcpy((void *) &buff[i*schur_m], (void *) &local_schur[i*m + m - schur_m], (size_t) schur_m * sizeof(*buff));
    
    /* printf("%d sending %d elements\n", rank, BB_cols*schur_m); */
    /* for ( i = 0; i < BB_cols; i++) { */
    /*   for (j = 0; j < schur_m; j++) {  */
    /* 	printf("%f\t", buff[i*schur_m + j]); */
    /*   } */
    /*   printf("\n"); */
    /* } */

    MPI_Send(buff, BB_cols*schur_m, MPI_DOUBLE, 0, 0, comm);
  } else {
    int BB_global   = col_blocks->nb_BB_cols;
    int nb_blocks   = col_blocks->nb_blocks;
    int *BB_indices = *col_blocks->BB_index;
    int *local_schur_m = col_blocks->local_schur_m;
    int max_nnz = 0, local_n, local_m, local_nnz;
    int block, current_row = 0;
    MPI_Status status;
    
    for (block = 1; block < nb_blocks; block++) {
      local_n = col_blocks->BB_size[block];
      local_m = local_schur_m[block+1] - local_schur_m[block];
      local_nnz = local_n * local_m;
      
      max_nnz = max_nnz > local_nnz ? max_nnz : local_nnz;
    }
    buff = malloc((size_t) max_nnz * sizeof(*buff));
    
    local_m = m - n + BB_cols;
    /* first habdle my block */
    for ( i = 0; i < BB_cols; i++)
      memcpy((void *) &global_schur[BB_indices[i]*BB_global], (void *) &local_schur[i*m + m - local_m], (size_t) local_m * sizeof(*global_schur));
    current_row += local_m;

    for (block = 1; block < nb_blocks; block++) {
      local_n = col_blocks->BB_size[block];
      local_m = local_schur_m[block+1] - local_schur_m[block]; 
      BB_indices = col_blocks->BB_index[block];

      /* printf("recieving fro %d  %d elements\n", block, local_n * local_m); */
      MPI_Recv(buff, local_n * local_m, MPI_DOUBLE, block, 0, comm, &status);

      for ( i = 0; i < local_n; i++)
	memcpy((void *) &global_schur[BB_indices[i]*BB_global + current_row], (void *) &buff[i*local_m],
	       (size_t) local_m * sizeof(*global_schur));
      current_row += local_m;
    }
    if (current_row != BB_global) 
      printf("KO with BB_global = %d and current_row %d\n", BB_global, current_row);
  }

  free(buff);
}

