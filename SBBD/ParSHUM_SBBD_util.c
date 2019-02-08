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
#define  BLOCK_SIZES   4
#define  BLOCK_COLS    0
#define  BLOCK_BB_COLS 1
#define  BLOCK_ROWS    2
#define  BLOCK_NNZ     3

static ParSHUM_matrix
ParSHUM_get_block(ParSHUM_schur_matrix matrix, row_block row_blocks,
		  col_block col_blocks, int block);


void 
ParSHUM_get_col_blocks(ParSHUM_schur_matrix A, col_block col_blocks, row_block row_blocks)
{
  int n = A->n, nb_blocks = row_blocks->nb_blocks, block, i, j, k;
  int *tmp; 
  
  col_blocks->nb_blocks = nb_blocks;
  col_blocks->nb_BB_cols = n;
  col_blocks->n = n;
  col_blocks->perms      = malloc((size_t) n * sizeof(*col_blocks->perms));
  col_blocks->invr_perms = malloc((size_t) n * sizeof(*col_blocks->invr_perms));
  col_blocks->sizes      = calloc((size_t) (nb_blocks + 2 ), sizeof(*col_blocks->sizes));
  col_blocks->nnz        = calloc((size_t) nb_blocks, sizeof(*col_blocks->nnz));
  col_blocks->BB_size    = calloc((size_t) nb_blocks, sizeof(*col_blocks->BB_size));
  tmp = malloc((size_t) col_blocks->n * sizeof(*tmp));

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
		    col_blocks->BB_size[block]++;
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
		col_blocks->BB_size[block]++;
		tmp[col] = block;
	      } else {
		col_blocks->perms[col_block_size] = col;
		col_blocks->invr_perms[col] = col_block_size++;
	      }
	    }
	}
      col_blocks->sizes[block+1] = col_block_size;
    }
  col_blocks->nb_BB_cols = n - col_blocks->nb_BB_cols ;
  col_blocks->sizes[nb_blocks+1] = n;
  
  if (col_blocks->sizes[nb_blocks] + col_blocks->nb_BB_cols != n)
    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"col_block sizes are not correct");

  for (i = 0; i < n - col_blocks->nb_BB_cols; i++) 
    if(tmp[col_blocks->perms[i]] != ParSHUM_UNUSED_PIVOT)
      ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"tmp is wrong before");
  
  for (i = n - col_blocks->nb_BB_cols; i < n; i++) 
    if(tmp[col_blocks->perms[i]] == ParSHUM_UNUSED_PIVOT || tmp[col_blocks->perms[i]] >= nb_blocks)
      ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"tmp is wrong before");

  free(tmp);
}

ParSHUM_matrix
ParSUM_Zoltan_distribute(ParSHUM_schur_matrix matrix, row_block row_blocks,
			 col_block col_blocks, ParSHUM_MPI_info MPI_info)
{
  int rank = MPI_info->rank;
  MPI_Comm comm = MPI_info->world;
  int n, nb_BB, m, nnz;
  ParSHUM_matrix A;

  if (rank == 0) {
    int block, nb_blocks = col_blocks->nb_blocks;

    A = ParSHUM_get_block(matrix, row_blocks, col_blocks, 0);
    
    for (block = 1; block < nb_blocks; block++) {
      ParSHUM_matrix block_matrix;  
      int block_sizes[BLOCK_SIZES];
      block_sizes[BLOCK_COLS]    = col_blocks->sizes[block+1] - col_blocks->sizes[block];
      block_sizes[BLOCK_BB_COLS] = col_blocks->BB_size[block];
      block_sizes[BLOCK_ROWS]    = row_blocks->sizes[block+1] - row_blocks->sizes[block];
      block_sizes[BLOCK_NNZ]     = col_blocks->nnz[block];
      MPI_Send(block_sizes, BLOCK_SIZES, MPI_INT, block, 0, comm);
      
      block_matrix = ParSHUM_get_block(matrix, row_blocks, col_blocks, block);

      MPI_Send(block_matrix->row,       block_matrix->nnz, MPI_INT,    block, 0, comm);
      MPI_Send(block_matrix->val,       block_matrix->nnz, MPI_DOUBLE, block, 0, comm);
      MPI_Send(block_matrix->col_ptr, block_matrix->n + 1, MPI_LONG,   block, 0, comm);

      ParSHUM_matrix_destroy(block_matrix);
    }
  } else {
    int my_sizes[BLOCK_SIZES];
    MPI_Status status;
    MPI_Recv(my_sizes, 4, MPI_INT, 0, 0, comm, &status);
    n      = my_sizes[BLOCK_COLS];
    nb_BB  = my_sizes[BLOCK_BB_COLS];
    m      = my_sizes[BLOCK_ROWS];
    nnz    = my_sizes[BLOCK_NNZ];

    A = ParSHUM_matrix_create();
    A->n   = n + nb_BB;
    A->m   = m;
    A->nnz = nnz;
    ParSHUM_matrix_allocate(A, A->n, A->m, A->nnz, 1.0, ParSHUM_CSC_matrix);
    MPI_Recv(A->row,     A->nnz,   MPI_INT,    0, 0, comm, &status);
    MPI_Recv(A->val,     A->nnz,   MPI_DOUBLE, 0, 0, comm, &status);
    MPI_Recv(A->col_ptr, A->n + 1, MPI_LONG,   0, 0, comm, &status);
  }

  return A;
}

static ParSHUM_matrix 
ParSHUM_get_block(ParSHUM_schur_matrix matrix, row_block row_blocks,
		  col_block col_blocks, int block)
{
  int i, j, k, l, start, end, n;

  int nnz = col_blocks->nnz[block], *rows;
  int global_BB = col_blocks->nb_BB_cols, local_BB = col_blocks->BB_size[block];
  int *col_perms = col_blocks->perms;
  int *invr_row_perms = row_blocks->invr_perms;
  int start_block = row_blocks->sizes[block];
  int end_block   = row_blocks->sizes[block+1];
  double *vals;
  long *col_ptr; 
  ParSHUM_matrix self = ParSHUM_matrix_create();  

  start = col_blocks->sizes[block];
  end   = col_blocks->sizes[block+1];
  n     = end - start;

  self->n = col_blocks->sizes[block+1] - col_blocks->sizes[block] + local_BB;
  self->m = n;
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

  start = col_blocks->n - global_BB;
  end   = col_blocks->n;
  for( j = start; j < end; j++)  {
    CSC_struct *CSC  = &matrix->CSC[col_perms[j]];
    int *CSC_rows    = CSC->row;
    double *CSC_vals = CSC->val;
    int nb_elem      = CSC->nb_elem;
    int BB_col       = 0;

    for( k = 0; k < nb_elem; k++) {
      int local_row = invr_row_perms[CSC_rows[k]];
      if (local_row >= start_block && local_row < end_block) {
	BB_col = 1;
	break;
      }
    }
    if (BB_col) {
      long col_index = col_ptr[l];
      for( ; k < nb_elem; k++) {
	int local_row = invr_row_perms[CSC_rows[k]];
	if (local_row >= start_block && local_row < end_block) {
	  if (col_index + 1 > (long) nnz) 
	    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"block has more entries then expected");
	  rows[col_index  ] = local_row - start_block;
	  vals[col_index++] = CSC_vals[k];
	}
      }
      if (++l > local_BB + n)  
	ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"block has more columns then expected");
      col_ptr[l] = col_index;
    }
  }

  if (l != n + local_BB) 
    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"did not found all the columns in the block matrix");
  if(nnz != col_ptr[l]) 
    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"did not found all the enteries in the block matrix");

  return self;
}

void
ParSHUM_Zoltan_print_stats(ParSHUM_schur_matrix A, row_block row_blocks, col_block col_blocks)
{
  int i, j, m = A->m, n = A->n ;
  int nb_blocks = col_blocks->nb_blocks;
  int min_n = INT_MAX, max_n = 0, min_nnz = INT_MAX, max_nnz = 0;
  int min_m = INT_MAX, max_m = 0;
  int BB_n, BB_nnz = 0;
  double avg_n, std_n = 0, avg_nnz = 0, std_nnz = 0, total_nnz = 0, block_nnz[n];
  double avg_m, std_m = 0;

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
  printf("%d\t%e\t%e\t%d\t%d\t%e\t%e\t%d\t%d\t%e\t%e\t%d\t%d\t%d\t%d\n", nb_blocks, avg_m,   std_m,   max_m,   min_n,
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

  check_vlaid_perms(row_blocks->perms, m, m);
  check_vlaid_perms(col_blocks->perms, n, n);
  check_perms_and_invr_perms(row_blocks->perms, row_blocks->invr_perms, m, "row_block");
  check_perms_and_invr_perms(col_blocks->perms, col_blocks->invr_perms, n, "col_block");
  
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
