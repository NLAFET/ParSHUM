#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "ParSHUM_auxiliary.h"
#include "ParSHUM_Zoltan.h"
#include "ParSHUM_schur_matrix.h"

static int ParSHUM_get_number_of_vertices(void *data, int *ierr);
static void ParSHUM_get_vertex_list(void *data, int sizeGID, int sizeLID,
				    ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
				    int wgt_dim, float *obj_wgts, int *ierr);
static void  ParSHUM_get_hypergraph_size(void *data, int *num_lists, int *num_nonzeroes,
					 int *format, int *ierr);
static void  ParSHUM_get_hypergraph(void *data, int sizeGID, int num_edges, int num_nonzeroes,
				    int format, ZOLTAN_ID_PTR edgeGID, int *vtxPtr,
				    ZOLTAN_ID_PTR vtxGID, int *ierr);

typedef struct _row_block {
  int nb_blocks;
  int n;
  
  int *perms;
  int *invr_perms;
  int *sizes;
} *row_block;

typedef struct _col_block {
  int nb_blocks;
  int n;
  int nb_BB_cols;
  
  int *perms;
  int *invr_perms;
  int *sizes;
  int *nnz;
  int *BB_size;
} *col_block;

struct _Zoltan_Hypergraph  {
  struct Zoltan_Struct *Zoltan;  /* Zoltan's internal structure */

  int rank;
  int MPI_size;
  MPI_Comm world;
  
  int numVertices;            /* number of all vertices */ 
  int numMyVertices;          /* number of vertices that I own */
  ZOLTAN_ID_TYPE *vtxGID;     /* global ID of these vertices */
  
  int numHEdges;              /* number of hyperedges */
  int numMyHEdges;            /* number of my hyperedges */
  int numAllNbors;            /* number of vertices in my hyperedges */
  ZOLTAN_ID_TYPE *edgeGID;    /* global ID of each of my hyperedges */
  int *nborIndex;             /* index into nborGID array of edge's vertices */
  ZOLTAN_ID_TYPE *nborGID;    /* Vertices of edge edgeGID[i] begin at nborGID[nborIndex[i]] */
  
  /* Data needed for the call of Zoltan_LB_Partition */
  int changes, numGidEntries, numLidEntries, numImport, numExport;
  ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
  int *importProcs, *importToPart, *exportProcs, *exportToPart;
};


float
ParSHUM_Zoltan_init(void)
{
  float version;
  int ret = Zoltan_Initialize(0, NULL, &version);

  if (ret != ZOLTAN_OK)
    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"Zoltan failed to initialize");

  return version;
}

Zoltan_Hypergraph 
ParSHUM_Zoltan_create(MPI_Comm comm)
{
  Zoltan_Hypergraph self = calloc((size_t) 1, sizeof(*self));     
  self->world = comm;
  MPI_Comm_rank(self->world, &self->rank);
  MPI_Comm_size(self->world, &self->MPI_size);

  self->Zoltan = Zoltan_Create(comm);
  
  /* General parameters */
  Zoltan_Set_Param(self->Zoltan, "DEBUG_LEVEL", "0");
  Zoltan_Set_Param(self->Zoltan, "LB_METHOD", "HYPERGRAPH");   /* partitioning method */
  Zoltan_Set_Param(self->Zoltan, "HYPERGRAPH_PACKAGE", "PHG"); /* version of method */
  Zoltan_Set_Param(self->Zoltan, "NUM_GID_ENTRIES", "1");      /* global IDs are integers */
  Zoltan_Set_Param(self->Zoltan, "NUM_LID_ENTRIES", "1");      /* local IDs are integers */
  Zoltan_Set_Param(self->Zoltan, "RETURN_LISTS", "ALL");       /* export AND import lists */
  Zoltan_Set_Param(self->Zoltan, "OBJ_WEIGHT_DIM", "0");       /* use Zoltan default vertex weights */
  Zoltan_Set_Param(self->Zoltan, "EDGE_WEIGHT_DIM", "0");      /* use Zoltan default hyperedge weights */
  Zoltan_Set_Param(self->Zoltan, "RETURN_LISTS", "PART");      /* controls what is returned in the  output
								  after the LB_partition call */
  Zoltan_Set_Param(self->Zoltan, "LB_APPROACH", "PARTITION");    /* approach used (PARTITION, REPARTITION, REFINE) */

  return self;
}


int
ParSHUM_Zoltan_init_distrubtion(Zoltan_Hypergraph self, 
				ParSHUM_schur_matrix matrix)
{
  MPI_Comm comm = self->world;
  int rank = self->rank, MPI_size = self->MPI_size;
  int i, j, k;
  MPI_Status status;
  
  /* fields from send_metaData: 
     0:m , 1:my_m, 2:first_m, 3:n, 4:my_n, 5:first_n, 6:my_nnz 
  */
  if (rank == 0) {
    int send_metaData[MPI_size][7];
    int n = matrix->n, m = matrix->m;
    int block_size = (m + MPI_size - 1) / MPI_size;
    int start, end;
    
    for(i = 0, start = 0; i < MPI_size; i++)
      {
	int this_block_size;
	send_metaData[i][0] = m;
	send_metaData[i][2] = start;

	if (start + block_size <= m) 
	  this_block_size = block_size;
	else 
	  this_block_size = m - start;
	send_metaData[i][1] = this_block_size;
	start += this_block_size;
      }

    block_size = (n + MPI_size - 1) / MPI_size;
    for(i = 0, start = 0; i < MPI_size; i++)
      {
	int this_block_size;
	int block_nnz = 0;
	send_metaData[i][3] = n;
	send_metaData[i][5] = start;
	
	if (start + block_size <= n) 
	  this_block_size = block_size;
	else 
	  this_block_size = n - start;
       send_metaData[i][4] = this_block_size;
       end = start + this_block_size;
       for( j = start; j < end; j++)
	 block_nnz += matrix->CSC[j].nb_elem;
       send_metaData[i][6] = block_nnz;
       
       start += this_block_size;
      }
    
    self->numVertices   = send_metaData[0][0];
    self->numMyVertices = send_metaData[0][1];

    self->numHEdges     = send_metaData[0][3];
    self->numMyHEdges   = send_metaData[0][4];
    self->numAllNbors   = send_metaData[0][6];

    if (self->numMyVertices) {
      self->vtxGID = malloc((size_t) self->numMyVertices * sizeof(*self->vtxGID));
      for( i = 0; i < self->numMyVertices; i++)
	self->vtxGID[i] = i;
    }

    if ( self->numMyHEdges) {
      self->edgeGID   = malloc((size_t) self->numMyHEdges * sizeof(*self->edgeGID));
      self->nborIndex = malloc((size_t) (self->numMyHEdges + 1) * sizeof(*self->nborIndex));
      
      for( i = 0; i < self->numMyHEdges; i++)
	self->edgeGID[i] = i;

      *self->nborIndex = 0;
      for( i = 0; i < self->numMyHEdges; i++)
	self->nborIndex[i+1] = matrix->CSC[i].nb_elem + self->nborIndex[i];

      if ( self->numAllNbors) {
	int nb = 0;
	self->nborGID = malloc((size_t) self->numAllNbors * sizeof(*self->nborGID));
	for( i = 0; i < self->numMyHEdges; i++) { 
	  int size = matrix->CSC[i].nb_elem;
	  memcpy((void *) &self->nborGID[nb], (void *) matrix->CSC[i].row, size * sizeof(*self->nborGID));
	  nb += size;
	}
      }
    }
    
    for(i = 1; i < MPI_size; i++) { 
      ZOLTAN_ID_TYPE *zoltan_buff;
      int *int_buff, nb=0;

      MPI_Send(send_metaData[i], 7, MPI_INT, i, 0, comm);

      if (send_metaData[i][4] > 0) {
	int_buff  = malloc((size_t) (send_metaData[i][4] + 1) * sizeof(*int_buff)) ;
	*int_buff = 0;
	for(j = 0, k = send_metaData[i][5]; j < send_metaData[i][4] ; j++, k++)
	  int_buff[j+1] = int_buff[j] + matrix->CSC[k].nb_elem;

	MPI_Send(int_buff, send_metaData[i][4] + 1, MPI_INT, i, 0, comm); 
	free(int_buff);

	if (send_metaData[i][6] > 0) {
	  zoltan_buff = malloc((size_t) send_metaData[i][6] * sizeof(*zoltan_buff));
	  
	  for(j = 0, k = send_metaData[i][5]; j < send_metaData[i][4]; j++, k++) {
	    int size = matrix->CSC[k].nb_elem;
	    memcpy((void *) &zoltan_buff[nb], (void *) matrix->CSC[k].row, size * sizeof(*zoltan_buff));
	    nb+=size;
	  }
	
	  MPI_Send(zoltan_buff, send_metaData[i][6], ZOLTAN_ID_MPI_TYPE, i, 0, comm); 
	  free(zoltan_buff);
	}
      }
    }
  } else {
    int metaData[7];
    
    MPI_Recv(metaData, 7, MPI_INT, 0, 0, comm, &status);
    
    self->numVertices   = metaData[0];
    self->numMyVertices = metaData[1];

    self->numHEdges     = metaData[3];
    self->numMyHEdges   = metaData[4];
    self->numAllNbors   = metaData[6];
    
    if ( self->numMyVertices > 0 ) {
      int start = metaData[2];
      self->vtxGID = malloc((size_t) self->numMyVertices * sizeof(*self->vtxGID));
      for( i = 0; i < self->numMyVertices; i++)
	self->vtxGID[i] = start + i;
    }
    
    if ( self->numMyHEdges > 0 ) {
      int start = metaData[5];
      self->edgeGID = malloc((size_t) self->numMyHEdges * sizeof(*self->edgeGID));
      self->nborIndex = malloc((size_t) (self->numMyHEdges + 1) * sizeof(*self->nborIndex));
      
      for( i = 0; i < self->numMyHEdges; i++)
	self->edgeGID[i] = start + i;

      MPI_Recv(self->nborIndex, self->numMyHEdges + 1, MPI_INT, 0, 0, comm, &status);

      if ( self->numAllNbors > 0 ) {
	self->nborGID = malloc((size_t) self->numAllNbors * sizeof(*self->nborGID));
	MPI_Recv(self->nborGID, self->numAllNbors, ZOLTAN_ID_MPI_TYPE, 0, 0, comm, &status);
      }
    }
  }
  
  return 0;
}

static row_block
ParSHUM_Zoltan_get_row_blocks(Zoltan_Hypergraph hypergraph)
{
  row_block self = calloc(1, sizeof(*self));
  int *parts, *belongs, i;
  int m = hypergraph->numHEdges;
  int rank = hypergraph->rank, MPI_size = hypergraph->MPI_size;
  MPI_Comm comm = hypergraph->world;

  self->nb_blocks = MPI_size;
  self->n         = m;
  if (!rank) {
    self->perms        = malloc((size_t) m * sizeof(*self->perms));
    self->invr_perms   = malloc((size_t) m * sizeof(*self->invr_perms));
    self->sizes        = calloc((size_t) (self->nb_blocks + 2), sizeof(*self->sizes));
    belongs            = malloc((size_t) m * sizeof(*belongs));
    int_array_memset(self->perms, ParSHUM_UNUSED_PIVOT, m);
  }
  parts = malloc((size_t) m * sizeof(*parts));
  int_array_memset(parts, ParSHUM_UNUSED_PIVOT, m);
  
  for( i = 0; i < hypergraph->numExport; i++)
    parts[hypergraph->exportGlobalGids[i]] = hypergraph->exportToPart[i];  

  MPI_Reduce(parts, belongs, m, MPI_INT, MPI_MAX, 0, comm);
  free(parts);

  if (!rank) {
    for( i = 0; i < self->n; i++)
      self->sizes[belongs[i] + 2]++;

    for( i = 1; i <= self->nb_blocks; i++)
      self->sizes[i+1] += self->sizes[i];

    for( i = 0; i < self->n; i++)
      self->perms[self->sizes[belongs[i] + 1]++] = i;

    for( i = 0; i < self->n; i++)
      self->invr_perms[self->perms[i]] = i;
    free(belongs);
  }
  
  return self;
}



static col_block
ParSHUM_Zoltan_get_col_blocks(Zoltan_Hypergraph hypergraph, 
			      ParSHUM_schur_matrix S, 
			      row_block row_blocks)
{
  int n = S->n, nb_blocks = row_blocks->nb_blocks, block, i, j, k;
  col_block col_blocks = calloc(1, sizeof(*col_blocks));
  MPI_Comm comm = hypergraph->world;
  int *tmp; 
  
  col_blocks->nb_blocks = nb_blocks;
  col_blocks->nb_BB_cols = n;
  col_blocks->n = n;
  col_blocks->perms      = malloc((size_t) col_blocks->n * sizeof(*col_blocks->perms));
  col_blocks->invr_perms = malloc((size_t) col_blocks->n * sizeof(*col_blocks->invr_perms));
  col_blocks->sizes      = calloc((size_t) ( col_blocks->nb_blocks + 2 ), sizeof(*col_blocks->sizes));
  col_blocks->nnz        = calloc((size_t) col_blocks->nb_blocks, sizeof(*col_blocks->nnz));
  col_blocks->BB_size    = calloc((size_t) col_blocks->nb_blocks, sizeof(*col_blocks->BB_size));
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
	  CSR_struct *CSR = &S->CSR[row];
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
	      CSC_struct *CSC = &S->CSC[col];
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
      if (block) {
	int local_sizes[4];
	local_sizes[0] = col_blocks->sizes[block+1] - col_blocks->sizes[block];
	local_sizes[1] = col_blocks->BB_size[block];
	local_sizes[2] = row_blocks->sizes[block+1] - row_blocks->sizes[block];
	local_sizes[3] = col_blocks->nnz[block];
	MPI_Send(local_sizes, 4, MPI_INT, block, 0, comm);
      }
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

  return col_blocks;
}

ParSHUM_matrix 
ParSUM_get_block(ParSHUM_schur_matrix matrix, row_block row_blocks,
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
ParSUM_Zoltan_distribute(Zoltan_Hypergraph hypergraph, 
			 ParSHUM_schur_matrix matrix,
			 row_block row_blocks,
			 col_block col_blocks)
{
  int rank = hypergraph->rank, MPI_size = hypergraph->MPI_size;
  MPI_Comm comm = hypergraph->world;
  int n, nb_BB, m, nnz;
  ParSHUM_matrix A;

  if (rank == 0) {
    int block, nb_blocks = col_blocks->nb_blocks, BB_cols = col_blocks->nb_BB_cols;
    

    A = ParSUM_get_block(matrix, row_blocks, col_blocks, 0);
    
    for (block = 1; block < nb_blocks; block++) {
      ParSHUM_matrix block_matrix;  
      
      block_matrix = ParSUM_get_block(matrix, row_blocks, col_blocks, block);
      ParSHUM_matrix_destroy(block_matrix);
    }
    ParSHUM_matrix_destroy(A);
  } else {
    /* int local_sizes[4]; */
    /* MPI_Status status; */
    /* MPI_Recv(local_sizes, 4, MPI_INT, 0, 0, comm, &status); */
    /* n     = local_sizes[0]; */
    /* nb_BB = local_sizes[1]; */
    /* m     = local_sizes[2]; */
    /* nnz   = local_sizes[3]; */
    /* A->n   = n + nb_BB; */
    /* A->m   = m; */
    /* A->nnz = nnz; */
    /* ParSHUM_matrix_allocate(A, A->n, A->m, A->nnz, 1.0, ParSHUM_CSC_matrix); */
  }
}

void
ParSHUM_Zoltan_print_stats(Zoltan_Hypergraph self,
			   ParSHUM_schur_matrix S, 
			   row_block row_blocks,
			   col_block col_blocks)
{
  int i, j, m = self->numVertices, n = self->numHEdges; 
  int nb_blocks = col_blocks->nb_blocks;
  int min_n = INT_MAX, max_n = 0, min_nnz = INT_MAX, max_nnz = 0;
  int min_m = INT_MAX, max_m = 0;
  int BB_n, BB_nnz = 0;
  double avg_n, std_n = 0, avg_nnz = 0, std_nnz = 0, total_nnz = 0, block_nnz[n];
  double avg_m, std_m = 0;

  int  MPI_size = self->MPI_size;

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
      block_nnz[i] += S->CSC[col_blocks->perms[j]].nb_elem;
    
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

  /* BB_n = col_blocks->sizes[nb_blocks+1] - col_blocks->sizes[nb_blocks];  */
  BB_n = col_blocks->nb_BB_cols;
  for( i = col_blocks->sizes[nb_blocks]; i < col_blocks->sizes[nb_blocks+1]; i++) 
    BB_nnz += S->CSC[col_blocks->perms[i]].nb_elem;

  std_nnz /= nb_blocks;
  std_nnz  = sqrt(std_nnz);
  std_n   /= nb_blocks;
  std_n    = sqrt(std_n);
  if (MPI_size == 1)
    printf("#blocks\tavg_m\t\tstd_m\t\tmax_m\tmin_m\tavg_n\t\tstd_n\t\tmax_n\tmin_n\tavg_nnz\t\tstd_nnz\t\tmax_nnz\tmin_nnz\tBB_n\tBB_nnz\n");
  printf("%d\t%e\t%e\t%d\t%d\t%e\t%e\t%d\t%d\t%e\t%e\t%d\t%d\t%d\t%d\n", MPI_size, avg_m,   std_m,   max_m,   min_n,
	                                                                           avg_n,   std_n,   max_n,   min_n,
	                                                                           avg_nnz, std_nnz, max_nnz,min_nnz,
	                                                                           BB_n,    BB_nnz);
}

static void 
ParSHUM_Zoltan_check_blocks(ParSHUM_schur_matrix S, row_block row_blocks, col_block col_blocks)
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
  	CSC_struct *CSC = &S->CSC[col];
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
      CSC_struct *CSC    = &S->CSC[col];
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
ParSHUM_Zoltan_partition(Zoltan_Hypergraph self, ParSHUM_schur_matrix A)
{
  int ret;
  row_block row_blocks;
  col_block col_blocks;
  int  rank = self->rank;

   /* Application defined query functions */
  Zoltan_Set_Num_Obj_Fn(self->Zoltan,    ParSHUM_get_number_of_vertices, self);
  Zoltan_Set_Obj_List_Fn(self->Zoltan,   ParSHUM_get_vertex_list,        self);
  Zoltan_Set_HG_Size_CS_Fn(self->Zoltan, ParSHUM_get_hypergraph_size,    self);
  Zoltan_Set_HG_CS_Fn(self->Zoltan,      ParSHUM_get_hypergraph,         self);
  
  ret = Zoltan_LB_Partition(self->Zoltan,             /* input (all remaining fields are output) */
			    &self->changes,           /* 1 if partitioning was changed, 0 otherwise */ 
			    &self->numGidEntries,     /* Number of integers used for a global ID */
			    &self->numLidEntries,     /* Number of integers used for a local ID */
			    &self->numImport,         /* Number of vertices to be sent to me */
			    &self->importGlobalGids,  /* Global IDs of vertices to be sent to me */
			    &self->importLocalGids,   /* Local IDs of vertices to be sent to me */
			    &self->importProcs,       /* Process rank for source of each incoming vertex */
			    &self->importToPart,      /* New partition for each incoming vertex */
			    &self->numExport,         /* Number of vertices I must send to other processes*/
			    &self->exportGlobalGids,  /* Global IDs of the vertices I must send */
			    &self->exportLocalGids,   /* Local IDs of the vertices I must send */
			    &self->exportProcs,       /* Process to which I send each of the vertices */
			    &self->exportToPart);     /* Partition to which each vertex will belong */
  

  if (ret != ZOLTAN_OK)
    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"Zoltan failed the paritioning");
  
  row_blocks = ParSHUM_Zoltan_get_row_blocks(self);
  if (!rank) {
    col_blocks = ParSHUM_Zoltan_get_col_blocks(self, A, row_blocks);
    ParSHUM_Zoltan_check_blocks(A, row_blocks, col_blocks);
    ParSHUM_Zoltan_print_stats(self, A, row_blocks, col_blocks);
  }
  ParSUM_Zoltan_distribute(self, A, row_blocks, col_blocks);
}

void
ParSHUM_Zoltan_print_distribution(Zoltan_Hypergraph self)
{
  int  rank = self->rank, MPI_size = self->MPI_size, i;
  MPI_Comm comm = self->world;

  for( i = 0; i < MPI_size; i++)
    {
      if ( i == rank) {
	char mess[2048];
	snprintf(mess, 2048, "on process %d, there are %d rows:", i, self->numMyVertices);
 	print_int_array((int *) self->vtxGID, self->numMyVertices, mess);
	snprintf(mess, 2048, "There are %d cols:", self->numMyHEdges);
 	print_int_array((int *) self->edgeGID, self->numMyHEdges, mess);
	snprintf(mess, 2048, "and their col_ptr is:");
 	print_int_array((int *) self->nborIndex, self->numMyHEdges + 1, mess);
	snprintf(mess, 2048, "Having %d entries:",  self->numAllNbors);
 	print_int_array((int *) self->nborGID, self->numAllNbors, mess);
	fprintf(stdout, "************************************************\n\n");
	fflush(stdout);
      }
      MPI_Barrier(comm);
    }
}


void
ParSHUM_Zoltan_destroy(Zoltan_Hypergraph self)
{
  Zoltan_LB_Free_Part(&self->importGlobalGids, &self->importLocalGids,
		      &self->importProcs, &self->importToPart);
  Zoltan_LB_Free_Part(&self->exportGlobalGids, &self->exportLocalGids,
		      &self->exportProcs, &self->exportToPart);

  Zoltan_Destroy(&self->Zoltan);

  if (self->numMyVertices)
    free(self->vtxGID);
  if (self->numMyHEdges) {
    free(self->edgeGID);
    free(self->nborIndex);
    if (self->numAllNbors)
      free(self->nborGID);
  }

  free(self);
}

static int 
ParSHUM_get_number_of_vertices(void *data, int *ierr)
{
  Zoltan_Hypergraph hg = (Zoltan_Hypergraph) data;
  *ierr = ZOLTAN_OK;
  return hg->numMyVertices;
}

static void 
ParSHUM_get_vertex_list(void *data, int sizeGID, int sizeLID,
			ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
			int wgt_dim, float *obj_wgts, int *ierr)
{
  int i;

  Zoltan_Hypergraph hg= (Zoltan_Hypergraph) data;
  *ierr = ZOLTAN_OK;

  for (i=0; i<hg->numMyVertices; i++){
    globalID[i] = hg->vtxGID[i];
    localID[i] = i;
  }
}

static void 
ParSHUM_get_hypergraph_size(void *data, int *num_lists, int *num_nonzeroes,
			    int *format, int *ierr)
{
  Zoltan_Hypergraph hg = (Zoltan_Hypergraph) data;
  *ierr = ZOLTAN_OK;

  *num_lists = hg->numMyHEdges;
  *num_nonzeroes = hg->numAllNbors;

  *format = ZOLTAN_COMPRESSED_EDGE;
}

static void 
ParSHUM_get_hypergraph(void *data, int sizeGID, int num_edges, int num_nonzeroes,
		       int format, ZOLTAN_ID_PTR edgeGID, int *vtxPtr,
		       ZOLTAN_ID_PTR vtxGID, int *ierr)
{
  int i;

  Zoltan_Hypergraph hg = (Zoltan_Hypergraph) data;
  *ierr = ZOLTAN_OK;

  if ( (num_edges != hg->numMyHEdges) || (num_nonzeroes != hg->numAllNbors) ||
       (format != ZOLTAN_COMPRESSED_EDGE)) {
    *ierr = ZOLTAN_FATAL;
    return;
  }

  for (i=0; i < num_edges; i++){
    edgeGID[i] = hg->edgeGID[i];
    vtxPtr[i] = hg->nborIndex[i];
  }

  for (i=0; i < num_nonzeroes; i++){
    vtxGID[i] = hg->nborGID[i];
  }
}
