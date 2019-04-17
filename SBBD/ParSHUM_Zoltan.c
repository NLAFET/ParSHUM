#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "ParSHUM_auxiliary.h"
#include "ParSHUM_Zoltan.h"
#include "ParSHUM_SBBD_util.h"
#include "ParSHUM_schur_matrix.h"

#define  BLOCK_SIZES   4
#define  BLOCK_COLS    0
#define  BLOCK_BB_COLS 1
#define  BLOCK_ROWS    2
#define  BLOCK_NNZ     3

static int
ParSHUM_Zoltan_get_number_of_vertices(void *data, int *ierr);

static void
ParSHUM_Zoltan_get_vertex_list(void *data, int sizeGID, int sizeLID,
			       ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
			       int wgt_dim, float *obj_wgts, int *ierr);

static void  
ParSHUM_Zoltan_get_hypergraph_size(void *data, int *num_lists, int *num_nonzeroes,
				   int *format, int *ierr);

static void  
ParSHUM_Zoltan_get_hypergraph(void *data, int sizeGID, int num_edges, int num_nonzeroes,
			      int format, ZOLTAN_ID_PTR edgeGID, int *vtxPtr,
			      ZOLTAN_ID_PTR vtxGID, int *ierr);


struct _Zoltan_Hypergraph  {
  struct Zoltan_Struct *Zoltan;  /* Zoltan's internal structure */
  ParSHUM_MPI_info MPI_info;
  float version;

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


Zoltan_Hypergraph 
ParSHUM_Zoltan_create(ParSHUM_MPI_info MPI_info)
{
  Zoltan_Hypergraph self = calloc((size_t) 1, sizeof(*self));     
  
  if (Zoltan_Initialize(0, NULL, &self->version) != ZOLTAN_OK)
    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"Zoltan failed to initialize");
  self->MPI_info = MPI_info;

  self->Zoltan = Zoltan_Create(self->MPI_info->world);

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
ParSHUM_Zoltan_register_data(Zoltan_Hypergraph self, ParSHUM_schur_matrix matrix)
{
  MPI_Comm comm = self->MPI_info->world;
  int rank = self->MPI_info->rank, MPI_size = self->MPI_info->MPI_size;
  int i, j, k;
  MPI_Status status;

  /* fields for send_metaData: 
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

void
ParSHUM_Zoltan_partition(Zoltan_Hypergraph self, ParSHUM_schur_matrix A)
{
  int ret;

   /* Application defined query functions */
  Zoltan_Set_Num_Obj_Fn(self->Zoltan,    ParSHUM_Zoltan_get_number_of_vertices, self);
  Zoltan_Set_Obj_List_Fn(self->Zoltan,   ParSHUM_Zoltan_get_vertex_list,        self);
  Zoltan_Set_HG_Size_CS_Fn(self->Zoltan, ParSHUM_Zoltan_get_hypergraph_size,    self);
  Zoltan_Set_HG_CS_Fn(self->Zoltan,      ParSHUM_Zoltan_get_hypergraph,         self);
  
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
}

ParSHUM_matrix
ParSUM_Zoltan_distribute(ParSHUM_schur_matrix matrix, row_block row_blocks,
			 col_block col_blocks, ParSHUM_solver solver, ParSHUM_MPI_info MPI_info)
{
  int rank = MPI_info->rank;
  MPI_Comm comm = MPI_info->world;
  int n, nb_BB, m, nnz;
  ParSHUM_matrix A;

  if (rank == 0) {
    int block, nb_blocks = col_blocks->nb_blocks;

    A = ParSHUM_get_block(matrix, row_blocks, col_blocks, 0);
    A->nnz = A->col_ptr[A->n];
    solver->BB_cols = *col_blocks->BB_size;
    
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
    A->nnz = A->col_ptr[A->n];
    solver->BB_cols = nb_BB;
  }

  return A;
}

void
ParSHUM_Zoltan_get_row_blocks(Zoltan_Hypergraph hypergraph, row_block row_blocks)
{
  int *parts, *belongs, i;
  int m = hypergraph->numHEdges;
  int rank = hypergraph->MPI_info->rank, MPI_size = hypergraph->MPI_info->MPI_size;
  MPI_Comm comm = hypergraph->MPI_info->world;

  row_blocks->nb_blocks = MPI_size;
  row_blocks->n         = m;
  if (!rank) {
    row_blocks->perms        = malloc((size_t) m * sizeof(*row_blocks->perms));
    row_blocks->invr_perms   = malloc((size_t) m * sizeof(*row_blocks->invr_perms));
    row_blocks->sizes        = calloc((size_t) (row_blocks->nb_blocks + 2), sizeof(*row_blocks->sizes));
    belongs            = malloc((size_t) m * sizeof(*belongs));
    int_array_memset(row_blocks->perms, ParSHUM_UNUSED_PIVOT, m);
  }
  parts = malloc((size_t) m * sizeof(*parts));
  int_array_memset(parts, ParSHUM_UNUSED_PIVOT, m);
  
  for( i = 0; i < hypergraph->numExport; i++)
    parts[hypergraph->exportGlobalGids[i]] = hypergraph->exportToPart[i];  

  MPI_Reduce(parts, belongs, m, MPI_INT, MPI_MAX, 0, comm);
  free(parts);

  if (!rank) {
    for( i = 0; i < row_blocks->n; i++)
      row_blocks->sizes[belongs[i] + 2]++;

    for( i = 1; i <= row_blocks->nb_blocks; i++)
      row_blocks->sizes[i+1] += row_blocks->sizes[i];

    for( i = 0; i < row_blocks->n; i++)
      row_blocks->perms[row_blocks->sizes[belongs[i] + 1]++] = i;

    for( i = 0; i < row_blocks->n; i++)
      row_blocks->invr_perms[row_blocks->perms[i]] = i;
    free(belongs);
  }
}

void
ParSHUM_Zoltan_print_distribution(Zoltan_Hypergraph self)
{
  int  rank = self->MPI_info->rank, MPI_size = self->MPI_info->MPI_size, i;
  MPI_Comm comm = self->MPI_info->world;

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
ParSHUM_Zoltan_get_number_of_vertices(void *data, int *ierr)
{
  Zoltan_Hypergraph hg = (Zoltan_Hypergraph) data;
  *ierr = ZOLTAN_OK;
  return hg->numMyVertices;
}

static void 
ParSHUM_Zoltan_get_vertex_list(void *data, int sizeGID, int sizeLID,
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
ParSHUM_Zoltan_get_hypergraph_size(void *data, int *num_lists, int *num_nonzeroes,
				   int *format, int *ierr)
{
  Zoltan_Hypergraph hg = (Zoltan_Hypergraph) data;
  *ierr = ZOLTAN_OK;

  *num_lists = hg->numMyHEdges;
  *num_nonzeroes = hg->numAllNbors;

  *format = ZOLTAN_COMPRESSED_EDGE;
}

static void 
ParSHUM_Zoltan_get_hypergraph(void *data, int sizeGID, int num_edges, int num_nonzeroes,
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
