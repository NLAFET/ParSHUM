#include <zoltan.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "ParSHUM_solver.h" 
#include "ParSHUM_matrix.h" 
#include "ParSHUM_dense.h"
#include "ParSHUM_enum.h"
#include "ParSHUM_pivot_list.h" 
#include "ParSHUM_auxiliary.h"

static char *global_fname="hypergraph.txt";

/* Structure to hold distributed hypergraph */

typedef struct _HGRAPH_DATA{
  int numMyVertices;  /* number of vertices that I own initially */
  ZOLTAN_ID_TYPE *vtxGID;        /* global ID of these vertices */

  int numMyHEdges;    /* number of my hyperedges */
  int numAllNbors; /* number of vertices in my hyperedges */
  ZOLTAN_ID_TYPE *edgeGID;       /* global ID of each of my hyperedges */
  int *nborIndex;     /* index into nborGID array of edge's vertices */
  ZOLTAN_ID_TYPE *nborGID;  /* Vertices of edge edgeGID[i] begin at nborGID[nborIndex[i]] */
} HGRAPH_DATA;


typedef struct _block_perms {
  int *perms;
  int *invr_perms;
  int *sizes;
  int n;
  int nb_blocks;
} *block_perms;
block_perms get_row_blocks(ZOLTAN_ID_TYPE *GIDs, int n, int rank,
			   int *init_block, int init_block_size, int nb_parts);
block_perms get_col_blocks(ParSHUM_schur_matrix A, block_perms row_perms);
block_perms get_col_blocks2(ParSHUM_schur_matrix A, block_perms row_perms);


/* 4 application defined query functions.  If we were going to define
 * a weight for each hyperedge, we would need to define 2 more query functions.
 */

static int get_number_of_vertices(void *data, int *ierr);
static void get_vertex_list(void *data, int sizeGID, int sizeLID,
			    ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
			    int wgt_dim, float *obj_wgts, int *ierr);
static void get_hypergraph_size(void *data, int *num_lists, int *num_nonzeroes,
                                int *format, int *ierr);
static void get_hypergraph(void *data, int sizeGID, int num_edges, int num_nonzeroes,
                           int format, ZOLTAN_ID_PTR edgeGID, int *vtxPtr,
                           ZOLTAN_ID_PTR vtxGID, int *ierr);

/* Functions to read hypergraph in from file, distribute it, view it, handle errors */

/* static int get_next_line(FILE *fp, char *buf, int bufsize); */
/* static int get_line_ints(char *buf, int bufsize, int *vals); */
/* static void input_file_error(int numProcs, int tag, int startProc); */
/* static void showHypergraph(int myProc, int numProc, int numIDs, ZOLTAN_ID_TYPE *GIDs, int *parts, int n); */
static void read_input_file(int myRank, int numProcs, char *fname, HGRAPH_DATA *data, ParSHUM_matrix A);

static HGRAPH_DATA global_hg; 

int 
main(int argc, char *argv[])
{
  int i, rc;
  float ver;
  struct Zoltan_Struct *zz;
  int changes, numGidEntries, numLidEntries, numImport, numExport;
  int myRank, numProcs;
  ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids;
  int *importProcs, *importToPart, *exportProcs, *exportToPart;
  int *parts;
  HGRAPH_DATA hg;

  ParSHUM_matrix A = ParSHUM_matrix_create();
  char *file_ext = strrchr(argv[1], '.');


  /******************************************************************
  ** Initialize MPI and Zoltan
  ******************************************************************/
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  rc = Zoltan_Initialize(argc, argv, &ver);

  if (rc != ZOLTAN_OK){
    printf("sorry...\n");
    MPI_Finalize();
    exit(0);
  }

  if (!strcmp(file_ext, ".mtl"))
    ParSHUM_read_mtl_file(A, argv[1]);
#ifdef HAVE_SPRAL
  else  if (!strcmp(file_ext, ".rb"))
    ParSHUM_read_rutherford_boeing(A, argv[1]);
#endif
  else
    ParSHUM_fatal_error(__FUNCTION__, __FILE__, __LINE__,"unsuported matrix file");
  
  /******************************************************************
  ** Read hypergraph from input file and distribute it 
  ******************************************************************/
  read_input_file(myRank, numProcs, global_fname, &hg, A);

  /******************************************************************
  ** Create a Zoltan library structure for this instance of load
  ** balancing.  Set the parameters and query functions that will
  ** govern the library's calculation.  See the Zoltan User's
  ** Guide for the definition of these and many other parameters.
  ******************************************************************/
  zz = Zoltan_Create(MPI_COMM_WORLD);

  /* General parameters */
  Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
  Zoltan_Set_Param(zz, "LB_METHOD", "HYPERGRAPH");   /* partitioning method */
  Zoltan_Set_Param(zz, "HYPERGRAPH_PACKAGE", "PHG"); /* version of method */
  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");/* global IDs are integers */
  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");/* local IDs are integers */
  Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL"); /* export AND import lists */
  Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0"); /* use Zoltan default vertex weights */
  Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", "0");/* use Zoltan default hyperedge weights */
  /* NEW PARM */
  Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");/* controls what is returned in the  output 
						  arrays after the LB_partition call */

  /* PHG parameters  - see the Zoltan User's Guide for many more
   *   (The "REPARTITION" approach asks Zoltan to create a partitioning that is
   *    better but is not too far from the current partitioning, rather than partitioning 
   *    from scratch.  It may be faster but of lower quality that LB_APPROACH=PARTITION.)
  */
  Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");

  /* Application defined query functions */
  Zoltan_Set_Num_Obj_Fn(zz, get_number_of_vertices, &hg);
  Zoltan_Set_Obj_List_Fn(zz, get_vertex_list, &hg);
  Zoltan_Set_HG_Size_CS_Fn(zz, get_hypergraph_size, &hg);
  Zoltan_Set_HG_CS_Fn(zz, get_hypergraph, &hg);
  
  /******************************************************************
  ** Zoltan can now partition the vertices of hypergraph.
  ** In this simple example, we assume the number of partitions is
  ** equal to the number of processes.  Process rank 0 will own
  ** partition 0, process rank 1 will own partition 1, and so on.
  ******************************************************************/
  rc = Zoltan_LB_Partition(zz,                 /* input (all remaining fields are output) */
			   &changes,           /* 1 if partitioning was changed, 0 otherwise */ 
			   &numGidEntries,     /* Number of integers used for a global ID */
			   &numLidEntries,     /* Number of integers used for a local ID */
			   &numImport,         /* Number of vertices to be sent to me */
			   &importGlobalGids,  /* Global IDs of vertices to be sent to me */
			   &importLocalGids,   /* Local IDs of vertices to be sent to me */
			   &importProcs,       /* Process rank for source of each incoming vertex */
			   &importToPart,      /* New partition for each incoming vertex */
			   &numExport,         /* Number of vertices I must send to other processes*/
			   &exportGlobalGids,  /* Global IDs of the vertices I must send */
			   &exportLocalGids,   /* Local IDs of the vertices I must send */
			   &exportProcs,       /* Process to which I send each of the vertices */
			   &exportToPart);     /* Partition to which each vertex will belong */
  
  if (rc != ZOLTAN_OK){
    printf("sorry...\n");
    MPI_Finalize();
    Zoltan_Destroy(&zz);
    exit(0);
  }


  for( i = 0; i < numProcs; i++)  {
    if ( i == myRank) {
      char mess[2048];
      snprintf(mess, 2048, "on process %d exportGlobalGids", i);
      print_int_array((int *) exportGlobalGids, numExport, mess);
      snprintf(mess, 2048, "on process %d exportLocalGids", i);
      print_int_array((int *) exportLocalGids, numExport, mess);
      snprintf(mess, 2048, "on process %d exportProcs", i);
      print_int_array((int *) exportProcs, numExport, mess);
      snprintf(mess, 2048, "on process %d exportToPart", i);
      print_int_array((int *) exportToPart, numExport, mess);
      printf("************************************************\n\n");
      snprintf(mess, 2048, "on process %d importGlobalGids", i);
      print_int_array((int *) importGlobalGids, numImport, mess);
      snprintf(mess, 2048, "on process %d importLocalGids", i);
      print_int_array((int *) importLocalGids, numImport, mess);
      snprintf(mess, 2048, "on process %d importProcs", i);
      print_int_array((int *) importProcs, numImport, mess);
      snprintf(mess, 2048, "on process %d importToPart", i);
      print_int_array((int *) importToPart, numImport, mess);
      printf("************************************************\n");
      printf("************************************************\n");
      printf("************************************************\n\n");
      fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  /******************************************************************
  ** Visualize the hypergraph partitioning before and after calling Zoltan.
  ******************************************************************/
  parts = (int *)malloc(sizeof(int) * hg.numMyVertices);
  int_array_memset(parts, ParSHUM_UNUSED_PIVOT, hg.numMyVertices);

  for( i = 0; i < numProcs; i++)  {
    if ( i == myRank) {
      char mess[2048];
      snprintf(mess, 2048, "on process %d parts after calloc", i);
      print_int_array((int *) parts, numExport, mess);
      fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  for (i=0; i < hg.numMyVertices; i++){
    parts[i] = myRank;
  }
  for( i = 0; i < numProcs; i++)  {
    if ( i == myRank) {
      char mess[2048];
      snprintf(mess, 2048, "on process %d parts in the begining", i);
      print_int_array((int *) parts, numExport, mess);
      fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  
  if (myRank== 0){
    printf("\nHypergraph partition before calling Zoltan\n");
  }

  /* showHypergraph(myRank, numProcs, hg.numMyVertices, hg.vtxGID, parts, A->n); */

  for (i=0; i < numExport; i++){
    parts[exportLocalGids[i]] = exportToPart[i];
  }

  for( i = 0; i < numProcs; i++)  {
    if ( i == myRank) {
      char mess[2048];
      snprintf(mess, 2048, "on process %d parts", i);
      print_int_array((int *) parts, numExport, mess);
      fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  
  if (myRank == 0){
    printf("Graph partition after calling Zoltan\n");
  }

  /* showHypergraph(myRank, numProcs, hg.numMyVertices, hg.vtxGID, parts, A->n); */

  block_perms  row_perms = get_row_blocks(hg.vtxGID, A->n, myRank, parts, numExport, numProcs);
  if (myRank == 0){

    /* ParSHUM_schur_matrix S = ParSHUM_schur_matrix_create(); */
    /* ParSHUM_schur_matrix_allocate(S, A->n, A->n, A->nnz, 0, 0, 0, 0.0, 0.0); */
    /* ParSHUM_schur_matrix_copy(A, S, 0.2); */
    ParSHUM_matrix_print_as_dense(A, "not permuted");

    int *tmp = malloc((size_t) A->n * sizeof(*tmp));
    int o;
    for(o = 0; o < A->n; o++)
      tmp[o] = o;
    ParSHUM_matrix C = ParSHUM_matrix_permute(A, tmp, row_perms->invr_perms);
    ParSHUM_matrix_print_as_dense(C, "row permuted");

    ParSHUM_schur_matrix S = ParSHUM_schur_matrix_create();
    ParSHUM_schur_matrix_allocate(S, C->n, C->n, C->nnz, 0, 0, 0, 0.0, 0.0);
    ParSHUM_schur_matrix_copy(C, S, 0.2);

    block_perms  col_perms = get_col_blocks2(S, row_perms);

    print_int_array(row_perms->perms, A->n, "row perms");
    print_int_array(col_perms->perms, A->n, "col perms");
    print_int_array(col_perms->invr_perms, A->n, "invr_col perms");

    ParSHUM_matrix B = ParSHUM_matrix_permute(A, col_perms->perms, row_perms->invr_perms);
    ParSHUM_matrix_print_as_dense(B, "all permuted");
  }

  /******************************************************************
  ** Free the arrays allocated by Zoltan_LB_Partition, and free
  ** the storage allocated for the Zoltan structure.
  ******************************************************************/
  Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids,
                      &importProcs, &importToPart);
  Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids,
                      &exportProcs, &exportToPart);

  Zoltan_Destroy(&zz);

  /**********************
  ** all done ***********
  **********************/
  MPI_Finalize();

  if (hg.numMyVertices > 0){
    free(parts);
    free(hg.vtxGID);
  }
  if (hg.numMyHEdges > 0){
    free(hg.edgeGID);
    free(hg.nborIndex);
    if (hg.numAllNbors > 0){
      free(hg.nborGID);
    }
  }

  ParSHUM_matrix_destroy(A);

  return 0;
}

/* Invr_perms from the block_perms is not the inverse permutation, but for each elemnt 
   it tells you to which partition that row belongs. */
block_perms
get_row_blocks(ZOLTAN_ID_TYPE *GIDs, int n, int rank,
	       int *init_block, int init_block_size, int nb_parts)
{
  int i;
  block_perms perms = calloc(1, sizeof(*perms));
  int *my_row_block = calloc((size_t) n, sizeof(*my_row_block));
  int *all_row_block;
  if( rank == 0 ) {
    perms->n = n;
    perms->nb_blocks = nb_parts;
    perms->perms = calloc((size_t) n, sizeof(*perms->perms));
    perms->invr_perms = calloc((size_t) n, sizeof(*perms->invr_perms));
    perms->sizes = calloc((size_t) nb_parts + 1, sizeof(*perms->sizes));
    all_row_block = malloc((size_t) n * sizeof(*all_row_block));
    int_array_memset(all_row_block, ParSHUM_UNUSED_PIVOT, n);
  }
  
  int_array_memset(my_row_block, ParSHUM_UNUSED_PIVOT, n);

  for(i = 0; i < init_block_size; i++)
    my_row_block[GIDs[i]] = init_block[i];

  MPI_Reduce(my_row_block, all_row_block, n, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

  free(my_row_block);
  if (rank > 0) {
    return NULL;
  }
  print_int_array((int *) all_row_block, n, "after MPI_Reduce");

  int *sizes = perms->sizes;
  int *row_perms = perms->perms;
  int *invr_row_perms = perms->invr_perms;

  for(i = 0; i < n; i++)
    sizes[all_row_block[i]+1]++;
  print_int_array((int *) sizes, nb_parts+1, "row block sizes first");

  for(i = 1; i <= nb_parts; i++)
    sizes[i] += sizes[i-1];
  print_int_array((int *) sizes, nb_parts+1, "row block sizes");

  for(i = 0; i < n; i++) 
    row_perms[sizes[all_row_block[i]]++] = i;
  print_int_array((int *) sizes, nb_parts+1, "row block sizes after");

  printf("nb parts = %d \n", nb_parts);
  for(i = nb_parts - 1; i > 0; i--)
    sizes[i] -= sizes[i] - sizes[i-1];
  *sizes = 0;
  print_int_array((int *) sizes, nb_parts+1, "row block sizes after all");

  for(i = 0; i < n; i++) 
    invr_row_perms[row_perms[i]] = i;

  return perms;
}

int
get_block(int row, int nb_blocks, int *blocks)
{
  int my_block, found = 1;
		
  for ( my_block = 0; my_block < nb_blocks; my_block++)
    if ( row >= blocks[my_block] &&  row < blocks[my_block + 1])  {
      found = 1;
      break;
    }
  
  if (found)
    return my_block;
  else
    return -1;
}

block_perms
get_col_blocks2(ParSHUM_schur_matrix A, block_perms row_perms)
{
  int n = A->n, BB_cols = A->n, affiliated_cols = 0, row, i, j;
  int nb_blocks = row_perms->nb_blocks;
  int *row_block_sizes = row_perms->sizes;
  block_perms col_perm = calloc(1, sizeof(*col_perm));

  col_perm->n = n;
  col_perm->nb_blocks  = nb_blocks;
  col_perm->perms      = calloc((size_t) n, sizeof(*col_perm->perms));
  col_perm->invr_perms = malloc((size_t) n * sizeof(*col_perm->invr_perms));
  int *col_perms = col_perm->perms;
  int *invr_col_perms = col_perm->invr_perms;

  int_array_memset(invr_col_perms, ParSHUM_UNUSED_PIVOT, n);
  for( row = 0; row < n; row++)
    {
      int my_block = get_block(row, nb_blocks, row_block_sizes);
      printf("my blocks is %d\n", my_block);
      CSR_struct *CSR = &A->CSR[row];
      int *cols       = CSR->col;
      int row_nb_elem     = CSR->nb_elem;
      for ( i = 0; i < row_nb_elem; i++)
	{
	  int col = cols[i];
	  if (invr_col_perms[col] != ParSHUM_UNUSED_PIVOT) 
	    continue;
	  CSC_struct *CSC = &A->CSC[col];
	  int *rows       = CSC->row;
	  int col_nb_elem     = CSC->nb_elem;
	  int BB_col = 0;
	  for ( j = 0; j < col_nb_elem; j++) {
	    int res = get_block(rows[j], nb_blocks, row_block_sizes);
	    if (res < 0)
	      printf("BIG KO\n");
	    else if ( my_block !=  res)  {
	      BB_col = 1;
	      break;
	    }
	  }
	  if (BB_col) {
	    col_perms[--BB_cols] = col;
	    invr_col_perms[col] = BB_cols;
	  } else {
	    col_perms[affiliated_cols] = col;
	    invr_col_perms[col] = affiliated_cols++;
	  }
	}
    }
  printf("there are %d BB cols, %d cols in the ii blocks and n = %d\n", n - BB_cols, affiliated_cols, n);

  return col_perm;
}


block_perms
get_col_blocks(ParSHUM_schur_matrix A, block_perms row_perms)
{
  int n = A->n, BB_cols = A->n, affiliated_cols = 0, i, j, k;
  int nb_blocks = row_perms->nb_blocks;
  int *row_perm = row_perms->perms;
  int *row_block_affiliation = row_perms->invr_perms;
  block_perms col_perm = calloc(1, sizeof(*col_perm));

  col_perm->n = n;
  col_perm->nb_blocks  = nb_blocks;
  col_perm->perms      = calloc((size_t) n, sizeof(*col_perm->perms));
  col_perm->invr_perms = malloc((size_t) n * sizeof(*col_perm->invr_perms));
  /* col_perms->sizes      = calloc((size_t) nb_blocks + 1, sizeof(*col_perms->sizes)); */
  int *col_perms = col_perm->perms;
  int *invr_col_perms = col_perm->invr_perms;
  /* int *col_sizes = col_perms->sizes; */
  int_array_memset(invr_col_perms, ParSHUM_UNUSED_PIVOT, n);

  for( i = 0; i < n; i++)
    {
      int my_row      = row_perm[i];
      int my_block    = row_block_affiliation[i];
      CSR_struct *CSR = &A->CSR[my_row];
      int *cols       = CSR->col;
      int row_nb_elem     = CSR->nb_elem;

      for ( j = 0; j < row_nb_elem; j++)
	{
	  int col = cols[j];
	  if (invr_col_perms[col] != ParSHUM_UNUSED_PIVOT) 
	    continue;
	  CSC_struct *CSC = &A->CSC[col];
	  int *rows = CSC->row;
	  int col_nb_elem = CSC->nb_elem;
	  int BB_col = 0;
	  for ( k = 0; k < col_nb_elem; k++)
	    {
	      int row = rows[k];
	      if (row_block_affiliation[row] != my_block)  {
		BB_col = 1;
		break;
	      }
	    }
	  if (BB_col) {
	    col_perms[--BB_cols] = col;
	    invr_col_perms[col] = BB_cols;
	  } else {
	    col_perms[affiliated_cols] = col;
	    invr_col_perms[col] = affiliated_cols++;
	  }
	}
    }
  printf("BB cols %d and affiliated cols %d\n", n - BB_cols, affiliated_cols);
  return  col_perm;
}


/* Application defined query functions */
static int get_number_of_vertices(void *data, int *ierr)
{
  HGRAPH_DATA *hg = (HGRAPH_DATA *)data;
  *ierr = ZOLTAN_OK;
  return hg->numMyVertices;
}

static void get_vertex_list(void *data, int sizeGID, int sizeLID,
			    ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
			    int wgt_dim, float *obj_wgts, int *ierr)
{
int i;

  HGRAPH_DATA *hg= (HGRAPH_DATA *)data;
  *ierr = ZOLTAN_OK;

  /* In this example, return the IDs of our vertices, but no weights.
   * Zoltan will assume equally weighted vertices.
   */

  for (i=0; i<hg->numMyVertices; i++){
    globalID[i] = hg->vtxGID[i];
    localID[i] = i;
  }
}

static void get_hypergraph_size(void *data, int *num_lists, int *num_nonzeroes,
                                int *format, int *ierr)
{
  HGRAPH_DATA *hg = (HGRAPH_DATA *)data;
  *ierr = ZOLTAN_OK;

  *num_lists = hg->numMyHEdges;
  *num_nonzeroes = hg->numAllNbors;

  /* We will provide compressed hyperedge (row) format.  The alternative is
   * is compressed vertex (column) format: ZOLTAN_COMPRESSED_VERTEX.
   */

  *format = ZOLTAN_COMPRESSED_EDGE;

  return;
}

static void get_hypergraph(void *data, int sizeGID, int num_edges, int num_nonzeroes,
                           int format, ZOLTAN_ID_PTR edgeGID, int *vtxPtr,
                           ZOLTAN_ID_PTR vtxGID, int *ierr)
{
int i;

  HGRAPH_DATA *hg = (HGRAPH_DATA *)data;
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

  return;
}


void read_input_file(int myRank, int numProcs,
		     char *fname, HGRAPH_DATA *hg,
		     ParSHUM_matrix A)
{
int numGlobalVertices, numGlobalEdges, numGlobalNZ;
int num, count, nnbors, ack=0;
int to=-1, from, remaining;
int i, j;
int send_count[3];
ZOLTAN_ID_TYPE *idx;
unsigned int id;
MPI_Status status;
int ack_tag = 5, count_tag = 10, id_tag = 15;
HGRAPH_DATA *send_hg;

  if (myRank == 0){

    numGlobalVertices = A->n;
    global_hg.numMyVertices = numGlobalVertices = numGlobalVertices;
    global_hg.vtxGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * numGlobalVertices);

    for (i=0; i < numGlobalVertices; i++){

      global_hg.vtxGID[i] = (ZOLTAN_ID_TYPE) i;
    }
    numGlobalEdges = A->n;
    global_hg.numMyHEdges = numGlobalEdges;
    global_hg.edgeGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * numGlobalEdges);
    global_hg.nborIndex = (int *)malloc(sizeof(int) * (numGlobalEdges + 1));

    numGlobalNZ = A->nnz;
    global_hg.nborGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * numGlobalNZ);

    /* Get the list of vertices in each hyperedge  */
    global_hg.nborIndex[0] = 0;

    for (i=0; i < numGlobalEdges; i++){
      int row_end  =  A->col_ptr[i+1];
      int row_start =  A->col_ptr[i];

      id = i;
      nnbors = row_end - row_start;

      global_hg.edgeGID[i] = (ZOLTAN_ID_TYPE)id;

      for (j=0; j < nnbors; j++){
        global_hg.nborGID[global_hg.nborIndex[i] + j] = (ZOLTAN_ID_TYPE)A->row[row_start + j];
      }

      global_hg.nborIndex[i+1] = global_hg.nborIndex[i] + nnbors;
    }


    /* Create a sub graph for each process */
    send_hg = (HGRAPH_DATA *)calloc(sizeof(HGRAPH_DATA) , numProcs);

    /* 
     * Divide the vertices across the processes
     */
    remaining = numGlobalVertices;
    count = (numGlobalVertices + numProcs  -1) / numProcs;
    idx = global_hg.vtxGID;

    for (i=0; i < numProcs; i++){

      if (remaining == 0) count = 0;
      if (count > remaining) count = remaining;

      send_hg[i].numMyVertices = count;

      if (count){

        send_hg[i].vtxGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * count);

        for (j=0; j < count; j++){
          send_hg[i].vtxGID[j] = *idx++;
        }
      }

      remaining -= count;
    }

    /*
     * Assign hyperedges to processes, and create a sub-hypergraph for each process.
     */
    remaining = numGlobalEdges;
    count = (numGlobalVertices + numProcs  -1) / numProcs;
    from = 0;

    for (i=0; i < numProcs; i++){

      if (remaining == 0) count = 0;
      if (count > remaining) count = remaining;

      send_hg[i].numMyHEdges = count;
      send_hg[i].numAllNbors = 0;

      if (count > 0){
  
        to = from + count;
  
        nnbors = global_hg.nborIndex[to] - global_hg.nborIndex[from];

        send_hg[i].numAllNbors = nnbors;
  
        send_hg[i].edgeGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * count);
        memcpy(send_hg[i].edgeGID, global_hg.edgeGID + from, sizeof(ZOLTAN_ID_TYPE) * count);

        send_hg[i].nborIndex = (int *)malloc(sizeof(int) * (count + 1));
        send_hg[i].nborIndex[0] = 0;

        if (nnbors > 0){

          num = global_hg.nborIndex[from];

          for (j=1; j <= count; j++){
            send_hg[i].nborIndex[j] = global_hg.nborIndex[from+j] - num;
          }
        
          send_hg[i].nborGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * nnbors);
          memcpy(send_hg[i].nborGID, 
                 global_hg.nborGID + global_hg.nborIndex[from], 
                 sizeof(ZOLTAN_ID_TYPE) * nnbors);
        }
      }

      remaining -= count;
      from = to;
    }

    /* Send each process its hyperedges and the vertices in its partition */
    *hg = send_hg[0];

    for (i=1; i < numProcs; i++){
      send_count[0] = send_hg[i].numMyVertices;
      send_count[1] = send_hg[i].numMyHEdges;
      send_count[2] = send_hg[i].numAllNbors;

      MPI_Send(send_count, 3, MPI_INT, i, count_tag, MPI_COMM_WORLD);
      MPI_Recv(&ack, 1, MPI_INT, i, ack_tag, MPI_COMM_WORLD, &status);

      if (send_count[0] > 0){
        MPI_Send(send_hg[i].vtxGID, send_count[0], ZOLTAN_ID_MPI_TYPE, i, id_tag, MPI_COMM_WORLD);
        free(send_hg[i].vtxGID);
      }

      if (send_count[1] > 0){
        MPI_Send(send_hg[i].edgeGID, send_count[1], ZOLTAN_ID_MPI_TYPE, i, id_tag + 1, MPI_COMM_WORLD);
        free(send_hg[i].edgeGID);

        MPI_Send(send_hg[i].nborIndex, send_count[1] + 1, MPI_INT, i, id_tag + 2, MPI_COMM_WORLD);
        free(send_hg[i].nborIndex);

        if (send_count[2] > 0){
          MPI_Send(send_hg[i].nborGID, send_count[2], ZOLTAN_ID_MPI_TYPE, i, id_tag + 3, MPI_COMM_WORLD);
          free(send_hg[i].nborGID);
        }
      }
    }

    free(send_hg);

    /* signal all procs it is OK to go on */
    ack = 0;
    for (i=1; i < numProcs; i++){
      MPI_Send(&ack, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    }
  }
  else{

    MPI_Recv(send_count, 3, MPI_INT, 0, count_tag, MPI_COMM_WORLD, &status);

    if (send_count[0] < 0){
      MPI_Finalize();
      exit(1);
    }

    ack = 0;

    memset(hg, 0, sizeof(HGRAPH_DATA));

    hg->numMyVertices = send_count[0];
    hg->numMyHEdges = send_count[1];
    hg->numAllNbors = send_count[2];

    if (send_count[0] > 0){
      hg->vtxGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * send_count[0]);
    }

    if (send_count[1] > 0){
      hg->edgeGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * send_count[1]);
      hg->nborIndex = (int *)malloc(sizeof(int) * (send_count[1] + 1));

      if (send_count[2] > 0){
        hg->nborGID = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * send_count[2]);
      }
    }

    MPI_Send(&ack, 1, MPI_INT, 0, ack_tag, MPI_COMM_WORLD);

    if (send_count[0] > 0){
      MPI_Recv(hg->vtxGID,send_count[0], ZOLTAN_ID_MPI_TYPE, 0, id_tag, MPI_COMM_WORLD, &status);

      if (send_count[1] > 0){
        MPI_Recv(hg->edgeGID,send_count[1], ZOLTAN_ID_MPI_TYPE, 0, id_tag + 1, MPI_COMM_WORLD, &status);
        MPI_Recv(hg->nborIndex,send_count[1] + 1, MPI_INT, 0, id_tag + 2, MPI_COMM_WORLD, &status);

        if (send_count[2] > 0){
          MPI_Recv(hg->nborGID,send_count[2], ZOLTAN_ID_MPI_TYPE, 0, id_tag + 3, MPI_COMM_WORLD, &status);
        }
      }
    }

    /* ok to go on? */
    MPI_Recv(&ack, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    if (ack < 0){
      MPI_Finalize();
      exit(1);
    }
  }

}
