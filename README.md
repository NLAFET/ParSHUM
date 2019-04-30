# The ParSHUM library 
The ParSHUM library implements parallel sparse a
LU factorization. The configuration and installation process for the library is explained
next, followed by a description of how to use the driver that is provided by the library and
the interface of the library. The last section explains how to obtain optimal performance
when using this solver.

# Getting sources
ParSHUM is developed at STFC in the context of the NLAFET project, and its source
codes are available on GitHub. The source codes can be obtained with the following git
command:
```console 
$ git clone https://github.com/NLAFET/ParSHUM.git
```
# Configuration
This library relies on cmake for the building of the library. The ParSHUM library depends
only on the MKL library for the BLAS and LAPACK routines. The CFLAGS and LDFLAGS 
for MKL must be provided by passing as arguments −DMKL_LDFLAGS and
−DMKL_CFLAGS to the cmake script. If these two argument are not given, then the
cmake commands fails. These can be obtained from the [MKL Link Advisor](https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor). Additionally, the SPRAL library could be provided in order to be able to read matrices stored in
Rutherford Boeing format. This can be done by providing the install root directory of
SPRAL to the −DSPRAL_DIR cmake’s argument.

# Installing
Once the library has been downloaded into a directory, issue the following commands in
that directory:
```console 
$ mkdir build &&  cd build
$ cmake -DMKL_LDFLAGS=<MKL LDFLAGS> -DMKL_CFLAGS=<MKL CFLAGS> -DSPRAL_DIR=<spral installation directory>
```
To build the library use:
```console 
$ make 
```

# Using the test driver
Once the library has been built a test driver called ParSHUM_simple is created in the bin
directory. This driver takes a large set of parameters:
* --matrix (string) =>     specifies the path to the matrix file. Our driver is able to read classic IJV files
(suffixed by "*.ijv" or "*.mtl"). Harwell-Boeing matrices (suffixed by "*.hb" or "*.rb")
can be read if the library is compiled with the −DSP RAL_DIR option.
* --marko_threshold (real) =>
Markowitz threshold of acceptability for a pivot. This threshold should be larger
than one. The default value is 4.
* --threshold_value (real) =>
threshold value of acceptability for a pivot. This threshold should be between zero
and one. The default value is 10 −4 .
* --schur_density_tolerance (real) =>
maximum density allowed of the Schur complement. Once the Schur complement
reaches this density, we switch to dense factorization. This parameter should be
between zero and one. The default value is 0.2.
* --nb_previous_steps (integer) =>
number of previous steps for which we keep track of how many pivots have been
found. This parameter should be at least one. The default value is 5.
* --min_pivot_per_steps (integer) =>
minimum number of pivots in the last nb_previous_steps steps. If the number
of pivots is less than this value, the solver switches to a dense factorization. This
parameter should be at least nb_previous_steps. The default value is ten times
nb_previous_steps.
* --max_dense_schur (real) =>
the maximum allowed size for the dense factorization. If the Schur complement is
larger than this value, then the dense factorization is not performed and an error is
returned by the solver. The default value is 20,000.
* --nb_threads (integer) =>
the number of threads used. This parameter should be at least 1. The default value
is 1.
* --extra_space (real) =>
the factor of the memory initially allocated for the solver in comparison with the
size of the input matrix. The default value is 3.
* --trace =>
activates a creation of a trace of the execution of the algorithm in Pajé format. By
default this option is deactivated. 
* --verbosity (integer) =>
Controls the verbosity level.
If the value is smaller and equal then 0, no otput is produced. If it is  
1 a general output is printed. If it is larger then  one, 
 extra information about each is step is given. Additionaly, plot, data and
fig directories are created and a data file is writen, prefixed by the parameters of
the algorithm, for each run. A script is created in the plot directory that if compiled
with gnuplot will create a graph of the execution time in the f ig directory. 

# Using the ParSHUM interface
In this subsection we will explain how to use the interface of the ParSHUM library. The
source code for the test driver presented previously is presented next: 

```C
int
main(int argc, char **argv)
{
  ParSHUM_solver solver;
  ParSHUM_vector X, B;
  /* Create the solver */
  solver = ParSHUM_solver_create();   

  /*  Parse the arguments */
  ParSHUM_solver_parse_args(solver, argc, argv);
  /*  Read the matrix */
  ParSHUM_solver_read_matrix(solver);

  /*  Initialize the vectors */
  X = ParSHUM_vector_create(solver->A->n);
  B = ParSHUM_vector_create(solver->A->n);
  ParSHUM_vector_read_file(solver, X);

  /* copy the vector B in X */
  ParSHUM_vector_copy(B, X);
  /*  Initialize the solver */
  ParSHUM_solver_init(solver);

  /*  Perform the factorization */
  ParSHUM_solver_factorize(solver);

  /*  Perform the solve operation */
  ParSHUM_solver_solve(solver, B);

  /*  Compute the norms */
  ParSHUM_solver_compute_norms(solver, X, B);

  /*  Finalize the solver */
  ParSHUM_solver_finalize(solver);

  /*  Free all the data */
  ParSHUM_vector_destroy(X);
  ParSHUM_vector_destroy(B);
  ParSHUM_solver_destroy(solver);

  return 0;
}
```
First the solver structure is created by calling ParSHUM_solver_create in line 7. This will
only allocate the basic data for the solver. Then on line 10, the arguments are processed
and the solver’s parameters are set to the values provided. The matrix is then read from
the file followed by the initialisation of the two vectors B and X. Once everything is in
place, the function ParSHUM_solver_init is called which will allocate all the internal
data that is needed for the factorization. The factorization is performed by the function
ParSHUM_solver_factorize on line 25, followed by the solve operation on line 28 and
3the computation of the norms on line 31. The ParSHUM_solver_finalize needs to be
called before the data is freed. This function finalises the PLASMA library, creating the
trace file if demanded and prints the output of the solver. Finally, all the data that has
been used is freed (lines 37-39). 

Alternatively, is possible for the user to control the input
matrix and parameters for the solver directly in the code and we now explain how it can
be done. In order to change the input matrix and parameters, we use the ParSHUM_solver
structure. The only fields that
should be modified by the user are the input matrix A and the exe_parms
field. The ParSHUM_solver is presented next: 
```C
typedef struct _ParSHUM_solver {
  ParSHUM_matrix   A;  /* The input matrix */

  /* The factors */
  ParSHUM_L_matrix L;   
  ParSHUM_matrix   D;
  ParSHUM_U_matrix U;

  /* The Schur complement */
  ParSHUM_schur_matrix S;

  ParSHUM_exe_parms exe_parms;
  ParSHUM_verbose verbose;
  .....
} * ParSHUM_solver; 
```
In the first example the input matrix was loaded by calling the
function ParSHUM_solver_read_matrix. The other option is to set it as the input matrix
in the ParSHUM_solver structure. We present the ParSHUM_matrix structure and a code  example 
on how to assign a CSC matrix as the input matrix:
```C
typedef  enum _ParSHUM_matrix_type
{
  ParSHUM_CSC_matrix,
  ParSHUM_CSR_matrix,
  ParSHUM_Diag_matrix,
  ParSHUM_Rutherford_matrix
} ParSHUM_matrix_type;

typedef struct _ParSHUM_matrix {
  ParSHUM_matrix_type type;
  int n;

  long allocated;
  long nz;

  int *row;
  long *col_ptr;
  double *val;

  int *col;
  long *row_ptr;

  /* This is used by the SPRAL matrix driver since it is a Fortran code (see http://www.test-numerical.rl.ac.uk/spral/doc/sphinx/C/rutherford_boeing.html) */
  void *handle;
} * ParSHUM_matrix;

   .....
   /* An example for assigning a CSC matrix as an input matrix. 
   The following code should replace the line 12 in the previous example. */

   ParSHUM_matrix matrix = calloc(1, sizeof(struct _ParSHUM_matrix);
   matrix->type = ParSHUM_CSC_matrix;
   matrix->n = n;
   matrix->allocated = nz;
   matrix->nz = nz;
   /* If the input matrix is a CSR matrix, the field type should be assigned with ParSHUM_CSR_matrix and instead of assigning arrays to the row and col_ptr fields, the col and row_ptr fields should be assigned. */
   matrix->row = row;
   matrix->col_ptr = col_ptr;
   solver->A = matrix;
   .....
```
In the first example, the solver is parametrized
by calling the function ParSHUM_solver_parse_args. The other option is to
modify directly the exe_parms field of the solver structure. The exe_parms field is of
type ParSHUM_exe_parms and an example of how to use this structure to change the
threshold value and Markowitz threshold. All the other parameters
could be changed in the same way.
```C
typedef struct ParSHUM_exe_parms {
  double value_tol;
  double marko_tol;
  double extra_space;
  double density_tolerance;

  char *matrix_file;
  char *RHS_file;

  int min_pivot_per_steps;
  int nb_threads;
  int nb_previous_pivots;
  int max_dense_schur;
  int trace;
} *ParSHUM_exe_parms;

   .....
   /* A code example for changing the parameters of the solver */ 
   /* The following code should replace the call the ParSHUM_solver_parse_args (line 10). */
   solver->exe_parms->value_threshold = 0.001;
   solver->exe_parms->marko_threshold = 16;
   .....
```
Once the matrix and the parameters are set, the solver should be initialised by calling
the ParSHUM_solver_init. From this point the ParSHUM_solver structure should not
be accessed or modified directly by the user.

# Getting optimal performance with ParSHUM
We strongly recommend that on NUMA machines, ParSHUM should be binded to a single
NUMA node. This algorithm is highly dependent on the bandwidth of the machine and
not on the computational power To do so, one way is to use the numactl command. For
example:
```console
$ numactl --cpunodebind=0 --membind=0 ./bin/ParSHUM_simple <arguments>
```
will bind the process to the cores and the memory of NUMA node 0.
Although internally in ParSHUM we bind the threads through the proc_bind option of
the omp pragma, we have observed that this argument is ignored for some OpenMP imple-
mentations. This could be also achieved through the OMP_PROC_BIN D environment
variable. Additionally, 10% to 20% could be gained by setting the OMP_WAIT_POLICY
environment variable to ACTIVE. Here is an example for an optimal run:
```console
$ export OMP_WAIT_POLICY=ACTIVE
$ export OMP_PROC_BIND=close
$ numactl --cpunodebind=0 --membind=0 ./bin/ParSHUM_simple <arguments>
```
