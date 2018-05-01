#ifndef _ParSHUM_PIVOT_LIST_H
#define _ParSHUM_PIVOT_LIST_H

#include "ParSHUM_schur_matrix.h"
#include "ParSHUM_solver.h"

typedef struct _ParSHUM_pivot_list *ParSHUM_pivot_list;
typedef struct _ParSHUM_pivot_set  *ParSHUM_pivot_set;
typedef struct _ParSHUM_pivot_cell *ParSHUM_pivot_cell;

struct _ParSHUM_pivot_list 
{
  ParSHUM_pivot_set sets;
  int nb_elem;
};

struct _ParSHUM_pivot_set 
{
  ParSHUM_pivot_set next;
  ParSHUM_pivot_cell cells;
  int nb_elem;
  int *rows_count;
  int *cols_count;
  int base;
};

struct _ParSHUM_pivot_cell
{
  int row; 
  int col;
  int marko;
  ParSHUM_pivot_cell next;
};

// PIVOT LIST
ParSHUM_pivot_list  ParSHUM_pivot_list_create();
ParSHUM_pivot_set   ParSHUM_pivot_set_create(ParSHUM_solver solver, int n, int m);
ParSHUM_pivot_cell  ParSHUM_pivot_cell_create(int row, int col, int marko);
ParSHUM_pivot_list  ParSHUM_pivot_list_insert_set(ParSHUM_pivot_list self, ParSHUM_pivot_set set);
ParSHUM_pivot_set   get_next_merging_set(ParSHUM_pivot_list self);
ParSHUM_pivot_set   merge_sorted_sets(ParSHUM_pivot_set self);

ParSHUM_pivot_set   merge_to_larger_set(ParSHUM_pivot_set self, ParSHUM_schur_matrix matrix, ParSHUM_solver solver);
ParSHUM_pivot_cell  add_cell_to_sorted_set(ParSHUM_pivot_set set, ParSHUM_pivot_cell cell,
					   ParSHUM_schur_matrix matrix);
void                ParSHUM_pivot_list_destroy(ParSHUM_pivot_list self, ParSHUM_solver solver);
void                ParSHUM_pivot_set_destroy(ParSHUM_pivot_set self, ParSHUM_solver solver);
void                print_pivot_list(ParSHUM_pivot_list self, char *mess);
void                ParSHUM_pivot_cell_destroy(ParSHUM_pivot_cell self);

// OPS
ParSHUM_pivot_list  get_possible_pivots(ParSHUM_solver solver, ParSHUM_schur_matrix matrix, int *random_col,
					ParSHUM_pivot_candidates candidates, int nb_threads,
					double value_tol, double marko_tol, int nb_init_blocks);
ParSHUM_pivot_list  merge_pivot_sets(ParSHUM_pivot_list self, ParSHUM_schur_matrix matrix, ParSHUM_solver solver);


#endif // ParSHUM_PIVOT_LIST_H
