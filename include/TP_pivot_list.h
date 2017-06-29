#ifndef _TP_PIVOT_LIST_H
#define _TP_PIVOT_LIST_H

#include "TP_schur_matrix.h"
#include "TP_solver.h"

typedef struct _TP_pivot_list *TP_pivot_list;
typedef struct _TP_pivot_set  *TP_pivot_set;
typedef struct _TP_pivot_cell *TP_pivot_cell;

struct _TP_pivot_list 
{
  TP_pivot_set sets;
  int nb_elem;
};

struct _TP_pivot_set 
{
  TP_pivot_set next;
  TP_pivot_cell cells;
  int nb_elem;
  int *rows_count;
  int *cols_count;
  int base;
};

struct _TP_pivot_cell
{
  int row; 
  int col;
  int marko;
  TP_pivot_cell next;
};

// PIVOT LIST
TP_pivot_list  TP_pivot_list_create();
TP_pivot_set   TP_pivot_set_create(TP_solver solver, int n, int m);
TP_pivot_cell  TP_pivot_cell_create(int row, int col, int marko);
/* TODO:clean up */
/* TP_pivot_list  TP_pivot_list_insert_new_set(TP_pivot_list self, TP_schur_matrix matrix, int row, int col, int marko); */
TP_pivot_list  TP_pivot_list_insert_set(TP_pivot_list self, TP_pivot_set set);
TP_pivot_set   get_next_merging_set(TP_pivot_list self);
TP_pivot_set   merge_sorted_sets(TP_pivot_set self);
/* TODO:clean up */
/* TP_pivot_set   get_independent_pivots(TP_pivot_set candidates, TP_schur_matrix matrix); */
TP_pivot_set   merge_to_larger_set(TP_pivot_set self, TP_schur_matrix matrix, TP_solver solver);
TP_pivot_cell  add_cell_to_sorted_set(TP_pivot_set set, TP_pivot_cell cell, TP_schur_matrix matrix);
void           TP_pivot_list_destroy(TP_pivot_list self, TP_solver solver);
void           TP_pivot_set_destroy(TP_pivot_set self, TP_solver solver);
void           print_pivot_list(TP_pivot_list self, char *mess);
void           TP_pivot_cell_destroy(TP_pivot_cell self);



// OPS
TP_pivot_list  get_possible_pivots(TP_solver solver, TP_schur_matrix matrix, int *random_col,
				   TP_pivot_candidates candidates, int nb_threads,
				   double value_tol, double marko_tol, int nb_init_blocks);
TP_pivot_list  merge_pivot_sets(TP_pivot_list self, TP_schur_matrix matrix, TP_solver solver);


#endif // TP_PIVOT_LIST_H
