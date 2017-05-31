#ifndef _TP_PIVOT_LIST_H
#define _TP_PIVOT_LIST_H

#include "TP_schur_matrix.h"

typedef struct _TP_pivot_list *TP_pivot_list;
typedef struct _TP_pivot_set  *TP_pivot_set;
typedef struct _TP_pivot_cell *TP_pivot_cell;



struct _TP_pivot_list 
{
  TP_pivot_set first;
  TP_pivot_set midle;
  TP_pivot_set last;
  int nb_elem;
};

struct _TP_pivot_set 
{
  TP_pivot_set next;
  TP_pivot_cell cells;
  int nb_elem;
  int *rows_count;
  int *cols_count;
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
TP_pivot_set   TP_pivot_set_create(int n, int m);
TP_pivot_cell  TP_pivot_cell_create(int row, int col, int marko);
TP_pivot_list  TP_pivot_list_insert_new_set(TP_pivot_list self, TP_schur_matrix matrix, int row, int col, int marko);
TP_pivot_list  TP_pivot_list_insert_set(TP_pivot_list self, TP_pivot_set set);
TP_pivot_set   get_next_merging_set(TP_pivot_list self);
TP_pivot_set   merge_sorted_sets(TP_pivot_set self);
TP_pivot_set   get_independent_pivots(TP_pivot_set candidates, TP_schur_matrix matrix);
TP_pivot_set   merge_to_larger_set(TP_pivot_set self, TP_schur_matrix matrix);
void           add_cell_to_sorted_set(TP_pivot_set set, TP_pivot_cell cell, TP_schur_matrix matrix);
void           TP_pivot_list_destroy(TP_pivot_list self);
void           TP_pivot_set_destroy(TP_pivot_set self);
void           print_pivot_list(TP_pivot_list self, char *mess);



// OPS
TP_pivot_list  get_possible_pivots(TP_schur_matrix matrix, int *random_col, double value_tol, double marko_tol, int nb_init_blocks);
TP_pivot_list  merge_pivot_sets(TP_pivot_list self, TP_schur_matrix matrix);


#endif // TP_PIVOT_LIST_H
