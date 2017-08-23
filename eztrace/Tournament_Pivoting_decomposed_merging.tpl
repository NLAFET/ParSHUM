BEGIN_MODULE
NAME Tournament_Pivoting
DESC "blabla"
ID 37
TP_pivot_list  get_possible_pivots(TP_solver solver, TP_schur_matrix matrix, int *random_col, TP_pivot_candidates candidates, int nb_threads, double value_tol, double marko_tol, int nb_init_blocks)
TP_pivot_cell  get_independent_cells(TP_pivot_cell cell, int *cols_count, int *rows_count, int base)
TP_pivot_cell  TP_cells_merge(TP_pivot_cell first, TP_pivot_cell second)
int            update_both(TP_pivot_cell cell, TP_schur_matrix matrix, int *cols_count, int *rows_count, int base)
void TP_schur_matrix_update_LD(TP_schur_matrix S, TP_matrix L, TP_matrix D, int *row_perm, int *col_perm, int nb_pivots)
void TP_schur_matrix_update_U(TP_schur_matrix S, TP_U_matrix U, TP_matrix L, int nb_pivots, int *row_perm, TP_U_struct *U_struct, int U_new_n, int U_new_nnz)
void TP_schur_matrix_update_S(TP_schur_matrix S, TP_matrix L, TP_U_matrix U, TP_U_struct *U_struct, int U_new_n, int *invr_row_perm)
END_MODULE