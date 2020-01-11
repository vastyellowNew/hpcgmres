void APPLY_GIVENS_ROTATION(REAL *column_vec, REAL *givens_c_vec, 
                           REAL *givens_s_vec, int applied_row) {
  REAL tmp1 = givens_c_vec[applied_row] * column_vec[applied_row]
              - givens_s_vec[applied_row] * column_vec[applied_row+1];
  REAL tmp2 = givens_s_vec[applied_row] * column_vec[applied_row]
              + givens_c_vec[applied_row] * column_vec[applied_row+1];
  column_vec[applied_row]   = tmp1;
  column_vec[applied_row+1] = tmp2;
}