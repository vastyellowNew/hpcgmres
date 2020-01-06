void APPLY_GIVENS_ROTATION(REAL *column_vec, REAL * givens_c_vec, 
                           REAL * givens_s_vec, int applied_index) {
  REAL tmp1 = givens_c_vec[applied_index] * column_vec[applied_index] 
                - givens_s_vec[applied_index] * column_vec[applied_index+1];
  REAL tmp2 = givens_s_vec[applied_index] * column_vec[applied_index] 
                + givens_c_vec[applied_index] * column_vec[applied_index+1];
  column_vec[applied_index] = tmp1;
  column_vec[applied_index+1] = tmp2;
}