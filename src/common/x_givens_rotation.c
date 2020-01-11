void APPLY_GIVENS_ROTATION_TO_MAT(struct STRMAT *mat, 
                                  struct STRVEC *givens_c_vec, 
                                  struct STRVEC *givens_s_vec, 
                                  int applied_column, int applied_row) {
  REAL tmp1 = VECEX1(givens_c_vec, applied_row) 
                  * MATEX1(mat, applied_column, applied_row)
              - VECEX1(givens_s_vec, applied_row) 
                  * MATEX1(mat, applied_column, applied_row+1);
  REAL tmp2 = VECEX1(givens_s_vec, applied_row) 
                  * MATEX1(mat, applied_column, applied_row) 
              + VECEX1(givens_c_vec, applied_row) 
                  * MATEX1(mat, applied_column, applied_row+1);
  MATIN1(tmp1, mat, applied_column, applied_row+1);
  MATIN1(tmp2, mat, applied_column, applied_row);
}


void APPLY_GIVENS_ROTATION_TO_VEC(struct STRVEC *vec, 
                                  struct STRVEC *givens_c_vec, 
                                  struct STRVEC *givens_s_vec, 
                                  int applied_index) {
  REAL tmp1 = VECEX1(givens_c_vec, applied_index) 
                  * VECEX1(vec, applied_index)
              - VECEX1(givens_s_vec, applied_index) 
                  * VECEX1(vec, applied_index+1);
  REAL tmp2 = VECEX1(givens_s_vec, applied_index) 
                  * VECEX1(vec, applied_index) 
              + VECEX1(givens_c_vec, applied_index) 
                  * VECEX1(vec, applied_index+1);
  VECIN1(tmp1, vec, applied_index);
  VECIN1(tmp2, vec, applied_index+1);
}


void APPLY_GIVENS_ROTATION(REAL *column_vec, REAL *givens_c_vec, 
                           REAL *givens_s_vec, int applied_row) {
  REAL tmp1 = givens_c_vec[applied_row] * column_vec[applied_row]
              - givens_s_vec[applied_row] * column_vec[applied_row+1];
  REAL tmp2 = givens_s_vec[applied_row] * column_vec[applied_row]
              + givens_c_vec[applied_row] * column_vec[applied_row+1];
  column_vec[applied_row]   = tmp1;
  column_vec[applied_row+1] = tmp2;
}