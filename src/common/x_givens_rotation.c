void MFGMRES_GIVENS_ROTATION_MAT(struct STRMAT *mat, 
                                 struct STRVEC *givens_c_vec, 
                                 struct STRVEC *givens_s_vec, int column_index, 
                                 int row_index) {
  REAL tmp1 = VECEX1(givens_c_vec, row_index) * MATEX1(mat, column_index, 
                                                       row_index) 
              - VECEX1(givens_s_vec, row_index) 
                  * MATEX1(mat, column_index, row_index+1);
  REAL tmp2 = VECEX1(givens_s_vec, row_index) 
                  * MATEX1(mat, column_index, row_index) 
              + VECEX1(givens_c_vec, row_index) 
                  * MATEX1(mat, column_index, row_index+1);
  MATIN1(tmp1, mat, column_index, row_index);
  MATIN1(tmp2, mat, column_index, row_index+1);
}


void MFGMRES_GIVENS_ROTATION_VEC(struct STRVEC *vec, 
                                 struct STRVEC *givens_c_vec, 
                                 struct STRVEC *givens_s_vec, int row_index) {
  REAL tmp1 = VECEX1(givens_c_vec, row_index) * VECEX1(vec, row_index) 
              - VECEX1(givens_s_vec, row_index) * VECEX1(vec, row_index+1);
  REAL tmp2 = VECEX1(givens_s_vec, row_index) * VECEX1(vec, row_index) 
              + VECEX1(givens_c_vec, row_index) * VECEX1(vec, row_index+1);
  VECIN1(tmp1, vec, row_index);
  VECIN1(tmp2, vec, row_index+1);
}