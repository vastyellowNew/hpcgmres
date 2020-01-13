#ifndef CGMRES_S_GIVENS_ROTATION_H_
#define CGMRES_S_GIVENS_ROTATION_H_


#ifdef __cplusplus
extern "C" {
#endif

void s_apply_givens_rotation(float *column_vec, float *givens_c_vec, 
                             float *givens_s_vec, int applied_row);

#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // CGMRES_S_GIVENS_ROTATION_H_