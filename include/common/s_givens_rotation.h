#ifndef HPCGMRES_S_GIVENS_ROTATION_H_
#define HPCGMRES_S_GIVENS_ROTATION_H_

void s_apply_givens_rotation(float *column_vec, float *givens_c_vec, 
                             float *givens_s_vec, int applied_index);

#endif // HPCGMRES_S_GIVENS_ROTATION_H_