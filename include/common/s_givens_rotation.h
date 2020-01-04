#ifndef HPCGMRES_S_GIVENS_ROTATION_H_
#define HPCGMRES_S_GIVENS_ROTATION_H_


#include <cmath.h>
#include <blasfeo.h>


void s_givens_rotation_mat(struct blasfeo_smat *mat, 
                           struct blasfeo_svec givens_c_vec, 
                           struct blasfeo_svec givens_s_vec, int column_index, 
                           int row_index);

void s_givens_rotation_vec(struct blasfeo_svec *vec, 
                           struct blasfeo_svec givens_c_vec, 
                           struct blasfeo_svec givens_s_vec, int row_index);


#endif // HPCGMRES_S_GIVENS_ROTATION_H_