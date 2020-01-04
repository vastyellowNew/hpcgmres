#ifndef HPCGMRES_D_GIVENS_ROTATION_H_
#define HPCGMRES_D_GIVENS_ROTATION_H_


#include <cmath.h>
#include <blasfeo.h>


void d_givens_rotation_mat(struct blasfeo_dmat *mat, 
                           struct blasfeo_dvec givens_c_vec, 
                           struct blasfeo_dvec givens_s_vec, int column_index, 
                           int row_index);

void d_givens_rotation_vec(struct blasfeo_dvec *vec, 
                           struct blasfeo_dvec givens_c_vec, 
                           struct blasfeo_dvec givens_s_vec, int row_index);


#endif // HPCGMRES_D_GIVENS_ROTATION_H_