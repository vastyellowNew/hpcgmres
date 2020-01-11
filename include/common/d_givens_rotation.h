#ifndef HPCGMRES_D_GIVENS_ROTATION_H_
#define HPCGMRES_D_GIVENS_ROTATION_H_


#include <blasfeo.h>

#ifdef __cplusplus
extern "C" {
#endif


void d_apply_givens_rotation_to_mat(struct blasfeo_dmat *mat, 
                                    struct blasfeo_dvec *givens_c_vec, 
                                    struct blasfeo_dvec *givens_s_vec, 
                                    int applied_column, int applied_row);

void d_apply_givens_rotation_to_vec(struct blasfeo_dvec *vec, 
                                    struct blasfeo_dvec *givens_c_vec, 
                                    struct blasfeo_dvec *givens_s_vec, 
                                    int applied_row);

void d_apply_givens_rotation(double *column_vec, double *givens_c_vec, 
                             double *givens_s_vec, int applied_row);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // HPCGMRES_D_GIVENS_ROTATION_H_