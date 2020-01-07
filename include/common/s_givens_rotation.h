#ifndef HPCGMRES_S_GIVENS_ROTATION_H_
#define HPCGMRES_S_GIVENS_ROTATION_H_


#include <blasfeo.h>

#ifdef __cplusplus
extern "C" {
#endif


void s_apply_givens_rotation_to_mat(struct blasfeo_smat *mat, 
                                    struct blasfeo_svec *givens_c_vec, 
                                    struct blasfeo_svec *givens_s_vec, 
                                    int applied_column, int applied_row);

void s_apply_givens_rotation_to_vec(struct blasfeo_svec *mat, 
                                    struct blasfeo_svec *givens_c_vec, 
                                    struct blasfeo_svec *givens_s_vec, 
                                    int applied_row);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // HPCGMRES_S_GIVENS_ROTATION_H_