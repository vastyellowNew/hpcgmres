#include "d_givens_rotation.h"


#define REAL double

#define APPLY_GIVENS_ROTATION_TO_VEC d_apply_givens_rotation_to_vec
#define APPLY_GIVENS_ROTATION_TO_MAT d_apply_givens_rotation_to_mat

#define STRVEC blasfeo_dvec
#define STRMAT blasfeo_dmat
#define VECIN1 blasfeo_dvecin1
#define VECEX1 blasfeo_dvecex1
#define MATIN1 blasfeo_dgein1 
#define MATEX1 blasfeo_dgeex1


#include "x_givens_rotation.c"