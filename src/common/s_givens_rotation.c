#include "s_givens_rotation.h"


#define REAL float

#define APPLY_GIVENS_ROTATION_TO_VEC s_apply_givens_rotation_to_vec
#define APPLY_GIVENS_ROTATION_TO_MAT s_apply_givens_rotation_to_mat

#define STRVEC blasfeo_svec
#define STRMAT blasfeo_smat
#define VECIN1 blasfeo_svecin1
#define VECEX1 blasfeo_svecex1
#define MATIN1 blasfeo_sgein1 
#define MATEX1 blasfeo_sgeex1


#include "x_givens_rotation.c"