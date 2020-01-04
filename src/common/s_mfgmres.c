#include <stdlib.h>
#include <stdio.h>
#include <float.h>

#include <blasfeo.h>

#include "s_mfgmres.h"


#define REAL float
#define MACHINE_EPSILON FLT_EPSILON

#define MFGMRES_GIVENS_ROTATION_MAT s_mfgmres_givens_rotation_mat
#define MFGMRES_GIVENS_ROTATION_VEC s_mfgmres_givens_rotation_vec

#define STRVEC blasfeo_svec
#define STRMAT blasfeo_smat
#define SIZE_STRMAT blasfeo_memsize_smat
#define SIZE_STRVEC blasfeo_memsize_svec
#define CREATE_STRMAT blasfeo_memsize_smat
#define CREATE_STRVEC blasfeo_memsize_svec

#define VECIN1 blasfeo_svecin1
#define VECEX1 blasfeo_svecex1
#define VECCSE blasfeo_svecse
#define MATIN1 blasfeo_sgein1 
#define MATEX1 blasfeo_sgeex1
#define MATCSE blasfeo_sgese
#define VECSC blasfeo_svecsc
#define VECCAD blasfeo_saxpby
#define VECCPSC blasfeo_sveccpsc
#define VECDOT blasfeo_sdot
#define TRSV_UNN blasfeo_strsv_unn


#include "x_mfgmres.c"