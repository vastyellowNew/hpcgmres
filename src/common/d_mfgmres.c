#include <stdlib.h>
#include <stdio.h>
#include <float.h>

#include <blasfeo.h>

#include "d_mfgmres.h"


#define REAL double
#define MACHINE_EPSILON DBL_EPSILON

#define MFGMRES_GIVENS_ROTATION_MAT d_mfgmres_givens_rotation_mat
#define MFGMRES_GIVENS_ROTATION_VEC d_mfgmres_givens_rotation_vec

#define STRVEC blasfeo_dvec
#define STRMAT blasfeo_dmat
#define SIZE_STRMAT blasfeo_memsize_dmat
#define SIZE_STRVEC blasfeo_memsize_dvec
#define CREATE_STRMAT blasfeo_memsize_dmat
#define CREATE_STRVEC blasfeo_memsize_dvec

#define VECIN1 blasfeo_dvecin1
#define VECEX1 blasfeo_dvecex1
#define VECCSE blasfeo_dvecse
#define MATIN1 blasfeo_dgein1 
#define MATEX1 blasfeo_dgeex1
#define MATCSE blasfeo_dgese
#define VECSC blasfeo_dvecsc
#define VECCAD blasfeo_daxpby
#define VECCPSC blasfeo_dveccpsc
#define VECDOT blasfeo_ddot
#define TRSV_UNN blasfeo_dtrsv_unn


#include "x_mfgmres.c"