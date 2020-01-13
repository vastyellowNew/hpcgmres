#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "common/s_memory_manager.h"
#include "common/s_linear_algebra.h"
#include "common/s_givens_rotation.h"


#define REAL float
#define MACHINE_EPSILON FLT_EPSILON

#define APPLY_GIVENS_ROTATION s_apply_givens_rotation

#define ALLOCATE_MAT allocate_smat
#define ALLOCATE_VEC allocate_svec
#define FREE_MAT free_smat
#define FREE_VEC free_svec

#define VECSET cgmres_svecset
#define VECNRM2 cgmres_svecnrm2
#define VECMCP cgmres_svecmcp
#define VECDOT cgmres_svecdot
#define VECMAD cgmres_svecmadd
#define VECMUL cgmres_svecmul


#include "x_mfgmres.c"