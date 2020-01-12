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

#define VECSET hpcgmres_svecset
#define VECNRM2 hpcgmres_svecnrm2
#define VECMCP hpcgmres_svecmcp
#define VECDOT hpcgmres_svecdot
#define VECMAD hpcgmres_svecmadd
#define VECMUL hpcgmres_svecmul


#include "x_mfgmres.c"