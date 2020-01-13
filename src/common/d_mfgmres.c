#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "common/d_memory_manager.h"
#include "common/d_linear_algebra.h"
#include "common/d_givens_rotation.h"


#define REAL double
#define MACHINE_EPSILON DBL_EPSILON

#define APPLY_GIVENS_ROTATION d_apply_givens_rotation

#define ALLOCATE_MAT allocate_dmat
#define ALLOCATE_VEC allocate_dvec
#define FREE_MAT free_dmat
#define FREE_VEC free_dvec

#define VECSET cgmres_dvecset
#define VECNRM2 cgmres_dvecnrm2
#define VECMCP cgmres_dvecmcp
#define VECDOT cgmres_dvecdot
#define VECMAD cgmres_dvecmadd
#define VECMUL cgmres_dvecmul


#include "x_mfgmres.c"