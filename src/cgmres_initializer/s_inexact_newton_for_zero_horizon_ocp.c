#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "cgmres_initializer/s_inexact_newton_for_zero_horizon_ocp.h"
#include "common/s_memory_manager.h"
#include "common/s_linear_algebra.h"


#define REAL float

#define INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_STRSIZE s_inexact_newton_for_zero_horizon_ocp_strsize
#define INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_MEMSIZE s_inexact_newton_for_zero_horizon_ocp_memsize
#define INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_CREATE s_inexact_newton_for_zero_horizon_ocp_create
#define INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_DELETE s_inexact_newton_for_zero_horizon_ocp_delete
#define INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_COMPUTE_B s_inexact_newton_for_zero_horizon_ocp_compute_b
#define INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_COMPUTE_AX s_inexact_newton_for_zero_horizon_ocp_compute_ax
#define INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_GET_ERROR_NORM s_inexact_newton_for_zero_horizon_ocp_get_error_norm
#define INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_GET_TERMINAL_COST_DERIVATIVE s_inexact_newton_for_zero_horizon_ocp_get_terminal_cost_derivative
#define INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_DIMX s_inexact_newton_for_zero_horizon_ocp_dimx
#define INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_DIMU s_inexact_newton_for_zero_horizon_ocp_dimu
#define INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_DIMC s_inexact_newton_for_zero_horizon_ocp_dimc

#define INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_MFGMRES_ARGS s_inexact_newton_for_zero_horizon_ocp_mfgmres_args
#define INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP s_inexact_newton_for_zero_horizon_ocp

#define ZERO_HORIZON_OCP_STRSIZE s_zero_horizon_ocp_strsize 
#define ZERO_HORIZON_OCP_MEMSIZE s_zero_horizon_ocp_memsize
#define ZERO_HORIZON_OCP_CREATE s_zero_horizon_ocp_create
#define ZERO_HORIZON_OCP_DELETE s_zero_horizon_ocp_delete
#define ZERO_HORIZON_OCP_COMPUTE_OPTIMALITY_RESIDUAL s_zero_horizon_ocp_compute_optimality_residual
#define ZERO_HORIZON_OCP_COMPUTE_TERMINAL_COST_DERIVATIVE s_zero_horizon_ocp_compute_terminal_cost_derivative
#define ZERO_HORIZON_OCP_DIMX s_zero_horizon_ocp_dimx
#define ZERO_HORIZON_OCP_DIMU s_zero_horizon_ocp_dimu
#define ZERO_HORIZON_OCP_DIMC s_zero_horizon_ocp_dimc

#define ALLOCATE_VEC allocate_svec
#define FREE_VEC free_svec

#define VECNRM2 hpcgmres_svecnrm2
#define AXPY hpcgmres_saxpy
#define AXPBY hpcgmres_saxpby


#include "x_inexact_newton_for_zero_horizon_ocp.c"