#include <stdlib.h>
#include <math.h>

#include "d_single_shooting_continuation.h"
#include "d_single_shooting_ocp.h"
#include "d_linear_algebra.h"
#include "d_memory_manager.h"


#define REAL double

#define SINGLE_SHOOTING_CONTINUATION_STRSIZE d_single_shooting_continuation_strsize
#define SINGLE_SHOOTING_CONTINUATION_MEMSIZE d_single_shooting_continuation_memsize
#define SINGLE_SHOOTING_CONTINUATION_CREATE d_single_shooting_continuation_create
#define SINGLE_SHOOTING_CONTINUATION_DELETE d_single_shooting_continuation_delete
#define SINGLE_SHOOTING_CONTINUATION_INTEGRATE_SOLUTION d_single_shooting_continuation_integrate_solution
#define SINGLE_SHOOTING_CONTINUATION_COMPUTE_B d_single_shooting_continuation_compute_b
#define SINGLE_SHOOTING_CONTINUATION_COMPUTE_AX d_single_shooting_continuation_compute_ax
#define SINGLE_SHOOTING_CONTINUATION_GET_ERROR_NORM d_single_shooting_continuation_get_error_norm
#define SINGLE_SHOOTING_CONTINUATION_RESET_HORIZON_LENGTH d_single_shooting_continuation_reset_horizon_length

#define SINGLE_SHOOTING_CONTINUATION_MFGMRES_ARGS d_single_shooting_continuation_mfgmres_args
#define SINGLE_SHOOTING_CONTINUATION d_single_shooting_continuation

#define SINGLE_SHOOTING_OCP_MEMSIZE d_single_shooting_ocp_memsize
#define SINGLE_SHOOTING_OCP_CREATE d_single_shooting_ocp_create
#define SINGLE_SHOOTING_OCP_DELETE d_single_shooting_ocp_delete
#define SINGLE_SHOOTING_OCP_COMPUTE_OPTIMALITY_RESIDUAL d_single_shooting_ocp_compute_optimality_residual
#define SINGLE_SHOOTING_OCP_PREDICT_STATE_FROM_SOLUTION d_single_shooting_ocp_predict_state_from_solution
#define SINGLE_SHOOTING_OCP_RESET_HORIZON_LENGTH d_single_shooting_ocp_reset_horizon_length
#define SINGLE_SHOOTING_OCP_DIMX d_single_shooting_ocp_dimx
#define SINGLE_SHOOTING_OCP_DIMU d_single_shooting_ocp_dimu
#define SINGLE_SHOOTING_OCP_DIMC d_single_shooting_ocp_dimc

#define VECMAD hpcgmres_dvecmadd
#define AXPY hpcgmres_daxpy
#define AXPBY hpcgmres_daxpby
#define VECNRM2 hpcgmres_dvecnrm2

#define ALLOCATE_VEC allocate_dvec
#define FREE_VEC free_dvec


#include "x_single_shooting_continuation.c"