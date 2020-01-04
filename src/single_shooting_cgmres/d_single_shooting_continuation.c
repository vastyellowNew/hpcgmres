#include <stdlib.h>
#include <stdio.h>
#include <cmath.h>

#include <blasfeo.h>

#include "d_single_shooting_continuation.h"


#define REAL double

#define SINGLE_SHOOTING_CONTINUATION d_single_shooting_continuation
#define SINGLE_SHOOTING_CONTINUATION_STRSIZE d_single_shooting_continuation_strsize
#define SINGLE_SHOOTING_CONTINUATION_MEMSIZE d_single_shooting_continuation_memsize
#define SINGLE_SHOOTING_CONTINUATION_CREATE d_single_shooting_continuation_create
#define SINGLE_SHOOTING_CONTINUATION_INTEGRATE_SOLUTION d_single_shooting_continuation_integrate_solution
#define SINGLE_SHOOTING_CONTINUATION_COMPUTE_B d_single_shooting_continuation_compute_b
#define SINGLE_SHOOTING_CONTINUATION_COMPUTE_AX d_single_shooting_continuation_compute_ax
#define SINGLE_SHOOTING_CONTINUATION_GET_ERROR_NORM d_single_shooting_continuation_get_error_norm
#define SINGLE_SHOOTING_CONTINUATION_RESET_HORIZON_LENGTH d_single_shooting_continuation_reset_horizon_length

#define SINGLE_SHOOTING_OCP_MEMSIZE d_single_shooting_ocp_memsize
#define SINGLE_SHOOTING_OCP_CREATE d_single_shooting_ocp_create
#define SINGLE_SHOOTING_OCP_COMPUTE_OPTIMALITY_RESIDUAL d_single_shooting_ocp_compute_optimality_residual
#define SINGLE_SHOOTING_OCP_PREDICT_STATE_FROM_SOLUTION d_single_shooting_ocp_predict_state_from_solution

#define STRVEC blasfeo_dvec
#define SIZE_STRVEC blasfeo_memsize_dvec
#define CREATE_STRVEC blasfeo_memsize_dvec

#define VECCSE blasfeo_dvecse
#define VECAD blasfeo_dvecad
#define VECDOT blasfeo_ddot
#define AXPY blasfeo_daxpy
#define AXBPY blasfeo_daxbpy


#include "x_single_shooting_continuation.c"