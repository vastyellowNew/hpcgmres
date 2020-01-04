#include <stdlib.h>
#include <stdio.h>

#include <blasfeo.h>

#include "d_zero_horizon_ocp.h"


#define REAL double

#define ZERO_HORIZON_OCP d_zero_horizon_ocp
#define ZERO_HORIZON_OCP_STRSIZE d_zero_horizon_ocp_strsize
#define ZERO_HORIZON_OCP_MEMSIZE d_zero_horizon_ocp_memsize
#define ZERO_HORIZON_OCP_CREATE d_zero_horizon_ocp_create
#define ZERO_HORIZON_OCP_COMPUTE_OPTIMALITY_RESIDUAL d_zero_horizon_ocp_compute_optimality_residual
#define ZERO_HORIZON_OCP_COMPUTE_TERMINAL_COST_DERIVATIVE d_zero_horizon_ocp_compute_terminal_cost_derivative

#define NMPC_MODEL d_nmpc_model
#define NMPC_MODEL_CREATE d_nmpc_model_create
#define NMPC_MODEL_PHIX d_nmpc_model_phix
#define NMPC_MODEL_HU d_nmpc_model_hu

#define STRVEC blasfeo_dvec
#define SIZE_STRVEC blasfeo_memsize_dvec
#define CREATE_STRVEC blasfeo_memsize_dvec

#define VECCSE blasfeo_dvecse


#include "x_zero_horizon_ocp.c"