#include <stdio.h>

#include <blasfeo.h>

#include "s_zero_horizon_ocp.h"


#define REAL float

#define ZERO_HORIZON_OCP s_zero_horizon_ocp
#define ZERO_HORIZON_OCP_STRSIZE s_zero_horizon_ocp_strsize
#define ZERO_HORIZON_OCP_MEMSIZE s_zero_horizon_ocp_memsize
#define ZERO_HORIZON_OCP_CREATE s_zero_horizon_ocp_create
#define ZERO_HORIZON_OCP_COMPUTE_OPTIMALITY_RESIDUAL s_zero_horizon_ocp_compute_optimality_residual
#define ZERO_HORIZON_OCP_COMPUTE_TERMINAL_COST_DERIVATIVE s_zero_horizon_ocp_compute_terminal_cost_derivative

#define NMPC_MODEL s_nmpc_model
#define NMPC_MODEL_CREATE s_nmpc_model_create
#define NMPC_MODEL_PHIX s_nmpc_model_phix
#define NMPC_MODEL_HU s_nmpc_model_hu
#define NMPC_MODEL_DIMX s_nmpc_model_dimx
#define NMPC_MODEL_DIMU s_nmpc_model_dimu
#define NMPC_MODEL_DIMC s_nmpc_model_dimc

#define STRVEC blasfeo_svec
#define SIZE_STRVEC blasfeo_memsize_svec
#define CREATE_STRVEC blasfeo_create_svec
#define VECSE blasfeo_svecse


#include "x_zero_horizon_ocp.c"