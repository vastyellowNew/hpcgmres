#include "s_mfgmres_for_cgmres_initializer.h"


#define LINEAR_PROBLEM s_inexact_newton_for_zero_horizon_ocp
#define LINEAR_PROBLEM_ARGS s_inexact_newton_for_zero_horizon_ocp_mfgmres_args


#include "../common/s_mfgmres.c"