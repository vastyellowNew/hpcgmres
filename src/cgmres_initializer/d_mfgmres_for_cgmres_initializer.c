#include "d_mfgmres_for_cgmres_initializer.h"


#define LINEAR_PROBLEM d_inexact_newton_for_zero_horizon_ocp
#define LINEAR_PROBLEM_ARGS d_inexact_newton_for_zero_horizon_ocp_mfgmres_args


#include "../common/d_mfgmres.c"