#include "cgmres_initializer/s_mfgmres_for_cgmres_initializer.h"


#define MFGMRES_STRSIZE s_mfgmres_for_cgmres_initializer_strsize
#define MFGMRES_MEMSIZE s_mfgmres_for_cgmres_initializer_memsize
#define MFGMRES_CREATE s_mfgmres_for_cgmres_initializer_create
#define MFGMRES_DELETE s_mfgmres_for_cgmres_initializer_delete
#define MFGMRES_SOLVE_LINEAR_PROBLEM s_mfgmres_for_cgmres_initializer_solve_linear_problem
#define MFGMRES s_mfgmres_for_cgmres_initializer

#define LINEAR_PROBLEM_COMPUTE_AX s_inexact_newton_for_zero_horizon_ocp_compute_ax
#define LINEAR_PROBLEM_COMPUTE_B s_inexact_newton_for_zero_horizon_ocp_compute_b
#define LINEAR_PROBLEM_ARGS s_inexact_newton_for_zero_horizon_ocp_mfgmres_args
#define LINEAR_PROBLEM s_inexact_newton_for_zero_horizon_ocp


#include "../common/s_mfgmres.c"