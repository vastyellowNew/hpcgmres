#include "cgmres_initializer/d_mfgmres_for_cgmres_initializer.h"


#define MFGMRES_STRSIZE d_mfgmres_for_cgmres_initializer_strsize
#define MFGMRES_MEMSIZE d_mfgmres_for_cgmres_initializer_memsize
#define MFGMRES_CREATE d_mfgmres_for_cgmres_initializer_create
#define MFGMRES_DELETE d_mfgmres_for_cgmres_initializer_delete
#define MFGMRES_SOLVE_LINEAR_PROBLEM d_mfgmres_for_cgmres_initializer_solve_linear_problem
#define MFGMRES d_mfgmres_for_cgmres_initializer

#define LINEAR_PROBLEM_COMPUTE_AX d_inexact_newton_for_zero_horizon_ocp_compute_ax
#define LINEAR_PROBLEM_COMPUTE_B d_inexact_newton_for_zero_horizon_ocp_compute_b
#define LINEAR_PROBLEM_ARGS d_inexact_newton_for_zero_horizon_ocp_mfgmres_args
#define LINEAR_PROBLEM d_inexact_newton_for_zero_horizon_ocp


#include "../common/d_mfgmres.c"