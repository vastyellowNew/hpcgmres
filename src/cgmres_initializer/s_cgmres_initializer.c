#include "s_cgmres_initializer.h"

#include "s_mfgmres_for_cgmres_initializer.h"
#include "s_inexact_newton_for_zero_horizon_ocp.h"
#include "s_memory_manager.h"
#include "s_linear_algebra.h"


#define REAL float

#define CGMRES_INITIALIZER_STRSIZE s_cgmres_initializer_strsize
#define CGMRES_INITIALIZER_MEMSIZE s_cgmres_initializer_memsize
#define CGMRES_INITIALIZER_CREATE s_cgmres_initializer_create
#define CGMRES_INITIALIZER_DELETE s_cgmres_initializer_delete
#define CGMRES_INITIALIZER_SET_TERMINATION_CRITERIONS s_cgmres_initializer_set_termination_criterions
#define CGMRES_INITIALIZER_SET_INITIAL_GUESS_SOLUTION s_cgmres_initializer_set_initial_guess_solution
#define CGMRES_INITIALIZER_COMPUTE_INITIAL_SOLUTION s_cgmres_initializer_compute_initial_solution
#define CGMRES_INITIALIZER_GET_TERMINAL_COST_DERIVATIVE s_cgmres_initializer_get_terminal_cost_derivative
#define CGMRES_INITIALIZER s_cgmres_initializer

#define MFGMRES_FOR_CGMRES_INITIALIZER_MEMSIZE s_mfgmres_for_cgmres_initializer_memsize
#define MFGMRES_FOR_CGMRES_INITIALIZER_CREATE s_mfgmres_for_cgmres_initializer_create
#define MFGMRES_FOR_CGMRES_INITIALIZER_DELETE s_mfgmres_for_cgmres_initializer_delete
#define MFGMRES_FOR_CGMRES_INITIALIZER_SOLVE_LINEAR_PROBLEM s_mfgmres_for_cgmres_initializer_solve_linear_problem

#define INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_MEMSIZE s_inexact_newton_for_zero_horizon_ocp_memsize
#define INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_CREATE s_inexact_newton_for_zero_horizon_ocp_create
#define INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_DELETE s_inexact_newton_for_zero_horizon_ocp_delete
#define INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_GET_ERROR_NORM s_inexact_newton_for_zero_horizon_ocp_get_error_norm
#define INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_GET_TERMINAL_COST_DERIVATIVE s_inexact_newton_for_zero_horizon_ocp_get_terminal_cost_derivative
#define INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_DIMX s_inexact_newton_for_zero_horizon_ocp_dimx
#define INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_DIMU s_inexact_newton_for_zero_horizon_ocp_dimu
#define INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_DIMC s_inexact_newton_for_zero_horizon_ocp_dimc

#define ALLOCATE_VEC allocate_svec
#define FREE_VEC free_svec

#define VECCP hpcgmres_sveccp
#define VECAD hpcgmres_svecadd


#include "x_cgmres_initializer.c"