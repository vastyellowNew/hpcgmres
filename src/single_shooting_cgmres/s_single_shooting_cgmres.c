#include "s_single_shooting_cgmres.h"

#include "s_mfgmres_for_single_shooting_cgmres.h"
#include "s_single_shooting_continuation.h"
#include "s_single_shooting_continuation_mfgmres_args.h"
#include "s_linear_algebra.h"
#include "s_memory_manager.h"


#define REAL float

#define SINGLE_SHOOTING_CGMRES_STRSIZE s_single_shooting_cgmres_strsize
#define SINGLE_SHOOTING_CGMRES_MEMSIZE s_single_shooting_cgmres_memsize
#define SINGLE_SHOOTING_CGMRES_CREATE s_single_shooting_cgmres_create
#define SINGLE_SHOOTING_CGMRES_DELETE s_single_shooting_cgmres_delete
#define SINGLE_SHOOTING_CGMRES_SET_INITIALIZATION_PARAMETERS s_single_shooting_cgmres_set_initialization_parameters
#define SINGLE_SHOOTING_CGMRES_INITIALIZE_SOLUTION s_single_shooting_cgmres_initialize_solution
#define SINGLE_SHOOTING_CGMRES_UPDATE_CONTROL_INPUT s_single_shooting_cgmres_update_control_input
#define SINGLE_SHOOTING_CGMRES_GET_CONTROL_INPUT s_single_shooting_cgmres_get_control_input
#define SINGLE_SHOOTING_CGMRES_GET_ERROR_NORM s_single_shooting_cgmres_get_error_norm
#define SINGLE_SHOOTING_CGMRES s_single_shooting_cgmres

#define SINGLE_SHOOTING_CONTINUATION_MEMSIZE s_single_shooting_continuation_memsize
#define SINGLE_SHOOTING_CONTINUATION_CREATE s_single_shooting_continuation_create
#define SINGLE_SHOOTING_CONTINUATION_DELETE s_single_shooting_continuation_delete
#define SINGLE_SHOOTING_CONTINUATION_INTEGRATE_SOLUTION s_single_shooting_continuation_integrate_solution
#define SINGLE_SHOOTING_CONTINUATION_GET_ERROR_NORM s_single_shooting_continuation_get_error_norm
#define SINGLE_SHOOTING_CONTINUATION_RESET_HORIZON_LENGTH s_single_shooting_continuation_reset_horizon_length
#define SINGLE_SHOOTING_CONTINUATION_DIMX s_single_shooting_continuation_dimx
#define SINGLE_SHOOTING_CONTINUATION_DIMU s_single_shooting_continuation_dimu
#define SINGLE_SHOOTING_CONTINUATION_DIMC s_single_shooting_continuation_dimc

#define MFGMRES_FOR_SINGLE_SHOOTING_CGMRES_MEMSIZE s_mfgmres_for_single_shooting_cgmres_memsize
#define MFGMRES_FOR_SINGLE_SHOOTING_CGMRES_CREATE s_mfgmres_for_single_shooting_cgmres_create
#define MFGMRES_FOR_SINGLE_SHOOTING_CGMRES_DELETE s_mfgmres_for_single_shooting_cgmres_delete
#define MFGMRES_FOR_SINGLE_SHOOTING_CGMRES_SOLVE_LINEAR_PROBLEM s_mfgmres_for_single_shooting_cgmres_solve_linear_problem

#define CGMRES_INITIALIZER_MEMSIZE s_cgmres_initializer_memsize
#define CGMRES_INITIALIZER_CREATE s_cgmres_initializer_create
#define CGMRES_INITIALIZER_DELETE s_cgmres_initializer_delete
#define CGMRES_INITIALIZER_SET_TERMINATION_CRITERIONS s_cgmres_initializer_set_termination_criterions
#define CGMRES_INITIALIZER_SET_INITIAL_GUESS_SOLUTION s_cgmres_initializer_set_initial_guess_solution
#define CGMRES_INITIALIZER_COMPUTE_INITIAL_SOLUTION s_cgmres_initializer_compute_initial_solution

#define VECCP hpcgmres_sveccp

#define ALLOCATE_VEC allocate_svec
#define FREE_VEC free_svec


#include "x_single_shooting_cgmres.c"