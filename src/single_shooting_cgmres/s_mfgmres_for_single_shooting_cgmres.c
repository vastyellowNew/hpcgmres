#include "single_shooting_cgmres/s_mfgmres_for_single_shooting_cgmres.h"


#define MFGMRES_STRSIZE s_mfgmres_for_single_shooting_cgmres_strsize
#define MFGMRES_MEMSIZE s_mfgmres_for_single_shooting_cgmres_memsize
#define MFGMRES_CREATE s_mfgmres_for_single_shooting_cgmres_create
#define MFGMRES_DELETE s_mfgmres_for_single_shooting_cgmres_delete
#define MFGMRES_SOLVE_LINEAR_PROBLEM s_mfgmres_for_single_shooting_cgmres_solve_linear_problem
#define MFGMRES s_mfgmres_for_single_shooting_cgmres

#define LINEAR_PROBLEM_COMPUTE_AX s_single_shooting_continuation_compute_ax
#define LINEAR_PROBLEM_COMPUTE_B s_single_shooting_continuation_compute_b
#define LINEAR_PROBLEM_ARGS s_single_shooting_continuation_mfgmres_args
#define LINEAR_PROBLEM s_single_shooting_continuation


#include "../common/s_mfgmres.c"