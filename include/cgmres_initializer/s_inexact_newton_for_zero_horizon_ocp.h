#ifndef HPCGMRES_S_INEXACT_NEWTON_FOR_ZERO_HORIZON_COP_H_
#define HPCGMRES_S_INEXACT_NEWTON_FOR_ZERO_HORIZON_COP_H_


#include <blasfeo.h>

#include "s_zero_horizon_ocp.h"
#include "s_inexact_newton_for_zero_horizon_ocp_mfgmres_args.h"


#ifdef __cplusplus
extern "C" {
#endif


struct s_inexact_newton_for_zero_horizon_ocp {
  struct s_zero_horizon_ocp ocp;
  struct blasfeo_svec *incremented_solution;
  struct blasfeo_svec *optimality_residual;
  struct blasfeo_svec *optimality_residual1;
  float finite_difference_increment; 
  int dim_solution;
  int memsize; // memory size in bytes
};

int s_inexact_newton_for_zero_horizon_ocp_strsize();

int s_inexact_newton_for_zero_horizon_ocp_memsize();

void s_inexact_newton_for_zero_horizon_ocp_create(
    struct s_inexact_newton_for_zero_horizon_ocp *newton, 
    float finite_difference_increment, void *memory);

void s_inexact_newton_for_zero_horizon_ocp_compute_b(
    struct s_inexact_newton_for_zero_horizon_ocp *newton, 
    struct s_inexact_newton_for_zero_horizon_ocp_mfgmres_args *args,
    struct blasfeo_svec *direction, struct blasfeo_svec *b);

void s_inexact_newton_for_zero_horizon_ocp_compute_ax(
    struct s_inexact_newton_for_zero_horizon_ocp *newton, 
    struct s_inexact_newton_for_zero_horizon_ocp_mfgmres_args *args,
    struct blasfeo_svec *direction, struct blasfeo_svec *ax);

float s_inexact_newton_for_zero_horizon_ocp_get_error_norm(
    struct s_inexact_newton_for_zero_horizon_ocp *newton, float initial_time, 
    struct blasfeo_svec *initial_state, struct blasfeo_svec *initial_solution);

void s_inexact_newton_for_zero_horizon_ocp_get_terminal_cost_derivative(
    struct s_inexact_newton_for_zero_horizon_ocp *newton, float initial_time, 
    struct blasfeo_svec *initial_state, 
    struct blasfeo_svec *terminal_cost_derivative);

int s_inexact_newton_for_zero_horizon_ocp_dimx();

int s_inexact_newton_for_zero_horizon_ocp_dimu();

int s_inexact_newton_for_zero_horizon_ocp_dimc();


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // HPCGMRES_S_INEXACT_NEWTON_FOR_ZERO_HORIZON_COP_H_