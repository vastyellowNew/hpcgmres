#ifndef CGMRES_S_INEXACT_NEWTON_FOR_ZERO_HORIZON_COP_H_
#define CGMRES_S_INEXACT_NEWTON_FOR_ZERO_HORIZON_COP_H_


#include "cgmres_initializer/s_zero_horizon_ocp.h"
#include "cgmres_initializer/s_inexact_newton_for_zero_horizon_ocp_mfgmres_args.h"


#ifdef __cplusplus
extern "C" {
#endif


struct s_inexact_newton_for_zero_horizon_ocp {
  struct s_zero_horizon_ocp ocp;
  float *incremented_solution;
  float *optimality_residual;
  float *optimality_residual1;
  float finite_difference_increment; 
  int dim_solution;
  int memsize; // memory size in bytes
};

int s_inexact_newton_for_zero_horizon_ocp_strsize();

int s_inexact_newton_for_zero_horizon_ocp_memsize();

void s_inexact_newton_for_zero_horizon_ocp_create(
    struct s_inexact_newton_for_zero_horizon_ocp *newton, 
    float finite_difference_increment);

void s_inexact_newton_for_zero_horizon_ocp_delete(
    struct s_inexact_newton_for_zero_horizon_ocp *newton);

void s_inexact_newton_for_zero_horizon_ocp_compute_b(
    struct s_inexact_newton_for_zero_horizon_ocp *newton, 
    struct s_inexact_newton_for_zero_horizon_ocp_mfgmres_args *args,
    float *direction, float *b);

void s_inexact_newton_for_zero_horizon_ocp_compute_ax(
    struct s_inexact_newton_for_zero_horizon_ocp *newton, 
    struct s_inexact_newton_for_zero_horizon_ocp_mfgmres_args *args,
    float *direction, float *ax);

float s_inexact_newton_for_zero_horizon_ocp_get_error_norm(
    struct s_inexact_newton_for_zero_horizon_ocp *newton, float initial_time, 
    float *initial_state, float *initial_solution);

void s_inexact_newton_for_zero_horizon_ocp_get_terminal_cost_derivative(
    struct s_inexact_newton_for_zero_horizon_ocp *newton, float initial_time, 
    float *initial_state, float *terminal_cost_derivative);

int s_inexact_newton_for_zero_horizon_ocp_dimx();

int s_inexact_newton_for_zero_horizon_ocp_dimu();

int s_inexact_newton_for_zero_horizon_ocp_dimc();


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // CGMRES_S_INEXACT_NEWTON_FOR_ZERO_HORIZON_COP_H_