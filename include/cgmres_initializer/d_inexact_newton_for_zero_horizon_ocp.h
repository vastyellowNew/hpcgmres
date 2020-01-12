#ifndef HPCGMRES_D_INEXACT_NEWTON_FOR_ZERO_HORIZON_COP_H_
#define HPCGMRES_D_INEXACT_NEWTON_FOR_ZERO_HORIZON_COP_H_


#include "cgmres_initializer/d_zero_horizon_ocp.h"
#include "cgmres_initializer/d_inexact_newton_for_zero_horizon_ocp_mfgmres_args.h"


#ifdef __cplusplus
extern "C" {
#endif


struct d_inexact_newton_for_zero_horizon_ocp {
  struct d_zero_horizon_ocp ocp;
  double *incremented_solution;
  double *optimality_residual;
  double *optimality_residual1;
  double finite_difference_increment; 
  int dim_solution;
  int memsize; // memory size in bytes
};

int d_inexact_newton_for_zero_horizon_ocp_strsize();

int d_inexact_newton_for_zero_horizon_ocp_memsize();

void d_inexact_newton_for_zero_horizon_ocp_create(
    struct d_inexact_newton_for_zero_horizon_ocp *newton, 
    double finite_difference_increment);

void d_inexact_newton_for_zero_horizon_ocp_delete(
    struct d_inexact_newton_for_zero_horizon_ocp *newton);

void d_inexact_newton_for_zero_horizon_ocp_compute_b(
    struct d_inexact_newton_for_zero_horizon_ocp *newton, 
    struct d_inexact_newton_for_zero_horizon_ocp_mfgmres_args *args,
    double *direction, double *b);

void d_inexact_newton_for_zero_horizon_ocp_compute_ax(
    struct d_inexact_newton_for_zero_horizon_ocp *newton, 
    struct d_inexact_newton_for_zero_horizon_ocp_mfgmres_args *args,
    double *direction, double *ax);

double d_inexact_newton_for_zero_horizon_ocp_get_error_norm(
    struct d_inexact_newton_for_zero_horizon_ocp *newton, double initial_time, 
    double *initial_state, double *initial_solution);

void d_inexact_newton_for_zero_horizon_ocp_get_terminal_cost_derivative(
    struct d_inexact_newton_for_zero_horizon_ocp *newton, double initial_time, 
    double *initial_state, double *terminal_cost_derivative);

int d_inexact_newton_for_zero_horizon_ocp_dimx();

int d_inexact_newton_for_zero_horizon_ocp_dimu();

int d_inexact_newton_for_zero_horizon_ocp_dimc();


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // HPCGMRES_D_INEXACT_NEWTON_FOR_ZERO_HORIZON_COP_H_