#ifndef HPCGMRES_D_INEXACT_NEWTON_FOR_ZERO_HORIZON_COP_H_
#define HPCGMRES_D_INEXACT_NEWTON_FOR_ZERO_HORIZON_COP_H_


#include <cmath>
#include <blasfeo.h>

#include "d_zero_horizon_ocp.h"
#include "d_inexact_newton_for_zero_horizon_ocp_mfgmres_args.h"


#ifdef __cplusplus
extern "C" {
#endif


struct d_inexact_newton_for_zero_horizon_ocp {
  struct d_zero_horizon_ocp *ocp;
  struct blasfeo_dvec *incremented_solution;
  struct blasfeo_dvec *optimality_residual;
  struct blasfeo_dvec *optimality_residual1;
  double finite_difference_increment; 
  int dim_solution;
  int memsize; // memory size in bytes
};

int d_inexact_newton_for_zero_horizon_ocp_strsize();

int d_inexact_newton_for_zero_horizon_ocp_memsize(
    struct d_inexact_newton_for_zero_horizon_ocp *newton);

void d_inexact_newton_for_zero_horizon_ocp_create(
    struct d_inexact_newton_for_zero_horizon_ocp *newton, 
    double finite_difference_increment, void *memory);

void d_inexact_newton_for_zero_horizon_ocp_compute_b(
    struct d_inexact_newton_for_zero_horizon_ocp *newton, 
    struct d_inexact_newton_for_zero_horizon_ocp_mfgmres_args *args,
    struct blasfeo_dvec *direction, struct blasfeo_dvec *b);

void d_inexact_newton_for_zero_horizon_ocp_compute_ax(
    struct d_inexact_newton_for_zero_horizon_ocp *newton, 
    struct d_inexact_newton_for_zero_horizon_ocp_mfgmres_args *args,
    struct blasfeo_dvec *direction, struct blasfeo_dvec *ax);

double d_inexact_newton_for_zero_horizon_ocp_get_error_norm(
    struct d_inexact_newton_for_zero_horizon_ocp *newton, double initial_time, 
    struct blasfeo_dvec *initial_state, struct blasfeo_dvec *initial_solution);

void d_inexact_newton_for_zero_horizon_ocp_get_terminal_cost_derivative(
    struct d_inexact_newton_for_zero_horizon_ocp *newton, double initial_time, 
    struct blasfeo_dvec *initial_state, 
    struct blasfeo_dvec *terminal_cost_derivative);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // HPCGMRES_D_INEXACT_NEWTON_FOR_ZERO_HORIZON_COP_H_