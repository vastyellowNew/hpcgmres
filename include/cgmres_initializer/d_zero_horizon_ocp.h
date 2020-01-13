#ifndef CGMRES_D_ZERO_HORIZON_OCP_H_
#define CGMRES_D_ZERO_HORIZON_OCP_H_


#include "d_nmpc_model.h"


#ifdef __cplusplus
extern "C" {
#endif


struct d_zero_horizon_ocp {
  struct d_nmpc_model model;
  double *lmd_vec;
  int dimx;
  int dimu;
  int dimc;
  int dimuc;
  int dim_solution;
  int memsize; // memory size in bytes
};

int d_zero_horizon_ocp_strsize();

int d_zero_horizon_ocp_memsize();

void d_zero_horizon_ocp_create(struct d_zero_horizon_ocp *ocp);

void d_zero_horizon_ocp_delete(struct d_zero_horizon_ocp *ocp);

void d_zero_horizon_ocp_compute_optimality_residual(
    struct d_zero_horizon_ocp *ocp, double current_time, double *current_state, 
    double *solution, double *optimality_residual);

void d_zero_horizon_ocp_compute_terminal_cost_derivative(
    struct d_zero_horizon_ocp *ocp, double current_time, double *current_state, 
    double *terminal_cost_derivative);

int d_zero_horizon_ocp_dimx();

int d_zero_horizon_ocp_dimu();

int d_zero_horizon_ocp_dimc();


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // CGMRES_D_ZERO_HORIZON_OCP_H_