#ifndef HPCGMRES_S_ZERO_HORIZON_OCP_H_
#define HPCGMRES_S_ZERO_HORIZON_OCP_H_


#include "s_nmpc_model.h"


#ifdef __cplusplus
extern "C" {
#endif


struct s_zero_horizon_ocp {
  struct s_nmpc_model model;
  float *lmd_vec;
  int dimx;
  int dimu;
  int dimc;
  int dimuc;
  int dim_solution;
  int memsize; // memory size in bytes
};

int s_zero_horizon_ocp_strsize();

int s_zero_horizon_ocp_memsize();

void s_zero_horizon_ocp_create(struct s_zero_horizon_ocp *ocp);

void s_zero_horizon_ocp_delete(struct s_zero_horizon_ocp *ocp);

void s_zero_horizon_ocp_compute_optimality_residual(
    struct s_zero_horizon_ocp *ocp, float current_time, float *current_state, 
    float *solution, float *optimality_residual);

void s_zero_horizon_ocp_compute_terminal_cost_derivative(
    struct s_zero_horizon_ocp *ocp, float current_time, float *current_state, 
    float *terminal_cost_derivative);

int s_zero_horizon_ocp_dimx();

int s_zero_horizon_ocp_dimu();

int s_zero_horizon_ocp_dimc();


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // HPCGMRES_S_ZERO_HORIZON_OCP_H_