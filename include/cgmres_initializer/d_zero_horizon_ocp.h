#ifndef HPCGMRES_D_ZERO_HORIZON_OCP_H_
#define HPCGMRES_D_ZERO_HORIZON_OCP_H_


#include <blasfeo.h>

#include "d_nmpc_model.h"


#ifdef __cplusplus
extern "C" {
#endif


struct d_zero_horizon_ocp {
  struct d_nmpc_model *model;
  struct blasfeo_dvec *lmb_vec;
  int dimx;
  int dimu;
  int dimc;
  int dimuc;
  int dim_solution;
  int memsize; // memory size in bytes
};

// returns the size of the strage
int d_zero_horizon_ocp_strsize();

// returns the size of the strage
int d_zero_horizon_ocp_memsize(struct d_zero_horizon_ocp *ocp);

void d_zero_horizon_ocp_create(struct d_zero_horizon_ocp *ocp, void *memory);

void d_zero_horizon_ocp_compute_optimality_residual(
    struct d_zero_horizon_ocp *ocp, double current_time, 
    struct blasfeo_dvec *current_state, struct blasfeo_dvec *solution, 
    struct blasfeo_dvec *optimality_residual);

void d_zero_horizon_ocp_compute_terminal_cost_derivative(
    struct d_zero_horizon_ocp *ocp, double current_time, 
    struct blasfeo_dvec *current_state, 
    struct blasfeo_dvec *terminal_cost_derivative);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // HPCGMRES_D_ZERO_HORIZON_OCP_H_