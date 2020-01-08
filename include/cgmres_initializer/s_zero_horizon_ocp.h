#ifndef HPCGMRES_S_ZERO_HORIZON_OCP_H_
#define HPCGMRES_S_ZERO_HORIZON_OCP_H_


#include <blasfeo.h>

#include "s_nmpc_model.h"


#ifdef __cplusplus
extern "C" {
#endif


struct s_zero_horizon_ocp {
  struct s_nmpc_model model;
  struct blasfeo_svec *lmd_vec;
  int dimx;
  int dimu;
  int dimc;
  int dimuc;
  int dim_solution;
  int memsize; // memory size in bytes
};

// returns the size of the strage
int s_zero_horizon_ocp_strsize();

// returns the size of the strage
int s_zero_horizon_ocp_memsize();

void s_zero_horizon_ocp_create(struct s_zero_horizon_ocp *ocp, void *memory);

void s_zero_horizon_ocp_compute_optimality_residual(
    struct s_zero_horizon_ocp *ocp, float current_time, 
    struct blasfeo_svec *current_state, struct blasfeo_svec *solution, 
    struct blasfeo_svec *optimality_residual);

void s_zero_horizon_ocp_compute_terminal_cost_derivative(
    struct s_zero_horizon_ocp *ocp, float current_time, 
    struct blasfeo_svec *current_state, 
    struct blasfeo_svec *terminal_cost_derivative);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // HPCGMRES_S_ZERO_HORIZON_OCP_H_