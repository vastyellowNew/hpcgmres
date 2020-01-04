#ifndef CGMRES_S_SINGLE_SHOOTING_OCP_H_
#define CGMRES_S_SINGLE_SHOOTING_OCP_H_


#include <blasfeo.h>

#include "s_nmpc_model.h"
#include "s_time_varying_smooth_horizon.h"


#ifdef __cplusplus
extern "C" {
#endif


struct s_single_shooting_ocp {
  struct s_nmpc_model *model;
  struct s_time_varying_smooth_horizon *horizon;
  struct blasfeo_svec *dx_vec;
  struct blasfeo_svec *x_sequence_vec;
  struct blasfeo_svec *lmd_sequence_vec;
  int dimx;
  int dimu;
  int dimc;
  int dimuc;
  int N;
  int dim_solution;
  int memsize; // memory size in bytes
};

int s_single_shooting_ocp_strsize();

int s_single_shooting_ocp_memsize(struct s_single_shooting_ocp *ocp, int N);


void s_single_shooting_ocp_create(struct s_single_shooting_ocp *ocp,
                                  float T_f, float alpha, float initial_time, 
                                  int N, void *memory);

void s_single_shooting_ocp_compute_optimality_residual(
    struct s_single_shooting_ocp *ocp, float current_time, 
    struct blasfeo_svec *current_state, struct blasfeo_svec *solution, 
    struct blasfeo_svec *optimality_residual);

void s_single_shooting_ocp_predict_state_from_solution(
    struct s_single_shooting_ocp *ocp, float current_time, 
    struct blasfeo_svec *current_state, struct blasfeo_svec *solution, 
    float prediction_length, struct blasfeo_svec *predicted_state);

void s_single_shooting_ocp_reset_horizon_length(
    struct s_single_shooting_ocp *ocp, float initial_time);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // CGMRES_s_single_SHOOTING_OCP_H_