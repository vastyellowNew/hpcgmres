#ifndef CGMRES_S_SINGLE_SHOOTING_OCP_H_
#define CGMRES_S_SINGLE_SHOOTING_OCP_H_


#include "s_nmpc_model.h"
#include "common/s_time_varying_smooth_horizon.h"


#ifdef __cplusplus
extern "C" {
#endif


struct s_single_shooting_ocp {
  struct s_nmpc_model model;
  struct s_time_varying_smooth_horizon horizon;
  float *dx;
  float **x_sequence;
  float **lmd_sequence;
  int dimx;
  int dimu;
  int dimc;
  int dimuc;
  int N;
  int dim_solution;
  int memsize; // memory size in bytes
};

int s_single_shooting_ocp_strsize();

int s_single_shooting_ocp_memsize(int N);

void s_single_shooting_ocp_create(struct s_single_shooting_ocp *ocp, float T_f, 
                                  float alpha, float initial_time, int N);

void s_single_shooting_ocp_delete(struct s_single_shooting_ocp *ocp);

void s_single_shooting_ocp_compute_optimality_residual(
    struct s_single_shooting_ocp *ocp, float current_time, float *current_state, 
    float *solution, float *optimality_residual);

void s_single_shooting_ocp_predict_state_from_solution(
    struct s_single_shooting_ocp *ocp, float current_time, float *current_state, 
    float *solution, float prediction_length, float *predictes_state);

void s_single_shooting_ocp_reset_horizon_length(
    struct s_single_shooting_ocp *ocp, float initial_time);

int s_single_shooting_ocp_dimx();

int s_single_shooting_ocp_dimu();

int s_single_shooting_ocp_dimc();


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // CGMRES_S_SINGLE_SHOOTING_OCP_H_