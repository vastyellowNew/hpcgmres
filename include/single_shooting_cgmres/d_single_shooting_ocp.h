#ifndef CGMRES_D_SINGLE_SHOOTING_OCP_H_
#define CGMRES_D_SINGLE_SHOOTING_OCP_H_


#include "d_nmpc_model.h"
#include "common/d_time_varying_smooth_horizon.h"


#ifdef __cplusplus
extern "C" {
#endif


struct d_single_shooting_ocp {
  struct d_nmpc_model model;
  struct d_time_varying_smooth_horizon horizon;
  double *dx;
  double **x_sequence;
  double **lmd_sequence;
  int dimx;
  int dimu;
  int dimc;
  int dimuc;
  int N;
  int dim_solution;
  int memsize; // memory size in bytes
};

int d_single_shooting_ocp_strsize();

int d_single_shooting_ocp_memsize(int N);

void d_single_shooting_ocp_create(struct d_single_shooting_ocp *ocp, double T_f, 
                                  double alpha, double initial_time, int N);

void d_single_shooting_ocp_delete(struct d_single_shooting_ocp *ocp);

void d_single_shooting_ocp_compute_optimality_residual(
    struct d_single_shooting_ocp *ocp, double current_time, 
    double *current_state, double *solution, double *optimality_residual);

void d_single_shooting_ocp_predict_state_from_solution(
    struct d_single_shooting_ocp *ocp, double current_time, 
    double *current_state, double *solution, double prediction_length, 
    double *predicted_state);

void d_single_shooting_ocp_reset_horizon_length(
    struct d_single_shooting_ocp *ocp, double initial_time);

int d_single_shooting_ocp_dimx();

int d_single_shooting_ocp_dimu();

int d_single_shooting_ocp_dimc();


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // CGMRES_D_SINGLE_SHOOTING_OCP_H_