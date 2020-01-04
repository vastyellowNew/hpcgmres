#ifndef CGMRES_D_SINGLE_SHOOTING_OCP_H_
#define CGMRES_D_SINGLE_SHOOTING_OCP_H_


#include <blasfeo.h>

#include "d_nmpc_model.h"
#include "d_time_varying_smooth_horizon.h"


#ifdef __cplusplus
extern "C" {
#endif


struct d_single_shooting_ocp {
  struct d_nmpc_model *model;
  struct d_time_varying_smooth_horizon *horizon;
  struct blasfeo_dvec *dx_vec;
  struct blasfeo_dvec *x_sequence_vec;
  struct blasfeo_dvec *lmd_sequence_vec;
  int dimx;
  int dimu;
  int dimc;
  int dimuc;
  int N;
  int dim_solution;
  int memsize; // memory size in bytes
};

int d_single_shooting_ocp_strsize();

int d_single_shooting_ocp_memsize(struct d_single_shooting_ocp *ocp, int N);


void d_single_shooting_ocp_create(struct d_single_shooting_ocp *ocp,
                                  double T_f, double alpha, double initial_time, 
                                  int N, void *memory);

void d_single_shooting_ocp_compute_optimality_residual(
    struct d_single_shooting_ocp *ocp, double current_time, 
    struct blasfeo_dvec *current_state, struct blasfeo_dvec *solution, 
    struct blasfeo_dvec *optimality_residual);

void d_single_shooting_ocp_predict_state_from_solution(
    struct d_single_shooting_ocp *ocp, double current_time, 
    struct blasfeo_dvec *current_state, struct blasfeo_dvec *solution, 
    double prediction_length, struct blasfeo_dvec *predicted_state);

void d_single_shooting_ocp_reset_horizon_length(
    struct d_single_shooting_ocp *ocp, double initial_time);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // CGMRES_D_SINGLE_SHOOTING_OCP_H_