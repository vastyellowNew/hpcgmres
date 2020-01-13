#ifndef CGMRES_S_SINGLE_SHOOTING_CONTINUATION_H_
#define CGMRES_S_SINGLE_SHOOTING_CONTINUATION_H_


#include "single_shooting_cgmres/s_single_shooting_ocp.h"
#include "single_shooting_cgmres/s_single_shooting_continuation_mfgmres_args.h"


#ifdef __cplusplus
extern "C" {
#endif


struct s_single_shooting_continuation {
  struct s_single_shooting_ocp ocp;
  float *incremented_state;
  float *incremented_solution;
  float *optimality_residual;
  float *optimality_residual1;
  float *optimality_residual2;
  float finite_difference_increment; 
  float zeta;
  float incremented_time;
  int dimx;
  int dimu;
  int dimc;
  int N;
  int dim_solution;
  int memsize; // memory size in bytes
};

int s_single_shooting_continuation_strsize();

int s_single_shooting_continuation_memsize(int N);

void s_single_shooting_continuation_create(
    struct s_single_shooting_continuation *continuation, float T_f, 
    float alpha, float initial_time, int N, 
    float finite_difference_increment, float zeta);

void s_single_shooting_continuation_delete(
    struct s_single_shooting_continuation *continuation);

void s_single_shooting_continuation_integrate_solution(
    struct s_single_shooting_continuation *continuation, float *solution,
    float *solution_update, float integration_length);

void s_single_shooting_continuation_compute_b(
    struct s_single_shooting_continuation *continuation, 
    struct s_single_shooting_continuation_mfgmres_args *args, 
    float *direction, float *b);

void s_single_shooting_continuation_compute_ax(
    struct s_single_shooting_continuation *continuation, 
    struct s_single_shooting_continuation_mfgmres_args *args, 
    float *direction, float *ax);

float s_single_shooting_continuation_get_error_norm(
    struct s_single_shooting_continuation *continuation, float current_time, 
    float *current_state, float *current_solution);

void s_single_shooting_continuation_reset_horizon_length(
    struct s_single_shooting_continuation *continuation, float initial_time);

int s_single_shooting_continuation_dimx();

int s_single_shooting_continuation_dimu();

int s_single_shooting_continuation_dimc();


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // CGMRES_S_SINGLE_SHOOTING_CONTINUATION_H_