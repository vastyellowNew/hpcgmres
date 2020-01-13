#ifndef CGMRES_D_SINGLE_SHOOTING_CONTINUATION_H_
#define CGMRES_D_SINGLE_SHOOTING_CONTINUATION_H_


#include "single_shooting_cgmres/d_single_shooting_ocp.h"
#include "single_shooting_cgmres/d_single_shooting_continuation_mfgmres_args.h"


#ifdef __cplusplus
extern "C" {
#endif


struct d_single_shooting_continuation {
  struct d_single_shooting_ocp ocp;
  double *incremented_state;
  double *incremented_solution;
  double *optimality_residual;
  double *optimality_residual1;
  double *optimality_residual2;
  double finite_difference_increment; 
  double zeta;
  double incremented_time;
  int dimx;
  int dimu;
  int dimc;
  int N;
  int dim_solution;
  int memsize; // memory size in bytes
};

int d_single_shooting_continuation_strsize();

int d_single_shooting_continuation_memsize(int N);

void d_single_shooting_continuation_create(
    struct d_single_shooting_continuation *continuation, double T_f, 
    double alpha, double initial_time, int N, 
    double finite_difference_increment, double zeta);

void d_single_shooting_continuation_delete(
    struct d_single_shooting_continuation *continuation);

void d_single_shooting_continuation_integrate_solution(
    struct d_single_shooting_continuation *continuation, double *solution,
    double *solution_update, double integration_length);

void d_single_shooting_continuation_compute_b(
    struct d_single_shooting_continuation *continuation, 
    struct d_single_shooting_continuation_mfgmres_args *args, 
    double *direction, double *b);

void d_single_shooting_continuation_compute_ax(
    struct d_single_shooting_continuation *continuation, 
    struct d_single_shooting_continuation_mfgmres_args *args, 
    double *direction, double *ax);

double d_single_shooting_continuation_get_error_norm(
    struct d_single_shooting_continuation *continuation, double current_time, 
    double *current_state, double *current_solution);

void d_single_shooting_continuation_reset_horizon_length(
    struct d_single_shooting_continuation *continuation, double initial_time);

int d_single_shooting_continuation_dimx();

int d_single_shooting_continuation_dimu();

int d_single_shooting_continuation_dimc();


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // CGMRES_D_SINGLE_SHOOTING_CONTINUATION_H_