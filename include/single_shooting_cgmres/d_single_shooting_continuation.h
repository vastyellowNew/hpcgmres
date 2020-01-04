#ifndef HPCGMRES_D_SINGLE_SHOOTING_CONTINUATION_H_
#define HPCGMRES_D_SINGLE_SHOOTING_CONTINUATION_H_


#include <blasfeo.h>

#include "d_single_shooting_ocp.h"
#include "d_single_shooting_continuation_mfgmres_args.h"


#ifdef __cplusplus
extern "C" {
#endif


struct d_single_shooting_continuation {
  struct d_single_shooting_ocp *ocp;
  struct blasfeo_dvec *incremented_state;
  struct blasfeo_dvec *incremented_solution;
  struct blasfeo_dvec *optimality_residual;
  struct blasfeo_dvec *optimality_residual1;
  struct blasfeo_dvec *optimality_residual2;
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

int d_single_shooting_continuation_memsize(
    struct d_single_shooting_continuation *continuation, int N);

void d_single_shooting_continuation_create(
    struct d_single_shooting_continuation *continuation, double T_f, 
    double alpha, double initial_time, int N, 
    double finite_difference_increment, double zeta, void *memory);

void d_single_shooting_continuation_integrate_solution(
    struct d_single_shooting_ocp *ocp, struct blasfeo_dvec *solution,
    struct blasfeo_dvec *solution_update_vec, double integration_length);

void d_single_shooting_continuation_compute_b(
    struct d_single_shooting_continuation *continuation, 
    struct d_single_shooting_continuation_mfgmres_args *args, 
    struct blasfeo_dvec *direction, struct blasfeo_dvec *b);

void d_single_shooting_continuation_compute_ax(
    struct d_single_shooting_continuation *continuation, 
    struct d_single_shooting_continuation_mfgmres_args *args, 
    struct blasfeo_dvec *direction, struct blasfeo_dvec *ax);

double d_single_shooting_continuation_get_error_norm(
    struct d_single_shooting_continuation *continuation, double current_time, 
    struct blasfeo_dvec *current_state, struct blasfeo_dvec *current_solution);

void d_single_shooting_continuation_reset_horizon_length(
    struct d_single_shooting_continuation *continuation, double initial_time);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // HPCGMRES_D_SINGLE_SHOOTING_CONTINUATION_H_