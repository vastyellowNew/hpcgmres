#ifndef HPCGMRES_S_SINGLE_SHOOTING_CONTINUATION_H_
#define HPCGMRES_S_SINGLE_SHOOTING_CONTINUATION_H_


#include <blasfeo.h>

#include "s_single_shooting_ocp.h"
#include "s_single_shooting_continuation_mfgmres_args.h"


#ifdef __cplusplus
extern "C" {
#endif


struct s_single_shooting_continuation {
  struct s_single_shooting_ocp *ocp;
  struct blasfeo_svec *incremented_state;
  struct blasfeo_svec *incremented_solution;
  struct blasfeo_svec *optimality_residual;
  struct blasfeo_svec *optimality_residual1;
  struct blasfeo_svec *optimality_residual2;
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

int s_single_shooting_continuation_memsize(
    struct s_single_shooting_continuation *continuation, int N);

void s_single_shooting_continuation_create(
    struct s_single_shooting_continuation *continuation, float T_f, 
    float alpha, float initial_time, int N, 
    float finite_difference_increment, float zeta, void *memory);

void s_single_shooting_continuation_integrate_solution(
    struct s_single_shooting_ocp *ocp, struct blasfeo_svec *solution,
    struct blasfeo_svec *solution_update_vec, float integration_length);

void s_single_shooting_continuation_compute_b(
    struct s_single_shooting_continuation *continuation, 
    struct s_single_shooting_continuation_mfgmres_args *args, 
    struct blasfeo_svec *direction, struct blasfeo_svec *b);

void s_single_shooting_continuation_compute_ax(
    struct s_single_shooting_continuation *continuation, 
    struct s_single_shooting_continuation_mfgmres_args *args, 
    struct blasfeo_svec *direction, struct blasfeo_svec *ax);

float s_single_shooting_continuation_get_error_norm(
    struct s_single_shooting_continuation *continuation, float current_time, 
    struct blasfeo_svec *current_state, struct blasfeo_svec *current_solution);

void s_single_shooting_continuation_reset_horizon_length(
    struct s_single_shooting_continuation *continuation, float initial_time);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // HPCGMRES_S_SINGLE_SHOOTING_CONTINUATION_H_