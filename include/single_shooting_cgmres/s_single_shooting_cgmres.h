#ifndef CGMRES_S_SINGLE_SHOOTING_CGMRES_H_
#define CGMRES_S_SINGLE_SHOOTING_CGMRES_H_


#include "single_shooting_cgmres/s_mfgmres_for_single_shooting_cgmres.h"
#include "single_shooting_cgmres/s_single_shooting_continuation.h"
#include "single_shooting_cgmres/s_single_shooting_continuation_mfgmres_args.h"
#include "cgmres_initializer/s_cgmres_initializer.h"


#ifdef __cplusplus
extern "C" {
#endif


struct s_single_shooting_cgmres {
  struct s_mfgmres_for_single_shooting_cgmres mfgmres;
  struct s_single_shooting_continuation continuation;
  struct s_single_shooting_continuation_mfgmres_args mfgmres_args;
  struct s_cgmres_initializer initializer;
  float *solution;
  float *solution_update;
  float *initial_solution;
  int dimu;
  int dimc;
  int dimuc;
  int N;
  int memsize; // memory size in bytes
};


int s_single_shooting_cgmres_strsize();

int s_single_shooting_cgmres_memsize(int N, int kmax);

void s_single_shooting_cgmres_create(struct s_single_shooting_cgmres *cgmres, 
                                     float T_f, float alpha, 
                                     float initial_time, int N,
                                     float finite_difference_increment,
                                     float zeta, int kmax);

void s_single_shooting_cgmres_delete(struct s_single_shooting_cgmres *cgmres);

void s_single_shooting_cgmres_set_initialization_parameters(
    struct s_single_shooting_cgmres *cgmres, float newton_residual_tolerance, 
    int max_newton_iteration, float *initial_guess_solution);

void s_single_shooting_cgmres_initialize_solution(
    struct s_single_shooting_cgmres *cgmres, float initial_time, 
    float *initial_state);

void s_single_shooting_cgmres_update_control_input(
    struct s_single_shooting_cgmres *cgmres, float current_time,
    float *current_state, float sampling_period, float *control_input);

void s_single_shooting_cgmres_get_control_input(
    struct s_single_shooting_cgmres *cgmres, float *control_input);

float s_single_shooting_cgmres_get_error_norm(
    struct s_single_shooting_cgmres *cgmres, float current_time,
    float *current_state);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // CGMRES_S_SINGLE_SHOOTING_CGMRES_H_