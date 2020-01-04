#ifndef HPCGMRES_S_SINGLE_SHOOTING_CGMRES_H_
#define HPCGMRES_S_SINGLE_SHOOTING_CGMRES_H_


#include <blasfeo.h>

#include "s_single_shooting_continuation.h"
#include "s_single_shooting_continuation_mfgmres_args.h"
#include "s_mfgmres_for_single_shooting_cgmres.h"
#include "s_cgmres_initializer.h"


#ifdef __cplusplus
extern "C" {
#endif


struct s_single_shooting_cgmres {
  struct s_single_shooting_continuation *continuation;
  struct s_single_shooting_continuation_mfgmres_args *mfgmres_args;
  struct s_mfgmres_for_single_shooting_cgmres *mfgmres_args;
  struct s_cgmres_initializer *initializer;
  struct blasfeo_svec *solution;
  struct blasfeo_svec *solution_update;
  struct blasfeo_svec *initial_solution;
  int dimu;
  int dimc;
  int dimuc;
  int N;
  int memsize; // memory size in bytes
};

int s_single_shooting_cgmres_strsize();

int s_single_shooting_cgmres_memsize(struct d_single_shooting_cgmres *cgmres,
                                     int N, int kmax);

void s_single_shooting_cgmres_create(struct d_single_shooting_cgmres *cgmres, 
                                     float T_f, float alpha, float initial_time, 
                                     int N, float finite_difference_increment,
                                     float zeta, int kmax);

void s_single_shooting_cgmres_set_initialization_parameters(
    struct s_single_shooting_cgmres *cgmres, 
    struct blasfeo_svec *initial_guess_solution,
    float newton_residual_tolerance, int max_newton_iteration);

void s_single_shooting_cgmres_initialize_solution(
    struct s_single_shooting_cgmres *cgmres, float initial_time, 
    struct blasfeo_svec *initial_state);

void s_single_shooting_cgmres_update_control_input(
    struct s_single_shooting_cgmres *cgmres, float current_time,
    struct blasfeo_svec *current_state, float sampling_period,
    struct blasfeo_svec *control_input);

void s_single_shooting_cgmres_get_control_input(
    struct s_single_shooting_cgmres *cgmres, 
    struct blasfeo_svec *control_input);

float s_single_shooting_cgmres_get_error_norm(
    struct s_single_shooting_cgmres *cgmres, float current_time,
    struct blasfeo_svec *current_state);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // HPCGMRES_S_SINGLE_SHOOTING_CGMRES_H_