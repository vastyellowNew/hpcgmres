#ifndef HPCGMRES_D_SINGLE_SHOOTING_CGMRES_H_
#define HPCGMRES_D_SINGLE_SHOOTING_CGMRES_H_


#include <blasfeo.h>

#include "d_single_shooting_continuation.h"
#include "d_single_shooting_continuation_mfgmres_args.h"
#include "d_mfgmres_for_single_shooting_cgmres.h"
#include "d_cgmres_initializer.h"


#ifdef __cplusplus
extern "C" {
#endif


struct d_single_shooting_cgmres {
  struct d_single_shooting_continuation *continuation;
  struct d_single_shooting_continuation_mfgmres_args *mfgmres_args;
  struct d_mfgmres_for_single_shooting_cgmres *mfgmres_args;
  struct d_cgmres_initializer *initializer;
  struct blasfeo_dvec *solution;
  struct blasfeo_dvec *solution_update;
  struct blasfeo_dvec *initial_solution;
  int dimu;
  int dimc;
  int dimuc;
  int N;
  int memsize; // memory size in bytes
};

int d_single_shooting_cgmres_strsize();

int d_single_shooting_cgmres_memsize(struct d_single_shooting_cgmres *cgmres,
                                     int N, int kmax);

void d_single_shooting_cgmres_create(struct d_single_shooting_cgmres *cgmres, 
                                     double T_f, double alpha, 
                                     double initial_time, int N, 
                                     double finite_difference_increment,
                                     double zeta, int kmax);

void d_single_shooting_cgmres_set_initialization_parameters(
    struct d_single_shooting_cgmres *cgmres, 
    struct blasfeo_dvec *initial_guess_solution,
    double newton_residual_tolerance, int max_newton_iteration);

void d_single_shooting_cgmres_initialize_solution(
    struct d_single_shooting_cgmres *cgmres, double initial_time, 
    struct blasfeo_dvec *initial_state);

void d_single_shooting_cgmres_update_control_input(
    struct d_single_shooting_cgmres *cgmres, double current_time,
    struct blasfeo_dvec *current_state, double sampling_period,
    struct blasfeo_dvec *control_input);

void d_single_shooting_cgmres_get_control_input(
    struct d_single_shooting_cgmres *cgmres, 
    struct blasfeo_dvec *control_input);

double d_single_shooting_cgmres_get_error_norm(
    struct d_single_shooting_cgmres *cgmres, double current_time,
    struct blasfeo_dvec *current_state);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // HPCGMRES_D_SINGLE_SHOOTING_CGMRES_H_