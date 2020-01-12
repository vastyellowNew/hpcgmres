#ifndef HPCGMRES_S_CGMRES_INITIALIZER_H_
#define HPCGMRES_S_CGMRES_INITIALIZER_H_


#include "cgmres_initializer/s_mfgmres_for_cgmres_initializer.h"
#include "cgmres_initializer/s_inexact_newton_for_zero_horizon_ocp.h"
#include "cgmres_initializer/s_inexact_newton_for_zero_horizon_ocp_mfgmres_args.h"


#ifdef __cplusplus
extern "C" {
#endif


struct s_cgmres_initializer {
  struct s_mfgmres_for_cgmres_initializer mfgmres;
  struct s_inexact_newton_for_zero_horizon_ocp newton;
  struct s_inexact_newton_for_zero_horizon_ocp_mfgmres_args mfgmres_args;
  float *initial_guess_solution; 
  float *solution_update;
  float newton_residual_tolerance;
  int max_newton_iteration;
  int dim_solution;
  int memsize; // memory size in bytes
};

int s_cgmres_initializer_strsize();

int s_cgmres_initializer_memsize(int kmax);

void s_cgmres_initializer_create(struct s_cgmres_initializer *initializer, 
                                 float finite_difference_increment, int kmax);

void s_cgmres_initializer_delete(struct s_cgmres_initializer *initializer);

void s_cgmres_initializer_set_termination_criterions(
    struct s_cgmres_initializer *initializer, 
    float newton_residual_tolerance, int max_newton_iteration);

void s_cgmres_initializer_set_initial_guess_solution(
    struct s_cgmres_initializer *initializer, float *initial_guess_solution);

void s_cgmres_initializer_compute_initial_solution(
    struct s_cgmres_initializer *initializer, float initial_time, 
    float *initial_state, float *initial_solution);

void s_cgmres_initializer_get_terminal_cost_derivative(
    struct s_cgmres_initializer *initializer, float initial_time, 
    float *initial_state, float *terminal_cost_derivative);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // HPCGMRES_S_CGMRES_INITIALIZER_H_