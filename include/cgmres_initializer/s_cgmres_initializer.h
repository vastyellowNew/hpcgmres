#ifndef HPCGMRES_S_CGMRES_INITIALIZER_H_
#define HPCGMRES_S_CGMRES_INITIALIZER_H_


#include <blasfeo.h>

#include "s_inexact_newton_for_zero_horizon_ocp_mfgmres_args.h"
#include "s_inexact_newton_for_zero_horizon_ocp.h"
#include "s_mfgmres_for_cgmres_initializer.h"


#ifdef __cplusplus
extern "C" {
#endif


struct s_cgmres_initializer {
  struct s_inexact_newton_for_zero_horizon_ocp *newton;
  struct s_inexact_newton_for_zero_horizon_ocp_mfgmres_args *mfgmres_args;
  struct s_mfgmres_for_cgmres_initializer *mfgmres;
  struct blasfeo_svec *initial_guess_solution; 
  struct blasfeo_svec *solution_update;
  double newton_residual_tolerance;
  int max_newton_iteration;
  int dim_solution;
  int memsize; // memory size in bytes
};

int s_cgmres_initializer_strsize();

int s_cgmres_initializer_memsize(struct s_cgmres_initializer *initializer);

void s_cgmres_initializer_create(struct s_cgmres_initializer *initializer, 
                                 float finite_difference_increment, int kmax,
                                 void *memory);

void s_cgmres_initializer_set_termination_criterions(
    struct s_cgmres_initializer *initializer, 
    float newton_residual_tolerance, int max_newton_iteration);

void s_cgmres_initializer_set_initial_guess_solution(
    struct s_cgmres_initializer *initializer, 
    struct blasfeo_svec *initial_guess_solution);

void s_cgmres_initializer_compute_initial_solution(
    struct s_cgmres_initializer *initializer, float initial_time, 
    struct blasfeo_svec *initial_state, struct blasfeo_svec *initial_solution);

void s_cgmres_initializer_get_terminal_cost_derivative(
    struct s_cgmres_initializer *initializer, float initial_time, 
    struct blasfeo_svec *initial_state, 
    struct blasfeo_svec *terminal_cost_derivative);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // HPCGMRES_S_CGMRES_INITIALIZER_H_