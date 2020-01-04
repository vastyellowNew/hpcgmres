#ifndef HPCGMRES_D_CGMRES_INITIALIZER_H_
#define HPCGMRES_D_CGMRES_INITIALIZER_H_


#include <blasfeo.h>

#include "d_inexact_newton_for_zero_horizon_ocp.h"
#include "d_inexact_newton_for_zero_horizon_ocp_mfgmres_args.h"
#include "d_mfgmres_for_cgmres_initializer.h"


#ifdef __cplusplus
extern "C" {
#endif


struct d_cgmres_initializer {
  struct d_inexact_newton_for_zero_horizon_ocp *newton;
  struct d_inexact_newton_for_zero_horizon_ocp_mfgmres_args *mfgmres_args;
  struct d_mfgmres_for_cgmres_initializer *mfgmres;
  struct blasfeo_dvec *initial_guess_solution; 
  struct blasfeo_dvec *solution_update;
  double newton_residual_tolerance;
  int max_newton_iteration;
  int dim_solution;
  int memsize; // memory size in bytes
};

int d_cgmres_initializer_strsize();

int d_cgmres_initializer_memsize(struct d_cgmres_initializer *initializer);

void d_cgmres_initializer_create(struct d_cgmres_initializer *initializer, 
                                 double finite_difference_increment, int kmax,
                                 void *memory);

void d_cgmres_initializer_set_termination_criterions(
    struct d_cgmres_initializer *initializer, 
    double newton_residual_tolerance, int max_newton_iteration);

void d_cgmres_initializer_set_initial_guess_solution(
    struct d_cgmres_initializer *initializer, 
    struct blasfeo_dvec *initial_guess_solution);

void d_cgmres_initializer_compute_initial_solution(
    struct d_cgmres_initializer *initializer, double initial_time, 
    struct blasfeo_dvec *initial_state, struct blasfeo_dvec *initial_solution);

void d_cgmres_initializer_get_terminal_cost_derivative(
    struct d_cgmres_initializer *initializer, double initial_time, 
    struct blasfeo_dvec *initial_state, 
    struct blasfeo_dvec *terminal_cost_derivative);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // HPCGMRES_D_CGMRES_INITIALIZER_H_