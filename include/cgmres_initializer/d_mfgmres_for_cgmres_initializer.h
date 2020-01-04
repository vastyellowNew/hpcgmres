#ifndef HPCGMRES_D_MFGMRES_FOR_CGMRES_INITIALIZER_H_
#define HPCGMRES_D_MFGMRES_FOR_CGMRES_INITIALIZER_H_


#include <blasfeo.h>

#include "d_inexact_newton_for_zero_horizon_ocp.h"
#include "d_inexact_newton_for_zero_horizon_ocp_mfgmres_args.h"


#ifdef __cplusplus
extern "C" {
#endif


struct d_mfgmres_for_cgmres_initializer {
  struct blasfeo_dmat *hessenberg_mat; // hessenberg matrix 
  struct blasfeo_dvec *basis_vec; // basis of the Krylov subspace
  struct blasfeo_dvec *b_vec; // strage of the residual
  struct blasfeo_dvec *givens_c_vec; // vector for Givens rotation
  struct blasfeo_dvec *givens_s_vec; // vector for Givens rotation
  struct blasfeo_dvec *g_vec; // temporal strage
  int dim_linear_problem; // dimension of the linear problem
  int kmax; // maximum dimension of the Krylov subspace
  int memsize; // memory size in bytes
};

// returns the size of the strage
int d_mfgmres_for_cgmres_initializer_strsize();

// returns the size of the memory
int d_mfgmres_for_cgmres_initializer_memsize(int dim_linear_problem, int kmax);

// constructor
void d_mfgmres_for_cgmres_initializer_create(
    struct d_mfgmres_for_cgmres_initializer *mfgmres, int dim_linear_problem, 
    int kmax, void *memory);

void d_mfgmres_for_cgmres_initializer_solve_linear_problem(
    struct d_mfgmres_for_cgmres_initializer *mfgmres, 
    struct d_inexact_newton_for_zero_horizon_ocp *newton, 
    struct d_inexact_newton_for_zero_horizon_ocp_mfgmres_args *linear_problem_args, 
    struct blasfeo_dvec *solution_vec);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // HPCGMRES_D_MFGMRES_FOR_CGMRES_INITIALIZER_H_