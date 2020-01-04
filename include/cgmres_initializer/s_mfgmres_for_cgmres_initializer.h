#ifndef HPCGMRES_S_MFGMRES_FOR_CGMRES_INITIALIZER_H_
#define HPCGMRES_S_MFGMRES_FOR_CGMRES_INITIALIZER_H_


#include <blasfeo.h>

#include "s_inexact_newton_for_zero_horizon_ocp.h"
#include "s_inexact_newton_for_zero_horizon_ocp_mfgmres_args.h"


#ifdef __cplusplus
extern "C" {
#endif


struct s_mfgmres_for_cgmres_initializer {
  struct blasfeo_smat *hessenberg_mat; // hessenberg matrix 
  struct blasfeo_svec *basis_vec; // basis of the Krylov subspace
  struct blasfeo_svec *b_vec; // strage of the residual
  struct blasfeo_svec *givens_c_vec; // vector for Givens rotation
  struct blasfeo_svec *givens_s_vec; // vector for Givens rotation
  struct blasfeo_svec *g_vec; // temporal strage
  int dim_linear_problem; // dimension of the linear problem
  int kmax; // maximum dimension of the Krylov subspace
  int memsize; // memory size in bytes
};

// returns the size of the strage
int s_mfgmres_for_cgmres_initializer_strsize();

// returns the size of the memory
int s_mfgmres_for_cgmres_initializer_memsize(int dim_linear_problem, int kmax);

// constructor
void s_mfgmres_for_cgmres_initializer_create(
    struct s_mfgmres_for_cgmres_initializer *mfgmres, int dim_linear_problem, 
    int kmax, void *memory);

void s_mfgmres_for_cgmres_initializer_solve_linear_problem(
    struct s_mfgmres_for_cgmres_initializer *mfgmres, 
    struct s_inexact_newton_for_zero_horizon_ocp *newton, 
    struct s_inexact_newton_for_zero_horizon_ocp_mfgmres_args *linear_problem_args, 
    struct blasfeo_svec *solution_vec);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // HPCGMRES_S_MFGMRES_FOR_CGMRES_INITIALIZER_H_