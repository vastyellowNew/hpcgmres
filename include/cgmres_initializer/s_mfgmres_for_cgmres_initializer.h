#ifndef HPCGMRES_S_MFGMRES_FOR_CGMRES_INITIALIZER_H_
#define HPCGMRES_S_MFGMRES_FOR_CGMRES_INITIALIZER_H_


#include "cgmres_initializer/s_inexact_newton_for_zero_horizon_ocp.h"
#include "cgmres_initializer/s_inexact_newton_for_zero_horizon_ocp_mfgmres_args.h"


#ifdef __cplusplus
extern "C" {
#endif


struct s_mfgmres_for_cgmres_initializer {
  float **hessenberg_mat; // hessenberg matrix 
  float **basis_mat; // basis of the Krylov subspace
  float *b_vec; // strage of the residual
  float *givens_c_vec; // vector for Givens rotation
  float *givens_s_vec; // vector for Givens rotation
  float *g_vec; // temporal strage
  int dim_linear_problem; // dimension of the linear problem
  int kmax; // maximum dimension of the Krylov subspace
  int memsize; // memory size in bytes
};

int s_mfgmres_for_cgmres_initializer_strsize();

int s_mfgmres_for_cgmres_initializer_memsize(int dim_linear_problem, int kmax);

void s_mfgmres_for_cgmres_initializer_create(
    struct s_mfgmres_for_cgmres_initializer *mfgmres, int dim_linear_problem, 
    int kmax);

void s_mfgmres_for_cgmres_initializer_delete(
    struct s_mfgmres_for_cgmres_initializer *mfgmres);

void s_mfgmres_for_cgmres_initializer_solve_linear_problem(
    struct s_mfgmres_for_cgmres_initializer *mfgmres, 
    struct s_inexact_newton_for_zero_horizon_ocp *newton, 
    struct s_inexact_newton_for_zero_horizon_ocp_mfgmres_args *linear_problem_args, 
    float *solution_vec);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // HPCGMRES_S_MFGMRES_FOR_CGMRES_INITIALIZER_H_