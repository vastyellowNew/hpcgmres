#ifndef HPCGMRES_D_MFGMRES_FOR_CGMRES_INITIALIZER_H_
#define HPCGMRES_D_MFGMRES_FOR_CGMRES_INITIALIZER_H_


#include "d_inexact_newton_for_zero_horizon_ocp.h"
#include "d_inexact_newton_for_zero_horizon_ocp_mfgmres_args.h"


#ifdef __cplusplus
extern "C" {
#endif


struct d_mfgmres_for_cgmres_initializer {
  double **hessenberg_mat; // hessenberg matrix 
  double **basis_mat; // basis of the Krylov subspace
  double *b_vec; // strage of the residual
  double *givens_c_vec; // vector for Givens rotation
  double *givens_s_vec; // vector for Givens rotation
  double *g_vec; // temporal strage
  int dim_linear_problem; // dimension of the linear problem
  int kmax; // maximum dimension of the Krylov subspace
  int memsize; // memory size in bytes
};


int d_mfgmres_for_cgmres_initializer_strsize();

int d_mfgmres_for_cgmres_initializer_memsize(int dim_linear_problem, int kmax);

void d_mfgmres_for_cgmres_initializer_create(
    struct d_mfgmres_for_cgmres_initializer *mfgmres, int dim_linear_problem, 
    int kmax);

void d_mfgmres_for_cgmres_initializer_delete(
    struct d_mfgmres_for_cgmres_initializer *mfgmres);

void d_mfgmres_for_cgmres_initializer_solve_linear_problem(
    struct d_mfgmres_for_cgmres_initializer *mfgmres, 
    struct d_inexact_newton_for_zero_horizon_ocp *newton, 
    struct d_inexact_newton_for_zero_horizon_ocp_mfgmres_args *linear_problem_args, 
    double *solution_vec);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // HPCGMRES_D_MFGMRES_FOR_CGMRES_INITIALIZER_H_