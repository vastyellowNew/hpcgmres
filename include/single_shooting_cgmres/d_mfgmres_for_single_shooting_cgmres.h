#ifndef HPCGMRES_D_MFGMRES_FOR_SINGLE_SHOOTING_CGMRES_H_
#define HPCGMRES_D_MFGMRES_FOR_SINGLE_SHOOTING_CGMRES_H_


#include "single_shooting_cgmres/d_single_shooting_continuation.h"
#include "single_shooting_cgmres/d_single_shooting_continuation_mfgmres_args.h"


#ifdef __cplusplus
extern "C" {
#endif


struct d_mfgmres_for_single_shooting_cgmres {
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


int d_mfgmres_for_single_shooting_cgmres_strsize();

int d_mfgmres_for_single_shooting_cgmres_memsize(int dim_linear_problem, 
                                                 int kmax);

void d_mfgmres_for_single_shooting_cgmres_create(
    struct d_mfgmres_for_single_shooting_cgmres *mfgmres,
    int dim_linear_problem, int kmax);

void d_mfgmres_for_single_shooting_cgmres_delete(
    struct d_mfgmres_for_single_shooting_cgmres *mfgmres);

void d_mfgmres_for_single_shooting_cgmres_solve_linear_problem(
    struct d_mfgmres_for_single_shooting_cgmres *mfgmres,
    struct d_single_shooting_continuation *continuation,
    struct d_single_shooting_continuation_mfgmres_args *linear_problem_args,
    double *solution_vec);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // HPCGMRES_D_MFGMRES_FOR_SINGLE_SHOOTING_CGMRES_H_