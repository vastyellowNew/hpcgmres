#ifndef HPCGMRES_S_MFGMRES_FOR_SINGLE_SHOOTING_CGMRES_H_
#define HPCGMRES_S_MFGMRES_FOR_SINGLE_SHOOTING_CGMRES_H_


#include "s_single_shooting_continuation.h"
#include "s_single_shooting_continuation_mfgmres_args.h"


#ifdef __cplusplus
extern "C" {
#endif


struct s_mfgmres_for_single_shooting_cgmres {
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


int s_mfgmres_for_single_shooting_cgmres_strsize();

int s_mfgmres_for_single_shooting_cgmres_memsize(int dim_linear_problem, 
                                                 int kmax);

void s_mfgmres_for_single_shooting_cgmres_create(
    struct s_mfgmres_for_single_shooting_cgmres *mfgmres,
    int dim_linear_problem, int kmax);

void s_mfgmres_for_single_shooting_cgmres_delete(
    struct s_mfgmres_for_single_shooting_cgmres *mfgmres);

void s_mfgmres_for_single_shooting_cgmres_solve_linear_problem(
    struct s_mfgmres_for_single_shooting_cgmres *mfgmres,
    struct s_single_shooting_continuation *continuation,
    struct s_single_shooting_continuation_mfgmres_args *linear_problem_args,
    float *solution_vec);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // HPCGMRES_S_MFGMRES_FOR_SINGLE_SHOOTING_CGMRES_H_