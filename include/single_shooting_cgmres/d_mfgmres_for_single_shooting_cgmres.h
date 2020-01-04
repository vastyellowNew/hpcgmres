#ifndef HPCGMRES_D_MFGMRES_FOR_SINGLE_SHOOTING_CGMRES_H_
#define HPCGMRES_D_MFGMRES_FOR_SINGLE_SHOOTING_CGMRES_H_


#include <blasfeo.h>

#include "d_single_shooting_continuation.h"
#include "d_single_shooting_continuation_mfgmres_args.h"


#ifdef __cplusplus
extern "C" {
#endif


struct d_mfgmres_for_single_shooting_cgmres {
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
int d_mfgmres_for_single_shooting_cgmres_strsize();

// returns the size of the memory
int d_mfgmres_for_single_shooting_cgmres_memsize(int dim_linear_problem, 
                                                 int kmax);

// constructor
void d_mfgmres_for_single_shooting_cgmres_create(
    struct d_mfgmres_for_single_shooting_cgmres *mfgmres,
    int dim_linear_problem, int kmax, void *memory);

void d_mfgmres_for_single_shooting_cgmres_solve_linear_problem(
    struct d_mfgmres_for_single_shooting_cgmres *mfgmres,
    struct d_single_shooting_continuation *continuation,
    struct d_single_shooting_continuation_mfgmres_args *linear_problem_args,
    struct blasfeo_dvec *solution_vec);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // HPCGMRES_D_MFGMRES_FOR_SINGLE_SHOOTING_CGMRES_H_