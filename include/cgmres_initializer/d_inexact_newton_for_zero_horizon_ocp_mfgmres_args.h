#ifndef CGMRES_D_INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_MFGMRES_ARGS_H_
#define CGMRES_D_INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_MFGMRES_ARGS_H_


#ifdef __cplusplus
extern "C" {
#endif


struct d_inexact_newton_for_zero_horizon_ocp_mfgmres_args {
  double current_time;
  double *current_state_ptr;
  double *current_solution_ptr;
};


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // CGMRES_D_INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_MFGMRES_ARGS_H_