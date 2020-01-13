#ifndef CGMRES_S_INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_MFGMRES_ARGS_H_
#define CGMRES_S_INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_MFGMRES_ARGS_H_


#ifdef __cplusplus
extern "C" {
#endif


struct s_inexact_newton_for_zero_horizon_ocp_mfgmres_args {
  float current_time;
  float *current_state_ptr;
  float *current_solution_ptr;
};


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // CGMRES_S_INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_MFGMRES_ARGS_H_