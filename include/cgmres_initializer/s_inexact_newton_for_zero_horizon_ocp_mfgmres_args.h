#ifndef HPCGMRES_S_INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_MFGMRES_ARGS_H_
#define HPCGMRES_S_INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_MFGMRES_ARGS_H_


#include <blasfeo.h>


#ifdef __cplusplus
extern "C" {
#endif


struct s_inexact_newton_for_zero_horizon_ocp_mfgmres_args {
  float current_time;
  struct blasfeo_svec *current_state_ptr;
  struct blasfeo_svec *current_solution_ptr;
};


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // HPCGMRES_S_INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_MFGMRES_ARGS_H_