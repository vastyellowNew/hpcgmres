#ifndef HPCGMRES_D_SINGLE_SHOOTING_CONTINUATION_MFGMRES_ARGS_H_
#define HPCGMRES_D_SINGLE_SHOOTING_CONTINUATION_MFGMRES_ARGS_H_


#include <blasfeo.h>


#ifdef __cplusplus
extern "C" {
#endif


struct d_single_shooting_continuation_mfgmres_args {
  double current_time;
  struct blasfeo_dvec *current_state_ptr;
  struct blasfeo_dvec *current_solution_ptr;
};


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // HPCGMRES_D_SINGLE_SHOOTING_CONTINUATION_MFGMRES_ARGS_H_