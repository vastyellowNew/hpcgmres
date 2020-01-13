#ifndef CGMRES_D_SINGLE_SHOOTING_CONTINUATION_MFGMRES_ARGS_H_
#define CGMRES_D_SINGLE_SHOOTING_CONTINUATION_MFGMRES_ARGS_H_


#ifdef __cplusplus
extern "C" {
#endif


struct d_single_shooting_continuation_mfgmres_args {
  double current_time;
  double *current_state_ptr;
  double *current_solution_ptr;
};


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // CGMRES_D_SINGLE_SHOOTING_CONTINUATION_MFGMRES_ARGS_H_