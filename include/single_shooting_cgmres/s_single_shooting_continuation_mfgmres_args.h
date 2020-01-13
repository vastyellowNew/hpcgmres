#ifndef CGMRES_S_SINGLE_SHOOTING_CONTINUATION_MFGMRES_ARGS_H_
#define CGMRES_S_SINGLE_SHOOTING_CONTINUATION_MFGMRES_ARGS_H_


#ifdef __cplusplus
extern "C" {
#endif


struct s_single_shooting_continuation_mfgmres_args {
  float current_time;
  float *current_state_ptr;
  float *current_solution_ptr;
};


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // CGMRES_S_SINGLE_SHOOTING_CONTINUATION_MFGMRES_ARGS_H_