#ifndef HPCGMRES_S_SINGLE_SHOOTING_CONTINUATION_MFGMRES_ARGS_H_
#define HPCGMRES_S_SINGLE_SHOOTING_CONTINUATION_MFGMRES_ARGS_H_


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


#endif // HPCGMRES_S_SINGLE_SHOOTING_CONTINUATION_MFGMRES_ARGS_H_