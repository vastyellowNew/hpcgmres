#ifndef HPCGMRES_S_TIME_VARYING_SMOOTH_HORIZON_H_
#define HPCGMRES_S_TIME_VARYING_SMOOTH_HORIZON_H_


#include <math.h>


#ifdef __cplusplus
extern "C" {
#endif


struct s_time_varying_smooth_horizon {
  float T_f;
  float alpha;
  float initial_time;
};

int s_time_varying_smooth_horizon_strsize();

void s_time_varying_smooth_horizon_create(
    struct s_time_varying_smooth_horizon *horizon, float T_f, float alpha, 
    float initial_time);

float s_time_varying_smooth_horizon_get_length(
    struct s_time_varying_smooth_horizon *horizon, float current_time);

void s_time_varying_smooth_horizon_reset_length(
    struct s_time_varying_smooth_horizon *horizon, float initial_time);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // HPCGMRES_S_TIME_VARYING_SMOOTH_HORIZON_H_