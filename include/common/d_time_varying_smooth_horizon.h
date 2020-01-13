#ifndef CGMRES_D_TIME_VARYING_SMOOTH_HORIZON_H_
#define CGMRES_D_TIME_VARYING_SMOOTH_HORIZON_H_


#ifdef __cplusplus
extern "C" {
#endif


struct d_time_varying_smooth_horizon {
  double T_f;
  double alpha;
  double initial_time;
};

int d_time_varying_smooth_horizon_strsize();

void d_time_varying_smooth_horizon_create(
    struct d_time_varying_smooth_horizon *horizon, double T_f, double alpha, 
    double initial_time);

double d_time_varying_smooth_horizon_get_length(
    struct d_time_varying_smooth_horizon *horizon, double current_time);

void d_time_varying_smooth_horizon_reset_length(
    struct d_time_varying_smooth_horizon *horizon, double initial_time);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // CGMRES_D_TIME_VARYING_SMOOTH_HORIZON_H_