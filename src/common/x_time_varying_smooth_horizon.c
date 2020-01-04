int TIME_VARYING_SMOOTH_HORIZON_STRSIZE() {
  return sizeof(struct TIME_VARYING_SMOOTH_HORIZON);
}


void TIME_VARYING_SMOOTH_HORIZON_CREATE(
    struct TIME_VARYING_SMOOTH_HORIZON *horizon, REAL T_f, REAL alpha, 
    REAL initial_time) {
  horizon->T_f = T_f;
  horizon->alpha = alpha;
  horizon->initial_time = initial_time;
}


REAL TIME_VARYING_SMOOTH_HORIZON_GET_LENGTH(
    struct TIME_VARYING_SMOOTH_HORIZON *horizon, REAL current_time) {
  return horizon->T_f 
      * (1.0-exp(-horizon->alpha*(current_time-horizon->initial_time)));
}


void TIME_VARYING_SMOOTH_HORIZON_RESET_LENGTH(
    struct TIME_VARYING_SMOOTH_HORIZON *horizon, REAL initial_time) {
  horizon->initial_time = initial_time;
}