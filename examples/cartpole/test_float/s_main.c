#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "hpcgmres.h"
#include "d_nmpc_model.h"


int main() {
  struct s_single_shooting_cgmres cgmres;
  struct s_nmpc_model model;
  float T_f = 2.0;
  float alpha = 1.0;
  float initial_time = 0;
  int N = 100;
  float finite_difference_increment = 1.0e-03;
  float zeta = 1000;
  int kmax = 10;
  int dimx = 4;
  float state[4] = {0, 0, 0, 0};

  s_single_shooting_cgmres_create(&cgmres, T_f, alpha, initial_time, N, 
                                  finite_difference_increment, zeta, kmax);
  s_nmpc_model_create(&model);
  float initial_guess[3] = {0.1, 0.1, 0.1};
  s_single_shooting_cgmres_set_initialization_parameters(&cgmres, 1.0e-08, 100, 
                                                         initial_guess);
  s_single_shooting_cgmres_initialize_solution(&cgmres, initial_time, state);

  float control_input[1], next_state[4], dx[4];
  float sampling_time = 0.001;
  float simulation_time = 10;
  float current_time = 0;
  s_single_shooting_cgmres_get_control_input(&cgmres, control_input);
  printf("start %lf [s] simulation.\n", simulation_time);
  clock_t start_time, end_time, total_time;
  int count = 0;
  total_time = 0;
  for (; current_time < simulation_time; current_time+=sampling_time, ++count) {
    float error_norm = s_single_shooting_cgmres_get_error_norm(&cgmres, 
                                                                current_time, 
                                                                state);
    printf("error norm = %lf\n", error_norm);
    s_nmpc_model_f(&model, current_time, state, control_input, dx);
    int i;
    for (i=0; i<dimx; ++i) {
      next_state[i] = state[i] + sampling_time * dx[i];
    }
    start_time = clock();
    s_single_shooting_cgmres_update_control_input(&cgmres, current_time, 
                                                  state, sampling_time,
                                                  control_input);
    end_time = clock();
    total_time += end_time - start_time;
    for (i=0; i<dimx; ++i) {
      state[i] = next_state[i];
    }
  }
  double total_cpu_time = ((double)(total_time*10000/CLOCKS_PER_SEC));
  total_cpu_time /= 10000;
  double ave_time = total_cpu_time / count;
  printf("end simulation. average computational time = %lf[ms].\n", 1000*ave_time);

  return 0;
}