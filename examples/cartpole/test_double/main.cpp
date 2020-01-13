#include <string>
#include <sys/stat.h>

extern "C" {
#include "cgmres.h"
}
#include "d_simulator.hpp"


int main() {
  struct d_single_shooting_cgmres cgmres;
  double T_f = 2.0;
  double alpha = 1.0;
  double initial_time = 0;
  int N = 100;
  double finite_difference_increment = 1.0e-08;
  double zeta = 1000;
  int kmax = 10;
  double state[4] = {0, 0, 0, 0};

  d_single_shooting_cgmres_create(&cgmres, T_f, alpha, initial_time, N, 
                                  finite_difference_increment, zeta, kmax);
  double initial_guess[3] = {0.01, 10, 0.01};
  d_single_shooting_cgmres_set_initialization_parameters(&cgmres, 1.0e-08, 100, 
                                                         initial_guess);
  double sampling_time = 0.001;
  double simulation_time = 10;
  double current_time = 0;

  std::string save_dir_name("simulation_result");
  int mkdir_err = mkdir(save_dir_name.c_str(), 0755);
  simulation(&cgmres, state, initial_time, simulation_time, sampling_time, 
              "simulation_result", "d_cart_pole");

  return 0;
}