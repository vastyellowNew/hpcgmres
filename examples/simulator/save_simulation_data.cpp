#include "save_simulation_data.hpp"


void saveData(const int dim_state, const int dim_control_input, 
              std::ofstream& state_data, std::ofstream& control_input_data, 
              std::ofstream& error_data, const double time_param, 
              const double* state_vec, const double* control_input_vec,
              const double error_norm) {
  for (int i=0; i<dim_state; i++) {
    state_data << state_vec[i] << " ";
  }
  state_data << "\n";
  for (int i=0; i<dim_control_input; i++) {
    control_input_data << control_input_vec[i] << " ";
  }
  control_input_data << "\n";
  error_data << error_norm << "\n";
}