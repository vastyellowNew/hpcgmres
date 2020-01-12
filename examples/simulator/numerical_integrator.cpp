#include "numerical_integrator.hpp"


NumericalIntegrator::NumericalIntegrator() {
  d_nmpc_model_create(&model_);
}

void NumericalIntegrator::euler(double current_time, 
                                double* current_state_vec, 
                                double* control_input_vec, 
                                double integration_length, 
                                double* integrated_state) {
  double dx_vec_[d_nmpc_model_dimx()];
  d_nmpc_model_f(&model_, current_time, current_state_vec, control_input_vec, 
                 dx_vec_);
  for (int i=0; i<d_nmpc_model_dimx(); i++) {
    integrated_state[i] = current_state_vec[i] + integration_length*dx_vec_[i];
  }
}

void NumericalIntegrator::rungeKuttaGill(double current_time, 
                                         double* current_state_vec, 
                                         double* control_input_vec, 
                                         double integration_length, 
                                         double* integrated_state) {
  double k1_vec[d_nmpc_model_dimx()],  k2_vec[d_nmpc_model_dimx()], 
      k3_vec[d_nmpc_model_dimx()], k4_vec[d_nmpc_model_dimx()], 
      tmp_vec[d_nmpc_model_dimx()];

  d_nmpc_model_f(&model_, current_time, current_state_vec, control_input_vec, 
                 k1_vec);
  for (int i=0; i<d_nmpc_model_dimx(); i++) {
      tmp_vec[i] = current_state_vec[i] + 0.5*integration_length*k1_vec[i];
  }

  d_nmpc_model_f(&model_, current_time+0.5*integration_length, tmp_vec, 
                 control_input_vec, k2_vec);
  for (int i=0; i<d_nmpc_model_dimx(); i++) {
    tmp_vec[i] = current_state_vec[i] 
        + integration_length*0.5*(std::sqrt(2)-1)*k1_vec[i] 
        + integration_length*(1-(1/std::sqrt(2)))*k2_vec[i];
  }

  d_nmpc_model_f(&model_, current_time+0.5*integration_length, tmp_vec, 
                 control_input_vec, k3_vec);
  for (int i=0; i<d_nmpc_model_dimx(); i++) {
    tmp_vec[i] = current_state_vec[i] 
        - integration_length*0.5*std::sqrt(2)*k2_vec[i] 
        + integration_length*(1+(1/std::sqrt(2)))*k3_vec[i];
  }

  d_nmpc_model_f(&model_, current_time+integration_length, tmp_vec, 
                 control_input_vec, k4_vec);
  for (int i=0; i<d_nmpc_model_dimx(); i++) {
    integrated_state[i] = current_state_vec[i] 
        + (integration_length/6)
        * (k1_vec[i]+(2-std::sqrt(2))*k2_vec[i]
            +(2+std::sqrt(2))*k3_vec[i]+k4_vec[i]);
  }
}