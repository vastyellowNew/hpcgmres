#ifndef NUMERICAL_INTEGRATOR_H
#define NUMERICAL_INTEGRATOR_H


#include <cmath>

extern "C" {
#include "d_nmpc_model.h"
}


// Supports numerical integration of the state equation of the system described 
// in nmpc_model.hpp for numerical simnulations.
class NumericalIntegrator {
public:
  NumericalIntegrator();

  // Euler method for the state equation.
  void euler(double current_time, double* current_state_vec, 
             double* control_input_vec, double integration_length,
             double* integrated_state);

  // The four-step Runge-Kutta-Gill method for the state equation.
  void rungeKuttaGill(double current_time, 
                      double* current_state_vec, 
                      double* control_input_vec, 
                      double integration_length, 
                      double* integrated_state);

private:
  struct d_nmpc_model model_;
};


#endif // NUMERICAL_INTEGRATOR_H