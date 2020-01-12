#ifndef S_CGMRES_SIMULATOR_H
#define S_CGMRES_SIMULATOR_H

#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include "numerical_integrator.hpp"
#include "save_simulation_data.hpp"

extern "C" {
#include "d_nmpc_model.h"
#include "hpcgmres.h"
}


void simulation(struct s_single_shooting_cgmres *cgmres, 
                const double* initial_state_vec, 
                const double start_time, const double end_time, 
                const double sampling_period, const std::string save_dir, 
                const std::string savefile_name);

#endif // S_CGMRES_SIMULATOR_H