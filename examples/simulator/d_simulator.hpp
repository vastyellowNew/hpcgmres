#ifndef D_CGMRES_SIMULATOR_H
#define D_CGMRES_SIMULATOR_H

#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include "numerical_integrator.hpp"
#include "save_simulation_data.hpp"

extern "C" {
#include "d_nmpc_model.h"
#include "cgmres.h"
}


void simulation(struct d_single_shooting_cgmres *cgmres, 
                const double* initial_state_vec, 
                const double start_time, const double end_time, 
                const double sampling_period, const std::string save_dir, 
                const std::string savefile_name);

#endif // D_CGMRES_SIMULATOR_H