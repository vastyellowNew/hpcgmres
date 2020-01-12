#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>

#include <gtest/gtest.h>

extern "C" {
#include "d_single_shooting_cgmres.h"
#include "d_mfgmres_for_single_shooting_cgmres.h"
#include "d_single_shooting_continuation.h"
#include "d_single_shooting_continuation_mfgmres_args.h"
#include "d_cgmres_initializer.h"
#include "d_memory_manager.h"
#include "d_linear_algebra.h"
#include "d_time_varying_smooth_horizon.h"
}


class d_single_shooting_cgmres_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand(time(NULL));
    T_f = 1.0;
    alpha = 1.0;
    initial_time = 0.0;
    finite_diff = 1.0e-08;
    zeta = 1000;
    current_time = fabs((double)rand()/RAND_MAX);
    N = 50;
    kmax = 100;
    dimx = d_single_shooting_continuation_dimx();
    dimu = d_single_shooting_continuation_dimu();
    dimc = d_single_shooting_continuation_dimc();
    dimuc = dimu + dimc;
    dim_solution = N * dimuc;
    d_single_shooting_cgmres_create(&cgmres, T_f, alpha, initial_time, N,
                                    finite_diff, zeta, kmax);
    d_single_shooting_continuation_create(&continuation, T_f, alpha, 
                                          initial_time, N, finite_diff, zeta);
    d_mfgmres_for_single_shooting_cgmres_create(&mfgmres, dim_solution, kmax);
    d_cgmres_initializer_create(&initializer, finite_diff, kmax);

    solution = allocate_dvec(dim_solution);
    solution_update = allocate_dvec(dim_solution);
    initial_solution = allocate_dvec(dimuc);
    state = allocate_dvec(dimx);
    for (int i=0; i<dimx; ++i) {
      state[i] = (double)rand()/RAND_MAX;
    }
    sampling_period = 0.001;
  }

  virtual void TearDown() {
    d_single_shooting_cgmres_delete(&cgmres);
    free_dvec(solution);
    free_dvec(solution_update);
    free_dvec(initial_solution);
    free_dvec(state);
  }

  struct d_single_shooting_cgmres cgmres;
  struct d_single_shooting_continuation continuation;
  struct d_single_shooting_continuation_mfgmres_args mfgmres_args;
  struct d_mfgmres_for_single_shooting_cgmres mfgmres;
  struct d_cgmres_initializer initializer;
  double T_f, alpha, initial_time, finite_diff, zeta, current_time;
  int N, kmax, dimx, dimu, dimc, dimuc, dim_solution;
  double *solution, *solution_update, *initial_solution, *state;
  double sampling_period;
};


TEST_F(d_single_shooting_cgmres_test, memsize) {
  int expect_size = 0;
  expect_size += d_mfgmres_for_single_shooting_cgmres_memsize(dim_solution, kmax);
  expect_size += d_single_shooting_continuation_memsize(N);
  expect_size += d_cgmres_initializer_memsize(kmax);
  expect_size += 2*dim_solution*sizeof(double); // solution, solution_update
  expect_size += 1*(dimu+dimc)*sizeof(double); // initial_solution
  expect_size = (expect_size+63)/64*64; // make multiple of typical cache line size
  expect_size += 64; // align to typical cache line size
  EXPECT_EQ(expect_size, d_single_shooting_cgmres_memsize(N, kmax));
}


TEST_F(d_single_shooting_cgmres_test, init_parameters) {
  double tol = 0.01;
  int max_itr = 50;
  double *initial_guess = allocate_dvec(dimuc);
  for (int i=0; i<dimuc; ++i) {
    initial_guess[i] = (double)rand()/RAND_MAX;
  }
  d_single_shooting_cgmres_set_initialization_parameters(&cgmres, tol, max_itr, 
                                                         initial_guess);
  d_cgmres_initializer_set_termination_criterions(&initializer, tol, max_itr);
  d_cgmres_initializer_set_initial_guess_solution(&initializer, initial_guess);
  EXPECT_EQ(cgmres.initializer.newton_residual_tolerance, 
            initializer.newton_residual_tolerance);
  EXPECT_EQ(cgmres.initializer.max_newton_iteration, 
            initializer.max_newton_iteration);
  for (int i=0; i<dimuc; ++i) {
    EXPECT_EQ(cgmres.initializer.initial_guess_solution[i], 
              initializer.initial_guess_solution[i]);
  }
  free_dvec(initial_guess);
}

TEST_F(d_single_shooting_cgmres_test, initialize_solution) {
  double tol = 1.0e-08;
  int max_itr = 100;
  double *initial_guess = allocate_dvec(dimuc);
  double *initial_solution_ref = allocate_dvec(dimuc);
  for (int i=0; i<dimuc; ++i) {
    initial_guess[i] = 0.1;
  }
  d_single_shooting_cgmres_set_initialization_parameters(&cgmres, tol, max_itr, 
                                                         initial_guess);
  d_cgmres_initializer_set_termination_criterions(&initializer, tol, max_itr);
  d_cgmres_initializer_set_initial_guess_solution(&initializer, initial_guess);
  d_single_shooting_cgmres_initialize_solution(&cgmres, current_time, state);
  d_cgmres_initializer_compute_initial_solution(&initializer, current_time, 
                                                state, initial_solution_ref);
  for (int i=0; i<N; ++i) {
    for (int j=0; j<dimuc; ++j) {
      EXPECT_EQ(cgmres.solution[i*dimuc+j], initial_solution_ref[j]);
    }
  }
  double horizon_length 
    = d_time_varying_smooth_horizon_get_length(&cgmres.continuation.ocp.horizon, 
                                               current_time);
  EXPECT_EQ(horizon_length, 0); 
  free_dvec(initial_guess);
  free_dvec(initial_solution_ref);
}


TEST_F(d_single_shooting_cgmres_test, update_control_input) {
  double *solution_ref = allocate_dvec(dim_solution);
  double *update_ref = allocate_dvec(dim_solution);
  double *control_input = allocate_dvec(dimu);
  for (int i=0; i<dim_solution; ++i) {
    cgmres.solution[i] = fabs((double)rand()/RAND_MAX);
    solution_ref[i] = cgmres.solution[i];
  }
  for (int i=0; i<dim_solution; ++i) {
    ASSERT_EQ(cgmres.solution_update[i], update_ref[i]);
  }
  for (int i=0; i<dim_solution; ++i) {
    ASSERT_EQ(cgmres.solution[i], solution_ref[i]);
  }
  d_single_shooting_cgmres_update_control_input(&cgmres, current_time, state,
                                                sampling_period, control_input);

  mfgmres_args.current_time = current_time;
  mfgmres_args.current_state_ptr = state;
  mfgmres_args.current_solution_ptr = solution_ref;
  d_single_shooting_continuation_reset_horizon_length(&continuation, 
                                                      initial_time);
  d_mfgmres_for_single_shooting_cgmres_solve_linear_problem(&mfgmres, 
                                                            &continuation, 
                                                            &mfgmres_args, 
                                                            update_ref);
  d_single_shooting_continuation_integrate_solution(&continuation, solution_ref,
                                                    update_ref, 
                                                    sampling_period);
  for (int i=0; i<dimu; ++i) {
    EXPECT_EQ(cgmres.solution[i], control_input[i]);
    control_input[i] += 100;
  }
  d_single_shooting_cgmres_get_control_input(&cgmres, control_input);
  for (int i=0; i<dimu; ++i) {
    EXPECT_EQ(cgmres.solution[i], control_input[i]);
  }

  for (int i=0; i<dim_solution; ++i) {
    EXPECT_EQ(cgmres.solution_update[i], update_ref[i]);
  }
  for (int i=0; i<dim_solution; ++i) {
    EXPECT_EQ(cgmres.solution[i], solution_ref[i]);
  }

  free_dvec(solution_ref);
  free_dvec(update_ref);
  free_dvec(control_input);
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}