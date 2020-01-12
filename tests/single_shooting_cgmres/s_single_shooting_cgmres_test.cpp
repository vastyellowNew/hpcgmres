#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>

#include <gtest/gtest.h>

extern "C" {
#include "s_single_shooting_cgmres.h"
#include "s_mfgmres_for_single_shooting_cgmres.h"
#include "s_single_shooting_continuation.h"
#include "s_single_shooting_continuation_mfgmres_args.h"
#include "s_cgmres_initializer.h"
#include "s_memory_manager.h"
#include "s_linear_algebra.h"
#include "s_time_varying_smooth_horizon.h"
}


class s_single_shooting_cgmres_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand(time(NULL));
    T_f = 1.0;
    alpha = 1.0;
    initial_time = 0.0;
    finite_diff = 1.0e-03;
    zeta = 1000;
    current_time = fabs((float)rand()/RAND_MAX);
    N = 50;
    kmax = 100;
    dimx = s_single_shooting_continuation_dimx();
    dimu = s_single_shooting_continuation_dimu();
    dimc = s_single_shooting_continuation_dimc();
    dimuc = dimu + dimc;
    dim_solution = N * dimuc;
    s_single_shooting_cgmres_create(&cgmres, T_f, alpha, initial_time, N,
                                    finite_diff, zeta, kmax);
    s_single_shooting_continuation_create(&continuation, T_f, alpha, 
                                          initial_time, N, finite_diff, zeta);
    s_mfgmres_for_single_shooting_cgmres_create(&mfgmres, dim_solution, kmax);
    s_cgmres_initializer_create(&initializer, finite_diff, kmax);

    solution = allocate_svec(dim_solution);
    solution_update = allocate_svec(dim_solution);
    initial_solution = allocate_svec(dimuc);
    state = allocate_svec(dimx);
    for (int i=0; i<dimx; ++i) {
      state[i] = (float)rand()/RAND_MAX;
    }
    sampling_period = 0.001;
  }

  virtual void TearDown() {
    s_single_shooting_cgmres_delete(&cgmres);
    s_single_shooting_continuation_delete(&continuation);
    s_mfgmres_for_single_shooting_cgmres_delete(&mfgmres);
    s_cgmres_initializer_delete(&initializer);
    free_svec(solution);
    free_svec(solution_update);
    free_svec(initial_solution);
    free_svec(state);
  }

  struct s_single_shooting_cgmres cgmres;
  struct s_single_shooting_continuation continuation;
  struct s_single_shooting_continuation_mfgmres_args mfgmres_args;
  struct s_mfgmres_for_single_shooting_cgmres mfgmres;
  struct s_cgmres_initializer initializer;
  float T_f, alpha, initial_time, finite_diff, zeta, current_time;
  int N, kmax, dimx, dimu, dimc, dimuc, dim_solution;
  float *solution, *solution_update, *initial_solution, *state;
  float sampling_period;
};


TEST_F(s_single_shooting_cgmres_test, memsize) {
  int expect_size = 0;
  expect_size += s_mfgmres_for_single_shooting_cgmres_memsize(dim_solution, kmax);
  expect_size += s_single_shooting_continuation_memsize(N);
  expect_size += s_cgmres_initializer_memsize(kmax);
  expect_size += 2*dim_solution*sizeof(float); // solution, solution_update
  expect_size += 1*(dimu+dimc)*sizeof(float); // initial_solution
  expect_size = (expect_size+63)/64*64; // make multiple of typical cache line size
  expect_size += 64; // align to typical cache line size
  EXPECT_EQ(expect_size, s_single_shooting_cgmres_memsize(N, kmax));
}


TEST_F(s_single_shooting_cgmres_test, init_parameters) {
  float tol = 0.01;
  int max_itr = 50;
  float *initial_guess = allocate_svec(dimuc);
  for (int i=0; i<dimuc; ++i) {
    initial_guess[i] = (float)rand()/RAND_MAX;
  }
  s_single_shooting_cgmres_set_initialization_parameters(&cgmres, tol, max_itr, 
                                                         initial_guess);
  s_cgmres_initializer_set_termination_criterions(&initializer, tol, max_itr);
  s_cgmres_initializer_set_initial_guess_solution(&initializer, initial_guess);
  EXPECT_EQ(cgmres.initializer.newton_residual_tolerance, 
            initializer.newton_residual_tolerance);
  EXPECT_EQ(cgmres.initializer.max_newton_iteration, 
            initializer.max_newton_iteration);
  for (int i=0; i<dimuc; ++i) {
    EXPECT_EQ(cgmres.initializer.initial_guess_solution[i], 
              initializer.initial_guess_solution[i]);
  }
  free_svec(initial_guess);
}

TEST_F(s_single_shooting_cgmres_test, initialize_solution) {
  float tol = 1.0e-08;
  int max_itr = 100;
  float *initial_guess = allocate_svec(dimuc);
  float *initial_solution_ref = allocate_svec(dimuc);
  for (int i=0; i<dimuc; ++i) {
    initial_guess[i] = 0.1;
  }
  s_single_shooting_cgmres_set_initialization_parameters(&cgmres, tol, max_itr, 
                                                         initial_guess);
  s_cgmres_initializer_set_termination_criterions(&initializer, tol, max_itr);
  s_cgmres_initializer_set_initial_guess_solution(&initializer, initial_guess);
  s_single_shooting_cgmres_initialize_solution(&cgmres, current_time, state);
  s_cgmres_initializer_compute_initial_solution(&initializer, current_time, 
                                                state, initial_solution_ref);
  for (int i=0; i<N; ++i) {
    for (int j=0; j<dimuc; ++j) {
      EXPECT_EQ(cgmres.solution[i*dimuc+j], initial_solution_ref[j]);
    }
  }
  float horizon_length 
    = s_time_varying_smooth_horizon_get_length(&cgmres.continuation.ocp.horizon, 
                                               current_time);
  EXPECT_EQ(horizon_length, 0); 
  free_svec(initial_guess);
  free_svec(initial_solution_ref);
}


TEST_F(s_single_shooting_cgmres_test, update_control_input) {
  float *solution_ref = allocate_svec(dim_solution);
  float *update_ref = allocate_svec(dim_solution);
  float *control_input = allocate_svec(dimu);
  for (int i=0; i<dim_solution; ++i) {
    cgmres.solution[i] = fabs((float)rand()/RAND_MAX);
    solution_ref[i] = cgmres.solution[i];
  }
  for (int i=0; i<dim_solution; ++i) {
    ASSERT_EQ(cgmres.solution_update[i], update_ref[i]);
  }
  for (int i=0; i<dim_solution; ++i) {
    ASSERT_EQ(cgmres.solution[i], solution_ref[i]);
  }
  s_single_shooting_cgmres_update_control_input(&cgmres, current_time, state,
                                                sampling_period, control_input);

  mfgmres_args.current_time = current_time;
  mfgmres_args.current_state_ptr = state;
  mfgmres_args.current_solution_ptr = solution_ref;
  s_single_shooting_continuation_reset_horizon_length(&continuation, 
                                                      initial_time);
  s_mfgmres_for_single_shooting_cgmres_solve_linear_problem(&mfgmres, 
                                                            &continuation, 
                                                            &mfgmres_args, 
                                                            update_ref);
  s_single_shooting_continuation_integrate_solution(&continuation, solution_ref,
                                                    update_ref, 
                                                    sampling_period);
  for (int i=0; i<dimu; ++i) {
    EXPECT_EQ(cgmres.solution[i], control_input[i]);
    control_input[i] += 100;
  }
  s_single_shooting_cgmres_get_control_input(&cgmres, control_input);
  for (int i=0; i<dimu; ++i) {
    EXPECT_EQ(cgmres.solution[i], control_input[i]);
  }

  for (int i=0; i<dim_solution; ++i) {
    EXPECT_EQ(cgmres.solution_update[i], update_ref[i]);
  }
  for (int i=0; i<dim_solution; ++i) {
    EXPECT_EQ(cgmres.solution[i], solution_ref[i]);
  }

  free_svec(solution_ref);
  free_svec(update_ref);
  free_svec(control_input);
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}