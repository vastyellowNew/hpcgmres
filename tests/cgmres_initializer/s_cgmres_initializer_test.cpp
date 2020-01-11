#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <gtest/gtest.h>

extern "C" {
#include "s_cgmres_initializer.h"
#include "s_memory_manager.h"
#include "s_linear_algebra.h"
#include "s_inexact_newton_for_zero_horizon_ocp.h"
#include "s_inexact_newton_for_zero_horizon_ocp_mfgmres_args.h"
#include "s_mfgmres_for_cgmres_initializer.h"
}


class s_cgmres_initializer_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand(time(NULL));
    finite_diff = 1.0e-02;
    delta = finite_diff;
    dimx = s_inexact_newton_for_zero_horizon_ocp_dimx();
    dimu = s_inexact_newton_for_zero_horizon_ocp_dimu();
    dimc = s_inexact_newton_for_zero_horizon_ocp_dimc();
    dimuc = dimu + dimc;
    dim_solution = dimuc;
    kmax = dim_solution;
    s_cgmres_initializer_create(&initializer, finite_diff, kmax);
    s_mfgmres_for_cgmres_initializer_create(&mfgmres, dim_solution, kmax);
    s_inexact_newton_for_zero_horizon_ocp_create(&newton, finite_diff);
    state = allocate_svec(dimx);
    solution = allocate_svec(dim_solution);
    solution_ref = allocate_svec(dim_solution);
    update = allocate_svec(dim_solution);
    current_time = (float)rand()/RAND_MAX;
    for (int i=0; i<dimx; ++i) {
      state[i] = (float)rand()/RAND_MAX; 
    }
    for (int i=0; i<dim_solution; ++i) {
      solution[i] = (float)rand()/RAND_MAX;
    }
  }

  virtual void TearDown() {
    s_cgmres_initializer_delete(&initializer);
    s_mfgmres_for_cgmres_initializer_delete(&mfgmres);
    s_inexact_newton_for_zero_horizon_ocp_delete(&newton);
    free_svec(state);
    free_svec(solution);
    free_svec(solution_ref);
    free_svec(update);
  }

  struct s_cgmres_initializer initializer;
  struct s_mfgmres_for_cgmres_initializer mfgmres;
  struct s_inexact_newton_for_zero_horizon_ocp newton;
  struct s_inexact_newton_for_zero_horizon_ocp_mfgmres_args args;
  int dimx, dimu, dimc, dimuc, dim_solution, kmax;
  float finite_diff, current_time, delta;
  float *state, *solution, *update, *solution_ref;
};


TEST_F(s_cgmres_initializer_test, memsize) {
  int expectes_size = 0;
  expectes_size += s_mfgmres_for_cgmres_initializer_memsize(dimuc, kmax);
  expectes_size += s_inexact_newton_for_zero_horizon_ocp_memsize();
  expectes_size += 2*dimuc*sizeof(float); 
  expectes_size = (expectes_size+63)/64*64; 
  expectes_size += 64; 
  EXPECT_EQ(expectes_size, s_cgmres_initializer_memsize(kmax));
}


TEST_F(s_cgmres_initializer_test, parameters) {
  EXPECT_EQ((float)1.0e-08, initializer.newton_residual_tolerance);
  EXPECT_EQ(50, initializer.max_newton_iteration);

  float tol = 0.5;
  float max_n = 100;
  s_cgmres_initializer_set_termination_criterions(&initializer, tol, max_n);
  EXPECT_EQ(tol, initializer.newton_residual_tolerance);
  EXPECT_EQ(max_n, initializer.max_newton_iteration);

  for (int i=0; i<dim_solution; ++i) {
    solution[i] = (float)rand()/RAND_MAX;
  }
  s_cgmres_initializer_set_initial_guess_solution(&initializer, solution);
  for (int i=0; i<dim_solution; ++i) {
    EXPECT_EQ(solution[i], initializer.initial_guess_solution[i]);
  }
}


TEST_F(s_cgmres_initializer_test, compute_initial_solution) {
  s_cgmres_initializer_set_initial_guess_solution(&initializer, solution);
  s_cgmres_initializer_compute_initial_solution(&initializer, current_time,
                                                state, solution_ref);
  args.current_time = current_time;
  args.current_state_ptr = state;
  args.current_solution_ptr = solution;
  int num_itr = 0;
  float optimality_error 
      = s_inexact_newton_for_zero_horizon_ocp_get_error_norm(&newton, 
                                                             current_time, 
                                                             state, solution);
  while (optimality_error > 1.0e-08 && num_itr < 50) {
    s_mfgmres_for_cgmres_initializer_solve_linear_problem(&mfgmres, &newton, 
                                                          &args, update);
    hpcgmres_svecadd(dim_solution, update, solution);
    optimality_error 
        = s_inexact_newton_for_zero_horizon_ocp_get_error_norm(&newton, 
                                                               current_time, 
                                                               state, solution);
    ++num_itr;
  }
  for (int i=0; i<dim_solution; ++i) {
    EXPECT_EQ(solution[i], solution_ref[i]);
  }
}


TEST_F(s_cgmres_initializer_test, get_terminal_cost_derivative) {
  float *phix_ref = allocate_svec(dimx);
  float *phix = allocate_svec(dimx);
  s_cgmres_initializer_get_terminal_cost_derivative(&initializer, current_time,
                                                    state, phix_ref);

  s_inexact_newton_for_zero_horizon_ocp_get_terminal_cost_derivative(
      &newton, current_time, state, phix);
  for (int i=0; i<dimx; ++i) {
    EXPECT_EQ(phix_ref[i], phix[i]);
  }
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}