#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <gtest/gtest.h>

extern "C" {
#include "d_cgmres_initializer.h"
#include "d_memory_manager.h"
#include "d_linear_algebra.h"
#include "d_inexact_newton_for_zero_horizon_ocp.h"
#include "d_inexact_newton_for_zero_horizon_ocp_mfgmres_args.h"
#include "d_mfgmres_for_cgmres_initializer.h"
}


class d_cgmres_initializer_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand(time(NULL));
    finite_diff = 1.0e-08;
    dimx = d_inexact_newton_for_zero_horizon_ocp_dimx();
    dimu = d_inexact_newton_for_zero_horizon_ocp_dimu();
    dimc = d_inexact_newton_for_zero_horizon_ocp_dimc();
    dimuc = dimu + dimc;
    dim_solution = dimuc;
    kmax = dim_solution;
    d_cgmres_initializer_create(&initializer, finite_diff, kmax);
    d_mfgmres_for_cgmres_initializer_create(&mfgmres, dim_solution, kmax);
    d_inexact_newton_for_zero_horizon_ocp_create(&newton, finite_diff);
    state = allocate_dvec(dimx);
    solution = allocate_dvec(dim_solution);
    solution_ref = allocate_dvec(dim_solution);
    update = allocate_dvec(dim_solution);
    current_time = (double)rand()/RAND_MAX;
    for (int i=0; i<dimx; ++i) {
      state[i] = (double)rand()/RAND_MAX; 
    }
    for (int i=0; i<dim_solution; ++i) {
      solution[i] = (double)rand()/RAND_MAX;
    }
  }

  virtual void TearDown() {
    d_cgmres_initializer_delete(&initializer);
    d_mfgmres_for_cgmres_initializer_delete(&mfgmres);
    d_inexact_newton_for_zero_horizon_ocp_delete(&newton);
    free_dvec(state);
    free_dvec(solution);
    free_dvec(solution_ref);
    free_dvec(update);
  }

  struct d_cgmres_initializer initializer;
  struct d_mfgmres_for_cgmres_initializer mfgmres;
  struct d_inexact_newton_for_zero_horizon_ocp newton;
  struct d_inexact_newton_for_zero_horizon_ocp_mfgmres_args args;
  int dimx, dimu, dimc, dimuc, dim_solution, kmax;
  double finite_diff, current_time;
  double *state, *solution, *update, *solution_ref;
};


TEST_F(d_cgmres_initializer_test, memsize) {
  int expected_size = 0;
  expected_size += d_mfgmres_for_cgmres_initializer_memsize(dimuc, kmax);
  expected_size += d_inexact_newton_for_zero_horizon_ocp_memsize();
  expected_size += 2*dimuc*sizeof(double); 
  expected_size = (expected_size+63)/64*64; 
  expected_size += 64; 
  EXPECT_EQ(expected_size, d_cgmres_initializer_memsize(kmax));
}


TEST_F(d_cgmres_initializer_test, parameters) {
  EXPECT_EQ(1.0e-08, initializer.newton_residual_tolerance);
  EXPECT_EQ(50, initializer.max_newton_iteration);

  double tol = 0.5;
  double max_n = 100;
  d_cgmres_initializer_set_termination_criterions(&initializer, tol, max_n);
  EXPECT_EQ(tol, initializer.newton_residual_tolerance);
  EXPECT_EQ(max_n, initializer.max_newton_iteration);

  for (int i=0; i<dim_solution; ++i) {
    solution[i] = (double)rand()/RAND_MAX;
  }
  d_cgmres_initializer_set_initial_guess_solution(&initializer, solution);
  for (int i=0; i<dim_solution; ++i) {
    EXPECT_EQ(solution[i], initializer.initial_guess_solution[i]);
  }
}


TEST_F(d_cgmres_initializer_test, compute_initial_solution) {
  d_cgmres_initializer_set_initial_guess_solution(&initializer, solution);
  d_cgmres_initializer_compute_initial_solution(&initializer, current_time,
                                                state, solution_ref);
  args.current_time = current_time;
  args.current_state_ptr = state;
  args.current_solution_ptr = solution;
  int num_itr = 0;
  double optimality_error 
      = d_inexact_newton_for_zero_horizon_ocp_get_error_norm(&newton, 
                                                             current_time, 
                                                             state, solution);
  while (optimality_error > 1.0e-08 && num_itr < 50) {
    d_mfgmres_for_cgmres_initializer_solve_linear_problem(&mfgmres, &newton, 
                                                          &args, update);
    hpcgmres_dvecadd(dim_solution, update, solution);
    optimality_error 
        = d_inexact_newton_for_zero_horizon_ocp_get_error_norm(&newton, 
                                                               current_time, 
                                                               state, solution);
    ++num_itr;
  }
  for (int i=0; i<dim_solution; ++i) {
    EXPECT_EQ(solution[i], solution_ref[i]);
  }
}


TEST_F(d_cgmres_initializer_test, get_terminal_cost_derivative) {
  double *phix_ref = allocate_dvec(dimx);
  double *phix = allocate_dvec(dimx);
  d_cgmres_initializer_get_terminal_cost_derivative(&initializer, current_time,
                                                    state, phix_ref);

  d_inexact_newton_for_zero_horizon_ocp_get_terminal_cost_derivative(
      &newton, current_time, state, phix);
  for (int i=0; i<dimx; ++i) {
    EXPECT_EQ(phix_ref[i], phix[i]);
  }
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}