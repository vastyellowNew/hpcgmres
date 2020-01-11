#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <gtest/gtest.h>

extern "C" {
#include "s_memory_manager.h"
#include "s_linear_algebra.h"
#include "s_zero_horizon_ocp.h"
#include "s_inexact_newton_for_zero_horizon_ocp.h"
#include "s_inexact_newton_for_zero_horizon_ocp_mfgmres_args.h"
}


class s_inexact_newton_for_zero_horizon_ocp_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand(time(NULL));
    finite_diff = 1.0e-08;
    s_inexact_newton_for_zero_horizon_ocp_create(&newton, finite_diff);
    dimx = s_inexact_newton_for_zero_horizon_ocp_dimx();
    dimu = s_inexact_newton_for_zero_horizon_ocp_dimu();
    dimc = s_inexact_newton_for_zero_horizon_ocp_dimc();
    dim_solution = dimu + dimc;
    state = allocate_svec(dimx);
    solution = allocate_svec(dim_solution);
    inc_solution = allocate_svec(dim_solution);
    direction = allocate_svec(dim_solution);
    b = allocate_svec(dim_solution);
    ax = allocate_svec(dim_solution);
    b_ref = allocate_svec(dim_solution);
    ax_ref = allocate_svec(dim_solution);
    for (int i=0; i<dimx; ++i) {
      state[i] = (float)rand()/RAND_MAX;
    }
    for (int i=0; i<dim_solution; ++i) {
      solution[i] = (float)rand()/RAND_MAX;
      direction[i] = (float)rand()/RAND_MAX;
    }
  }

  virtual void TearDown() {
    s_inexact_newton_for_zero_horizon_ocp_delete(&newton);
    free_svec(state);
    free_svec(solution);
    free_svec(inc_solution);
    free_svec(direction);
    free_svec(b);
    free_svec(ax);
    free_svec(b_ref);
    free_svec(ax_ref);
  }

  struct s_inexact_newton_for_zero_horizon_ocp newton;
  struct s_inexact_newton_for_zero_horizon_ocp_mfgmres_args args;
  int dimx, dimu, dimc, dim_solution;
  float finite_diff, current_time;
  float *state, *solution, *inc_solution, *direction, *b, *ax, *b_ref, *ax_ref;
};


TEST_F(s_inexact_newton_for_zero_horizon_ocp_test, memsize) {
  int expected_size = 0;
  expected_size += 1*s_zero_horizon_ocp_memsize();
  expected_size += 3*(dimu+dimc)*sizeof(float); // incremented_solution, optimality_residual, optimality_residual1;
  expected_size = (expected_size+63)/64*64; // make multiple of typical cache line size
  expected_size += 64; // align to typical cache line size
  EXPECT_EQ(expected_size, s_inexact_newton_for_zero_horizon_ocp_memsize());
}


TEST_F(s_inexact_newton_for_zero_horizon_ocp_test, b_and_ax) {
  args.current_time = current_time;
  args.current_state_ptr = state;
  args.current_solution_ptr = solution;
  float *opt_res = allocate_svec(dim_solution);
  float *opt_res1 = allocate_svec(dim_solution);
  hpcgmres_saxpy(dim_solution, finite_diff, direction, solution, inc_solution);
  s_zero_horizon_ocp_compute_optimality_residual(&newton.ocp, current_time, 
                                                 state, solution, opt_res);
  s_zero_horizon_ocp_compute_optimality_residual(&newton.ocp, current_time, 
                                                 state, inc_solution, opt_res1);
  hpcgmres_saxpby(dim_solution, 1/finite_diff-1, opt_res, -1/finite_diff, 
                  opt_res1, b_ref);
  s_inexact_newton_for_zero_horizon_ocp_compute_b(&newton, &args, direction, b);
  for (int i=0; i<dim_solution; ++i) {
    EXPECT_EQ(b_ref[i], b[i]);
  }

  hpcgmres_saxpby(dim_solution, -1/finite_diff, opt_res, 1/finite_diff, 
                  opt_res1, ax_ref);
  s_inexact_newton_for_zero_horizon_ocp_compute_ax(&newton, &args, direction, ax);
  for (int i=0; i<dim_solution; ++i) {
    EXPECT_EQ(ax_ref[i], ax[i]);
  }
  free_svec(opt_res1);
  free_svec(opt_res);
}


TEST_F(s_inexact_newton_for_zero_horizon_ocp_test, error_norm) {
  float *opt_res = allocate_svec(dim_solution);
  s_zero_horizon_ocp_compute_optimality_residual(&newton.ocp, current_time, 
                                                 state, solution, opt_res);
  float err_norm_ref = sqrt(hpcgmres_svecnrm2(dim_solution, opt_res));
  float err_norm = s_inexact_newton_for_zero_horizon_ocp_get_error_norm(&newton, current_time, state, solution);
  EXPECT_EQ(err_norm_ref, err_norm);
  free_svec(opt_res);
}


TEST_F(s_inexact_newton_for_zero_horizon_ocp_test, terminal_cost_derivative) {
  float *phix = allocate_svec(dimx);
  float *phix_ref = allocate_svec(dimx);
  s_zero_horizon_ocp_compute_terminal_cost_derivative(&newton.ocp, current_time, 
                                                      state, phix_ref);
  s_inexact_newton_for_zero_horizon_ocp_get_terminal_cost_derivative(&newton, current_time, state, phix);
  for (int i=0; i<dimx; ++i) {
    EXPECT_EQ(phix_ref[i], phix[i]);
  }
  free_svec(phix);
  free_svec(phix_ref);
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}