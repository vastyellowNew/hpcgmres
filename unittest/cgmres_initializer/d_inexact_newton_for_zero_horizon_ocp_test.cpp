#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <gtest/gtest.h>

extern "C" {
#include "d_memory_manager.h"
#include "d_linear_algebra.h"
#include "d_zero_horizon_ocp.h"
#include "d_inexact_newton_for_zero_horizon_ocp.h"
#include "d_inexact_newton_for_zero_horizon_ocp_mfgmres_args.h"
}


class d_inexact_newton_for_zero_horizon_ocp_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand(time(NULL));
    finite_diff = 1.0e-08;
    d_inexact_newton_for_zero_horizon_ocp_create(&newton, finite_diff);
    dimx = d_inexact_newton_for_zero_horizon_ocp_dimx();
    dimu = d_inexact_newton_for_zero_horizon_ocp_dimu();
    dimc = d_inexact_newton_for_zero_horizon_ocp_dimc();
    dim_solution = dimu + dimc;
    state = allocate_dvec(dimx);
    solution = allocate_dvec(dim_solution);
    inc_solution = allocate_dvec(dim_solution);
    direction = allocate_dvec(dim_solution);
    b = allocate_dvec(dim_solution);
    ax = allocate_dvec(dim_solution);
    b_ref = allocate_dvec(dim_solution);
    ax_ref = allocate_dvec(dim_solution);
    for (int i=0; i<dimx; ++i) {
      state[i] = (double)rand()/RAND_MAX;
    }
    for (int i=0; i<dim_solution; ++i) {
      solution[i] = (double)rand()/RAND_MAX;
      direction[i] = (double)rand()/RAND_MAX;
    }
  }

  virtual void TearDown() {
    d_inexact_newton_for_zero_horizon_ocp_delete(&newton);
    free_dvec(state);
    free_dvec(solution);
    free_dvec(inc_solution);
    free_dvec(direction);
    free_dvec(b);
    free_dvec(ax);
    free_dvec(b_ref);
    free_dvec(ax_ref);
  }

  struct d_inexact_newton_for_zero_horizon_ocp newton;
  struct d_inexact_newton_for_zero_horizon_ocp_mfgmres_args args;
  int dimx, dimu, dimc, dim_solution;
  double finite_diff, current_time;
  double *state, *solution, *inc_solution, *direction, *b, *ax, *b_ref, *ax_ref;
};


TEST_F(d_inexact_newton_for_zero_horizon_ocp_test, memsize) {
  int expected_size = 0;
  expected_size += 1*d_zero_horizon_ocp_memsize();
  expected_size += 3*(dimu+dimc)*sizeof(double); // incremented_solution, optimality_residual, optimality_residual1;
  expected_size = (expected_size+63)/64*64; // make multiple of typical cache line size
  expected_size += 64; // align to typical cache line size
  EXPECT_EQ(expected_size, d_inexact_newton_for_zero_horizon_ocp_memsize());
}


TEST_F(d_inexact_newton_for_zero_horizon_ocp_test, b_and_ax) {
  args.current_time = current_time;
  args.current_state_ptr = state;
  args.current_solution_ptr = solution;
  double *opt_res = allocate_dvec(dim_solution);
  double *opt_res1 = allocate_dvec(dim_solution);
  hpcgmres_daxpy(dim_solution, finite_diff, direction, solution, inc_solution);
  d_zero_horizon_ocp_compute_optimality_residual(&newton.ocp, current_time, 
                                                 state, solution, opt_res);
  d_zero_horizon_ocp_compute_optimality_residual(&newton.ocp, current_time, 
                                                 state, inc_solution, opt_res1);
  hpcgmres_daxpby(dim_solution, 1/finite_diff-1, opt_res, -1/finite_diff, 
                  opt_res1, b_ref);
  d_inexact_newton_for_zero_horizon_ocp_compute_b(&newton, &args, direction, b);
  for (int i=0; i<dim_solution; ++i) {
    EXPECT_EQ(b_ref[i], b[i]);
  }

  hpcgmres_daxpby(dim_solution, -1/finite_diff, opt_res, 1/finite_diff, 
                  opt_res1, ax_ref);
  d_inexact_newton_for_zero_horizon_ocp_compute_ax(&newton, &args, direction, ax);
  for (int i=0; i<dim_solution; ++i) {
    EXPECT_EQ(ax_ref[i], ax[i]);
  }
  free_dvec(opt_res1);
  free_dvec(opt_res);
}


TEST_F(d_inexact_newton_for_zero_horizon_ocp_test, error_norm) {
  double *opt_res = allocate_dvec(dim_solution);
  d_zero_horizon_ocp_compute_optimality_residual(&newton.ocp, current_time, 
                                                 state, solution, opt_res);
  double err_norm_ref = sqrt(hpcgmres_dvecnrm2(dim_solution, opt_res));
  double err_norm = d_inexact_newton_for_zero_horizon_ocp_get_error_norm(&newton, current_time, state, solution);
  EXPECT_EQ(err_norm_ref, err_norm);
  free_dvec(opt_res);
}


TEST_F(d_inexact_newton_for_zero_horizon_ocp_test, terminal_cost_derivative) {
  double *phix = allocate_dvec(dimx);
  double *phix_ref = allocate_dvec(dimx);
  d_zero_horizon_ocp_compute_terminal_cost_derivative(&newton.ocp, current_time, 
                                                      state, phix_ref);
  d_inexact_newton_for_zero_horizon_ocp_get_terminal_cost_derivative(&newton, current_time, state, phix);
  for (int i=0; i<dimx; ++i) {
    EXPECT_EQ(phix_ref[i], phix[i]);
  }
  free_dvec(phix);
  free_dvec(phix_ref);
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}