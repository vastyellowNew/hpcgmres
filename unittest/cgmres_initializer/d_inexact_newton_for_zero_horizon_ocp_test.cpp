#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <blasfeo.h>

#include <gtest/gtest.h>

extern "C" {
#include "d_zero_horizon_ocp.h"
#include "d_inexact_newton_for_zero_horizon_ocp.h"
#include "d_inexact_newton_for_zero_horizon_ocp_mfgmres_args.h"
}


class d_inexact_newton_for_zero_horizon_ocp_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    int memsize = d_inexact_newton_for_zero_horizon_ocp_memsize();
    memory = malloc(memsize);
    finite_diff = 1.0e-08;
    d_inexact_newton_for_zero_horizon_ocp_create(&newton, finite_diff, memory);
    dimx = d_inexact_newton_for_zero_horizon_ocp_dimx();
    dimu = d_inexact_newton_for_zero_horizon_ocp_dimu();
    dimc = d_inexact_newton_for_zero_horizon_ocp_dimc();
    dim_solution = dimu + dimc;
    blasfeo_allocate_dvec(dimx, &state);
    blasfeo_allocate_dvec(dim_solution, &solution);
    blasfeo_allocate_dvec(dim_solution, &direction);
    blasfeo_allocate_dvec(dim_solution, &increment_solution);
    blasfeo_allocate_dvec(dim_solution, &opt_res);
    blasfeo_allocate_dvec(dim_solution, &opt_res1);
    blasfeo_allocate_dvec(dim_solution, &b);
    blasfeo_allocate_dvec(dim_solution, &b_ref);
    blasfeo_allocate_dvec(dim_solution, &ax);
    blasfeo_allocate_dvec(dim_solution, &ax_ref);
    blasfeo_allocate_dvec(dimx, &phix);
    blasfeo_allocate_dvec(dimx, &phix_ref);

    current_time = (double)rand()/RAND_MAX;
    for (int i=0; i<dimx; ++i) {
      blasfeo_dvecin1((double)rand()/RAND_MAX, &state, i);
    }
    for (int i=0; i<dim_solution; ++i) {
      blasfeo_dvecin1((double)rand()/RAND_MAX, &solution, i);
      blasfeo_dvecin1((double)rand()/RAND_MAX, &direction, i);
    }
  }

  virtual void TearDown() {
    free(memory);
    blasfeo_free_dvec(&state);
    blasfeo_free_dvec(&solution);
    blasfeo_free_dvec(&direction);
    blasfeo_free_dvec(&increment_solution);
    blasfeo_free_dvec(&opt_res);
    blasfeo_free_dvec(&opt_res1);
    blasfeo_free_dvec(&b);
    blasfeo_free_dvec(&b_ref);
    blasfeo_free_dvec(&ax);
    blasfeo_free_dvec(&ax_ref);
    blasfeo_free_dvec(&phix);
    blasfeo_free_dvec(&phix_ref);
  }

  void *memory;
  double finite_diff, current_time;
  int dimx, dimu, dimc, dim_solution;
  struct d_inexact_newton_for_zero_horizon_ocp newton;
  struct d_inexact_newton_for_zero_horizon_ocp_mfgmres_args args;
  struct blasfeo_dvec state, solution, direction, increment_solution, opt_res, 
                      opt_res1, b, b_ref, ax, ax_ref, phix, phix_ref;
};


TEST_F(d_inexact_newton_for_zero_horizon_ocp_test, memsize) {
  int expected_size = 0;
  expected_size += 1*d_zero_horizon_ocp_memsize();
  expected_size += 3*sizeof(struct blasfeo_dvec*); // incremented_solution, optimality_residual, optimality_residual1;
  int dimu = d_zero_horizon_ocp_dimu();
  int dimc = d_zero_horizon_ocp_dimc();
  expected_size += 1*blasfeo_memsize_dvec(dimu+dimc); // incremented_solution
  expected_size += 1*blasfeo_memsize_dvec(dimu+dimc); // optimality_residual
  expected_size += 1*blasfeo_memsize_dvec(dimu+dimc); // optimality_residual1
  expected_size = (expected_size+63)/64*64; // make multiple of typical cache line size
  expected_size += 64; // align to typical cache line size
  EXPECT_EQ(expected_size, d_inexact_newton_for_zero_horizon_ocp_memsize());
}


TEST_F(d_inexact_newton_for_zero_horizon_ocp_test, b_and_ax) {
  args.current_time = current_time;
  args.current_state_ptr = &state;
  args.current_solution_ptr = &solution;
  blasfeo_daxpy(dim_solution, finite_diff, &direction, 0, &solution, 0, 
                &increment_solution, 0);
  d_zero_horizon_ocp_compute_optimality_residual(&newton.ocp, current_time, 
                                                 &state, &solution, &opt_res);
  d_zero_horizon_ocp_compute_optimality_residual(&newton.ocp, current_time, 
                                                 &state, &increment_solution, 
                                                 &opt_res1);
  blasfeo_daxpby(dim_solution, 1/finite_diff-1, &opt_res, 0, -1/finite_diff, 
                 &opt_res1, 0, &b_ref, 0);
  d_inexact_newton_for_zero_horizon_ocp_compute_b(&newton, &args, &direction, 
                                                  &b);
  for (int i=0; i<dim_solution; ++i) {
    EXPECT_EQ(blasfeo_dvecex1(&b_ref, i), blasfeo_dvecex1(&b, i));
  }
    for (int i=0; i<dim_solution; ++i) {
      blasfeo_dvecin1((double)rand()/RAND_MAX, &direction, i);
    }
  blasfeo_daxpy(dim_solution, finite_diff, &direction, 0, &solution, 0, 
                &increment_solution, 0);
  d_zero_horizon_ocp_compute_optimality_residual(&newton.ocp, current_time, 
                                                 &state, &increment_solution, 
                                                 &opt_res1);
  blasfeo_daxpby(dim_solution, 1/finite_diff-1, &opt_res, 0, 1/finite_diff, 
                 &opt_res1, 0, &ax_ref, 0);
  d_inexact_newton_for_zero_horizon_ocp_compute_ax(&newton, &args, &direction, 
                                                  &ax);
  for (int i=0; i<dim_solution; ++i) {
    EXPECT_EQ(blasfeo_dvecex1(&ax_ref, i), blasfeo_dvecex1(&ax, i));
  }
}


TEST_F(d_inexact_newton_for_zero_horizon_ocp_test, error_norm) {
  d_zero_horizon_ocp_compute_optimality_residual(&newton.ocp, current_time, 
                                                 &state, &solution, &opt_res);
  double err_norm_ref = sqrt(blasfeo_ddot(dim_solution, &opt_res, 0, &opt_res, 0));
  double err_norm = d_inexact_newton_for_zero_horizon_ocp_get_error_norm(&newton, current_time, &state, &solution);
  EXPECT_EQ(err_norm_ref, err_norm);
}


TEST_F(d_inexact_newton_for_zero_horizon_ocp_test, terminal_cost_derivative) {
  d_zero_horizon_ocp_compute_terminal_cost_derivative(&newton.ocp, current_time, 
                                                      &state, &phix_ref);
  d_inexact_newton_for_zero_horizon_ocp_get_terminal_cost_derivative(&newton, current_time, &state, &phix);
  for (int i=0; i<dimx; ++i) {
    EXPECT_EQ(blasfeo_dvecex1(&phix_ref, i), blasfeo_dvecex1(&phix, i));
  }
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}