#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>

#include <gtest/gtest.h>

extern "C" {
#include "d_single_shooting_continuation.h"
#include "d_single_shooting_continuation_mfgmres_args.h"
#include "d_single_shooting_ocp.h"
#include "d_memory_manager.h"
#include "d_linear_algebra.h"
}


class d_single_shooting_continuation_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand(time(NULL));
    T_f = 1.0;
    alpha = 1.0;
    initial_time = 0.0;
    N = 50;
    finite_diff = 1.0e-08;
    zeta = 1000;
    d_single_shooting_continuation_create(&continuation, T_f, alpha, 
                                          initial_time, N, finite_diff, zeta);
    d_single_shooting_ocp_create(&ocp, T_f, alpha, initial_time, N);
    dimx = d_single_shooting_ocp_dimx();
    dimu = d_single_shooting_ocp_dimu();
    dimc = d_single_shooting_ocp_dimc();
    dimuc = dimu + dimc;
    dim_solution = N * dimuc;
    x0 = allocate_dvec(dimx);
    x1 = allocate_dvec(dimx);
    solution = allocate_dvec(dim_solution);
    direction = allocate_dvec(dim_solution);
    inc_solution = allocate_dvec(dim_solution);
    b = allocate_dvec(dim_solution);
    b_ref = allocate_dvec(dim_solution);
    ax = allocate_dvec(dim_solution);
    ax_ref = allocate_dvec(dim_solution);
    opt = allocate_dvec(dim_solution);
    opt1 = allocate_dvec(dim_solution);
    opt2 = allocate_dvec(dim_solution);
    for (int i=0; i<dimx; ++i) {
      x0[i] = (double)rand()/RAND_MAX;
    }
    for (int i=0; i<dim_solution; ++i) {
      solution[i] = (double)rand()/RAND_MAX;
    }
    for (int i=0; i<dim_solution; ++i) {
      direction[i] = (double)rand()/RAND_MAX;
    }
  }

  virtual void TearDown() {
    d_single_shooting_continuation_delete(&continuation);
    d_single_shooting_ocp_delete(&ocp);
    free_dvec(x0);
    free_dvec(x1);
    free_dvec(solution);
    free_dvec(direction);
    free_dvec(inc_solution);
    free_dvec(b);
    free_dvec(b_ref);
    free_dvec(ax);
    free_dvec(ax_ref);
    free_dvec(opt);
    free_dvec(opt1);
    free_dvec(opt2);
  }

  struct d_single_shooting_continuation continuation;
  struct d_single_shooting_continuation_mfgmres_args args;
  struct d_single_shooting_ocp ocp;
  double T_f, alpha, initial_time, finite_diff, zeta;
  int N, dimx, dimu, dimc, dimuc, dim_solution;
  double *x0, *x1, *solution, *direction, *inc_solution, *b, *b_ref, 
         *ax, *ax_ref, *opt, *opt1, *opt2;
};


TEST_F(d_single_shooting_continuation_test, memsize) {
  int expect_size = 0;
  expect_size += d_single_shooting_ocp_memsize(N);
  expect_size += dimx*sizeof(double); 
  expect_size += 4*dim_solution*sizeof(double); 
  expect_size = (expect_size+63)/64*64; // make multiple of typical cache line size
  expect_size += 64; // align to typical cache line size
  EXPECT_EQ(expect_size, d_single_shooting_continuation_memsize(N));
}


TEST_F(d_single_shooting_continuation_test, integrate_solution) {
  double length = fabs((double)rand()/RAND_MAX);
  hpcgmres_daxpy(dim_solution, length, direction, solution, inc_solution);
  d_single_shooting_continuation_integrate_solution(&continuation, solution, 
                                                    direction, length);
  for (int i=0; i<dim_solution; ++i) {
    EXPECT_EQ(inc_solution[i], solution[i]);
  }
}


TEST_F(d_single_shooting_continuation_test, compute_b_and_ax) {
  double current_time = fabs((double)rand()/RAND_MAX);
  args.current_time = current_time;
  args.current_state_ptr = x0;
  args.current_solution_ptr = solution;
  d_single_shooting_continuation_compute_b(&continuation, &args, direction, b);
  d_single_shooting_continuation_compute_ax(&continuation, &args, direction, ax);

  double incremented_time = current_time + finite_diff;
  d_single_shooting_ocp_predict_state_from_solution(&ocp, current_time, x0, 
                                                    solution, finite_diff, x1);
  hpcgmres_daxpy(dim_solution, finite_diff, direction, solution, inc_solution);
  d_single_shooting_ocp_compute_optimality_residual(&ocp, current_time, x0, 
                                                    solution, opt);
  d_single_shooting_ocp_compute_optimality_residual(&ocp, incremented_time, x1, 
                                                    solution, opt1);
  d_single_shooting_ocp_compute_optimality_residual(&ocp, incremented_time, x1, 
                                                    inc_solution, opt2);
  hpcgmres_daxpby(dim_solution, 1/finite_diff-zeta, opt, -1/finite_diff, opt2, 
                  b_ref); 
  for (int i=0; i<dim_solution; ++i) {
    EXPECT_EQ(b_ref[i], b[i]);
  }

  hpcgmres_daxpy(dim_solution, finite_diff, direction, solution, inc_solution);
  d_single_shooting_ocp_compute_optimality_residual(&ocp, incremented_time, x1, 
                                                    inc_solution, opt2);
  hpcgmres_daxpby(dim_solution, 1/finite_diff, opt2, -1/finite_diff, opt1, 
                  ax_ref); 
  for (int i=0; i<dim_solution; ++i) {
    EXPECT_EQ(ax_ref[i], ax[i]);
  }
}


TEST_F(d_single_shooting_continuation_test, err_norm) {
  double current_time = fabs((double)rand()/RAND_MAX);

  d_single_shooting_ocp_compute_optimality_residual(&ocp, current_time, x0, 
                                                    solution, opt);
  double err_ref = sqrt(hpcgmres_dvecnrm2(dim_solution, opt));
  double err  = d_single_shooting_continuation_get_error_norm(&continuation, 
                                                              current_time, x0, 
                                                              solution);
  EXPECT_EQ(err_ref, err);
}



TEST_F(d_single_shooting_continuation_test, reset_horizon) {
  double current_time = fabs((double)rand()/RAND_MAX);
  d_single_shooting_continuation_reset_horizon_length(&continuation, 
                                                      current_time);
  double horizon_length = d_time_varying_smooth_horizon_get_length(&continuation.ocp.horizon, 
                                                                   current_time);
  EXPECT_EQ(horizon_length, 0);
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}