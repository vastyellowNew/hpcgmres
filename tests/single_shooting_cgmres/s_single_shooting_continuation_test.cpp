#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>

#include <gtest/gtest.h>

extern "C" {
#include "s_single_shooting_continuation.h"
#include "s_single_shooting_continuation_mfgmres_args.h"
#include "s_single_shooting_ocp.h"
#include "s_memory_manager.h"
#include "s_linear_algebra.h"
}


class s_single_shooting_continuation_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand(time(NULL));
    T_f = 1.0;
    alpha = 1.0;
    initial_time = 0.0;
    N = 50;
    finite_diff = 1.0e-08;
    zeta = 1000;
    s_single_shooting_continuation_create(&continuation, T_f, alpha, 
                                          initial_time, N, finite_diff, zeta);
    s_single_shooting_ocp_create(&ocp, T_f, alpha, initial_time, N);
    dimx = s_single_shooting_ocp_dimx();
    dimu = s_single_shooting_ocp_dimu();
    dimc = s_single_shooting_ocp_dimc();
    dimuc = dimu + dimc;
    dim_solution = N * dimuc;
    x0 = allocate_svec(dimx);
    x1 = allocate_svec(dimx);
    solution = allocate_svec(dim_solution);
    direction = allocate_svec(dim_solution);
    inc_solution = allocate_svec(dim_solution);
    b = allocate_svec(dim_solution);
    b_ref = allocate_svec(dim_solution);
    ax = allocate_svec(dim_solution);
    ax_ref = allocate_svec(dim_solution);
    opt = allocate_svec(dim_solution);
    opt1 = allocate_svec(dim_solution);
    opt2 = allocate_svec(dim_solution);
    for (int i=0; i<dimx; ++i) {
      x0[i] = (float)rand()/RAND_MAX;
    }
    for (int i=0; i<dim_solution; ++i) {
      solution[i] = (float)rand()/RAND_MAX;
    }
    for (int i=0; i<dim_solution; ++i) {
      direction[i] = (float)rand()/RAND_MAX;
    }
  }

  virtual void TearDown() {
    s_single_shooting_continuation_delete(&continuation);
    s_single_shooting_ocp_delete(&ocp);
    free_svec(x0);
    free_svec(x1);
    free_svec(solution);
    free_svec(direction);
    free_svec(inc_solution);
    free_svec(b);
    free_svec(b_ref);
    free_svec(ax);
    free_svec(ax_ref);
    free_svec(opt);
    free_svec(opt1);
    free_svec(opt2);
  }

  struct s_single_shooting_continuation continuation;
  struct s_single_shooting_continuation_mfgmres_args args;
  struct s_single_shooting_ocp ocp;
  float T_f, alpha, initial_time, finite_diff, zeta;
  int N, dimx, dimu, dimc, dimuc, dim_solution;
  float *x0, *x1, *solution, *direction, *inc_solution, *b, *b_ref, 
         *ax, *ax_ref, *opt, *opt1, *opt2;
};


TEST_F(s_single_shooting_continuation_test, memsize) {
  int expect_size = 0;
  expect_size += s_single_shooting_ocp_memsize(N);
  expect_size += dimx*sizeof(float); 
  expect_size += 4*dim_solution*sizeof(float); 
  expect_size = (expect_size+63)/64*64; // make multiple of typical cache line size
  expect_size += 64; // align to typical cache line size
  EXPECT_EQ(expect_size, s_single_shooting_continuation_memsize(N));
}


TEST_F(s_single_shooting_continuation_test, integrate_solution) {
  float length = fabs((float)rand()/RAND_MAX);
  hpcgmres_saxpy(dim_solution, length, direction, solution, inc_solution);
  s_single_shooting_continuation_integrate_solution(&continuation, solution, 
                                                    direction, length);
  for (int i=0; i<dim_solution; ++i) {
    EXPECT_EQ(inc_solution[i], solution[i]);
  }
}


TEST_F(s_single_shooting_continuation_test, compute_b_ans_ax) {
  float current_time = fabs((float)rand()/RAND_MAX);
  args.current_time = current_time;
  args.current_state_ptr = x0;
  args.current_solution_ptr = solution;
  s_single_shooting_continuation_compute_b(&continuation, &args, direction, b);
  s_single_shooting_continuation_compute_ax(&continuation, &args, direction, ax);

  float incrementes_time = current_time + finite_diff;
  s_single_shooting_ocp_predict_state_from_solution(&ocp, current_time, x0, 
                                                    solution, finite_diff, x1);
  hpcgmres_saxpy(dim_solution, finite_diff, direction, solution, inc_solution);
  s_single_shooting_ocp_compute_optimality_residual(&ocp, current_time, x0, 
                                                    solution, opt);
  s_single_shooting_ocp_compute_optimality_residual(&ocp, incrementes_time, x1, 
                                                    solution, opt1);
  s_single_shooting_ocp_compute_optimality_residual(&ocp, incrementes_time, x1, 
                                                    inc_solution, opt2);
  hpcgmres_saxpby(dim_solution, 1/finite_diff-zeta, opt, -1/finite_diff, opt2, 
                  b_ref); 
  for (int i=0; i<dim_solution; ++i) {
    EXPECT_EQ(b_ref[i], b[i]);
  }

  hpcgmres_saxpy(dim_solution, finite_diff, direction, solution, inc_solution);
  s_single_shooting_ocp_compute_optimality_residual(&ocp, incrementes_time, x1, 
                                                    inc_solution, opt2);
  hpcgmres_saxpby(dim_solution, 1/finite_diff, opt2, -1/finite_diff, opt1, 
                  ax_ref); 
  for (int i=0; i<dim_solution; ++i) {
    EXPECT_EQ(ax_ref[i], ax[i]);
  }
}


TEST_F(s_single_shooting_continuation_test, err_norm) {
  float current_time = fabs((float)rand()/RAND_MAX);

  s_single_shooting_ocp_compute_optimality_residual(&ocp, current_time, x0, 
                                                    solution, opt);
  float err_ref = sqrt(hpcgmres_svecnrm2(dim_solution, opt));
  float err  = s_single_shooting_continuation_get_error_norm(&continuation, 
                                                              current_time, x0, 
                                                              solution);
  EXPECT_EQ(err_ref, err);
}



TEST_F(s_single_shooting_continuation_test, reset_horizon) {
  float current_time = fabs((float)rand()/RAND_MAX);
  s_single_shooting_continuation_reset_horizon_length(&continuation, 
                                                      current_time);
  float horizon_length = s_time_varying_smooth_horizon_get_length(&continuation.ocp.horizon, 
                                                                   current_time);
  EXPECT_EQ(horizon_length, 0);
}


TEST_F(s_single_shooting_continuation_test, dim) {
  EXPECT_EQ(dimx, s_single_shooting_continuation_dimx());
  EXPECT_EQ(dimu, s_single_shooting_continuation_dimu());
  EXPECT_EQ(dimc, s_single_shooting_continuation_dimc());
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}