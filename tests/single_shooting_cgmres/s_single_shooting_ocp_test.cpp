#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>

#include <gtest/gtest.h>

extern "C" {
#include "s_single_shooting_ocp.h"
#include "s_nmpc_model.h"
#include "s_memory_manager.h"
#include "s_linear_algebra.h"
#include "s_time_varying_smooth_horizon.h"
}


class s_single_shooting_ocp_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand(time(NULL));
    T_f = 1.0;
    alpha = 1.0;
    initial_time = 0.0;
    N = 50;
    s_single_shooting_ocp_create(&ocp, T_f, alpha, initial_time, N);
    s_time_varying_smooth_horizon_create(&horizon, T_f, alpha, initial_time);
    s_nmpc_model_create(&model);
    dimx = s_nmpc_model_dimx();
    dimu = s_nmpc_model_dimu();
    dimc = s_nmpc_model_dimc();
    dimuc = dimu + dimc;
    dim_solution = N * dimuc;
    x_seq = allocate_smat(N, dimx);
    lmd_seq = allocate_smat(N, dimx);
    dx = allocate_svec(dimx);
    x0 = allocate_svec(dimx);
    solution = allocate_svec(dim_solution);
    opt_res = allocate_svec(dim_solution);
    opt_res_ref = allocate_svec(dim_solution);
    for (int i=0; i<dimx; ++i) {
      x0[i] = (float)rand()/RAND_MAX;
    }
    for (int i=0; i<dim_solution; ++i) {
      solution[i] = (float)rand()/RAND_MAX;
    }
  }

  virtual void TearDown() {
    s_single_shooting_ocp_delete(&ocp);
    free_smat(x_seq);
    free_smat(lmd_seq);
    free_svec(dx);
    free_svec(x0);
    free_svec(solution);
    free_svec(opt_res);
    free_svec(opt_res_ref);
  }

  struct s_single_shooting_ocp ocp;
  struct s_time_varying_smooth_horizon horizon;
  struct s_nmpc_model model;
  float T_f, alpha, initial_time;
  int N, dimx, dimu, dimc, dimuc, dim_solution;
  float **x_seq, **lmd_seq, *dx, *x0, *solution, *opt_res, *opt_res_ref;
};


TEST_F(s_single_shooting_ocp_test, memsize) {
  int expect_size = 0;
  expect_size += 2*N*sizeof(float*); 
  expect_size += 2*N*s_nmpc_model_dimx()*sizeof(float); 
  expect_size += s_nmpc_model_dimx()*sizeof(float); 
  expect_size = (expect_size+63)/64*64; // make multiple of typical cache line size
  expect_size += 64; // align to typical cache line size
  EXPECT_EQ(expect_size, s_single_shooting_ocp_memsize(N));
}


TEST_F(s_single_shooting_ocp_test, dim) {
  EXPECT_EQ(4, s_single_shooting_ocp_dimx());
  EXPECT_EQ(2, s_single_shooting_ocp_dimu());
  EXPECT_EQ(1, s_single_shooting_ocp_dimc());
  EXPECT_EQ(4, ocp.model.dimx);
  EXPECT_EQ(2, ocp.model.dimu);
  EXPECT_EQ(1, ocp.model.dimc);
}


TEST_F(s_single_shooting_ocp_test, compute_optimality_residual) {
  float current_time = fabs((float)rand()/RAND_MAX);
  float horizon_length = s_time_varying_smooth_horizon_get_length(&horizon, 
                                                                   current_time);
  float delta_tau = horizon_length / N;
  s_nmpc_model_f(&model, current_time, x0, solution, dx);
  hpcgmres_saxpy(dimx, delta_tau, dx, x0, x_seq[0]);
  float tau = current_time + delta_tau;
  for (int ii=1; ii<N; ++ii, tau+=delta_tau) {
    s_nmpc_model_f(&model, tau, x_seq[ii-1], &(solution[ii*dimuc]), dx);
    hpcgmres_saxpy(dimx, delta_tau, dx, x_seq[ii-1], x_seq[ii]);
  }
  s_nmpc_model_phix(&model, tau, x_seq[N-1], lmd_seq[N-1]);
  for (int ii=N-1; ii>=1; --ii, tau-=delta_tau) {
    s_nmpc_model_hx(&model, tau, x_seq[ii-1], &(solution[ii*dimuc]), 
                    lmd_seq[ii], dx);
    hpcgmres_saxpy(dimx, delta_tau, dx, lmd_seq[ii], lmd_seq[ii-1]);
  }
  tau = current_time;
  s_nmpc_model_hu(&model, tau, x0, solution, lmd_seq[0], opt_res);
  for (int ii=1; ii<N; ++ii, tau+=delta_tau) {
    s_nmpc_model_hu(&model, tau, x_seq[ii-1], &(solution[ii*dimuc]), lmd_seq[ii],
                    &(opt_res[ii*dimuc]));
  }

  s_single_shooting_ocp_compute_optimality_residual(&ocp, current_time, x0, 
                                                    solution, opt_res_ref);
  for (int i=0; i<dim_solution; ++i) {
    EXPECT_EQ(opt_res[i], opt_res_ref[i]);
  }
}


TEST_F(s_single_shooting_ocp_test, predict_state) {
  float current_time = fabs((float)rand()/RAND_MAX);
  float prediction_length = fabs((float)rand()/RAND_MAX);
  s_nmpc_model_f(&model, current_time, x0, solution, dx);
  hpcgmres_saxpy(dimx, prediction_length, dx, x0, x_seq[0]);
  s_single_shooting_ocp_predict_state_from_solution(&ocp, current_time, x0,
                                                    solution, prediction_length, 
                                                    x_seq[1]);
  for (int i=0; i<dimx; ++i) {
    EXPECT_EQ(x_seq[0][i], x_seq[1][i]);
  }
}


TEST_F(s_single_shooting_ocp_test, reset_horizon) {
  float current_time = fabs((float)rand()/RAND_MAX);
  s_single_shooting_ocp_reset_horizon_length(&ocp, current_time);
  float horizon_length = s_time_varying_smooth_horizon_get_length(&ocp.horizon, 
                                                                   current_time);
  EXPECT_EQ(horizon_length, 0);
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}