
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <blasfeo.h>
#include <iostream>

#include <gtest/gtest.h>

extern "C" {
#include "s_zero_horizon_ocp.h"
#include "s_nmpc_model.h"
}


class s_zero_horizon_ocp_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    memory = malloc(s_zero_horizon_ocp_memsize());
    s_zero_horizon_ocp_create(&ocp, memory);
    blasfeo_allocate_svec(ocp.model.dimx, &state);
    blasfeo_allocate_svec(ocp.model.dimu+ocp.model.dimc, &sol);
    blasfeo_allocate_svec(ocp.model.dimu+ocp.model.dimc, &opt_res);
    blasfeo_allocate_svec(ocp.model.dimx, &terminal_cost_der);
  }

  virtual void TearDown() {
    free(memory);
    blasfeo_free_svec(&state);
    blasfeo_free_svec(&sol);
    blasfeo_free_svec(&opt_res);
    blasfeo_free_svec(&terminal_cost_der);
  }

  void *memory;
  struct s_zero_horizon_ocp ocp;
  struct blasfeo_svec state, sol, opt_res, terminal_cost_der;
};


TEST_F(s_zero_horizon_ocp_test, memsize) {
  int expect_size = 0;
  expect_size += 1*sizeof(struct blasfeo_svec*);
  expect_size += 1*blasfeo_memsize_svec(s_nmpc_model_dimx()); 
  expect_size = (expect_size+63)/64*64; 
  expect_size += 64; 
  EXPECT_EQ(expect_size, s_zero_horizon_ocp_memsize());
}


TEST_F(s_zero_horizon_ocp_test, dim) {
  EXPECT_EQ(4, ocp.model.dimx);
  EXPECT_EQ(2, ocp.model.dimu);
  EXPECT_EQ(1, ocp.model.dimc);
}


TEST_F(s_zero_horizon_ocp_test, compute_optimality_residual) {
  float current_time;
  float x[ocp.model.dimx];
  float u[ocp.model.dimu+ocp.model.dimu];
  float lmd[ocp.model.dimx];
  float hu[ocp.model.dimu+ocp.model.dimu];
  current_time = (float)rand()/RAND_MAX;
  for (int i=0; i<ocp.model.dimx; ++i) {
    x[i] = (float)rand()/RAND_MAX;
    blasfeo_svecin1(x[i], &state, i);
  }
  for (int i=0; i<ocp.model.dimu+ocp.model.dimu; ++i) {
    u[i] = (float)rand()/RAND_MAX;
    blasfeo_svecin1(u[i], &sol, i);
  }
  s_nmpc_model_phix(&ocp.model, current_time, x, lmd);
  s_nmpc_model_hu(&ocp.model, current_time, x, u, lmd, hu);
  s_zero_horizon_ocp_compute_optimality_residual(&ocp, current_time, &state, 
                                                 &sol, &opt_res);
  for (int i=0; i<ocp.model.dimu+ocp.model.dimc; ++i) {
    EXPECT_EQ(hu[i], blasfeo_svecex1(&opt_res, i));
  }
}


TEST_F(s_zero_horizon_ocp_test, terminal_cost_derivative) {
  float current_time;
  float x[ocp.model.dimx];
  float u[ocp.model.dimu+ocp.model.dimu];
  float lmd[ocp.model.dimx];
  current_time = (float)rand()/RAND_MAX;
  for (int i=0; i<ocp.model.dimx; ++i) {
    x[i] = (float)rand()/RAND_MAX;
    blasfeo_svecin1(x[i], &state, i);
  }
  s_nmpc_model_phix(&ocp.model, current_time, x, lmd);
  s_zero_horizon_ocp_compute_terminal_cost_derivative(&ocp, current_time, 
                                                      &state, &terminal_cost_der);
  for (int i=0; i<ocp.model.dimu+ocp.model.dimc; ++i) {
    EXPECT_EQ(lmd[i], blasfeo_svecex1(&terminal_cost_der, i));
  }
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}