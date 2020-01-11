#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>

#include <gtest/gtest.h>

extern "C" {
#include "s_zero_horizon_ocp.h"
#include "s_nmpc_model.h"
}


class s_zero_horizon_ocp_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    s_zero_horizon_ocp_create(&ocp);
  }

  virtual void TearDown() {
    s_zero_horizon_ocp_delete(&ocp);
  }

  struct s_zero_horizon_ocp ocp;
};


TEST_F(s_zero_horizon_ocp_test, memsize) {
  int expect_size = 0;
  expect_size += s_nmpc_model_dimx()*sizeof(float); 
  expect_size = (expect_size+63)/64*64; 
  expect_size += 64; 
  EXPECT_EQ(expect_size, s_zero_horizon_ocp_memsize());
}


TEST_F(s_zero_horizon_ocp_test, dim) {
  EXPECT_EQ(4, s_zero_horizon_ocp_dimx());
  EXPECT_EQ(2, s_zero_horizon_ocp_dimu());
  EXPECT_EQ(1, s_zero_horizon_ocp_dimc());
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
  float hu_ref[ocp.model.dimu+ocp.model.dimu];

  current_time = (float)rand()/RAND_MAX;
  for (int i=0; i<ocp.model.dimx; ++i) {
    x[i] = (float)rand()/RAND_MAX;
  }
  for (int i=0; i<ocp.model.dimu+ocp.model.dimu; ++i) {
    u[i] = (float)rand()/RAND_MAX;
  }
  s_nmpc_model_phix(&ocp.model, current_time, x, lmd);
  s_nmpc_model_hu(&ocp.model, current_time, x, u, lmd, hu_ref);
  s_zero_horizon_ocp_compute_optimality_residual(&ocp, current_time, x, u, hu);
  for (int i=0; i<ocp.model.dimu+ocp.model.dimc; ++i) {
    EXPECT_EQ(hu[i], hu_ref[i]);
  }
}


TEST_F(s_zero_horizon_ocp_test, terminal_cost_derivative) {
  float current_time;
  float x[ocp.model.dimx];
  float u[ocp.model.dimu+ocp.model.dimu];
  float phix[ocp.model.dimx];
  float phix_ref[ocp.model.dimx];
  current_time = (float)rand()/RAND_MAX;
  for (int i=0; i<ocp.model.dimx; ++i) {
    x[i] = (float)rand()/RAND_MAX;
  }
  s_nmpc_model_phix(&ocp.model, current_time, x, phix_ref);
  s_zero_horizon_ocp_compute_terminal_cost_derivative(&ocp, current_time, x,
                                                      phix);
  for (int i=0; i<ocp.model.dimu+ocp.model.dimc; ++i) {
    EXPECT_EQ(phix[i], phix_ref[i]);
  }
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}