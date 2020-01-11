#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>

#include <gtest/gtest.h>

extern "C" {
#include "d_zero_horizon_ocp.h"
#include "d_nmpc_model.h"
}


class d_zero_horizon_ocp_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    d_zero_horizon_ocp_create(&ocp);
  }

  virtual void TearDown() {
    d_zero_horizon_ocp_delete(&ocp);
  }

  struct d_zero_horizon_ocp ocp;
};


TEST_F(d_zero_horizon_ocp_test, memsize) {
  int expect_size = 0;
  expect_size += d_nmpc_model_dimx()*sizeof(double); 
  expect_size = (expect_size+63)/64*64; 
  expect_size += 64; 
  EXPECT_EQ(expect_size, d_zero_horizon_ocp_memsize());
}


TEST_F(d_zero_horizon_ocp_test, dim) {
  EXPECT_EQ(4, d_zero_horizon_ocp_dimx());
  EXPECT_EQ(2, d_zero_horizon_ocp_dimu());
  EXPECT_EQ(1, d_zero_horizon_ocp_dimc());
  EXPECT_EQ(4, ocp.model.dimx);
  EXPECT_EQ(2, ocp.model.dimu);
  EXPECT_EQ(1, ocp.model.dimc);
}


TEST_F(d_zero_horizon_ocp_test, compute_optimality_residual) {
  double current_time;
  double x[ocp.model.dimx];
  double u[ocp.model.dimu+ocp.model.dimu];
  double lmd[ocp.model.dimx];
  double hu[ocp.model.dimu+ocp.model.dimu];
  double hu_ref[ocp.model.dimu+ocp.model.dimu];

  current_time = (double)rand()/RAND_MAX;
  for (int i=0; i<ocp.model.dimx; ++i) {
    x[i] = (double)rand()/RAND_MAX;
  }
  for (int i=0; i<ocp.model.dimu+ocp.model.dimu; ++i) {
    u[i] = (double)rand()/RAND_MAX;
  }
  d_nmpc_model_phix(&ocp.model, current_time, x, lmd);
  d_nmpc_model_hu(&ocp.model, current_time, x, u, lmd, hu_ref);
  d_zero_horizon_ocp_compute_optimality_residual(&ocp, current_time, x, u, hu);
  for (int i=0; i<ocp.model.dimu+ocp.model.dimc; ++i) {
    EXPECT_EQ(hu[i], hu_ref[i]);
  }
}


TEST_F(d_zero_horizon_ocp_test, terminal_cost_derivative) {
  double current_time;
  double x[ocp.model.dimx];
  double u[ocp.model.dimu+ocp.model.dimu];
  double phix[ocp.model.dimx];
  double phix_ref[ocp.model.dimx];
  current_time = (double)rand()/RAND_MAX;
  for (int i=0; i<ocp.model.dimx; ++i) {
    x[i] = (double)rand()/RAND_MAX;
  }
  d_nmpc_model_phix(&ocp.model, current_time, x, phix_ref);
  d_zero_horizon_ocp_compute_terminal_cost_derivative(&ocp, current_time, x,
                                                      phix);
  for (int i=0; i<ocp.model.dimu+ocp.model.dimc; ++i) {
    EXPECT_EQ(phix[i], phix_ref[i]);
  }
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}