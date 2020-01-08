#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <blasfeo.h>
#include <iostream>

#include <gtest/gtest.h>

extern "C" {
#include "d_zero_horizon_ocp.h"
#include "d_nmpc_model.h"
}


class d_zero_horizon_ocp_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    memory = malloc(d_zero_horizon_ocp_memsize());
    d_zero_horizon_ocp_create(&ocp, memory);
    blasfeo_allocate_dvec(ocp.model.dimx, &state);
    blasfeo_allocate_dvec(ocp.model.dimu+ocp.model.dimc, &sol);
    blasfeo_allocate_dvec(ocp.model.dimu+ocp.model.dimc, &opt_res);
    blasfeo_allocate_dvec(ocp.model.dimx, &terminal_cost_der);
  }

  virtual void TearDown() {
    free(memory);
    blasfeo_free_dvec(&state);
    blasfeo_free_dvec(&sol);
    blasfeo_free_dvec(&opt_res);
    blasfeo_free_dvec(&terminal_cost_der);
  }

  void *memory;
  struct d_zero_horizon_ocp ocp;
  struct blasfeo_dvec state, sol, opt_res, terminal_cost_der;
};


TEST_F(d_zero_horizon_ocp_test, memsize) {
  int expect_size = 0;
  expect_size += 1*sizeof(struct blasfeo_dvec*);
  expect_size += 1*blasfeo_memsize_dvec(d_nmpc_model_dimx()); 
  expect_size = (expect_size+63)/64*64; 
  expect_size += 64; 
  EXPECT_EQ(expect_size, d_zero_horizon_ocp_memsize());
}


TEST_F(d_zero_horizon_ocp_test, dim) {
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
  current_time = (double)rand()/RAND_MAX;
  for (int i=0; i<ocp.model.dimx; ++i) {
    x[i] = (double)rand()/RAND_MAX;
    blasfeo_dvecin1(x[i], &state, i);
  }
  for (int i=0; i<ocp.model.dimu+ocp.model.dimu; ++i) {
    u[i] = (double)rand()/RAND_MAX;
    blasfeo_dvecin1(u[i], &sol, i);
  }
  d_nmpc_model_phix(&ocp.model, current_time, x, lmd);
  d_nmpc_model_hu(&ocp.model, current_time, x, u, lmd, hu);
  d_zero_horizon_ocp_compute_optimality_residual(&ocp, current_time, &state, 
                                                 &sol, &opt_res);
  for (int i=0; i<ocp.model.dimu+ocp.model.dimc; ++i) {
    EXPECT_EQ(hu[i], blasfeo_dvecex1(&opt_res, i));
  }
}


TEST_F(d_zero_horizon_ocp_test, terminal_cost_derivative) {
  double current_time;
  double x[ocp.model.dimx];
  double u[ocp.model.dimu+ocp.model.dimu];
  double lmd[ocp.model.dimx];
  current_time = (double)rand()/RAND_MAX;
  for (int i=0; i<ocp.model.dimx; ++i) {
    x[i] = (double)rand()/RAND_MAX;
    blasfeo_dvecin1(x[i], &state, i);
  }
  d_nmpc_model_phix(&ocp.model, current_time, x, lmd);
  d_zero_horizon_ocp_compute_terminal_cost_derivative(&ocp, current_time, 
                                                      &state, &terminal_cost_der);
  for (int i=0; i<ocp.model.dimu+ocp.model.dimc; ++i) {
    EXPECT_EQ(lmd[i], blasfeo_dvecex1(&terminal_cost_der, i));
  }
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}