
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <blasfeo.h>
#include <iostream>

#include <gtest/gtest.h>

extern "C" {
#include "s_inexact_newton_for_zero_horizon_ocp_mfgmres_args.h"
}


class s_inexact_newton_for_zero_horizon_ocp_mfgmres_args_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    dim = 10;
    blasfeo_allocate_svec(dim, &state);
    blasfeo_allocate_svec(dim, &sol);
  }

  virtual void TearDown() {
    blasfeo_free_svec(&state);
    blasfeo_free_svec(&sol);
  }

  int dim;
  float time_param;
  struct blasfeo_svec state, sol;
  struct s_inexact_newton_for_zero_horizon_ocp_mfgmres_args args;
};


TEST_F(s_inexact_newton_for_zero_horizon_ocp_mfgmres_args_test, ptr) {
  time_param = (float)rand()/RAND_MAX;
  for (int i=0; i<dim; ++i) {
    blasfeo_svecin1((float)rand()/RAND_MAX, &state, i);
  }
  for (int i=0; i<dim; ++i) {
    blasfeo_svecin1((float)rand()/RAND_MAX, &sol, i);
  }
  args.current_time = time_param;
  args.current_state_ptr = &state;
  args.current_solution_ptr = &sol;

  EXPECT_EQ(time_param, args.current_time);
  for (int i=0; i<dim; ++i) {
    EXPECT_EQ(blasfeo_svecex1(&state, i), 
              blasfeo_svecex1(args.current_state_ptr, i));
  }
  for (int i=0; i<dim; ++i) {
    EXPECT_EQ(blasfeo_svecex1(&sol, i), 
              blasfeo_svecex1(args.current_solution_ptr, i));
  }
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}