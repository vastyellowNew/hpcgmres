#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>

#include <gtest/gtest.h>

extern "C" {
#include "d_inexact_newton_for_zero_horizon_ocp_mfgmres_args.h"
#include "d_memory_manager.h"
}


class d_inexact_newton_for_zero_horizon_ocp_mfgmres_args_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    dim = 50;
    state = allocate_dvec(dim);
    state_ref = allocate_dvec(dim);
    solution = allocate_dvec(dim);
    solution_ref = allocate_dvec(dim);
  }

  virtual void TearDown() {
    free_dvec(state);
    free_dvec(state_ref);
    free_dvec(solution);
    free_dvec(solution_ref);
  }

  int dim;
  double time_param;
  double *state, *state_ref, *solution, *solution_ref;
  struct d_inexact_newton_for_zero_horizon_ocp_mfgmres_args args;
};


TEST_F(d_inexact_newton_for_zero_horizon_ocp_mfgmres_args_test, ptr) {
  time_param = (double)rand()/RAND_MAX;
  for (int i=0; i<dim; ++i) {
    state[i] = (double)rand()/RAND_MAX;
    state_ref[i] = state[i];
  }
  for (int i=0; i<dim; ++i) {
    solution[i] = (double)rand()/RAND_MAX;
    solution_ref[i] = solution[i];
  }
  args.current_time = time_param;
  args.current_state_ptr = state;
  args.current_solution_ptr = solution;

  EXPECT_EQ(time_param, args.current_time);
  for (int i=0; i<dim; ++i) {
    EXPECT_EQ(state_ref[i], args.current_state_ptr[i]);
  }
  for (int i=0; i<dim; ++i) {
    EXPECT_EQ(solution_ref[i], args.current_solution_ptr[i]);
  }
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}