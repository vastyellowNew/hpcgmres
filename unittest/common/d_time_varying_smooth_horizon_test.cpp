#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <gtest/gtest.h>

extern "C" {
#include "d_time_varying_smooth_horizon.h"
}


class d_time_varying_smooth_horizon_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    T_f = (double)rand()/RAND_MAX;
    alpha = (double)rand()/RAND_MAX;
    initial_time = (double)rand()/RAND_MAX;
    d_time_varying_smooth_horizon_create(&horizon, T_f, alpha, initial_time);
  }

  virtual void TearDown() {
  }

  struct d_time_varying_smooth_horizon horizon;
  double T_f, alpha, initial_time;
};


TEST_F(d_time_varying_smooth_horizon_test, strsize) {
  EXPECT_EQ(d_time_varying_smooth_horizon_strsize(), 3*sizeof(double));
}


TEST_F(d_time_varying_smooth_horizon_test, get_length) {
  double current_time = (double)rand()/RAND_MAX + initial_time;
  double expect_length = T_f * (1.0 - exp(-alpha*(current_time - initial_time)));
  EXPECT_EQ(d_time_varying_smooth_horizon_get_length(&horizon, current_time), 
            expect_length);
}


TEST_F(d_time_varying_smooth_horizon_test, reset_length) {
  double current_time = (double)rand()/RAND_MAX + initial_time;
  d_time_varying_smooth_horizon_reset_length(&horizon, current_time);
  EXPECT_EQ(d_time_varying_smooth_horizon_get_length(&horizon, current_time), 0);
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}