#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>

#include <gtest/gtest.h>

extern "C" {
#include "s_time_varying_smooth_horizon.h"
}


class s_time_varying_smooth_horizon_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    T_f = (float)rand()/RAND_MAX;
    alpha = (float)rand()/RAND_MAX;
    initial_time = (float)rand()/RAND_MAX;
    s_time_varying_smooth_horizon_create(&horizon, T_f, alpha, initial_time);
  }

  virtual void TearDown() {
  }

  struct s_time_varying_smooth_horizon horizon;
  float T_f, alpha, initial_time;
};


TEST_F(s_time_varying_smooth_horizon_test, get_length) {
  float current_time = (float)rand()/RAND_MAX + initial_time;
  float expect_length = T_f * (1.0 - exp(-alpha*(current_time - initial_time)));
  float diff = s_time_varying_smooth_horizon_get_length(&horizon, current_time)
               - expect_length;
  EXPECT_TRUE(fabs(diff) < FLT_EPSILON);
}


TEST_F(s_time_varying_smooth_horizon_test, reset_length) {
  float current_time = (float)rand()/RAND_MAX + initial_time;
  s_time_varying_smooth_horizon_reset_length(&horizon, current_time);
  EXPECT_EQ(s_time_varying_smooth_horizon_get_length(&horizon, current_time), 0);
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}