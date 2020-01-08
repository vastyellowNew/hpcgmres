#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <gtest/gtest.h>

extern "C" {
#include "s_nmpc_model.h"
}


class s_nmpc_model_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    s_nmpc_model_create(&model);
  }

  struct s_nmpc_model model;
};


TEST_F(s_nmpc_model_test, dim) {
  EXPECT_EQ(4, s_nmpc_model_dimx());
  EXPECT_EQ(2, s_nmpc_model_dimu());
  EXPECT_EQ(1, s_nmpc_model_dimc());
  EXPECT_EQ(4, model.dimx);
  EXPECT_EQ(2, model.dimu);
  EXPECT_EQ(1, model.dimc);
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}