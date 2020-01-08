#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <gtest/gtest.h>

extern "C" {
#include "d_nmpc_model.h"
}


class d_nmpc_model_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    d_nmpc_model_create(&model);
  }

  struct d_nmpc_model model;
};


TEST_F(d_nmpc_model_test, dim) {
  EXPECT_EQ(4, d_nmpc_model_dimx());
  EXPECT_EQ(2, d_nmpc_model_dimu());
  EXPECT_EQ(1, d_nmpc_model_dimc());
  EXPECT_EQ(4, model.dimx);
  EXPECT_EQ(2, model.dimu);
  EXPECT_EQ(1, model.dimc);
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}