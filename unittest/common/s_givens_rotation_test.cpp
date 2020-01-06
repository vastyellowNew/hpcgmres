#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <gtest/gtest.h>

extern "C" {
#include "s_givens_rotation.h"
}


class s_givens_rotation_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    for (int i=0; i<dim; ++i) {
      vec[i] = 0.0;
      givens_c_vec[i] = 0.0;
      givens_s_vec[i] = 0.0;
    }
  }

  virtual void TearDown() {
  }

  static constexpr int dim = 50;
  float vec[dim];
  float givens_c_vec[50];
  float givens_s_vec[50];
};


TEST_F(s_givens_rotation_test, mat) {
  for (int i=0; i<dim; ++i) {
    vec[i] = (float)rand()/RAND_MAX;
  }
  for (int i=0; i<dim-1; ++i) {
    float tmp1 = givens_c_vec[i] * vec[i] - givens_s_vec[i] * vec[i+1];
    float tmp2 = givens_s_vec[i] * vec[i] - givens_c_vec[i] * vec[i+1];
    s_apply_givens_rotation(vec, givens_c_vec, givens_s_vec, i);
    EXPECT_EQ(tmp1, vec[i]);
    EXPECT_EQ(tmp2, vec[i+1]);
  }
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}