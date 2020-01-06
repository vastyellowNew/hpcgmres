#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <gtest/gtest.h>

extern "C" {
#include "d_givens_rotation.h"
}


class d_givens_rotation_test : public ::testing::Test {
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
  double vec[dim];
  double givens_c_vec[50];
  double givens_s_vec[50];
};


TEST_F(d_givens_rotation_test, mat) {
  for (int i=0; i<dim; ++i) {
    vec[i] = (double)rand()/RAND_MAX;
  }
  for (int i=0; i<dim-1; ++i) {
    double tmp1 = givens_c_vec[i] * vec[i] - givens_s_vec[i] * vec[i+1];
    double tmp2 = givens_s_vec[i] * vec[i] - givens_c_vec[i] * vec[i+1];
    d_apply_givens_rotation(vec, givens_c_vec, givens_s_vec, i);
    EXPECT_EQ(tmp1, vec[i]);
    EXPECT_EQ(tmp2, vec[i+1]);
  }
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}