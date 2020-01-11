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
    srand(time(NULL));
    dim = 50;
  }

  virtual void TearDown() {
  }

  int dim;
};


TEST_F(d_givens_rotation_test, test_vec) {
  double dvec[dim], dgivens_c_vec[dim], dgivens_s_vec[dim];
  for (int i=0; i<dim; ++i) {
    dvec[i] = (double)rand()/RAND_MAX;
    dgivens_c_vec[i] = 0.0;
    dgivens_s_vec[i] = 0.0;
  }
  for (int i=0; i<dim-1; ++i) {
    double tmp1 = dgivens_c_vec[i] * dvec[i]
                    - dgivens_s_vec[i] * dvec[i+1];
    double tmp2 = dgivens_s_vec[i] * dvec[i]
                    + dgivens_c_vec[i] * dvec[i+1];
    dvec[i]   = tmp1;
    dvec[i+1] = tmp2;
    d_apply_givens_rotation(dvec, dgivens_c_vec, dgivens_s_vec, i);
    EXPECT_EQ(tmp1, dvec[i]);
    EXPECT_EQ(tmp2, dvec[i+1]);
  }
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}