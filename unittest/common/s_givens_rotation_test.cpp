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
    srand(time(NULL));
    dim = 50;
  }

  virtual void TearDown() {
  }

  int dim;
};


TEST_F(s_givens_rotation_test, test_vec) {
  float svec[dim], sgivens_c_vec[dim], sgivens_s_vec[dim];
  for (int i=0; i<dim; ++i) {
    svec[i] = (float)rand()/RAND_MAX;
    sgivens_c_vec[i] = 0.0;
    sgivens_s_vec[i] = 0.0;
  }
  for (int i=0; i<dim-1; ++i) {
    float tmp1 = sgivens_c_vec[i] * svec[i]
                    - sgivens_s_vec[i] * svec[i+1];
    float tmp2 = sgivens_s_vec[i] * svec[i]
                    + sgivens_c_vec[i] * svec[i+1];
    svec[i]   = tmp1;
    svec[i+1] = tmp2;
    s_apply_givens_rotation(svec, sgivens_c_vec, sgivens_s_vec, i);
    EXPECT_EQ(tmp1, svec[i]);
    EXPECT_EQ(tmp2, svec[i+1]);
  }
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}