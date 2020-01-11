#include <gtest/gtest.h>

extern "C" {
#include "s_memory_manager.h"
}


class s_memory_manager : public ::testing::Test {
protected:
  virtual void SetUp() {
    dim = 50;
  }

  virtual void TearDown() {
  }

  int dim;
};


TEST_F(s_memory_manager, allocate_vec) {
  float *vec = allocate_svec(dim);
  for (int i=0; i<dim; ++i) {
    EXPECT_EQ(vec[i], 0);
  }
  free_svec(vec);
}


TEST_F(s_memory_manager, allocate_mat) {
  float **mat = allocate_smat(dim, dim);
  for (int i=0; i<dim; ++i) {
    for (int j=0; j<dim; ++j) {
      EXPECT_EQ(mat[i][j], 0);
    }
  }
  free_smat(mat);
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}