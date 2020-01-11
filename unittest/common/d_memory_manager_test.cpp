#include <gtest/gtest.h>

extern "C" {
#include "d_memory_manager.h"
}


class d_memory_manager : public ::testing::Test {
protected:
  virtual void SetUp() {
    dim = 50;
  }

  virtual void TearDown() {
  }

  int dim;
};


TEST_F(d_memory_manager, allocate_vec) {
  double *vec = allocate_dvec(dim);
  for (int i=0; i<dim; ++i) {
    EXPECT_EQ(vec[i], 0);
  }
  free_dvec(vec);
}


TEST_F(d_memory_manager, allocate_mat) {
  double **mat = allocate_dmat(dim, dim);
  for (int i=0; i<dim; ++i) {
    for (int j=0; j<dim; ++j) {
      EXPECT_EQ(mat[i][j], 0);
    }
  }
  free_dmat(mat);
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}