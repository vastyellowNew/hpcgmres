#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <blasfeo.h>

#include <gtest/gtest.h>

extern "C" {
#include "d_givens_rotation.h"
}


class d_givens_rotation_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    dim = 50;
    blasfeo_allocate_dmat(dim, dim, &mat);
    blasfeo_allocate_dvec(dim, &vec);
    blasfeo_allocate_dvec(dim, &givens_c_vec);
    blasfeo_allocate_dvec(dim, &givens_s_vec);
  }

  virtual void TearDown() {
    blasfeo_free_dmat(&mat);
    blasfeo_free_dvec(&vec);
    blasfeo_free_dvec(&givens_c_vec);
    blasfeo_free_dvec(&givens_s_vec);
  }

  int dim;
  struct blasfeo_dmat mat;
  struct blasfeo_dvec vec, givens_c_vec, givens_s_vec;
};


TEST_F(d_givens_rotation_test, test_mat) {
  for (int i=0; i<dim; ++i) {
    for (int j=0; j<dim; ++j) {
      blasfeo_dgein1((double)rand()/RAND_MAX, &mat, i, j);
    }
  }
  blasfeo_dvecse(dim, 0.0, &givens_c_vec, 0);
  blasfeo_dvecse(dim, 0.0, &givens_s_vec, 0);

  for (int i=0; i<dim-1; ++i) {
    for (int j=0; j<dim-1; ++j) {
      double tmp1 = blasfeo_dvecex1(&givens_c_vec, j) * blasfeo_dgeex1(&mat, i, j)
                    - blasfeo_dvecex1(&givens_s_vec, j) * blasfeo_dgeex1(&mat, i, j+1);
      double tmp2 = blasfeo_dvecex1(&givens_s_vec, j) * blasfeo_dgeex1(&mat, i, j) 
                    + blasfeo_dvecex1(&givens_c_vec, j) * blasfeo_dgeex1(&mat, i, j+1);
      d_apply_givens_rotation_to_mat(&mat, &givens_c_vec, &givens_s_vec, i, j);
      EXPECT_EQ(tmp1, blasfeo_dgeex1(&mat, i, j));
      EXPECT_EQ(tmp2, blasfeo_dgeex1(&mat, i, j+1));
    }
  }
}


TEST_F(d_givens_rotation_test, test_vec) {
  for (int i=0; i<dim; ++i) {
    blasfeo_dvecin1((double)rand()/RAND_MAX, &vec, i);
  }
  blasfeo_dvecse(dim, 0.0, &givens_c_vec, 0);
  blasfeo_dvecse(dim, 0.0, &givens_s_vec, 0);

  for (int i=0; i<dim-1; ++i) {
    double tmp1 = blasfeo_dvecex1(&givens_c_vec, i) * blasfeo_dvecex1(&vec, i)
                  - blasfeo_dvecex1(&givens_s_vec, i) * blasfeo_dvecex1(&vec, i+1);
    double tmp2 = blasfeo_dvecex1(&givens_s_vec, i) * blasfeo_dvecex1(&vec, i) 
                  + blasfeo_dvecex1(&givens_c_vec, i) * blasfeo_dvecex1(&vec, i+1);
    d_apply_givens_rotation_to_vec(&vec, &givens_c_vec, &givens_s_vec, i);
    EXPECT_EQ(tmp1, blasfeo_dvecex1(&vec, i));
    EXPECT_EQ(tmp2, blasfeo_dvecex1(&vec, i+1));
  }
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}