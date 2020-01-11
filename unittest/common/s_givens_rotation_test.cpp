#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <blasfeo.h>

#include <gtest/gtest.h>

extern "C" {
#include "s_givens_rotation.h"
}


class s_givens_rotation_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    dim = 50;
    blasfeo_allocate_smat(dim, dim, &mat);
    blasfeo_allocate_svec(dim, &vec);
    blasfeo_allocate_svec(dim, &givens_c_vec);
    blasfeo_allocate_svec(dim, &givens_s_vec);
  }

  virtual void TearDown() {
    blasfeo_free_smat(&mat);
    blasfeo_free_svec(&vec);
    blasfeo_free_svec(&givens_c_vec);
    blasfeo_free_svec(&givens_s_vec);
  }

  int dim;
  struct blasfeo_smat mat;
  struct blasfeo_svec vec, givens_c_vec, givens_s_vec;
};


TEST_F(s_givens_rotation_test, blasfeo_mat) {
  for (int i=0; i<dim; ++i) {
    for (int j=0; j<dim; ++j) {
      blasfeo_sgein1((float)rand()/RAND_MAX, &mat, i, j);
    }
  }
  blasfeo_svecse(dim, 0.0, &givens_c_vec, 0);
  blasfeo_svecse(dim, 0.0, &givens_s_vec, 0);

  for (int i=0; i<dim-1; ++i) {
    for (int j=0; j<dim-1; ++j) {
      float tmp1 = blasfeo_svecex1(&givens_c_vec, j) * blasfeo_sgeex1(&mat, i, j)
                    - blasfeo_svecex1(&givens_s_vec, j) * blasfeo_sgeex1(&mat, i, j+1);
      float tmp2 = blasfeo_svecex1(&givens_s_vec, j) * blasfeo_sgeex1(&mat, i, j) 
                    + blasfeo_svecex1(&givens_c_vec, j) * blasfeo_sgeex1(&mat, i, j+1);
      s_apply_givens_rotation_to_mat(&mat, &givens_c_vec, &givens_s_vec, i, j);
      EXPECT_EQ(tmp1, blasfeo_sgeex1(&mat, i, j));
      EXPECT_EQ(tmp2, blasfeo_sgeex1(&mat, i, j+1));
    }
  }
}


TEST_F(s_givens_rotation_test, blasfeo_vec) {
  for (int i=0; i<dim; ++i) {
    blasfeo_svecin1((float)rand()/RAND_MAX, &vec, i);
  }
  blasfeo_svecse(dim, 0.0, &givens_c_vec, 0);
  blasfeo_svecse(dim, 0.0, &givens_s_vec, 0);

  for (int i=0; i<dim-1; ++i) {
    float tmp1 = blasfeo_svecex1(&givens_c_vec, i) * blasfeo_svecex1(&vec, i)
                  - blasfeo_svecex1(&givens_s_vec, i) * blasfeo_svecex1(&vec, i+1);
    float tmp2 = blasfeo_svecex1(&givens_s_vec, i) * blasfeo_svecex1(&vec, i) 
                  + blasfeo_svecex1(&givens_c_vec, i) * blasfeo_svecex1(&vec, i+1);
    s_apply_givens_rotation_to_vec(&vec, &givens_c_vec, &givens_s_vec, i);
    EXPECT_EQ(tmp1, blasfeo_svecex1(&vec, i));
    EXPECT_EQ(tmp2, blasfeo_svecex1(&vec, i+1));
  }
}


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