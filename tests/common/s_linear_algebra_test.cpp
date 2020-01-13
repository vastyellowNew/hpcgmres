#include <stdlib.h>
#include <time.h>

#include <gtest/gtest.h>

extern "C" {
#include "s_memory_manager.h"
#include "s_linear_algebra.h"
}


class s_linear_algebra_test: public ::testing::Test {
protected:
  virtual void SetUp() {
    srand(time(NULL));
    dim = 50;
    x = allocate_svec(dim);
    x_ref = allocate_svec(dim);
    y = allocate_svec(dim);
    y_ref = allocate_svec(dim);
    res = allocate_svec(dim);
    res_ref = allocate_svec(dim);
    a = (float)rand()/RAND_MAX;
    b = (float)rand()/RAND_MAX;
  }

  virtual void TearDown() {
    free_svec(x);
    free_svec(x_ref);
    free_svec(y);
    free_svec(y_ref);
    free_svec(res);
    free_svec(res_ref);
  }

  int dim;
  float *x, *x_ref, *y, *y_ref, *res, *res_ref;
  float a, b;
};


TEST_F(s_linear_algebra_test, svecset) {
  for (int i=0; i<dim; ++i) {
    x_ref[i] = a;
  }
  cgmres_svecset(dim, a, x);
  for (int i=0; i<dim; ++i) {
    EXPECT_EQ(x_ref[i], x[i]);
  }
}


TEST_F(s_linear_algebra_test, sveccp) {
  for (int i=0; i<dim; ++i) {
    x[i] = (float)rand()/RAND_MAX;
  }
  for (int i=0; i<dim; ++i) {
    y_ref[i] = x[i];
  }
  cgmres_sveccp(dim, x, y);
  for (int i=0; i<dim; ++i) {
    EXPECT_EQ(y_ref[i], y[i]);
  }
}


TEST_F(s_linear_algebra_test, svecmcp) {
  for (int i=0; i<dim; ++i) {
    x_ref[i] = (float)rand()/RAND_MAX;
  }
  for (int i=0; i<dim; ++i) {
    y_ref[i] = b*x[i];
  }
  cgmres_svecmcp(dim, b, x, y);
  for (int i=0; i<dim; ++i) {
    EXPECT_EQ(y_ref[i], y[i]);
  }
}


TEST_F(s_linear_algebra_test, svecdot) {
  for (int i=0; i<dim; ++i) {
    x[i] = (float)rand()/RAND_MAX;
  }
  for (int i=0; i<dim; ++i) {
    y[i] = (float)rand()/RAND_MAX;
  }
  float dot_ref = 0;
  for (int i=0; i<dim; ++i) {
    dot_ref += x[i] * y[i];
  }
  float dot = cgmres_svecdot(dim, x, y);
  EXPECT_EQ(dot_ref, dot);
}


TEST_F(s_linear_algebra_test, svecnrm2) {
  for (int i=0; i<dim; ++i) {
    x[i] = (float)rand()/RAND_MAX;
  }
  float nrm_ref = 0;
  for (int i=0; i<dim; ++i) {
    nrm_ref += x[i] * x[i];
  }
  float nrm = cgmres_svecnrm2(dim, x);
  EXPECT_EQ(nrm_ref, nrm);
}


TEST_F(s_linear_algebra_test, svecmul) {
  for (int i=0; i<dim; ++i) {
    x[i] = (float)rand()/RAND_MAX;
  }
  for (int i=0; i<dim; ++i) {
    x_ref[i] = b * x[i];
  }
  cgmres_svecmul(dim, b, x);
  for (int i=0; i<dim; ++i) {
    EXPECT_EQ(x_ref[i], x[i]);
  }
}


TEST_F(s_linear_algebra_test, svecadd) {
  for (int i=0; i<dim; ++i) {
    x[i] = (float)rand()/RAND_MAX;
  }
  for (int i=0; i<dim; ++i) {
    y[i] = (float)rand()/RAND_MAX;
  }
  for (int i=0; i<dim; ++i) {
    y_ref[i] = x[i] + y[i];
  }
  cgmres_svecadd(dim, x, y);
  for (int i=0; i<dim; ++i) {
    EXPECT_EQ(y_ref[i], y[i]);
  }
}


TEST_F(s_linear_algebra_test, svecmadd) {
  for (int i=0; i<dim; ++i) {
    x[i] = (float)rand()/RAND_MAX;
  }
  for (int i=0; i<dim; ++i) {
    y[i] = (float)rand()/RAND_MAX;
  }
  for (int i=0; i<dim; ++i) {
    y_ref[i] = a*x[i] + y[i];
  }
  cgmres_svecmadd(dim, a, x, y);
  for (int i=0; i<dim; ++i) {
    EXPECT_EQ(y_ref[i], y[i]);
  }
}


TEST_F(s_linear_algebra_test, saxpy) {
  for (int i=0; i<dim; ++i) {
    x[i] = (float)rand()/RAND_MAX;
  }
  for (int i=0; i<dim; ++i) {
    y[i] = (float)rand()/RAND_MAX;
  }
  for (int i=0; i<dim; ++i) {
    res_ref[i] = a*x[i] + y[i];
  }
  cgmres_saxpy(dim, a, x, y, res);
  for (int i=0; i<dim; ++i) {
    EXPECT_EQ(res_ref[i], res[i]);
  }
}


TEST_F(s_linear_algebra_test, saxpby) {
  for (int i=0; i<dim; ++i) {
    x[i] = (float)rand()/RAND_MAX;
  }
  for (int i=0; i<dim; ++i) {
    y[i] = (float)rand()/RAND_MAX;
  }
  for (int i=0; i<dim; ++i) {
    res_ref[i] = a*x[i] + b*y[i];
  }
  cgmres_saxpby(dim, a, x, b, y, res);
  for (int i=0; i<dim; ++i) {
    EXPECT_EQ(res_ref[i], res[i]);
  }
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}