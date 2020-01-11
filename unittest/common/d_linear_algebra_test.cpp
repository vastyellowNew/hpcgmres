#include <stdlib.h>
#include <time.h>

#include <gtest/gtest.h>

extern "C" {
#include "d_memory_manager.h"
#include "d_linear_algebra.h"
}


class d_linear_algebra_test: public ::testing::Test {
protected:
  virtual void SetUp() {
    dim = 50;
    x = allocate_dvec(dim);
    x_ref = allocate_dvec(dim);
    y = allocate_dvec(dim);
    y_ref = allocate_dvec(dim);
    res = allocate_dvec(dim);
    res_ref = allocate_dvec(dim);
    a = (double)rand()/RAND_MAX;
    b = (double)rand()/RAND_MAX;
  }

  virtual void TearDown() {
    free_dvec(x);
    free_dvec(x_ref);
    free_dvec(y);
    free_dvec(y_ref);
    free_dvec(res);
    free_dvec(res_ref);
  }

  int dim;
  double *x, *x_ref, *y, *y_ref, *res, *res_ref;
  double a, b;
};


TEST_F(d_linear_algebra_test, dvecset) {
  for (int i=0; i<dim; ++i) {
    x_ref[i] = a;
  }
  hpcgmres_dvecset(dim, a, x);
  for (int i=0; i<dim; ++i) {
    EXPECT_EQ(x_ref[i], x[i]);
  }
}


TEST_F(d_linear_algebra_test, dveccp) {
  hpcgmres_dvecset(dim, a, x);
  for (int i=0; i<dim; ++i) {
    y_ref[i] = x[i];
  }
  hpcgmres_dveccp(dim, x, y);
  for (int i=0; i<dim; ++i) {
    EXPECT_EQ(y_ref[i], y[i]);
  }
}


TEST_F(d_linear_algebra_test, dvecmcp) {
  hpcgmres_dvecset(dim, a, x);
  for (int i=0; i<dim; ++i) {
    y_ref[i] = b*x[i];
  }
  hpcgmres_dvecmcp(dim, b, x, y);
  for (int i=0; i<dim; ++i) {
    EXPECT_EQ(y_ref[i], y[i]);
  }
}


TEST_F(d_linear_algebra_test, dvecdot) {
  hpcgmres_dvecset(dim, a, x);
  hpcgmres_dvecset(dim, b, y);
  double dot_ref = 0;
  for (int i=0; i<dim; ++i) {
    dot_ref += x[i] * y[i];
  }
  double dot = hpcgmres_dvecdot(dim, x, y);
  EXPECT_EQ(dot_ref, dot);
}


TEST_F(d_linear_algebra_test, dvecnrm2) {
  hpcgmres_dvecset(dim, a, x);
  double nrm_ref = 0;
  for (int i=0; i<dim; ++i) {
    nrm_ref += x[i] * x[i];
  }
  double nrm = hpcgmres_dvecnrm2(dim, x);
  EXPECT_EQ(nrm_ref, nrm);
}


TEST_F(d_linear_algebra_test, dvecmul) {
  hpcgmres_dvecset(dim, a, x);
  for (int i=0; i<dim; ++i) {
    x_ref[i] = b * x[i];
  }
  hpcgmres_dvecmul(dim, b, x);
  for (int i=0; i<dim; ++i) {
    EXPECT_EQ(x_ref[i], x[i]);
  }
}


TEST_F(d_linear_algebra_test, dvecadd) {
  hpcgmres_dvecset(dim, a, x);
  hpcgmres_dvecset(dim, b, y);
  for (int i=0; i<dim; ++i) {
    y_ref[i] = x[i] + y[i];
  }
  hpcgmres_dvecadd(dim, x, y);
  for (int i=0; i<dim; ++i) {
    EXPECT_EQ(y_ref[i], y[i]);
  }
}


TEST_F(d_linear_algebra_test, dvecmadd) {
  hpcgmres_dvecset(dim, a, x);
  hpcgmres_dvecset(dim, b, y);
  for (int i=0; i<dim; ++i) {
    y_ref[i] = a*x[i] + y[i];
  }
  hpcgmres_dvecmadd(dim, a, x, y);
  for (int i=0; i<dim; ++i) {
    EXPECT_EQ(y_ref[i], y[i]);
  }
}


TEST_F(d_linear_algebra_test, daxpy) {
  hpcgmres_dvecset(dim, a, x);
  hpcgmres_dvecset(dim, b, y);
  for (int i=0; i<dim; ++i) {
    res_ref[i] = a*x[i] + y[i];
  }
  hpcgmres_daxpy(dim, a, x, y, res);
  for (int i=0; i<dim; ++i) {
    EXPECT_EQ(res_ref[i], res[i]);
  }
}


TEST_F(d_linear_algebra_test, daxpby) {
  hpcgmres_dvecset(dim, a, x);
  hpcgmres_dvecset(dim, b, y);
  for (int i=0; i<dim; ++i) {
    res_ref[i] = a*x[i] + b*y[i];
  }
  hpcgmres_daxpby(dim, a, x, b, y, res);
  for (int i=0; i<dim; ++i) {
    EXPECT_EQ(res_ref[i], res[i]);
  }
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}