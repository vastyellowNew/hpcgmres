#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <gtest/gtest.h>

extern "C" {
#include "s_nmpc_model.h"
}


class s_nmpc_model_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    s_nmpc_model_create(&model);
    m_c = 2;
    m_p = 0.2;
    l = 0.5;
    g = 9.80665;
    u_min = -15;
    u_max = 15;
    dummy_weight = 0.1;

    q[0] = 2.5; 
    q[1] = 10; 
    q[2] = 0.01; 
    q[3] = 0.01; 

    r[0] = 1.0; 
    r[1] = 0.01; 

    q_terminal[0] = 2.5; 
    q_terminal[1] = 10; 
    q_terminal[2] = 0.01; 
    q_terminal[3] = 0.01; 

    x_ref[0] = 0; 
    x_ref[1] = M_PI; 
    x_ref[2] = 0; 
    x_ref[3] = 0; 

    t = (float)rand()/RAND_MAX;
    for (int i=0; i<4; ++i) {
      x[i] = (float)rand()/RAND_MAX;
    }
    for (int i=0; i<3; ++i) {
      u[i] = (float)rand()/RAND_MAX;
    }
    for (int i=0; i<4; ++i) {
      lmd[i] = (float)rand()/RAND_MAX;
    }
  }

  struct s_nmpc_model model;
  float t, x[4], u[3], lmd[4];
  float m_c, m_p, l, g, u_min, u_max, dummy_weight;
  float q[4], r[2], q_terminal[4], x_ref[4];
};


TEST_F(s_nmpc_model_test, dim) {
  EXPECT_EQ(4, s_nmpc_model_dimx());
  EXPECT_EQ(2, s_nmpc_model_dimu());
  EXPECT_EQ(1, s_nmpc_model_dimc());
  EXPECT_EQ(4, model.dimx);
  EXPECT_EQ(2, model.dimu);
  EXPECT_EQ(1, model.dimc);
}


TEST_F(s_nmpc_model_test, f) {
  float f_ref[4], f[4];
  float x0 = sin(x[1]);
  float x1 = 1.0/(m_c + m_p*pow(x0, 2));
  float x2 = cos(x[1]);
  float x3 = l*pow(x[1], 2);
  float x4 = m_p*x0;
  f_ref[0] = x[2];
  f_ref[1] = x[3];
  f_ref[2] = x1*(u[0] + x4*(g*x2 + x3));
  f_ref[3] = x1*(-g*x0*(m_c + m_p) - u[0]*x2 - x2*x3*x4)/l;
  s_nmpc_model_f(&model, t, x, u, f);
  for (int i=0; i<model.dimx; ++i) {
    EXPECT_EQ(f_ref[i], f[i]);
  }
}


TEST_F(s_nmpc_model_test, phix) {
  float phix_ref[4], phix[4];
  phix_ref[0] = (1.0/2.0)*q_terminal[0]*(2*x[0] - 2*x_ref[0]);
  phix_ref[1] = (1.0/2.0)*q_terminal[1]*(2*x[1] - 2*x_ref[1]);
  phix_ref[2] = (1.0/2.0)*q_terminal[2]*(2*x[2] - 2*x_ref[2]);
  phix_ref[3] = (1.0/2.0)*q_terminal[3]*(2*x[3] - 2*x_ref[3]);
  s_nmpc_model_phix(&model, t, x, phix);
  for (int i=0; i<model.dimx; ++i) {
    EXPECT_EQ(phix_ref[i], phix[i]);
  }
}


TEST_F(s_nmpc_model_test, hx) {
  float hx_ref[4], hx[4];
  float x0 = 2*x[1];
  float x1 = sin(x[1]);
  float x2 = cos(x[1]);
  float x3 = g*x2;
  float x4 = l*pow(x[1], 2);
  float x5 = m_p*(x3 + x4);
  float x6 = m_p*pow(x1, 2);
  float x7 = m_c + x6;
  float x8 = m_p*x1;
  float x9 = x2*x8;
  float x10 = 2*x9/pow(x7, 2);
  float x11 = 1.0/x7;
  float x12 = l*x0;
  float x13 = g*x1;
  float x14 = m_c + m_p;
  float x15 = lmd[3]/l;
  hx_ref[0] = (1.0/2.0)*q[0]*(2*x[0] - 2*x_ref[0]);
  hx_ref[1] = -lmd[2]*x10*(u[0] + x1*x5) + lmd[2]*x11*(x2*x5 + x8*(x12 - x13)) + (1.0/2.0)*q[1]*(x0 - 2*x_ref[1]) - x10*x15*(-u[0]*x2 - x13*x14 - x4*x9) + x11*x15*(-m_p*pow(x2, 2)*x4 + u[0]*x1 - x12*x9 - x14*x3 + x4*x6);
  hx_ref[2] = lmd[0] + (1.0/2.0)*q[2]*(2*x[2] - 2*x_ref[2]);
  hx_ref[3] = lmd[1] + (1.0/2.0)*q[3]*(2*x[3] - 2*x_ref[3]);
  s_nmpc_model_hx(&model, t, x, u, lmd, hx);
  for (int i=0; i<model.dimx; ++i) {
    EXPECT_EQ(hx_ref[i], hx[i]);
  }
}


TEST_F(s_nmpc_model_test, hu) {
  float hu_ref[3], hu[3];
  float x0 = 2*u[2];
  float x1 = 1.0/(m_c + m_p*pow(sin(x[1]), 2));
  hu_ref[0] = lmd[2]*x1 + r[0]*u[0] + u[0]*x0 - lmd[3]*x1*cos(x[1])/l;
  hu_ref[1] = -dummy_weight + u[1]*x0;
  hu_ref[2] = pow(u[0], 2) + pow(u[1], 2) - 1.0/4.0*pow(u_max - u_min, 2);
  s_nmpc_model_hu(&model, t, x, u, lmd, hu);
  for (int i=0; i<model.dimu+model.dimc; ++i) {
    EXPECT_EQ(hu_ref[i], hu[i]);
  }
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}