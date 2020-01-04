#ifndef HPCGMRES_S_NMPC_MODEL_H_
#define HPCGMRES_S_NMPC_MODEL_H_


#include <cmath.h>
#include <cstdlib.h>


struct s_nmpc_model {
  int dimx = 4;
  int dimu = 2;
  int dimc = 1;

  float m_c = 2;
  float m_p = 0.2;
  float l = 0.5;
  float g = 9.80665;
  float u_min = -15;
  float u_max = 15;
  float dummy_weight = 0.1;

  float q[4] = {2.5, 10, 0.01, 0.01};
  float r[2] = {1, 0.01};
  float q_terminal[4] = {2.5, 10, 0.01, 0.01};
  float x_ref[4] = {0, M_PI, 0, 0};
}

int s_nmpc_model_strsize();

void s_nmpc_model_f(struct s_nmpc_model *model, float t, float *x, float *u, 
                    float *f);

void s_nmpc_model_phix(struct s_nmpc_model *model, float t, float *x, 
                       float *phix);

void s_nmpc_model_hx(struct s_nmpc_model *model, float t, float *x, float *u, 
                     float *lmd, float *hx);

void s_nmpc_model_hu(struct s_nmpc_model *model, float t, float *x, float *u, 
                     float *lmd, float *hu);


#endif // HPCGMRES_S_NMPC_MODEL_H_