#ifndef HPCGMRES_D_NMPC_MODEL_H_
#define HPCGMRES_D_NMPC_MODEL_H_


#include <cmath.h>
#include <cstdlib.h>


struct d_nmpc_model {
  int dimx = 4;
  int dimu = 2;
  int dimc = 1;

  double m_c = 2;
  double m_p = 0.2;
  double l = 0.5;
  double g = 9.80665;
  double u_min = -15;
  double u_max = 15;
  double dummy_weight = 0.1;

  double q[4] = {2.5, 10, 0.01, 0.01};
  double r[2] = {1, 0.01};
  double q_terminal[4] = {2.5, 10, 0.01, 0.01};
  double x_ref[4] = {0, M_PI, 0, 0};
}

int d_nmpc_model_strsize();

void d_nmpc_model_f(struct d_nmpc_model *model, double t, double *x, double *u, 
                    double *f);

void d_nmpc_model_phix(struct d_nmpc_model *model, double t, double *x, 
                       double *phix);

void d_nmpc_model_hx(struct d_nmpc_model *model, double t, double *x, 
                     double *u, double *lmd, double *hx);

void d_nmpc_model_hu(struct d_nmpc_model *model, double t, double *x, 
                     double *u, double *lmd, double *hu);


#endif // HPCGMRES_D_NMPC_MODEL_H_