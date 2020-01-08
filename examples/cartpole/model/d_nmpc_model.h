#ifndef HPCGMRES_D_NMPC_MODEL_H_
#define HPCGMRES_D_NMPC_MODEL_H_


#ifdef __cplusplus
extern "C" {
#endif

struct d_nmpc_model {
  int dimx, dimu, dimc;
  double m_c, m_p, l, g, u_min, u_max, dummy_weight;
  double q[4], r[2], q_terminal[4], x_ref[4];
};

int d_nmpc_model_strsize();

void d_nmpc_model_create(struct d_nmpc_model *model);

void d_nmpc_model_f(struct d_nmpc_model *model, double t, double *x, double *u, 
                    double *f);

void d_nmpc_model_phix(struct d_nmpc_model *model, double t, double *x, 
                       double *phix);

void d_nmpc_model_hx(struct d_nmpc_model *model, double t, double *x, 
                     double *u, double *lmd, double *hx);

void d_nmpc_model_hu(struct d_nmpc_model *model, double t, double *x, 
                     double *u, double *lmd, double *hu);

int d_nmpc_model_dimx();
int d_nmpc_model_dimu();
int d_nmpc_model_dimc();

#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // HPCGMRES_D_NMPC_MODEL_H_