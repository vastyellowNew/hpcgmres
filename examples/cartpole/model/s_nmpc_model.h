#ifndef HPCGMRES_S_NMPC_MODEL_H_
#define HPCGMRES_S_NMPC_MODEL_H_


#ifdef __cplusplus
extern "C" {
#endif

struct s_nmpc_model {
  int dimx, dimu, dimc;
  float m_c, m_p, l, g, u_min, u_max, dummy_weight;
  float q[4], r[2], q_terminal[4], x_ref[4];
};

int s_nmpc_model_strsize();

void s_nmpc_model_create(struct s_nmpc_model *model);

void s_nmpc_model_f(struct s_nmpc_model *model, float t, float *x, float *u, 
                    float *f);

void s_nmpc_model_phix(struct s_nmpc_model *model, float t, float *x, 
                       float *phix);

void s_nmpc_model_hx(struct s_nmpc_model *model, float t, float *x, float *u, 
                     float *lmd, float *hx);

void s_nmpc_model_hu(struct s_nmpc_model *model, float t, float *x, float *u, 
                     float *lmd, float *hu);

int s_nmpc_model_dimx();
int s_nmpc_model_dimu();
int s_nmpc_model_dimc();

#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // HPCGMRES_S_NMPC_MODEL_H_