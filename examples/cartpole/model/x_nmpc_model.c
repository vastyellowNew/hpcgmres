int NMPC_MODEL_STRSIZE() {
  return sizeof(struct NMPC_MODEL);
}


void NMPC_MODEL_F(struct NMPC_MODEL *model, REAL t, REAL *x, REAL *u, REAL *f) {
  REAL x0 = sin(x[1]);
  REAL x1 = 1.0/(model->m_c + model->m_p*pow(x0, 2));
  REAL x2 = cos(x[1]);
  REAL x3 = model->l*pow(x[1], 2);
  REAL x4 = model->m_p*x0;
  f[0] = x[2];
  f[1] = x[3];
  f[2] = x1*(u[0] + x4*(model->g*x2 + x3));
  f[3] = x1*(-model->g*x0*(model->m_c + model->m_p) - u[0]*x2 - x2*x3*x4)/model->l;
}


void NMPC_MODEL_PHIX(struct NMPC_MODEL *model, REAL t, REAL *x, REAL *phix) {
  phix[0] = (1.0/2.0)*model->q_terminal[0]*(2*x[0] - 2*model->x_ref[0]);
  phix[1] = (1.0/2.0)*model->q_terminal[1]*(2*x[1] - 2*model->x_ref[1]);
  phix[2] = (1.0/2.0)*model->q_terminal[2]*(2*x[2] - 2*model->x_ref[2]);
  phix[3] = (1.0/2.0)*model->q_terminal[3]*(2*x[3] - 2*model->x_ref[3]);
}


void NMPC_MODEL_HX(struct NMPC_MODEL *model, REAL t, REAL *x, REAL *u, 
                   REAL *lmd, REAL *hx) {
  REAL x0 = 2*x[1];
  REAL x1 = sin(x[1]);
  REAL x2 = cos(x[1]);
  REAL x3 = model->g*x2;
  REAL x4 = model->l*pow(x[1], 2);
  REAL x5 = model->m_p*(x3 + x4);
  REAL x6 = model->m_p*pow(x1, 2);
  REAL x7 = model->m_c + x6;
  REAL x8 = model->m_p*x1;
  REAL x9 = x2*x8;
  REAL x10 = 2*x9/pow(x7, 2);
  REAL x11 = 1.0/x7;
  REAL x12 = model->l*x0;
  REAL x13 = model->g*x1;
  REAL x14 = model->m_c + model->m_p;
  REAL x15 = lmd[3]/model->l;
  hx[0] = (1.0/2.0)*model->q[0]*(2*x[0] - 2*model->x_ref[0]);
  hx[1] = -lmd[2]*x10*(u[0] + x1*x5) + lmd[2]*x11*(x2*x5 + x8*(x12 - x13)) + (1.0/2.0)*model->q[1]*(x0 - 2*model->x_ref[1]) - x10*x15*(-u[0]*x2 - x13*x14 - x4*x9) + x11*x15*(-model->m_p*pow(x2, 2)*x4 + u[0]*x1 - x12*x9 - x14*x3 + x4*x6);
  hx[2] = lmd[0] + (1.0/2.0)*model->q[2]*(2*x[2] - 2*model->x_ref[2]);
  hx[3] = lmd[1] + (1.0/2.0)*model->q[3]*(2*x[3] - 2*model->x_ref[3]);
}


void NMPC_MODEL_HU(struct NMPC_MODEL *model, REAL t, REAL *x, REAL *u, 
                   REAL *lmd, REAL *hu) {
  REAL x0 = 2*u[2];
  REAL x1 = 1.0/(model->m_c + model->m_p*pow(sin(x[1]), 2));
  hu[0] = lmd[2]*x1 + model->r[0]*u[0] + u[0]*x0 - lmd[3]*x1*cos(x[1])/model->l;
  hu[1] = -model->dummy_weight + u[1]*x0;
  hu[2] = pow(u[0], 2) + pow(u[1], 2) - 1.0/4.0*pow(model->u_max - model->u_min, 2);
}