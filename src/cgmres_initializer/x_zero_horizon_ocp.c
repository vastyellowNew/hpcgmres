int ZERO_HORIZON_OCP_STRSIZE() {
  return sizeof(struct ZERO_HORIZON_OCP);
}


int ZERO_HORIZON_OCP_MEMSIZE() {
  int size = 0;
  size += NMPC_MODEL_DIMX()*sizeof(REAL); // lmd_vec
  size = (size+63)/64*64; // make multiple of typical cache line size
  size += 64; // align to typical cache line size
  return size;
}


void ZERO_HORIZON_OCP_CREATE(struct ZERO_HORIZON_OCP *ocp) {
  NMPC_MODEL_CREATE(&ocp->model);
  ocp->lmd_vec = ALLOCATE_VEC(NMPC_MODEL_DIMX());
  ocp->dimx = ocp->model.dimx;
  ocp->dimu = ocp->model.dimu;
  ocp->dimc = ocp->model.dimc;
  ocp->dimuc = ocp->model.dimu + ocp->model.dimc;
  ocp->dim_solution = ocp->model.dimu + ocp->model.dimc;
  ocp->memsize = ZERO_HORIZON_OCP_MEMSIZE();
}


void ZERO_HORIZON_OCP_DELETE(struct ZERO_HORIZON_OCP *ocp) {
  FREE_VEC(ocp->lmd_vec);
}


void ZERO_HORIZON_OCP_COMPUTE_OPTIMALITY_RESIDUAL(
    struct ZERO_HORIZON_OCP *ocp, REAL current_time, REAL *current_state, 
    REAL *solution, REAL *optimality_residual) {
  NMPC_MODEL_PHIX(&ocp->model, current_time, current_state, ocp->lmd_vec);
  NMPC_MODEL_HU(&ocp->model, current_time, current_state, solution, 
                ocp->lmd_vec, optimality_residual);
}


void ZERO_HORIZON_OCP_COMPUTE_TERMINAL_COST_DERIVATIVE(
    struct ZERO_HORIZON_OCP *ocp, REAL current_time, REAL *current_state, 
    REAL *terminal_cost_derivative) {
  NMPC_MODEL_PHIX(&ocp->model, current_time, current_state, 
                  terminal_cost_derivative);
}


int ZERO_HORIZON_OCP_DIMX() {
  return NMPC_MODEL_DIMX();
}

int ZERO_HORIZON_OCP_DIMU() {
  return NMPC_MODEL_DIMU();
}

int ZERO_HORIZON_OCP_DIMC() {
  return NMPC_MODEL_DIMC();
}