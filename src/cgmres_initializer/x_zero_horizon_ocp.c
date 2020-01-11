int ZERO_HORIZON_OCP_STRSIZE() {
  return sizeof(struct ZERO_HORIZON_OCP);
}


int ZERO_HORIZON_OCP_MEMSIZE() {
  int size = 0;
  size += 1*SIZE_STRVEC(NMPC_MODEL_DIMX()); // lmd_vec
  size += 1*sizeof(struct STRVEC*); // lmd_vec
  size = (size+63)/64*64; // make multiple of typical cache line size
  size += 64; // align to typical cache line size
  return size;
}


void ZERO_HORIZON_OCP_CREATE(struct ZERO_HORIZON_OCP *ocp, void *memory) {
  // zero memory (to avoid corrupted memory like e.g. NaN)
  int memsize = ZERO_HORIZON_OCP_MEMSIZE();
  int memsize_m8 = memsize/8; // sizeof(double) is 8
  //	int memsize_r8 = memsize - 8*memsize_m8;
  double *double_ptr = memory;
  int ii;
  for(ii=0; ii<memsize_m8-7; ii+=8) {
    double_ptr[ii+0] = 0.0;
    double_ptr[ii+1] = 0.0;
    double_ptr[ii+2] = 0.0;
    double_ptr[ii+3] = 0.0;
    double_ptr[ii+4] = 0.0;
    double_ptr[ii+5] = 0.0;
    double_ptr[ii+6] = 0.0;
    double_ptr[ii+7] = 0.0;
    }
  // XXX exploit that it is multiple of 64 bytes !!!!!
  // for(; ii<memsize_m8; ii++) {
  // double_ptr[ii] = 0.0;
  // }
  // char *char_ptr = (char *) (&double_ptr[ii]);
  // for(ii=0; ii<memsize_r8; ii++) {
  // char_ptr[ii] = 0;
  // }

  // vector struct stuff
  struct STRVEC *sv_ptr = (struct STRVEC *) memory;
  ocp->lmd_vec = sv_ptr;
  sv_ptr += 1;

  // align to typical cache line size
  size_t s_ptr = (size_t) sv_ptr;
  s_ptr = (s_ptr+63)/64*64;

  // void stuff
  char *c_ptr = (char *) s_ptr;

  NMPC_MODEL_CREATE(&ocp->model);
  CREATE_STRVEC(ocp->model.dimx, ocp->lmd_vec, c_ptr);
  c_ptr += ocp->lmd_vec->memsize;
  VECSE(ocp->model.dimx, 0.0, ocp->lmd_vec, 0);

  ocp->dimx = ocp->model.dimx;
  ocp->dimu = ocp->model.dimu;
  ocp->dimc = ocp->model.dimc;
  ocp->dimuc = ocp->model.dimu + ocp->model.dimc;
  ocp->dim_solution = ocp->model.dimu + ocp->model.dimc;
  ocp->memsize = ZERO_HORIZON_OCP_MEMSIZE();

#if defined(RUNTIME_CHECKS)
  if(c_ptr > ((char *) mem) + ocp->memsize) {
    printf("\nerror: ZERO_HORIZON_OCP_CREATE: outside memory bounds!\n\n");
    exit(1);
  }
#endif
}


void ZERO_HORIZON_OCP_COMPUTE_OPTIMALITY_RESIDUAL(
    struct ZERO_HORIZON_OCP *ocp, REAL current_time, 
    struct STRVEC *current_state, struct STRVEC *solution, 
    struct STRVEC *optimality_residual) {
  NMPC_MODEL_PHIX(&ocp->model, current_time, current_state->pa, 
                  ocp->lmd_vec->pa);
  NMPC_MODEL_HU(&ocp->model, current_time, current_state->pa, solution->pa, 
                ocp->lmd_vec->pa, optimality_residual->pa);
}


void ZERO_HORIZON_OCP_COMPUTE_TERMINAL_COST_DERIVATIVE(
    struct ZERO_HORIZON_OCP *ocp, REAL current_time, 
    struct STRVEC *current_state, struct STRVEC *terminal_cost_derivative) {
  NMPC_MODEL_PHIX(&ocp->model, current_time, current_state->pa, 
                  terminal_cost_derivative->pa);
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
