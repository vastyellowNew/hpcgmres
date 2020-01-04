int SINGLE_SHOOTING_OCP_STRSIZE() {
  return sizeof(struct SINGLE_SHOOTING_OCP);
}


int SINGLE_SHOOTING_OCP_MEMSIZE(struct SINGLE_SHOOTING_OCP *ocp, int N) {
  int size = 0;

  size += TIME_VARYING_SMOOTH_HORIZON_MEMSIZE()

  // size of pointers and int
  size += 1*sizeof(struct STRMAT); // hessenberg_mat
  size += 1*(kmax+1)*sizeof(struct STRVEC); // basis_vec
  size += 4*sizeof(struct STRVEC); // b_vec, givens_c_vec, givens_s_vec, g_vec

  // size of matrices and vectors
  size += 1*SIZE_STRVEC(ocp->model->dimx); // dx_vec
  size += 2*SIZE_STRVEC(N*ocp->model->dimx); // x_sequence_vec, lmd_sequence_vec

  size = (size+63)/64*64; // make multiple of typical cache line size
  size += 64; // align to typical cache line size

  return size;
}


void SINGLE_SHOOTING_OCP_CREATE(struct SINGLE_SHOOTING_OCP *ocp, REAL T_f, 
                                REAL alpha, REAL initial_time, int N, 
                                void *memory) {
  TIME_VARYING_SMOOTH_HORIZON_CREATE(ocp->horizon, T_f, alpha, initial_time);

  // zero memory (to avoid corrupted memory like e.g. NaN)
  int memsize = SINGLE_SHOOTING_OCP_MEMSIZE(ocp, N);
  int memsize_m8 = memsize/8; // sizeof(double) is 8
  //	int memsize_r8 = memsize - 8*memsize_m8;
  double *double_ptr = mem;
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
  mfgmres->dx_vec = sv_ptr;
  sv_ptr += 1;
  mfgmres->x_sequence_vec = sv_ptr;
  sv_ptr += N;
  mfgmres->lmd_sequence_vec = sv_ptr;
  sv_ptr += N;

  // align to typical cache line size
  size_t s_ptr = (size_t) sv_ptr;
  s_ptr = (s_ptr+63)/64*64;

  // void stuff
  char *c_ptr = (char *) s_ptr;

  int dimx = ocp->model->dimx;

  CREATE_STRVEC(dimx, ocp->dx_vec, c_ptr);
  c_ptr += ocp->dx_vec->memsize;
  VECCSE(dimx, 0.0, ocp->dx_vec, 0);

  char *tmp_ptr = c_ptr;
  c_ptr += SIZE_STRVEC(N*dimx);
  for (ii=0; ii<N; ++ii) {
    CREATE_STRVEC(dimx, ocp->x_sequence_vec+ii, tmp_ptr);
    tmp_ptr += dimx*sizeof(REAL);
    VECCSE(dimx, 0.0, ocp->x_sequence_vec+ii, 0);
  }

  tmp_ptr = c_ptr;
  c_ptr += SIZE_STRVEC(N*dimx);
  for (ii=0; ii<N; ++ii) {
    CREATE_STRVEC(dimx, ocp->lmd_sequence_vec+ii, tmp_ptr);
    tmp_ptr += dimx*sizeof(REAL);
    VECCSE(dimx, 0.0, ocp->lmd_sequence_vec+ii, 0);
  }

  ocp->dimx = ocp->model->dimx;
  ocp->dimu = ocp->model->dimu;
  ocp->dimc = ocp->model->dimc;
  ocp->dimuc = ocp->model->dimu + ocp->model->dimc;
  ocp->N = N;
  ocp->dim_solution = N*(ocp->model->dimu+ocp->model->dimc);
  ocp->memsize = SINGLE_SHOOTING_OCP_MEMSIZE(ocp, N);

#if defined(RUNTIME_CHECKS)
  if(c_ptr > ((char *) mem) + ocp->memsize) {
    printf("\nerror: SINGLE_SHOOTING_OCP_CREATE: outside memory bounds!\n\n");
    exit(1);
  }
#endif
}


void SINGLE_SHOOTING_OCP_COMPUTE_OPTIMALITY_RESIDUAL( 
    struct SINGLE_SHOOTING_OCP *ocp, REAL current_time, 
    struct STRVEC *current_state, struct STRVEC *solution, 
    struct STRVEC *optimality_residual) {
  struct STRVEC *x_sequence_vec = ocp->x_sequence_vec;
  struct STRVEC *lmd_sequence_vec = ocp->lmd_sequence_vec;
  struct STRVEC *dx_vec = ocp->dx_vec;
  REAL horizon_length = TIME_VARYING_SMOOTH_HORIZON_GET_LENGTH(ocp->horizon, 
                                                               current_time);
  REAL delta_tau = horizon_length / ocp->N;

  // Compute the state trajectory over the horizon on the basis of the 
  // time, solution_vec and the state_vec.
  NMPC_MODEL_F(ocp->model, current_time, current_state->pa, solution->pa, 
               dx_vec->pa);
  AXPY(ocp->dimx, delta_tau, dx_vec, current_state, 0, x_sequence_vec, 0);
  REAL tau = time + delta_tau;
  for (int i=1; i<ocp->N; ++i, tau+=delta_tau) {
    NMPC_MODEL_F(ocp->model, tau, (x_sequence_vec+(i-1))->pa, 
                 solution->pa+i*ocp->dimuc, dx_vec->pa);
    AXPY(ocp->dimx, delta_tau, dx_vec, x_sequence_vec, (i-1)*ocp->dimx, 
         x_sequence_vec, i*ocp->dimx);
  }

  // Compute the Lagrange multiplier over the horizon on the basis of 
  // time, solution_vec and the state_vec.
  NMPC_MODEL_PHIX(ocp->model, tau, (x_sequence_vec+(ocp->N-1))->pa, 
                  (lmd_sequence_vec+(ocp->N-1))->pa);
  for (int i=ocp->N-1; i>=1; --i, tau-=delta_tau) {
    NMPC_MODEL_HX(ocp->model, tau, (x_sequence_vec+(i-1))->pa, 
                  solution->pa+i*ocp->dimuc, (lmd_sequence_vec+(i-1))->pa, 
                  dx_vec->pa);
    AXPY(ocp->dimx, delta_tau, dx_vec, lmd_sequence_vec, i*ocp->dimx, 
         lmd_sequence_vec, (i-1)*ocp->dimx);
  }

  // Compute the erros in optimality over the horizon on the basis of the 
  // control_input_vec and the state_vec.
  tau = current_time,;
  NMPC_MODEL_HU(ocp->model, tau, current_state->pa, solution->pa, 
                lmd_sequence_vec->pa, optimality_residual->pa);
  for (int i=1; i<ocp->N; ++i, tau+=delta_tau) {
    NMPC_MODEL_HU(ocp->model, tau, (x_sequence_vec+(i-1))->pa, 
                  solution->pa+i*ocp->dimuc, (lmd_sequence_vec+(i-1))->pa, 
                  optimality_residual->pa+i*ocp->dimuc);
  }
}


void SINGLE_SHOOTING_OCP_PREDICT_STATE_FROM_SOLUTION( 
    struct SINGLE_SHOOTING_OCP *ocp, REAL current_time, 
    struct STRVEC *current_state, struct STRVEC *solution, 
    REAL prediction_length, struct STRVEC *predicted_state) {
  NMPC_MODEL_F(ocp->model, current_time, current_state->pa, solution->pa, 
               ocp->dx_vec->pa);
  AXPY(ocp->dimx, prediction_length, ocp->dx_vec, current_state, 
       predicted_state);
}


void SINGLE_SHOOTING_OCP_RESET_HORIZON_LENGTH( 
    struct SINGLE_SHOOTING_OCP *ocp, REAL initial_time) {
  TIME_VARYING_SMOOTH_HORIZON_RESET_LENGTH(ocp->horizon, initial_time);
}