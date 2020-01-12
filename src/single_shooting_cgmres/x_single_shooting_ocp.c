int SINGLE_SHOOTING_OCP_STRSIZE() {
  return sizeof(struct SINGLE_SHOOTING_OCP);
}


int SINGLE_SHOOTING_OCP_MEMSIZE(int N) {
  int size = 0;

  size += 2*N*sizeof(REAL*); // x_sequence, lmd_sequence
  size += 2*N*NMPC_MODEL_DIMX()*sizeof(REAL); // x_sequence, lmd_sequence
  size += NMPC_MODEL_DIMX()*sizeof(REAL); // dx

  size = (size+63)/64*64; // make multiple of typical cache line size
  size += 64; // align to typical cache line size

  return size;
}


void SINGLE_SHOOTING_OCP_CREATE(struct SINGLE_SHOOTING_OCP *ocp, REAL T_f, 
                                REAL alpha, REAL initial_time, int N) {
  NMPC_MODEL_CREATE(&ocp->model);
  TIME_VARYING_SMOOTH_HORIZON_CREATE(&ocp->horizon, T_f, alpha, initial_time);
  ocp->x_sequence = ALLOCATE_MAT(N, NMPC_MODEL_DIMX());
  ocp->lmd_sequence = ALLOCATE_MAT(N, NMPC_MODEL_DIMX());
  ocp->dx = ALLOCATE_VEC(NMPC_MODEL_DIMX());
  ocp->dimx = NMPC_MODEL_DIMX();
  ocp->dimu = NMPC_MODEL_DIMU();
  ocp->dimc = NMPC_MODEL_DIMC();
  ocp->dimuc = ocp->dimu + ocp->dimc;
  ocp->N = N;
  ocp->dim_solution = N * (ocp->dimu+ocp->dimc);
  ocp->memsize = SINGLE_SHOOTING_OCP_MEMSIZE(N);
}


void SINGLE_SHOOTING_OCP_DELETE(struct SINGLE_SHOOTING_OCP *ocp) {
  FREE_MAT(ocp->x_sequence);
  FREE_MAT(ocp->lmd_sequence);
  FREE_VEC(ocp->dx);
}


void SINGLE_SHOOTING_OCP_COMPUTE_OPTIMALITY_RESIDUAL( 
    struct SINGLE_SHOOTING_OCP *ocp, REAL current_time, REAL *current_state, 
    REAL *solution, REAL *optimality_residual) {
  REAL horizon_length = TIME_VARYING_SMOOTH_HORIZON_GET_LENGTH(&ocp->horizon, 
                                                               current_time);
  REAL delta_tau = horizon_length / ocp->N;

  // Compute the state trajectory over the horizon on the basis of the 
  // time, solution_vec and the state_vec.
  NMPC_MODEL_F(&ocp->model, current_time, current_state, solution, ocp->dx);
  AXPY(ocp->dimx, delta_tau, ocp->dx, current_state, ocp->x_sequence[0]);
  REAL tau = current_time + delta_tau;
  int ii;
  for (ii=1; ii<ocp->N; ++ii, tau+=delta_tau) {
    NMPC_MODEL_F(&ocp->model, tau, ocp->x_sequence[ii-1], 
                 &(solution[ii*ocp->dimuc]), ocp->dx);
    AXPY(ocp->dimx, delta_tau, ocp->dx, ocp->x_sequence[ii-1], 
         ocp->x_sequence[ii]);
  }

  // Compute the Lagrange multiplier over the horizon on the basis of 
  // time, solution_vec and the state_vec.
  NMPC_MODEL_PHIX(&ocp->model, tau, ocp->x_sequence[ocp->N-1], 
                  ocp->lmd_sequence[ocp->N-1]);
  tau = horizon_length - delta_tau;
  for (ii=ocp->N-1; ii>=1; --ii, tau-=delta_tau) {
    NMPC_MODEL_HX(&ocp->model, tau, ocp->x_sequence[ii-1], 
                  &(solution[ii*ocp->dimuc]), ocp->lmd_sequence[ii],
                  ocp->dx);
    AXPY(ocp->dimx, delta_tau, ocp->dx, ocp->lmd_sequence[ii], 
         ocp->lmd_sequence[ii-1]);
  }

  // Compute the erros in optimality over the horizon on the basis of the 
  // control_input_vec and the state_vec.
  NMPC_MODEL_HU(&ocp->model, current_time, current_state, solution, 
                ocp->lmd_sequence[0], optimality_residual);
  tau = current_time + delta_tau;
  for (ii=1; ii<ocp->N; ++ii, tau+=delta_tau) {
    NMPC_MODEL_HU(&ocp->model, tau, ocp->x_sequence[ii-1], 
                  &(solution[ii*ocp->dimuc]), ocp->lmd_sequence[ii],
                  &(optimality_residual[ii*ocp->dimuc]));
  }
}


void SINGLE_SHOOTING_OCP_PREDICT_STATE_FROM_SOLUTION( 
    struct SINGLE_SHOOTING_OCP *ocp, REAL current_time, REAL *current_state, 
    REAL *solution, REAL prediction_length, REAL *predicted_state) {
  NMPC_MODEL_F(&ocp->model, current_time, current_state, solution, ocp->dx);
  AXPY(ocp->dimx, prediction_length, ocp->dx, current_state, predicted_state);
}


void SINGLE_SHOOTING_OCP_RESET_HORIZON_LENGTH( 
    struct SINGLE_SHOOTING_OCP *ocp, REAL initial_time) {
  TIME_VARYING_SMOOTH_HORIZON_RESET_LENGTH(&ocp->horizon, initial_time);
}


int SINGLE_SHOOTING_OCP_DIMX() {
  return NMPC_MODEL_DIMX();
}


int SINGLE_SHOOTING_OCP_DIMU() {
  return NMPC_MODEL_DIMU();
}


int SINGLE_SHOOTING_OCP_DIMC() {
  return NMPC_MODEL_DIMC();
}