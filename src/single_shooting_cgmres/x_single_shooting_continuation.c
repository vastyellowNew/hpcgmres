int SINGLE_SHOOTING_CONTINUATION_STRSIZE() {
  return sizeof(struct SINGLE_SHOOTING_CONTINUATION);
}


int SINGLE_SHOOTING_CONTINUATION_MEMSIZE(int N) {
  int dimx = SINGLE_SHOOTING_OCP_DIMX();
  int dimu = SINGLE_SHOOTING_OCP_DIMU();
  int dimc = SINGLE_SHOOTING_OCP_DIMC();
  int dim_solution = N * (dimu+dimc);

  int size = 0;

  size += SINGLE_SHOOTING_OCP_MEMSIZE(N);
  size += dimx*sizeof(REAL); 
  size += 4*dim_solution*sizeof(REAL); 

  size = (size+63)/64*64; // make multiple of typical cache line size
  size += 64; // align to typical cache line size

  return size;
}


void SINGLE_SHOOTING_CONTINUATION_CREATE(
    struct SINGLE_SHOOTING_CONTINUATION *continuation, REAL T_f, REAL alpha, 
    REAL initial_time, int N, REAL finite_difference_increment, REAL zeta) {
  int dimx = SINGLE_SHOOTING_OCP_DIMX();
  int dimu = SINGLE_SHOOTING_OCP_DIMU();
  int dimc = SINGLE_SHOOTING_OCP_DIMC();
  int dim_solution = N * (dimu+dimc);

  SINGLE_SHOOTING_OCP_CREATE(&continuation->ocp, T_f, alpha, initial_time, N); 
  continuation->incremented_state = ALLOCATE_VEC(dimx);
  continuation->incremented_solution = ALLOCATE_VEC(dim_solution);
  continuation->optimality_residual = ALLOCATE_VEC(dim_solution);
  continuation->optimality_residual1 = ALLOCATE_VEC(dim_solution);
  continuation->optimality_residual2 = ALLOCATE_VEC(dim_solution);
  continuation->finite_difference_increment = finite_difference_increment; 
  continuation->zeta = zeta;
  continuation->incremented_time = 0;
  continuation->dimx = dimx;
  continuation->dimu = dimu;
  continuation->dimc = dimc;
  continuation->N = N;
  continuation->dim_solution = dim_solution;
  continuation->memsize = SINGLE_SHOOTING_CONTINUATION_MEMSIZE(N);
}


void SINGLE_SHOOTING_CONTINUATION_DELETE(
    struct SINGLE_SHOOTING_CONTINUATION *continuation) {
  SINGLE_SHOOTING_OCP_DELETE(&continuation->ocp);
  FREE_VEC(continuation->incremented_state);
  FREE_VEC(continuation->incremented_solution);
  FREE_VEC(continuation->optimality_residual);
  FREE_VEC(continuation->optimality_residual1);
  FREE_VEC(continuation->optimality_residual2);
}


void SINGLE_SHOOTING_CONTINUATION_INTEGRATE_SOLUTION(
    struct SINGLE_SHOOTING_CONTINUATION *continuation, REAL *solution, 
    REAL *solution_update, REAL integration_length) {
  VECMAD(continuation->dim_solution, integration_length, solution_update, 
         solution);
}


void SINGLE_SHOOTING_CONTINUATION_COMPUTE_B(
    struct SINGLE_SHOOTING_CONTINUATION *continuation, 
    struct SINGLE_SHOOTING_CONTINUATION_MFGMRES_ARGS *args, REAL *direction, 
    REAL *b) {
  continuation->incremented_time 
      = args->current_time + continuation->finite_difference_increment;
  SINGLE_SHOOTING_OCP_PREDICT_STATE_FROM_SOLUTION(
      &continuation->ocp, args->current_time, args->current_state_ptr, 
      args->current_solution_ptr, continuation->finite_difference_increment,
      continuation->incremented_state);
  AXPY(continuation->dim_solution, continuation->finite_difference_increment, 
       direction, args->current_solution_ptr, 
       continuation->incremented_solution);
  SINGLE_SHOOTING_OCP_COMPUTE_OPTIMALITY_RESIDUAL(
      &continuation->ocp, args->current_time, args->current_state_ptr, 
      args->current_solution_ptr, continuation->optimality_residual);
  SINGLE_SHOOTING_OCP_COMPUTE_OPTIMALITY_RESIDUAL(
      &continuation->ocp, continuation->incremented_time,
      continuation->incremented_state, args->current_solution_ptr, 
      continuation->optimality_residual1);
  SINGLE_SHOOTING_OCP_COMPUTE_OPTIMALITY_RESIDUAL(
      &continuation->ocp, continuation->incremented_time,
      continuation->incremented_state, continuation->incremented_solution,  
      continuation->optimality_residual2);
  AXPBY(continuation->dim_solution, 
        1/continuation->finite_difference_increment-continuation->zeta, 
        continuation->optimality_residual, 
        -1/continuation->finite_difference_increment, 
        continuation->optimality_residual2, b); 
}


void SINGLE_SHOOTING_CONTINUATION_COMPUTE_AX(
    struct SINGLE_SHOOTING_CONTINUATION *continuation, 
    struct SINGLE_SHOOTING_CONTINUATION_MFGMRES_ARGS *args, REAL *direction, 
    REAL *ax) {
  AXPY(continuation->dim_solution, continuation->finite_difference_increment, 
       direction, args->current_solution_ptr, 
       continuation->incremented_solution);
  SINGLE_SHOOTING_OCP_COMPUTE_OPTIMALITY_RESIDUAL(
      &continuation->ocp, continuation->incremented_time,
      continuation->incremented_state, continuation->incremented_solution,  
      continuation->optimality_residual2);
  AXPBY(continuation->dim_solution, 1/continuation->finite_difference_increment, 
        continuation->optimality_residual2,  
        -1/continuation->finite_difference_increment, 
        continuation->optimality_residual1, ax); 
}


REAL SINGLE_SHOOTING_CONTINUATION_GET_ERROR_NORM(
    struct SINGLE_SHOOTING_CONTINUATION *continuation, REAL current_time, 
    REAL *current_state, REAL *current_solution) {
  SINGLE_SHOOTING_OCP_COMPUTE_OPTIMALITY_RESIDUAL(
      &continuation->ocp, current_time, current_state, current_solution,
      continuation->optimality_residual);
  return sqrt(VECNRM2(continuation->dim_solution, 
                      continuation->optimality_residual));
}


void SINGLE_SHOOTING_CONTINUATION_RESET_HORIZON_LENGTH(
    struct SINGLE_SHOOTING_CONTINUATION *continuation, REAL initial_time) {
  SINGLE_SHOOTING_OCP_RESET_HORIZON_LENGTH(&continuation->ocp, initial_time);
}


int SINGLE_SHOOTING_CONTINUATION_DIMX() {
  return SINGLE_SHOOTING_OCP_DIMX();
}


int SINGLE_SHOOTING_CONTINUATION_DIMU() {
  return SINGLE_SHOOTING_OCP_DIMU();
}


int SINGLE_SHOOTING_CONTINUATION_DIMC() {
  return SINGLE_SHOOTING_OCP_DIMC();
}