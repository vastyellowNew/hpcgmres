int SINGLE_SHOOTING_CONTINUATION_STRSIZE() {
  return sizeof(struct SINGLE_SHOOTING_CONTINUATION);
}


int SINGLE_SHOOTING_CONTINUATION_MEMSIZE(
    struct SINGLE_SHOOTING_CONTINUATION *continuation, int N) {
  int size = 0;

  size += SINGLE_SHOOTING_OCP_MEMSIZE(continuation->ocp, N);

  // size of pointers and int
  size += 5*sizeof(struct STRVEC); 

  int dimx = continuation->ocp->model->dimu;
  int dimu = continuation->ocp->model->dimu;
  int dimc = continuation->ocp->model->dimc;
  int dim_solution = N * (dimu+dimc);

  // size of matrices and vectors
  size += 1*SIZE_STRVEC(dimx); 
  size += 4*SIZE_STRVEC(dim_solution); 

  size = (size+63)/64*64; // make multiple of typical cache line size
  size += 64; // align to typical cache line size

  return size;
}


void SINGLE_SHOOTING_CONTINUATION_CREATE(
    struct SINGLE_SHOOTING_CONTINUATION *continuation, REAL T_f, REAL alpha, 
    REAL initial_time, int N, REAL finite_difference_increment, REAL zeta, 
    void *memory) {

  SINGLE_SHOOTING_OCP_CREATE(continuation->ocp, T_f, alpha, initial_time, N, 
                             memory);

  // zero memory (to avoid corrupted memory like e.g. NaN)
  int memsize = SINGLE_SHOOTING_OCP_MEMSIZE(continuation, N);
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
  continuation->incremented_state = sv_ptr;
  sv_ptr += 1;
  continuation->incremented_solution = sv_ptr;
  sv_ptr += 1;
  continuation->optimality_residual = sv_ptr;
  sv_ptr += 1;
  continuation->optimality_residual1 = sv_ptr;
  sv_ptr += 1;
  continuation->optimality_residual2 = sv_ptr;
  sv_ptr += 1;

  // align to typical cache line size
  size_t s_ptr = (size_t) sv_ptr;
  s_ptr = (s_ptr+63)/64*64;

  // void stuff
  char *c_ptr = (char *) s_ptr;

  int dimx = continuation->ocp->model->dimx;
  int dimu = continuation->ocp->model->dimu;
  int dimc = continuation->ocp->model->dimc;
  int dimuc = dimu + dimc;
  int dim_solution = N * dimuc;

  CREATE_STRVEC(dim_solution, continuation->incremented_state, c_ptr);
  c_ptr += continuation->incremented_state->memsize;
  VECCSE(dimx, 0.0, continuation->incremented_state, 0);

  CREATE_STRVEC(dim_solution, continuation->incremented_solution, c_ptr);
  c_ptr += continuation->incremented_solution->memsize;
  VECCSE(dim_solution, 0.0, continuation->incremented_solution, 0);

  CREATE_STRVEC(dim_solution, continuation->optimality_residual, c_ptr);
  c_ptr += continuation->optimality_residual->memsize;
  VECCSE(dim_solution, 0.0, continuation->optimality_residual, 0);

  CREATE_STRVEC(dim_solution, continuation->optimality_residual1, c_ptr);
  c_ptr += continuation->optimality_residual1->memsize;
  VECCSE(dim_solution, 0.0, continuation->optimality_residual1, 0);

  CREATE_STRVEC(dim_solution, continuation->optimality_residual2, c_ptr);
  c_ptr += continuation->optimality_residual2->memsize;
  VECCSE(dim_solution, 0.0, continuation->optimality_residual2, 0);

  continuation->dimx = continuation->ocp->model->dimx;
  continuation->dimu = continuation->ocp->model->dimu;
  continuation->dimc = continuation->ocp->model->dimc;
  continuation->N = N;
  continuation->dim_solution = dim_solution;
  continuation->finite_difference_increment = finite_difference_increment;
  continuation->zeta = zeta;
  continuation->incremented_time = 0.0;
  continuation->memsize = SINGLE_SHOOTING_CONTINUATION_MEMSIZE(continuation, N);

#if defined(RUNTIME_CHECKS)
  if(c_ptr > ((char *) mem) + ocp->memsize) {
    printf("\nerror: SINGLE_SHOOTING_OCP_CREATE: outside memory bounds!\n\n");
    exit(1);
  }
#endif
}


void SINGLE_SHOOTING_CONTINUATION_INTEGRATE_SOLUTION(
    struct SINGLE_SHOOTING_CONTINUATION *continuation, struct STRVEC *solution, 
    struct STRVEC *solution_update, REAL integration_length) {
  VECAD(continuation->dim_solution, integration_length, solution_update, 0, 
        solution, 0)
}


void SINGLE_SHOOTING_CONTINUATION_COMPUTE_B(
    struct SINGLE_SHOOTING_CONTINUATION *continuation, 
    struct SINGLE_SHOOTING_CONTINUATION_MFGMRES_ARGS *args, 
    struct STRVEC *direction, struct STRVEC *b) {
  continuation->incremented_time 
      = args->current_time + continuation->finite_difference_increment;
  SINGLE_SHOOTING_OCP_PREDICT_STATE_FROM_SOLUTION(
      continuation->ocp, args->current_time, args->current_state_ptr, 
      args->current_solution_ptr, continuation->finite_difference_increment,
      continuation->incremented_state);
  AXPY(continuation->dim_solution, continuation->finite_difference_increment, 
       direction, args->current_solution_ptr, 0, 
      continuation->incremented_solution, 0);
  SINGLE_SHOOTING_OCP_COMPUTE_OPTIMALITY_RESIDUAL(
      continuation->ocp, args->current_time, args->current_state_ptr, 
      args->current_solution_ptr, continuation->optimality_residual);
  SINGLE_SHOOTING_OCP_COMPUTE_OPTIMALITY_RESIDUAL(
      continuation->ocp, continuation->incremented_time,
      continuation->incremented_state, args->current_solution_ptr, 
      continuation->optimality_residual1);
  SINGLE_SHOOTING_OCP_COMPUTE_OPTIMALITY_RESIDUAL(
      continuation->ocp, continuation->incremented_time,
      continuation->incremented_state, continuation->incremented_solution,  
      continuation->optimality_residual2);
  AXBPY(continuation->dim_solution, 
        1/continuation->finite_difference_increment-continuation->zeta, 
        continuation->optimality_residual, 0,
        -1/continuation->finite_difference_increment, 
        continuation->optimality_residual2, 0, b, 0); 
}


void SINGLE_SHOOTING_CONTINUATION_COMPUTE_AX(
    struct SINGLE_SHOOTING_CONTINUATION *continuation, 
    struct SINGLE_SHOOTING_CONTINUATION_MFGMRES_ARGS *args, 
    struct STRVEC *direction, struct STRVEC *ax) {
  AXPY(continuation->dim_solution, continuation->finite_difference_increment, 
       direction, args->current_solution_ptr, 0, 
       continuation->incremented_solution, 0);
  SINGLE_SHOOTING_OCP_COMPUTE_OPTIMALITY_RESIDUAL(
      continuation->ocp, continuation->incremented_time,
      continuation->incremented_state, continuation->incremented_solution,  
      continuation->optimality_residual2);
  AXBPY(continuation->dim_solution, 1/continuation->finite_difference_increment, 
        continuation->optimality_residual2, 0, 
        -1/continuation->finite_difference_increment, 
        continuation->optimality_residual1, 0, ax, 0); 
}


REAL SINGLE_SHOOTING_CONTINUATION_GET_ERROR_NORM(
    struct SINGLE_SHOOTING_CONTINUATION *continuation, REAL current_time, 
    struct STRVEC *current_state, struct STRVEC *current_solution) {
  SINGLE_SHOOTING_OCP_COMPUTE_OPTIMALITY_RESIDUAL(
      continuation->ocp, current_time, current_state, current_solution,
      continuation->optimality_residual);
  return sqrt(VECDOT(continuation->dim_solution, 
                     continuation->optimality_residual, 0,
                     continuation->optimality_residual, 0));
}


void SINGLE_SHOOTING_CONTINUATION_RESET_HORIZON_LENGTH(
    struct SINGLE_SHOOTING_CONTINUATION *continuation, REAL initial_time) {
  SINGLE_SHOOTING_OCP_RESET_HORIZON_LENGTH(continuation->ocp, initial_time);
}