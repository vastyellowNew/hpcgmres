int SINGLE_SHOOTING_CGMRES_STRSIZE() {
  return sizeof(struct SINGLE_SHOOTING_CGMRES);
}


int SINGLE_SHOOTING_CGMRES_MEMSIZE(int N, int kmax) {
  int dimx = SINGLE_SHOOTING_CONTINUATION_DIMX();
  int dimu = SINGLE_SHOOTING_CONTINUATION_DIMU();
  int dimc = SINGLE_SHOOTING_CONTINUATION_DIMC();
  int dim_solution = N * (dimu+dimc);
  int size = 0;
  size += MFGMRES_FOR_SINGLE_SHOOTING_CGMRES_MEMSIZE(dim_solution, kmax);
  size += SINGLE_SHOOTING_CONTINUATION_MEMSIZE(N);
  size += CGMRES_INITIALIZER_MEMSIZE(kmax);
  size += 2*dim_solution*sizeof(REAL); // solution, solution_update
  size += 1*(dimu+dimc)*sizeof(REAL); // initial_solution
  size = (size+63)/64*64; // make multiple of typical cache line size
  size += 64; // align to typical cache line size
  return size;
}


void SINGLE_SHOOTING_CGMRES_CREATE(struct SINGLE_SHOOTING_CGMRES *cgmres, 
                                   REAL T_f, REAL alpha, REAL initial_time, 
                                   int N, REAL finite_difference_increment,
                                   REAL zeta, int kmax) {
  int dimx = SINGLE_SHOOTING_CONTINUATION_DIMX();
  int dimu = SINGLE_SHOOTING_CONTINUATION_DIMU();
  int dimc = SINGLE_SHOOTING_CONTINUATION_DIMC();
  int dim_solution = N * (dimu+dimc);
  SINGLE_SHOOTING_CONTINUATION_CREATE(&cgmres->continuation, T_f, alpha, 
                                      initial_time, N, 
                                      finite_difference_increment, zeta);
  MFGMRES_FOR_SINGLE_SHOOTING_CGMRES_CREATE(&cgmres->mfgmres, dim_solution, 
                                            kmax);
  CGMRES_INITIALIZER_CREATE(&cgmres->initializer, finite_difference_increment, 
                            kmax);
  cgmres->solution = ALLOCATE_VEC(dim_solution);
  cgmres->solution_update = ALLOCATE_VEC(dim_solution);
  cgmres->initial_solution = ALLOCATE_VEC(dimu+dimc);
  cgmres->dimu = dimu;
  cgmres->dimc = dimc;
  cgmres->dimuc = dimu + dimc;
  cgmres->N = N;
  cgmres->memsize = SINGLE_SHOOTING_CGMRES_MEMSIZE(N, kmax);
}


void SINGLE_SHOOTING_CGMRES_DELETE(struct SINGLE_SHOOTING_CGMRES *cgmres) {
  SINGLE_SHOOTING_CONTINUATION_DELETE(&cgmres->continuation);
  MFGMRES_FOR_SINGLE_SHOOTING_CGMRES_DELETE(&cgmres->mfgmres);
  CGMRES_INITIALIZER_DELETE(&cgmres->initializer);
  FREE_VEC(cgmres->solution);
  FREE_VEC(cgmres->solution_update);
  FREE_VEC(cgmres->initial_solution);
}


void SINGLE_SHOOTING_CGMRES_SET_INITIALIZATION_PARAMETERS(
    struct SINGLE_SHOOTING_CGMRES *cgmres, REAL newton_residual_tolerance, 
    int max_newton_iteration, REAL *initial_guess_solution) {
  CGMRES_INITIALIZER_SET_TERMINATION_CRITERIONS(&cgmres->initializer, 
                                                newton_residual_tolerance, 
                                                max_newton_iteration);
  CGMRES_INITIALIZER_SET_INITIAL_GUESS_SOLUTION(&cgmres->initializer, 
                                                initial_guess_solution);
}


void SINGLE_SHOOTING_CGMRES_INITIALIZE_SOLUTION(
    struct SINGLE_SHOOTING_CGMRES *cgmres, REAL initial_time, 
    REAL *initial_state) {
  CGMRES_INITIALIZER_COMPUTE_INITIAL_SOLUTION(&cgmres->initializer, 
                                              initial_time, initial_state, 
                                              cgmres->initial_solution);
  for (int i=0; i<cgmres->N; ++i) {
    VECCP(cgmres->dimuc, cgmres->initial_solution, 
          &(cgmres->solution[i*cgmres->dimuc]));
  }
  SINGLE_SHOOTING_CONTINUATION_RESET_HORIZON_LENGTH(&cgmres->continuation, 
                                                    initial_time);
}


void SINGLE_SHOOTING_CGMRES_UPDATE_CONTROL_INPUT(
    struct SINGLE_SHOOTING_CGMRES *cgmres, REAL current_time,
    REAL *current_state, REAL sampling_period, REAL *control_input) {
  cgmres->mfgmres_args.current_time = current_time;
  cgmres->mfgmres_args.current_state_ptr = current_state;
  cgmres->mfgmres_args.current_solution_ptr = cgmres->solution;
  MFGMRES_FOR_SINGLE_SHOOTING_CGMRES_SOLVE_LINEAR_PROBLEM(
      &cgmres->mfgmres, &cgmres->continuation, &cgmres->mfgmres_args, 
      cgmres->solution_update);
  SINGLE_SHOOTING_CONTINUATION_INTEGRATE_SOLUTION(&cgmres->continuation, 
                                                  cgmres->solution, 
                                                  cgmres->solution_update,
                                                  sampling_period);
  VECCP(cgmres->dimu, cgmres->solution, control_input);
}


void SINGLE_SHOOTING_CGMRES_GET_CONTROL_INPUT(
    struct SINGLE_SHOOTING_CGMRES *cgmres, REAL *control_input) {
  VECCP(cgmres->dimu, cgmres->solution, control_input);
}


REAL SINGLE_SHOOTING_CGMRES_GET_ERROR_NORM(
    struct SINGLE_SHOOTING_CGMRES *cgmres, REAL current_time,
    REAL *current_state) {
  return SINGLE_SHOOTING_CONTINUATION_GET_ERROR_NORM(&cgmres->continuation, 
                                                     current_time, 
                                                     current_state,
                                                     cgmres->solution);
}