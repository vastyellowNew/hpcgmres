int SINGLE_SHOOTING_CGMRES_STRSIZE() {
  return sizeof(struct SINGLE_SHOOTING_CGMRES);
}


int SINGLE_SHOOTING_CGMRES_MEMSIZE(struct SINGLE_SHOOTING_CGMRES *cgmres, 
                                   int N, int kmax) {
  int size = 0;

  int dimu = cgmres->continuation->ocp->dimu;
  int dimc = cgmres->continuation->ocp->dimc;
  int dim_solution = N * (dimu+dimc);
  size += SINGLE_SHOOTING_CONTINUATION_MEMSIZE(cgmres->continuation, N);
  size += MFGMRES_MEMSIZE(dim_solution, kmax);
  size += CGMRES_INITIALIZER_MEMSIZE(cgmres->initializer, kmax);

  // size of pointers and int
  size += 3*sizeof(struct STRVEC); 

  // size of matrices and vectors
  size += 2*SIZE_STRVEC(dim_solution); // solution, solution_update
  size += 1*SIZE_STRVEC(dimu+dimc); // initial_solution

  size = (size+63)/64*64; // make multiple of typical cache line size
  size += 64; // align to typical cache line size

  return size;
}


void SINGLE_SHOOTING_CGMRES_CREATE(struct SINGLE_SHOOTING_CGMRES *cgmres, 
                                   REAL T_f, double alpha, REAL initial_time, 
                                   int N, REAL finite_difference_increment,
                                   REAL zeta, int kmax, void *memory) {
  int dimu = cgmres->continuation->ocp->dimu;
  int dimc = cgmres->continuation->ocp->dimc;
  int dimuc = dimu + dimc;
  int dim_solution = N * dimuc;

  SINGLE_SHOOTING_CONTINUATION_CREATE(cgmres->continuation, T_f, alpha, 
                                      initial_time, N, 
                                      finite_difference_increment, zeta,
                                      memory);
  MFGMRES_CREATE(cgmres->mfgmres, dim_solution, kmax, memory);
  CGMRES_INITIALIZER_CREATE(cgmres->initializer, finite_difference_increment, 
                            kmax, memory);

  // zero memory (to avoid corrupted memory like e.g. NaN)
  int memsize =  SINGLE_SHOOTING_CGMRES_MEMSIZE(cgmres, N, kmax);
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
  cgmres->solution = sv_ptr;
  sv_ptr += 1;
  cgmres->solution_update = sv_ptr;
  sv_ptr += 1;
  cgmres->initial_solution = sv_ptr;
  sv_ptr += 1;

  // align to typical cache line size
  size_t s_ptr = (size_t) sv_ptr;
  s_ptr = (s_ptr+63)/64*64;

  // void stuff
  char *c_ptr = (char *) s_ptr;

  CREATE_STRVEC(dim_solution, cgmres->solution, c_ptr);
  c_ptr += cgmres->solution->memsize;
  VECCSE(dim_solution, 0.0, cgmres->solution, 0);

  CREATE_STRVEC(dim_solution, cgmres->solution_update, c_ptr);
  c_ptr += cgmres->solution_update->memsize;
  VECCSE(dim_solution, 0.0, cgmres->solution_update, 0);

  CREATE_STRVEC(dimuc, cgmres->initial_solution, c_ptr);
  c_ptr += cgmres->initial_solution->memsize;
  VECCSE(dimuc, 0.0, cgmres->initial_solution, 0);

  cgmres->dimu = dimu;
  cgmres->dimc = dimc;
  cgmres->dimuc = dimuc;
  cgmres->N = N;
  cgmres->memsize = SINGLE_SHOOTING_CGMRES_MEMSIZE(cgmres, N, kmax);

#if defined(RUNTIME_CHECKS)
  if(c_ptr > ((char *) mem) + mfgmres->memsize) {
    printf("\nerror: CGMRES_INITIALIZER_CREATE: outside memory bounds!\n\n");
    exit(1);
  }
#endif
}


void SINGLE_SHOOTING_CGMRES_SET_INITIALIZATION_PARAMETERS(
    struct SINGLE_SHOOTING_CGMRES *cgmres, 
    struct STRVEC *initial_guess_solution, REAL newton_residual_tolerance, 
    int max_newton_iteration) {
  CGMRES_INITIALIZER_SET_INITIAL_GUESS_SOLUTION(cgmres->initializer, 
                                                initial_guess_solution);
  CGMRES_INITIALIZER_SET_TERMINATION_CRITERIONS(cgmres->initializer, 
                                                newton_residual_tolerance, 
                                                max_newton_iteration);
}


void SINGLE_SHOOTING_CGMRES_INITIALIZE_SOLUTION(
    struct SINGLE_SHOOTING_CGMRES *cgmres, REAL initial_time, 
    struct blasfeo_dvec *initial_state) {
  CGMRES_INITIALIZER_COMPUTE_INITIAL_SOLUTION(cgmres->initializer, initial_time, 
                                              initial_state, 
                                              cgmres->initial_solution);
  for (int i=0; i<cgmres->N; ++i) {
    VECCP(cgmres->dimuc, cgmres->initial_solution, 0, cgmres->solution, 
          i*cgmres->dimuc);
  }
  SINGLE_SHOOTING_CONTINUATION_RESET_HORIZON_LENGTH(cgmres->continuation, 
                                                    initial_time)
}


void SINGLE_SHOOTING_CGMRES_UPDATE_CONTROL_INPUT(
    struct SINGLE_SHOOTING_CGMRES *cgmres, REAL current_time,
    struct STRVEC *current_state, REAL sampling_period,
    struct STRVEC *control_input) {
  cgmres->mfgmres_args->current_time = current_time;
  cgmres->mfgmres_args->current_state_ptr = current_state;
  cgmres->mfgmres_args->current_solution_ptr = cgmres->solution;
  MFGMRES_SOLVE_LINEAR_PROBLEM(cgmres->mfgmres, cgmres->continuation, 
                               cgmres->mfgmres_args, cgmres->solution_update);
  SINGLE_SHOOTING_CONTINUATION_INTEGRATE_SOLUTION(cgmres->continuation, 
                                                  cgmres->solution, 
                                                  cgmres->solution_update,
                                                  sampling_period);
  VECCP(cgmres->dimu, cgmres->solution, 0, control_input, 0);
}


void SINGLE_SHOOTING_CGMRES_GET_CONTROL_INPUT(
    struct SINGLE_SHOOTING_CGMRES *cgmres, struct blasfeo_dvec *control_input) {
  VECCP(cgmres->dimu, cgmres->solution, 0, control_input, 0);
}


REAL SINGLE_SHOOTING_CGMRES_GET_ERROR_NORM(
    struct SINGLE_SHOOTING_CGMRES *cgmres, REAL current_time,
    struct blasfeo_dvec *current_state) {
  return SINGLE_SHOOTING_CONTINUATION_GET_ERROR_NORM(cgmres->continuation, 
                                                     current_time, 
                                                     current_state,
                                                     cgmres->solution);
}