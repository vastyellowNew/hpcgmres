int CGMRES_INITIALIZER_STRSIZE() {
  return sizeof(struct CGMRES_INITIALIZER);
}


int CGMRES_INITIALIZER_MEMSIZE(struct CGMRES_INITIALIZER *initializer, 
                               int kmax) {
  int dimu = initializer->newton->ocp->model->dimu;
  int dimc = initializer->newton->ocp->model->dimc;
  int dimuc = dimu + dimc;

  int size = 0;

  size += INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_MEMSIZE(initializer->newton);
  size += MFGMRES_MEMSIZE(dimu+dimc, kmax);

  // size of pointers and int
  size += 2*sizeof(struct STRVEC); // 

  // size of matrices and vectors
  size += 2*SIZE_STRVEC(dimuc); // initial_guess_solution, solution_update

  size = (size+63)/64*64; // make multiple of typical cache line size
  size += 64; // align to typical cache line size

  return size;
}


void CGMRES_INITIALIZER_CREATE(struct CGMRES_INITIALIZER *initializer, 
                               REAL finite_difference_increment, int kmax,
                               void *memory) {
  int dimu = initializer->newton->ocp->model->dimu;
  int dimc = initializer->newton->ocp->model->dimc;
  int dimuc = dimu + dimc;

  INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_CREATE(initializer->newton, 
                                             finite_difference_increment, 
                                             memory);
  MFGMRES_CREATE(initializer->mfgmres, dimuc, kmax, memory);

  // zero memory (to avoid corrupted memory like e.g. NaN)
  int memsize =  CGMRES_INITIALIZER_MEMSIZE(initializer->newton->ocp);
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
  initializer->initial_guess_solution = sv_ptr;
  sv_ptr += 1;
  initializer->solution_update = sv_ptr;
  sv_ptr += 1;

  // align to typical cache line size
  size_t s_ptr = (size_t) sv_ptr;
  s_ptr = (s_ptr+63)/64*64;

  // void stuff
  char *c_ptr = (char *) s_ptr;

  int dim_solution = initializer->newton->ocp->dim_solution;
  CREATE_STRVEC(dim_solution, initializer->initial_guess_solution, c_ptr);
  c_ptr += initializer->initial_guess_solution->memsize;
  REAL initial_solution_scalar = 0.001;
  VECCSE(dim_solution, initial_solution_scalar, 
         initializer->initial_guess_solution, 0);

  CREATE_STRVEC(dim_solution, initializer->solution_update, c_ptr);
  c_ptr += initializer->solution_update->memsize;
  VECCSE(dim_solution, 0.0, initializer->solution_update, 0);

  initializer->newton_residual_tolerance = 1.0e-08;
  initializer->max_newton_iteration = 50;

  initializer->dim_solution = dimuc;
  initializer->memsize = CGMRES_INITIALIZER_MEMSIZE(initializer, kmax);

#if defined(RUNTIME_CHECKS)
  if(c_ptr > ((char *) mem) + mfgmres->memsize) {
    printf("\nerror: CGMRES_INITIALIZER_CREATE: outside memory bounds!\n\n");
    exit(1);
  }
#endif
}


void CGMRES_INITIALIZER_SET_TERMINATION_CRITERIONS(
    struct CGMRES_INITIALIZER *initializer, REAL newton_residual_tolerance, 
    int max_newton_iteration) {
  initializer->newton_residual_tolerance = newton_residual_tolerance;
  initializer->max_newton_iteration = max_newton_iteration;
}


void CGMRES_INITIALIZER_SET_INITIAL_GUESS_SOLUTION(
    struct CGMRES_INITIALIZER *initializer, 
    struct STRVEC *initial_guess_solution) {
  VECCP(initializer->dim_solution, initial_guess_solution, 0, 
        initializer->initial_guess_solution, 0);
}


void CGMRES_INITIALIZER_COMPUTE_INITIAL_SOLUTION(
    struct CGMRES_INITIALIZER *initializer, REAL initial_time, 
    struct STRVEC *initial_state, struct STRVEC *initial_solution) {
  VECCP(initializer->dim_solution, initializer->initial_guess_solution, 0, 
        initial_solution, 0);
  initializer->mfgmres_args->current_time = initial_time;
  initializer->mfgmres_args->current_state_ptr = initial_state;
  initializer->mfgmres_args->current_solution_ptr = initial_solution;

  int num_itr = 0;
  REAL optimality_error = INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_GET_ERROR_NORM(
      initializer->newton, initial_time, initial_state, initial_solution);

  while (optimality_error > initializer->newton_residual_tolerance 
         && num_itr < initializer->max_newton_iteration) {
    MFGMRES_SOLVE_LINEAR_PROBLEM(initializer->mfgmres, initializer->newton, 
                                 initializer->mfgmres_args, 
                                 initializer->solution_update);
    VECAD(initializer->dim_solution, 1.0, initializer->solution_update, 0,
          initial_solution, 0)
    optimality_error = INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_GET_ERROR_NORM(
        initializer->newton, initial_time, initial_state, initial_solution);
    ++num_itr;
  }
}


void CGMRES_INITIALIZER_GET_TERMINAL_COST_DERIVATIVE(
    struct CGMRES_INITIALIZER *initializer, REAL initial_time, 
    struct STRVEC *initial_state, struct STRVEC *terminal_cost_derivative) [
  INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_GET_TERMINAL_COST_DERIVATIVE(
    initializer->newton, initial_time, initial_state, terminal_cost_derivative);
]