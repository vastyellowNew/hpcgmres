int CGMRES_INITIALIZER_STRSIZE() {
  return sizeof(struct CGMRES_INITIALIZER);
}


int CGMRES_INITIALIZER_MEMSIZE(int kmax) { 
  int dimu = INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_DIMU();
  int dimc = INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_DIMC();
  int dimuc = dimu + dimc;

  int size = 0;
  size += MFGMRES_FOR_CGMRES_INITIALIZER_MEMSIZE(dimuc, kmax);
  size += INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_MEMSIZE();
  size += 2*dimuc*sizeof(REAL); // initial_guess_solution, solution_update;
  size = (size+63)/64*64; // make multiple of typical cache line size
  size += 64; // align to typical cache line size

  return size;
}


void CGMRES_INITIALIZER_CREATE(struct CGMRES_INITIALIZER *initializer, 
                               REAL finite_difference_increment, int kmax) {
  int dimx = INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_DIMX();
  int dimu = INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_DIMU();
  int dimc = INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_DIMC();
  int dimuc = dimu + dimc;
  MFGMRES_FOR_CGMRES_INITIALIZER_CREATE(&initializer->mfgmres, dimuc, kmax);
  INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_CREATE(&initializer->newton, 
                                             finite_difference_increment);

  initializer->initial_guess_solution = ALLOCATE_VEC(dimuc);
  initializer->solution_update = ALLOCATE_VEC(dimuc);

  initializer->newton_residual_tolerance = 1.0e-08;

  initializer->max_newton_iteration = 50;
  initializer->dim_solution = dimuc;
  initializer->memsize = CGMRES_INITIALIZER_MEMSIZE(kmax);
}


void CGMRES_INITIALIZER_DELETE(struct CGMRES_INITIALIZER *initializer) {
  MFGMRES_FOR_CGMRES_INITIALIZER_DELETE(&initializer->mfgmres);
  INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_DELETE(&initializer->newton);
  FREE_VEC(initializer->initial_guess_solution);
  FREE_VEC(initializer->solution_update);
}


void CGMRES_INITIALIZER_SET_TERMINATION_CRITERIONS(
    struct CGMRES_INITIALIZER *initializer, REAL newton_residual_tolerance, 
    int max_newton_iteration) {
  initializer->newton_residual_tolerance = newton_residual_tolerance;
  initializer->max_newton_iteration = max_newton_iteration;
}


void CGMRES_INITIALIZER_SET_INITIAL_GUESS_SOLUTION(
    struct CGMRES_INITIALIZER *initializer, REAL *initial_guess_solution) {
  VECCP(initializer->dim_solution, initial_guess_solution,
        initializer->initial_guess_solution);
}


void CGMRES_INITIALIZER_COMPUTE_INITIAL_SOLUTION(
    struct CGMRES_INITIALIZER *initializer, REAL initial_time, 
    REAL *initial_state, REAL *initial_solution) {
  VECCP(initializer->dim_solution, initializer->initial_guess_solution,  
        initial_solution);
  initializer->mfgmres_args.current_time = initial_time;
  initializer->mfgmres_args.current_state_ptr = initial_state;
  initializer->mfgmres_args.current_solution_ptr = initial_solution;

  int num_itr = 0;
  REAL optimality_error = INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_GET_ERROR_NORM(
      &initializer->newton, initial_time, initial_state, initial_solution);

  while (optimality_error > initializer->newton_residual_tolerance 
         && num_itr < initializer->max_newton_iteration) {
    MFGMRES_FOR_CGMRES_INITIALIZER_SOLVE_LINEAR_PROBLEM(
        &initializer->mfgmres, &initializer->newton, &initializer->mfgmres_args, 
        initializer->solution_update);
    VECAD(initializer->dim_solution, initializer->solution_update,
          initial_solution);
    optimality_error = INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_GET_ERROR_NORM(
        &initializer->newton, initial_time, initial_state, initial_solution);
    ++num_itr;
  }
}


void CGMRES_INITIALIZER_GET_TERMINAL_COST_DERIVATIVE(
    struct CGMRES_INITIALIZER *initializer, REAL initial_time, 
    REAL *initial_state, REAL *terminal_cost_derivative) {
  INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_GET_TERMINAL_COST_DERIVATIVE(
    &initializer->newton, initial_time, initial_state, 
    terminal_cost_derivative);
} 