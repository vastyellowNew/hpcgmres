int INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_STRSIZE() {
  return sizeof(struct INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP);
}


int INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_MEMSIZE(
    struct INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP *newton) {
  int size = 0;

  size += 1*ZERO_HORIZON_OCP_MEMSIZE(newton->ocp);

  // size of pointers and int
  size += 3*sizeof(struct STRVEC); // hessenberg_mat

  int dimu = newton->ocp->model->dimu;
  int dimc = newton->ocp->model->dimc;

  // size of matrices and vectors
  size += 1*SIZE_STRVEC(dimu+dimc); // incremented_solution
  size += 1*SIZE_STRVEC(dimu+dimc); // optimality_residual
  size += 1*SIZE_STRVEC(dimu+dimc); // optimality_residual1

  size = (size+63)/64*64; // make multiple of typical cache line size
  size += 64; // align to typical cache line size

  return size;
}


void INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_CREATE(
    struct INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP *newton, 
    REAL finite_difference_increment, void *memory) {
  ZERO_HORIZON_OCP_CREATE(newton->ocp, memory);

  // zero memory (to avoid corrupted memory like e.g. NaN)
  int memsize =  INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_MEMSIZE(newton);
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
  newton->incremented_solution = sv_ptr;
  sv_ptr += 1;
  newton->optimality_residual = sv_ptr;
  sv_ptr += 1;
  newton->optimality_residual1 = sv_ptr;
  sv_ptr += 1;

  // align to typical cache line size
  size_t s_ptr = (size_t) sv_ptr;
  s_ptr = (s_ptr+63)/64*64;

  // void stuff
  char *c_ptr = (char *) s_ptr;

  CREATE_STRVEC(newton->ocp->dim_solution, newton->incremented_solution, c_ptr);
  c_ptr += newton->incremented_solution->memsize;
  VECCSE(newton->ocp->dim_solution, 0.0, newton->incremented_solution, 0);

  CREATE_STRVEC(newton->ocp->dim_solution, newton->optimality_residual, c_ptr);
  c_ptr += newton->optimality_residual->memsize;
  VECCSE(newton->ocp->dim_solution, 0.0, newton->optimality_residual, 0);

  CREATE_STRVEC(newton->ocp->dim_solution, newton->optimality_residual1, c_ptr);
  c_ptr += newton->optimality_residual1->memsize;
  VECCSE(newton->ocp->dim_solution, 0.0, newton->optimality_residual1, 0);

  newton->finite_difference_increment = finite_difference_increment;
  newton->dim_solution = newton->ocp->dim_solution;
  newton->memsize = INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_MEMSIZE(newton->ocp);

#if defined(RUNTIME_CHECKS)
  if(c_ptr > ((char *) mem) + mfgmres->memsize) {
    printf("\nerror: INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_CREATE: outside memory 
            bounds!\n\n");
    exit(1);
  }
#endif
}


void INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_COMPUTE_B(
    struct INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP *newton, 
    struct INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_MFGMRES_ARGS *args,
    struct STRVEC *direction, struct STRVEC *b) {
  AXPY(newton->ocp->dim_solution, newton->finite_difference_increment, 
       direction, args->current_solution_ptr, 0, newton->incremented_solution, 
       0);
  ZERO_HORIZON_OCP_COMPUTE_OPTIMALITY_RESIDUAL(newton->ocp, args->current_time, 
                                               args->current_state_ptr, 
                                               args->current_solution_ptr, 
                                               newton->optimality_residual);
  ZERO_HORIZON_OCP_COMPUTE_OPTIMALITY_RESIDUAL(newton->ocp, args->current_time, 
                                               args->current_state_ptr, 
                                               newton->incremented_solution, 
                                               newton->optimality_residual1);
  AXBPY(newton->ocp->dim_solution, 1/newton->finite_difference_increment-1, 
        newton->optimality_residual, 1/newton->finite_difference_increment, 
        newton->optimality_residual1, 0, b, 0); 
}


void INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_COMPUTE_AX(
    struct INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP *newton, 
    struct INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_MFGMRES_ARGS *args,
    struct STRVEC *direction, struct STRVEC *ax) {
  AXPY(newton->ocp->dim_solution, newton->finite_difference_increment, 
       direction, args->solution_ptr, 0, newton->incremented_solution, 0);
  ZERO_HORIZON_OCP_COMPUTE_OPTIMALITY_RESIDUAL(newton->ocp, args->time_param, 
                                               args->state_ptr, 
                                               newton->incremented_solution, 
                                               newton->optimality_residual1);
  AXBPY(newton->ocp->dim_solution, -1/newton->finite_difference_increment, 
        newton->optimality_residual, 1/newton->finite_difference_increment, 
        newton->optimality_residual1, 0, ax, 0); 
}


REAL INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_GET_ERROR_NORM(
    struct INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP *newton, REAL initial_time, 
    struct STRVEC *initial_state, struct STRVEC *initial_solution) {
  ZERO_HORIZON_OCP_COMPUTE_OPTIMALITY_RESIDUAL(newton->ocp, initial_time, 
                                               initial_state, initial_solution, 
                                               newton->optimality_residual);
  return sqrt(VECDOT(newton->ocp->dim_solution, newton->optimality_residual, 0,
                     newton->optimality_residual, 0));
}


void INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_GET_TERMINAL_COST_DERIVATIVE(
    struct INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP *newton, REAL initial_time, 
    struct STRVEC *initial_state, struct STRVEC *terminal_cost_derivative) {
  ZERO_HORIZON_OCP_COMPUTE_TERMINAL_COST_DERIVATIVE(newton->ocp, initial_time, 
                                                    initial_state, 
                                                    terminal_cost_derivative);
}