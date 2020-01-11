int INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_STRSIZE() {
  return sizeof(struct INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP);
}


int INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_MEMSIZE() {
  int dimuc = ZERO_HORIZON_OCP_DIMU() + ZERO_HORIZON_OCP_DIMC();
  int size = 0;
  size += 1*ZERO_HORIZON_OCP_MEMSIZE();
  size += 3*dimuc*sizeof(REAL); // incremented_solution, optimality_residual, optimality_residual1;
  size = (size+63)/64*64; // make multiple of typical cache line size
  size += 64; // align to typical cache line size
  return size;
}


void INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_CREATE(
    struct INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP *newton, 
    REAL finite_difference_increment) {
  ZERO_HORIZON_OCP_CREATE(&newton->ocp);
  int dimuc = ZERO_HORIZON_OCP_DIMU() + ZERO_HORIZON_OCP_DIMC();
  newton->incremented_solution = ALLOCATE_VEC(dimuc);
  newton->optimality_residual = ALLOCATE_VEC(dimuc);
  newton->optimality_residual1 = ALLOCATE_VEC(dimuc);
  newton->finite_difference_increment = finite_difference_increment; 
  newton->dim_solution = dimuc;
  newton->memsize = INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_MEMSIZE();
}


void INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_DELETE(
    struct INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP *newton) {
  FREE_VEC(newton->incremented_solution);
  FREE_VEC(newton->optimality_residual);
  FREE_VEC(newton->optimality_residual1);
  ZERO_HORIZON_OCP_DELETE(&newton->ocp);
}


void INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_COMPUTE_B(
    struct INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP *newton, 
    struct INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_MFGMRES_ARGS *args,
    REAL *direction, REAL *b) {
  AXPY(newton->dim_solution, newton->finite_difference_increment, 
       direction, args->current_solution_ptr, newton->incremented_solution);
  ZERO_HORIZON_OCP_COMPUTE_OPTIMALITY_RESIDUAL(&newton->ocp, args->current_time, 
                                               args->current_state_ptr, 
                                               args->current_solution_ptr, 
                                               newton->optimality_residual);
  ZERO_HORIZON_OCP_COMPUTE_OPTIMALITY_RESIDUAL(&newton->ocp, args->current_time, 
                                               args->current_state_ptr, 
                                               newton->incremented_solution, 
                                               newton->optimality_residual1);
  AXPBY(newton->dim_solution, 1/newton->finite_difference_increment-1, 
        newton->optimality_residual, -1/newton->finite_difference_increment, 
        newton->optimality_residual1, b); 
}


void INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_COMPUTE_AX(
    struct INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP *newton, 
    struct INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_MFGMRES_ARGS *args,
    REAL *direction, REAL *ax) {
  AXPY(newton->dim_solution, newton->finite_difference_increment, 
       direction, args->current_solution_ptr, newton->incremented_solution);
  ZERO_HORIZON_OCP_COMPUTE_OPTIMALITY_RESIDUAL(&newton->ocp, args->current_time, 
                                               args->current_state_ptr, 
                                               newton->incremented_solution, 
                                               newton->optimality_residual1);
  AXPBY(newton->dim_solution, -1/newton->finite_difference_increment, 
        newton->optimality_residual, 1/newton->finite_difference_increment, 
        newton->optimality_residual1, ax); 
}


REAL INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_GET_ERROR_NORM(
    struct INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP *newton, REAL initial_time, 
    REAL *initial_state, REAL *initial_solution) {
  ZERO_HORIZON_OCP_COMPUTE_OPTIMALITY_RESIDUAL(&newton->ocp, initial_time, 
                                               initial_state, initial_solution, 
                                               newton->optimality_residual);
  return sqrt(VECNRM2(newton->dim_solution, newton->optimality_residual));
}


void INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_GET_TERMINAL_COST_DERIVATIVE(
    struct INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP *newton, REAL initial_time, 
    REAL *initial_state, REAL *terminal_cost_derivative) {
  ZERO_HORIZON_OCP_COMPUTE_TERMINAL_COST_DERIVATIVE(&newton->ocp, initial_time, 
                                                    initial_state, 
                                                    terminal_cost_derivative);
}


int INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_DIMX() {
  return ZERO_HORIZON_OCP_DIMX();
}


int INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_DIMU() {
  return ZERO_HORIZON_OCP_DIMU();
}


int INEXACT_NEWTON_FOR_ZERO_HORIZON_OCP_DIMC() {
  return ZERO_HORIZON_OCP_DIMC();
}