#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <gtest/gtest.h>

extern "C" {
#include "d_memory_manager.h"
#include "d_linear_algebra.h"
#include "d_inexact_newton_for_zero_horizon_ocp.h"
#include "d_inexact_newton_for_zero_horizon_ocp_mfgmres_args.h"
#include "d_mfgmres_for_cgmres_initializer.h"
}


class d_mfgmres_for_cgmres_initializer_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand(time(NULL));
    finite_diff = 1.0e-08;
    dimx = d_inexact_newton_for_zero_horizon_ocp_dimx();
    dimu = d_inexact_newton_for_zero_horizon_ocp_dimu();
    dimc = d_inexact_newton_for_zero_horizon_ocp_dimc();
    dim_solution = dimu + dimc;
    kmax = dim_solution;
    d_mfgmres_for_cgmres_initializer_create(&mfgmres, dim_solution, kmax);
    d_inexact_newton_for_zero_horizon_ocp_create(&newton, finite_diff);
    state = allocate_dvec(dimx);
    solution = allocate_dvec(dim_solution);
    update = allocate_dvec(dim_solution);
    current_time = (double)rand()/RAND_MAX;
    for (int i=0; i<dimx; ++i) {
      state[i] = (double)rand()/RAND_MAX;
    }
    for (int i=0; i<dim_solution; ++i) {
      solution[i] = (double)rand()/RAND_MAX;
    }
  }

  virtual void TearDown() {
    d_mfgmres_for_cgmres_initializer_delete(&mfgmres);
    d_inexact_newton_for_zero_horizon_ocp_delete(&newton);
    free_dvec(state);
    free_dvec(solution);
    free_dvec(update);
  }

  struct d_mfgmres_for_cgmres_initializer mfgmres;
  struct d_inexact_newton_for_zero_horizon_ocp newton;
  struct d_inexact_newton_for_zero_horizon_ocp_mfgmres_args args;
  int dimx, dimu, dimc, dim_solution, kmax;
  double finite_diff, current_time;
  double *state, *solution, *update;
};


TEST_F(d_mfgmres_for_cgmres_initializer_test, memsize) {
  int expected_size = 0;
  expected_size += (kmax+1)*sizeof(double*); 
  expected_size += (kmax+1)*(kmax+1)*sizeof(double); 
  expected_size += (kmax+1)*sizeof(double*); 
  expected_size += (kmax+1)*(dim_solution)*sizeof(double); 
  expected_size += dim_solution*sizeof(double); 
  expected_size += 3*(kmax+1)*sizeof(double); 
  expected_size = (expected_size+63)/64*64; 
  expected_size += 64; 
  EXPECT_EQ(expected_size, d_mfgmres_for_cgmres_initializer_memsize(dim_solution, kmax));
}


TEST_F(d_mfgmres_for_cgmres_initializer_test, solve_linear_problem) {
  args.current_time = current_time;
  args.current_state_ptr = state;
  args.current_solution_ptr = solution;
  int imax = 100;
  for (int i=0; i<imax; ++i) {
    d_mfgmres_for_cgmres_initializer_solve_linear_problem(&mfgmres, &newton, 
                                                          &args, update);
    hpcgmres_dvecadd(dim_solution, update, solution);
  }
  double delta = 1.0e-05;
  double err_norm 
      = d_inexact_newton_for_zero_horizon_ocp_get_error_norm(&newton, 
                                                             current_time, 
                                                             state, solution);
  EXPECT_TRUE(err_norm < delta);
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}