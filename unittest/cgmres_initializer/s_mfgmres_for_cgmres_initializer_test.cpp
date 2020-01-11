#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <gtest/gtest.h>

extern "C" {
#include "s_memory_manager.h"
#include "s_linear_algebra.h"
#include "s_inexact_newton_for_zero_horizon_ocp.h"
#include "s_inexact_newton_for_zero_horizon_ocp_mfgmres_args.h"
#include "s_mfgmres_for_cgmres_initializer.h"
}


class s_mfgmres_for_cgmres_initializer_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    finite_diff = 1.0e-06;
    dimx = s_inexact_newton_for_zero_horizon_ocp_dimx();
    dimu = s_inexact_newton_for_zero_horizon_ocp_dimu();
    dimc = s_inexact_newton_for_zero_horizon_ocp_dimc();
    dim_solution = dimu + dimc;
    kmax = dim_solution;
    s_mfgmres_for_cgmres_initializer_create(&mfgmres, dim_solution, kmax);
    s_inexact_newton_for_zero_horizon_ocp_create(&newton, finite_diff);
    state = allocate_svec(dimx);
    solution = allocate_svec(dim_solution);
    update = allocate_svec(dim_solution);
    current_time = (float)rand()/RAND_MAX;
    for (int i=0; i<dimx; ++i) {
      state[i] = 0.0; 
    }
    for (int i=0; i<dim_solution; ++i) {
      solution[i] = 0.0;
    }
  }

  virtual void TearDown() {
    s_mfgmres_for_cgmres_initializer_delete(&mfgmres);
    s_inexact_newton_for_zero_horizon_ocp_delete(&newton);
    free_svec(state);
    free_svec(solution);
    free_svec(update);
  }

  struct s_mfgmres_for_cgmres_initializer mfgmres;
  struct s_inexact_newton_for_zero_horizon_ocp newton;
  struct s_inexact_newton_for_zero_horizon_ocp_mfgmres_args args;
  int dimx, dimu, dimc, dim_solution, kmax;
  float finite_diff, current_time;
  float *state, *solution, *update;
};


TEST_F(s_mfgmres_for_cgmres_initializer_test, memsize) {
  int expectes_size = 0;
  expectes_size += (kmax+1)*sizeof(float*); 
  expectes_size += (kmax+1)*(kmax+1)*sizeof(float); 
  expectes_size += (kmax+1)*sizeof(float*); 
  expectes_size += (kmax+1)*(dim_solution)*sizeof(float); 
  expectes_size += dim_solution*sizeof(float); 
  expectes_size += 3*(kmax+1)*sizeof(float); 
  expectes_size = (expectes_size+63)/64*64; 
  expectes_size += 64; 
  EXPECT_EQ(expectes_size, s_mfgmres_for_cgmres_initializer_memsize(dim_solution, kmax));
}


TEST_F(s_mfgmres_for_cgmres_initializer_test, solve_linear_problem) {
  args.current_time = current_time;
  args.current_state_ptr = state;
  args.current_solution_ptr = solution;
  int imax = 100;
  for (int i=0; i<imax; ++i) {
    std::cout << s_inexact_newton_for_zero_horizon_ocp_get_error_norm(&newton, 
                                                             current_time, 
                                                             state, solution) << std::endl;
    s_mfgmres_for_cgmres_initializer_solve_linear_problem(&mfgmres, &newton, 
                                                          &args, update);
    hpcgmres_svecadd(dim_solution, update, solution);
  }
  float delta = 1.0e-05;
  float err_norm 
      = s_inexact_newton_for_zero_horizon_ocp_get_error_norm(&newton, 
                                                             current_time, 
                                                             state, solution);
  EXPECT_TRUE(err_norm < delta);
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}