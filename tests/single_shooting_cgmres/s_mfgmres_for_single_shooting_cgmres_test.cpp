#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <gtest/gtest.h>

extern "C" {
#include "s_memory_manager.h"
#include "s_linear_algebra.h"
#include "s_single_shooting_continuation.h"
#include "s_single_shooting_continuation_mfgmres_args.h"
#include "s_mfgmres_for_single_shooting_cgmres.h"
}


class s_mfgmres_for_single_shooting_cgmres_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand(time(NULL));
    finite_diff = 1.0e-03;
    N = 50;
    dimx = 4;
    dimu = 2;
    dimc = 1;
    dim_solution = N*(dimu+dimc);
    kmax = dim_solution;
    s_mfgmres_for_single_shooting_cgmres_create(&mfgmres, dim_solution, kmax);
  }

  virtual void TearDown() {
    s_mfgmres_for_single_shooting_cgmres_delete(&mfgmres);
  }

  struct s_mfgmres_for_single_shooting_cgmres mfgmres;
  int dimx, dimu, dimc, N, dim_solution, kmax;
  float finite_diff, current_time;
  float *state, *solution, *update;
};


TEST_F(s_mfgmres_for_single_shooting_cgmres_test, memsize) {
  int expectes_size = 0;
  expectes_size += (kmax+1)*sizeof(float*); 
  expectes_size += (kmax+1)*(kmax+1)*sizeof(float); 
  expectes_size += (kmax+1)*sizeof(float*); 
  expectes_size += (kmax+1)*(dim_solution)*sizeof(float); 
  expectes_size += dim_solution*sizeof(float); 
  expectes_size += 3*(kmax+1)*sizeof(float); 
  expectes_size = (expectes_size+63)/64*64; 
  expectes_size += 64; 
  EXPECT_EQ(expectes_size, s_mfgmres_for_single_shooting_cgmres_memsize(dim_solution, kmax));
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}