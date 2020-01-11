#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <gtest/gtest.h>

extern "C" {
#include "d_memory_manager.h"
#include "d_linear_algebra.h"
#include "d_single_shooting_continuation.h"
#include "d_single_shooting_continuation_mfgmres_args.h"
#include "d_mfgmres_for_single_shooting_cgmres.h"
}


class d_mfgmres_for_single_shooting_cgmres_test : public ::testing::Test {
protected:
  virtual void SetUp() {
    srand(time(NULL));
    finite_diff = 1.0e-08;
    N = 50;
    dimx = 4;
    dimu = 2;
    dimc = 1;
    dim_solution = N*(dimu+dimc);
    kmax = dim_solution;
    d_mfgmres_for_single_shooting_cgmres_create(&mfgmres, dim_solution, kmax);
  }

  virtual void TearDown() {
    d_mfgmres_for_single_shooting_cgmres_delete(&mfgmres);
  }

  struct d_mfgmres_for_single_shooting_cgmres mfgmres;
  int dimx, dimu, dimc, N, dim_solution, kmax;
  double finite_diff, current_time;
  double *state, *solution, *update;
};


TEST_F(d_mfgmres_for_single_shooting_cgmres_test, memsize) {
  int expected_size = 0;
  expected_size += (kmax+1)*sizeof(double*); 
  expected_size += (kmax+1)*(kmax+1)*sizeof(double); 
  expected_size += (kmax+1)*sizeof(double*); 
  expected_size += (kmax+1)*(dim_solution)*sizeof(double); 
  expected_size += dim_solution*sizeof(double); 
  expected_size += 3*(kmax+1)*sizeof(double); 
  expected_size = (expected_size+63)/64*64; 
  expected_size += 64; 
  EXPECT_EQ(expected_size, d_mfgmres_for_single_shooting_cgmres_memsize(dim_solution, kmax));
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}