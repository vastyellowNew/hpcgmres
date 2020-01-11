int MFGMRES_STRSIZE() {
  return sizeof(struct MFGMRES);
}


int MFGMRES_MEMSIZE(int dim_linear_problem, int kmax) {
  int size = 0;

  // size of pointers and int
  size += (kmax+1)*sizeof(REAL*); // hessenberg_mat
  size += (kmax+1)*(kmax+1)*sizeof(REAL); // hessenberg_mat

  size += (kmax+1)*sizeof(REAL*); // basis_mat
  size += (kmax+1)*(dim_linear_problem)*sizeof(REAL); // basis_mat

  size += dim_linear_problem*sizeof(REAL); // b_vec 
  size += 3*(kmax+1)*sizeof(REAL); // givens_c_vec, givens_s_vec, g_vec

  size = (size+63)/64*64; // make multiple of typical cache line size
  size += 64; // align to typical cache line size

  return size;
}


void MFGMRES_CREATE(struct MFGMRES *mfgmres, int dim_linear_problem, 
                    int kmax) {
  // ALLOCATE_MAT(kmax+1, kmax+1, mfgmres->hessenberg_mat);
  // ALLOCATE_MAT(kmax+1, dim_linear_problem, mfgmres->basis_mat);
  // ALLOCATE_VEC(dim_linear_problem, mfgmres->b_vec);
  // ALLOCATE_VEC(kmax+1, mfgmres->givens_c_vec);
  // ALLOCATE_VEC(kmax+1, mfgmres->givens_s_vec);
  // ALLOCATE_VEC(kmax+1, mfgmres->g_vec);
  mfgmres->hessenberg_mat = ALLOCATE_MAT(kmax+1, kmax+1);
  mfgmres->basis_mat = ALLOCATE_MAT(kmax+1, dim_linear_problem);
  mfgmres->b_vec = ALLOCATE_VEC(dim_linear_problem);
  mfgmres->givens_c_vec = ALLOCATE_VEC(kmax+1);
  mfgmres->givens_s_vec = ALLOCATE_VEC(kmax+1);
  mfgmres->g_vec = ALLOCATE_VEC(kmax+1);
}


void MFGMRES_DELETE(struct MFGMRES *mfgmres, int dim_linear_problem, 
                    int kmax) {
  FREE_MAT(mfgmres->hessenberg_mat);
  FREE_MAT(mfgmres->basis_mat);
  FREE_VEC(mfgmres->b_vec);
  FREE_VEC(mfgmres->givens_c_vec);
  FREE_VEC(mfgmres->givens_s_vec);
  FREE_VEC(mfgmres->g_vec);
}


void MFGMRES_SOLVE_LINEAR_PROBLEM(
    struct MFGMRES *mfgmres, struct LINEAR_PROBLEM *linear_problem, 
    struct LINEAR_PROBLEM_ARGS *linear_problem_args, REAL *solution) {
  int dim_linear_problem = mfgmres->dim_linear_problem;
  int kmax = mfgmres->kmax;
  VECSET(kmax+1, 0., mfgmres->givens_c_vec);
  VECSET(kmax+1, 0., mfgmres->givens_s_vec);
  VECSET(kmax+1, 0., mfgmres->g_vec);
  // Generates the initial basis of the Krylov subspace.
  LINEAR_PROBLEM_COMPUTE_B(linear_problem, linear_problem_args, solution, 
                           mfgmres->b_vec);
  // g_vec[0] = sqrt(b_vec^2)
  mfgmres->g_vec[0] = sqrt(VECNRM2(dim_linear_problem, mfgmres->b_vec));
  // basis_mat[0] = b_vec / g_vec[0]
  VECMCP(dim_linear_problem, 1/mfgmres->g_vec[0], mfgmres->b_vec, basis_mat[0]);

  // k : the dimension of the Krylov subspace at the current iteration.
  int k;
  for (k=0; k<kmax; ++k) {
    LINEAR_PROBLEM_COMPUTE_AX(linear_problem, linear_problem_args, basis_mat[k], 
                              basis_mat[k+1]);
    for (int j=0; j<=k; ++j) {
      mfgmres->hessenberg_mat[k][j] = VECDOT(dim_linear_problem,  
                                             mfgmres->basis_mat[k+1], 
                                             mfgmres->basis_mat[j]);
      VECMAD(dim_linear_problem, -mfgmres->hessenberg_mat[k][j], 
             mfgmres->basis_mat[j], mfgmres->basis_mat[k+1]);
    }
    mfgmres->hessenberg_mat[k][k+1] = sqrt(VECNRM2(dim_linear_problem, 
                                                    mfgmres->basis_mat[k+1]));
    if (abs(mfgmres->hessenberg_mat[k][k+1])) < MACHINE_EPSILON)) {
#if defined(RUNTIME_CHECKS)
      pinrtf("The modified Gram-Schmidt breakdown at k = %d.\n", k);
#endif
      break;
    }
    // basis_mat[k+1] = basis_mat[k+1] / hessenberg_mat[k][k+1];
    VECMUL(dim_linear_problem, 1/mfgmres->hessenberg_mat[k][k+1], 
           mfgmres->basis_mat[k+1]);
    // Givens Rotation for QR factrization of the least squares problem.
    for (int j=0; j<k; ++j) {
      MFGMRES_GIVENS_ROTATION(mfgmres->hessenberg_mat[k], mfgmres->givens_c_vec, 
                              mfgmres->givens_s_vec, j);
    }
    REAL nu 
        = sqrt(mfgmres->hessenberg_mat[k][k]*mfgmres->hessenberg_mat[k][k]
               +mfgmres->hessenberg_mat[k][k+1]*mfgmres->hessenberg_mat[k][k+1]);

      mfgmres->givens_c_vec_[k] = mfgmres->hessenberg_mat_[k][k] / nu;
      mfgmres->givens_s_vec_[k] = - mfgmres->hessenberg_mat_[k][k+1] / nu;
      mfgmres->hessenberg_mat_[k][k] = mfgmres->givens_c_vec_[k] 
                                           * mfgmres->hessenberg_mat_[k][k] 
                                       - mfgmres->givens_s_vec_[k] 
                                           * mfgmres->hessenberg_mat_[k][k+1];
      MFGMRES_GIVENS_ROTATION(mfgmres->g_vec, mfgmres->givens_c_vec, 
                              mfgmres->givens_s_vec, k);
      mfgmres->hessenberg_mat_[k][k+1] = 0;
  }

  // Computes solution_vec by solving hessenberg_mat * y = g_vec.
  int ii;
  REAL tmp;
  for (ii=k-1: ii>=0; --ii) {
    tmp = mfgmres->g_vec[i];
    int j;
    for (j=ii+1; j<k; ++j) {
      tmp -= mfgmres->hessenberg_mat[j][ii] * mfgmres->givens_c_vec[j];
    }
    mfgmres->givens_c_vec[ii] = tmp / mfgmres->hessenberg_mat[ii][ii];
  }
  for (ii=0; ii<dim_linear_problem; ++ii) {
    tmp = 0;
    for (int j=0; j<k; ++j) {
      tmp += mfgmres->basis_mat[j][ii] * mfgmres->givens_c_vec[j];
    }
    mfgmres->solution_vec[ii] += tmp;
  }
}


    for (int i=k-1; i>=0; --i) {
      double tmp = g_vec_[i];
      for (int j=i+1; j<k; ++j) {
        tmp -= hessenberg_mat_[j][i] * givens_c_vec_[j];
      }
      givens_c_vec_[i] = tmp / hessenberg_mat_[i][i];
    }
    for (int i=0; i<dim_linear_problem_; ++i) {
      double tmp = 0;
      for (int j=0; j<k; ++j) { 
        tmp += basis_mat_[j][i] * givens_c_vec_[j];
      }
      solution_vec[i] += tmp;
    }
  }
