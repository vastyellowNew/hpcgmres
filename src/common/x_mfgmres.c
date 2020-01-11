int MFGMRES_STRSIZE() {
  return sizeof(struct MFGMRES);
}


int MFGMRES_MEMSIZE(int dim_linear_problem, int kmax) {
  int size = 0;

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


void MFGMRES_CREATE(struct MFGMRES *mfgmres, int dim_linear_problem, int kmax) {
  mfgmres->hessenberg_mat = ALLOCATE_MAT(kmax+1, kmax+1);
  mfgmres->basis_mat = ALLOCATE_MAT(kmax+1, dim_linear_problem);
  mfgmres->b_vec = ALLOCATE_VEC(dim_linear_problem);
  mfgmres->givens_c_vec = ALLOCATE_VEC(kmax+1);
  mfgmres->givens_s_vec = ALLOCATE_VEC(kmax+1);
  mfgmres->g_vec = ALLOCATE_VEC(kmax+1);
  mfgmres->dim_linear_problem = dim_linear_problem;
  if (kmax > dim_linear_problem) {
    mfgmres->kmax = dim_linear_problem;
  } 
  else {
    mfgmres->kmax = kmax;
  }
  mfgmres->memsize = MFGMRES_MEMSIZE(dim_linear_problem, kmax);
}


void MFGMRES_DELETE(struct MFGMRES *mfgmres) {
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
  mfgmres->g_vec[0] = sqrt(VECNRM2(dim_linear_problem, mfgmres->b_vec));
  VECMCP(dim_linear_problem, 1/mfgmres->g_vec[0], mfgmres->b_vec, 
         mfgmres->basis_mat[0]);

  // kk : the dimension of the Krylov subspace at the current iteration.
  int kk, jj;
  REAL nu;
  for (kk=0; kk<kmax; ++kk) {
    LINEAR_PROBLEM_COMPUTE_AX(linear_problem, linear_problem_args, 
                              mfgmres->basis_mat[kk], mfgmres->basis_mat[kk+1]);
    for (jj=0; jj<=kk; ++jj) {
      mfgmres->hessenberg_mat[kk][jj] = VECDOT(dim_linear_problem,  
                                               mfgmres->basis_mat[kk+1], 
                                               mfgmres->basis_mat[jj]);
      VECMAD(dim_linear_problem, -mfgmres->hessenberg_mat[kk][jj], 
             mfgmres->basis_mat[jj], mfgmres->basis_mat[kk+1]);
    }
    mfgmres->hessenberg_mat[kk][kk+1] = sqrt(VECNRM2(dim_linear_problem, 
                                                     mfgmres->basis_mat[kk+1]));
    if (fabs(mfgmres->hessenberg_mat[kk][kk+1]) < MACHINE_EPSILON) {
#if defined(RUNTIME_CHECKS)
      pinrtf("The modified Gram-Schmidt breakdown at k = %d.\n", k);
#endif
      break;
    }
    VECMUL(dim_linear_problem, 1/mfgmres->hessenberg_mat[kk][kk+1], 
           mfgmres->basis_mat[kk+1]);
    // Givens Rotation for QR factrization of the least squares problem.
    for (jj=0; jj<kk; ++jj) {
      APPLY_GIVENS_ROTATION(mfgmres->hessenberg_mat[kk], 
                            mfgmres->givens_c_vec, mfgmres->givens_s_vec, jj);
    }
    nu = sqrt(mfgmres->hessenberg_mat[kk][kk]*mfgmres->hessenberg_mat[kk][kk]
              +mfgmres->hessenberg_mat[kk][kk+1]*mfgmres->hessenberg_mat[kk][kk+1]);

    mfgmres->givens_c_vec[kk] = mfgmres->hessenberg_mat[kk][kk]/nu;
    mfgmres->givens_s_vec[kk] = - mfgmres->hessenberg_mat[kk][kk+1]/nu;
    mfgmres->hessenberg_mat[kk][kk] = mfgmres->givens_c_vec[kk] 
                                          * mfgmres->hessenberg_mat[kk][kk] 
                                      - mfgmres->givens_s_vec[kk] 
                                          * mfgmres->hessenberg_mat[kk][kk+1];
    mfgmres->hessenberg_mat[kk][kk+1] = 0;
    APPLY_GIVENS_ROTATION(mfgmres->g_vec, mfgmres->givens_c_vec, 
                          mfgmres->givens_s_vec, kk);
  }
  // Computes solution by solving hessenberg_mat * y = g_vec.
  int ii;
  REAL tmp;
  for (ii=kk-1; ii>=0; --ii) {
    tmp = mfgmres->g_vec[ii];
    for (jj=ii+1; jj<kk; ++jj) {
      tmp -= mfgmres->hessenberg_mat[jj][ii] * mfgmres->givens_c_vec[jj];
    }
    mfgmres->givens_c_vec[ii] = tmp / mfgmres->hessenberg_mat[ii][ii];
  }
  for (ii=0; ii<dim_linear_problem; ++ii) {
    tmp = 0;
    for (jj=0; jj<kk; ++jj) {
      tmp += mfgmres->basis_mat[jj][ii] * mfgmres->givens_c_vec[jj];
    }
    solution[ii] += tmp;
  }
}