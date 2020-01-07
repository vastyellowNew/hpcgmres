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
  size += 3*(kmax+1) * sizeof(REAL); // givens_c_vec, givens_s_vec, g_vec

  size = (size+63)/64*64; // make multiple of typical cache line size
  size += 64; // align to typical cache line size

  return size;
}


void MFGMRES_CREATE(struct MFGMRES *mfgmres, int dim_linear_problem, 
                    int kmax) {
  ALLOCATE_MAT(kmax+1, kmax+1, mfgmres->hessenberg_mat);
  ALLOCATE_MAT(kmax+1, dim_linear_problem, mfgmres->basis_mat);
  ALLOCATE_VEC(dim_linear_problem, mfgmres->b_vec);
  ALLOCATE_VEC(kmax+1, mfgmres->givens_c_vec);
  ALLOCATE_VEC(kmax+1, mfgmres->givens_s_vec);
  ALLOCATE_VEC(kmax+1, mfgmres->g_vec);
}


void MFGMRES_DELETE(struct MFGMRES *mfgmres, int dim_linear_problem, 
                    int kmax) {
  FREE_MAT()
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
  REAL **hessenberg_mat = mfgmres->hessenberg_mat;
  REAL **basis_vec = mfgmres->basis_vec;
  REAL *givens_c_vec = mfgmres->givens_c_vec;
  REAL *givens_s_vec = mfgmres->givens_s_vec;
  REAL *g_vec = mfgmres->g_vec;

  VECCSE(kmax+1, 0., givens_c_vec, 0);
  VECCSE(kmax+1, 0., givens_s_vec, 0);
  VECCSE(kmax+1, 0., g_vec, 0);
  // Generates the initial basis of the Krylov subspace.
  LINEAR_PROBLEM_COMPUTE_B(linear_problem, linear_problem_args, solution, b_vec);
  // g_vec[0] = sqrt(b_vec^2)
  VECIN1(sqrt(VECDOT(dim_linear_problem, b_vec, 0, b_vec, 0)), g_vec, 0);
  // basis_mat[0] = b_vec / g_vec[0]
  VECCPSC(dim_linear_problem, 1/VECEX1(g_vec, 0), b_vec, 0, basis_vec, 0);

  // k : the dimension of the Krylov subspace at the current iteration.
  int k;
  for (k=0; k<kmax_; ++k) {
    LINEAR_PROBLEM_COMPUTE_AX(linear_problem, linear_problem_args, basis_vec+k, 
                              basis_vec+k+1);
    for (int j=0; j<=k; ++j) {
      MATIN1(VECDOT(dim_linear_problem, basis_vec+(k+1), 
                    (k+1)*dim_linear_problem, basis_vec, j*dim_linear_problem), 
             hessenberg_mat, k, j);
      VECCAD(dim_linear_problem, -MATEX1(hessenberg_mat, k, j), basis_vec, 
             j*dim_linear_problem, basis_vec, (k+1)*dim_linear_problem);
      MATIN1(sqrt(VECDOT(dim_linear_problem, 
                  basis_vec, (k+1)*dim_linear_problem, 
                  basis_vec, (k+1)*dim_linear_problem)), 
              hessenberg_mat, k, k+1);
    if (abs(MATEX1(hessenberg_mat, k, k+1)) < MACHINE_EPSILON)) {
#if defined(RUNTIME_CHECKS)
      pinrtf("The modified Gram-Schmidt breakdown at k = %d.\n", k);
#endif
      break;
    }
    // basis_mat[k+1] = basis_mat[k+1] / hessenberg_mat[k][k+1];
    VECSC(dim_linear_problem, 1/MATEX1(hessenberg_mat, k, k+1), basis_vec, 
          (k+1)*dim_linear_problem);
    // Givens Rotation for QR factrization of the least squares problem.
    for (int j=0; j<k; ++j) {
      MFGMRES_GIVENS_ROTATION_MAT(hessenberg_mat, givens_c_vec, givens_s_vec, 
                                  k, j);
    }
    REAL nu 
        = sqrt(MATEX1(hessenberg_mat, k, k)*MATEX1(hessenberg_mat, k, k)
               +MATEX1(hessenberg_mat, k, k+1)*MATEX1(hessenberg_matm k, k+1));
      VECIN1(MATEX1(hessenberg_mat, k, k)/nu, givens_c_vec, k);
      VECIN1(-MATEX1(hessenberg_mat, k, k+1)/nu, givens_s_vec, k);
      MATIN1(VECEX1(givens_c_vec, k)*MATEX1(hessenberg_mat, k, k)
              -VECEX1(givens_s_vec, k)*MATEX1(hessenberg_mat, k, k+1),
              hessenberg_mat, k, k);
      MATN1(0, hessenberg_mat, k, k+1);
      MFGMRES_GIVENS_ROTATION_VEC(g_vec, givens_c_vec, givens_s_vec, k);
  }

  // Computes solution_vec by solving hessenberg_mat * y = g_vec.
  TRSV_UNN(k, hessenberg_mat, 0, 0, g_vec, 0, solution_vec, 0);
}