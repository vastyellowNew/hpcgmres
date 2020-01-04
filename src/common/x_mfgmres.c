int MFGMRES_STRSIZE() {
  return sizeof(struct MFGMRES);
}


int MFGMRES_MEMSIZE(int dim_linear_problem, int kmax) {
  int size = 0;

  // size of pointers and int
  size += 1*sizeof(struct STRMAT); // hessenberg_mat
  size += 1*(kmax+1)*sizeof(struct STRVEC); // basis_vec
  size += 4*sizeof(struct STRVEC); // b_vec, givens_c_vec, givens_s_vec, g_vec

  // size of matrices and vectors
  size += 1*SIZE_STRMAT(kmax+1, kmax+1); // hessenberg_mat
  size += 1*SIZE_STRVEC((kmax+1)*dim_linear_problem); // basis_mat
  size += 1*SIZE_STRVEC(dim_linear_problem); // b_vec 
  size += 3*SIZE_STRVEC(kmax+1); // givens_c_vec, givens_s_vec, g_vec

  size = (size+63)/64*64; // make multiple of typical cache line size
  size += 64; // align to typical cache line size

  return size;
}


void MFGMRES_CREATE(struct MFGMRES *mfgmres, int dim_linear_problem, 
                    int kmax, void *memory) {
  // zero memory (to avoid corrupted memory like e.g. NaN)
  int memsize = MFGMRES_MEMSIZE(dim_linear_problem, kmax);
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

  // matrix struct stuff
  struct STRMAT *sm_ptr = (struct STRMAT *) memory;
  mfgmres->hessenberg_mat = sm_ptr;
  sm_ptr += 1;

  // vector struct stuff
  struct STRVEC *sv_ptr = (struct STRVEC *) sm_ptr;
  mfgmres->basis_vec = sv_ptr;
  sv_ptr += kmax+1;
  mfgmres->b_vec = sv_ptr;
  sv_ptr += 1;
  mfgmres->givens_c_vec = sv_ptr;
  sv_ptr += 1;
  mfgmres->givens_s_vec = sv_ptr;
  sv_ptr += 1;
  mfgmres->g_vec = sv_ptr;
  sv_ptr += 1;

  // align to typical cache line size
  size_t s_ptr = (size_t) sv_ptr;
  s_ptr = (s_ptr+63)/64*64;

  // void stuff
  char *c_ptr = (char *) s_ptr;
  char *tmp_ptr;

  CREATE_STRMAT(kmax+1, kmax+1, mfgmres->hessenberg_mat, c_ptr);
  c_ptr += mfgmres->hessenberg_mat->memsize;
  MATCSE(kmax+1, kmax+1, 0.0, mfgmres->hessenberg_mat, 0, 0);

  tmp_ptr = c_ptr;
  c_ptr += SIZE_STRVEC((kmax+1)*dim_linear_problem);
  for (ii=0; ii<kmax+1; ++ii) {
    CREATE_STRVEC(dim_linear_problem, mfgmres->basis_vec+ii, tmp_ptr);
    tmp_ptr += dim_linear_problem*sizeof(REAL);
    VECCSE(dim_linear_problem, 0.0, mfgmres->basis_vec+ii, 0);
  }

  CREATE_STRVEC(dim_linear_problem, mfgmres->b_vec, c_ptr);
  c_ptr += mfgmres->b_vec->memsize;
  VECCSE(dim_linear_problem, 0.0, mfgmres->b_vec, 0);

  CREATE_STRVEC(kmax+1, mfgmres->givens_c_vec, c_ptr);
  c_ptr += mfgmres->givens_c_vec->memsize;
  VECCSE(kmax+1, 0.0, mfgmres->givens_c_vec, 0);

  CREATE_STRVEC(kmax+1, mfgmres->givens_s_vec, c_ptr);
  c_ptr += mfgmres->givens_s_vec->memsize;
  VECCSE(kmax+1, 0.0, mfgmres->givens_s_vec, 0);

  CREATE_STRVEC(kmax+1, mfgmres->g_vec, c_ptr);
  c_ptr += mfgmres->g_vec->memsize;
  VECCSE(kmax+1, 0.0, mfgmres->g_vec, 0);

  mfgmres->dim_linear_problem = dim_linear_problem;
  mfgmres->kmax = kmax;
  mfgmres->memsize = MFGMRES_MEMSIZE(dim_linear_problem, kmax);

#if defined(RUNTIME_CHECKS)
  if(c_ptr > ((char *) mem) + mfgmres->memsize) {
    printf("\nerror: MFGMRES_CREATE: outside memory bounds!\n\n");
    exit(1);
  }
#endif
}


void MFGMRES_SOLVE_LINEAR_PROBLEM(
    struct MFGMRES *mfgmres, struct LINEAR_PROBLEM *linear_problem, 
    struct LINEAR_PROBLEM_ARGS *linear_problem_args,
    struct STRVEC *solution_vec) {

  int dim_linear_problem = mfgmres->dim->dim_linear_problem;
  int kmax = mfgmres->dim->kmax;
  struct STRMAT *hessenberg_mat = mfgmres->hessenberg_mat;
  struct STRVEC *basis_vec = mfgmres->basis_vec;
  struct STRVEC *givens_c_vec = mfgmres->givens_c_vec;
  struct STRVEC *givens_s_vec = mfgmres->givens_s_vec;
  struct STRVEC *g_vec = mfgmres->g_vec;

  VECCSE(kmax+1, 0., givens_c_vec, 0);
  VECCSE(kmax+1, 0., givens_s_vec, 0);
  VECCSE(kmax+1, 0., g_vec, 0);
  // Generates the initial basis of the Krylov subspace.
  B_FUNC(linear_problem, linear_problem_args, solution_vec, b_vec);
  // g_vec[0] = sqrt(b_vec^2)
  VECIN1(sqrt(VECDOT(dim_linear_problem, b_vec, 0, b_vec, 0)), g_vec, 0);
  // basis_mat[0] = b_vec / g_vec[0]
  VECCPSC(dim_linear_problem, 1/VECEX1(g_vec, 0), b_vec, 0, basis_vec, 0);

  // k : the dimension of the Krylov subspace at the current iteration.
  int k;
  for (k=0; k<kmax_; ++k) {
    AX_FUNC(linear_problem, linear_problem_args, basis_vec+k, basis_vec+(k+1));
    for (int j=0; j<=k; ++j) {
      MATIN1(VECDOT(dim_linear_problem, basis_vec, (k+1)*dim_linear_problem, 
                      basis_vec, j*dim_linear_problem), hessenberg_mat, k, j);
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