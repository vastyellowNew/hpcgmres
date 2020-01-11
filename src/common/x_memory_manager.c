REAL* ALLOCATE_VEC(int dim) {
  REAL *vec = (REAL *)malloc((size_t)(dim*sizeof(REAL)));
  if (!vec) {
    printf("allocation failuer in ALLOCATE_VEC.\n");
    exit(1);
  }
  int i;
  for (i=0; i<dim; ++i) {
    vec[i] = 0.0;
  }
  return vec;
}


REAL** ALLOCATE_MAT(int dim_column, int dim_row) {
  REAL **mat = (REAL **)malloc((size_t)(dim_column*sizeof(REAL*)));
  if (!mat) {
    printf("allocation failuer in ALLOCATE_MAT.\n");
    exit(1);
  }
  mat[0] = (REAL *)malloc((size_t)(dim_column*dim_row*sizeof(REAL)));
  if (!mat[0]) {
    printf("allocation failuer in ALLOCATE_MAT.\n");
    exit(1);
  }
  int i;
  for (i=0; i<dim_column-1; ++i) {
    mat[i+1] = mat[i] + dim_row;
  }
  for (i=0; i<dim_column*dim_row; ++i) {
    mat[0][i] = 0.0;
  }
  return mat;
}


void FREE_VEC(REAL *vec) {
  free(vec);
}


void FREE_MAT(REAL **mat) {
  free(mat[0]);
  free(mat);
}