#ifndef HPCGMRES_S_LINEAR_ALGEBRA_H_
#define HPCGMRES_S_LINEAR_ALGEBRA_H_


#ifdef __cplusplus
extern "C" {
#endif


// Set value in x whose dimension is dim.
void hpcgmres_svecset(int dim, float value, float *x);

// Copy x to y whose dimensions are dim.
void hpcgmres_sveccp(int dim, float *x, float *y);

// Copy ax to y whose dimensions are dim.
void hpcgmres_svecmcp(int dim, float a, float *x, float *y);

// Return the dot products of x and y.
float hpcgmres_svecdot(int dim, float *x, float *y);

// Return the squared norm of x.
float hpcgmres_svecnrm2(int dim, float *x);

// Multyply the scalar a to each element of x.
void hpcgmres_svecmul(int dim, float a, float *x);

// Add x to y elementwise.
void hpcgmres_svecadd(int dim, float *x, float *y);

// Add ax to y elementwise.
void hpcgmres_svecmadd(int dim, float a, float *x, float *y);

// Set ax + y to result.
void hpcgmres_saxpy(int dim, float a, float *x, float *y, float *result);

// Set ax + by to result.
void hpcgmres_saxpby(int dim, float a, float *x, float b, float *y, 
                     float *result);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // HPCGMRES_S_LINEAR_ALGEBRA_H_