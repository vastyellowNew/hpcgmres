#ifndef CGMRES_S_LINEAR_ALGEBRA_H_
#define CGMRES_S_LINEAR_ALGEBRA_H_


#ifdef __cplusplus
extern "C" {
#endif


// Set value in x whose dimension is dim.
void cgmres_svecset(int dim, float value, float *x);

// Copy x to y whose dimensions are dim.
void cgmres_sveccp(int dim, float *x, float *y);

// Copy ax to y whose dimensions are dim.
void cgmres_svecmcp(int dim, float a, float *x, float *y);

// Return the dot products of x and y.
float cgmres_svecdot(int dim, float *x, float *y);

// Return the squared norm of x.
float cgmres_svecnrm2(int dim, float *x);

// Multyply the scalar a to each element of x.
void cgmres_svecmul(int dim, float a, float *x);

// Add x to y elementwise.
void cgmres_svecadd(int dim, float *x, float *y);

// Add ax to y elementwise.
void cgmres_svecmadd(int dim, float a, float *x, float *y);

// Set ax + y to result.
void cgmres_saxpy(int dim, float a, float *x, float *y, float *result);

// Set ax + by to result.
void cgmres_saxpby(int dim, float a, float *x, float b, float *y, 
                   float *result);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // CGMRES_S_LINEAR_ALGEBRA_H_