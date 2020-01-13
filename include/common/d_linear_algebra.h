#ifndef CGMRES_D_LINEAR_ALGEBRA_H_
#define CGMRES_D_LINEAR_ALGEBRA_H_


#ifdef __cplusplus
extern "C" {
#endif


// Set value in x whose dimension is dim.
void cgmres_dvecset(int dim, double value, double *x);

// Copy x to y whose dimensions are dim.
void cgmres_dveccp(int dim, double *x, double *y);

// Copy ax to y whose dimensions are dim.
void cgmres_dvecmcp(int dim, double a, double *x, double *y);

// Return the dot products of x and y.
double cgmres_dvecdot(int dim, double *x, double *y);

// Return the squared norm of x.
double cgmres_dvecnrm2(int dim, double *x);

// Multyply the scalar a to each element of x.
void cgmres_dvecmul(int dim, double a, double *x);

// Add x to y elementwise.
void cgmres_dvecadd(int dim, double *x, double *y);

// Add ax to y elementwise.
void cgmres_dvecmadd(int dim, double a, double *x, double *y);

// Set ax + y to result.
void cgmres_daxpy(int dim, double a, double *x, double *y, double *result);

// Set ax + by to result.
void cgmres_daxpby(int dim, double a, double *x, double b, double *y, 
                   double *result);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // CGMRES_D_LINEAR_ALGEBRA_H_