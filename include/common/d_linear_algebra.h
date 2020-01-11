#ifndef HPCGMRES_D_LINEAR_ALGEBRA_H_
#define HPCGMRES_D_LINEAR_ALGEBRA_H_


#ifdef __cplusplus
extern "C" {
#endif


// Set value in x whose dimension is dim.
void hpcgmres_dvecset(int dim, double value, double *x);

// Copy x to y whose dimensions are dim.
void hpcgmres_dveccp(int dim, double *x, double *y);

// Copy ax to y whose dimensions are dim.
void hpcgmres_dvecmcp(int dim, double a, double *x, double *y);

// Return the dot products of x and y.
double hpcgmres_dvecdot(int dim, double *x, double *y);

// Return the squared norm of x.
double hpcgmres_dvecnrm2(int dim, double *x);

// Multyply the scalar a to each element of x.
void hpcgmres_dvecmul(int dim, double a, double *x);

// Add x to y elementwise.
void hpcgmres_dvecadd(int dim, double *x, double *y);

// Add ax to y elementwise.
void hpcgmres_dvecmadd(int dim, double a, double *x, double *y);

// Set ax + y to result.
void hpcgmres_daxpy(int dim, double a, double *x, double *y, double *result);

// Set ax + by to result.
void hpcgmres_daxpby(int dim, double a, double *x, double b, double *y, 
                     double *result);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // HPCGMRES_D_LINEAR_ALGEBRA_H_