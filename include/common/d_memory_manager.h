#ifndef CGMRES_D_MEMORY_MANAGER_H_
#define CGMRES_D_MEMORY_MANAGER_H_


#ifdef __cplusplus
extern "C" {
#endif


double* allocate_dvec(int dim);
double** allocate_dmat(int dim_column, int dim_row);
void free_dvec(double *dvec);
void free_dmat(double **dmat);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // CGMRES_D_MEMORY_MANAGER_H_