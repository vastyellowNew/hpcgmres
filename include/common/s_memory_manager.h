#ifndef CGMRES_S_MEMORY_MANAGER_H_
#define CGMRES_S_MEMORY_MANAGER_H_


#ifdef __cplusplus
extern "C" {
#endif


float* allocate_svec(int dim);
float** allocate_smat(int dim_column, int dim_row);
void free_svec(float *dvec);
void free_smat(float **dmat);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif // CGMRES_S_MEMORY_MANAGER_H_