#include "d_memory_manager.h"

#include <stdio.h>
#include <stdlib.h>


#define REAL double 
#define ALLOCATE_VEC allocate_dvec
#define ALLOCATE_MAT allocate_dmat
#define FREE_VEC free_dvec
#define FREE_MAT free_dmat


#include "x_memory_manager.c"