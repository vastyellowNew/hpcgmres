#include "common/s_memory_manager.h"

#include <stdio.h>
#include <stdlib.h>


#define REAL float
#define ALLOCATE_VEC allocate_svec
#define ALLOCATE_MAT allocate_smat
#define FREE_VEC free_svec
#define FREE_MAT free_smat


#include "x_memory_manager.c"