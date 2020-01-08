#include "d_nmpc_model.h"

#include <math.h>


#define REAL double
#define NMPC_MODEL d_nmpc_model
#define NMPC_MODEL_STRSIZE d_nmpc_model_strsize
#define NMPC_MODEL_CREATE d_nmpc_model_create
#define NMPC_MODEL_F d_nmpc_model_f
#define NMPC_MODEL_PHIX d_nmpc_model_phix
#define NMPC_MODEL_HX d_nmpc_model_hx
#define NMPC_MODEL_HU d_nmpc_model_hu
#define NMPC_MODEL_DIMX d_nmpc_model_dimx
#define NMPC_MODEL_DIMU d_nmpc_model_dimu
#define NMPC_MODEL_DIMC d_nmpc_model_dimc


#include "x_nmpc_model.c"