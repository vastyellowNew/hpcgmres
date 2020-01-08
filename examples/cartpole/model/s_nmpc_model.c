#include "s_nmpc_model.h"

#include <math.h>


#define REAL float
#define NMPC_MODEL s_nmpc_model
#define NMPC_MODEL_STRSIZE s_nmpc_model_strsize
#define NMPC_MODEL_CREATE s_nmpc_model_create
#define NMPC_MODEL_F s_nmpc_model_f
#define NMPC_MODEL_PHIX s_nmpc_model_phix
#define NMPC_MODEL_HX s_nmpc_model_hx
#define NMPC_MODEL_HU s_nmpc_model_hu
#define NMPC_MODEL_DIMX s_nmpc_model_dimx
#define NMPC_MODEL_DIMU s_nmpc_model_dimu
#define NMPC_MODEL_DIMC s_nmpc_model_dimc


#include "x_nmpc_model.c"