#include "s_linear_algebra.h"


void hpcgmres_svecset(int dim, float value, float *x) {
  int i = 0;
  for ( ; i<dim-3; i+=4) {
    x[i  ] = value;
    x[i+1] = value;
    x[i+2] = value;
    x[i+3] = value;
  }
  for ( ; i<dim; i++) {
    x[i  ] = value;
  }
}


void hpcgmres_sveccp(int dim, float *x, float *y) {
  int i = 0;
  for ( ; i<dim-3; i+=4) {
    y[i  ] = x[i  ];
    y[i+1] = x[i+1];
    y[i+2] = x[i+2];
    y[i+3] = x[i+3];
  }
  for ( ; i<dim; i++) {
    y[i] = x[i];
  }
}


void hpcgmres_svecmcp(int dim, float a, float *x, float *y) {
  int i = 0;
  for ( ; i<dim-3; i+=4) {
    y[i  ] = a*x[i  ];
    y[i+1] = a*x[i+1];
    y[i+2] = a*x[i+2];
    y[i+3] = a*x[i+3];
  }
  for ( ; i<dim; i++) {
    y[i] = a*x[i];
  }
}


float hpcgmres_svecdot(int dim, float *x, float *y) {
  int i = 0;
  float dot = 0.0;
  for ( ; i<dim-3; i+=4) {
    dot += x[i  ] * y[i  ];
    dot += x[i+1] * y[i+1];
    dot += x[i+2] * y[i+2];
    dot += x[i+3] * y[i+3];
  }
  for ( ; i<dim; i++) {
    dot += x[i] * y[i];
  }
  return dot;
}


float hpcgmres_svecnrm2(int dim, float *x) {
  int i = 0;
  float nrm = 0.0;
  for ( ; i<dim-3; i+=4) {
    nrm += x[i  ] * x[i  ];
    nrm += x[i+1] * x[i+1];
    nrm += x[i+2] * x[i+2];
    nrm += x[i+3] * x[i+3];
  }
  for ( ; i<dim; i++) {
    nrm += x[i] * x[i];
  }
  return nrm;
}


void hpcgmres_svecmul(int dim, float a, float *x) {
  int i = 0;
  for ( ; i<dim-3; i+=4) {
    x[i  ] *= a;
    x[i+1] *= a;
    x[i+2] *= a;
    x[i+3] *= a;
  }
  for ( ; i<dim; i++) {
    x[i] *= a;
  }
}


void hpcgmres_svecadd(int dim, float *x, float *y) {
  int i = 0;
  for ( ; i<dim-3; i+=4) {
    y[i  ] += x[i  ];
    y[i+1] += x[i+1];
    y[i+2] += x[i+2];
    y[i+3] += x[i+3];
  }
  for ( ; i<dim; i++) {
    y[i] += x[i];
  }
}


void hpcgmres_svecmadd(int dim, float a, float *x, float *y) {
  int i = 0;
  for ( ; i<dim-3; i+=4) {
    y[i  ] += a*x[i  ];
    y[i+1] += a*x[i+1];
    y[i+2] += a*x[i+2];
    y[i+3] += a*x[i+3];
  }
  for ( ; i<dim; i++) {
    y[i] += a*x[i];
  }
}


void hpcgmres_saxpy(int dim, float a, float *x, float *y, float *result) {
  int i = 0;
  for ( ; i<dim-3; i+=4) {
    result[i  ] = a*x[i  ] + y[i  ];
    result[i+1] = a*x[i+1] + y[i+1];
    result[i+2] = a*x[i+2] + y[i+2];
    result[i+3] = a*x[i+3] + y[i+3];
  }
  for ( ; i<dim; i++) {
    result[i] = a*x[i] + y[i];
  }
}


void hpcgmres_saxpby(int dim, float a, float *x, float b, float *y, 
                     float *result) {
  int i = 0;
  for ( ; i<dim-3; i+=4) {
    result[i  ] = a*x[i  ] + b*y[i  ];
    result[i+1] = a*x[i+1] + b*y[i+1];
    result[i+2] = a*x[i+2] + b*y[i+2];
    result[i+3] = a*x[i+3] + b*y[i+3];
  }
  for ( ; i<dim; i++) {
    result[i] = a*x[i] + b*y[i];
  }
}