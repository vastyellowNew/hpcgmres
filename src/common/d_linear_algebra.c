#include "d_linear_algebra.h"


void hpcgmres_dvecset(int dim, double value, double *x) {
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


void hpcgmres_dveccp(int dim, double *x, double *y) {
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


void hpcgmres_dvecmcp(int dim, double a, double *x, double *y) {
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


double hpcgmres_dvecdot(int dim, double *x, double *y) {
  int i = 0;
  double dot = 0.0;
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


double hpcgmres_dvecnrm2(int dim, double *x) {
  int i = 0;
  double nrm = 0.0;
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


void hpcgmres_dvecmul(int dim, double a, double *x) {
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


void hpcgmres_dvecadd(int dim, double *x, double *y) {
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


void hpcgmres_dvecmadd(int dim, double a, double *x, double *y) {
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


void hpcgmres_daxpy(int dim, double a, double *x, double *y, double *result) {
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


void hpcgmres_daxpby(int dim, double a, double *x, double b, double *y, 
                     double *result) {
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