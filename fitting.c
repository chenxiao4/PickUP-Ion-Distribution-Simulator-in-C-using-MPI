#include "commonincludes.h"



void linfit(double *x, double *y, int n,double *fit){

  double c0, c1, cov00, cov01, cov11, sumsq;
  gsl_fit_linear(x, 1, y, 1, n, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
  
  fit[0] = c1;
  fit[1] = c0;
  
  //printf ("# best fit: Y = %g + %g X\n", c0, c1);

}


void linfit_w(double *x, double *y,double *err, int n,double *fit){

  int i;
  double c0, c1, cov00, cov01, cov11, chisq;
  double *w;

  w = (double *)malloc(n*sizeof(double));
  
  for (i = 0; i < n; i++)
    w[i] = 1. / err[i];

  gsl_fit_wlinear (x, 1, w, 1, y, 1, n, 
                   &c0, &c1, &cov00, &cov01, &cov11, 
                   &chisq);


  fit[0] = c1;
  fit[1] = c0;

}
