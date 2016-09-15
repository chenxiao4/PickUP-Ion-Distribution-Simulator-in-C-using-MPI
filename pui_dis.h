#ifndef PUI_DIS_H
#define PUI_DIS_H

#include "parameter.h"


//#define EPS 1.E-6
//#define JMAX 30



double dupwind(double r,PARA_t *input);
double ColdGas(double rho, double theta, PARA_t *input, int io);
double v_injection(double r, PARA_t *input);
double v_injection_long(double rho, double theta, PARA_t *input);
double vs76(double w, double cl_index, PARA_t *input);
double dlnr(double r, PARA_t *input);
double trapzd(double a, double b, int n, PARA_t *input);
double qsimp(double a, double b, PARA_t *input);
double dupwind_with_ei(double r, PARA_t *input);
double upwind_density(double r, PARA_t *input);
void init_uwden(double *x, double *y, PARA_t *input);
double interp_spline(double xi,double *x, double *y);
double interpol(double t, double *x, double *y);


#endif
