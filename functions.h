#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "parameter.h"

//#define EPS 1.E-6
//#define CA 1.E-4
//#define CB 10.
//#define MAXIT 100
//#define EULER 0.5772156649
//#define FPMIN 1.0e-30




double **cmatrix(int n, int m);
void dmatrix(double **u);

void Froot(double (*fn) (double *, double, double, PARA_t *), 
	    double *x, double vi, double cl_index, PARA_t *input);

double fw(double *x, double vi,double cl_index,PARA_t *input);
double tran_w_to_sw(double vi, double cl_index, PARA_t *input);
double expint(int n, double x);
void SCtoGSE(double **u, int t,PARA_t *input);
void coord_tran(double *p, double **mat, double *q);
void carsph(double *car, double *sph, int flag, int io);
double vec_len(double *vec);
void vec_min(double *a, double *b, double *c);



#endif
