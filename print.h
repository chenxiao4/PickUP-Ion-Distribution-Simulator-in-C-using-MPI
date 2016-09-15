#ifndef PRINT_H
#define PRINT_H

#include "parameter.h"




void print1D(double *arr, int dim, char *fname);
void print2D(double *x, double *y, int dim, char *fname);
void print3D(DATA_t *output, int n, char *fname);
char *printFname(char *str, int year);

#endif
