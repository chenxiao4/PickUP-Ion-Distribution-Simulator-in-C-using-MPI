#ifndef E_IMP_H
#define E_IMP_H

#include "parameter.h"

//#define IN_IND 0.394 // r = [0.3,1]
//#define OUT_IND 0.68 // r = [,0.3]
//#define FAR_IND 2.2





double voronvo(double r, PARA_t *input, int flag);
double ei_rate(double r, PARA_t *input);
double rucinski(double r, PARA_t *input, int flag);
void Tradial(double *ptr, double vsw);
double ei_rate_loss(double r, PARA_t *input);



#endif
