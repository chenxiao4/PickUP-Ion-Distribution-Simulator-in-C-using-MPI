#include "commonincludes.h"



double **cmatrix(int n, int m)
{
  int i;

  double **u = malloc(n * sizeof(*u));
  u[0] = malloc(n * m * sizeof(*u[0]));
  for (i =1; i<n; i++)
    {
      u[i] = u[i-1] + m;
    }
  return u;
}



void dmatrix(double **u)
{
  free(u[0]);
  free(u);

}



double fw(double *x, double vi,double cl_index,PARA_t *input){

  double res, r;
  r = pow((*x),cl_index);

  res = (*x) - vi / v_injection(r,input);
  //  printf("fw[%lf] = %lf\n",*x,res);
  return res;
}





void Froot(double (*fn) (double *, double, double, PARA_t *), double *x, double vi, double cl_index, PARA_t *input){
  
  double x1, x2, yy,y1,y2;


  x1 = CA;
  x2 = CB;
  y1 = fw(&x1,vi,cl_index,input);
  y2 = fw(&x2,vi,cl_index,input);

  if (y2 < 0.){
    x2 *= 1.5;
    y2 = fw(&x2,vi,cl_index,input);
  }
  

  if (fabs(y1) < EPS){
    *x = x1;
    return;
  }
  /*
  if (fabs(y2) < EPS){
    *x = x2;
    return;
  }
  */
  /*
  while ((y1 < 0 && y2 < 0) || (y1 > 0 && y2 > 0)){
    x1 *= 1.E-2;
    x2 *= 2.;
    y1 = fw(&x1,vi,cl_index,input);
    y2 = fw(&x2,vi,cl_index,input);

    if (fabs(y1) < EPS){
      *x = x1;
      return;
    }
    if (fabs(y2) < EPS){
      *x = x2;
      return;
    }
    //message("W Transformation correction.");
    // printf("y1[%lE] = %lE, y2[%lf] = %lf\n",x1,y1,x2,y2);
  }//while
  */
 
  if ((y1 < 0 && y2 < 0) || (y1 > 0 && y2 > 0)){
    fprintf(stdout,"g[%lf] is %lE\n",x1,y1);
    err_exit("W transformation is fail due to the failure of solve the equation");
  }


  *x = (x1 + x2) /2.;
  yy = fw(x,vi,cl_index,input);
  if (fabs(yy) < EPS){
    return;
  }
  

  while (fabs(yy) > EPS){
    
    if (yy > 0.){    
      x2 = *x;
    } else {
      x1 = *x;    
    }
    *x = (x1 + x2) / 2.;
    yy = fw(x,vi,cl_index,input);
  }



}



double tran_w_to_sw(double vi, double cl_index, PARA_t *input){

  double res, temp;

  Froot(fw,&res,vi,cl_index,input);
  
  temp = fw(&res,vi,cl_index,input);

  if (fabs(temp) > 1.E-6)    
      err_exit("PUI Speed Transformation Fail.");
  
  // printf("w is : %lf\n",res);
  
  return res;
}



double expint(int n, double x){
	
  int i,ii,nm1;
  double a,b,c,d,del,fact,h,psi,ans;

  nm1=n-1;
  if (n < 0 || x < 0.0 || (x==0.0 && (n==0 || n==1))){
    printf("bad arguments in expint");
    exit(EXIT_FAILURE);
  } else {
    if (n == 0) ans=exp(-x)/x;
    else {
      if (x == 0.0) ans=1.0/nm1;

      else {
	if (x > 1.0) {
	  b=x+n;
	  c=1.0/FPMIN;
	  d=1.0/b;
	  h=d;
	  for (i=1;i<=MAXIT;i++) {
	    a = -i*(nm1+i);
	    b += 2.0;
	    d=1.0/(a*d+b);
	    c=b+a/c;
	    del=c*d;
	    h *= del;
	    if (fabs(del-1.0) < EPS) {
	      ans=h*exp(-x);
	      return ans;
	    }
	  }
	  printf("continued fraction failed in expint");
	  exit(EXIT_FAILURE);
	} else {
	  ans = (nm1!=0 ? 1.0/nm1 : -log(x)-EULER);
	  fact=1.0;
	  for (i=1;i<=MAXIT;i++) {
	    fact *= -x/i;
	    if (i != nm1) del = -fact/(i-nm1);
	    else {
	      psi = -EULER;
	      for (ii=1;ii<=nm1;ii++) psi += 1.0/ii;
	      del=fact*(-log(x)+psi);
	    }
	    ans += del;
	    if (fabs(del) < fabs(ans)*EPS) return ans;
	  }
	  printf("series failed in expint");
	  exit(EXIT_FAILURE);
	}
      }
    }
  }
  return ans;
}



void SCtoGSE(double **u, int t,PARA_t *input){

  double omeg, thet, phi, pphi, lon, time;

  omeg= PI / 6.;
  thet = input->SC.theta;
  lon = input->SC.phi;
  
  pphi = lon + PI / 2.;
  time = t / 6.;

  phi = omeg * time;

  u[0][0] = cos(pphi)*cos(phi) - cos(thet)*sin(pphi)*sin(phi);
  u[0][1] = -sin(phi)*cos(phi) - cos(thet)*cos(phi)*sin(pphi);
  u[0][2] = sin(thet)*sin(pphi);
  u[1][0] = cos(thet)*sin(phi)*cos(pphi) + sin(pphi)*cos(phi);
  u[1][1] = cos(thet)*cos(pphi)*cos(phi) - sin(pphi)*sin(phi);
  u[1][2] = -sin(thet)*cos(pphi);
  u[2][0] = sin(thet)*sin(phi);
  u[2][1] = sin(thet)*cos(phi);
  u[2][2] = cos(thet);

}





void coord_tran(double *p, double **mat, double *q){

  int i;

  for (i = 0; i < 3; i++)
    p[i] = mat[i][0]*q[0] + mat[i][1]*q[1] + mat[i][2]*q[2];

}


void carsph(double *car, double *sph, int flag, int io){

  //flag = 1: sph = [r,long,lat];

  double F, X1;
  
  F = 1.;
  if (flag == 1)
    F = DEG_TO_RAD;

  if (io == 0) {
    car[0] = sph[0]*sin(sph[2]*F)*cos(sph[1]*F);
    car[1] = sph[0]*sin(sph[2]*F)*sin(sph[1]*F);
    car[2] = sph[0]*cos(sph[2]*F);
  } else if (io == 1){
    X1 = car[0];
    if (X1 == 0.)
      X1 = 1.E-30;
    sph[0] = sqrt(vec_len(car));
    sph[2] = acos(car[2] / sph[0]) / F;
    sph[1] = atan2(car[1],X1) / F;
  } else {
    err_exit("No Such Mode Found in CARSPH().");
  }
}


double vec_len(double *vec){

  double l;

  l = sqrt(pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2));
  
  return l;
}


void vec_min(double *a, double *b, double *c){
  
  int i;

  for(i = 0; i < 3; i++)
    c[i] = a[i] - b[i];
 
}






