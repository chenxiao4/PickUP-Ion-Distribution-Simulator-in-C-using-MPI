#include "commonincludes.h"

extern double R_AU[IP_DIM];
extern double D_ISM[IP_DIM];

double dupwind(double r,PARA_t *input){

  double source;
  double rho, vr, vi;
  double part, index, index1;
  double loss;
  
  loss = input->ism.loss_rate;// + EILOSSAU;
 

  vi = input->ism.v;
  rho = r * RAU;
  vr = sqrt(pow(vi,2.)+2.*GM/rho);
  
  part = pow(RAU,2.)/GM;
  index1 = loss * pow(RAU,2.)/GM*vi;
  index = -loss * part * vr;

  source = input->ism.density * exp(index + index1);

  return source;
}



double ColdGas(double rho, double theta, PARA_t *input, int io){
  //note here theta is in rad
  /* cold gas model as a function of r and theta 
     theta ~ (0,180)
  */
 
  //io = 0 for density output
  //io = 1 for density ratio of direct orbit 

  double r, res;
  double c,lam,loss;
  double si, rsin, rout, b1, b2, db, db1, db2;
  double indir, dir;

  assert(theta != PI);

  loss = input->ism.loss_rate;
  r = rho * RAU;
  c = VISM_INF * VISM_INF / GM;
  lam = RAU * RAU / VISM_INF *loss;
  
  if (theta > 0. && theta < PI){

    si = sin(theta);
    rsin = r * si / 2.;
    rout = sqrt(rsin * rsin + 1./c * r * (1. - cos(theta)));
    b1 = rout + rsin;
    b2 = rout - rsin;

    db = (rsin*si + 1./c*(1. - cos(theta)))/rout/2.;
    db1 = db + si/2.;
    db2 = db - si/2.;

    indir = NINF * db2 *exp(-lam/b2*(2.*PI - theta)) / si;
    dir = NINF * db1 * exp(-lam/b1*theta) / si;

  } else if (theta == 0.){

    indir = 0.;
    dir = NINF * 0.5 * (1. + sqrt(1.+2./r/c)-sqrt(1./r/c/(r*c+2.))) * \
      exp(-2.*lam/r/(1.+sqrt(1.+2./r/c)));  
  }


    if (io == 0)
      res = indir + dir;
    else 
      res = dir / (indir + dir);
  return res;
}




double v_injection_long(double rho, double theta, PARA_t *input){

  double r, c, ecen, a, p1, p2;
  double vth2, vrt, vt1, vt2, perih, vr, vr2;
  double rt1,rt2;
  double vin;

  if (theta > PI) theta = 2*PI - theta;
 
  r = rho * RAU;
  vth2 = K * TINF / mhe * 1.e-6;
  c = VISM_INF * VISM_INF / GM;
  ecen = 2. * vth2 * c;
  vrt = sqrt(VISM_INF * VISM_INF + 2. * GM / r);

  if (theta != PI && theta != 0.) {
    a = 1. + 4./r/ecen/(1. + cos(theta));
    p1 = 0.5 * r * sin(theta) * (1+sqrt(a));
    p2 = 0.5 * r * sin(theta) * (sqrt(a)-1.);
    vt1 = p1 * VISM_INF /r;
    vt2 = p2 * VISM_INF /r;
    perih = acos(-1./sqrt(1.+pow(ecen,2.)*pow(p1,2.)));
    vr2 = sqrt(vrt * vrt - vt2 * vt2);

    if (theta == perih)
      vr = 0.;
    else if (theta < perih)
      vr = sqrt(vrt * vrt - vt1 * vt1);
    else 
      vr = -sqrt(vrt * vrt - vt1 * vt1);
  
      rt1 = ColdGas(rho,theta,input,1);
      rt2 = 1. - rt1;
      vin = sqrt(pow((input->sw.v + vr*rt1 - vr2*rt2),2.) + \
		 pow((vt1*rt1 - vt2*rt2),2.));

  } else if (theta == PI) {

    vin = v_injection_long(rho,0.99*PI,input);

  } else {

    vin = input->sw.v + vrt;

  }
  return vin;
}





double v_injection(double r, PARA_t *input){

  //injection speed in the upwind direction

  double vism, rho, vr ,vin;
  
  vism = input->ism.v;
  rho = r * RAU;
  
  vr = sqrt(pow(vism,2.) + 2.*GM/rho); 
  vin = input->sw.v + vr;
  return vin;
}





double dlnr(double r, PARA_t *input){

  double dln, tloss, densw, vism;

  tloss = input->ism.loss_rate;
  densw = input->sw.density;
  vism = input->ism.v;

  dln = -(tloss * pow(r,-2.) + ei_rate_loss(r,input)) * RAU 
    /sqrt(pow(vism,2.) + 2.*GM/RAU/r);
  
  return dln;

}


double trapzd(double a, double b, int n, PARA_t *input){

  int it, j;
  double del, sum, tnm, x;
  static double s;


  if (n == 1) {
    s = 0.5 * (b-a) * (dlnr(a,input) + dlnr(b,input));
  } else {
    for (it = 1, j = 1; j < n-1; j++)
      it <<= 1;
    tnm = it;
    del = (b-a) / tnm;
    x = a + 0.5*del;
    
    for (sum = 0., j = 1; j <= it; j++, x+=del)
      sum += dlnr(x,input);
    s = 0.5 * (s + (b-a)*sum/tnm);
  }//else

  return s;
}


double qsimp(double a, double b, PARA_t *input){
  
  int j;
  double s;
  double os,ost,st;

  ost = -1.E-30;
  os = 1.E-30;

  for(j = 1; j <= JMAX; j++){
    
    st = trapzd(a,b,j,input);
    s = (4.*st - ost) / 3.;

    if (fabs(s-os) < EPS*fabs(os))
      return s;    
    
    os = s;
    ost = st;

  }
  
  fprintf(stdout,"Too many steps in qsimp.Parameter a = %lf, b = %lf",a,b);
  return EXIT_FAILURE;
}


double dupwind_with_ei(double r, PARA_t *input){

  double dens,c0,c,r0,r1,r2,c1,d1,d0,dhe;

  r1 = 100.;
  r0 = 1.;
  r2 = 0.1;
  dhe = input->ism.density;
  
  c0 = qsimp(r1,r0,input);
  d0 = exp(-c0);
  if (r == 1.) {
    dens = d0;
    return dens*dhe;
  }

  c1 = qsimp(r0,r2,input);
  d1 = exp(-c1) * d0;
  if (r == 0.1) {
    dens = d1;
    return dens*dhe;
  }
 
  if (r > 0.1 && r < 1.){
    c = qsimp(r0,r,input);
    dens = exp(-c) * d0;
    return dens*dhe;

  }

  if (r < 0.1){
    c = qsimp(r2,r,input);
    dens = exp(-c) * d1;
    return dens*dhe;

  }

  if (r > 1.){
    c = qsimp(r1,r,input);
    dens = exp(-c);
    return dens*dhe;
  }

  exit(EXIT_SUCCESS);
}



double upwind_density(double r, PARA_t *input){

  double dens;

  if (EIRATE == 0){
    dens = dupwind(r,input);
  } else if (EIRATE == 1){
    dens = dupwind_with_ei(r,input);
  } else {
    err_exit("EIRATE is not initialized correctly.");
  }


  return dens;


}




void init_uwden(double *x, double *y,PARA_t *input){

  int i;
  double foo;

  foo = (IP_DIM-1.) / log(UPPER/LOWER);
  
  for (i = 0; i < IP_DIM; i++){

    x[i] = LOWER * exp( i / foo);

    y[i] = upwind_density(x[i],input);

    //   printf("Initi density : %lf\n",x[i]);

  }

}



double interp_spline(double xi,double *x, double *y){

  double res;



  gsl_interp_accel *acc 
    = gsl_interp_accel_alloc ();
  gsl_spline *spline 
    = gsl_spline_alloc (gsl_interp_cspline, IP_DIM);

  gsl_spline_init (spline, x, y, IP_DIM);
  res = gsl_spline_eval (spline, xi, acc);

  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);


  return res;
}




double interpol(double t,double *x, double *y)
{ int i,j,k,m,l;
  int n = IP_DIM;
  double z,h,b[8];
  z=0.0;
  if (n<1) return(z);
  if (n==1) { z=y[0]; return(z);}
  if (n<=8) { k=0; m=n;}
  else if (t<x[4]) { k=0; m=8;}
  else if (t>x[n-5]) { k=n-8; m=8;}
  else
    { k=1; j=n;
      while (j-k!=1)
	{ i=(k+j)/2;
	  if (t<x[i-1]) j=i;
	  else k=i;
	}
      k=k-4; m=8;
    }
  b[0]=y[k];
  for (i=2;i<=m;i++)
    { h=y[i+k-1]; l=0; j=1;
      while ((l==0)&&(j<=i-1))
	{ if (fabs(h-b[j-1])+1.0==1.0) l=1;
	  else h=(x[i+k-1]-x[j+k-1])/(h-b[j-1]);
	  j=j+1;
	}
      b[i-1]=h;
      if (l!=0) b[i-1]=1.0e+35;
    }
  z=b[m-1];
  for (i=m-1;i>=1;i--) z=b[i-1]+(t-x[i+k-1])/z;
  return(z);
}






double vs76(double w, double cl_index, PARA_t *input){

  double psd,norm, ionrate;
  double v_max, neutral, r;
  int flag;
  
  flag = EIRATE;

  if (w <= 0.)
    err_exit("W is not PROPER");


  r = pow(w,cl_index);
  v_max = v_injection(r,input);

  if (flag == OFF){
    //    neutral = dupwind(r,input); 
    neutral = upwind_density(r,input);
    psd = cl_index / 4. / PI * RAU / input->sw.v / pow(v_max,3.) *
      input->ism.ion_rate * neutral * pow(w,cl_index - 3.)*1.E15;

  } else if (flag == ON) { 
    //    neutral = dupwind_with_ei(r,input);
    //neutral = upwind_density(r,input);

    //neutral = interp_spline(r,R_AU,D_ISM);
    neutral = interpol(r,R_AU,D_ISM);
   
    norm = cl_index / 4. / PI * RAU / input->sw.v / pow(v_max,3.) *
      neutral * 1.E15;

    ionrate = (input->ism.ion_rate + ei_rate(r,input) * pow(r,2)) *
      pow(w,cl_index - 3.);

    psd = norm * ionrate;
  } else {
    fprintf(stderr,"Option in VS76 is not correct.\n");
    return EXIT_FAILURE;
  }

  return psd;
  
}
