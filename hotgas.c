#include "commonincludes.h"


double mbsi(int n,double x, int iop)
  { int i,m;
    double t,y,p,b0,b1,q;
    static double a[7]={ 1.0,3.5156229,3.0899424,1.2067492,
                         0.2659732,0.0360768,0.0045813};
    static double b[7]={ 0.5,0.87890594,0.51498869,
              0.15084934,0.02658773,0.00301532,0.00032411};
    static double c[9]={ 0.39894228,0.01328592,0.00225319,
                        -0.00157565,0.00916281,-0.02057706,
                         0.02635537,-0.01647633,0.00392377};
    static double d[9]={ 0.39894228,-0.03988024,-0.00362018,
                        0.00163801,-0.01031555,0.02282967,
                        -0.02895312,0.01787654,-0.00420059};
    if (n<0) n=-n;
    t=fabs(x);
    if (n!=1)
      { if (t<3.75)
          { y=(x/3.75)*(x/3.75); p=a[6];
            for (i=5; i>=0; i--)
              p=p*y+a[i];
          }
        else
          { y=3.75/t; p=c[8];
            for (i=7; i>=0; i--)
              p=p*y+c[i];
            p=p*exp(t)/sqrt(t);
          }
      }
    if (n==0) return(p);
    q=p;
    if (t<3.75)
      { y=(x/3.75)*(x/3.75); p=b[6];
        for (i=5; i>=0; i--) p=p*y+b[i];
        p=p*t;
      }
    else
      { y=3.75/t; p=d[8];
        for (i=7; i>=0; i--) p=p*y+d[i];
        p=p*exp(t)/sqrt(t);
      }
    if (x<0.0) p=-p;
    if (n==1) return(p);
    if (x==0.0) return(0.0);
    y=2.0/t; t=0.0; b1=1.0; b0=0.0;
    m=n+(int)sqrt(40.0*n);
    m=2*m;
    for (i=m; i>0; i--)
      { p=b0+i*y*b1; b0=b1; b1=p;
        if (fabs(b1)>1.0e+10)
          { t=t*1.0e-10; b0=b0*1.0e-10;
            b1=b1*1.0e-10;
          }
        if (i==n) t=b0;
      }
    p=t*q/b1;
    if ((x<0.0)&&(n%2==1)) p=-p;

    
    if (iop == 2) p = exp(-fabs(x))*p;
    return(p);
  }


double hedist(double r, double theta, double x1, double x2, PARA_t *input)

{
  double tdist = 0., rr, Q2,v0;
  double r01,vbsin,vbcos,factor;
  double calph,salph,C,calphm,A;
  double sgn,sign,p1,p,yy,zz,E,vmomen,val,r0,r1;
  double theta1,theta0,thpri,thprif,term1,term2;
  double phloss,sinal,sinal1,cosal,bes0;
  double F2,F1,F0;
  int nstep,iop;
  double loss;

  loss = input->ism.loss_rate;

  
  /**************************************/
  rr = r * RAU;
  Q2 = GM/rr;
  r01 = rr/Q2;

  vbsin = -fm1 * sin(theta);
  vbcos = -fm1 * cos(theta);

  v0 = x1 * scale * vt;
  factor = (1 - scale*scale)*x1*x1;

  calph = cos(x2);
  salph = sin(x2);
  C = v0*v0/GM;
  
  calphm = 1. - calph;
  A = (calphm  == 0.)? 1.e36 : (1. + 4./rr/C/calphm);
  sgn = 1.;
  p1 = 0.5 * rr * salph;

  /****************************************/
  for (nstep = 1; nstep < 3; nstep++)
    {
      p = (x2 == 0.)? (sgn*p1 + sqrt(p1*p1 + rr*(1.+calph)/C)) : (sgn*p1 + sqrt(A)*p1);

      E = (sqrt(A) + sgn*2. + 1./sqrt(A))/4.;//function nj*
      yy = fabs(p*v0/rr);
      vmomen = rr*yy;
      val = v0*yy/Q2;
      r0 = r01*yy*yy;
      r1 = GM/v0/v0;
  
      term1 = sqrt(r0/r1 + 2.*r0/rr - r0*r0/rr/rr);
      term2 = r0/rr - 1.;

      theta1 = (atan(term1/term2) < 0.)? (M_PI + atan(term1/term2)) : atan(term1/term2);
      theta0 = (atan(-val) < 0.)? (M_PI + atan(-val)) : atan(-val);
      thpri = M_PI - sgn*x2;

      sign = -1.;
      if (sgn == -1.) sign = 1.;
      if (sgn == 1. && thpri > theta0) sign = 1.;
      if (sgn == 1. && x2 == M_PI) sign = -1.;

      zz = sign * sqrt(v0*v0 + 2.*Q2 - yy*yy);
      thprif = theta0 + sign*theta1;
  
      phloss = (yy == 0.)? (-loss*RAU*RAU*2./rr/(v0*(1.+sqrt(1.+2./rr/C)))) : (-loss*RAU*RAU*thprif/vmomen);

      sinal = yy*(v0-zz)/(v0*(v0-zz)+Q2);
      F2 = vbsin*v0*sinal;

      iop = (fabs(F2) <= 80.)? 1 : 2;
      bes0 = 2.*M_PI*mbsi(0,F2,iop);
  
      sinal1 = (v0-zz)*(v0-zz)/(v0*(v0-zz)+Q2);
      cosal = 1. - sinal1;
  
      F1 = -vbcos*v0*cosal + phloss + factor;
  
      if (F1 < -77.)
	{
	  F0 = 0.;
	}
      else
	{
	  if (iop == 1) F0 = cf0*exp(F1)*bes0;
	  if (iop == 2) F0 = cf0*exp(fabs(F2))*exp(F1)*bes0;
	  if (F0 < 1.e-35) F0 = 0.;
	}
  
      tdist += F0*E*x1*salph;
      sgn = -1.;
    }
    
  return (tdist*N0);
}





/*************************************************/
/******/double nhedr(double r, double theta, PARA_t *input)/****/
{
 
  double dist,sum1,sum2;
  double x1,x2;
  double *fvinal;
  int i,j;
  
//  printf("test config:%g\n",config.HeDeninf);
  fvinal = (double *)malloc(N*sizeof(double));
  

  sum1 = 0.;
  for (i=0; i<N; i++)
    {
      sum2 = 0.;
      x1 = xvinf[i];

      for (j=0; j<N; j++)
	{
	  x2 = M_PI/2. + M_PI*root[j]/2.;
	  fvinal[j] = hedist(r,theta,x1,x2,input);
	  sum2 += weight[j]*fvinal[j];
	}
      
      sum2 = sum2 * M_PI/2.;
      sum1 += wvinf[i]*sum2;
    }
  
  dist = sum1;
  //printf("test inside nhedr:%g,%g,%g\n",r,theta,dist);
  free(fvinal);
  return (dist);


}
