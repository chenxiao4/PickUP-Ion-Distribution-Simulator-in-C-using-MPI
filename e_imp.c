#include "commonincludes.h"


void Tradial(double *ptr, double vsw){

  if (ISOTHERMAL){

    ptr[0] = 0;
    ptr[1] = 0;
    return;

  }


  if (vsw <= 400.){
    
    ptr[0] = 0.527;
    ptr[1] = 0.65;
  } else if (vsw > 400. && vsw <= 500.){
    
    ptr[0] = 0.394;
    ptr[1] = 0.68;
  } else if (vsw > 500 && vsw <= 600){

    ptr[0] = 0.2;
    ptr[1] = 0.767;

  } else if (vsw > 600. && vsw <= 700){

    ptr[0] = 0.226;
    ptr[1] = 0.805;
  } else if (vsw > 700 && vsw <= 800){

    ptr[0] = 0.296;
    ptr[1] = 0.812;

  } else {

    ptr[0] = 0.389;
    ptr[1] = 0.825;
  }


}



double voronvo(double r, PARA_t *input, int flag){
  
  // 0: loss rate
  // 1: ion rate

  double voronvo,densw = 0.;
  double de,p,a,x,k,kb;
  double den0core,den0hal,t0core,t0hal;
  double coreind,halind;
  double dencore,denhal,tcore = 0,thal = 0.;
  double tec,teh,uc,uh,ratc,rath;
  double tind[2];

  de = 24.6;
  p = 0.;
  a = 0.175E-7;
  x = 0.18;
  k = 0.35;
  kb = 8.617385E-5;


  if (flag == 1){//ion rate 
    Tradial(tind, input->sw.v);

    densw = input->sw.density;

    //t0core = 1.21E5;
    //t0hal = 9.42E5;

    t0core = TCOREAU;
    t0hal = THALOAU;

  } else if (flag == 0){//loss

    //Tradial(tind, 437.);
    Tradial(tind, input->sw.v);
    densw = 5.;

    t0core = TCOREAU;
    t0hal = THALOAU;

  } else {

    printf("In file %s : %d: flag = %d not found\n",__FILE__,__LINE__,flag);

  }



    
    
  den0core = (1. + 2.*0.04) / (1. + 0.05) * densw;
  den0hal = 0.05 * den0core;


  if (r >= 0.3 && r <= 1.) {
    coreind = tind[0];
    halind = coreind;
    tcore = t0core * pow(r,-coreind);
    thal = t0hal * pow(r,-halind);

  } else if (r < 0.3){
    coreind = tind[1];
    halind = coreind;
    tcore = t0core * pow(0.3,-tind[0]) * pow((r/0.3),-coreind);
    thal = t0hal * pow(0.3,-tind[0]) * pow((r/0.3),-halind);

  } else {

    printf("In file %s : %d: r = %f not found\n",__FILE__,__LINE__,r);

  }
    
    
  //printf("r[%lf]--- Tc[%.2lf],Th[%.2lf]\n",r,tcore,thal);
    
  dencore = den0core * pow(r, -2.);
  denhal = den0hal * pow(r, -2.);

  tec = tcore * kb;
  teh = thal * kb;
  uc = de / tec;
  uh = de / teh;
  
  ratc = dencore * a * (1.+p*sqrt(uc)) * pow(uc,k) * exp(-uc) / (x+uc);
  rath = denhal * a * (1.+p*sqrt(uh)) * pow(uh,k) * exp(-uh) / (x+uh);
  voronvo = ratc + rath;

  //printf("r[%lf]--- ratc[%.2lf],rath[%.2lf]\n",r,ratc,rath);
    
  return voronvo;
}






double rucinski(double r, PARA_t *input, int flag){

  double res,densw;
  double den0core,den0hal,t0core,t0hal;
  double coreind,halind;
  double dencore,denhal,tcore = 0.,thal = 0.;
  double kb,kb1,me,HeP,conv;
  double PIKTC,PIKTH,A1,C1,C2,G1,G2,ALPH1,ALPH2,PHI1,PHI2;
  double tind[2];

  Tradial(tind, input->sw.v);

  kb1 = 1.380658E-16;
  kb = 8.617385E-5;
  me = 9.10938987E-28;
  densw = input->sw.density;

  den0core = (1. + 2.*0.04) / (1. + 0.05) * densw;
  den0hal = 0.05 * den0core;

  if (flag == 1){

    t0core = TCOREAU; // Should be loaded form observation
    t0hal = THALOAU; // EI_t structure contains this data

  } 

  if (flag == 0){

    t0core = TCOREAU; // Should be loaded form observation
    t0hal = THALOAU; // EI_t structure contains this data

  }
  
  HeP = 2.459E1;
  conv = 1.60217733E-12;

  dencore = den0core * pow(r,-2.);
  denhal = den0hal * pow(r,-2);

  if (r >= 0.3 && r <= 1.){
    coreind = tind[0];
    halind = tind[0];
    tcore = t0core * pow(r,-coreind);
    thal = t0hal * pow(r,-halind);

  } else if (r < 0.3){
    coreind = tind[1];
    halind = tind[1]; 
    tcore = t0core * pow(0.3,-tind[0]) * pow((r/0.3),-coreind);
    thal = t0hal * pow(0.3,-tind[0]) * pow((r/0.3),-halind);
  } else {

    printf("In file %s : %d: flag = %f not found\n",__FILE__,__LINE__,r);
    
  }

  PIKTC = PI * kb1 * tcore;
  PIKTH = PI * kb1 * thal;

  A1 = 8. * PI / pow(me,2.);
  C1 = dencore * pow((me/2./PIKTC),1.5);
  C2 = denhal * pow((me/2./PIKTH),1.5);
  G1 = 0.75 * exp(0.46) * C1;
  G2 = 0.75 * exp(0.46) * C2;

  ALPH1 = HeP / kb / tcore;
  ALPH2 = HeP / kb / thal;

  PHI1 = ALPH1 + 0.46;
  PHI2 = ALPH2 + 0.46;
  

  res = 8.E-14*A1*(C1/ALPH1*expint(1,ALPH1)-G1/PHI1*expint(1,PHI1)+
		   C2/ALPH2*expint(1,ALPH2)-G2/PHI2*expint(1,PHI2));
  return res * pow(conv,2.);
}





double ei_rate(double r, PARA_t *input){

  double rate = 0., r0;
  char flag;

  flag = EIMODE;
  r0 = 1.;


  if (flag == 'S'){

    return voronvo(r0,input,1) * pow(r,-EICINDEX);

  }
 
  if(flag == 'V') {
    if(r > 1.) {
      rate = voronvo(r0,input,1) * pow(r,-FAR_IND);
    } else {
      rate = voronvo(r,input,1);
    }
  } else if (flag == 'R'){
    if(r > 1.) {
      rate =rucinski(r0,input,1) * pow(r,-FAR_IND);
    } else {
      rate = rucinski(r,input,1);
    }
  } else if (flag == (char)0) {
    //   err_exit("No Such EI MODE Found.");
    rate = 0.;
  } else {
    err_exit("No Such EI MODE Found.");
  }
  return rate;
}


double ei_rate_loss(double r, PARA_t *input){

  double rate = 0., r0;
  double temp;
  char flag;

  flag = EIMODE;
  r0 = 1.;


  if (flag == 'S'){

    return EILOSSAU * pow(r,-EICINDEX);

  }
 
  if(flag == 'V') {

    temp = voronvo(r0,input,0);
    if(r > 1.) {
      rate = EILOSSAU * pow(r,-FAR_IND);
    } else {
      rate = voronvo(r,input,0) * EILOSSAU / temp;
    }

  } else if (flag == 'R'){

    temp = rucinski(r0,input,0);
    if(r > 1.) {
      rate = EILOSSAU * pow(r,-FAR_IND);
    } else {
      rate = rucinski(r,input,0) * EILOSSAU / temp;
    }

  } else if (flag == (char)0) {
    //   err_exit("No Such EI MODE Found.");
    rate = 0.;
  } else {
    err_exit("No Such EI MODE Found.");
  }
  return rate;
 
}
