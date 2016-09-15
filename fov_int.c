#include "commonincludes.h"


void fov_integration(int iter,double cl_index,DATA_t *output, PARA_t *input){

  double dvol,bas;
  double dphi,dthet,DGF,res,dvolu;
  int i,ie,it,jp,jt;
  double eq[CC_STEP], eql[CC_STEP], equ[CC_STEP];
  double phise[CC_TIME], phisl[CC_TIME], phisu[CC_TIME];
  double x,y,z,zz,z1;
  double fl,dfl,efl;
  double w1,w2,wori,dvi,deqi,vion,w_v;
  double phii,thet;
  double car[3],sph[3],GSE[3],SWF[3];
  double swsph[] = {input->sw.v,input->sw.phi,input->sw.theta}; 
  double swcar[3];
  double sf1;
  int IS;

  IS = 0;
  dvol = 0.1;
  bas = 1.0744;
  res = 0.064;
  DGF = 0.02;
  dphi = 10./input->grid.AN;
  dthet = 66./input->grid.AE;
  


  //initialize E-step
  for (i = 0; i < CC_STEP; i++){

    eq[i] = pow(bas,(i+1))*0.4271;   
    eql[i] = eq[i] * (1. - res / 2.);
    equ[i] = eq[i] * (1. + res / 2.);
  }

  
  for (i = 0; i < CC_TIME; i++){
    phise[i] = 2.5 + 5. * i;
    phisl[i] = phise[i] - 5.;
    phisu[i] = phise[i] + 5.;
  }

  //step through E-step and FOV

  for (ie = 0; ie < CC_STEP; ie++){
    
    w1 = sqrt(eql[ie] / MQ) * KEV_TO_KMS;//km/s not cm/s
    w2 = sqrt(equ[ie] / MQ) * KEV_TO_KMS;
    wori = sqrt(w1*w2);

    output->v[ie] = wori / input->sw.v;
    dvi = w2 - w1;
    deqi = dvi * wori;


    x = 0.;
    y = 0.;
    z = 0.;
    zz = 0.;
    z1 = 0.;
    
    double **SATGSE;
    SATGSE = cmatrix(3,3);
    
    for (it = 0; it < CC_TIME; it++){

      SCtoGSE(SATGSE,it+1,input);
      
      
      fl = 0.;
      efl = 0.;
      dfl = 0.;


      for (jp = 0; jp < input->grid.AN; jp++){

	phii = -2.5 + dphi * (jp + 0.5);
	
	for (jt = 0; jt < input->grid.AE; jt++){

	  thet = dthet * (jt + 0.5); 
	  dvolu = deqi * (cos((thet - dthet / 2.)*DEG_TO_RAD) - 
			  cos((thet + dthet /2.)*DEG_TO_RAD)) * 
	    dphi * DEG_TO_RAD * wori;

	  sph[0] = wori;
	  
	  if (phii < 0.){
	    sph[1] = 180. + phii;
	  } else {
	    sph[1] = phii - 180.;
	  }

	  sph[2] = 180. - thet;
	  
	  carsph(car,sph,1,0);
	  
	  //check("Step 1");
	  
	  coord_tran(GSE,SATGSE,car);

	  // check("Step 2");

	  carsph(swcar,swsph,0,0);
	  
	  //check("Step 3");
	  
	  vec_min(GSE,swcar,SWF);

	  // check("Step 4");

	  vion = vec_len(SWF);

	  //check("Step 5");

	  w_v = tran_w_to_sw(vion,cl_index,input);
	  //w_v = vion / input->sw.v;

	  //	  check("Step 6");
	  //fprintf(stdout,"ie[%d],jt[%d],jp[%d],it[%d]\n",ie,jt,jp,it);
	  
	  if (w_v < 1.) {

	    sf1 = vs76(w_v,cl_index,input);
	    //printf("f[%lf] = %lf\n",w_v,sf1);
	    dfl += dvolu * sf1 * wori;
	  }//if w_v < 1 
	}//jt
      }//jp
      
      zz += dfl;

      if ((it+1) % 9){
	IS++;
	x = 0.;
	y = 0.;
      }//if it+1 % 9
    }//it

    if (IS == 8)
      IS = 0;

    output->psd[iter][ie] = 1./ 8. * 378.4 * zz * 
      DGF /pow((eq[ie] / MQ),2) * 1.E-10;


    dmatrix(SATGSE);

    fprintf(stdout,"E-Step %d is Finished, PSD[%d] = %lf\n",ie,ie,output->psd[iter][ie]);

  }//ie

}
