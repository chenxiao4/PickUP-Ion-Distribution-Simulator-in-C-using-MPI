#ifndef PARAMETER_H
#define PARAMETER_H


#define CC_STEP 76
#define CC_SECTOR 8
#define CC_TIME 72

typedef struct{

  int start_year;
  int end_year;
  int size;
  double cooling_index[11];
}CONF_t;



typedef struct{

  int year;
  double **psd;
  double v[CC_STEP];

} DATA_t;




typedef struct{

  double v;
  double theta;
  double phi;
  double density;

} SW_t;


typedef struct {

  double b;
  double theta;
  double phi;

} IMF_t;


typedef struct {

  double v;
  double density;
  double ion_rate;
  double loss_rate;
  double cl_index;

} ISM_t;


typedef struct {

  double theta;
  double phi;

} SC_t;

typedef struct{

  double TcoreAU;
  double ThalAU;

} EI_t;


typedef struct {

  int AM;
  int AN;
  int AE;

} GRID_t;


typedef struct {

  SW_t sw;
  IMF_t imf;
  ISM_t ism;
  SC_t SC;
  EI_t eimp;
  GRID_t grid;

} PARA_t;




//OPTION
#define OPTIONS ":hspc:"
#define WHITE_SPACE " \n\r\t\v\f"


#define PI 3.14159265
#define DEG_TO_RAD PI/180.
#define RAD_TO_DEG 180./PI


//pui_dis
#define EPS 1.E-6
#define JMAX 30

//neutral gas parameter
#define NINF 1.5E-2
#define VISM_INF 23.2
#define TINF 6300.
#define RAU 1.496e8
#define GM 1.32743e11


//functions
#define EPS 1.E-6
#define CA 1.E-6
#define CB 10.
#define MAXIT 100
#define EULER 0.5772156649
#define FPMIN 1.0e-30

//e_imp
#define IN_IND 0.394 // r = [0.3,1]
#define OUT_IND 0.68 // r = [,0.3]
#define FAR_IND 2.2

//FOV Integration
#define KEV_TO_KMS 4.4E2
#define MQ 4



//define interpolate parameters when EIRATE is on

#define IP_DIM 100
#define LOWER 1.E-3
#define UPPER 1.1


//define output name
/* if you need to change the input data file name, goto menu.h */
#define NAME_PSD "SWICS_PER_"
#define NAME_W "SWICS_WPER_"


/*****************define some modes**********************************/
/********************************************************************/

//debug
#define OFF      0
#define ON       1
#define DEBUG    OFF

#define OUT      ON
#define EIRATE   ON //0:off,1:in
#define EICINDEX 2.0
#define EIMODE   'V' //V = Voronvo, R=Rucinski, S=define a cooling rate 
#define GASMODE  'C' //C = Cold, H = Hot
#define EILOSSAU  1.5E-8 // 1/s
#define TCOREAU   1.5E5 //  K  used to calculate EI loss
#define THALOAU   7.E5  //  K  used to calculate EI loss
#define ISOTHERMAL OFF
/********************************************************************/
#endif
