#include "commonincludes.h"




void print1D(double *arr, int dim, char *fname){

  FILE *fp;
  int i;

  fp = fopen(fname,"w");

  for(i = 0; i < dim; i++)
    fprintf(fp,"%lE\n",arr[i]);

  fclose(fp);

}


void print2D(double *x, double *y, int dim, char *fname){

  FILE *fp;
  int i;
  
  fp = fopen(fname,"w");

  for(i = 0; i < dim; i++)
    fprintf(fp,"%lf %lE\n",x[i],y[i]);

  fclose(fp);
}



void print3D(DATA_t *output, int n, char *fname){


  FILE *fp;
  int i, it = n, j;
  
  fp = fopen(fname,"w");

  if (n != 1) {
    for(i = 0; i < CC_STEP; i++){
      for (j = 0; j < it; j++){
	if (j <= it - 2){
	  fprintf(fp,"%lf ",output->psd[j][i]);
	} else {
	  fprintf(fp,"%lf\n",output->psd[j][i]);
	}
      }
    }
  } else {
    for(i = 0; i < CC_STEP; i++)
      fprintf(fp,"%lf\n",output->psd[0][i]);
  }
  fclose(fp);

}


char *printFname(char *str, int year){

  char *buf = malloc(50*sizeof(char));

  sprintf(buf,"%s%d.dat",str,year);
  return buf;
}
