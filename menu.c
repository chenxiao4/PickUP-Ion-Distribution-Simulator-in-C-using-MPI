#include "commonincludes.h"



char *input_fname(int year){

  char str[5];
  char *fname;

  
  fname = (char *)malloc((strlen(HFILE)+strlen(EXT)+5)*sizeof(char));
  sprintf(str,"%d",year);
  
  strcpy(fname,HFILE);
  strcat(fname,str);
  strcat(fname,EXT);
  
  if (strlen(fname) != strlen(HFILE)+strlen(EXT)+4)
    err_say("File name load failed.");

  return fname;
  
}




int trim_line(char *buf) {

   int len;
   
   len = strlen(buf);
   if (len > 0 && buf[len-1] == '\n')
     buf[--len] = '\0';
   return len;
}




void read_input(int year, PARA_t *input){

  FILE *fp;
  char *buf, *fname;
  double temp, tmp;
  int i;

  fname = input_fname(year);

  if ((buf = (char *)malloc(MAXLINE*sizeof(char))) == NULL)
    err_say("Cannot malloc the memory.");

  if ((fp  = fopen(fname,"r")) == NULL){
    perror(fname);
    err_say("Cannot open input file.");
  }


  fprintf(stdout,"\n%20s************************************************\n"," ");
  fprintf(stdout,"%50s\n","INPUT FILE");
  // read solar wind data
  fgets(buf,MAXLINE,fp);
  trim_line(buf);
  fprintf(stdout," %s\n",buf);
  
  fscanf(fp,"%lf %lf %lf %lf\n",&input->sw.v,&input->sw.phi,
	 &input->sw.theta,&input->sw.density);
  
  fprintf(stdout," %lf %lf %lf %lf\n",input->sw.v,input->sw.phi,
	  input->sw.theta,input->sw.density);


  //read IMF data
  fgets(buf,MAXLINE,fp);
  trim_line(buf);
  fprintf(stdout," %s\n",buf);
  
  fscanf(fp,"%lf %lf %lf\n",&input->imf.b,&input->imf.phi,
	 &input->imf.theta);
  fprintf(stdout," %lf %lf %lf\n",input->imf.b,input->imf.phi,
	  input->imf.theta);


  //read ISM data
  for (i = 0; i < 2; i++){   
    fgets(buf,MAXLINE,fp);
    trim_line(buf);
    fprintf(stdout," %s\n",buf);
  }

  fscanf(fp,"%lE %lf %lE %lE %lE\n",&input->ism.density,&input->ism.v,
	 &input->ism.loss_rate,&input->ism.ion_rate,&temp);

  input->ism.density = NINF;
  input->ism.v = VISM_INF;

  fprintf(stdout," %lE %lf %lE %lE %lE\n",input->ism.density,input->ism.v,
	  input->ism.loss_rate,input->ism.ion_rate,temp);



  for (i = 0; i < 8; i++){   
    fgets(buf,MAXLINE,fp);
    trim_line(buf);
    fprintf(stdout," %s\n",buf);
  }

  fscanf(fp,"%lf %lf %lf\n",&input->SC.theta, &input->SC.phi, &temp);
  fprintf(stdout," %lf %lf %lf\n",input->SC.theta,input->SC.phi,temp);

  for (i = 0; i < 3; i++){   
    fgets(buf,MAXLINE,fp);
    trim_line(buf);
    fprintf(stdout," %s\n",buf);
  }  

  fscanf(fp,"%d %d %d\n",&input->grid.AM, &input->grid.AN,&input->grid.AE);
  fprintf(stdout," %10d %10d %10d\n",input->grid.AM,input->grid.AN,
	  input->grid.AE);


  for (i = 0; i < 1; i++){   
    fgets(buf,MAXLINE,fp);
    trim_line(buf);
    fprintf(stdout," %s\n",buf);
  } 

  fscanf(fp,"%lf %lf %lf\n",&input->ism.cl_index,&tmp,&temp);
  fprintf(stdout," %lf %lf %lf\n",input->ism.cl_index,tmp,temp);

  fprintf(stdout,"\n%20s************************************************\n"," ");


  fclose(fp);
  free(buf);
}
