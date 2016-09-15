#include "commonincludes.h"

double R_AU[IP_DIM];
double D_ISM[IP_DIM];


int main(int argc, char *argv[]){

  /*
   * starting the normal mode
   *
   */


  if (!DEBUG) {
    
    int YST, YED;
    int opts[] = {0,0};
    int i,year;
    int counter = 0;

    double perc;

    PARA_t input;
    DATA_t *output;
    CONF_t init;

    int procs, pid;

    /*
     * MPI Initialization
     *
     */

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&procs);
    MPI_Comm_rank(MPI_COMM_WORLD,&pid);

    //initialization
    if (pid == 0){
      fprintf(stdout,"Having %d procs running in the analysis\n",procs);
      set_options(argc,argv,opts);
      init_config(opts,&init);
    }

    // root bcast the init structure
    Root_Bcat_data(&init,0);

    /*
    printf("Proc %d: \n",pid);
    printf("start: %d, end: %d, size: %d\n",init.start_year,init.end_year,init.size);
    printf("coolind index: ");
    for (int i = 0; i < 11; i++)
      printf("%f ",init.cooling_index[i]);
    printf("\n");
    */

    /*
    #pragma omp parallel
    {
        int thnum = omp_get_thread_num();
        int thtotal = omp_get_num_threads();
        printf("parallel: %d out of %d from proc %d out of %d\n",thnum,thtotal,pid,procs);

    }
    */

    YST = init.start_year;
    YED = init.end_year;  

    output = (DATA_t *)malloc(sizeof(DATA_t));
    output->psd = cmatrix(init.size,CC_STEP);
  

    //Simulation the psd in SWICS
    year = YST + pid;
   
    read_input(year, &input);
    output->year = year;
    
    if (EIRATE == ON)
      init_uwden(R_AU, D_ISM, &input);


    /*
    omp_set_num_threads(init.size);   
#pragma omp parallel
    {
      int tid = omp_get_thread_num();
      //int thtotal = omp_get_num_threads();
      //printf("parallel: %d out of %d from proc %d out of %d\n",tid,thtotal,pid,procs);
      */
    

    //#pragma omp parallel for private(i)
      for (i = 0; i < init.size; i++){
      
      double temp;
      temp = init.cooling_index[i];

      fov_integration(i,temp,output,&input);  

      //#pragma omp critical      
      counter++;
      perc = 100. * counter * 1.0 /init.size;
      
      fprintf(stdout,"COOLING INDEX is: %lf\n",temp);
      fprintf(stdout,"YEAR : %d completed: %lf%%. Whole Simulation completed: %lf%%\n",year,100. * (i+1.)/init.size,perc);

    }//i
      // }

    /*
     * output simulated results
     */
    //#pragma omp barrier
    if (OUT){
      print3D(output,init.size, printFname(NAME_PSD,year));
      print1D(output->v,CC_STEP, printFname(NAME_W,year));  
    }


    //free memory
    dmatrix(output->psd);   
  }// end of debug mode




  /*
   * starting the debug mode
   *
   */



  if (DEBUG) {

    int DIM = 100;
    PARA_t input;
    //    DATA_t output;
    //  double dhe;
    double arr[DIM], arr1[DIM];//,arr2[DIM];
    int i;

    //  double x;
  
    //  output.v = (double *)malloc(CC_STEP*sizeof(double));
    // output.psd = cmatrix(2,CC_STEP);

    //  init_uwden(R_AU, D_ISM,&input);
  
    read_input(1998, &input);

    //init_uwden(R_AU, D_ISM,&input);
    //    fov_integration(0,1.5,&output,&input);

  
    for (i = 0; i < DIM; i++){
      arr1[i] = 0.01*(i+1.);
      //arr[i] = upwind_density(arr1[i],&input);
      arr[i] = ei_rate(arr1[i], &input);
      //arr[i] = vs76(arr1[i], 1.5, &input);
      //arr[i] = ColdGas(arr1[i],0.,&input,0);
      //arr2[i] = input.ism.ion_rate * pow(arr1[i],-2);
      // arr2[i] = ei_rate_loss(arr1[i],&input);
    }
  
    print1D(arr,DIM,"eirate.txt");
    //print1D(D_ISM,IP_DIM,"updei");
    //print2D(R_AU,D_ISM,IP_DIM,"updei.txt");
    //print2D(arr1,arr,DIM,"eirate.txt");



    // dmatrix(output.psd);
    //  free(output.v);


    //Aitken(fn,&x);
    // Froot(fw,&x,500.,1.5,&input);
    // x =  tran_w_to_sw(500.,1.5,&input);
    //fprintf(stdout,"x is : %lf\n",x);
    // DATA_t mat;
    //  mat.data = create_matrix(10,10);
    // destroy_matrix(mat.data);

    //fprintf(stdout,"test expint(1,%lf) = %lf\n",0.5,expint(1,0.5));

  
    // dhe = dupwind(200.,&input);

    //dhe = v_injection(1., &input);


    // for (i = 0; i < DIM; i++){
    //   arr1[i] = (i+1)*0.01;
    //   arr[i] = ei_rate((i+1)*0.01,&input);
    //arr[i] = vs76((i+1)*0.01,1.5,&input,1);// * pow((i+1)*0.01,-2);
    // arr[i] = input.ism.ion_rate * pow((i+1)*0.01,-2);
    //arr[i] = fw(arr1 + i,500.,1.5,&input);
    //   printf("dendity[%d = %lf] is : %lf\n",i,i*0.01+0.01,arr[i]);
    // }

    // test(arr,DIM,"eio");

  }//debug



  MPI_Finalize();
  return EXIT_SUCCESS;
}









