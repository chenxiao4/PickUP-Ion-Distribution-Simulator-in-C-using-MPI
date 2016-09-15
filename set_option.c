#include "commonincludes.h"




void set_options(int argc, char **argv, int *opts){

  int c;
  
  while ((c = getopt(argc,argv,OPTIONS)) != -1)
    {
      switch(c) {
      case 'h'://print help info
	printf("\nUsage: ./pui [-hsp] [-c num]\n\n");
	printf(" -h: this is the help information.\n");
	printf(" -s: Do the simulation for some year.\n");
	printf(" -p: Calculate the predicted PSD for some cooling index in some year.\n");
	break;
	exit(EXIT_SUCCESS);
      case 's':     /* Select input line only if pattern string starts
		       at the beginning of the line.*/
	opts[0] = 1;
	break;
      case 'p':    /* Select input line only if pattern string ends
		      at the end of the line. */
	opts[1] = 1;
	break;
      case ':':
	fprintf(stderr,"Option -%c requires an argument.\n",optopt);
	break;
	exit(EXIT_SUCCESS);
      case '?':
	fprintf(stderr,"Illegal or Unknownn option '%c' \n",optopt);
	break;
	exit(EXIT_SUCCESS);
      }
    }

}



void init_config(int *opts, CONF_t *init){

  int pflag,sflag;
  int st,ed,i;
  double temp;
  
  pflag = opts[1];
  sflag = opts[0];


  if(pflag == 1 && sflag == 0) {

    fprintf(stdout,"Enter the YEAR: ");
    fscanf(stdin,"%d",&st);
    if (st < 1990 || st > 2020)
      err_exit("No data for this YEAR");

    fprintf(stdout,"Enter the COOLING INDEX: ");
    fscanf(stdin,"%lf",&temp);
    if(temp < 0.)
      err_exit("COOLING INDEX is not proper");

    init->start_year = st;
    init->end_year = st;
    init->size = 1;
    init->cooling_index[0] = temp;

  } else if(sflag == 1 && pflag == 0){

    fprintf(stdout,"Enter the Start YEAR: \n");
    fscanf(stdin,"%d",&st);
    if (st < 1990 || st > 2020)
      err_exit("No data for this YEAR");

    fprintf(stdout,"Enter the Start YEAR: \n");
    fscanf(stdin,"%d",&ed);
    if (st > ed)
      err_exit("END YEAR is not proper");

    init->size = 11;
    init->start_year = st;
    init->end_year = ed;


    for (i = 0; i < 11; i++)
      
      init->cooling_index[i] = 1. + 0.1*i;
   
  } else {
    
    err_exit("Wrong Usage Initialization.");

  }
   
}



void Root_Bcat_data(CONF_t *init, int root){

  MPI_Datatype mystruct;
  int blocklens[4];
  MPI_Aint indices[4];
  MPI_Datatype old_types[4];


  /*
   * value of each type
   */

  for (int i = 0; i < 3; i++){
    blocklens[i] = 1;
    old_types[i] = MPI_INT;
  }
  blocklens[3] = 11;
  old_types[3] = MPI_DOUBLE;

  /*
   * location of each element
   */
  MPI_Get_address( &init->start_year, &indices[0] );
  MPI_Get_address( &init->end_year, &indices[1] );
  MPI_Get_address( &init->size, &indices[2] );
  MPI_Get_address( init->cooling_index, &indices[3] );

  /* Make relative */
  for (int i = 1; i < 4; i++)
    indices[i] = indices[i] - indices[0];
  indices[0] = 0;

  MPI_Type_create_struct( 4, blocklens, indices, old_types, &mystruct );
  MPI_Type_commit( &mystruct );

  MPI_Bcast(init, 1, mystruct, root, MPI_COMM_WORLD );
  MPI_Type_free( &mystruct );

}


