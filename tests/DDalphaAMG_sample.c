/*
 * Copyright (C) 2016, Matthias Rottmann, Artur Strebel, Simon Heybrock, Simone Bacchio, Bjoern Leder.
 * 
 * This file is part of the DDalphaAMG solver library.
 * 
 * The DDalphaAMG solver library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * The DDalphaAMG solver library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * 
 * You should have received a copy of the GNU General Public License
 * along with the DDalphaAMG solver library. If not, see http://www.gnu.org/licenses/.
 * 
 */

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <mpi.h>
#include <omp.h>

//library to include
#include "DDalphaAMG.h"

enum { T, Z, Y, X };
#define printf0(...) do{if(rank==0)printf(__VA_ARGS__);}while(0)

int rank, max_setup=0, conf_format;
float residual=1e-9;
MPI_Comm comm_cart;
DDalphaAMG_init init;
DDalphaAMG_parameters params;
DDalphaAMG_status status;
char * conf_file = "../conf/8x8x8x8b6.0000id3n1";
char * options = "c:f:i:L:p:B:l:k:w:u:t:hr:C:K:V:m:s:S:v";
/*
 * Setting standard values for DDalphaAMG_init
 */
void standard_init() {
  init.global_lattice[T] = 8;
  init.global_lattice[Z] = 8;
  init.global_lattice[Y] = 8;
  init.global_lattice[X] = 8;

  init.kappa = 0.142857143;
  init.mu = 0.1;
  init.csw = 1;

  init.bc = 1;

  init.number_of_levels = 3;

  init.procs[T] = 2;
  init.procs[Z] = 1;
  init.procs[Y] = 1;
  init.procs[X] = 1;

  init.number_openmp_threads = 2;

  init.block_lattice[T]=2;
  init.block_lattice[Z]=2;
  init.block_lattice[Y]=2;
  init.block_lattice[X]=2;

  init.comm_cart = MPI_COMM_WORLD;
  init.Cart_coords = NULL;
  init.Cart_rank = NULL;
  init.init_file = NULL;
}

/*
 * Printing implemented parameters
 */
void help( char * arg0 ) {
  static int printed = 0;
  if(!printed) {
    printf0("\n\n");
    printf0("Usage: %s -c <conf> [<option(s)>]\n", arg0);
    printf0("   -c PATH      Configuration to load\n");
    printf0("   -f #         Configuration format (0 -> DDalphaAMG, 1 -> Lime)\n");
    printf0("   -i PATH      Input file (optional)\n");
    printf0("   -L T Y X Z   Lattice size in each direction\n");
    printf0("   -p T Y X Z   Processors in each direction\n");
    printf0("   -B T Y X Z   Block size in each direction on first level.\n");
    printf0("   -k #         kappa for the configuration\n");
    printf0("   -w #         c_sw for the configuration\n");
    printf0("   -u #         mu for the configuration\n");
    printf0("   -t #         Number of OpenMp threads\n");
    printf0("   -r #         Relative residual\n");
    printf0("   -K #         K-cycle tolerance\n");
    printf0("   -C #         Tolerance on coarsest grid\n");
    printf0("   -l #         Number of levels, l (from 1 to 4)\n");
    printf0("   -V 1 [2] [3] Basis vectors between each level (l-1)\n");
    printf0("   -m 2 [3] [4] Factor for mu on coarse levels\n");
    printf0("   -s 1 [2] [3] Setup iterations on each level (l-1)\n");
    printf0("   -S 1         Test additional setup iterations\n");
    printf0("   -v           Verbose\n");
  }
  printf0("\n\n");
  printed++;
} 

/*
 * Printing implemented parameters
 */
void read_init_arg(int argc, char *argv[] ) {

  int opt, p=0, fail=0, mu;
  optind = 0;
  while ((opt = getopt(argc, argv, options)) != -1) {
    switch (opt) {
    case 'c':
      conf_file = optarg;
      break;
    case 'f':
      conf_format = atoi(optarg);
      break;
    case 'i':
      init.init_file = optarg;
      break;
    case 'L':
      optind--;
      mu=0;
      for ( ; optind < argc && *argv[optind] != '-'; optind++){
	if(mu > 3) {
	  printf0("Error: too many arguments in -L.\n");
	  p++;
	  fail++;
	  break;
	}
	init.global_lattice[mu] = atoi(argv[optind]);
	mu++;        
      }
      if(mu < 4) {
	printf0("Warning: too few arguments in -L.\n");
	p++;
	fail++;
      }
      break;
    case 'p':
      optind--;
      mu=0;
      for ( ; optind < argc && *argv[optind] != '-'; optind++){
	if(mu > 3) {
	  printf0("Error: too many arguments in -p.\n");
	  p++;
	  fail++;
	  break;
	}
	init.procs[mu] = atoi(argv[optind]);
	mu++;        
      }
      if(mu < 4) {
	printf0("Warning: too few arguments in -p.\n");
	p++;
      }
      break;
    case 'B':
      optind--;
      mu=0;
      for ( ; optind < argc && *argv[optind] != '-'; optind++){
	if(mu > 3) {
	  printf0("Error: too many arguments in -B.\n");
	  p++;
	  fail++;
	  break;
	}
	init.block_lattice[mu] = atoi(argv[optind]);
	mu++;        
      }
      if(mu < 4) {
	printf0("Warning: too few arguments in -B.\n");
	p++;
      }
      break;
    case 'l':
      init.number_of_levels = atoi(optarg);
      break;
    case 'k':
      init.kappa = atof(optarg);
      break;
    case 'w':
      init.csw = atof(optarg);
      break;
    case 'u':
      init.mu = atof(optarg);
      break;
    case 't':
      init.number_openmp_threads = atoi(optarg);
      break;
    case '?':
    case 'h':
      p++;
      break;
    default: 
      break;
    }
  }

  if(conf_file == NULL) {
    printf0("Error: configuration file is missing (use -c PATH).\n");
    p++;
    fail++;
  }
  
  if(p) {
    help(argv[0]);
    MPI_Abort(MPI_COMM_WORLD,0);
  }
  if(fail)
    MPI_Abort(MPI_COMM_WORLD,0);
}

void read_params_arg(int argc, char *argv[] ) {

  int opt, p=0, fail=0, mu;
  optind = 0;
  while ((opt = getopt(argc, argv, options)) != -1) {
    switch (opt) {
    case 'r':
      residual = atof(optarg);
      break;
    case 'K':
      params.kcycle_tolerance = atof(optarg);
      break;
    case 'C':
      params.coarse_tolerance = atof(optarg);
      break;
    case 'V':
      optind--;
      mu=0;
      for ( ; optind < argc && *argv[optind] != '-'; optind++){
	if(mu > 2) {
	  printf0("Error: too many arguments in -V.\n");
	  p++;
	  fail++;
	  break;
	}
	params.mg_basis_vectors[mu] = atoi(argv[optind]);
	mu++;        
      }
      break;
    case 'm':
      optind--;
      params.mu_factor[0]=1;
      mu=1;
      for ( ; optind < argc && *argv[optind] != '-'; optind++){
	if(mu > 3) {
	  printf0("Error: too many arguments in -m.\n");
	  p++;
	  fail++;
	  break;
	}
	params.mu_factor[mu] = atof(argv[optind]);
	mu++;        
      }
      break;
    case 's':
      optind--;
      mu=0;
      for ( ; optind < argc && *argv[optind] != '-'; optind++){
	if(mu > 2) {
	  printf0("Error: too many arguments in -s.\n");
	  p++;
	  fail++;
	  break;
	}
	params.setup_iterations[mu] = atoi(argv[optind]);
	mu++;        
      }
      break;
    case 'S':
      max_setup=atoi(optarg);
      break;
    case 'v':
      params.print = 1;
      break;
    case '?':
    default: 
      break;
    }
  }

  if(conf_file == NULL) {
    printf0("Error: configuration file is missing (use -c PATH).\n");
    p++;
    fail++;
  }
  
  if(p) {
    help(argv[0]);
    MPI_Abort(MPI_COMM_WORLD,0);
  }
  if(fail)
    MPI_Abort(MPI_COMM_WORLD,0);
}

int main( int argc, char *argv[] ) {

  
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  standard_init();
  read_init_arg(argc, argv);
  
  printf0("Running initialization...\n");
  DDalphaAMG_initialize( &init, &params, &status );
  printf0("Initialized %d levels in %.2f sec\n", status.success, status.time);

  read_params_arg(argc, argv);

  comm_cart =  DDalphaAMG_get_communicator();
  MPI_Comm_rank( comm_cart, &rank );

  printf0("Running updating\n");
  DDalphaAMG_update_parameters( &params, &status );
  if (status.success)
    printf0("Updating time %.2f sec\n", status.time);

  /*
   * Reading the configuration. In plaq, it returns the plaquette value
   *  if provided in the configuration header.
   */
  double *gauge_field;
  int vol = init.global_lattice[T] * init.global_lattice[X] * init.global_lattice[Y] *
    init.global_lattice[Z] / init.procs[T] / init.procs[X] / init.procs[Y] / init.procs[Z];
  gauge_field = (double *) malloc(18*4*vol*sizeof(double));

  printf0("Reading config.\n");
  DDalphaAMG_read_configuration( gauge_field, conf_file, conf_format, &status );
  printf0("Reading configuration time %.2f sec\n", status.time);
  printf0("Desired plaquette %.13lf\n", status.info);
      
  printf0("Setting config.\n");
  DDalphaAMG_set_configuration( gauge_field, &status );
  printf0("Setting configuration time %.2f sec\n", status.time);
  printf0("Computed plaquette %.13lf\n", status.info);
  

  printf0("Running setup\n");
  DDalphaAMG_setup( &status );
  printf0("Run %d setup iterations in %.2f sec (%.1f %% on coarse grid)\n", status.success,
	  status.time, 100.*(status.coarse_time/status.time));
  printf0("Total iterations on fine grid %d\n", status.iter_count);
  printf0("Total iterations on coarse grids %d\n", status.coarse_iter_count);
 
  /*
   * Defining the vector randomly.
   */
  double *vector_in, *vector_out;
  vector_in = (double *) malloc(24*vol*sizeof(double));
  vector_out = (double *) malloc(24*vol*sizeof(double));
  
  DDalphaAMG_define_vector_rand(vector_in);
  
  do {
    printf0("Running solver for up propagator\n");
    DDalphaAMG_solve( vector_out, vector_in, residual, &status );
    if (status.success)
      printf0("Converged with final relative residual %e\n", status.info);
    else
      printf0("ERROR: not converged with final relative residual %e\n", status.info);
    printf0("Solving time %.2f sec (%.1f %% on coarse grid)\n", status.time,
	    100.*(status.coarse_time/status.time));
    printf0("Total iterations on fine grid %d\n", status.iter_count);
    printf0("Total iterations on coarse grids %d\n", status.coarse_iter_count);
     
    
    if (params.mu != 0) {
      // Changing mu to down.
      DDalphaAMG_change_mu_sign(&status);
      printf0("Changing mu time %.2f sec \n", status.time);
      
      printf0("Running solver for down propagator\n");
      DDalphaAMG_solve( vector_out, vector_in, residual, &status );
      if (status.success) {
      printf0("Solving time %.2f sec (%.1f %% on coarse grid)\n", status.time,
	      100.*(status.coarse_time/status.time));
      printf0("Total iterations on fine grid %d\n", status.iter_count);
      printf0("Total iterations on coarse grids %d\n", status.coarse_iter_count);
      } 
      else 
	printf0("ERROR: not converged\n");
      
      DDalphaAMG_change_mu_sign(&status);
      printf0("Changing mu time %.2f sec \n", status.time);
    }

    if( params.setup_iterations[0]<max_setup){
      if (params.mu != 0) {
	DDalphaAMG_change_mu_sign(&status);
	printf0("Changing mu time %.2f sec \n", status.time);
      }
      
      params.setup_iterations[0]++;
      printf0("Updating setup\n");
      DDalphaAMG_update_setup( 1, &status );
      printf0("Setup updated in %.2f sec (%.1f %% on coarse grid)\n",
	      status.time, 100.*(status.coarse_time/status.time));
      printf0("Total iterations on fine grid %d\n", status.iter_count);
      printf0("Total iterations on coarse grids %d\n", status.coarse_iter_count);
    }
    else
      break;
  }while(0);
  
  //  free(vector_in);
  // free(vector_out);
  //free(gauge_field);
  //  DDalphaAMG_finalize();
  MPI_Finalize();
}
