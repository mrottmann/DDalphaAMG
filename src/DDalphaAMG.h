/*
 * Copyright (C) 2016, Simone Bacchio.
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
 
#ifndef DDalphaAMG_HEADER
  #define DDalphaAMG_HEADER
  #define MAX_MG_LEVELS 4

  #include "mpi.h"

  /*
   *   Welcome user of the DDalphaAMG library!
   *
   *   In the following we will give the main instructions useful for including the library
   * in your code. The procedure should be quite straightforward and we suggest an accurate
   * reading of this header and of the user documentation in doc/.
   *
   *   For further information or remarks please contact 
   *   
   *     Simone Bacchio at <s.bacchio@hpc-leap.eu>.
   *
   *   The library makes use the following three structures,
   */

  typedef struct init        DDalphaAMG_init;
  typedef struct parameters  DDalphaAMG_parameters;
  typedef struct status {
    int success;
    double info;
    double time;
    double coarse_time;
    int iter_count;
    int coarse_iter_count; } DDalphaAMG_status;

  /*
   * where all the details of variables in init and parameters are given after the function
   * declarations coming in the following. Briefly, the first is the minimal set of parameters
   * required for the initialization; and the second is the full set of parameters tunable
   * by the user. The third structure, status, is required by any function and it returns some
   * details on the execution status. For simplicity's sake, you can print the output in
   * message for having a full view of the informations. 
   *
   *   Looking now to the functions, first of all the software has to be initialized calling
   */

  /**
   ** Initialize and return the default set of parameters.
   **  mg_status.success = number of levels initialized
   **/
  void DDalphaAMG_initialize( DDalphaAMG_init *mg_init, DDalphaAMG_parameters *mg_params, DDalphaAMG_status *mg_status );

  /*
   * which requires an appropiate set of input in mg_init and it returns the default set of
   * parameters in mg_params; mg_status.success is equal to the number of level that has been
   * initialized. The full set of parameter can be always retrieved by a call of
   */

  /**
   ** Optional - Get parameters.
   **  mg_status.success = 1-> same paramameters of given mg_params; 0-> some params cheanged;
   **/
  void DDalphaAMG_get_parameters( DDalphaAMG_parameters *mg_params );

  /*
   * and mg_status.success==0 alerts that some parameters has been changed wrt the previous 
   * set given to the function (in message the full set is given, and the same is done in
   * DDalphaAMG_init). 
   *
   *   Before proceeding with the setup, you may want to change the default parameters calling
   */

  /**
   ** Optional - Update parameters.
   **  mg_status.success = 0 -> some parameter not accepted + warning, 1 -> parameters updated,
   **                      2 -> operators updated, 3 -> setup run again
   **/
  void DDalphaAMG_update_parameters( DDalphaAMG_parameters *mg_params, DDalphaAMG_status *mg_status );

  /*
   * and be aware that some parameters can be changed only at this point; as described later
   * in the parameters structure. Then, the configuration has to be loaded giving to the
   * function
   */

  /**
   ** Set the configuration. 
   **  mg_status.success = 0 -> configuration not loaded, 1 -> configuration loaded/updated
   **  mg_status.info = plaquette
   **/
  void DDalphaAMG_set_configuration( double *gauge_field, DDalphaAMG_status *mg_status );

  /*
   * the complete gauge field stored in an array of double with length ( 18 x 4 x V ), where
   * V is the number of sites in the volume, 4 is the number matrices per site and 18 is the
   * size of the matrices (complex, 3x3). The matrices MUST stored in row major format,
   * alternating the real and immaginary components. The order of the sites and links is not
   * fixed and can be adjusted mapping the indeces with the function conf_index_fct(...) in
   * the init structure. A good way to check that everything has been done correctly is to 
   * compare the average plaquette stored in mg_status.info with the desired one.
   *
   *   Finally, the setup can be run using
   */

  /**
   ** Setup. If called again, it destroys the previous setup.
   **  mg_status.success = number of setup iters performed
   **/
  void DDalphaAMG_setup( DDalphaAMG_status *mg_status );

  /*
   * and more setup iterations can be done calling
   */

  /**
   ** Optional - More setup iterations.
   **  mg_status.success = TOTAL number setup iters (since DDalphaAMG_setup has been called)
   **/
  void DDalphaAMG_update_setup( int iterations, DDalphaAMG_status *mg_status );

  /*
   *   Now the method is ready for computing the invers of Wilson Dirac operator and several
   * functions are given for different purposes:
   */

  /**
   ** Optional - Invert the operator against vector_in at required tolerance, tol:
   **   vector_out = D^{-1} * vector_in.
   **  mg_status.success = 0: not converged, 1: converged
   **  mg_status.info = final residual
   **/
  void DDalphaAMG_solve( double *vector_out, double *vector_in,
                         double tol, DDalphaAMG_status *mg_status );
  
  void DDalphaAMG_solve_doublet( double *vector1_out, double *vector1_in, 
                                 double *vector2_out, double *vector2_in,
                                 double tol, DDalphaAMG_status *mg_status );

  /**
   ** Optional - Solve squared operator performing two inversions: 
   **   vector_out = D_{d/u}^{-1} * \Gamma_5 * D_{u/d}^{-1} * \Gamma_5 * vector_in.
   **  mg_status.success = 0: not converged, 1: converged
   **  mg_status.info = final residual
   **/
  void DDalphaAMG_solve_squared( double *vector_out, double *vector_in,
                                 double tol, DDalphaAMG_status *mg_status );
  
  void DDalphaAMG_solve_doublet_squared( double *vector1_out, double *vector1_in, 
                                         double *vector2_out, double *vector2_in,
                                         double tol, DDalphaAMG_status *mg_status );

  /**
   ** Optional - Solve squared operator against the odd compoments performing two inversions: 
   **   vector_out = D_{d/u}^{-1} * \Gamma_5 * P_{odd} * D_{u/d}^{-1} * \Gamma_5 * vector_in.
   **  mg_status.success = 0: not converged, 1: converged
   **  mg_status.info = final residual
   **/
  void DDalphaAMG_solve_squared_odd( double *vector_out, double *vector_in,
                                    double tol, DDalphaAMG_status *mg_status );
  
  void DDalphaAMG_solve_doublet_squared_odd( double *vector1_out, double *vector1_in, 
                                             double *vector2_out, double *vector2_in,
                                             double tol, DDalphaAMG_status *mg_status );

  /**
   ** Optional - Solve squared operator against the even compoments performing two inversions:
   **   vector_out = D_{d/u}^{-1} * \Gamma_5 * P_{even} * D_{u/d}^{-1} * \Gamma_5 * vector_in.
   **  mg_status.success = 0: not converged, 1: converged
   **  mg_status.info = final residual
   **/
  void DDalphaAMG_solve_squared_even( double *vector_out, double *vector_in,
                                     double tol, DDalphaAMG_status *mg_status );
  
  void DDalphaAMG_solve_doublet_squared_even( double *vector1_out, double *vector1_in, 
                                              double *vector2_out, double *vector2_in,
                                              double tol, DDalphaAMG_status *mg_status );

  /**
   ** Optional - Apply the operator:
   **   vector_out = D * vector_in.
   **  mg_status.success = 1
   **/
  void DDalphaAMG_apply_operator( double *vector_out, double *vector_in,
                                  DDalphaAMG_status *mg_status );

  void DDalphaAMG_apply_operator_doublet( double *vector1_out, double *vector1_in,
                                          double *vector2_out, double *vector2_in, DDalphaAMG_status *mg_status );

  /**
   ** Optional - Apply a preconditioner step:
   **   vector_out = D_c^{-1} * vector_in.
   **  mg_status.success = 1
   **  mg_status.info = residual after preconditioning
   **/
  void DDalphaAMG_preconditioner( double *vector_out, double *vector_in,
                                  DDalphaAMG_status *mg_status );
  void DDalphaAMG_preconditioner_doublet( double *vector1_out, double *vector1_in,
                                          double *vector2_out, double *vector2_in, DDalphaAMG_status *mg_status );

  /*
   *  Concluding the following functions have to be call for freeing the memory and finalizing
   * the software. 
   */

  /**
   ** Free.
   **/
  void DDalphaAMG_free( void );

  /**
   ** Free if not freed and finalize.
   **/
  void DDalphaAMG_finalize( void );

  /*
   *  Some extra functions are
   */

  /**
   ** Extra - Change sign of mu and update all the operators.
   **/
  void DDalphaAMG_change_mu_sign( DDalphaAMG_status *mg_status );

  /**
   ** Extra - Returns the MPI communicator used by the library. 
   **/
  MPI_Comm DDalphaAMG_get_communicator( void );

  /**
   ** Extra - Read configuration and read/write vector
   **  Formats: 0 -> DDalphaAMG standard, 1 -> LIME format
   **  Note: configurations and vector read with this function do NOT require reordering,
   **    -> mg_params.conf_index_fct = NULL, mg_params.vector_index_fct = NULL;
   **/
  void DDalphaAMG_read_configuration( double *gauge_field, char *filename, int format,
                                      DDalphaAMG_status *mg_status );
  void DDalphaAMG_read_vector( double *vector_in, char *filename, int format,
                               DDalphaAMG_status *mg_status );
  void DDalphaAMG_write_vector( double *vector_out, char *filename, int format,
                                DDalphaAMG_status *mg_status );

  /**
   ** Extra - Define vector with constant or random components
   **/
  void DDalphaAMG_define_vector_const( double *vector, double re, double im );
  void DDalphaAMG_define_vector_rand( double *vector );

  /**
   ** Extra - Vector algebra: norm, saxpy
   **/
  double DDalphaAMG_vector_norm( double *vector );
  void DDalphaAMG_vector_saxpy( double *vector_out, double a, double *x, double *y );

  /**
   ** Extra - Full test routinle
   **/
  void DDalphaAMG_test_routine( DDalphaAMG_status *mg_status );

  /*
   *  And here the descriptions of variables in the stracture init and following stracture
   * parameters. Enjoy! 
   */
  struct init {
    
    /**
     ** MPI commuticator. It can be
     **
     ** (a) MPI_COMM_WORLD
     ** (b) A split of MPI_COMM_WORLD
     ** (c) A cartesian communicator with 4 dims, number of processes in each dimension equals
     **     to procs[4] and proper bondary conditions
     **/  
    MPI_Comm comm_cart;

    /**
     ** MPI cartesian functions.
     **   Custom functions can be constracted and they have to work as the respective
     ** MPI_Cart_* functions. In case of NULL, the standard MPI function will be used. In case
     ** of standard MPI functions and MPI communicator equals to 1 or 2, a cartesian comm will
     ** be defined and returned in mg_params.comm_cart.
     **/
    int (*Cart_rank)(MPI_Comm comm, const int coords[], int *rank);
    int (*Cart_coords)(MPI_Comm comm, int rank, int maxdims, int coords[]);
  
    /**
     **  Size of lattice and number of processes in each direction.
     **  The directions order is T, Z, Y, X.
     **/
    int global_lattice[4];
    int procs[4];

    /**
     **  Size of aggregation blocks. 
     **    It can be NULL and the initialization will produce the recommended size.
     **  The order of directions is T, Z, Y, X.
     **/
    int block_lattice[4];
    
    /**
     **  Boundary conditions:
     **    0: periodic,
     **    1: anti-periodic,
     **    2: twisted -> a phase proportional to M_PI * theta[i] will be all links in dir. i
     **  The order of directions is T, Z, Y, X.
     **/
    int bc;
    double theta[4];

    /**
     **  Number of levels for multigrid:
     **    from 1 (no MG) to 4 (max number of levels)
     **/
    int number_of_levels;
    
    /**
     **  Number of openmp threads:
     **    from 1 to omp_get_num_threads()
     **/
    int number_openmp_threads;
    
    /**
     **  Operator parameters
     **/
    double kappa;
    double mu;
    double csw;

    /**
     **  Parameter file-name.
     **   If NULL no file will be read.
     **   If a file is given, all the parameters except the previous ones are read.
     **/
    char * init_file;

    /**
     **  Random seed generator (srand is used)
     **   If NULL random numbers will be automatically produced.
     **   If defined, it has to be a vector of lenght MPI tasks. 
     **    The task with rank n will take the nth seed.
     **/
    unsigned int * rnd_seeds;

  };

  struct parameters {
                
    /**
     ** Method used by the solver:                                
     **         
     **  method = -1 - pure CGN (no AMG)         
     **  method =  0 - pure GMRES (no AMG)       
     **  method =  1 - FGMRES + additive Schwarz             
     **  method =  2 - FGMRES + red-black Schwarz                 
     **  method =  3 - FGMRES + 16 color Schwarz      
     **  method =  4 - FGMRES + GMRES                 
     **  method =  5 - FGMRES + biCGstab (no AMG) 
     **/  
    int method;

    /**
     ** Interpolation used by the multigrid:
     **
     **  interpolation = 0 - no interpolation
     **  interpolation = 1 - successive inv-iter
     **  interpolation = 2 - f-cycle inv-iter
     **  interpolation = 3 - f-cycle fgmres+amg-inv-it
     **/
    int interpolation;

    /**
     ** Mixed precision solver:                               
     **  
     **  0 : double precision on fine and coarse levels
     **  1 : double precision on fine level, single precision on coarse levels
     **  2 : mixed precision on fine level, single precision on coarse levels
     **/  
    int mixed_precision;
    
    /**
     ** Aggregation block lattice size.
     **
     **  The same size is used for aggregation and red-black Schwarz.
     **  For details on allowed sizes, look at the doc/.         
     **/
    int block_lattice[MAX_MG_LEVELS][4];

    /**
     ** Number of basis vectors used used between level i and i+1.
     **  Default 20 -> 28 -> 40
     **  If the numbers are DECREASED after the setup, you do NOT NEED to run the setup again
     **/
    int mg_basis_vectors[MAX_MG_LEVELS-1];
  
    /**
     ** Number of setup iterations on each level i.
     **  Default 6 -> 3 -> 2 -> 2
     **/
    int setup_iterations[MAX_MG_LEVELS];

    /**
     ** Smoother iterations.
     **  Default 4
     **/
    int smoother_iterations;

    /**
     ** Tolerance used for inversion during k-cycle.
     **  Default 1e-1
     **/
    double kcycle_tolerance;
    
    /**
     ** Tolerance used on the coarsest level for inversion.
     **  Default 1e-1
     **/
    double coarse_tolerance;
        
    /**
     ** Hopping parameter
     **/
    double kappa;

    /**
     ** Twisted mass parameter and shifts on even/odd sites
     **/
    double mu;
    double mu_odd_shift;
    double mu_even_shift;

    /**
     ** Twisted mass factor for the preconditioner on each level, l.
     **  Default 6 on the coarsest level 
     ** 
     **   -> mu_o[l] = (mu + mu_odd_shift)  * mu_factor[l]
     **   -> mu_e[l] = (mu + mu_even_shift) * mu_factor[l]
     **/
    double mu_factor[MAX_MG_LEVELS];  

    /**
     ** Twisted mass doublet parameter and shifts on even/odd sites
     **/
    double epsbar;
    double epsbar_ig5_odd_shift;
    double epsbar_ig5_even_shift;

    /**
     ** Twisted mass doublet factor for the preconditioner on each level, l.
     **  Default 6 on the coarsest level 
     ** 
     **   -> epsbar_o[l] = ( epsbar + i * gamma_5 * epsbar_ig5_odd_shift ) * epsbar_factor[l]
     **   -> epsbar_e[l] = ( epsbar + i * gamma_5 * epsbar_ig5_even_shift ) * epsbar_factor[l]
     **/
    double epsbar_factor[MAX_MG_LEVELS];  
    
    /**
     ** Function returning the index of a element at the corresponding
     **    position (t,z,y,x are local position w.r.t the process ).
     **
     ** In DDalphaAMG format:
     **   directions ordering: T, Z, Y, X
     **   link ordering: T, Z, Y, X
     **   SU3 matrices stored in row major format
     **
     ** If their value is NULL, the conf/vector is assumed in DDalphaAMG format.
     ** If the included library lime_io is used, their value has to be NULL
     **/  
    int (*conf_index_fct)(int t, int z, int y, int x, int mu);
    int (*vector_index_fct)(int t, int z, int y, int x );

    /**
     ** Printing level:
     **  -1: silent (errors or warnings)
     **   0: minimal //default
     **   1: verbose
     **/
    int print;
 
  };
 
#endif
