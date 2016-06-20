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

#ifndef DDalphaAMG_INTERFACE
  #define DDalphaAMG_INTERFACE

  #define STRINGLENGTH 500
  #define MAX_MG_LEVELS 4

  typedef struct DDalphaAMG_parameters {

    int number_of_levels;

    int global_lattice[MAX_MG_LEVELS][4];
    int local_lattice[MAX_MG_LEVELS][4];
    int block_lattice[MAX_MG_LEVELS][4];

    int mg_basis_vectors[MAX_MG_LEVELS];
    int setup_iterations[MAX_MG_LEVELS];
    int discard_setup_after;
    int update_setup_iterations[MAX_MG_LEVELS];
    int update_setup_after;

    int post_smooth_iterations[MAX_MG_LEVELS];
    int post_smooth_block_iterations[MAX_MG_LEVELS];

    int coarse_grid_iterations;
    int coarse_grid_maximum_number_of_restarts;
    double coarse_grid_tolerance;

    double solver_mass;
    double setup_mass;
    double c_sw;

  } DDalphaAMG_parameters;


  typedef struct {
    char param_file_path[STRINGLENGTH];
    int (*conf_index_fct)(int t, int z, int y, int x, int mu);
    int (*vector_index_fct)(int t, int z, int y, int x);
    int (*global_time)(int t);
    int bc; /* bc: 0 dirichlet, 1 periodic, 2 anti-periodic */
    double m0;
    double csw;
    double setup_m0;
    struct DDalphaAMG_parameters amg_params;
  } DDalphaAMG_par;

  /* SU3 matrices stored in row major format */
  void DDalphaAMG_init( DDalphaAMG_par p );
  void DDalphaAMG_init_external_threading( DDalphaAMG_par p, int n_core, int n_thread );

  double* DDalphaAMG_get_gauge_pointer( void );
  double* DDalphaAMG_get_clover_pointer( void );
  /* Notifies MG that the (gauge+clover) fields were/are being updated externally
   * (via the pointers obtained by DDalphaAMG_get_gauge_pointer() and
   * DDalphaAMG_get_clover_pointer()).
   * This is required to allow MG to determine when the MG setup vectors should
   * be updated or discarded.
   * When the functions DDalphaAMG_put_conf() and DDalphaAMG_put_clover() are used this is
   * done automatically, and calling this function is *not* necessary.
   */
  void DDalphaAMG_fields_updated( void );
    
  double DDalphaAMG_set_conf( double *gauge_field );
  /* does not use conf_index_fct, assumes that caller
   * provides gauge_field in proper data layout */

  void DDalphaAMG_update_parameters( const struct DDalphaAMG_parameters *amg_params );
  
  void DDalphaAMG_setup( int iterations, int *status );
  void DDalphaAMG_setup_external_threading( int iterations, int *status,
      int core, int thread, void *thread_barrier_data, void (*thread_barrier)(void *, int) );
  
  void DDalphaAMG_setup_update( int iterations, int *status );
  void DDalphaAMG_setup_update_external_threading( int iterations, int *status,
      int core, int thread, void *thread_barrier_data, void (*thread_barrier)(void *, int) );
  
  double DDalphaAMG_wilson_solve( double *vector_out, double *vector_in, double tol, double scale_even, double scale_odd, int *status );
  /* does not use vector_index_fct, assumes that caller
   * provides vector_out and vector_in in proper data layout */
  void DDalphaAMG_preconditioner( double *vector_out, double *vector_in, double scale_even, double scale_odd, int *status );
  /* threads spawned by user, each thread can call this function,
   * for now we assume that barrier among cores can simply use omp barrier,
   * barrier among hyperthreads is passed in extra arguments */
  void DDalphaAMG_preconditioner_external_threading( double *vector_out, double *vector_in, int *status,
      int core, int thread, void *thread_barrier_data, void (*thread_barrier)(void *, int));
  
  void DDalphaAMG_free( void );

#endif
