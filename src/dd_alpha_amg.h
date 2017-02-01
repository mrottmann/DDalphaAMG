/*
 * Copyright (C) 2016, Matthias Rottmann, Artur Strebel, Simon Heybrock, Simone Bacchio, Bjoern Leder, Issaku Kanamori.
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

#include "dd_alpha_amg_parameters.h"
  
  #define STRINGLENGTH 500

  typedef struct {
    char param_file_path[STRINGLENGTH];
    int (*conf_index_fct)(int t, int z, int y, int x, int mu);
    int (*vector_index_fct)(int t, int z, int y, int x);
    int (*global_time)(int t);
    int bc; /* bc: 0 dirichlet, 1 periodic, 2 anti-periodic */
    double m0;
    double csw;
    double setup_m0;
    struct dd_alpha_amg_parameters amg_params;
  } dd_alpha_amg_par;

  /* SU3 matrices stored in row major format */
  void dd_alpha_amg_init( dd_alpha_amg_par p );
  void dd_alpha_amg_init_external_threading( dd_alpha_amg_par p, int n_core, int n_thread );

  double* dd_alpha_amg_get_gauge_pointer( void );
  double* dd_alpha_amg_get_clover_pointer( void );
  /* Notifies MG that the (gauge+clover) fields were/are being updated externally
   * (via the pointers obtained by dd_alpha_amg_get_gauge_pointer() and
   * dd_alpha_amg_get_clover_pointer()).
   * This is required to allow MG to determine when the MG setup vectors should
   * be updated or discarded.
   * When the functions dd_alpha_amg_put_conf() and dd_alpha_amg_put_clover() are used this is
   * done automatically, and calling this function is *not* necessary.
   */
  void dd_alpha_amg_fields_updated( void );
    
  double dd_alpha_amg_set_conf( double *gauge_field );
  /* does not use conf_index_fct, assumes that caller
   * provides gauge_field in proper data layout */

  void dd_alpha_amg_update_parameters( const struct dd_alpha_amg_parameters *amg_params );
  
  void dd_alpha_amg_setup( int iterations, int *status );
  void dd_alpha_amg_setup_external_threading( int iterations, int *status,
      int core, int thread, void *thread_barrier_data, void (*thread_barrier)(void *, int) );
  
  void dd_alpha_amg_setup_update( int iterations, int *status );
  void dd_alpha_amg_setup_update_external_threading( int iterations, int *status,
      int core, int thread, void *thread_barrier_data, void (*thread_barrier)(void *, int) );
  
  double dd_alpha_amg_wilson_solve( double *vector_out, double *vector_in, double tol, double scale_even, double scale_odd, int *status );
  /* does not use vector_index_fct, assumes that caller
   * provides vector_out and vector_in in proper data layout
   * CAUTION: this interface is not designed for vector_in = 0
   *          users are supposed to take care of such cases themselves */
  void dd_alpha_amg_preconditioner( double *vector_out, double *vector_in, double scale_even, double scale_odd, int *status );
  /* threads spawned by user, each thread can call this function,
   * for now we assume that barrier among cores can simply use omp barrier,
   * barrier among hyperthreads is passed in extra arguments */
  void dd_alpha_amg_preconditioner_external_threading( double *vector_out, double *vector_in, int *status,
      int core, int thread, void *thread_barrier_data, void (*thread_barrier)(void *, int));
  
  void dd_alpha_amg_free( void );

#endif
