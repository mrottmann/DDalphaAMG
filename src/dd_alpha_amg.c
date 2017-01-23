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

#include "main.h"
#include "dd_alpha_amg.h"


#define NCORE 1
 
global_struct g;
static level_struct l;
static int (*get_global_time)(int t);
static struct common_thread_data *commonthreaddata;
static struct Thread **threading;
struct Thread *no_threading;

int get_thread_id( int core, int thread ) {
#ifdef FLAT_OMP
  return core+60*thread;
#else
  return threading[0]->n_thread*core+thread;
#endif
}

void update_threading_struct_with_barrier_info( int thread_id, void *thread_barrier_data, void (*thread_barrier)(void *, int) ) {
  threading[thread_id]->thread_barrier_data = thread_barrier_data;
  threading[thread_id]->thread_barrier      = thread_barrier;
}

void run_setup( struct Thread *threading ) {

  START_LOCKED_MASTER(threading)
  g.coarse_iter_count = 0;
  if ( g.setup_flag )
    method_free( &l );
  END_LOCKED_MASTER(threading)

  method_setup( NULL, &l, threading );
  START_LOCKED_MASTER(threading)
  g.setup_flag = 1;
  END_LOCKED_MASTER(threading)
  method_update( g.setup_iter[0], &l, threading );

  START_LOCKED_MASTER(threading)
  g.conf_flag = 0;
  END_LOCKED_MASTER(threading)
}

void run_setup_update( struct Thread *threading ) {

  if ( g.conf_flag == 1 ) {
    START_LOCKED_MASTER(threading)
    g.conf_flag = 0;
    END_LOCKED_MASTER(threading)
    if ( g.mixed_precision )
      re_setup_float( &l, threading );
    else
      re_setup_double( &l, threading );
  }

  START_LOCKED_MASTER(threading)
  g.coarse_iter_count = 0;
  END_LOCKED_MASTER(threading)
  method_update( g.setup_iter[0], &l, threading );
}

void run_dd_alpha_amg_setup_if_necessary( struct Thread *threading ) {
  if( g.mg_setup_status.gauge_updates_since_last_setup >= g.amg_params.discard_setup_after )
    run_setup( threading );
  else if( g.mg_setup_status.gauge_updates_since_last_setup_update >= g.amg_params.update_setup_after )
    run_setup_update( threading );

  if( g.mass_for_next_solve != l.dirac_shift )
    shift_update( (complex_double)g.mass_for_next_solve, &l, threading );
}

void dd_alpha_amg_init( dd_alpha_amg_par p ) {
  
  predefine_rank();
  l_init( &l );
  g_init( &l );
  g.conf_index_fct = p.conf_index_fct;
  g.vector_index_fct = p.vector_index_fct;
  g.bc = p.bc;
  get_global_time = p.global_time;
  lg_in( p.param_file_path, &l );
  g.csw = p.csw;
  g.setup_m0 = p.setup_m0;
  l.real_shift = p.m0;
  l.dirac_shift = l.real_shift;
  l.even_shift = l.dirac_shift;
  l.odd_shift = l.dirac_shift;
  cart_define( &l );
  data_layout_init( &l );
  operator_double_alloc( &(g.op_double), _ORDINARY, &l );
  operator_double_define( &(g.op_double), &l );
  MALLOC( g.odd_even_table, int, l.num_inner_lattice_sites );
  define_odd_even_table( &l );
  g.conf_flag = 0;
  g.setup_flag = 0;

  commonthreaddata = (struct common_thread_data *)malloc(sizeof(struct common_thread_data));
  init_common_thread_data(commonthreaddata);
  threading = (struct Thread **)malloc(sizeof(struct Thread *)*NCORE);
  for(int i=0; i<NCORE; i++) {
    threading[i] = (struct Thread *)malloc(sizeof(struct Thread));
  }
#pragma omp parallel num_threads(NCORE)
  {
    setup_threading(threading[omp_get_thread_num()], commonthreaddata, &l);
  }
  no_threading = (struct Thread *)malloc(sizeof(struct Thread));
  setup_no_threading(no_threading, &l);
}


void dd_alpha_amg_init_external_threading( dd_alpha_amg_par p, int n_core, int n_thread ) {

  predefine_rank();
  l_init( &l );
  g_init( &l );
  g.conf_index_fct = p.conf_index_fct;
  g.vector_index_fct = p.vector_index_fct;
  g.bc = p.bc;
  get_global_time = p.global_time;
  set_dd_alpha_amg_parameters( &(p.amg_params), &l );
  g.csw = p.csw;
  g.setup_m0 = p.setup_m0;
  l.real_shift = p.m0;
  l.dirac_shift = l.real_shift;
  l.even_shift = l.dirac_shift;
  l.odd_shift = l.dirac_shift;
  cart_define( &l );
  data_layout_init( &l );
  operator_double_alloc( &(g.op_double), _ORDINARY, &l );
  operator_double_define( &(g.op_double), &l );
  MALLOC( g.odd_even_table, int, l.num_inner_lattice_sites );
  define_odd_even_table( &l );
  g.conf_flag = 0;
  g.setup_flag = 0;

  commonthreaddata = (struct common_thread_data *)malloc(sizeof(struct common_thread_data));
  init_common_thread_data(commonthreaddata);
  threading = (struct Thread **)malloc(sizeof(struct Thread *)*n_core*n_thread);
  for(int i=0; i<n_core*n_thread; i++) {
    threading[i] = (struct Thread *)malloc(sizeof(struct Thread));
#ifdef FLAT_OMP
    setup_threading_external(threading[i], commonthreaddata, &l, n_core*n_thread, 1, i, 0);
#else
    setup_threading_external(threading[i], commonthreaddata, &l, n_core, n_thread, i/n_thread, i%n_thread);
#endif
  }
  no_threading = (struct Thread *)malloc(sizeof(struct Thread));
  setup_no_threading(no_threading, &l);
}


double *dd_alpha_amg_get_gauge_pointer() {
  return dirac_setup_get_gauge_pointer();
}
double *dd_alpha_amg_get_clover_pointer() {
  return dirac_setup_get_clover_pointer();
}
void dd_alpha_amg_fields_updated() {
  g.mg_setup_status.gauge_updates_since_last_setup++;
  g.mg_setup_status.gauge_updates_since_last_setup_update++;
}


double dd_alpha_amg_set_conf( double *gauge_field ) {
  
  int t, z, y, x, mu, i, j, k, *ll = l.local_lattice, ifail=0, tg;
  double *hopp = NULL,*clover = NULL;
  
  MALLOC( hopp, double, 2*3*l.inner_vector_size );
  
  if (g.bc>0) {
    for ( j=0, t=0; t<ll[T]; t++ )
      for ( z=0; z<ll[Z]; z++ )
        for ( y=0; y<ll[Y]; y++ )
          for ( x=0; x<ll[X]; x++ )
            for ( mu=0; mu<4; mu++ ) {
              i = g.conf_index_fct( t, z, y, x, mu );
              for ( k=0; k<18; k++, j++ )
                hopp[j] = gauge_field[i+k];
            }
  } else {
    MALLOC( clover, double, 2*3*l.inner_vector_size );
    for ( j=0, t=0; t<ll[T]; t++ ) {
      tg=get_global_time(t);
      for ( z=0; z<ll[Z]; z++ ) {
        for ( y=0; y<ll[Y]; y++ ) {
          for ( x=0; x<ll[X]; x++ ) {
            mu=0;
            if ((tg==0)||(tg>=g.process_grid[T]*ll[T]-2)) {
               i = g.conf_index_fct( t, z, y, x, mu );
               for ( k=0; k<18; k++, j++ ) {
                  hopp[j] = 0.0;
                  clover[j] = gauge_field[i+k];
                  if ((tg==g.process_grid[T]*ll[T]-1)&&(gauge_field[i+k]!=0.0))
                    ifail+=1;
               }
               mu=1;
            }
            for ( ; mu<4; mu++ ) {
              i = g.conf_index_fct( t, z, y, x, mu );
              for ( k=0; k<18; k++, j++ ) {
                hopp[j] = gauge_field[i+k];
                clover[j] = gauge_field[i+k];
              }
            }  
          }
        }
      }
    }
  }
  
  if (ifail)
    error0("Error in \"dd_alpha_amg_set_conf\": Gauge field does not fit expected boundary conditions.\n");
  
  if (g.bc>0)
    dirac_setup( (complex_double*) hopp, (complex_double*) hopp, &l );
  else
    dirac_setup( (complex_double*) hopp, (complex_double*) clover, &l );
  
  FREE( hopp, double, 2*3*l.inner_vector_size );
  if (g.bc==0)
    FREE( clover, double, 2*3*l.inner_vector_size );

  g.conf_flag = 1;
  return g.plaq;
}


void dd_alpha_amg_update_parameters( const struct dd_alpha_amg_parameters *amg_params ) {
  update_dd_alpha_amg_parameters( amg_params, &l );
}


void dd_alpha_amg_setup( int iterations, int *status ) {
  
  g.coarse_iter_count = 0;
  if ( g.setup_flag )
    method_free( &l );
#pragma omp parallel num_threads(threading[0]->n_core)
  {
    method_setup( NULL, &l, threading[omp_get_thread_num()] );
    START_LOCKED_MASTER(threading[omp_get_thread_num()])
    g.setup_flag = 1;
    END_LOCKED_MASTER(threading[omp_get_thread_num()])
    method_update( iterations, &l, threading[omp_get_thread_num()] );
  }
  status[1] = g.coarse_iter_count;
  status[0] = 1;
  g.conf_flag = 0;
}
void dd_alpha_amg_setup_external_threading( int iterations, int *status,
    int core, int thread, void *thread_barrier_data, void (*thread_barrier)(void *, int) ) {

  int thread_id = get_thread_id( core, thread );
  update_threading_struct_with_barrier_info( thread_id, thread_barrier_data, thread_barrier );

  run_setup( threading[thread_id] );

  status[1] = g.coarse_iter_count;
  status[0] = 1;
}


void dd_alpha_amg_setup_update( int iterations, int *status ) {

  if ( g.conf_flag == 1 ) {
    g.conf_flag = 0;
#pragma omp parallel num_threads(threading[0]->n_core)
    {
    if ( g.mixed_precision )
      re_setup_float( &l, threading[omp_get_thread_num()] );
    else
      re_setup_double( &l, threading[omp_get_thread_num()] );
    }
  }
  
  g.coarse_iter_count = 0;
#pragma omp parallel num_threads(threading[0]->n_core)
  {
  method_update( iterations, &l, threading[omp_get_thread_num()] );
  }
  
  status[0] = 1;
  status[1] = g.coarse_iter_count;
  status[0] = 1;
}
void dd_alpha_amg_setup_update_external_threading( int iterations, int *status,
    int core, int thread, void *thread_barrier_data, void (*thread_barrier)(void *, int) ) {

  int thread_id = get_thread_id( core, thread );
  update_threading_struct_with_barrier_info( thread_id, thread_barrier_data, thread_barrier );

  run_setup_update( threading[thread_id] );

  status[0] = 1;
  status[1] = g.coarse_iter_count;
}


double dd_alpha_amg_wilson_solve( double *vector_out, double *vector_in, double tol, double scale_even, double scale_odd, int *status ) {
  
  int t, z, y, x, i, j, k, clover_size, *ll = l.local_lattice;
  
  double *source = NULL, *solution = NULL;
  config_double clover_tmp = NULL;
  
  if ( g.csw != 0 )
    clover_size = 42*l.num_inner_lattice_sites;
  else
    clover_size = 12*l.num_inner_lattice_sites;
  
  MALLOC( source, double, 2*l.inner_vector_size );
  MALLOC( solution, double, 2*l.inner_vector_size );
  MALLOC( clover_tmp, complex_double, clover_size );
  
  g.coarse_iter_count = 0;
  g.iter_count = 0;
  g.p.tol = tol;
  g.p_MP.dp.tol = tol;
  
  for ( j=0, t=0; t<ll[T]; t++ )
    for ( z=0; z<ll[Z]; z++ )
      for ( y=0; y<ll[Y]; y++ )
        for ( x=0; x<ll[X]; x++ ) {
          i = g.vector_index_fct( t, z, y, x );
          for ( k=0; k<24; k++, j++ )
            source[j] = vector_in[i+k];
        }
  
  vector_double_copy( clover_tmp, g.op_double.clover, 0, clover_size, &l );  
  scale_clover( &(g.op_double), scale_even, scale_odd, &l );
  
  if ( g.mixed_precision ) {
    operator_updates_float( &l ); 
  } else {
    operator_updates_double( &l );
  }
  
#pragma omp parallel num_threads(threading[0]->n_core)
  {
  wilson_driver( (vector_double)solution, (vector_double)source, &l, threading[omp_get_thread_num()] );
  }
  
  vector_double_copy( g.op_double.clover, clover_tmp, 0, clover_size, &l );
  if ( g.mixed_precision ) {
    operator_updates_float( &l );
  } else {
    operator_updates_double( &l );
  }
  
  for ( j=0, t=0; t<ll[T]; t++ )
    for ( z=0; z<ll[Z]; z++ )
      for ( y=0; y<ll[Y]; y++ )
        for ( x=0; x<ll[X]; x++ ) {
          i = g.vector_index_fct( t, z, y, x );
          for ( k=0; k<24; k++, j++ )
            vector_out[i+k] = solution[j];
        }
  
  status[0] = g.iter_count;
  status[1] = g.coarse_iter_count;
  
  FREE( source, double, 2*l.inner_vector_size );
  FREE( solution, double, 2*l.inner_vector_size );
  FREE( clover_tmp, complex_double, clover_size );
  
  if ( g.norm_res > tol )
    status[0] = -1;
  
  return g.norm_res;
}


void dd_alpha_amg_free( void ) {
  
  finalize_common_thread_data(commonthreaddata);
  finalize_no_threading(no_threading);
  method_free( &l );
  method_finalize( &l );
}
