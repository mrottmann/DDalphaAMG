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
 
#include "main.h"
#include "DDalphaAMG.h"

global_struct g;
static level_struct l;
static int (*conf_index_fct)(int t, int z, int y, int x, int mu); 
static int (*vector_index_fct)(int t, int z, int y, int x); 
struct common_thread_data *commonthreaddata;
struct Thread **threading;
struct Thread *no_threading;

/*
 * Mirror of 
 *   method_init( &argc, &argv, &l );
 *   setup_threading( ... );
 */
void DDalphaAMG_initialize( DDalphaAMG_init *mg_init, DDalphaAMG_parameters *mg_params, DDalphaAMG_status *mg_status ) {
  
  double t0, t1;
  t0 = MPI_Wtime();

  /*
   * BEGIN: method_init( &argc, &argv, &l );
   */
  predefine_rank(mg_init->comm_cart);

  l_init( &l );
  g_init( &l );

  set_DDalphaAMG_parameters( mg_init, &l );
    
  conf_index_fct = NULL;
  vector_index_fct = NULL;

  int topo_type;
  MPI_Topo_test(mg_init->comm_cart, &topo_type);

  if(mg_init->Cart_coords==NULL)
    g.Cart_coords = MPI_Cart_coords;
  else
    g.Cart_coords = mg_init->Cart_coords;

  if(mg_init->Cart_rank==NULL)
    g.Cart_rank = MPI_Cart_rank;
  else
    g.Cart_rank = mg_init->Cart_rank;

  if( (g.Cart_rank==MPI_Cart_rank ||
       g.Cart_coords==MPI_Cart_coords) && topo_type!=MPI_CART ){
    warning0("Defining a new communicator\n");
    cart_define( mg_init->comm_cart, &l );
  } else
    cart_validate( mg_init->comm_cart, &l );
  
  if (mg_init->rnd_seeds == NULL)
    srand( time( 0 ) + 1000*g.my_rank );
  else
    srand( mg_init->rnd_seeds[g.my_rank] );

  data_layout_init( &l );
  operator_double_alloc( &(g.op_double), _ORDINARY, &l );
  operator_double_define( &(g.op_double), &l );
  MALLOC( g.odd_even_table, int, l.num_inner_lattice_sites );
  define_odd_even_table( &l );
  
  /*
   * BEGIN: setup_threading( ... );
   */
  no_threading = NULL;
  no_threading = (struct Thread *)malloc(sizeof(struct Thread));
  setup_no_threading(no_threading, &l);

  commonthreaddata = NULL;
  commonthreaddata = (struct common_thread_data *)malloc(sizeof(struct common_thread_data));
  init_common_thread_data(commonthreaddata);

  threading = NULL;
  MALLOC( threading, struct Thread *, g.num_openmp_processes);
  for(int i=0; i<g.num_openmp_processes; i++) {
    threading[i] = NULL;
    MALLOC( threading[i], struct Thread, 1);
  }
  THREADED(g.num_openmp_processes)
  setup_threading(threading[omp_get_thread_num()], commonthreaddata, &l);

  g.conf_flag = 0;
  g.setup_flag = 0;

  DDalphaAMG_get_parameters( mg_params );

  t1 = MPI_Wtime();

  mg_status->success = g.num_levels;// 1: OK, 2: re_setup done
  mg_status->time = t1-t0;
  
}

void DDalphaAMG_update_parameters( DDalphaAMG_parameters *mg_params, DDalphaAMG_status * mg_status ) {
  
  int i, j, re_setup=0, re_projs=0, re_dirac=0;
  level_struct *l_tmp;
  double t0, t1, m0;
  t0 = MPI_Wtime();
  g.coarse_time = 0;
  g.iter_count = 0;
  g.coarse_iter_count = 0;

  // int method;
  if ( mg_params->method != g.method ) {
    g.method = mg_params->method;
    if( g.setup_flag ) {
      //TODO: test which cases work and what to do for making the other working
      warning0("Change of method parameter after setup not guaranteed\n");
    }
  } 

  // int interpolation;
  if ( g.interpolation != mg_params->interpolation ) {
    g.interpolation = mg_params->interpolation;
    if( g.setup_flag ) {
      //TODO: test which cases work and what to do for making the other working
      warning0("Change of interpolation parameter after setup not guaranteed\n");
    }
  } 
  
  // int mixed_precision;
  if ( mg_params->mixed_precision != g.mixed_precision ) {
    g.mixed_precision = mg_params->mixed_precision;
#ifndef INIT_ONE_PREC //change between 1 and 2 allowed
    if( g.setup_flag && mg_params->mixed_precision * g.mixed_precision == 0 ) {
      warning0("Change from mixed_precision==0 to !=0 (or viceversa) needs a new setup.\n");
      re_setup++;
    }
#else
    warning0("Change of mixed_precision needs a new setup.\n");
    re_setup++;
#endif
  }

  // int block_lattice[MAX_MG_LEVELS][4];
  for ( i=0; i<g.num_levels; i++ )
    for ( j=0; j<4; j++ )
      if (g.block_lattice[i][j] != mg_params->block_lattice[i][j]) {
        g.block_lattice[i][j] = mg_params->block_lattice[i][j];
        parameter_update(&l);
        if (g.setup_flag) {
          warning0("Change of block_lattice needs a new setup.\n");
          re_setup++;
        }
      }
  
  // int mg_basis_vectors[MAX_MG_LEVELS-1];
  l_tmp=&l;
  for ( i=0; i<g.num_levels; i++ ){
    if ( mg_params->mg_basis_vectors[i] != g.num_eig_vect[i] ) {
      g.num_eig_vect[i] = mg_params->mg_basis_vectors[i];
      if( i==0 )
        parameter_update(&l);
      if( g.setup_flag ) {
        if ( mg_params->mg_basis_vectors[i] < g.num_eig_vect[i] )
          re_projs++; //TODO: check if this works
        else { //TODO just compute the extra vectors
          warning0("Increasing mg_basis_vectors needs a new setup.\n");
          re_setup++;
        }
      }
    }
    if( g.setup_flag )
      l_tmp = l_tmp->next_level;
  }
  
  // int setup_iterations[MAX_MG_LEVELS];
  l_tmp=&l;
  for ( i=0; i<g.num_levels; i++ ){
    if ( mg_params->setup_iterations[i] != g.setup_iter[i] ) {
      g.setup_iter[i] = mg_params->setup_iterations[i];
      if( (g.setup_flag && i>0) || (!g.setup_flag && i==0) ) 
        //after setup, l.setup_iter[i] is used as a counter for total number of setup iters
        l_tmp->setup_iter = mg_params->setup_iterations[i];
    }
    if( g.setup_flag )
      l_tmp = l_tmp->next_level;
  }

  // int smoother_iterations;
  l_tmp=&l;
  for ( i=0; i<g.num_levels; i++ ){
    if ( mg_params->smoother_iterations != g.post_smooth_iter[i] ) {
      l_tmp->post_smooth_iter = mg_params->smoother_iterations;
      g.post_smooth_iter[i] = mg_params->smoother_iterations;
    }
    if( g.setup_flag )
      l_tmp = l_tmp->next_level;
  }

  l_tmp=&l;
  for ( i=0; i<g.num_levels; i++ ){
    if (l_tmp->level > 0) {
      // double kcycle_tolerance;
      if ( mg_params->kcycle_tolerance != g.kcycle_tol ) {
        g.kcycle_tol = mg_params->kcycle_tolerance;
        if( g.setup_flag || i==0 ) {  
          if ( g.mixed_precision )
            l_tmp->p_float.tol = g.kcycle_tol;
          else
            l_tmp->p_float.tol = g.kcycle_tol;
        }
      }
    } else {
      // double coarse_tolerance;
      if ( mg_params->coarse_tolerance != g.coarse_tol ){
        g.coarse_tol = mg_params->coarse_tolerance;
        if (g.setup_flag && g.mixed_precision )
          l_tmp->p_float.tol = g.coarse_tol;
        else if(g.setup_flag)
          l_tmp->p_float.tol = g.coarse_tol;
      }
    }
    
    if( g.setup_flag )
      l_tmp = l_tmp->next_level;
    else
      break;
  }
  
  // double kappa;
  m0 = 1./(2.*mg_params->kappa)-4.; 
  if( g.m0 != m0 ){
    g.m0 = m0;
    THREADED(threading[0]->n_core)
      if ( g.setup_flag )
        m0_update( g.m0, &l, threading[omp_get_thread_num()] );
      else if ( g.conf_flag )
        m0_update_double( g.m0, &(g.op_double), &l, threading[omp_get_thread_num()] );
    re_dirac++;
  }
  
  // double mu;
  // double mu_odd_shift;
  // double mu_even_shift;
  // double mu_factor[MAX_MG_LEVELS];
#ifdef HAVE_TM
  int update_mu = 0;
  for ( i=0; i<g.num_levels; i++ )
    if (mg_params->mu_factor[i] != g.mu_factor[i] ) {
      g.mu_factor[i] = mg_params->mu_factor[i];
      update_mu = 1;
    }

  if( update_mu || mg_params->mu != g.mu || mg_params->mu_odd_shift != g.mu_odd_shift || 
      mg_params->mu_even_shift != g.mu_even_shift ){
    g.setup_mu = mg_params->mu;
    g.mu = mg_params->mu;
    g.mu_even_shift = mg_params->mu_even_shift;
    g.mu_odd_shift = mg_params->mu_odd_shift;
    THREADED(threading[0]->n_core)
      if ( g.setup_flag )
        tm_term_update( g.mu, &l, threading[omp_get_thread_num()] );
      else if ( g.conf_flag )
        tm_term_double_setup( g.mu, g.mu_even_shift, g.mu_odd_shift, &(g.op_double), &l, threading[omp_get_thread_num()] ); 
    re_dirac++;
  }
  
#else
  if ( mg_params->mu != 0 || mg_params->mu_odd_shift != 0 || mg_params->mu_even_shift != 0 )
    warning0("Parameters mu, mu_odd_shift, mu_even_shift not supported without defining HAVE_TM flag.");
#endif

#ifdef HAVE_TM1p1  
  int update_eps = 0;

  for ( i=0; i<g.num_levels; i++ )
    if (mg_params->epsbar_factor[i] != g.epsbar_factor[i] ) {
      g.epsbar_factor[i] = mg_params->epsbar_factor[i];
      update_eps = 1;
     }

  if( update_eps || mg_params->epsbar != g.epsbar || mg_params->epsbar_ig5_odd_shift != g.epsbar_ig5_odd_shift || mg_params->epsbar_ig5_even_shift != g.epsbar_ig5_even_shift ){
    g.epsbar = mg_params->epsbar;
    g.epsbar_ig5_even_shift = mg_params->epsbar_ig5_even_shift;
    g.epsbar_ig5_odd_shift = mg_params->epsbar_ig5_odd_shift;
    THREADED(threading[0]->n_core)
      if ( g.setup_flag )
        epsbar_term_update( &l, threading[omp_get_thread_num()] );
      else if ( g.conf_flag )
        epsbar_term_double_setup( g.epsbar, g.epsbar_ig5_even_shift, g.epsbar_ig5_odd_shift, &(g.op_double), &l, threading[omp_get_thread_num()] ); 
    re_dirac++;
  }
  
#else
  if ( mg_params->epsbar != 0 || mg_params->epsbar_ig5_odd_shift != 0 || mg_params->epsbar_ig5_even_shift != 0 )
    warning0("Parameters epsbar, epsbar_odd_shift, epsbar_even_shift not supported without defining HAVE_TM1p1 flag.");
#endif

  // int (*conf_index_fct)(int t, int z, int y, int x, int mu);
  // int (*vector_index_fct)(int t, int z, int y, int x );
  conf_index_fct = mg_params->conf_index_fct;
  vector_index_fct = mg_params->vector_index_fct;
  
  // int print;
  g.print = mg_params->print;
  
  // UPDATING
  if ( re_setup && g.setup_flag ){ // destroy and repeate setup
    DDalphaAMG_setup( mg_status ); // TODO handle status

  } else if ( re_projs && g.setup_flag ) { //project again the operators
    if ( g.mixed_precision )
      THREADED(threading[0]->n_core)
        re_setup_float( &l, threading[omp_get_thread_num()] ); 
    else
      THREADED(threading[0]->n_core)
        re_setup_double( &l, threading[omp_get_thread_num()] );
    
  } else if ( re_dirac && g.setup_flag ) { //update just the oddeven and vecorized operators
    THREADED(threading[0]->n_core)
      finalize_operator_update( &l, threading[omp_get_thread_num()]);
  }
  
  
  DDalphaAMG_get_parameters( mg_params );
  
  t1 = MPI_Wtime();
  
  mg_status->success = 1+re_setup;// 1: OK, 2: re_setup done
  mg_status->time = t1-t0;
  mg_status->coarse_time = g.coarse_time;
  mg_status->iter_count = g.iter_count;
  mg_status->coarse_iter_count = g.coarse_iter_count;
}

void DDalphaAMG_change_mu_sign( DDalphaAMG_status *mg_status ) {
  
  double t0, t1;
  t0 = MPI_Wtime();
  g.coarse_time = 0;
  g.iter_count = 0;
  g.coarse_iter_count = 0;
  mg_status->success = 0;
  mg_status->info = 0;  
  
  g.mu *= -1;
  g.mu_even_shift *= -1;
  g.mu_odd_shift *= -1;

  if (g.conf_flag && !g.setup_flag ) {
    
    THREADED(threading[0]->n_core)
    tm_term_double_setup( g.mu, g.mu_even_shift, g.mu_odd_shift, &(g.op_double), &l, threading[omp_get_thread_num()]);
    
  } else if (g.conf_flag && g.setup_flag )
    THREADED(threading[0]->n_core) {
      tm_term_update( g.mu, &l, threading[omp_get_thread_num()] );
      finalize_operator_update( &l, threading[omp_get_thread_num()] );
    }
  t1 = MPI_Wtime();
  
  mg_status->success = 1;// 1: OK, 2: re_setup done
  mg_status->time = t1-t0;
  mg_status->info = g.mu;
  mg_status->coarse_time = g.coarse_time;
  
}

void DDalphaAMG_change_epsbar_shift_sign( DDalphaAMG_status *mg_status ) {
  
  double t0, t1;
  t0 = MPI_Wtime();
  g.coarse_time = 0;
  g.iter_count = 0;
  g.coarse_iter_count = 0;
  mg_status->success = 0;
  mg_status->info = 0;  
  
  if ( g.epsbar_ig5_even_shift !=0 || g.epsbar_ig5_odd_shift !=0 ) {
    g.epsbar_ig5_even_shift *= -1;
    g.epsbar_ig5_odd_shift *= -1;

    if (g.conf_flag && !g.setup_flag ) {
      
      THREADED(threading[0]->n_core) {
        tm_term_double_setup( g.mu, g.mu_even_shift, g.mu_odd_shift, &(g.op_double),
                              &l, threading[omp_get_thread_num()]);
        epsbar_term_double_setup( g.epsbar, g.epsbar_ig5_even_shift, g.epsbar_ig5_odd_shift, &(g.op_double),
                                  &l, threading[omp_get_thread_num()]);
      }
    } else if (g.conf_flag && g.setup_flag )
      THREADED(threading[0]->n_core) {
        tm_term_update( g.mu, &l, threading[omp_get_thread_num()] );
        epsbar_term_update( &l, threading[omp_get_thread_num()] );
        finalize_operator_update( &l, threading[omp_get_thread_num()] );
      }
  }
  t1 = MPI_Wtime();
  
  mg_status->success = 1;// 1: OK, 2: re_setup done
  mg_status->time = t1-t0;
  mg_status->info = g.mu;
  mg_status->coarse_time = g.coarse_time;
  
}


void DDalphaAMG_set_configuration( double *gauge_field, DDalphaAMG_status *mg_status ) {
  
  int t, z, y, x, mu, i, j, k;
  double t0, t1;
  SU3_storage U;
  complex_double phase[4] = { _COMPLEX_double_ZERO, _COMPLEX_double_ZERO, _COMPLEX_double_ZERO, _COMPLEX_double_ZERO};
  int *ll=l.local_lattice, onb[4];
  for (i=0; i<4; i++)
    onb[i] = (g.my_coords[i]==g.process_grid[i]-1)?1:0;

  t0 = MPI_Wtime();
  mg_status->success = 0;
  mg_status->info = 0;  

  // START: dirac_setup
  if ( g.print > 0 ) printf0("%s\n", CLIFFORD_BASIS );
  if ( g.bc == _ANTIPERIODIC ) printf0("antiperiodic in time");
  else if ( g.bc == _TWISTED ) printf0("twisted (%.2f, %.2f, %.2f, %.2f)", g.twisted_bc[0], 
               g.twisted_bc[1], g.twisted_bc[2], g.twisted_bc[3]);
  else printf0("periodic in time");
  printf0(" boundary conditions \n");

  SU3_storage_alloc( &U, &l );
  
  if(g.bc == _ANTIPERIODIC && onb[T] ) {
    phase[Z] = 1; phase[Y] = 1; phase[X] = 1;
    for ( t=1, i=0, k=0; t<ll[T]+1; t++ ) {
      if (t<ll[T]) phase[T] = 1; 
      else phase[T] = -1;
      for ( z=1; z<ll[Z]+1; z++ )
        for ( y=1; y<ll[Y]+1; y++ )
          for ( x=1; x<ll[X]+1; x++ )
            for ( mu=0; mu<4; mu++ ) {
              if ( conf_index_fct != NULL )
                k = conf_index_fct( t-1, z-1, y-1, x-1, mu );
              for (j=0; j<9; j++, i++, k+=2) {
                g.op_double.D[i] = 0.5*phase[mu]*(gauge_field[k]+I*gauge_field[k+1]);
                U[t][z][y][x][mu][j] = phase[mu]*(gauge_field[k]+I*gauge_field[k+1]);
              }
            }     
    }
  }
  else if(g.bc == _TWISTED && ( onb[T] || onb[Z] || onb[Y] || onb[X] ))
    for ( t=1, i=0, k=0; t<ll[T]+1; t++ ) {
      if ( !onb[T] || t<ll[T] || g.twisted_bc[T]==0) phase[T] = 1; 
      else phase[T] = cexp(I*g.twisted_bc[T]);
      for ( z=1; z<ll[Z]+1; z++ ) {
        if ( !onb[Z] || z<ll[Z] || g.twisted_bc[Z]==0) phase[Z] = 1; 
        else phase[Z] = cexp(I*g.twisted_bc[Z]);
        for ( y=1; y<ll[Y]+1; y++ ) {
          if ( !onb[Y] || y<ll[Y] || g.twisted_bc[Y]==0) phase[Y] = 1; 
          else phase[Y] = cexp(I*g.twisted_bc[Y]);
          for ( x=1; x<ll[X]+1; x++ ) {
            if ( !onb[X] || x<ll[X] || g.twisted_bc[X]==0) phase[X] = 1; 
            else phase[X] = cexp(I*g.twisted_bc[X]);
            for ( mu=0; mu<4; mu++ ) {
              if ( conf_index_fct != NULL )
                k = conf_index_fct( t-1, z-1, y-1, x-1, mu );
              for (j=0; j<9; j++, i++, k+=2) {
                g.op_double.D[i] = 0.5*phase[mu]*(gauge_field[k]+I*gauge_field[k+1]);
                U[t][z][y][x][mu][j] = phase[mu]*(gauge_field[k]+I*gauge_field[k+1]);
              }
            }
          }
        }
      }
    }
  
  else
    for ( t=1, i=0, k=0; t<ll[T]+1; t++ )
      for ( z=1; z<ll[Z]+1; z++ )
        for ( y=1; y<ll[Y]+1; y++ )
          for ( x=1; x<ll[X]+1; x++ )
            for ( mu=0; mu<4; mu++ ) {
              if ( conf_index_fct != NULL )
                k = conf_index_fct( t-1, z-1, y-1, x-1, mu );
              for (j=0; j<9; j++, i++, k+=2) {
                g.op_double.D[i] = 0.5*(gauge_field[k]+I*gauge_field[k+1]);
                U[t][z][y][x][mu][j] = (gauge_field[k]+I*gauge_field[k+1]);
              }
            }
  
  SU3_ghost_update( &U, &l );
  if ( g.print > 0 ) printf0("Configuration stored...\n");
  
  compute_clover_term( U, &l );
  
  // calculate the plaquette
  g.plaq_clov = calc_plaq( U, &l );
  if(g.print > 0) printf0("average plaquette: %.13lf\n", g.plaq_clov);
    
  SU3_storage_free( &U, &l );
  //END: dirac_setup
  
  mg_status->success = 1;
  g.conf_flag = 1;
  mg_status->info = g.plaq;
  
  if ( g.setup_flag == 1 ) {
    if(l.s_double.op.clover != NULL)
      schwarz_double_setup( &(l.s_double), &(g.op_double), &l );
    if(l.s_float.op.clover != NULL)
      schwarz_float_setup( &(l.s_float), &(g.op_double), &l );
    THREADED(threading[0]->n_core)
    if ( g.mixed_precision ) 
      operator_updates_float( &l, threading[omp_get_thread_num()] );
    else
      operator_updates_double( &l, threading[omp_get_thread_num()] );
    
    mg_status->success++; //0: error, 1: OK, 2 re_setup done
  }    

  t1 = MPI_Wtime();
  mg_status->time = t1-t0;
  mg_status->coarse_time = 0;
  mg_status->iter_count = 0;
  mg_status->coarse_iter_count = 0;
  
}

void DDalphaAMG_setup( DDalphaAMG_status * mg_status ) {
  
  double t0, t1;
  t0 = MPI_Wtime();
  g.coarse_time = 0;
  g.iter_count = 0;
  g.coarse_iter_count = 0;
  mg_status->success = 0;
  mg_status->info = 0;
  
  if(g.conf_flag == 1) {
    if ( g.setup_flag )
      method_free( &l );
    THREADED(threading[0]->n_core)
    {
      method_setup( NULL, &l, threading[omp_get_thread_num()] );
      method_update( g.setup_iter[0], &l, threading[omp_get_thread_num()] );
    }
    g.setup_flag = 1;
        
    l.setup_iter = g.setup_iter[0];
    mg_status->success = l.setup_iter;
    
  }
  
  t1 = MPI_Wtime();
  mg_status->time = t1-t0;
  mg_status->coarse_time = g.coarse_time;
  mg_status->iter_count = g.iter_count;
  mg_status->coarse_iter_count = g.coarse_iter_count;
}

void DDalphaAMG_update_setup( int iterations, DDalphaAMG_status * mg_status ) {

  if(g.setup_flag) {
    double t0, t1;
    t0 = MPI_Wtime();
    g.coarse_time = 0;
    g.iter_count = 0;
    g.coarse_iter_count = 0;
    mg_status->success = 0;
    mg_status->info = 0;

    THREADED(threading[0]->n_core) 
    method_update( iterations, &l, threading[omp_get_thread_num()] );
    //    method_update( iterations, &l, no_threading );
    
    l.setup_iter += iterations;
    mg_status->success = l.setup_iter;
    t1 = MPI_Wtime();
    mg_status->time = t1-t0;
    mg_status->coarse_time = g.coarse_time;
    mg_status->iter_count = g.iter_count;
    mg_status->coarse_iter_count = g.coarse_iter_count;

  }
  else {
    g.setup_iter[0] = iterations;
    DDalphaAMG_setup( mg_status );
  }
}

enum {_SOLVE, _SOLVE_SQ, _SOLVE_SQ_ODD, _SOLVE_SQ_EVEN, _PRECOND, _OPERATOR};
static inline void DDalphaAMG_driver( double *vector1_out, double *vector1_in, double *vector2_out, double *vector2_in, double tol, DDalphaAMG_status *mg_status, int _TYPE ) {
  
  int t, z, y, x, i, j, k, mu, *ll = l.local_lattice, *gl=l.global_lattice, sl[4], precision_changed;
  complex_double twisted_bc, tmp1, tmp2;
  double phase[4] = {_COMPLEX_double_ZERO, _COMPLEX_double_ZERO, _COMPLEX_double_ZERO, _COMPLEX_double_ZERO}, vmin=1, vmax=EPS_float, vtmp;
  gmres_double_struct *p = g.mixed_precision==2?&(g.p_MP.dp):&(g.p);
  vector_double vb, rhs = p->b;
  vector_double vx, sol = p->x;
  DDalphaAMG_status tmp_status;

  double t0, t1;
  t0 = MPI_Wtime();
  g.coarse_time = 0;
  g.iter_count = 0;
  g.coarse_iter_count = 0;
  mg_status->success = 0;
  mg_status->info = 0;
  
  ASSERT(vector1_out!=NULL);
  ASSERT(vector1_in!=NULL);
#ifdef HAVE_TM1p1
  if(g.n_flavours==2) {
    ASSERT(vector2_out!=NULL);
    ASSERT(vector2_in!=NULL);
  }
#endif

  if(g.mixed_precision!=2)
    g.p.tol = tol;
  else
    g.p_MP.dp.tol = tol;

  for (i=0; i<4; i++)
    sl[i] = ll[i]*g.my_coords[i];
 
  for (t=0, j=0; t<ll[T]; t++) {
    if (g.bc==_TWISTED) phase[T] = g.twisted_bc[T]*((double)sl[T]+t)/(double)gl[T];
    for (z=0; z<ll[Z]; z++) {
      if (g.bc==_TWISTED) phase[Z] = phase[T] + g.twisted_bc[Z]*((double)sl[Z]+z)/(double)gl[Z];
      for (y=0; y<ll[Y]; y++) {
        if (g.bc==_TWISTED) phase[Y] = phase[Z] + g.twisted_bc[Y]*((double)sl[Y]+y)/(double)gl[Y];
        for (x=0; x<ll[X]; x++) {
          if (g.bc==_TWISTED) {
            phase[X] = phase[Y] + g.twisted_bc[X]*((double)sl[X]+x)/(double)gl[X];
            twisted_bc = cexp(I*phase[X]);
          } else
            twisted_bc = 1.;
          if(vector_index_fct!=NULL )
            i = vector_index_fct( t, z, y, x );
          else 
            i = 2*j;
          
#ifdef HAVE_TM1p1
          if(g.n_flavours==2) {
            for ( mu=0; mu<4; mu++ ) {
              for ( k=0; k<3; k++, j++ ) {
#ifndef BASIS4 
                rhs[j] = ((complex_double)vector1_in[i+2*(k+3*mu)] + I*(complex_double)vector1_in[i+2*(k+3*mu)+1]) * twisted_bc;
                rhs[j+6] = ((complex_double)vector2_in[i+2*(k+3*mu)] + I*(complex_double)vector2_in[i+2*(k+3*mu)+1]) * twisted_bc;

#else
                rhs[j] = ((complex_double)vector1_in[i+2*(k+3*(3-mu))] + I*(complex_double)vector1_in[i+2*(k+3*(3-mu))+1]) * twisted_bc;
                rhs[j+6] = ((complex_double)vector2_in[i+2*(k+3*(3-mu))] + I*(complex_double)vector2_in[i+2*(k+3*(3-mu))+1]) * twisted_bc;
#endif
                
#ifndef INIT_ONE_PREC
                if(g.mixed_precision==2) {
                  vtmp=cabs(rhs[j]);
                  if(vtmp > vmax)
                    vmax=vtmp;
                  if( vtmp > EPS_double && vtmp < vmin )
                    vmin=vtmp;
                  vtmp=cabs(rhs[j+6]);
                  if(vtmp > vmax)
                    vmax=vtmp;
                  if( vtmp > EPS_double && vtmp < vmin )
                    vmin=vtmp;
                }
              }
#endif
              if(mu%2)
                j+=6;
            }
          } else
#endif
            for ( mu=0; mu<4; mu++ )
              for ( k=0; k<3; k++, j++ ) {
#ifndef BASIS4 
                rhs[j] = ((complex_double)vector1_in[i+2*(k+3*mu)] + I*(complex_double)vector1_in[i+2*(k+3*mu)+1]) * twisted_bc;
#else
                rhs[j] = ((complex_double)vector1_in[i+2*(k+3*(3-mu))] + I*(complex_double)vector1_in[i+2*(k+3*(3-mu))+1]) * twisted_bc;
#endif
                
#ifndef INIT_ONE_PREC
                if(g.mixed_precision==2) {
                  vtmp=cabs(rhs[j]);
                  if(vtmp > vmax)
                    vmax=vtmp;
                  if( vtmp > EPS_double && vtmp < vmin )
                    vmin=vtmp;
                }
              }
#endif
        }
      }
    }
  }
   
#ifndef INIT_ONE_PREC
    
  double gvmin, gvmax;
  if(g.mixed_precision==2) {
    MPI_Allreduce(&vmin, &gvmin, 1, MPI_DOUBLE, MPI_MIN, g.comm_cart);
    MPI_Allreduce(&vmax, &gvmax, 1, MPI_DOUBLE, MPI_MAX, g.comm_cart);
  }
  
  //switching to double precision on the fine level
  if(g.mixed_precision==2 && gvmin/gvmax<EPS_float) {
    warning0("Changing solver precision on fine level due to rhs elements (min/max=%e)\n", vmin/vmax);
    precision_changed=1;
    g.mixed_precision=1;
    p = &(g.p);
    // storing pointer in x and b
    vb = p->b; 
    vx = p->x;
    p->b = g.p_MP.dp.b;
    p->x = g.p_MP.dp.x;
    p->tol = g.p_MP.dp.tol;
  } else precision_changed = 0;
#endif

  switch(_TYPE) {
    
  case _SOLVE :
    THREADED(threading[0]->n_core)
    if ( g.method == -1 ) {
      cgn_double( &(g.p), &l, threading[omp_get_thread_num()] );
    } else if ( g.mixed_precision == 2 ) {
      fgmres_MP( &(g.p_MP), &l, threading[omp_get_thread_num()] );
    } else {
      fgmres_double( &(g.p), &l, threading[omp_get_thread_num()] );
    }
    break;

  case _SOLVE_SQ :
    THREADED(threading[0]->n_core) {
#ifdef HAVE_TM1p1
      if(g.n_flavours==2) 
        // sol = (D_h^{-1})*g5*tau1*(D_h^{-1})*g5*tau1*rhs
        tau1_gamma5_double( rhs, rhs, &l, threading[omp_get_thread_num()] );
      else
#endif
        // sol = (D_d^{-1})*g5*(D_u^{-1})*g5*rhs
        gamma5_double( rhs, rhs, &l, threading[omp_get_thread_num()] );
      if ( g.method == -1 ) {
        cgn_double( &(g.p), &l, threading[omp_get_thread_num()] );
      } else if ( g.mixed_precision == 2 ) {
        fgmres_MP( &(g.p_MP), &l, threading[omp_get_thread_num()] );
      } else {
        fgmres_double( &(g.p), &l, threading[omp_get_thread_num()] );
      }
      
#ifdef HAVE_TM1p1
      if(g.n_flavours==2) 
        tau1_gamma5_double(rhs, sol, &l, threading[omp_get_thread_num()] );
      else
#endif
        gamma5_double(rhs, sol, &l, threading[omp_get_thread_num()] );
    }
#ifdef HAVE_TM1p1
    if(g.n_flavours==2) 
      DDalphaAMG_change_epsbar_shift_sign( &tmp_status );
    else
#endif
      DDalphaAMG_change_mu_sign( &tmp_status );
    THREADED(threading[0]->n_core) {
      if ( g.method == -1 ) {
        cgn_double( &(g.p), &l, threading[omp_get_thread_num()] );
      } else if ( g.mixed_precision == 2 ) {
        fgmres_MP( &(g.p_MP), &l, threading[omp_get_thread_num()] );
      } else {
        fgmres_double( &(g.p), &l, threading[omp_get_thread_num()] );
      }
    }
    // DDalphaAMG_change_mu_sign( &tmp_status );
    warning0("sign of mu changed during the inversion of squared operator\n");
    break;
    
  case _SOLVE_SQ_ODD :    
    THREADED(threading[0]->n_core)
    {
#ifdef HAVE_TM1p1
      if(g.n_flavours==2) 
        // sol = (D_h^{-1})*g5*tau1*(D_h^{-1})*g5*tau1*rhs
        tau1_gamma5_set_even_to_zero_double(rhs, rhs, &l, threading[omp_get_thread_num()]);
      else
#endif
        // sol = (D_d^{-1})*g5*(D_u^{-1})*g5*rhs
        gamma5_set_even_to_zero_double(rhs, rhs, &l, threading[omp_get_thread_num()]);
      if ( g.method == -1 ) {
        cgn_double( &(g.p), &l, threading[omp_get_thread_num()] );
      } else if ( g.mixed_precision == 2 ) {
        fgmres_MP( &(g.p_MP), &l, threading[omp_get_thread_num()] );
      } else {
        fgmres_double( &(g.p), &l, threading[omp_get_thread_num()] );
      }
#ifdef HAVE_TM1p1
      if(g.n_flavours==2) 
        tau1_gamma5_set_even_to_zero_double(rhs, sol, &l, threading[omp_get_thread_num()]);
      else
#endif
        gamma5_set_even_to_zero_double(rhs, sol, &l, threading[omp_get_thread_num()]);
    }
#ifdef HAVE_TM1p1
    if(g.n_flavours==2) 
      DDalphaAMG_change_epsbar_shift_sign( &tmp_status );
    else
#endif
      DDalphaAMG_change_mu_sign( &tmp_status );
    THREADED(threading[0]->n_core) {
      if ( g.method == -1 ) {
        cgn_double( &(g.p), &l, threading[omp_get_thread_num()] );
      } else if ( g.mixed_precision == 2 ) {
        fgmres_MP( &(g.p_MP), &l, threading[omp_get_thread_num()] );
      } else {
        fgmres_double( &(g.p), &l, threading[omp_get_thread_num()] );
      }
    }
    // DDalphaAMG_change_mu_sign( &tmp_status );
    warning0("sign of mu changed during the inversion of squared operator\n");
    break;
    
  case _SOLVE_SQ_EVEN :    
    THREADED(threading[0]->n_core)
    {
#ifdef HAVE_TM1p1
      if(g.n_flavours==2) 
        // sol = (D_h^{-1})*g5*tau1*(D_h^{-1})*g5*tau1*rhs
        tau1_gamma5_set_odd_to_zero_double(rhs, rhs, &l, threading[omp_get_thread_num()]);
      else
#endif
        // sol = (D_d^{-1})*g5*(D_u^{-1})*g5*rhs
        gamma5_set_odd_to_zero_double(rhs, rhs, &l, threading[omp_get_thread_num()]);
      if ( g.method == -1 ) {
        cgn_double( &(g.p), &l, threading[omp_get_thread_num()] );
      } else if ( g.mixed_precision == 2 ) {
        fgmres_MP( &(g.p_MP), &l, threading[omp_get_thread_num()] );
      } else {
        fgmres_double( &(g.p), &l, threading[omp_get_thread_num()] );
      }
#ifdef HAVE_TM1p1
      if(g.n_flavours==2) 
        tau1_gamma5_set_odd_to_zero_double(rhs, sol, &l, threading[omp_get_thread_num()]);
      else
#endif
        gamma5_set_odd_to_zero_double(rhs, sol, &l, threading[omp_get_thread_num()]);
    }
#ifdef HAVE_TM1p1
    if(g.n_flavours==2) 
      DDalphaAMG_change_epsbar_shift_sign( &tmp_status );
    else
#endif
      DDalphaAMG_change_mu_sign( &tmp_status );
    THREADED(threading[0]->n_core) {
      if ( g.method == -1 ) {
        cgn_double( &(g.p), &l, threading[omp_get_thread_num()] );
      } else if ( g.mixed_precision == 2 ) {
        fgmres_MP( &(g.p_MP), &l, threading[omp_get_thread_num()] );
      } else {
        fgmres_double( &(g.p), &l, threading[omp_get_thread_num()] );
      }
    }
    // DDalphaAMG_change_mu_sign( &tmp_status );
    warning0("sign of mu changed during the inversion of squared operator\n");
    break;

  case _PRECOND :
    THREADED(threading[0]->n_core)
    preconditioner( sol, NULL, rhs, _NO_RES, &l, threading[omp_get_thread_num()] );
    break;

  case _OPERATOR :
    THREADED(threading[0]->n_core)
    if ( g.mixed_precision == 2 ) {
      apply_operator_double( sol, rhs, &(g.p_MP.dp), &l, threading[omp_get_thread_num()] );
    } else {
      apply_operator_double( sol, rhs, &(g.p), &l, threading[omp_get_thread_num()] );
    }
    break;

  default :
    warning0("_TYPE not found in DDalphaAMG_driver. Returing vector in as vector out.");
    sol=rhs;
    break;
  }

  for (t=0, j=0; t<ll[T]; t++) {
    if (g.bc==_TWISTED) phase[T] = g.twisted_bc[T]*((double)sl[T]+t)/(double)gl[T];
    for (z=0; z<ll[Z]; z++) {
      if (g.bc==_TWISTED) phase[Z] = phase[T] + g.twisted_bc[Z]*((double)sl[Z]+z)/(double)gl[Z];
      for (y=0; y<ll[Y]; y++) {
        if (g.bc==_TWISTED) phase[Y] = phase[Z] + g.twisted_bc[Y]*((double)sl[Y]+y)/(double)gl[Y];
        for (x=0; x<ll[X]; x++) {
          if (g.bc==_TWISTED) {
            phase[X] = phase[Y] + g.twisted_bc[X]*((double)sl[X]+x)/(double)gl[X];
            twisted_bc = cexp(-I*phase[X]);
          } else
            twisted_bc = 1.;
          if(vector_index_fct!=NULL )
            i = vector_index_fct( t, z, y, x );
          else 
            i = 2*j;

#ifdef HAVE_TM1p1
          if(g.n_flavours==2) {
            for ( mu=0; mu<4; mu++ ) {
              for ( k=0; k<3; k++, j++ ) {
                tmp1 = sol[j] * twisted_bc;
                tmp2 = sol[j+6] * twisted_bc;
#ifndef BASIS4 
                vector1_out[i+2*(k+3*mu)]   = creal(tmp1);
                vector1_out[i+2*(k+3*mu)+1] = cimag(tmp1);
                vector2_out[i+2*(k+3*mu)]   = creal(tmp2);
                vector2_out[i+2*(k+3*mu)+1] = cimag(tmp2);
#else
                vector1_out[i+2*(k+3*(3-mu))]   = creal(tmp1);
                vector1_out[i+2*(k+3*(3-mu))+1] = cimag(tmp1);
                vector2_out[i+2*(k+3*(3-mu))]   = creal(tmp2);
                vector2_out[i+2*(k+3*(3-mu))+1] = cimag(tmp2);
#endif   
              }
              if(mu%2)
                j+=6;
            }
          } else
#endif
            for ( mu=0; mu<4; mu++ )
              for ( k=0; k<3; k++, j++ ) {
                tmp1 = sol[j] * twisted_bc;
#ifndef BASIS4 
                vector1_out[i+2*(k+3*mu)]   = creal(tmp1);
                vector1_out[i+2*(k+3*mu)+1] = cimag(tmp1);
#else
                vector1_out[i+2*(k+3*(3-mu))]   = creal(tmp1);
                vector1_out[i+2*(k+3*(3-mu))+1] = cimag(tmp1);
#endif   
              }
        }
      }
    }
  }
    
#ifndef INIT_ONE_PREC
  if (precision_changed) {
    g.mixed_precision=2;
    // recovering pointer from x and b
    p->b = vb; 
    p->x = vx;
  }
#endif
    
  if ( g.norm_res <= tol || _TYPE == _OPERATOR || _TYPE == _PRECOND )
    mg_status->success = 1;
  mg_status->info = g.norm_res;
  t1 = MPI_Wtime();
  mg_status->time = t1-t0;
  mg_status->coarse_time = g.coarse_time;
  mg_status->iter_count = g.iter_count;
  mg_status->coarse_iter_count = g.coarse_iter_count;
  
}

static inline void vector_copy( vector_double vector_out, vector_double vector_in )
{
  THREADED(threading[0]->n_core) {
    int start = threading[omp_get_thread_num()]->start_index[0], 
      end = threading[omp_get_thread_num()]->end_index[0];
    vector_double_copy( vector_out, vector_in, start, end, &l );
  }  
}

static inline void solver( )
{
  THREADED(threading[0]->n_core)
    if ( g.method == -1 ) {
      cgn_double( &(g.p), &l, threading[omp_get_thread_num()] );
    } else if ( g.mixed_precision == 2 ) {
      fgmres_MP( &(g.p_MP), &l, threading[omp_get_thread_num()] );
    } else {
      fgmres_double( &(g.p), &l, threading[omp_get_thread_num()] );
    }
}

static inline void correct_guess( vector_double guess, vector_double solution, vector_double solution2,
                                  double  even_dshift, double odd_dshift )
{
  // guess = D^{-1}*rhs - i*dshift*D^{-2}*rhs 
  THREADED(threading[0]->n_core) {
    int start = threading[omp_get_thread_num()]->start_index[0], 
      end = threading[omp_get_thread_num()]->end_index[0];
    if( odd_dshift == 0 || even_dshift == 0 || even_dshift == odd_dshift ) {
      double dshift = ( odd_dshift == 0 ) ? even_dshift:odd_dshift;
      printf0("correcting with dshift %le\n", dshift);
      vector_double_scale( guess, solution2, -I*dshift, g.p.v_start, g.p.v_end, &l );
      vector_double_plus( guess, guess, solution, start, end, &l );
    } else
      vector_double_copy( guess, solution, start, end, &l );
  }  
}


static inline void DDalphaAMG_ms_driver( double **vector1_out, double *vector1_in, 
                                         double **vector2_out, double *vector2_in, 
                                         double  *even_shifts, double *odd_shifts, int n_shifts,
                                         double tol, DDalphaAMG_status *mg_status, int _TYPE ) 
{
  int t, z, y, x, i, j, k, n, mu, *ll = l.local_lattice, *gl=l.global_lattice, sl[4], precision_changed;
  complex_double twisted_bc, tmp1, tmp2;
  double phase[4] = {_COMPLEX_double_ZERO, _COMPLEX_double_ZERO, _COMPLEX_double_ZERO, _COMPLEX_double_ZERO},
      vmin=1, vmax=EPS_float, vtmp;
  gmres_double_struct *p = g.mixed_precision==2?&(g.p_MP.dp):&(g.p);
  vector_double vb, rhs = p->b;
  vector_double vx, sol = p->x;
  vector_double source = NULL, solution = NULL, solution2 = NULL;
  DDalphaAMG_status tmp_status;

  double t0, t1;
  t0 = MPI_Wtime();
  g.coarse_time = 0;
  g.iter_count = 0;
  g.coarse_iter_count = 0;
  mg_status->success = 0;
  mg_status->info = 0;
  
  ASSERT(vector1_out!=NULL);
  ASSERT(vector1_in!=NULL);
#ifdef HAVE_TM1p1
  if(g.n_flavours==2) {
    ASSERT(vector2_out!=NULL);
    ASSERT(vector2_in!=NULL);
  }
#endif

  if(g.mixed_precision!=2)
    g.p.tol = tol;
  else
    g.p_MP.dp.tol = tol;

  for (i=0; i<4; i++)
    sl[i] = ll[i]*g.my_coords[i];
 
  for (t=0, j=0; t<ll[T]; t++) {
    if (g.bc==_TWISTED) phase[T] = g.twisted_bc[T]*((double)sl[T]+t)/(double)gl[T];
    for (z=0; z<ll[Z]; z++) {
      if (g.bc==_TWISTED) phase[Z] = phase[T] + g.twisted_bc[Z]*((double)sl[Z]+z)/(double)gl[Z];
      for (y=0; y<ll[Y]; y++) {
        if (g.bc==_TWISTED) phase[Y] = phase[Z] + g.twisted_bc[Y]*((double)sl[Y]+y)/(double)gl[Y];
        for (x=0; x<ll[X]; x++) {
          if (g.bc==_TWISTED) {
            phase[X] = phase[Y] + g.twisted_bc[X]*((double)sl[X]+x)/(double)gl[X];
            twisted_bc = cexp(I*phase[X]);
          } else
            twisted_bc = 1.;
          if(vector_index_fct!=NULL )
            i = vector_index_fct( t, z, y, x );
          else 
            i = 2*j;
          
#ifdef HAVE_TM1p1
          if(g.n_flavours==2) {
            for ( mu=0; mu<4; mu++ ) {
              for ( k=0; k<3; k++, j++ ) {
#ifndef BASIS4 
                rhs[j] = ((complex_double)vector1_in[i+2*(k+3*mu)] + I*(complex_double)vector1_in[i+2*(k+3*mu)+1]) * twisted_bc;
                rhs[j+6] = ((complex_double)vector2_in[i+2*(k+3*mu)] + I*(complex_double)vector2_in[i+2*(k+3*mu)+1]) * twisted_bc;

#else
                rhs[j] = ((complex_double)vector1_in[i+2*(k+3*(3-mu))] + I*(complex_double)vector1_in[i+2*(k+3*(3-mu))+1]) * twisted_bc;
                rhs[j+6] = ((complex_double)vector2_in[i+2*(k+3*(3-mu))] + I*(complex_double)vector2_in[i+2*(k+3*(3-mu))+1]) * twisted_bc;
#endif
                
#ifndef INIT_ONE_PREC
                if(g.mixed_precision==2) {
                  vtmp=cabs(rhs[j]);
                  if(vtmp > vmax)
                    vmax=vtmp;
                  if( vtmp > EPS_double && vtmp < vmin )
                    vmin=vtmp;
                  vtmp=cabs(rhs[j+6]);
                  if(vtmp > vmax)
                    vmax=vtmp;
                  if( vtmp > EPS_double && vtmp < vmin )
                    vmin=vtmp;
                }
              }
#endif
              if(mu%2)
                j+=6;
            }
          } else
#endif
            for ( mu=0; mu<4; mu++ )
              for ( k=0; k<3; k++, j++ ) {
#ifndef BASIS4 
                rhs[j] = ((complex_double)vector1_in[i+2*(k+3*mu)] + I*(complex_double)vector1_in[i+2*(k+3*mu)+1]) * twisted_bc;
#else
                rhs[j] = ((complex_double)vector1_in[i+2*(k+3*(3-mu))] + I*(complex_double)vector1_in[i+2*(k+3*(3-mu))+1]) * twisted_bc;
#endif
                
#ifndef INIT_ONE_PREC
                if( g.mixed_precision == 2 ) {
                  vtmp = cabs(rhs[j]);
                  if(vtmp > vmax)
                    vmax = vtmp;
                  if( vtmp > EPS_double && vtmp < vmin )
                    vmin = vtmp;
                }
              }
#endif
        }
      }
    }
  }
   
#ifndef INIT_ONE_PREC
    
  double gvmin, gvmax;
  if( g.mixed_precision == 2 ) {
    MPI_Allreduce(&vmin, &gvmin, 1, MPI_DOUBLE, MPI_MIN, g.comm_cart);
    MPI_Allreduce(&vmax, &gvmax, 1, MPI_DOUBLE, MPI_MAX, g.comm_cart);
  }
  
  //switching to double precision on the fine level
  if(g.mixed_precision==2 && gvmin/gvmax<EPS_float) {
    warning0("Changing solver precision on fine level due to rhs elements (min/max=%e)\n", vmin/vmax);
    precision_changed=1;
    g.mixed_precision=1;
    p = &(g.p);
    // storing pointer in x and b
    vb = p->b; 
    vx = p->x;
    p->b = g.p_MP.dp.b;
    p->x = g.p_MP.dp.x;
    p->tol = g.p_MP.dp.tol;
  } else precision_changed = 0;
#endif

  if ( n_shifts > 0 ) {
    ASSERT( even_shifts != NULL );
    ASSERT( odd_shifts != NULL );
  }
  if ( n_shifts > 1 ) {
    MALLOC( source, complex_double, l.inner_vector_size );
    MALLOC( solution, complex_double, l.inner_vector_size );
    if( _TYPE == _SOLVE_SQ || _TYPE == _SOLVE_SQ_ODD || _TYPE == _SOLVE_SQ_EVEN )
      MALLOC( solution2, complex_double, l.inner_vector_size );
  }
  
  for ( n = 0; n < n_shifts || n_shifts == 0; n++ ) {
    
    if ( n_shifts ) {
#ifdef HAVE_TM1p1
      if(g.n_flavours==2) {
        if( g.epsbar_ig5_even_shift != even_shifts[n] || g.epsbar_ig5_odd_shift != odd_shifts[n] ) {
          g.epsbar_ig5_even_shift = even_shifts[n];
          g.epsbar_ig5_odd_shift  =  odd_shifts[n];
          THREADED(threading[0]->n_core)
            epsbar_term_update( &l, threading[omp_get_thread_num()] );
          THREADED(threading[0]->n_core)
            finalize_operator_update( &l, threading[omp_get_thread_num()]);
        }        
      } else
#endif
        {
          if( g.mu_even_shift != even_shifts[n] || g.mu_odd_shift != odd_shifts[n] ) {
            g.mu_even_shift = even_shifts[n];
            g.mu_odd_shift  =  odd_shifts[n];
            THREADED(threading[0]->n_core)
              tm_term_update( g.mu, &l, threading[omp_get_thread_num()] );
            THREADED(threading[0]->n_core)
              finalize_operator_update( &l, threading[omp_get_thread_num()]);
          }
        }
    }

    switch(_TYPE) {
      
    case _SOLVE :
      if ( n ) {
        vector_copy( rhs, source );
        p->initial_guess_zero = 0;
      } else if ( n_shifts > 1 )
        vector_copy( source, rhs );
      
      solver( );
      break;

      
    case _SOLVE_SQ :
      if ( n ) {
        vector_copy( rhs, source );
        p->initial_guess_zero = 0;
      } else if ( n_shifts > 1 ) {
        THREADED(threading[0]->n_core) 
#ifdef HAVE_TM1p1
          if(g.n_flavours==2) 
            // sol = (D_h^{-1})*g5*tau1*(D_h^{-1})*g5*tau1*rhs
            tau1_gamma5_double( rhs, rhs, &l, threading[omp_get_thread_num()] );
          else
#endif
            // sol = (D_d^{-1})*g5*(D_u^{-1})*g5*rhs
            gamma5_double( rhs, rhs, &l, threading[omp_get_thread_num()] );
        vector_copy( source, rhs );
      }

      if( n )
        correct_guess( sol, solution, solution2, even_shifts[n]-even_shifts[n-1], odd_shifts[n]-odd_shifts[n-1]);
      solver( );
      if ( n < n_shifts-1 ) 
        vector_copy( solution, sol );
        
      THREADED(threading[0]->n_core) 
#ifdef HAVE_TM1p1
        if(g.n_flavours==2) 
          tau1_gamma5_double(rhs, sol, &l, threading[omp_get_thread_num()] );
        else
#endif
          gamma5_double(rhs, sol, &l, threading[omp_get_thread_num()] );

#ifdef HAVE_TM1p1
      if(g.n_flavours==2) 
        DDalphaAMG_change_epsbar_shift_sign( &tmp_status );
      else
#endif
        DDalphaAMG_change_mu_sign( &tmp_status );

      if( n )
        vector_copy( sol, solution2 );
      solver( );
      if ( n < n_shifts-1 ) 
        vector_copy( solution2, sol );
     
      // DDalphaAMG_change_mu_sign( &tmp_status );
      warning0("sign of mu changed during the inversion of squared operator\n");
      break;
      

    case _SOLVE_SQ_ODD :    
      if ( n ) {
        vector_copy( rhs, source );
        p->initial_guess_zero = 0;
      } else if ( n_shifts > 1 ) {
        THREADED(threading[0]->n_core)
#ifdef HAVE_TM1p1
          if(g.n_flavours==2) 
            // sol = (D_h^{-1})*g5*tau1*(D_h^{-1})*g5*tau1*rhs
            tau1_gamma5_set_even_to_zero_double(rhs, rhs, &l, threading[omp_get_thread_num()]);
          else
#endif
            // sol = (D_d^{-1})*g5*(D_u^{-1})*g5*rhs
            gamma5_set_even_to_zero_double(rhs, rhs, &l, threading[omp_get_thread_num()]);
       
        vector_copy( source, rhs );
      }

      if( n )
        correct_guess( sol, solution, solution2, even_shifts[n]-even_shifts[n-1], odd_shifts[n]-odd_shifts[n-1]);
      solver( );
      if ( n < n_shifts-1 ) 
        vector_copy( solution, sol );

      THREADED(threading[0]->n_core)
#ifdef HAVE_TM1p1
          if(g.n_flavours==2) 
            tau1_gamma5_set_even_to_zero_double(rhs, sol, &l, threading[omp_get_thread_num()]);
          else
#endif
            gamma5_set_even_to_zero_double(rhs, sol, &l, threading[omp_get_thread_num()]);

#ifdef HAVE_TM1p1
      if(g.n_flavours==2) 
        DDalphaAMG_change_epsbar_shift_sign( &tmp_status );
      else
#endif
        DDalphaAMG_change_mu_sign( &tmp_status );

      if( n )
        vector_copy( sol, solution2 );

      solver( );
      if ( n < n_shifts-1 ) 
        vector_copy( solution2, sol );

      // DDalphaAMG_change_mu_sign( &tmp_status );
      warning0("sign of mu changed during the inversion of squared operator\n");
      break;
      

    case _SOLVE_SQ_EVEN :    
      if ( n ) {
        vector_copy( rhs, source );
        p->initial_guess_zero = 0;
      } else if ( n_shifts > 1 ) {
        THREADED(threading[0]->n_core)
#ifdef HAVE_TM1p1
          if(g.n_flavours==2) 
            // sol = (D_h^{-1})*g5*tau1*(D_h^{-1})*g5*tau1*rhs
            tau1_gamma5_set_odd_to_zero_double(rhs, rhs, &l, threading[omp_get_thread_num()]);
          else
#endif
            // sol = (D_d^{-1})*g5*(D_u^{-1})*g5*rhs
            gamma5_set_odd_to_zero_double(rhs, rhs, &l, threading[omp_get_thread_num()]);

        vector_copy( source, rhs );
      }

      if( n )
        correct_guess( sol, solution, solution2, even_shifts[n]-even_shifts[n-1], odd_shifts[n]-odd_shifts[n-1]);
      solver( );
      if ( n < n_shifts-1 ) 
        vector_copy( solution, sol );

      THREADED(threading[0]->n_core)
#ifdef HAVE_TM1p1
        if(g.n_flavours==2) 
          tau1_gamma5_set_odd_to_zero_double(rhs, sol, &l, threading[omp_get_thread_num()]);
        else
#endif
          gamma5_set_odd_to_zero_double(rhs, sol, &l, threading[omp_get_thread_num()]);
       
#ifdef HAVE_TM1p1
      if(g.n_flavours==2) 
        DDalphaAMG_change_epsbar_shift_sign( &tmp_status );
      else
#endif
        DDalphaAMG_change_mu_sign( &tmp_status );
      
      if( n )
        vector_copy( sol, solution2 );
      solver( );
      if ( n < n_shifts-1 ) 
        vector_copy( solution2, sol );

      // DDalphaAMG_change_mu_sign( &tmp_status );
      warning0("sign of mu changed during the inversion of squared operator\n");
      break;
      

    case _PRECOND :
      THREADED(threading[0]->n_core)
        preconditioner( sol, NULL, rhs, _NO_RES, &l, threading[omp_get_thread_num()] );
      break;

      
    case _OPERATOR :
      THREADED(threading[0]->n_core)
        if ( g.mixed_precision == 2 ) {
          apply_operator_double( sol, rhs, &(g.p_MP.dp), &l, threading[omp_get_thread_num()] );
        } else {
          apply_operator_double( sol, rhs, &(g.p), &l, threading[omp_get_thread_num()] );
        }
      break;

      
    default :
      warning0("_TYPE not found in DDalphaAMG_driver. Returing vector in as vector out.");
      sol=rhs;
      break;
    }
    
    for (t=0, j=0; t<ll[T]; t++) {
      if (g.bc==_TWISTED) phase[T] = g.twisted_bc[T]*((double)sl[T]+t)/(double)gl[T];
      for (z=0; z<ll[Z]; z++) {
        if (g.bc==_TWISTED) phase[Z] = phase[T] + g.twisted_bc[Z]*((double)sl[Z]+z)/(double)gl[Z];
        for (y=0; y<ll[Y]; y++) {
          if (g.bc==_TWISTED) phase[Y] = phase[Z] + g.twisted_bc[Y]*((double)sl[Y]+y)/(double)gl[Y];
          for (x=0; x<ll[X]; x++) {
            if (g.bc==_TWISTED) {
              phase[X] = phase[Y] + g.twisted_bc[X]*((double)sl[X]+x)/(double)gl[X];
              twisted_bc = cexp(-I*phase[X]);
            } else
              twisted_bc = 1.;
            if(vector_index_fct!=NULL )
              i = vector_index_fct( t, z, y, x );
            else 
              i = 2*j;
            
#ifdef HAVE_TM1p1
            if(g.n_flavours==2) {
              for ( mu=0; mu<4; mu++ ) {
                for ( k=0; k<3; k++, j++ ) {
                  tmp1 = sol[j] * twisted_bc;
                  tmp2 = sol[j+6] * twisted_bc;
#ifndef BASIS4 
                  vector1_out[n][i+2*(k+3*mu)]   = creal(tmp1);
                  vector1_out[n][i+2*(k+3*mu)+1] = cimag(tmp1);
                  vector2_out[n][i+2*(k+3*mu)]   = creal(tmp2);
                  vector2_out[n][i+2*(k+3*mu)+1] = cimag(tmp2);
#else
                  vector1_out[n][i+2*(k+3*(3-mu))]   = creal(tmp1);
                  vector1_out[n][i+2*(k+3*(3-mu))+1] = cimag(tmp1);
                  vector2_out[n][i+2*(k+3*(3-mu))]   = creal(tmp2);
                  vector2_out[n][i+2*(k+3*(3-mu))+1] = cimag(tmp2);
#endif   
                }
                if(mu%2)
                  j+=6;
              }
            } else
#endif
              for ( mu=0; mu<4; mu++ )
                for ( k=0; k<3; k++, j++ ) {
                  tmp1 = sol[j] * twisted_bc;
#ifndef BASIS4 
                  vector1_out[n][i+2*(k+3*mu)]   = creal(tmp1);
                  vector1_out[n][i+2*(k+3*mu)+1] = cimag(tmp1);
#else
                  vector1_out[n][i+2*(k+3*(3-mu))]   = creal(tmp1);
                  vector1_out[n][i+2*(k+3*(3-mu))+1] = cimag(tmp1);
#endif   
                }
          }
        }
      }
    }
    
  }

  if ( n_shifts > 0 ) {
    FREE( source, complex_double, l.inner_vector_size );
    FREE( solution, complex_double, l.inner_vector_size );
    if( _TYPE == _SOLVE_SQ || _TYPE == _SOLVE_SQ_ODD || _TYPE == _SOLVE_SQ_EVEN )
      FREE( solution2, complex_double, l.inner_vector_size );
    p->initial_guess_zero = 1;
  }

  
#ifndef INIT_ONE_PREC
  if (precision_changed) {
    g.mixed_precision=2;
    // recovering pointer from x and b
    p->b = vb; 
    p->x = vx;
  }
#endif
    
  if ( g.norm_res <= tol || _TYPE == _OPERATOR || _TYPE == _PRECOND )
    mg_status->success = 1;
  mg_status->info = g.norm_res;
  t1 = MPI_Wtime();
  mg_status->time = t1-t0;
  mg_status->coarse_time = g.coarse_time;
  mg_status->iter_count = g.iter_count;
  mg_status->coarse_iter_count = g.coarse_iter_count;
  
}

static inline void set_n_flavours( int n) {

#ifdef HAVE_TM1p1
  THREADED(threading[0]->n_core)
    data_layout_n_flavours( n, &l, threading[omp_get_thread_num()] );
#else
  if( n==2 )
      error0("For DDalphaAMG_solve_doublet_*, HAVE_TM1p1 flag required\n");
#endif
    
}

void DDalphaAMG_solve( double *vector_out, double *vector_in, double tol, DDalphaAMG_status *mg_status )
{
  DDalphaAMG_driver( vector_out, vector_in, NULL, NULL, tol, mg_status, _SOLVE );
}

void DDalphaAMG_solve_doublet( double *vector1_out, double *vector1_in,
                               double *vector2_out, double *vector2_in,
                               double tol, DDalphaAMG_status *mg_status )
{
  set_n_flavours( 2 );
  DDalphaAMG_driver( vector1_out, vector1_in, vector2_out, vector2_in, tol, mg_status, _SOLVE );
  set_n_flavours( 1 );
}

void DDalphaAMG_solve_ms_doublet( double **vector1_out, double *vector1_in,
                                  double **vector2_out, double *vector2_in,
                                  double  *even_shifts, double *odd_shifts, int n_shifts,
                                  double tol, DDalphaAMG_status *mg_status )
{
  set_n_flavours( 2 );
  DDalphaAMG_ms_driver( vector1_out, vector1_in, vector2_out, vector2_in, even_shifts, odd_shifts, n_shifts,
                        tol, mg_status, _SOLVE );
  set_n_flavours( 1 );
}

void DDalphaAMG_solve_squared( double *vector_out, double *vector_in, double tol, DDalphaAMG_status *mg_status )
{
  DDalphaAMG_driver( vector_out, vector_in, NULL, NULL, tol, mg_status, _SOLVE_SQ );
}

void DDalphaAMG_solve_doublet_squared( double *vector1_out, double *vector1_in,
                                       double *vector2_out, double *vector2_in,
                                       double tol, DDalphaAMG_status *mg_status )
{
  set_n_flavours( 2 );
  DDalphaAMG_driver( vector1_out, vector1_in, vector2_out, vector2_in, tol, mg_status, _SOLVE_SQ );
  set_n_flavours( 1 );
}

void DDalphaAMG_solve_ms_doublet_squared( double **vector1_out, double *vector1_in,
                                          double **vector2_out, double *vector2_in,
                                          double  *even_shifts, double *odd_shifts, int n_shifts,
                                          double tol, DDalphaAMG_status *mg_status )
{
  set_n_flavours( 2 );
  DDalphaAMG_ms_driver( vector1_out, vector1_in, vector2_out, vector2_in, even_shifts, odd_shifts, n_shifts,
                        tol, mg_status, _SOLVE_SQ );
  set_n_flavours( 1 );
}

void DDalphaAMG_solve_squared_odd( double *vector_out, double *vector_in, double tol, DDalphaAMG_status *mg_status )
{
  DDalphaAMG_driver( vector_out, vector_in, NULL, NULL, tol, mg_status, _SOLVE_SQ_ODD );
}

void DDalphaAMG_solve_doublet_squared_odd( double *vector1_out, double *vector1_in,
                                           double *vector2_out, double *vector2_in,
                                           double tol, DDalphaAMG_status *mg_status )
{
  set_n_flavours( 2 );
  DDalphaAMG_driver( vector1_out, vector1_in, vector2_out, vector2_in, tol, mg_status, _SOLVE_SQ_ODD );
  set_n_flavours( 1 );
}

void DDalphaAMG_solve_ms_doublet_squared_odd( double **vector1_out, double *vector1_in,
                                              double **vector2_out, double *vector2_in,
                                              double  *even_shifts, double *odd_shifts, int n_shifts,
                                              double tol, DDalphaAMG_status *mg_status )
{
  set_n_flavours( 2 );
  DDalphaAMG_ms_driver( vector1_out, vector1_in, vector2_out, vector2_in, even_shifts, odd_shifts, n_shifts,
                        tol, mg_status, _SOLVE_SQ_ODD );
  set_n_flavours( 1 );
}

void DDalphaAMG_solve_squared_even( double *vector_out, double *vector_in, double tol, DDalphaAMG_status *mg_status )
{
  DDalphaAMG_driver( vector_out, vector_in, NULL, NULL, tol, mg_status, _SOLVE_SQ_EVEN );
}

void DDalphaAMG_solve_doublet_squared_even( double *vector1_out, double *vector1_in,
                                            double *vector2_out, double *vector2_in,
                                            double tol, DDalphaAMG_status *mg_status )
{
  set_n_flavours( 2 );
  DDalphaAMG_driver( vector1_out, vector1_in, vector2_out, vector2_in, tol, mg_status, _SOLVE_SQ_EVEN );
  set_n_flavours( 1 );
}

void DDalphaAMG_solve_ms_doublet_squared_even( double **vector1_out, double *vector1_in,
                                               double **vector2_out, double *vector2_in,
                                               double  *even_shifts, double *odd_shifts, int n_shifts,
                                               double tol, DDalphaAMG_status *mg_status )
{
  set_n_flavours( 2 );
  DDalphaAMG_ms_driver( vector1_out, vector1_in, vector2_out, vector2_in, even_shifts, odd_shifts, n_shifts,
                        tol, mg_status, _SOLVE_SQ_EVEN );
  set_n_flavours( 1 );
}

void DDalphaAMG_apply_operator( double *vector_out, double *vector_in, DDalphaAMG_status *mg_status ) {
  DDalphaAMG_driver( vector_out, vector_in, NULL, NULL, 0, mg_status, _OPERATOR );
}

void DDalphaAMG_apply_operator_doublet( double *vector1_out, double *vector1_in,
                                        double *vector2_out, double *vector2_in, DDalphaAMG_status *mg_status )
{
  set_n_flavours( 2 );
  DDalphaAMG_driver( vector1_out, vector1_in, vector2_out, vector2_in, 0, mg_status, _OPERATOR );
  set_n_flavours( 1 );
}

void DDalphaAMG_preconditioner( double *vector_out, double *vector_in, DDalphaAMG_status * mg_status ) {
  DDalphaAMG_driver( vector_out, vector_in, NULL, NULL, 0, mg_status, _PRECOND );
}

void DDalphaAMG_preconditioner_doublet( double *vector1_out, double *vector1_in,
                                        double *vector2_out, double *vector2_in, DDalphaAMG_status *mg_status )
{
  set_n_flavours( 2 );
  DDalphaAMG_driver( vector1_out, vector1_in, vector2_out, vector2_in, 0, mg_status, _PRECOND );
  set_n_flavours( 1 );
}

void DDalphaAMG_free( void ) {
  method_free( &l );
  g.setup_flag = 0;
}


void DDalphaAMG_finalize( void ) {

  finalize_common_thread_data(commonthreaddata);
  finalize_no_threading(no_threading);
  for(int i=0; i<g.num_openmp_processes; i++) {
    FREE( threading[i], struct Thread, 1);
  }
  FREE( threading, struct Thread *, g.num_openmp_processes);
  
  if (g.setup_flag)
    method_free( &l );
  method_finalize( &l );
  
}

MPI_Comm DDalphaAMG_get_communicator( void ){
  MPI_Comm comm;
  MPI_Comm_dup( g.comm_cart, &comm);
  return comm;
}

void DDalphaAMG_read_configuration( double *gauge_field, char *filename, int format, DDalphaAMG_status *mg_status ) {

  double plaq;

  if(format==1)
    lime_read_conf( gauge_field, filename, &plaq );
  else
    read_conf( gauge_field, filename, &plaq, &l );

}

void DDalphaAMG_read_vector( double *vector_in, char *filename, int format, DDalphaAMG_status *mg_status ){

  if(format==1)
    lime_read_vector( vector_in, filename );
  else
    vector_io( vector_in, filename, _READ, &l );

}

void DDalphaAMG_write_vector( double *vector_out, char *filename, int format, DDalphaAMG_status *mg_status ){

  if(format==1)
    lime_write_vector( vector_out, filename );
  else
    vector_io( vector_out, filename, _WRITE, &l );

}

void DDalphaAMG_define_vector_const( double *vector, double re, double im ) {

  THREADED(threading[0]->n_core)
  if(vector!=NULL){
    int start, end;
    compute_core_start_end( 0, l.inner_vector_size, &start, &end, &l, threading[omp_get_thread_num()]);
    vector_double_define( (vector_double) vector, re+I*im, start, end, &l );
  }
  else {
    warning0("Vector NULL when calling DDalphaAMG_define_vector_const!");
  }
}

void DDalphaAMG_define_vector_rand( double *vector ) {

  THREADED(threading[0]->n_core)
  if(vector!=NULL){
    int start, end;
    compute_core_start_end( 0, l.inner_vector_size, &start, &end, &l, threading[omp_get_thread_num()]);
    vector_double_define_random( (vector_double) vector, start, end, &l );
  }
  else {
    warning0("Vector NULL when calling DDalphaAMG_define_vector_const!");
  }

}

double DDalphaAMG_vector_norm( double *vector ) {

  double norm = 0;
  THREADED(threading[0]->n_core)
  if(vector!=NULL){
    norm = global_norm_double( (vector_double) vector, 0, l.inner_vector_size, &l, threading[omp_get_thread_num()] );
   }
  else {
    warning0("Vector NULL when calling DDalphaAMG_define_vector_const!");
  }

  return norm;
}

void DDalphaAMG_vector_saxpy( double *vector_out, double a, double *x, double *y ) {

  THREADED(threading[0]->n_core)
  if(vector_out!=NULL && x!=NULL && y!=NULL){
    int start, end;
    compute_core_start_end( 0, l.inner_vector_size, &start, &end, &l, threading[omp_get_thread_num()]);
    vector_double_saxpy( (vector_double) vector_out, (vector_double) x, (vector_double) y, a, start, end, &l );
  }
  else {
    warning0("Vector NULL when calling DDalphaAMG_define_vector_const!");
  }

}

void DDalphaAMG_test_routine( DDalphaAMG_status *mg_status ) {
  
  double t0, t1;
  t0 = MPI_Wtime();

  printf00("\n");
  THREADED(threading[0]->n_core)
  test_routine( &l, threading[omp_get_thread_num()]);

  if (g.test < 1e-5)
    mg_status->success = 1;
  else
    mg_status->success = 0;    
  mg_status->info = g.test; //highest error
  t1 = MPI_Wtime();
  mg_status->time = t1-t0;
  mg_status->coarse_time = 0;
  mg_status->iter_count = 0;
  mg_status->coarse_iter_count = 0;
}

void DDalphaAMG_get_parameters( DDalphaAMG_parameters *mg_params ){

  int i, j;
   
  mg_params->method = g.method;
  mg_params->interpolation = g.interpolation;
  mg_params->mixed_precision = g.mixed_precision;
  mg_params->kcycle_tolerance = g.kcycle_tol;
  mg_params->coarse_tolerance = g.coarse_tol;
  mg_params->smoother_iterations = g.post_smooth_iter[0];
  mg_params->conf_index_fct = conf_index_fct;
  mg_params->vector_index_fct = vector_index_fct;
  mg_params->kappa = 0.5/(g.m0 + 4.);
#ifdef HAVE_TM
  mg_params->mu = g.mu;
  mg_params->mu_odd_shift = g.mu_odd_shift;
  mg_params->mu_even_shift = g.mu_even_shift;
#else
  mg_params->mu = 0;
  mg_params->mu_odd_shift = 0;
  mg_params->mu_even_shift = 0;
#endif
#ifdef HAVE_TM1p1
  mg_params->epsbar = g.epsbar;
  mg_params->epsbar_ig5_odd_shift = g.epsbar_ig5_odd_shift;
  mg_params->epsbar_ig5_even_shift = g.epsbar_ig5_even_shift;
#else
  mg_params->epsbar = 0;
  mg_params->epsbar_ig5_odd_shift = 0;
  mg_params->epsbar_ig5_even_shift = 0;
#endif
  mg_params->print = g.print;
  
  for( i=0; i<g.num_levels; i++ ) {
    for( j=0; j<4; j++ )
      mg_params->block_lattice[i][j] = g.block_lattice[i][j];
    if( i<g.num_levels-1 )
      mg_params->mg_basis_vectors[i] = g.num_eig_vect[i];
    mg_params->setup_iterations[i] = g.setup_iter[i];
    mg_params->mu_factor[i] = g.mu_factor[i];
#ifdef HAVE_TM1p1
    mg_params->epsbar_factor[i] = g.epsbar_factor[i];
#else
    mg_params->epsbar_factor[i] = 1;
#endif
  }  
}
