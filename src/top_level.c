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

#include "main.h"

void rhs_define( vector_double rhs, level_struct *l, struct Thread *threading ) {
  
  // no hyperthreading here
  if(threading->thread != 0)
    return;

  if ( g.rhs == 0 ) {
    vector_double_define_real( rhs, 1, 0, l->inner_vector_size, l, threading );
    START_MASTER(threading)
    if ( g.print > 0 ) printf0("rhs = ones\n");
    END_MASTER(threading)
  } else if ( g.rhs == 1 )  {
    vector_double_define_zero( rhs, 0, l->inner_vector_size, l, threading );
    if ( g.my_rank == 0 ) {
      START_LOCKED_MASTER(threading)
      rhs[0] = 1.0;
      END_LOCKED_MASTER(threading)
    }
    START_MASTER(threading)
    if ( g.print > 0 ) printf0("rhs = first unit vector\n");
    END_MASTER(threading)
  } else if ( g.rhs == 2 ) {
    // this would yield different results if we threaded it, so we don't
    vector_double_define_random( rhs, 0, l->inner_vector_size, l, threading );
    START_MASTER(threading)
    if ( g.print > 0 ) printf0("rhs = random\n");
    END_MASTER(threading)
  } else if ( g.rhs == 3 ) {
    vector_double_define_zero( rhs, 0, l->inner_vector_size, l, threading );
  } else {
    ASSERT( g.rhs >= 0 && g.rhs <= 4 );
  }
    
}


int wilson_driver( vector_double solution, vector_double source, level_struct *l, struct Thread *threading ) {
  
  int iter = 0, start = threading->start_index[l->depth], end = threading->end_index[l->depth];
  
  vector_double rhs = (g.mixed_precision==2 && g.method >= 0)?g.p_MP.dp.b:g.p.b;
  vector_double sol = (g.mixed_precision==2 && g.method >= 0)?g.p_MP.dp.x:g.p.x;

#ifdef WILSON_BENCHMARK
  START_MASTER(threading)
  prof_init( l );
  END_MASTER(threading)
  double t = -MPI_Wtime();
  double t_min = 1000;
  for ( int i=0; i<100; i++ ) {
    double tmp_t = -MPI_Wtime();
#endif
  
  vector_double_copy( rhs, source, start, end, l );  
  if ( g.method == -1 ) {
    cgn_double( &(g.p), l, threading );
  } else if ( g.mixed_precision == 2 ) {
    iter = fgmres_MP( &(g.p_MP), l, threading );
  } else {
    iter = fgmres_double( &(g.p), l, threading );
  }
  vector_double_copy( solution, sol, start, end, l );
#ifdef WILSON_BENCHMARK
    tmp_t += MPI_Wtime();
    if ( tmp_t < t_min )
      t_min = tmp_t;
  }
  t +=MPI_Wtime();
  START_MASTER(threading)
  printf0("average over 100 solves: %lf seconds\n", t/100 );
  printf0("minimum out of 100 solves: %lf seconds\n", t_min );
  prof_print( l );
  END_MASTER(threading)
#endif
  
  return iter;
}


void solve( vector_double solution, vector_double source, level_struct *l, struct Thread *threading ) {
  
  if ( g.vt.evaluation ) {
    vector_double rhs = g.mixed_precision==2?g.p_MP.dp.b:g.p.b;
    // this would yield different results if we threaded it, so we don't
    vector_double_define_random( rhs, 0, l->inner_vector_size, l, threading );
    START_LOCKED_MASTER(threading)
    scan_var( &(g.vt), l );
    END_LOCKED_MASTER(threading)
  } else {
    wilson_driver( solution, source, l, threading );
  }
}


void solve_driver( level_struct *l, struct Thread *threading ) {
  
  vector_double solution = NULL, source = NULL;
  double minus_twisted_bc[4], norm;
 
  if(g.bc==2)
    for ( int i=0; i<4; i++ )
      minus_twisted_bc[i] = -1*g.twisted_bc[i];
  
#ifdef HAVE_TM1p1
  if( g.epsbar != 0 || g.epsbar_ig5_odd_shift != 0 || g.epsbar_ig5_odd_shift != 0 ) { 
    data_layout_n_flavours( 2, l, threading );
    printf0("inverting doublet operator\n");
  }
#endif
  PUBLIC_MALLOC( solution, complex_double, l->inner_vector_size );
  PUBLIC_MALLOC( source, complex_double, l->inner_vector_size );

  rhs_define( source, l, threading );

  if(g.bc==2)
    apply_twisted_bc_to_vector_double( source, source, g.twisted_bc, l);

  norm = global_norm_double( source, 0, l->inner_vector_size, l, threading );
  printf0("source vector norm: %le\n",norm);

#ifdef HAVE_TM1p1
  if( g.n_flavours == 1 )
#endif
#ifdef HAVE_TM
  if ( g.mu + g.mu_odd_shift != 0.0 || g.mu + g.mu_even_shift != 0.0 )
    if(g.downprop) {
      
      START_MASTER(threading)  
      printf0("\n\n+--------------------------- up ---------------------------+\n\n");
      END_MASTER(threading)

      solve( solution, source, l, threading );    
      
      if(g.bc==2)
     apply_twisted_bc_to_vector_double( solution, solution, minus_twisted_bc, l);
      
      START_LOCKED_MASTER(threading)  
      printf0("\n\n+-------------------------- down --------------------------+\n\n");
      g.mu*=-1;
      g.mu_odd_shift*=-1;
      g.mu_even_shift*=-1;
      END_LOCKED_MASTER(threading)
  
      tm_term_update( g.mu, l, threading );
      finalize_operator_update( l, threading );
    } 
#endif

  solve( solution, source, l, threading );

  if(g.bc==2)
    apply_twisted_bc_to_vector_double( solution, solution, minus_twisted_bc, l);

  norm = global_norm_double( solution, 0, l->inner_vector_size, l, threading );
  printf0("solution vector norm: %le\n",norm);

  PUBLIC_FREE( solution, complex_double, l->inner_vector_size );
  PUBLIC_FREE( source, complex_double, l->inner_vector_size );

#ifdef HAVE_TM1p1
  if( g.n_flavours == 2 ) 
    data_layout_n_flavours( 1, l, threading );
#endif
}

