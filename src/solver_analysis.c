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


void test_routine( level_struct *l, struct Thread *threading ) {

  if ( g.method >= 0 ) {
    START_MASTER(threading)
    g.test = 0;
    if ( l->depth == 0 ) {
#ifdef HAVE_TM1p1
      if( g.n_flavours==2 )
        printf0("\nRunning tests with D = TM doublet operator:\n");
      else
#endif
#ifdef HAVE_TM
        printf0("\nRunning tests with D = TM Wilson operator:\n");
#else
      printf0("\nRunning tests with D = Wilson operator:\n");
#endif
    }
    END_MASTER(threading)
    if ( g.mixed_precision ) {
      operator_float_test_routine( &(l->s_float.op), l, threading );
      if ( g.method > 0 && g.method < 4 ) schwarz_float_mvm_testfun( &(l->s_float), l, threading );
      if ( g.method > 0 && g.method < 4 && g.odd_even ) block_oddeven_float_test( l, threading );
    } else {
      operator_double_test_routine( &(l->s_double.op), l, threading );
      if ( g.method > 0 && g.method < 4 ) schwarz_double_mvm_testfun( &(l->s_double), l, threading );
      if ( g.method > 0 && g.method < 4 && g.odd_even ) block_oddeven_double_test( l, threading );
    }
    
    if ( g.interpolation && g.method > 0 ) {
      if ( g.mixed_precision )
        coarse_operator_float_test_routine( l, threading );
      else
        coarse_operator_double_test_routine( l, threading );
    }
    START_MASTER(threading)
    if (g.test < 1e-5)
      printf0("TESTS passed, highest error %e < 1e-5\n", g.test);
    else
      warning0("some TESTS not passed, highest error %e > 1e-5\n", g.test);
    printf0("\n");
    END_MASTER(threading)
  }


#ifdef HAVE_TM1p1
  if( g.n_flavours==1 &&
      (g.epsbar != 0 || g.epsbar_ig5_odd_shift != 0 || g.epsbar_ig5_odd_shift != 0) ) {
    
    if ( g.method >= 0 ) {
      START_MASTER(threading)
      g.test = 0;
      printf0("Running tests with D = TM doublet operator:\n");
      END_MASTER(threading)

      data_layout_n_flavours( 2, l, threading );
      
      if ( g.mixed_precision ) 
        two_flavours_test_float( &(l->s_float.op), l, threading ); 
      else 
        two_flavours_test_double( &(l->s_double.op), l, threading ); 
      
      if ( g.mixed_precision ) {
        operator_float_test_routine( &(l->s_float.op), l, threading );
        if ( g.method > 0 && g.method < 4 ) schwarz_float_mvm_testfun( &(l->s_float), l, threading );
        if ( g.method > 0 && g.method < 4 && g.odd_even ) block_oddeven_float_test( l, threading );
      } else {
        operator_double_test_routine( &(l->s_double.op), l, threading );
        if ( g.method > 0 && g.method < 4 ) schwarz_double_mvm_testfun( &(l->s_double), l, threading );
        if ( g.method > 0 && g.method < 4 && g.odd_even ) block_oddeven_double_test( l, threading );
      }
      
      if ( g.interpolation  && g.method > 0 ) {
        if ( g.mixed_precision )
          coarse_operator_float_test_routine( l, threading );
        else
          coarse_operator_double_test_routine( l, threading );
      }
     
      START_MASTER(threading)
      if (g.test < 1e-5)
        printf0("TESTS passed, highest error %e < 1e-5\n", g.test);
      else
        warning0("some TESTS not passed, highest error %e > 1e-5\n", g.test);
      printf0("\n");
      END_MASTER(threading)
      
      data_layout_n_flavours( 1, l, threading );
    }
  }
#endif

}


void prof_init( level_struct *l ) {
  if ( l->depth == 0 ) { g.coarse_time=0; g.coarse_iter_count=0; }
  prof_double_init( l );
  prof_float_init( l );
  if ( l->next_level != NULL )
    prof_init( l->next_level );
}


double prof_print( level_struct *l ) {
  double flop = 0;
#ifdef PROFILING
  
  if ( l != NULL && g.print > 0 ) {
    if ( l->depth == 0 ) printf0("\n+----------------------------------------------------------+\n");
    if ( l->depth == 0 ) printf0("| solver profiling                                         |\n");
    printf0("+----------------------------------------------------------+\n");
    printf0("| depth: %3d / level: %3d                time    ( count ) |\n", l->depth, l->level );
    printf0("+----------------------------------------------------------+\n");
    flop += prof_double_print( l );
    flop += prof_float_print( l );
    flop += prof_print( l->next_level );
    if ( l->depth == 0 ) {
      int *ll = l->local_lattice;
      printf0("+----------------------------------------------------------+\n");
      printf0("| flop/lattice site: %9.2le                             |\n", flop );
      printf0("| flop/s/MPIprocess: %9.2le                             |\n",
              (flop/g.total_time)*ll[T]*ll[Z]*ll[Y]*ll[X] );
      printf0("+----------------------------------------------------------+\n\n");
    }
  }
#endif
  return flop;
}
