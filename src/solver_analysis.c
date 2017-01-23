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


void test_routine( level_struct *l, struct Thread *threading ) {
  
  if ( g.method > 0 ) {
    if ( g.mixed_precision ) {
      operator_float_test_routine( &(l->s_float.op), l, threading );
      if ( g.method > 0 && g.method < 4 ) schwarz_float_mvm_testfun( &(l->s_float), l, threading );
      if ( g.method > 0 && g.method < 4 && g.odd_even ) block_oddeven_float_test( l, threading );
    } else {
      operator_double_test_routine( &(l->s_double.op), l, threading );
      if ( g.method > 0 && g.method < 4 ) schwarz_double_mvm_testfun( &(l->s_double), l, threading );
      if ( g.method > 0 && g.method < 4 && g.odd_even ) block_oddeven_double_test( l, threading );
    }
    
    if ( g.interpolation ) {
      if ( g.mixed_precision )
        coarse_operator_float_test_routine( l, threading );
      else
        coarse_operator_double_test_routine( l, threading );
    }
  }

  START_LOCKED_MASTER(threading)
  printf0("\n");
  prof_init( l );
  END_LOCKED_MASTER(threading)

  if ( g.restart > 0 )
    rhs_define( g.p.b, l, threading );
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
