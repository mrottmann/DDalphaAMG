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

void prof_PRECISION_init( level_struct *l ) {

/*********************************************************************************
* Initializes the profiling struct by specifying the name of every entry.
*********************************************************************************/
  
  if ( l != NULL ) {
    for ( int i=0; i<_NUM_PROF; i++ ) {
      l->prof_PRECISION.time[i] = 0.0;
      l->prof_PRECISION.count[i] = 0.0;
      l->prof_PRECISION.flop[i] = 0.0;
    }
    
    double level_ratio = 1;
    for ( int mu=0; mu<4; mu++ )
      level_ratio *= (double)g.global_lattice[l->depth][mu]/(double)g.global_lattice[0][mu];
    
    sprintf( l->prof_PRECISION.name[_GIP], "global inner product, PRECISION" );
    l->prof_PRECISION.flop[_GIP] = level_ratio*l->num_lattice_site_var*8.0;
    sprintf( l->prof_PRECISION.name[_PIP], "process inner product, PRECISION" );
    l->prof_PRECISION.flop[_PIP] = level_ratio*l->num_lattice_site_var*8.0;
    sprintf( l->prof_PRECISION.name[_LA2], "2 flop vector operations, PRECISION" );
    l->prof_PRECISION.flop[_LA2] = level_ratio*l->num_lattice_site_var*2.0;
    sprintf( l->prof_PRECISION.name[_LA6], "6 flop vector operations, PRECISION" );
    l->prof_PRECISION.flop[_LA6] = level_ratio*l->num_lattice_site_var*6.0;
    sprintf( l->prof_PRECISION.name[_LA8], "8 flop vector operations, PRECISION" );
    l->prof_PRECISION.flop[_LA8] = level_ratio*l->num_lattice_site_var*8.0;
    sprintf( l->prof_PRECISION.name[_LA], "other vector operations, PRECISION" );
    sprintf( l->prof_PRECISION.name[_GRAM_SCHMIDT], "Gram-Schmidt, PRECISION" );
    sprintf( l->prof_PRECISION.name[_GRAM_SCHMIDT_ON_AGGREGATES], "Gram-Schmidt on aggregates, PRECISION" );
    sprintf( l->prof_PRECISION.name[_CPY], "copy operations, PRECISION" );
    sprintf( l->prof_PRECISION.name[_SET], "set value operations, PRECISION" );
    sprintf( l->prof_PRECISION.name[_PR], "interpolation and restriction, PRECISION" );
    l->prof_PRECISION.flop[_PR] = level_ratio*l->num_lattice_site_var*8.0*(l->num_lattice_site_var/2);
    sprintf( l->prof_PRECISION.name[_SC], "self coupling, PRECISION" );
    l->prof_PRECISION.flop[_SC] = (l->depth==0)?552.0:level_ratio*SQUARE(l->num_lattice_site_var)*8.0;
    sprintf( l->prof_PRECISION.name[_NC], "neighbor coupling, PRECISION" );
    l->prof_PRECISION.flop[_NC] = (l->depth==0)?1368.0:level_ratio*8.0*SQUARE(l->num_lattice_site_var)*8.0;
    sprintf( l->prof_PRECISION.name[_SM], "smoother, PRECISION" );
    double ncflops = l->prof_PRECISION.flop[_SC];
    for ( int mu=0; mu<4; mu++ )
      ncflops += (l->prof_PRECISION.flop[_NC]/4.0)*((double)(l->block_lattice[mu]-1)/(double)l->block_lattice[mu]);
    l->prof_PRECISION.flop[_SM] = ncflops * (double)(g.odd_even?l->block_iter+1:l->block_iter);
    l->prof_PRECISION.flop[_SM] += (l->prof_PRECISION.flop[_NC] + l->prof_PRECISION.flop[_SC]);
    sprintf( l->prof_PRECISION.name[_OP_COMM], "operator comm init, PRECISION" );
    sprintf( l->prof_PRECISION.name[_OP_IDLE], "operator comm wait, PRECISION" );
    sprintf( l->prof_PRECISION.name[_ALLR], "allreduces, PRECISION" );
    sprintf( l->prof_PRECISION.name[_GD_COMM], "data re-distribution comm init, PRECISION" );
    sprintf( l->prof_PRECISION.name[_GD_IDLE], "data re-distribution comm wait, PRECISION" );
    sprintf( l->prof_PRECISION.name[_SM1], "smoother - pt 1, res no comm, PRECISION" );
    sprintf( l->prof_PRECISION.name[_SM2], "smoother - pt 2, solve no comm, PRECISION" );
    sprintf( l->prof_PRECISION.name[_SM3], "smoother - pt 3, res comm, PRECISION" );
    sprintf( l->prof_PRECISION.name[_SM4], "smoother - pt 4, solve comm, PRECISION" );
    
    sprintf( l->prof_PRECISION.name[_SMALL1], "Hessenberg: qr update PRECISION" );
    sprintf( l->prof_PRECISION.name[_SMALL2], "Hessenberg: bkwd subst PRECISION" );
  }
}


double prof_PRECISION_print( level_struct *l ) {
  double flop = 0;
  for ( int i=0; i<_NUM_PROF; i++ )
    if ( l->prof_PRECISION.count[i] > 0 ) {
      if ( l->prof_PRECISION.count[i] > 9999999 )
        printf0("| %37s: %8.2le(%7.1le) |\n", l->prof_PRECISION.name[i], l->prof_PRECISION.time[i], l->prof_PRECISION.count[i] );
      else
        printf0("| %37s: %8.2le(%7d) |\n", l->prof_PRECISION.name[i], l->prof_PRECISION.time[i], (int)l->prof_PRECISION.count[i] );
      flop += (double)l->prof_PRECISION.count[i] * l->prof_PRECISION.flop[i];
    }
  return flop;
}


void fine_level_PRECISION_alloc( level_struct *l ) {
  
  int n = 8;
  
  MALLOC( l->vbuf_PRECISION[0], complex_PRECISION, n*l->vector_size );  
  for ( int i=1; i<n; i++ )
    l->vbuf_PRECISION[i] = l->vbuf_PRECISION[0] + i*l->vector_size;
  MALLOC( l->p_PRECISION.b, complex_PRECISION, 2*l->inner_vector_size );
  l->p_PRECISION.x = l->p_PRECISION.b + l->inner_vector_size;
}


void fine_level_PRECISION_free( level_struct *l ) {
  
  int n = 8;
  
  FREE( l->vbuf_PRECISION[0], complex_PRECISION, n*l->vector_size );  
  for ( int i=1; i<n; i++ )
    l->vbuf_PRECISION[i] = NULL;
  FREE( l->p_PRECISION.b, complex_PRECISION, 2*l->inner_vector_size );
  l->p_PRECISION.x = NULL;
}


void next_level_PRECISION_setup( level_struct *l ) {
    
  prof_float_init( l->next_level );
  prof_double_init( l->next_level );
  gathering_PRECISION_next_level_init( &(l->next_level->gs_PRECISION), l );  
  gathering_PRECISION_setup( &(l->next_level->gs_PRECISION), l->next_level );
  
  if ( !l->idle ) {
    coarsening_index_table_PRECISION_alloc( &(l->is_PRECISION), l );
    coarsening_index_table_PRECISION_define( &(l->is_PRECISION), &(l->s_PRECISION), l );

    if ( l->level == 1 && !l->next_level->idle ) {
      fgmres_PRECISION_struct_alloc( g.coarse_iter, g.coarse_restart, l->next_level->vector_size, g.coarse_tol, 
                                     _COARSE_GMRES, _NOTHING, NULL,
                                     g.method==6?(g.odd_even?g5D_coarse_apply_schur_complement_PRECISION:g5D_apply_coarse_operator_PRECISION)
                                     :(g.odd_even?coarse_apply_schur_complement_PRECISION:apply_coarse_operator_PRECISION),
                                     &(l->next_level->p_PRECISION), l->next_level );
    } else {
      if ( g.kcycle ) {
        fgmres_PRECISION_struct_alloc( g.kcycle_restart, g.kcycle_max_restart, l->next_level->vector_size, g.kcycle_tol, 
                                       _K_CYCLE, _RIGHT, vcycle_PRECISION,
                                       g.method==6?g5D_apply_coarse_operator_PRECISION:apply_coarse_operator_PRECISION,
                                       &(l->next_level->p_PRECISION), l->next_level );
      } else {
        MALLOC( l->next_level->p_PRECISION.b, complex_PRECISION, 2*l->next_level->vector_size );
        l->next_level->p_PRECISION.x = l->next_level->p_PRECISION.b + l->next_level->vector_size;
        l->next_level->p_PRECISION.shift = 0;
        l->next_level->p_PRECISION.v_start = 0;
        l->next_level->p_PRECISION.v_end = l->inner_vector_size;
      }
    }

    int i, n = (l->next_level->level>0)?6:4;
    MALLOC( l->next_level->vbuf_PRECISION[0], complex_PRECISION, n*l->next_level->vector_size );
    for ( i=1; i<n; i++ )
      l->next_level->vbuf_PRECISION[i] = l->next_level->vbuf_PRECISION[0] + i*l->next_level->vector_size;
  }
}


void next_level_PRECISION_free( level_struct *l ) {
  
  coarse_grid_correction_PRECISION_free( l );
  
  if ( !l->idle ) {
    if ( ( l->level == 1 && !l->next_level->idle ) || g.kcycle ) {
      fgmres_PRECISION_struct_free( &(l->next_level->p_PRECISION), l->next_level );
    } else {
      FREE( l->next_level->p_PRECISION.b, complex_PRECISION, 2*l->next_level->vector_size );
    }
  
    int i, n = (l->next_level->level>0)?6:4;  
    for ( i=1; i<n; i++)
      l->next_level->vbuf_PRECISION[i] = NULL;
    FREE( l->next_level->vbuf_PRECISION[0], complex_PRECISION, n*l->next_level->vector_size );
    coarsening_index_table_PRECISION_free( &(l->is_PRECISION), l );
  }
  
  gathering_PRECISION_free( &(l->next_level->gs_PRECISION), l->next_level );
}


void level_PRECISION_init( level_struct *l ) {

  for ( int i=0; i<9; i++ )
    l->vbuf_PRECISION[i] = NULL;
  
  operator_PRECISION_init( &(l->op_PRECISION) );
  operator_PRECISION_init( &(l->oe_op_PRECISION) );
  schwarz_PRECISION_init( &(l->s_PRECISION), l );
  interpolation_PRECISION_struct_init( &(l->is_PRECISION) );
  fgmres_PRECISION_struct_init( &(l->p_PRECISION) );
  fgmres_PRECISION_struct_init( &(l->sp_PRECISION) );
}


void vcycle_timing_PRECISION( int n, level_struct *l, struct Thread *threading ) {
  
  ASSERT( g.mixed_precision );
  vector_PRECISION v1 = NULL, v2 = NULL;
  double t0=0, t1=0;
  PUBLIC_MALLOC( v1, complex_PRECISION, l->inner_vector_size );
  PUBLIC_MALLOC( v2, complex_PRECISION, l->inner_vector_size );

  START_LOCKED_MASTER(threading)
  vector_PRECISION_define_random( v2, 0, l->inner_vector_size, l );
  END_LOCKED_MASTER(threading)
  
  START_MASTER(threading)
  t0 = MPI_Wtime();
  END_MASTER(threading)
  for ( int i=0; i<n; i++ ) {
    vcycle_PRECISION( v1, NULL, v2, _NO_RES, l, threading );
  }
  START_MASTER(threading)
  t1 = MPI_Wtime();
  printf0("100 v-cycles: %le seconds\n", t1-t0 );
  END_MASTER(threading)

  START_LOCKED_MASTER(threading)
  PUBLIC_FREE( v1, complex_PRECISION, l->inner_vector_size );
  PUBLIC_FREE( v2, complex_PRECISION, l->inner_vector_size );
  END_LOCKED_MASTER(threading)
}
