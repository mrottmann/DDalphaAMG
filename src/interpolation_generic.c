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

#if ( !defined( SSE ) || !defined( INTERPOLATION_OPERATOR_LAYOUT_OPTIMIZED_PRECISION ) )

void interpolation_PRECISION_alloc( level_struct *l ) {
  
  int k, n = l->num_eig_vect;
  
  MALLOC( l->is_PRECISION.eigenvalues, complex_PRECISION, n );
  MALLOC( l->is_PRECISION.test_vector, complex_PRECISION*, n );
  MALLOC( l->is_PRECISION.interpolation, complex_PRECISION*, n );
  l->is_PRECISION.interpolation[0] = NULL;
  MALLOC_HUGEPAGES( l->is_PRECISION.interpolation[0], complex_PRECISION, n*l->vector_size, 64 );
  for ( k=1; k<n; k++ )
    l->is_PRECISION.interpolation[k] = l->is_PRECISION.interpolation[0] + k*l->vector_size;
  MALLOC( l->is_PRECISION.operator, complex_PRECISION, n*l->inner_vector_size );
  l->is_PRECISION.test_vector[0] = NULL;
  MALLOC_HUGEPAGES( l->is_PRECISION.test_vector[0], complex_PRECISION, n*l->inner_vector_size, 64 );
  for ( k=1; k<n; k++ ) {
    l->is_PRECISION.test_vector[k] = l->is_PRECISION.test_vector[0] + k*l->inner_vector_size;
  }
}


void interpolation_PRECISION_dummy_alloc( level_struct *l ) {
  
  MALLOC( l->is_PRECISION.test_vector, complex_PRECISION*, l->num_eig_vect );
  MALLOC( l->is_PRECISION.interpolation, complex_PRECISION*, l->num_eig_vect );  
}


void interpolation_PRECISION_dummy_free( level_struct *l ) {
  
  FREE( l->is_PRECISION.test_vector, complex_PRECISION*, l->num_eig_vect );
  FREE( l->is_PRECISION.interpolation, complex_PRECISION*, l->num_eig_vect );  
}


void interpolation_PRECISION_free( level_struct *l ) {
  
  int n = l->num_eig_vect;
  
  FREE_HUGEPAGES( l->is_PRECISION.test_vector[0], complex_PRECISION, n*l->inner_vector_size );
  FREE( l->is_PRECISION.eigenvalues, complex_PRECISION, n );
  FREE( l->is_PRECISION.test_vector, complex_PRECISION*, n );
  FREE_HUGEPAGES( l->is_PRECISION.interpolation[0], complex_PRECISION, n*l->vector_size );
  FREE( l->is_PRECISION.interpolation, complex_PRECISION*, n );
  FREE( l->is_PRECISION.operator, complex_PRECISION, n*l->inner_vector_size );

}


void define_interpolation_PRECISION_operator( complex_PRECISION **interpolation, level_struct *l, struct Thread *threading ) {
  
  int j, num_eig_vect = l->num_eig_vect;
  complex_PRECISION *operator = l->is_PRECISION.operator;

  int start = threading->start_index[l->depth];
  int end = threading->end_index[l->depth];
      
  SYNC_CORES(threading)
  operator += start*num_eig_vect;
  for ( int i=start; i<end; i++ )
    for ( j=0; j<num_eig_vect; j++ ) {
      *operator = interpolation[j][i];
      operator++;
    }
  SYNC_CORES(threading)
}


void interpolate_PRECISION( vector_PRECISION phi, vector_PRECISION phi_c, level_struct *l, struct Thread *threading ) {
  
  PROF_PRECISION_START( _PR, threading );
  int i, j, k, k1, k2, num_aggregates = l->is_PRECISION.num_agg, num_eig_vect = l->num_eig_vect, sign = 1,
      num_parent_eig_vect = l->num_lattice_site_var/2, aggregate_sites = l->num_inner_lattice_sites / num_aggregates;
  complex_PRECISION *operator = l->is_PRECISION.operator, *phi_pt = phi,
                    *phi_c_pt = l->next_level->gs_PRECISION.transfer_buffer;
                    
  START_LOCKED_MASTER(threading)
  vector_PRECISION_distribute( phi_c_pt, phi_c, l->next_level );
  END_LOCKED_MASTER(threading)
  SYNC_HYPERTHREADS(threading)
  
  for ( i=threading->n_thread*threading->core + threading->thread; i<num_aggregates; i+=threading->n_core*threading->n_thread ) {
    phi_pt   = phi + i*2*num_parent_eig_vect*aggregate_sites;
    phi_c_pt = l->next_level->gs_PRECISION.transfer_buffer + i*2*num_eig_vect;
    operator = l->is_PRECISION.operator + i*2*num_eig_vect*num_parent_eig_vect*aggregate_sites;
    for ( k=0; k<aggregate_sites; k++ ) {
      for ( k1=0; k1<2; k1++ ) {
        for ( k2=0; k2<num_parent_eig_vect; k2++ ) {
          for ( j=0; j<num_eig_vect; j++ ) {
            *phi_pt += phi_c_pt[j] * (*operator);
            operator++;
          }
          phi_pt++;
        }
        phi_c_pt += sign*num_eig_vect; sign*=-1;
      }
    }
  }
    
  PROF_PRECISION_STOP( _PR, 1, threading );

  SYNC_HYPERTHREADS(threading)
}


void interpolate3_PRECISION( vector_PRECISION phi, vector_PRECISION phi_c, level_struct *l, struct Thread *threading ) {
  
  PROF_PRECISION_START( _PR, threading );
  int i, j, k, k1, k2, num_aggregates = l->is_PRECISION.num_agg, num_eig_vect = l->num_eig_vect,
      num_parent_eig_vect = l->num_lattice_site_var/2, aggregate_sites = l->num_inner_lattice_sites / num_aggregates;
  complex_PRECISION *operator = l->is_PRECISION.operator, *phi_pt = phi,
                    *phi_c_pt = l->next_level->gs_PRECISION.transfer_buffer;
  
  START_LOCKED_MASTER(threading)
  vector_PRECISION_distribute( phi_c_pt, phi_c, l->next_level );
  END_LOCKED_MASTER(threading)
  SYNC_HYPERTHREADS(threading)
  
  for ( i=threading->n_thread*threading->core + threading->thread; i<num_aggregates; i+=threading->n_core*threading->n_thread ) {
    phi_pt   = phi + i*2*num_parent_eig_vect*aggregate_sites;
    phi_c_pt = l->next_level->gs_PRECISION.transfer_buffer + i*2*num_eig_vect;
    int sign = 1;
    operator = l->is_PRECISION.operator + i*2*num_eig_vect*num_parent_eig_vect*aggregate_sites;
    for ( k=0; k<aggregate_sites; k++ ) {
      for ( k1=0; k1<2; k1++ ) {
        for ( k2=0; k2<num_parent_eig_vect; k2++ ) {
          *phi_pt = phi_c_pt[0] * (*operator);
          operator++;
          for ( j=1; j<num_eig_vect; j++ ) {
            *phi_pt += phi_c_pt[j] * (*operator);
            operator++;
          }
          phi_pt++;
        }
        phi_c_pt += sign*num_eig_vect; sign*=-1;
      }
    }
  }
  PROF_PRECISION_STOP( _PR, 1, threading );

  SYNC_HYPERTHREADS(threading)
}


void restrict_PRECISION( vector_PRECISION phi_c, vector_PRECISION phi, level_struct *l, struct Thread *threading ) {
  
  SYNC_CORES(threading)
  SYNC_HYPERTHREADS(threading)

  PROF_PRECISION_START( _PR, threading );
  int i, j, k, k1, k2, num_aggregates = l->is_PRECISION.num_agg, num_eig_vect = l->num_eig_vect, sign = 1,
      num_parent_eig_vect = l->num_lattice_site_var/2, aggregate_sites = l->num_inner_lattice_sites / num_aggregates;
  complex_PRECISION *operator = l->is_PRECISION.operator, *phi_pt = phi,
                    *phi_c_pt = l->next_level->gs_PRECISION.transfer_buffer;

  for ( i=threading->n_thread*threading->core + threading->thread; i<num_aggregates; i+=threading->n_core*threading->n_thread ) {   
    phi_pt   = phi + i*2*num_parent_eig_vect*aggregate_sites;
    phi_c_pt = l->next_level->gs_PRECISION.transfer_buffer + i*2*num_eig_vect;
    operator = l->is_PRECISION.operator + i*2*num_eig_vect*num_parent_eig_vect*aggregate_sites;

    for ( j=0; j<2*num_eig_vect; j++ )
      phi_c_pt[j] = 0;
    
    for ( k=0; k<aggregate_sites; k++ ) {
      for ( k1=0; k1<2; k1++ ) {
        for ( k2=0; k2<num_parent_eig_vect; k2++ ) {
          for ( j=0; j<num_eig_vect; j++ ) {
            phi_c_pt[j] += *phi_pt * conj_PRECISION(*operator);
            operator++;
          }
          phi_pt++;
        }
        phi_c_pt += sign*num_eig_vect; sign*=-1;
      }
    }
  }
  
  SYNC_HYPERTHREADS(threading)
  START_LOCKED_MASTER(threading)
  vector_PRECISION_gather( phi_c, l->next_level->gs_PRECISION.transfer_buffer, l->next_level );
  END_LOCKED_MASTER(threading)
  PROF_PRECISION_STOP( _PR, 1, threading );
}

#endif
