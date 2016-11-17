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

#ifndef OPTIMIZED_INTERPOLATION_SETUP_PRECISION  
void gram_schmidt_on_aggregates_PRECISION( vector_PRECISION *V, const int num_vect, level_struct *l, struct Thread *threading ) 
#else
void gram_schmidt_on_aggregates_PRECISION( complex_PRECISION *V, const int num_vec, level_struct *l, struct Thread *threading )
#endif
{

  PROF_PRECISION_START( _GRAM_SCHMIDT_ON_AGGREGATES, threading );

#ifndef OPTIMIZED_INTERPOLATION_SETUP_PRECISION  
  SYNC_CORES(threading)
  SYNC_HYPERTHREADS(threading)
  int i, j, k, k1, k2, num_aggregates = l->s_PRECISION.num_aggregates,
      aggregate_size = l->inner_vector_size / num_aggregates, offset = l->num_lattice_site_var/2;
      
  complex_PRECISION alpha1, alpha2;
  vector_PRECISION v_pt1, v_pt2;
  PRECISION norm1, norm2;
      
  for ( j=threading->n_thread*threading->core+threading->thread; j<num_aggregates; j+=threading->n_thread*threading->n_core ) {
    for ( k1=0; k1<num_vect; k1++ ) {
      v_pt1 = V[k1] + j*aggregate_size;
      
      for ( k2=0; k2<k1; k2++ ) {
        v_pt2 = V[k2] + j*aggregate_size;
        alpha1 = 0; alpha2 = 0;
        // V[k1] -= <V[k2],V[k1]> V[k2] | 2*j-th and 2*j+1-st aggregate
        for ( i=0; i<aggregate_size; ) {
          for ( k=0; k<offset; k++, i++ )
            alpha1 += conj_PRECISION(v_pt2[i]) * v_pt1[i];
          for ( k=0; k<offset; k++, i++ )
            alpha2 += conj_PRECISION(v_pt2[i]) * v_pt1[i];
        }
        for ( i=0; i<aggregate_size; ) {
          for ( k=0; k<offset; k++, i++ )
            v_pt1[i] -=  alpha1 * v_pt2[i];
          for ( k=0; k<offset; k++, i++ )
            v_pt1[i] -=  alpha2 * v_pt2[i];
        }
      }
      
      norm1 = 0; norm2 = 0;
      // V[k1] = V[k1]/norm(V[k1]) | 2*j-th and 2*j+1-st aggregate    
      for ( i=0; i<aggregate_size; ) {
        for ( k=0; k<offset; k++, i++ )
          norm1 += NORM_SQUARE_PRECISION(v_pt1[i]);
        for ( k=0; k<offset; k++, i++ )
          norm2 += NORM_SQUARE_PRECISION(v_pt1[i]);
      }
      norm1 = 1/sqrt(norm1); norm2 = 1/sqrt(norm2);
      for ( i=0; i<aggregate_size; ) {
        for ( k=0; k<offset; k++, i++ )
          v_pt1[i] =  norm1 * creal_PRECISION(v_pt1[i]) + I*norm1* cimag_PRECISION(v_pt1[i]);
        for ( k=0; k<offset; k++, i++ )
          v_pt1[i] =  norm2 * creal_PRECISION(v_pt1[i]) + I*norm2* cimag_PRECISION(v_pt1[i]);
      }
    }
  }
  SYNC_HYPERTHREADS(threading)
  SYNC_CORES(threading)

#else

  // Here the ortogonalization is done on a operator layout.
  // the block version has some optimizations which are correct only on the fine grid
  if(l->depth == 0)
    block_gram_schmidt_PRECISION( V, num_vec, l, threading );  
  else
    aggregate_gram_schmidt_PRECISION( V, num_vec, l, threading );
#endif

  PROF_PRECISION_STOP( _GRAM_SCHMIDT_ON_AGGREGATES, 1, threading );
}


void gram_schmidt_PRECISION( vector_PRECISION *V, complex_PRECISION *buffer, const int begin, const int n, level_struct *l, struct Thread *threading ) {
  
  // NOTE: only thread safe, if "buffer" is the same buffer for all threads belonging to a common MPI process
  START_MASTER(threading)
  PROF_PRECISION_START( _LA );
  END_MASTER(threading)
  SYNC_CORES(threading)
  
  PRECISION beta;
  int i, j, start, end;
  
  compute_core_start_end_custom( 0, l->inner_vector_size, &start, &end, l, threading, l->num_lattice_site_var );
  
  for ( i=begin; i<n; i++ ) {
    
    complex_PRECISION tmp[i];
    process_multi_inner_product_PRECISION( i, tmp, V, V[i], 0, l->inner_vector_size, l, threading );
    SYNC_CORES(threading)
    START_MASTER(threading)
    for ( j=0; j<i; j++ ) {
      buffer[j] = tmp[j];
    }
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
    
    if ( i>0 ) {
      START_MASTER(threading)
      PROF_PRECISION_START( _ALLR );
      MPI_Allreduce( buffer, buffer+n, i, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
      PROF_PRECISION_STOP( _ALLR, 1 );
      END_MASTER(threading)
      SYNC_MASTER_TO_ALL(threading)
    }
    
    for( j=0; j<i; j++ ) {
      vector_PRECISION_saxpy( V[i], V[i], V[j], -(buffer+n)[j], 0, l->inner_vector_size, l, threading );
      SYNC_CORES(threading)
    }
    
    SYNC_CORES(threading)
      
    beta = global_norm_PRECISION( V[i], 0, l->inner_vector_size, l, threading );
    SYNC_MASTER_TO_ALL(threading)
    vector_PRECISION_real_scale( V[i], V[i], creal(1.0/beta), start, end, l );
    SYNC_CORES(threading)
  }
  
  START_MASTER(threading)
  PROF_PRECISION_STOP( _LA, 1 );
  END_MASTER(threading)
  SYNC_CORES(threading)
}
