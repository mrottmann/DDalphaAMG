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

#include "sse_float_intrinsic.h"
#include "sse_linalg.h"
#include "sse_linalg_PRECISION.h"

#ifndef OPTIMIZED_LINALG_PRECISION
complex_PRECISION global_inner_product_PRECISION( vector_PRECISION phi, vector_PRECISION psi, int start, int end, level_struct *l, struct Thread *threading ) {
  
  PROF_PRECISION_START( _GIP, threading );
  complex_PRECISION local_alpha = 0, global_alpha = 0;

  int thread_start;
  int thread_end;
  compute_core_start_end(start, end, &thread_start, &thread_end, l, threading);
  
  SYNC_CORES(threading)
  
  VECTOR_FOR( int i=thread_start, i<thread_end, local_alpha += conj_PRECISION(phi[i])*psi[i], i++, l );

  // sum over cores
  START_NO_HYPERTHREADS(threading)
  ((complex_PRECISION *)threading->workspace)[threading->core] = local_alpha;
  END_NO_HYPERTHREADS(threading)
  // master sums up all results
  SYNC_CORES(threading)
  START_MASTER(threading)
  for(int i=1; i<threading->n_core; i++)
    ((complex_PRECISION *)threading->workspace)[0] += ((complex_PRECISION *)threading->workspace)[i];
  local_alpha = ((complex_PRECISION *)threading->workspace)[0];
  END_MASTER(threading)
  
  if ( g.num_processes > 1 ) {
    START_MASTER(threading)
    PROF_PRECISION_START( _ALLR );
    MPI_Allreduce( &local_alpha, &global_alpha, 1, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
    PROF_PRECISION_STOP( _ALLR, 1 );
    ((complex_PRECISION *)threading->workspace)[0] = global_alpha;
    END_MASTER(threading)
    // all threads need the result of the norm
    SYNC_MASTER_TO_ALL(threading)
    global_alpha = ((complex_PRECISION *)threading->workspace)[0];
    PROF_PRECISION_STOP( _GIP, (double)(end-start)/(double)l->inner_vector_size, threading );
    return global_alpha;
  } else {
    // all threads need the result of the norm
    SYNC_MASTER_TO_ALL(threading)
    local_alpha = ((complex_PRECISION *)threading->workspace)[0];
    PROF_PRECISION_STOP( _GIP, (double)(end-start)/(double)l->inner_vector_size, threading );
    return local_alpha;
  }
}
#endif


complex_PRECISION process_inner_product_PRECISION( vector_PRECISION phi, vector_PRECISION psi, int start, int end, level_struct *l, struct Thread *threading ) {
  
  PROF_PRECISION_START( _PIP, threading );
  int i;
  complex_PRECISION local_alpha = 0;
  
  SYNC_CORES(threading)
  
  THREADED_VECTOR_FOR( i, start, end, local_alpha += conj_PRECISION(phi[i])*psi[i], i++, l, threading );

  START_NO_HYPERTHREADS(threading)
  ((complex_PRECISION *)threading->workspace)[threading->core] = local_alpha;
  END_NO_HYPERTHREADS(threading)
  // master sums up all results
  SYNC_CORES(threading)
  START_MASTER(threading)
  for(int i=1; i<threading->n_core; i++)
    ((complex_PRECISION *)threading->workspace)[0] += ((complex_PRECISION *)threading->workspace)[i];
  END_MASTER(threading)
  // all threads need the result of the norm
  SYNC_MASTER_TO_ALL(threading)
  local_alpha = ((complex_PRECISION *)threading->workspace)[0];

  PROF_PRECISION_STOP( _PIP, (double)(end-start)/(double)l->inner_vector_size, threading );

  return local_alpha;
}


#if !defined( OPTIMIZED_LINALG_PRECISION ) && !defined( __MIC__ )
void process_multi_inner_product_PRECISION( int count, complex_PRECISION *results, vector_PRECISION *phi, vector_PRECISION psi,
    int start, int end, level_struct *l, struct Thread *threading ) {

  PROF_PRECISION_START( _PIP, threading );
  int i;
  for(int c=0; c<count; c++)
    results[c] = 0.0;

  int thread_start;
  int thread_end;

  SYNC_CORES(threading)
  if ( l->depth == 0 ) {
    compute_core_start_end_custom(start, end, &thread_start, &thread_end, l, threading, 12);
    for(int c=0; c<count; c++)
      for ( i=thread_start; i<thread_end; )
        FOR12( results[c] += conj_PRECISION(phi[c][i])*psi[i]; i++; )
  } else {
#ifdef _M10TV
    compute_core_start_end_custom(start, end, &thread_start, &thread_end, l, threading, 20);
    for(int c=0; c<count; c++)
      for ( i=thread_start; i<thread_end; )
        FOR20( results[c] += conj_PRECISION(phi[c][i])*psi[i]; i++; )
#else
    compute_core_start_end_custom(start, end, &thread_start, &thread_end, l, threading, 2);
    for(int c=0; c<count; c++)
      for ( i=thread_start; i<thread_end; )
        FOR2( results[c] += conj_PRECISION(phi[c][i])*psi[i]; i++; )
#endif
  }

  START_NO_HYPERTHREADS(threading)
  ((complex_PRECISION **)threading->workspace)[threading->core] = results;
  END_NO_HYPERTHREADS(threading)
  // master sums up all results
  SYNC_CORES(threading)
  START_MASTER(threading)
  for(int c=0; c<count; c++)
    for(int i=1; i<threading->n_core; i++)
      ((complex_PRECISION **)threading->workspace)[0][c] += ((complex_PRECISION **)threading->workspace)[i][c];
  END_MASTER(threading)
  // all threads need the result of the norm
  SYNC_MASTER_TO_ALL(threading)
  for(int c=0; c<count; c++)
    results[c] = ((complex_PRECISION **)threading->workspace)[0][c];

  PROF_PRECISION_STOP( _PIP, (double)(end-start)/(double)l->inner_vector_size, threading );
}
#endif


complex_PRECISION local_xy_over_xx_PRECISION( vector_PRECISION phi, vector_PRECISION psi, int start, int end, level_struct *l  ) {
  
  complex_PRECISION numerator = 0.0; PRECISION denominator = 0.0;
  
  VECTOR_FOR( int i=start, i<end, numerator += conj_PRECISION(phi[i])*psi[i]; denominator += NORM_SQUARE_PRECISION(phi[i]), i++, l );
  
  if ( abs_PRECISION(denominator) < EPS_PRECISION ) {
    return 0.0;
  }
  
  return numerator/denominator;
}

#ifndef OPTIMIZED_LINALG_PRECISION
PRECISION global_norm_PRECISION( vector_PRECISION x, int start, int end, level_struct *l, struct Thread *threading ) {
  
  PROF_PRECISION_START( _GIP, threading );
  
  PRECISION local_alpha = 0, global_alpha = 0;

  int thread_start;
  int thread_end;
  compute_core_start_end(start, end, &thread_start, &thread_end, l, threading);
  
  SYNC_CORES(threading)
  
  VECTOR_FOR( int i=thread_start, i<thread_end, local_alpha += NORM_SQUARE_PRECISION(x[i]), i++, l );

  // sum over cores
  START_NO_HYPERTHREADS(threading)
  ((PRECISION *)threading->workspace)[threading->core] = local_alpha;
  END_NO_HYPERTHREADS(threading)
  // master sums up all results
  SYNC_CORES(threading)
  START_MASTER(threading)
  for(int i=1; i<threading->n_core; i++)
    ((PRECISION *)threading->workspace)[0] += ((PRECISION *)threading->workspace)[i];
  local_alpha = ((PRECISION *)threading->workspace)[0];
  END_MASTER(threading)

  if ( g.num_processes > 1 ) {
    START_MASTER(threading)
    PROF_PRECISION_START( _ALLR );
    MPI_Allreduce( &local_alpha, &global_alpha, 1, MPI_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
    PROF_PRECISION_STOP( _ALLR, 1 );
    ((PRECISION *)threading->workspace)[0] = global_alpha;
    END_MASTER(threading)
    // all threads need the result of the norm
    SYNC_MASTER_TO_ALL(threading)
    global_alpha = ((PRECISION *)threading->workspace)[0];
    PROF_PRECISION_STOP( _GIP, (double)(end-start)/(double)l->inner_vector_size, threading );
    return (PRECISION)sqrt((double)global_alpha);
  } else {
    // all threads need the result of the norm
    SYNC_MASTER_TO_ALL(threading)
    local_alpha = ((PRECISION *)threading->workspace)[0];
    PROF_PRECISION_STOP( _GIP, (double)(end-start)/(double)l->inner_vector_size, threading );
    return (PRECISION)sqrt((double)local_alpha);
  }
}
#endif

PRECISION process_norm_PRECISION( vector_PRECISION x, int start, int end, level_struct *l, struct Thread *threading ) {
     
  int i;
  PRECISION local_alpha = 0;
  PROF_PRECISION_START( _PIP, threading );
  
  SYNC_CORES(threading)
  
  THREADED_VECTOR_FOR( i, start, end, local_alpha += NORM_SQUARE_PRECISION(x[i]), i++, l, threading );

  START_NO_HYPERTHREADS(threading)
  ((PRECISION *)threading->workspace)[threading->core] = local_alpha;
  END_NO_HYPERTHREADS(threading)
  // master sums up all results
  SYNC_CORES(threading)
  START_MASTER(threading)
  for(int i=1; i<threading->n_core; i++)
    ((PRECISION *)threading->workspace)[0] += ((PRECISION *)threading->workspace)[i];
  END_MASTER(threading)
  // all threads need the result of the norm
  SYNC_MASTER_TO_ALL(threading)
  local_alpha = ((PRECISION *)threading->workspace)[0];

  PROF_PRECISION_STOP( _PIP, (double)(end-start)/(double)l->inner_vector_size, threading );

  return (PRECISION)sqrt((double)local_alpha);
}


void vector_PRECISION_plus( vector_PRECISION z, vector_PRECISION x, vector_PRECISION y, int start, int end, level_struct *l ) {
  
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
  PROF_PRECISION_START( _LA2 );
  
  VECTOR_FOR( int i=start, i<end, z[i] = x[i] + y[i], i++, l );
  
  if(thread == 0 && start != end)
  PROF_PRECISION_STOP( _LA2, (double)(end-start)/(double)l->inner_vector_size );
}


void vector_PRECISION_minus( vector_PRECISION z, vector_PRECISION x, vector_PRECISION y, int start, int end, level_struct *l ) {
  
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
  PROF_PRECISION_START( _LA2 );

  VECTOR_FOR( int i=start, i<end, z[i] = x[i] - y[i], i++, l );
  
  if(thread == 0 && start != end)
  PROF_PRECISION_STOP( _LA2, (double)(end-start)/(double)l->inner_vector_size );
}

#ifndef OPTIMIZED_LINALG_PRECISION
void vector_PRECISION_scale( vector_PRECISION z, vector_PRECISION x, complex_PRECISION alpha, int start, int end, level_struct *l ) {
  
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
  PROF_PRECISION_START( _LA6 );
  
  VECTOR_FOR( int i=start, i<end, z[i] = alpha*x[i], i++, l );
  
  if(thread == 0 && start != end)
  PROF_PRECISION_STOP( _LA6, (double)(end-start)/(double)l->inner_vector_size );
}
#endif


void vector_PRECISION_real_scale( vector_PRECISION z, vector_PRECISION x, complex_PRECISION alpha,
                                  int start, int end, level_struct *l ) {
  
  PRECISION *r_z = (PRECISION*)z, *r_x = (PRECISION*)x, r_alpha = creal_PRECISION(alpha);
  int r_start = 2*start, r_end = 2*end;
  
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
  PROF_PRECISION_START( _LA2 );
  
  REAL_VECTOR_FOR( int i=r_start, i<r_end, r_z[i] = r_alpha*r_x[i], i++, l );
  
  if(thread == 0 && start != end)
  PROF_PRECISION_STOP( _LA2, (double)(end-start)/(double)l->inner_vector_size );
}


void vector_PRECISION_copy( vector_PRECISION z, vector_PRECISION x, int start, int end, level_struct *l ) {
  
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
  PROF_PRECISION_START( _CPY );
  
  VECTOR_FOR( int i=start, i<end, z[i] = x[i], i++, l );
  
  if(thread == 0 && start != end)
  PROF_PRECISION_STOP( _CPY, (double)(end-start)/(double)l->inner_vector_size );
}

#ifndef OPTIMIZED_LINALG_PRECISION
void vector_PRECISION_saxpy( vector_PRECISION z, vector_PRECISION x, vector_PRECISION y, complex_PRECISION alpha, int start, int end, level_struct *l ) {
  
  int thread = omp_get_thread_num();
  if (thread == 0 && start != end )
  PROF_PRECISION_START( _LA8 );
  
  VECTOR_FOR( int i=start, i<end, z[i] = x[i] + alpha*y[i], i++, l );
  
  if( thread == 0 && start != end )
  PROF_PRECISION_STOP( _LA8, (double)(end-start)/(double)l->inner_vector_size );
}
#endif

#ifndef OPTIMIZED_LINALG_PRECISION
void vector_PRECISION_multi_saxpy( vector_PRECISION z, vector_PRECISION *V, complex_PRECISION *alpha,
                               int sign, int count, int start, int end, level_struct *l ) {
  
  int thread = omp_get_thread_num();
  if (thread == 0 && start != end )
  PROF_PRECISION_START( _LA8 );
  
  complex_PRECISION alpha_signed[count];
  for ( int c=0; c<count; c++ ) {
    alpha_signed[c] = sign*alpha[c];
  }
  
  for ( int c=0; c<count; c++ ) {
    for ( int i=start; i<end; ) {
      FOR12( z[i] += V[c][i]*alpha_signed[c]; i++; )
    }
  }
  
  if( thread == 0 && start != end )
  PROF_PRECISION_STOP( _LA8, (PRECISION)(count) );
}
#endif

void vector_PRECISION_projection( vector_PRECISION z, vector_PRECISION v, int k, vector_PRECISION *W, complex_PRECISION *diag, 
                                  int orthogonal, level_struct *l, Thread *threading ) {
  
  int j, start, end;
  
  compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
                    
  vector_PRECISION v_tmp = NULL, *W_tmp = NULL;
  complex_PRECISION ip[k], ip_buffer[2*k];      
  
  MALLOC( v_tmp, complex_PRECISION, l->inner_vector_size );
  vector_PRECISION_define(v_tmp, 0, 0, l->inner_vector_size, l );
  
  MALLOC( W_tmp, complex_PRECISION*, k );
  W_tmp[0] = NULL; 
  MALLOC( W_tmp[0], complex_PRECISION, k*l->inner_vector_size );
  for ( j = 1; j<k; j++ )
    W_tmp[j] = W_tmp[0]+j*l->inner_vector_size;
  
  for ( j=0; j<k; j++ ) {
   vector_PRECISION_scale( W_tmp[j], W[j], diag[j], 0, l->inner_vector_size, l );
  }
  process_multi_inner_product_PRECISION( k, ip, W_tmp, v, 0, l->inner_vector_size, l, threading );
  
  START_MASTER(threading)
  for ( j=0; j<k; j++ ) {
    ip_buffer[j] = ip[j];
  }
  MPI_Allreduce( ip_buffer, ip_buffer+k, k, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)  
  
  vector_PRECISION_multi_saxpy( v_tmp, W_tmp, ip_buffer+k, 1, k, 0, l->inner_vector_size, l );
   
  if (orthogonal) 
    vector_PRECISION_minus( z, v, v_tmp, 0, l->inner_vector_size, l );
  else
    vector_PRECISION_copy( z, v_tmp, 0, l->inner_vector_size, l );
  
  FREE( v_tmp, complex_PRECISION, l->inner_vector_size );
  FREE( W_tmp[0], complex_PRECISION, k*l->inner_vector_size );
  FREE( W_tmp, complex_PRECISION*, k );
}

void gram_schmidt_on_aggregates_PRECISION( vector_PRECISION *V, const int num_vect, level_struct *l, struct Thread *threading ) {
  
  PROF_PRECISION_START( _GRAM_SCHMIDT_ON_AGGREGATES, threading );
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
  PROF_PRECISION_STOP( _GRAM_SCHMIDT_ON_AGGREGATES, 1, threading );
}


void spinwise_PRECISION_skalarmultiply( vector_PRECISION eta1, vector_PRECISION eta2, vector_PRECISION phi, complex_PRECISION alpha,
                                        int start, int end, level_struct *l ) {
  
  PROF_PRECISION_START( _LA6 );  
  for ( int i=start; i<end; ) {
    FOR6( eta1[i] = alpha*phi[i]; eta2[i] = _COMPLEX_PRECISION_ZERO; i++; )
    FOR6( eta2[i] = alpha*phi[i]; eta1[i] = _COMPLEX_PRECISION_ZERO; i++; )
  }
  PROF_PRECISION_STOP( _LA6, 1 );
}


void set_boundary_PRECISION( vector_PRECISION phi, complex_PRECISION alpha, level_struct *l, struct Thread *threading ) {
  
  PROF_PRECISION_START( _SET, threading );
  int i;
  
  SYNC_CORES(threading)
  
  THREADED_VECTOR_FOR( i, l->inner_vector_size, l->vector_size, phi[i] = alpha, i++, l, threading );
  
  SYNC_CORES(threading)
  PROF_PRECISION_STOP( _SET, (double)(l->vector_size-l->inner_vector_size)/(double)l->inner_vector_size, threading );
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
      vector_PRECISION_saxpy( V[i], V[i], V[j], -(buffer+n)[j], start, end, l );
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


#if !defined( SSE ) || !defined( GRAM_SCHMIDT_VECTORIZED_PRECISION )
void setup_gram_schmidt_PRECISION_compute_dots(
    complex_PRECISION *thread_buffer, vector_PRECISION *V, int count, int offset,
    int start, int end, level_struct *l, struct Thread *threading) {

  int thread_start;
  int thread_end;
  int cache_block_size = 12*64;
  complex_PRECISION tmp[cache_block_size];

  for(int i=0; i<2*offset; i++)
    thread_buffer[i] = 0.0;

  SYNC_CORES(threading)
  compute_core_start_end_custom(start, end, &thread_start, &thread_end, l, threading, cache_block_size);
  
  for ( int i=thread_start; i<thread_end; i+=cache_block_size) {
    coarse_gamma5_PRECISION( tmp, V[count]+i, 0, cache_block_size, l );
    for ( int j=0; j<count; j++ ) {
      for ( int k=0; k<cache_block_size; k++) {
        thread_buffer[j]   += conj_PRECISION(V[j][i+k])*V[count][i+k];
        thread_buffer[j+offset] += conj_PRECISION(V[j][i+k])*tmp[k];
      }
    }
  }

  START_NO_HYPERTHREADS(threading)
  ((complex_PRECISION **)threading->workspace)[threading->core] = thread_buffer;
  END_NO_HYPERTHREADS(threading)
  // master sums up all results
  SYNC_CORES(threading)
  START_MASTER(threading)
  for(int i=1; i<threading->n_core; i++) {
    for(int j=0; j<count; j++) {
      ((complex_PRECISION **)threading->workspace)[0][j]        += ((complex_PRECISION **)threading->workspace)[i][j];
      ((complex_PRECISION **)threading->workspace)[0][j+offset] += ((complex_PRECISION **)threading->workspace)[i][j+offset];
    }
  }
  END_MASTER(threading)
  // only master needs the result in this case (it will be distributed later)
}
#endif


#if !defined( SSE ) || !defined( GRAM_SCHMIDT_VECTORIZED_PRECISION )
void setup_gram_schmidt_PRECISION_axpys(
    complex_PRECISION *thread_buffer, vector_PRECISION *V, int count, int offset,
    int start, int end, level_struct *l, struct Thread *threading) {
  
  int thread_start;
  int thread_end;
  int cache_block_size = 12*64;
  complex_PRECISION tmp[cache_block_size];

  compute_core_start_end_custom(start, end, &thread_start, &thread_end, l, threading, cache_block_size);

  for ( int i=thread_start; i<thread_end; i+=cache_block_size) {
    for ( int j=0; j<count; j++ ) {
      coarse_gamma5_PRECISION( tmp, V[j]+i, 0, cache_block_size, l );
      for ( int k=0; k<cache_block_size; k++) {
        V[count][i+k] -= thread_buffer[2*offset+j]*V[j][i+k];
        V[count][i+k] -= thread_buffer[3*offset+j]*tmp[k];
      }
    }
  }
}
#endif


void setup_gram_schmidt_PRECISION( vector_PRECISION *V, vector_PRECISION g5v,
                                   complex_PRECISION *buffer, const int n, level_struct *l,
                                   struct Thread *threading ) {
  
  PROF_PRECISION_START( _GRAM_SCHMIDT, threading );
  PRECISION beta;
  int i, j;
  int start = 0;
  int end = l->inner_vector_size;
  int thread_start = threading->start_index[l->depth];
  int thread_end   = threading->end_index[l->depth];

  complex_PRECISION thread_buffer[4*n];
  
  for ( i=0; i<4*n; i++ )
    thread_buffer[i] = 0;
  
  for ( i=0; i<n; i++ ) {
    
    if ( l->depth > 0 ) {
      coarse_gamma5_PRECISION( g5v, V[i], thread_start, thread_end, l );
      for ( j=0; j<i; j++ ) {
        thread_buffer[j] = process_inner_product_PRECISION( V[j], V[i], start, end, l, threading );
        thread_buffer[j+n] = process_inner_product_PRECISION( V[j], g5v, start, end, l, threading );
      }
    }
    else
      setup_gram_schmidt_PRECISION_compute_dots( thread_buffer, V, i, n, start, end, l, threading);
    
    
    START_LOCKED_MASTER(threading)
    if ( i>0 ) {
      PROF_PRECISION_START( _ALLR );
      MPI_Allreduce( thread_buffer, thread_buffer+2*n, 2*n, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
      PROF_PRECISION_STOP( _ALLR, 1 );
    }
    for ( j=2*n; j<4*n; j++ )
      ((complex_PRECISION *)(threading->workspace))[j] = thread_buffer[j];
    END_LOCKED_MASTER(threading)
    for ( j=2*n; j<4*n; j++ )
      thread_buffer[j] = ((complex_PRECISION *)(threading->workspace))[j];

    
    if ( l->depth > 0 ) {
      for( j=0; j<i; j++ ) {
        vector_PRECISION_saxpy( V[i], V[i], V[j], -(thread_buffer+2*n)[j], thread_start, thread_end, l );
        coarse_gamma5_PRECISION( g5v, V[j], thread_start, thread_end, l );
        vector_PRECISION_saxpy( V[i], V[i], g5v, -(thread_buffer+3*n)[j], thread_start, thread_end, l );
      }
    } else {
      setup_gram_schmidt_PRECISION_axpys( thread_buffer, V, i, n, start, end, l, threading);
    }
    
    beta = global_norm_PRECISION( V[i], start, end, l, threading );
    vector_PRECISION_real_scale( V[i], V[i], 1.0/beta, thread_start, thread_end, l );
  }
  PROF_PRECISION_STOP( _GRAM_SCHMIDT, (double)(end-start)/(double)l->inner_vector_size, threading );
}

