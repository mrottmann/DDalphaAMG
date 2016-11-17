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

complex_PRECISION global_inner_product_PRECISION( vector_PRECISION phi, vector_PRECISION psi, int start, int end,
                                                  level_struct *l, struct Thread *threading ) {
  
  PROF_PRECISION_START( _GIP, threading );
  complex_PRECISION local_alpha = 0, global_alpha = 0;

  int thread_start;
  int thread_end;
  
  SYNC_CORES(threading);
  
#ifndef OPTIMIZED_LINALG_PRECISION

  compute_core_start_end(start, end, &thread_start, &thread_end, l, threading);
  VECTOR_FOR( int i=thread_start, i<thread_end, local_alpha += conj_PRECISION(phi[i])*psi[i], i++, l );
  
#else
  compute_core_start_end_custom( start, end, &thread_start, &thread_end, l, threading, SIMD_LENGTH_PRECISION );
  mm_PRECISION alpha_re = mm_setzero_PRECISION();
  mm_PRECISION alpha_im = mm_setzero_PRECISION();
  
  for( int i=thread_start; i<thread_end; i+=SIMD_LENGTH_PRECISION ) {
    mm_PRECISION phi_re; mm_PRECISION phi_im;
    mm_PRECISION psi_re; mm_PRECISION psi_im;
    cload_PRECISION( (PRECISION*)(phi+i), &phi_re, &phi_im );
    cload_PRECISION( (PRECISION*)(psi+i), &psi_re, &psi_im );
    cfmadd_conj_PRECISION( phi_re, phi_im, psi_re, psi_im, &alpha_re, &alpha_im );
  }
  local_alpha = mm_reduce_add_PRECISION( alpha_re ) + I* mm_reduce_add_PRECISION( alpha_im );
#endif

  // sum over cores
  START_NO_HYPERTHREADS(threading)
    ((complex_PRECISION *)threading->workspace)[threading->core] = local_alpha;
  END_NO_HYPERTHREADS(threading);
  // master sums up all results
  SYNC_CORES(threading);
  MASTER(threading) {
    for(int i=1; i<threading->n_core; i++)
      ((complex_PRECISION *)threading->workspace)[0] += ((complex_PRECISION *)threading->workspace)[i];
    local_alpha = ((complex_PRECISION *)threading->workspace)[0];
  }
  
  if ( g.num_processes > 1 ) {
    MASTER(threading) {
      PROF_PRECISION_START( _ALLR );
      MPI_Allreduce( &local_alpha, &global_alpha, 1, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
      PROF_PRECISION_STOP( _ALLR, 1 );
      ((complex_PRECISION *)threading->workspace)[0] = global_alpha;
    }
    // all threads need the result of the norm
    SYNC_MASTER_TO_ALL(threading);
    global_alpha = ((complex_PRECISION *)threading->workspace)[0];
    PROF_PRECISION_STOP( _GIP, (double)(end-start)/(double)l->inner_vector_size, threading );
    return global_alpha;
  } else {
    // all threads need the result of the norm
    SYNC_MASTER_TO_ALL(threading);
    local_alpha = ((complex_PRECISION *)threading->workspace)[0];
    PROF_PRECISION_STOP( _GIP, (double)(end-start)/(double)l->inner_vector_size, threading );
    return local_alpha;
  }
}


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


void process_multi_inner_product_PRECISION( int count, complex_PRECISION *results, vector_PRECISION *phi, vector_PRECISION psi,
    int start, int end, level_struct *l, struct Thread *threading ) {

  PROF_PRECISION_START( _PIP, threading );
  int i;
  for(int c=0; c<count; c++)
    results[c] = 0.0;

  int thread_start;
  int thread_end;

  SYNC_CORES(threading)

#ifndef OPTIMIZED_LINALG_PRECISION  

  if ( l->depth == 0 ) {
    compute_core_start_end_custom(start, end, &thread_start, &thread_end, l, threading, 12);
    for(int c=0; c<count; c++)
      for ( i=thread_start; i<thread_end; )
        FOR12( results[c] += conj_PRECISION(phi[c][i])*psi[i]; i++; )
  } else {
    compute_core_start_end_custom(start, end, &thread_start, &thread_end, l, threading, 2);
    for(int c=0; c<count; c++)
      for ( i=thread_start; i<thread_end; )
        FOR2( results[c] += conj_PRECISION(phi[c][i])*psi[i]; i++; )
  }

#else
  compute_core_start_end_custom(start, end, &thread_start, &thread_end, l, threading, SIMD_LENGTH_PRECISION);
  for(int c=0; c<count; c++) {
    mm_PRECISION result_re = mm_setzero_PRECISION();
    mm_PRECISION result_im = mm_setzero_PRECISION();
    for ( i=thread_start; i<thread_end; i+=SIMD_LENGTH_PRECISION ) {
      mm_PRECISION phi_re; mm_PRECISION phi_im;
      mm_PRECISION pdi_re; mm_PRECISION pdi_im;
      
      // deinterleave complex numbers into 4 real parts and 4 imag parts        
      cload_PRECISION( (PRECISION*)(phi[c]+i), &phi_re, &phi_im );
      cload_PRECISION( (PRECISION*)(psi+i), &pdi_re, &pdi_im );
      
      cfmadd_conj_PRECISION(phi_re, phi_im, pdi_re, pdi_im, &result_re, &result_im);
    }
    results[c] += mm_reduce_add_PRECISION(result_re) + I*mm_reduce_add_PRECISION(result_im);
  }
#endif

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


complex_PRECISION local_xy_over_xx_PRECISION( vector_PRECISION phi, vector_PRECISION psi, int start, int end, level_struct *l  ) {
  
  complex_PRECISION numerator = 0.0; PRECISION denominator = 0.0;
  
  VECTOR_FOR( int i=start, i<end, numerator += conj_PRECISION(phi[i])*psi[i]; denominator += NORM_SQUARE_PRECISION(phi[i]), i++, l );
  
  if ( abs_PRECISION(denominator) < EPS_PRECISION ) {
    return 0.0;
  }
  
  return numerator/denominator;
}


PRECISION global_norm_PRECISION( vector_PRECISION x, int start, int end, level_struct *l, struct Thread *threading ) {
  
  PROF_PRECISION_START( _GIP, threading );
  
  PRECISION local_alpha = 0, global_alpha = 0;

  int thread_start;
  int thread_end;
  
  SYNC_CORES(threading);
  
#ifndef OPTIMIZED_LINALG_PRECISION
  
  compute_core_start_end(start, end, &thread_start, &thread_end, l, threading);
  VECTOR_FOR( int i=thread_start, i<thread_end, local_alpha += NORM_SQUARE_PRECISION(x[i]), i++, l );
  
#else
  compute_core_start_end_custom(start, end, &thread_start, &thread_end, l, threading, SIMD_LENGTH_PRECISION);
  mm_PRECISION alpha = mm_setzero_PRECISION(); 
  for( int i=thread_start; i<thread_end; i += SIMD_LENGTH_PRECISION/2 ) {
    mm_PRECISION phi = mm_load_PRECISION((PRECISION*)(x+i));
    alpha = mm_fmadd_PRECISION( phi, phi, alpha );
  }
  local_alpha = mm_reduce_add_PRECISION( alpha );
#endif

  // sum over cores
  START_NO_HYPERTHREADS(threading)
    ((PRECISION *)threading->workspace)[threading->core] = local_alpha;
  END_NO_HYPERTHREADS(threading);
  // master sums up all results
  SYNC_CORES(threading);
  MASTER(threading) {
    for(int i=1; i<threading->n_core; i++)
      ((PRECISION *)threading->workspace)[0] += ((PRECISION *)threading->workspace)[i];
    local_alpha = ((PRECISION *)threading->workspace)[0];
  }

  if ( g.num_processes > 1 ) {
    MASTER(threading) {
      PROF_PRECISION_START( _ALLR );
      MPI_Allreduce( &local_alpha, &global_alpha, 1, MPI_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
      PROF_PRECISION_STOP( _ALLR, 1 );
      ((PRECISION *)threading->workspace)[0] = global_alpha;
    }
    // all threads need the result of the norm
    SYNC_MASTER_TO_ALL(threading);
    global_alpha = ((PRECISION *)threading->workspace)[0];
    PROF_PRECISION_STOP( _GIP, (double)(end-start)/(double)l->inner_vector_size, threading );
    return (PRECISION)sqrt((double)global_alpha);
  } else {
    // all threads need the result of the norm
    SYNC_MASTER_TO_ALL(threading);
    local_alpha = ((PRECISION *)threading->workspace)[0];
    PROF_PRECISION_STOP( _GIP, (double)(end-start)/(double)l->inner_vector_size, threading );
    return (PRECISION)sqrt((double)local_alpha);
  }
}


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

// vector storage for PRECISION precision
void vector_PRECISION_define( vector_PRECISION phi, complex_PRECISION value, int start, int end, level_struct *l, Thread *threading ) {
  
  int i;
  PROF_PRECISION_START( _SET, threading );

  THREADED_VECTOR_FOR( i, start, end, phi[i] = value, i++, l, threading );

  PROF_PRECISION_STOP( _SET, 1, threading );
}

void vector_PRECISION_define_real( vector_PRECISION phi, PRECISION value, int start, int end, level_struct *l, Thread *threading ) {
  
  int i;
  PROF_PRECISION_START( _SET, threading );

  PRECISION *phi_pt = (PRECISION*) phi;
  THREADED_VECTOR_FOR( i, 2*start, 2*end, phi_pt[i] = value; phi_pt[i+1] = 0, i+=2, l, threading );

  PROF_PRECISION_STOP( _SET, 1, threading );
}

void vector_PRECISION_define_zero( vector_PRECISION phi, int start, int end, level_struct *l, Thread *threading ) {
  
  int i;
  PROF_PRECISION_START( _SET, threading );

  PRECISION *phi_pt = (PRECISION*) phi;
  THREADED_VECTOR_FOR( i, 2*start, 2*end, phi_pt[i] = phi_pt[i+1] = 0, i+=2, l, threading );

  PROF_PRECISION_STOP( _SET, 1, threading );
}


void vector_PRECISION_define_random( vector_PRECISION phi, int start, int end, level_struct *l, Thread *threading ) {
  
  int i;
  PROF_PRECISION_START( _SET, threading );
  
  // this would yield different results if we threaded it, so we don't
  START_LOCKED_MASTER(threading)
  VECTOR_FOR( i=start, i<end, phi[i] = (PRECISION)(((double)rand()/(double)RAND_MAX))-0.5 + ( (PRECISION)((double)rand()/(double)RAND_MAX)-0.5)*_Complex_I, i++, l );
  END_LOCKED_MASTER(threading)

  PROF_PRECISION_STOP( _SET, 1, threading );
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

void vector_PRECISION_scale( vector_PRECISION z, vector_PRECISION x, complex_PRECISION alpha, int start, int end, level_struct *l, struct Thread *threading ) {
  
  int thread_start, thread_end;
  PROF_PRECISION_START( _LA6, threading );
  
#ifndef OPTIMIZED_LINALG_PRECISION
  
  compute_core_start_end(start, end, &thread_start, &thread_end, l, threading);
  VECTOR_FOR( int i=start, i<end, z[i] = alpha*x[i], i++, l );
  
#else
  compute_core_start_end_custom(start, end, &thread_start, &thread_end, l, threading, SIMD_LENGTH_PRECISION);
  mm_PRECISION alpha_re = mm_set1_PRECISION( creal_PRECISION(alpha) );
  mm_PRECISION alpha_im = mm_set1_PRECISION( cimag_PRECISION(alpha) );
  
  for( int i=start; i<end; i+=SIMD_LENGTH_PRECISION ) {
    mm_PRECISION z_re, z_im, x_re, x_im;
    cload_PRECISION( (PRECISION*)(x+i), &x_re, &x_im );
    cmul_PRECISION( alpha_re, alpha_im, x_re, x_im, &z_re, &z_im );
    cstore_PRECISION( (PRECISION*)(z+i), z_re, z_im );
  }
#endif
  
  PROF_PRECISION_STOP( _LA6, (double)(end-start)/(double)l->inner_vector_size, threading );
}


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

void vector_PRECISION_saxpy( vector_PRECISION z, vector_PRECISION x, vector_PRECISION y, complex_PRECISION alpha,
                             int start, int end, level_struct *l, struct Thread *threading ) {
  
 int thread_start, thread_end;
 PROF_PRECISION_START( _LA8, threading );
  
#ifndef OPTIMIZED_LINALG_PRECISION
  
  compute_core_start_end(start, end, &thread_start, &thread_end, l, threading);
  VECTOR_FOR( int i=start, i<end, z[i] = x[i] + alpha*y[i], i++, l );
  
#else
  compute_core_start_end_custom(start, end, &thread_start, &thread_end, l, threading, SIMD_LENGTH_PRECISION);
  mm_PRECISION alpha_re = mm_set1_PRECISION( creal_PRECISION(alpha) );
  mm_PRECISION alpha_im = mm_set1_PRECISION( cimag_PRECISION(alpha) );
  
  for ( int i=start; i<end; i+=SIMD_LENGTH_PRECISION ) {
    mm_PRECISION x_re, x_im, y_re, y_im;
    cload_PRECISION( (PRECISION*)(x+i), &x_re, &x_im );
    cload_PRECISION( (PRECISION*)(y+i), &y_re, &y_im );
    cfmadd_PRECISION(alpha_re, alpha_im, y_re, y_im, &x_re, &x_im);
    cstore_PRECISION( (PRECISION*)(z+i), x_re, x_im );
  }
#endif

  PROF_PRECISION_STOP( _LA8, (double)(end-start)/(double)l->inner_vector_size, threading );
}

void vector_PRECISION_multi_saxpy( vector_PRECISION z, vector_PRECISION *V, complex_PRECISION *alpha, int sign, 
                                   int count, int start, int end, level_struct *l, struct Thread *threading ) {
  
 int thread_start, thread_end;
 PROF_PRECISION_START( _LA8, threading );
  
#ifndef OPTIMIZED_LINALG_PRECISION

  complex_PRECISION alpha_signed[count];
  for ( int c=0; c<count; c++ ) {
    alpha_signed[c] = sign*alpha[c];
  }
  
  compute_core_start_end_custom(start, end, &thread_start, &thread_end, l, threading, 12);
  for ( int c=0; c<count; c++ ) {
    for ( int i=thread_start; i<thread_end; ) {
      FOR12( z[i] += V[c][i]*alpha_signed[c]; i++; )
    }
  }  

#else

  compute_core_start_end_custom(start, end, &thread_start, &thread_end, l, threading, SIMD_LENGTH_PRECISION);
  for ( int c=0; c<count; c++ ) {
    mm_PRECISION alpha_re, alpha_im;
    alpha_re = mm_set1_PRECISION( sign*creal_PRECISION(alpha[c]) );
    alpha_im = mm_set1_PRECISION( sign*cimag_PRECISION(alpha[c]) );

    if ( abs_PRECISION(cimag_PRECISION(alpha[c])) < EPS_PRECISION ) {
      for ( int i=thread_start; i<thread_end; i+=SIMD_LENGTH_PRECISION/2 ) {
        mm_PRECISION z_re = mm_load_PRECISION( (PRECISION*)(z+i) );
        mm_PRECISION V_re = mm_load_PRECISION( (PRECISION*)(V[c]+i) );
        z_re = mm_fmadd_PRECISION( alpha_re, V_re, z_re );
        mm_store_PRECISION( (PRECISION*)(z+i), z_re );
      }
    } else {
      for ( int i=thread_start; i<thread_end; i+=SIMD_LENGTH_PRECISION ) {
        mm_PRECISION z_re, z_im, V_re, V_im; 
        cload_PRECISION( (PRECISION*)(z+i), &z_re, &z_im );
        cload_PRECISION( (PRECISION*)(V[c]+i), &V_re, &V_im );
        cfmadd_PRECISION( alpha_re, alpha_im, V_re, V_im, &z_re, &z_im );
        cstore_PRECISION( (PRECISION*)(z+i), z_re, z_im );
      }
    }
  }
#endif
  
  PROF_PRECISION_STOP( _LA8, (PRECISION)(count), threading );
}


void vector_PRECISION_projection( vector_PRECISION z, vector_PRECISION v, int k, vector_PRECISION *W, complex_PRECISION *diag, 
                                  int orthogonal, level_struct *l, Thread *threading ) {
  
  int j, start, end;
  
  compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
                    
  vector_PRECISION v_tmp = NULL, *W_tmp = NULL;
  complex_PRECISION ip[k], ip_buffer[2*k];      
  
  MALLOC( v_tmp, complex_PRECISION, l->inner_vector_size );
  vector_PRECISION_define_zero( v_tmp,  0, l->inner_vector_size, l, threading );
  
  MALLOC( W_tmp, complex_PRECISION*, k );
  W_tmp[0] = NULL; 
  MALLOC( W_tmp[0], complex_PRECISION, k*l->inner_vector_size );
  for ( j = 1; j<k; j++ )
    W_tmp[j] = W_tmp[0]+j*l->inner_vector_size;
  
  for ( j=0; j<k; j++ ) {
    vector_PRECISION_scale( W_tmp[j], W[j], diag[j], 0, l->inner_vector_size, l, threading );
  }
  process_multi_inner_product_PRECISION( k, ip, W_tmp, v, 0, l->inner_vector_size, l, threading );
  
  START_MASTER(threading)
  for ( j=0; j<k; j++ ) {
    ip_buffer[j] = ip[j];
  }
  MPI_Allreduce( ip_buffer, ip_buffer+k, k, MPI_COMPLEX_PRECISION, MPI_SUM, (l->depth==0)?g.comm_cart:l->gs_PRECISION.level_comm );
  END_MASTER(threading);
  SYNC_MASTER_TO_ALL(threading);
  
  vector_PRECISION_multi_saxpy( v_tmp, W_tmp, ip_buffer+k, 1, k, 0, l->inner_vector_size, l, threading );
   
  if (orthogonal) 
    vector_PRECISION_minus( z, v, v_tmp, 0, l->inner_vector_size, l );
  else
    vector_PRECISION_copy( z, v_tmp, 0, l->inner_vector_size, l );
  
  FREE( v_tmp, complex_PRECISION, l->inner_vector_size );
  FREE( W_tmp[0], complex_PRECISION, k*l->inner_vector_size );
  FREE( W_tmp, complex_PRECISION*, k );
}

void set_boundary_PRECISION( vector_PRECISION phi, complex_PRECISION alpha, level_struct *l, struct Thread *threading ) {
  
  PROF_PRECISION_START( _SET, threading );
  int i;
  
  SYNC_CORES(threading)
  
  THREADED_VECTOR_FOR( i, l->inner_vector_size, l->vector_size, phi[i] = alpha, i++, l, threading );
  
  SYNC_CORES(threading)
  PROF_PRECISION_STOP( _SET, (double)(l->vector_size-l->inner_vector_size)/(double)l->inner_vector_size, threading );
}
