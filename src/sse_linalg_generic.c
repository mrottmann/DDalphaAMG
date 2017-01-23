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

#ifdef SSE

#include "sse_complex_float_intrinsic.h"
#include "sse_float_intrinsic.h"
#include "sse_linalg.h"

void aggregate_gram_schmidt_PRECISION_vectorized( complex_PRECISION *V, const int num_vec, level_struct *l, struct Thread *threading ) {
    sse_aggregate_gram_schmidt_PRECISION( V, num_vec, l, threading );
}


void aggregate_gram_schmidt_block_PRECISION( PRECISION *V,
    int num_vec, int leading_dimension, level_struct *l, struct Thread *threading ) {
    sse_aggregate_gram_schmidt_block_PRECISION( V, num_vec, leading_dimension, l, threading );
}


void aggregate_block_gram_schmidt_PRECISION_vectorized( complex_PRECISION *V, const int num_vec, level_struct *l, struct Thread *threading ) {

  PROF_PRECISION_START( _GRAM_SCHMIDT_ON_AGGREGATES, threading );
  SYNC_CORES(threading)
  for ( int i=0; i<num_vec; i+=SIMD_LENGTH_PRECISION ) {

    int vecs = SIMD_LENGTH_PRECISION;
    if(num_vec-i < SIMD_LENGTH_PRECISION)
      vecs = num_vec-i;

    for ( int j=0; j<i; j+=SIMD_LENGTH_PRECISION )
      aggregate_orthogonalize_block_wrt_orthonormal_block_PRECISION( (PRECISION *)(V + i*l->vector_size),
          (PRECISION *)(V + j*l->vector_size), vecs, l, threading );
    aggregate_gram_schmidt_block_PRECISION( (PRECISION *)(V + i*l->vector_size), vecs, SIMD_LENGTH_PRECISION, l, threading );
  }
  SYNC_CORES(threading)
  PROF_PRECISION_STOP( _GRAM_SCHMIDT_ON_AGGREGATES, 1, threading );
}


void gram_schmidt_on_aggregates_PRECISION_vectorized( complex_PRECISION *V, const int num_vec, level_struct *l, struct Thread *threading ) {

  // the block version has some optimizations which are correct only on the fine grid
  if(l->depth == 0)
    aggregate_block_gram_schmidt_PRECISION_vectorized(V, num_vec, l, threading);
  else
    aggregate_gram_schmidt_PRECISION_vectorized(V, num_vec, l, threading);
}


void aggregate_orthogonalize_block_wrt_orthonormal_block_PRECISION( PRECISION *B, PRECISION *U, int num_vec, level_struct *l, struct Thread *threading ) {
  START_NO_HYPERTHREADS(threading)

  PRECISION *S = NULL;
  START_LOCKED_MASTER(threading)
  // factors 2 are for complex and spin01/23 aggregates
  MALLOC_HUGEPAGES(S, PRECISION, 2*2*l->s_PRECISION.num_aggregates*SIMD_LENGTH_PRECISION*SIMD_LENGTH_PRECISION, 64);
  ((PRECISION **)threading->workspace)[0] = S;
  END_LOCKED_MASTER(threading)
  S = ((PRECISION **)threading->workspace)[0];

  aggregate_block_dot_block_PRECISION(S, U, B, num_vec, SIMD_LENGTH_PRECISION, l , threading);
  aggregate_block_minus_block_times_dot_PRECISION(B, U, S, num_vec, SIMD_LENGTH_PRECISION, l , threading);

  START_LOCKED_MASTER(threading)
  FREE_HUGEPAGES(S, PRECISION, 2*2*l->s_PRECISION.num_aggregates*SIMD_LENGTH_PRECISION*SIMD_LENGTH_PRECISION);
  END_LOCKED_MASTER(threading)

  END_NO_HYPERTHREADS(threading)
}


void aggregate_block_dot_block_PRECISION( PRECISION *S, PRECISION *U, PRECISION *B,
    int num_vec, int leading_dimension, level_struct *l, struct Thread *threading ) {
    sse_aggregate_block_dot_block_PRECISION( S, U, B, num_vec, leading_dimension, l, threading );
}


void aggregate_block_minus_block_times_dot_PRECISION( PRECISION *B, PRECISION *U, PRECISION *S,
    int num_vec, int leading_dimension, level_struct *l, struct Thread *threading ) {
    sse_aggregate_block_minus_block_times_dot_PRECISION( B, U, S, num_vec, leading_dimension, l, threading );
}

#ifdef GRAM_SCHMIDT_VECTORIZED_PRECISION
void setup_gram_schmidt_PRECISION_compute_dots(
    complex_PRECISION *thread_buffer, vector_PRECISION *V, int count, int offset,
    int start, int end, level_struct *l, struct Thread *threading) {

  int thread_start;
  int thread_end;
  int cache_block_size = 12*16;

  for(int i=0; i<2*offset; i++)
    thread_buffer[i] = 0.0;

  SYNC_CORES(threading)
  compute_core_start_end_custom(start, end, &thread_start, &thread_end, l, threading, cache_block_size);

  __m128 dot_re[count];
  __m128 dot_im[count];
  __m128 dot_gamma5_re[count];
  __m128 dot_gamma5_im[count];

  for ( int j=0; j<count; j++) {
    dot_re[j] = _mm_setzero_ps();
    dot_im[j] = _mm_setzero_ps();
    dot_gamma5_re[j] = _mm_setzero_ps();
    dot_gamma5_im[j] = _mm_setzero_ps();
  }

  for ( int i=thread_start; i<thread_end; i+=cache_block_size) {
    for ( int j=0; j<count; j++ ) {
      for ( int k=0; k<cache_block_size; k+=12) {
        __m128 vj_re;
        __m128 vj_im;
        __m128 v_re;
        __m128 v_im;
        __m128 gamma5_v_re;
        __m128 gamma5_v_im;

        // gamma5 multiplies the first 6 out of 12 components with -1
        // SIMD_LENGTH is 4, so the pattern repeats after 12 elements = 3 cachelines
        // => can use 3 pre-defined +/-1 patterns
        __m128 gamma5[3];
        gamma5[0] = _mm_set_ps( -1.0,-1.0,-1.0,-1.0 );
        gamma5[1] = _mm_set_ps(  1.0, 1.0,-1.0,-1.0 );
        gamma5[2] = _mm_set_ps(  1.0, 1.0, 1.0, 1.0 );

        for(int m=0; m<3; m++) {
          
          sse_complex_deinterleaved_load( (float*)(V[j]+i+k+4*m), &vj_re, &vj_im  );
          sse_complex_deinterleaved_load( (float*)(V[count]+i+k+4*m), &v_re, &v_im  );

          gamma5_v_re = _mm_mul_ps(gamma5[m], v_re);
          gamma5_v_im = _mm_mul_ps(gamma5[m], v_im);

          cfmadd_conj(vj_re, vj_im, v_re, v_im, dot_re+j, dot_im+j);
          cfmadd_conj(vj_re, vj_im, gamma5_v_re, gamma5_v_im, dot_gamma5_re+j, dot_gamma5_im+j);
        }
      }
    }
  }
  for ( int j=0; j<count; j++ ) {
    thread_buffer[j]        = sse_reduce_add_ps(dot_re[j]) + I * sse_reduce_add_ps(dot_im[j]);
    thread_buffer[j+offset] = sse_reduce_add_ps(dot_gamma5_re[j]) + I * sse_reduce_add_ps(dot_gamma5_im[j]);
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

#ifdef GRAM_SCHMIDT_VECTORIZED_PRECISION
void setup_gram_schmidt_PRECISION_axpys(
    complex_PRECISION *thread_buffer, vector_PRECISION *V, int count, int offset,
    int start, int end, level_struct *l, struct Thread *threading) {

  int thread_start;
  int thread_end;
  int cache_block_size = 12*16;

  compute_core_start_end_custom(start, end, &thread_start, &thread_end, l, threading, cache_block_size);

  __m128 dot_re[count];
  __m128 dot_im[count];
  __m128 dot_gamma5_re[count];
  __m128 dot_gamma5_im[count];

  for ( int j=0; j<count; j++) {
    dot_re[j] = _mm_set1_ps(creal(thread_buffer[2*offset+j]));
    dot_im[j] = _mm_set1_ps(cimag(thread_buffer[2*offset+j]));
    dot_gamma5_re[j] = _mm_set1_ps(creal(thread_buffer[3*offset+j]));
    dot_gamma5_im[j] = _mm_set1_ps(cimag(thread_buffer[3*offset+j]));
  }

  for ( int i=thread_start; i<thread_end; i+=cache_block_size) {
    for ( int j=0; j<count; j++ ) {
      for ( int k=0; k<cache_block_size; k+=12) {
        __m128 vj_re;
        __m128 vj_im;
        __m128 gamma5_vj_re;
        __m128 gamma5_vj_im;
        __m128 v_re;
        __m128 v_im;

        // gamma5 multiplies the first 6 out of 12 components with -1
        // SIMD_LENGTH is 4, so the pattern repeats after 12 elements = 3 cachelines
        // => can use 3 pre-defined +/-1 patterns
        __m128 gamma5[3];
        gamma5[0] = _mm_set_ps( -1.0,-1.0,-1.0,-1.0 );
        gamma5[1] = _mm_set_ps(  1.0, 1.0,-1.0,-1.0 );
        gamma5[2] = _mm_set_ps(  1.0, 1.0, 1.0, 1.0 );

        for(int m=0; m<3; m++) {
          
          sse_complex_deinterleaved_load( (float*)(V[j]+i+k+4*m), &vj_re, &vj_im  );
          sse_complex_deinterleaved_load( (float*)(V[count]+i+k+4*m), &v_re, &v_im  );

          gamma5_vj_re = _mm_mul_ps(gamma5[m], vj_re);
          gamma5_vj_im = _mm_mul_ps(gamma5[m], vj_im);

          cfnmadd(vj_re, vj_im, dot_re[j], dot_im[j], &v_re, &v_im);
          cfnmadd(gamma5_vj_re, gamma5_vj_im, dot_gamma5_re[j], dot_gamma5_im[j], &v_re, &v_im);

          sse_complex_interleaved_store(v_re, v_im, (float*)(V[count]+i+k+4*m) ); 
        }
      }
    }
  }
}
#endif

#endif
