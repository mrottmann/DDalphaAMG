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

#ifndef LINALG_PRECISION_HEADER
  #define LINALG_PRECISION_HEADER
  
#ifdef _M10TV
  #define VECTOR_FOR( start, end, expression, update, l ) do{ \
    if ( l->depth == 0 ) { \
      for ( start; end; ) \
        FOR12( expression; update; ) \
    } else { \
      for ( start; end; ) \
        FOR20( expression; update; ) \
    } \
  } while(0)
#else
  #define VECTOR_FOR( start, end, expression, update, l ) do{ \
    if ( l->depth == 0 ) { \
      for ( start; end; ) \
        FOR12( expression; update; ) \
    } else { \
      for ( start; end; ) \
        FOR2( expression; update; ) \
    } \
  } while(0)
#endif


#ifdef _M10TV
  #define REAL_VECTOR_FOR( start, end, expression, update, l ) do{ \
    if ( l->depth == 0 ) { \
      for ( start; end; ) \
        FOR24( expression; update; ) \
    } else { \
      for ( start; end; ) \
        FOR40( expression; update; ) \
    } \
  } while(0)
#else
  #define REAL_VECTOR_FOR( start, end, expression, update, l ) do{ \
    if ( l->depth == 0 ) { \
      for ( start; end; ) \
        FOR24( expression; update; ) \
    } else { \
      for ( start; end; ) \
        FOR4( expression; update; ) \
    } \
  } while(0)
#endif


#ifdef _M10TV
  #define THREADED_VECTOR_FOR( i, start_index, end_index, expression, update, l, threading ) do{ \
    int thread_start, thread_end; \
    if ( l->depth == 0 ) { \
      compute_core_start_end_custom(start_index, end_index, &thread_start, &thread_end, l, threading, 12); \
      for ( i=thread_start; i<thread_end; ) \
        FOR12( expression; update; ) \
    } else { \
      compute_core_start_end_custom(start_index, end_index, &thread_start, &thread_end, l, threading, 20); \
      for ( i=thread_start; i<thread_end; ) \
        FOR20( expression; update; ) \
    } \
  } while(0)
#else
  #define THREADED_VECTOR_FOR( i, start_index, end_index, expression, update, l, threading ) do{ \
    int thread_start, thread_end; \
    if ( l->depth == 0 ) { \
      compute_core_start_end_custom(start_index, end_index, &thread_start, &thread_end, l, threading, 12); \
      for ( i=thread_start; i<thread_end; ) \
        FOR12( expression; update; ) \
    } else { \
      compute_core_start_end_custom(start_index, end_index, &thread_start, &thread_end, l, threading, 2); \
      for ( i=thread_start; i<thread_end; ) \
        FOR2( expression; update; ) \
    } \
  } while(0)
#endif
  

  struct Thread;

  complex_PRECISION global_inner_product_PRECISION( vector_PRECISION x, vector_PRECISION y, int start, int end, level_struct *l, struct Thread *threading );
  complex_PRECISION process_inner_product_PRECISION( vector_PRECISION phi, vector_PRECISION psi, int start, int end, level_struct *l, struct Thread *threading );

  void process_multi_inner_product_PRECISION( int count, complex_PRECISION *results, vector_PRECISION *phi, vector_PRECISION psi,
      int start, int end, level_struct *l, struct Thread *threading );

  PRECISION global_norm_PRECISION( vector_PRECISION phi, int start, int end, level_struct *l, struct Thread *threading );
  PRECISION process_norm_PRECISION( vector_PRECISION x, int start, int end, level_struct *l, struct Thread *threading );
  
  complex_PRECISION local_xy_over_xx_PRECISION( vector_PRECISION phi, vector_PRECISION psi, int start, int end, level_struct *l  );
  void vector_PRECISION_plus( vector_PRECISION z, vector_PRECISION x, vector_PRECISION y, int start, int end, level_struct *l ); // z := x + y
  void vector_PRECISION_minus( vector_PRECISION z, vector_PRECISION x, vector_PRECISION y, int start, int end, level_struct *l ); // z := x - y
  void vector_PRECISION_scale( vector_PRECISION z, vector_PRECISION x, complex_PRECISION alpha, int start, int end, level_struct *l ); // z := alpha*x
  void vector_PRECISION_real_scale( vector_PRECISION z, vector_PRECISION x, complex_PRECISION alpha,
                                    int start, int end, level_struct *l );
  void vector_PRECISION_saxpy( vector_PRECISION z, vector_PRECISION x, vector_PRECISION y, complex_PRECISION alpha, int start, int end, level_struct *l ); // z := x + alpha*y
  void vector_PRECISION_copy( vector_PRECISION z, vector_PRECISION x, int start, int end, level_struct *l ); // z := x
  void vector_PRECISION_projection( vector_PRECISION z, vector_PRECISION v, int k, vector_PRECISION *W, complex_PRECISION *diag, 
                                  int orthogonal, level_struct *l, Thread *threading );
  
  void gram_schmidt_on_aggregates_PRECISION( vector_PRECISION *V, const int num_vec, level_struct *l, struct Thread *threading );
  
  // Gram-Schmidt on a block of vectors, used by Block-Gram-Schmidt
  void aggregate_gram_schmidt_block_PRECISION( PRECISION *V,
      int num_vec, int leading_dimension, level_struct *l, struct Thread *threading );
  // used by Block-Gram-Schmidt
  void aggregate_orthogonalize_block_wrt_orthonormal_block_PRECISION( PRECISION *B, PRECISION *U,
      int num_vec, level_struct *l, struct Thread *threading );
  // used by Block-Gram-Schmidt
  void aggregate_block_dot_block_PRECISION( PRECISION *S, PRECISION *U, PRECISION *B,
      int num_vec, int leading_dimension, level_struct *l, struct Thread *threading );
  // used by Block-Gram-Schmidt
  void aggregate_block_minus_block_times_dot_PRECISION( PRECISION *B, PRECISION *U, PRECISION *S,
      int num_vec, int leading_dimension, level_struct *l, struct Thread *threading );

  void gram_schmidt_PRECISION( vector_PRECISION *V, complex_PRECISION *buffer, const int start, const int n, level_struct *l, struct Thread *threading );
  void setup_gram_schmidt_PRECISION( vector_PRECISION *V, vector_PRECISION g5v,
                                     complex_PRECISION *buffer, const int n, level_struct *l,
                                     struct Thread *threading );
  void spinwise_PRECISION_skalarmultiply( vector_PRECISION eta1, vector_PRECISION eta2,
                                          vector_PRECISION phi, complex_PRECISION alpha, int start, int end, level_struct *l );
  void set_boundary_PRECISION( vector_PRECISION phi, complex_PRECISION alpha, level_struct *l, struct Thread *threading );
  
#endif
