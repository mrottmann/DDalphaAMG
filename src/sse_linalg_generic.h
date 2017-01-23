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

#ifndef SSE_LINALG_PRECISION_HEADER
  #define SSE_LINALG_PRECISION_HEADER
  #ifdef SSE
  
  void gram_schmidt_on_aggregates_PRECISION_vectorized( complex_PRECISION *V, const int num_vec, level_struct *l, struct Thread *threading );
  // Block-Gram-Schmidt on aggregates
  void aggregate_block_gram_schmidt_PRECISION_vectorized( complex_PRECISION *V, const int num_vec, level_struct *l, struct Thread *threading );
  // Standard Gram-Schmidt on aggregates
  void aggregate_gram_schmidt_PRECISION_vectorized( complex_PRECISION *V, const int num_vec, level_struct *l, struct Thread *threading );
  
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
  
  void setup_gram_schmidt_PRECISION_compute_dots(
    complex_PRECISION *thread_buffer, vector_PRECISION *V, int count, int offset,
    int start, int end, level_struct *l, struct Thread *threading);
  
  void setup_gram_schmidt_PRECISION_axpys(
    complex_PRECISION *thread_buffer, vector_PRECISION *V, int count, int offset,
    int start, int end, level_struct *l, struct Thread *threading);
  
#endif
#endif