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

#ifndef FLOAT_INTRINSIC_SSE_H
#define FLOAT_INTRINSIC_SSE_H

#ifdef SSE
#include <xmmintrin.h>
#include <emmintrin.h>

// res = a*b + c
static inline __m128 sse_fmadd( __m128 a, __m128 b, __m128 c ) {
  __m128 res;
  res = _mm_mul_ps( a, b );
  res = _mm_add_ps( res, c );
  return res;
}

// res = -a*b + c
static inline __m128 sse_fnmadd( __m128 a, __m128 b, __m128 c ) {
  __m128 res;
  res = _mm_mul_ps( a, b );
  res = _mm_sub_ps( c, res );
  return res;
}

// res = a*b - c
static inline __m128 sse_fmsub( __m128 a, __m128 b, __m128 c ) {
  __m128 res;
  res = _mm_mul_ps( a, b );
  res = _mm_sub_ps( res, c );
  return res;
}

// res = -a*b - c
static inline __m128 sse_fnmsub( __m128 a, __m128 b, __m128 c ) {
  __m128 res; __m128 minus_a;
  minus_a = _mm_setzero_ps();
  minus_a = _mm_sub_ps( minus_a, a );
  res = _mm_mul_ps( minus_a, b );
  res = _mm_sub_ps( res, c );
  return res;
}

static inline void transpose_4_registers( __m128 *data)
{
   __m128 tmp[4];

   tmp[0] = _mm_unpacklo_ps( data[0], data[1] );
   tmp[1] = _mm_unpacklo_ps( data[2], data[3] );
   tmp[2] = _mm_unpackhi_ps( data[0], data[1] );
   tmp[3] = _mm_unpackhi_ps( data[2], data[3] );

   data[0] = _mm_movelh_ps( tmp[0], tmp[1] );
   data[1] = _mm_movehl_ps( tmp[1], tmp[0] );
   data[2] = _mm_movelh_ps( tmp[2], tmp[3] );
   data[3] = _mm_movehl_ps( tmp[3], tmp[2] );
}


static inline float sse_reduce_add_ps( __m128 data ) {
  float result;
  
  __m128 tmp;
  tmp = _mm_add_ps( data, _mm_movehl_ps( data, data ) );
  data = _mm_add_ss( tmp, _mm_shuffle_ps( tmp, tmp, 1 ) );
  _mm_store_ss( &result, data );
  
  return result;
}

#endif

#endif // FLOAT_INTRINSIC_SSE_H
