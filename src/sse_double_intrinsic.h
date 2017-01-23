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

#ifndef DOUBLE_INTRINSIC_SSE_H
#define DOUBLE_INTRINSIC_SSE_H

#ifdef SSE
#include <xmmintrin.h>
#include <emmintrin.h>

// res = a*b + c
static inline __m128d sse_fmadd_pd( __m128d a, __m128d b, __m128d c ) {
  __m128d res;
  res = _mm_mul_pd( a, b );
  res = _mm_add_pd( res, c );
  return res;
}

// res = -a*b + c
static inline __m128d sse_fnmadd_pd( __m128d a, __m128d b, __m128d c ) {
  __m128d res;
  res = _mm_mul_pd( a, b );
  res = _mm_sub_pd( c, res );
  return res;
}

// res = a*b - c
static inline __m128d sse_fmsub_pd( __m128d a, __m128d b, __m128d c ) {
  __m128d res;
  res = _mm_mul_pd( a, b );
  res = _mm_sub_pd( res, c );
  return res;
}

static inline double sse_reduce_add_pd( __m128d data ) {
  double result;
  data = _mm_add_pd( data, _mm_unpackhi_pd( data, data ) );
  _mm_store_sd( &result, data );
  return result;
}

#endif
#endif