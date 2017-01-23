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

#ifndef COMPLEX_SSE_INTRINSIC_H
#define COMPLEX_SSE_INTRINSIC_H

#ifdef SSE
#include "sse_float_intrinsic.h"

// c = a*b
static inline void cmul(
        __m128 a_real, __m128 a_imag,
        __m128 b_real, __m128 b_imag,
        __m128 *c_real, __m128 *c_imag)
{
    *c_real = _mm_mul_ps(a_imag, b_imag);
    *c_imag = _mm_mul_ps(a_imag, b_real);
    *c_real = sse_fmsub(a_real, b_real, *c_real);
    *c_imag = sse_fmadd(a_real, b_imag, *c_imag);
}

// c = -a*b
static inline void cnmul(
        __m128 a_real, __m128 a_imag,
        __m128 b_real, __m128 b_imag,
        __m128 *c_real, __m128 *c_imag)
{
    *c_real = _mm_mul_ps(a_imag, b_imag);
    *c_imag = _mm_mul_ps(a_imag, b_real);
    *c_real = sse_fnmadd(a_real, b_real, *c_real);
    *c_imag = sse_fnmsub(a_real, b_imag, *c_imag);
}

// c = conj(a)*b
static inline void cmul_conj(
        __m128 a_real, __m128 a_imag,
        __m128 b_real, __m128 b_imag,
        __m128 *c_real, __m128 *c_imag)
{
    *c_real = _mm_mul_ps(a_imag, b_imag);
    *c_imag = _mm_mul_ps(a_imag, b_real);
    *c_real = sse_fmadd(a_real, b_real, *c_real);
    *c_imag = sse_fmsub(a_real, b_imag, *c_imag);
}

// c = -conj(a)*b
static inline void cnmul_conj(
        __m128 a_real, __m128 a_imag,
        __m128 b_real, __m128 b_imag,
        __m128 *c_real, __m128 *c_imag)
{
    *c_real = _mm_mul_ps(a_imag, b_imag);
    *c_imag = _mm_mul_ps(a_imag, b_real);
    *c_real = sse_fnmsub(a_real, b_real, *c_real);
    *c_imag = sse_fnmadd(a_real, b_imag, *c_imag);
}

// c = a*b + c
static inline void cfmadd(
        __m128 a_real, __m128 a_imag,
        __m128 b_real, __m128 b_imag,
        __m128 *c_real, __m128 *c_imag)
{
  *c_real = sse_fmsub( a_imag, b_imag, *c_real );
  *c_imag = sse_fmadd( a_imag, b_real, *c_imag );
  *c_real = sse_fmsub( a_real, b_real, *c_real );
  *c_imag = sse_fmadd( a_real, b_imag, *c_imag );
}

// c = -a*b + c
static inline void masked_cfnmadd(
        __m128 a_real, __m128 a_imag,
        __m128 b_real, __m128 b_imag,
        __m128 *c_real, __m128 *c_imag, __m128 mask )
{
    __m128 minus_a_real; __m128 minus_a_imag;
    minus_a_real = _mm_setzero_ps();
    minus_a_imag = _mm_sub_ps( minus_a_real, a_imag );
    minus_a_real = _mm_sub_ps( minus_a_real, a_real );
    minus_a_real = _mm_and_ps( minus_a_real, mask );
    minus_a_imag = _mm_and_ps( minus_a_imag, mask );
    
    *c_real = sse_fmsub(minus_a_imag, b_imag, *c_real);
    *c_imag = sse_fmadd(minus_a_imag, b_real, *c_imag);
    *c_real = sse_fmsub(minus_a_real, b_real, *c_real);
    *c_imag = sse_fmadd(minus_a_real, b_imag, *c_imag);
}

// c = -a*b + c
static inline void cfnmadd(
        __m128 a_real, __m128 a_imag,
        __m128 b_real, __m128 b_imag,
        __m128 *c_real, __m128 *c_imag)
{
    __m128 minus_a_real; __m128 minus_a_imag;
    minus_a_real = _mm_setzero_ps();
    minus_a_imag = _mm_sub_ps( minus_a_real, a_imag );
    minus_a_real = _mm_sub_ps( minus_a_real, a_real );
    
    *c_real = sse_fmsub(minus_a_imag, b_imag, *c_real);
    *c_imag = sse_fmadd(minus_a_imag, b_real, *c_imag);
    *c_real = sse_fmsub(minus_a_real, b_real, *c_real);
    *c_imag = sse_fmadd(minus_a_real, b_imag, *c_imag);
}

// c = conj(a)*b + c
static inline void cfmadd_conj(
        __m128 a_real, __m128 a_imag,
        __m128 b_real, __m128 b_imag,
        __m128 *c_real, __m128 *c_imag)
{
    *c_real = sse_fmadd(a_imag, b_imag, *c_real);
    *c_imag = sse_fmsub(a_imag, b_real, *c_imag);
    *c_real = sse_fmadd(a_real, b_real, *c_real);
    *c_imag = sse_fmsub(a_real, b_imag, *c_imag);
}

// c = -conj(a)*b + c
static inline void cfnmadd_conj(
        __m128 a_real, __m128 a_imag,
        __m128 b_real, __m128 b_imag,
        __m128 *c_real, __m128 *c_imag)
{
    __m128 minus_a_real; __m128 minus_a_imag;
    minus_a_real = _mm_setzero_ps();
    minus_a_imag = _mm_sub_ps( minus_a_real, a_imag );
    minus_a_real = _mm_sub_ps( minus_a_real, a_real );
  
    *c_real = sse_fmadd(minus_a_imag, b_imag, *c_real);
    *c_imag = sse_fmsub(minus_a_imag, b_real, *c_imag);
    *c_real = sse_fmadd(minus_a_real, b_real, *c_real);
    *c_imag = sse_fmsub(minus_a_real, b_imag, *c_imag);
}

static inline void sse_complex_deinterleaved_load( float *data, __m128 *result_re, __m128 *result_im  ) {  
  *result_re = _mm_setr_ps( data[0], data[2], data[4], data[6] );
  *result_im = _mm_setr_ps( data[1], data[3], data[5], data[7] );
}


static inline void sse_complex_interleaved_store( __m128 data_re, __m128 data_im, float *result  ) { 
  _mm_storeu_ps( result,                   _mm_unpacklo_ps( data_re, data_im ) );
  _mm_storeu_ps( result+SIMD_LENGTH_float, _mm_unpackhi_ps( data_re, data_im ) );
}

#endif
#endif
