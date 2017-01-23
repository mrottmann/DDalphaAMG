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

#ifndef COMPLEX_SSE_INTRINSIC_DOUBLE_H
#define COMPLEX_SSE_INTRINSIC_DOUBLE_H

#ifdef SSE
#include "sse_double_intrinsic.h"

// c = a*b
static inline void cmul_pd(
        __m128d a_real, __m128d a_imag,
        __m128d b_real, __m128d b_imag,
        __m128d *c_real, __m128d *c_imag)
{
    *c_real = _mm_mul_pd(a_imag, b_imag);
    *c_imag = _mm_mul_pd(a_imag, b_real);
    *c_real = sse_fmsub_pd(a_real, b_real, *c_real);
    *c_imag = sse_fmadd_pd(a_real, b_imag, *c_imag);
}

// c = conj(a)*b
static inline void cmul_conj_pd(
        __m128d a_real, __m128d a_imag,
        __m128d b_real, __m128d b_imag,
        __m128d *c_real, __m128d *c_imag)
{
    *c_real = _mm_mul_pd(a_imag, b_imag);
    *c_imag = _mm_mul_pd(a_imag, b_real);
    *c_real = sse_fmadd_pd(a_real, b_real, *c_real);
    *c_imag = sse_fmsub_pd(a_real, b_imag, *c_imag);
}

// c = a*b + c
static inline void cfmadd_pd(
        __m128d a_real, __m128d a_imag,
        __m128d b_real, __m128d b_imag,
        __m128d *c_real, __m128d *c_imag)
{
  *c_real = sse_fmsub_pd( a_imag, b_imag, *c_real );
  *c_imag = sse_fmadd_pd( a_imag, b_real, *c_imag );
  *c_real = sse_fmsub_pd( a_real, b_real, *c_real );
  *c_imag = sse_fmadd_pd( a_real, b_imag, *c_imag );
}


// c = conj(a)*b + c
static inline void cfmadd_conj_pd(
        __m128d a_real, __m128d a_imag,
        __m128d b_real, __m128d b_imag,
        __m128d *c_real, __m128d *c_imag)
{
    *c_real = sse_fmadd_pd(a_imag, b_imag, *c_real);
    *c_imag = sse_fmsub_pd(a_imag, b_real, *c_imag);
    *c_real = sse_fmadd_pd(a_real, b_real, *c_real);
    *c_imag = sse_fmsub_pd(a_real, b_imag, *c_imag);
}


static inline void sse_complex_deinterleaved_load_pd( double *data, __m128d *result_re, __m128d *result_im  ) {  
  *result_re = _mm_setr_pd( data[0], data[2] );
  *result_im = _mm_setr_pd( data[1], data[3] );
}


static inline void sse_complex_interleaved_store_pd( __m128d data_re, __m128d data_im, double *result  ) { 
  _mm_storeu_pd( result,                   _mm_unpacklo_pd( data_re, data_im ) );
  _mm_storeu_pd( result+SIMD_LENGTH_double, _mm_unpackhi_pd( data_re, data_im ) );
}

#endif
#endif
