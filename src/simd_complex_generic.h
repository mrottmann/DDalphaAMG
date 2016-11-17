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

#ifndef SIMD_COMPLEX_PRECISION_HEADER
#define SIMD_COMPLEX_PRECISION_HEADER

// c = a*b
static inline void cmul_PRECISION(
        mm_PRECISION a_real, mm_PRECISION a_imag,
        mm_PRECISION b_real, mm_PRECISION b_imag,
        mm_PRECISION *c_real, mm_PRECISION *c_imag)
{
    *c_real = mm_mul_PRECISION(a_imag, b_imag);
    *c_imag = mm_mul_PRECISION(a_imag, b_real);
    *c_real = mm_fmsub_PRECISION(a_real, b_real, *c_real);
    *c_imag = mm_fmadd_PRECISION(a_real, b_imag, *c_imag);
}

// c = -a*b
static inline void cnmul_PRECISION(
        mm_PRECISION a_real, mm_PRECISION a_imag,
        mm_PRECISION b_real, mm_PRECISION b_imag,
        mm_PRECISION *c_real, mm_PRECISION *c_imag)
{
    *c_real = mm_mul_PRECISION(a_imag, b_imag);
    *c_imag = mm_mul_PRECISION(a_imag, b_real);
    *c_real = mm_fnmadd_PRECISION(a_real, b_real, *c_real);
    *c_imag = mm_fnmsub_PRECISION(a_real, b_imag, *c_imag);
}

// c = conj(a)*b
static inline void cmul_conj_PRECISION(
        mm_PRECISION a_real, mm_PRECISION a_imag,
        mm_PRECISION b_real, mm_PRECISION b_imag,
        mm_PRECISION *c_real, mm_PRECISION *c_imag)
{
    *c_real = mm_mul_PRECISION(a_imag, b_imag);
    *c_imag = mm_mul_PRECISION(a_imag, b_real);
    *c_real = mm_fmadd_PRECISION(a_real, b_real, *c_real);
    *c_imag = mm_fmsub_PRECISION(a_real, b_imag, *c_imag);
}

// c = -conj(a)*b
static inline void cnmul_conj_PRECISION(
        mm_PRECISION a_real, mm_PRECISION a_imag,
        mm_PRECISION b_real, mm_PRECISION b_imag,
        mm_PRECISION *c_real, mm_PRECISION *c_imag)
{
    *c_real = mm_mul_PRECISION(a_imag, b_imag);
    *c_imag = mm_mul_PRECISION(a_imag, b_real);
    *c_real = mm_fnmsub_PRECISION(a_real, b_real, *c_real);
    *c_imag = mm_fnmadd_PRECISION(a_real, b_imag, *c_imag);
}

// c = a*b + c
static inline void cfmadd_PRECISION(
        mm_PRECISION a_real, mm_PRECISION a_imag,
        mm_PRECISION b_real, mm_PRECISION b_imag,
        mm_PRECISION *c_real, mm_PRECISION *c_imag)
{
  *c_real = mm_fmsub_PRECISION( a_imag, b_imag, *c_real );
  *c_imag = mm_fmadd_PRECISION( a_imag, b_real, *c_imag );
  *c_real = mm_fmsub_PRECISION( a_real, b_real, *c_real );
  *c_imag = mm_fmadd_PRECISION( a_real, b_imag, *c_imag );
}

// c = -a*b + c
static inline void masked_cfnmadd_PRECISION(
        mm_PRECISION a_real, mm_PRECISION a_imag,
        mm_PRECISION b_real, mm_PRECISION b_imag,
        mm_PRECISION *c_real, mm_PRECISION *c_imag, mm_PRECISION mask )
{
    mm_PRECISION minus_a_real; mm_PRECISION minus_a_imag;
    minus_a_real = mm_setzero_PRECISION();
    minus_a_imag = mm_sub_PRECISION( minus_a_real, a_imag );
    minus_a_real = mm_sub_PRECISION( minus_a_real, a_real );
    minus_a_real = mm_and_PRECISION( minus_a_real, mask );
    minus_a_imag = mm_and_PRECISION( minus_a_imag, mask );
    
    *c_real = mm_fmsub_PRECISION(minus_a_imag, b_imag, *c_real);
    *c_imag = mm_fmadd_PRECISION(minus_a_imag, b_real, *c_imag);
    *c_real = mm_fmsub_PRECISION(minus_a_real, b_real, *c_real);
    *c_imag = mm_fmadd_PRECISION(minus_a_real, b_imag, *c_imag);
}

// c = -a*b + c
static inline void cfnmadd_PRECISION(
        mm_PRECISION a_real, mm_PRECISION a_imag,
        mm_PRECISION b_real, mm_PRECISION b_imag,
        mm_PRECISION *c_real, mm_PRECISION *c_imag)
{
    mm_PRECISION minus_a_real; mm_PRECISION minus_a_imag;
    minus_a_real = mm_setzero_PRECISION();
    minus_a_imag = mm_sub_PRECISION( minus_a_real, a_imag );
    minus_a_real = mm_sub_PRECISION( minus_a_real, a_real );
    
    *c_real = mm_fmsub_PRECISION(minus_a_imag, b_imag, *c_real);
    *c_imag = mm_fmadd_PRECISION(minus_a_imag, b_real, *c_imag);
    *c_real = mm_fmsub_PRECISION(minus_a_real, b_real, *c_real);
    *c_imag = mm_fmadd_PRECISION(minus_a_real, b_imag, *c_imag);
}

// c = conj(a)*b + c
static inline void cfmadd_conj_PRECISION(
        mm_PRECISION a_real, mm_PRECISION a_imag,
        mm_PRECISION b_real, mm_PRECISION b_imag,
        mm_PRECISION *c_real, mm_PRECISION *c_imag)
{
    *c_real = mm_fmadd_PRECISION(a_imag, b_imag, *c_real);
    *c_imag = mm_fmsub_PRECISION(a_imag, b_real, *c_imag);
    *c_real = mm_fmadd_PRECISION(a_real, b_real, *c_real);
    *c_imag = mm_fmsub_PRECISION(a_real, b_imag, *c_imag);
}

// c = -conj(a)*b + c
static inline void cfnmadd_conj_PRECISION(
        mm_PRECISION a_real, mm_PRECISION a_imag,
        mm_PRECISION b_real, mm_PRECISION b_imag,
        mm_PRECISION *c_real, mm_PRECISION *c_imag)
{
    mm_PRECISION minus_a_real; mm_PRECISION minus_a_imag;
    minus_a_real = mm_setzero_PRECISION();
    minus_a_imag = mm_sub_PRECISION( minus_a_real, a_imag );
    minus_a_real = mm_sub_PRECISION( minus_a_real, a_real );
  
    *c_real = mm_fmadd_PRECISION(minus_a_imag, b_imag, *c_real);
    *c_imag = mm_fmsub_PRECISION(minus_a_imag, b_real, *c_imag);
    *c_real = mm_fmadd_PRECISION(minus_a_real, b_real, *c_real);
    *c_imag = mm_fmsub_PRECISION(minus_a_real, b_imag, *c_imag);
}

static inline void cload_PRECISION( PRECISION *data, mm_PRECISION *result_re, mm_PRECISION *result_im  ) {  
  *result_re = mm_seteven_PRECISION( data );
  *result_im = mm_setodd_PRECISION( data );
}


static inline void cstore_PRECISION( PRECISION *result, mm_PRECISION data_re, mm_PRECISION data_im ) { 
  mm_store_PRECISION( result,                       mm_unpacklo_PRECISION( data_re, data_im ) );
  mm_store_PRECISION( result+SIMD_LENGTH_PRECISION, mm_unpackhi_PRECISION( data_re, data_im ) );
}

#endif
