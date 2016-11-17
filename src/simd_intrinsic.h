/*
 * Copyright (C) 2016 Simone Bacchio.
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

#ifndef SIMD_INTRINSIC_HEADER
#define SIMD_INTRINSIC_HEADER

#if   defined AVX

#include "immintrin.h"
#include "xmmintrin.h"
#include "emmintrin.h"
#include "pmmintrin.h"

#define SIMD               _AVX
#define SIMD_LENGTH_float  8
#define SIMD_LENGTH_double 4

#define mm_float  __m256
#define mm_double __m256d

#define mm_mul_float  _mm256_mul_ps
#define mm_mul_double _mm256_mul_pd
#define mm_add_float  _mm256_add_ps
#define mm_add_double _mm256_add_pd
#define mm_sub_float  _mm256_sub_ps
#define mm_sub_double _mm256_sub_pd
#define mm_and_float  _mm256_and_ps
#define mm_and_double _mm256_and_pd

#define mm_setzero_float   _mm256_setzero_ps
#define mm_setzero_double  _mm256_setzero_pd
#define mm_setr_float      _mm256_setr_ps
#define mm_setr_double     _mm256_setr_pd
#define mm_set1_float      _mm256_set1_ps
#define mm_set1_double     _mm256_set1_pd
#define mm_load_float      _mm256_load_ps
#define mm_load_double     _mm256_load_pd
#define mm_unpacklo_float  _mm256_unpacklo_ps
#define mm_unpacklo_double _mm256_unpacklo_pd
#define mm_unpackhi_float  _mm256_unpackhi_ps
#define mm_unpackhi_double _mm256_unpackhi_pd
#define mm_store_float    _mm256_store_ps
#define mm_store_double   _mm256_store_pd

static inline mm_float mm_seteven_float( float *data ) {
  return mm_setr_float( data[0], data[2], data[4], data[6], data[8], data[10], data[12], data[14] );
}
static inline mm_double mm_seteven_double( double *data ) {
  return mm_setr_double( data[0], data[2], data[4], data[6] );
}

static inline float mm_reduce_add_float( mm_float v) {
  __m128 vlow  = _mm256_castps256_ps128(v);
  __m128 vhigh = _mm256_extractf128_ps(v, 1); // high 128
  vlow         = _mm_add_ps(vlow, vhigh);     // add the low 128
  // same of SSE
  __m128 shuf  = _mm_movehdup_ps(v);          // broadcast elements 3,1 to 2,0
  __m128 sums  = _mm_add_ps(v, shuf);
  shuf         = _mm_movehl_ps(shuf, sums);   // high half -> low half
  sums         = _mm_add_ss(sums, shuf);
  return       _mm_cvtss_f32(sums);
}
static inline double mm_reduce_add_double( mm_double v ) {
  __m128d vlow  = _mm256_castpd256_pd128(v);
  __m128d vhigh = _mm256_extractf128_pd(v, 1);
  vlow          = _mm_add_pd(vlow, vhigh);
  // same of SSE
  double tmp;
  _mm_storeh_pd(&tmp, vlow);        // store the high half
  return _mm_cvtsd_f64(vlow) + tmp; // cast the low half and sum
}

static inline void mm_transpose_float( mm_float *data ) {
  mm_float __t0, __t1, __t2, __t3, __t4, __t5, __t6, __t7;
  mm_float __tt0, __tt1, __tt2, __tt3, __tt4, __tt5, __tt6, __tt7;
  __t0 = _mm256_unpacklo_ps(data[0], data[1]);
  __t1 = _mm256_unpackhi_ps(data[0], data[1]);
  __t2 = _mm256_unpacklo_ps(data[2], data[3]);
  __t3 = _mm256_unpackhi_ps(data[2], data[3]);
  __t4 = _mm256_unpacklo_ps(data[4], data[5]);
  __t5 = _mm256_unpackhi_ps(data[4], data[5]);
  __t6 = _mm256_unpacklo_ps(data[6], data[7]);
  __t7 = _mm256_unpackhi_ps(data[6], data[7]);
  __tt0 = _mm256_shuffle_ps(__t0,__t2,_MM_SHUFFLE(1,0,1,0));
  __tt1 = _mm256_shuffle_ps(__t0,__t2,_MM_SHUFFLE(3,2,3,2));
  __tt2 = _mm256_shuffle_ps(__t1,__t3,_MM_SHUFFLE(1,0,1,0));
  __tt3 = _mm256_shuffle_ps(__t1,__t3,_MM_SHUFFLE(3,2,3,2));
  __tt4 = _mm256_shuffle_ps(__t4,__t6,_MM_SHUFFLE(1,0,1,0));
  __tt5 = _mm256_shuffle_ps(__t4,__t6,_MM_SHUFFLE(3,2,3,2));
  __tt6 = _mm256_shuffle_ps(__t5,__t7,_MM_SHUFFLE(1,0,1,0));
  __tt7 = _mm256_shuffle_ps(__t5,__t7,_MM_SHUFFLE(3,2,3,2));
  data[0] = _mm256_permute2f128_ps(__tt0, __tt4, 0x20);
  data[1] = _mm256_permute2f128_ps(__tt1, __tt5, 0x20);
  data[2] = _mm256_permute2f128_ps(__tt2, __tt6, 0x20);
  data[3] = _mm256_permute2f128_ps(__tt3, __tt7, 0x20);
  data[4] = _mm256_permute2f128_ps(__tt0, __tt4, 0x31);
  data[5] = _mm256_permute2f128_ps(__tt1, __tt5, 0x31);
  data[6] = _mm256_permute2f128_ps(__tt2, __tt6, 0x31);
  data[7] = _mm256_permute2f128_ps(__tt3, __tt7, 0x31);
}
static inline void mm_transpose_double( mm_double *data)
{
   mm_double tmp[4];

   tmp[0] = _mm256_unpacklo_pd( data[0], data[1] );
   tmp[1] = _mm256_unpacklo_pd( data[2], data[3] );
   tmp[2] = _mm256_unpackhi_pd( data[0], data[1] );
   tmp[3] = _mm256_unpackhi_pd( data[2], data[3] );
   //TODO
   data[0] = _mm256_movelh_pd( tmp[0], tmp[1] );
   data[1] = _mm256_movehl_pd( tmp[1], tmp[0] );
   data[2] = _mm256_movelh_pd( tmp[2], tmp[3] );
   data[3] = _mm256_movehl_pd( tmp[3], tmp[2] );
}
#elif defined SSE

#include "xmmintrin.h"
#include "emmintrin.h"
#include "pmmintrin.h"

#define SIMD               _SSE
#define SIMD_LENGTH_float  4
#define SIMD_LENGTH_double 2

#define mm_float  __m128
#define mm_double __m128d

#define mm_mul_float  _mm_mul_ps
#define mm_mul_double _mm_mul_pd
#define mm_add_float  _mm_add_ps
#define mm_add_double _mm_add_pd
#define mm_sub_float  _mm_sub_ps
#define mm_sub_double _mm_sub_pd
#define mm_and_float  _mm_and_ps
#define mm_and_double _mm_and_pd

#define mm_setzero_float   _mm_setzero_ps
#define mm_setzero_double  _mm_setzero_pd
#define mm_setr_float      _mm_setr_ps
#define mm_setr_double     _mm_setr_pd
#define mm_set1_float      _mm_set1_ps
#define mm_set1_double     _mm_set1_pd
#define mm_load_float      _mm_load_ps
#define mm_load_double     _mm_load_pd
#define mm_unpacklo_float  _mm_unpacklo_ps
#define mm_unpacklo_double _mm_unpacklo_pd
#define mm_unpackhi_float  _mm_unpackhi_ps
#define mm_unpackhi_double _mm_unpackhi_pd
#define mm_store_float    _mm_store_ps
#define mm_store_double   _mm_store_pd

static inline mm_float mm_seteven_float( float *data ) {
  return mm_setr_float( data[0], data[2], data[4], data[6] );
}
static inline mm_double mm_seteven_double( double *data ) {
  return mm_setr_double( data[0], data[2] );
}

static inline float mm_reduce_add_float( mm_float v ) {
  mm_float shuf = _mm_movehdup_ps(v);        // broadcast elements 3,1 to 2,0
  mm_float sums = _mm_add_ps(v, shuf);
  shuf          = _mm_movehl_ps(shuf, sums); // high half -> low half
  sums          = _mm_add_ss(sums, shuf);
  return        _mm_cvtss_f32(sums);
}
static inline double mm_reduce_add_double( mm_double v ) {
  double tmp;
  _mm_storeh_pd(&tmp, v);        // store the high half
  return _mm_cvtsd_f64(v) + tmp; // cast the low half and sum
}
#endif 

static inline mm_float mm_setodd_float( float *data ) {
  return mm_seteven_float( data+1 );
}
static inline mm_double mm_setodd_double( double *data ) {
  return mm_seteven_double( data+1 );
}

static inline void mm_transpose_float( mm_float *data ) {
  _MM_TRANSPOSE4_PS(data[0],data[1],data[2],data[3]);
}
static inline void mm_transpose_double( mm_double *data ) {
  double tmp01, tmp10 = _mm_cvtsd_f64(data[1]);
  _mm_storeh_pd(&tmp01, data[0]);
  _mm_loadl_pd(data[1], &tmp01);
  _mm_loadh_pd(data[0], &tmp10);
}

#ifdef FMA

#if   defined AVX

#define mm_fmadd_float   _mm256_fmadd_ps
#define mm_fmadd_double  _mm256_fmadd_pd
#define mm_fnmadd_float  _mm256_fnmadd_ps
#define mm_fnmadd_double _mm256_fnmadd_pd
#define mm_fmsub_float   _mm256_fmsub_ps
#define mm_fmsub_double  _mm256_fmsub_pd
#define mm_fnmsub_float  _mm256_fnmsub_ps
#define mm_fnmsub_double _mm256_fnmsub_pd

#elif defined SSE

#include "immintrin.h"

#define mm_fmadd_float   _mm_fmadd_ps
#define mm_fmadd_double  _mm_fmadd_pd
#define mm_fnmadd_float  _mm_fnmadd_ps
#define mm_fnmadd_double _mm_fnmadd_pd
#define mm_fmsub_float   _mm_fmsub_ps
#define mm_fmsub_double  _mm_fmsub_pd
#define mm_fnmsub_float  _mm_fnmsub_ps
#define mm_fnmsub_double _mm_fnmsub_pd

#endif

#else

// a*b + c
static inline mm_double mm_fmadd_double( mm_double a, mm_double b, mm_double c ) {
  return mm_add_double( mm_mul_double( a, b ), c );
}
static inline mm_float mm_fmadd_float( mm_float a, mm_float b, mm_float c ) {
  return mm_add_float( mm_mul_float( a, b ), c );
}

// -a*b + c
static inline mm_double mm_fnmadd_double( mm_double a, mm_double b, mm_double c ) {
  return mm_sub_double( c, mm_mul_double( a, b ) );
}
static inline mm_float mm_fnmadd_float( mm_float a, mm_float b, mm_float c ) {
  return mm_sub_float( c, mm_mul_float( a, b ) );
}

// a*b - c
static inline mm_double mm_fmsub_double( mm_double a, mm_double b, mm_double c ) {
  return mm_sub_double( mm_mul_double( a, b ), c );
}
static inline mm_float mm_fmsub_float( mm_float a, mm_float b, mm_float c ) {
  return mm_sub_float( mm_mul_float( a, b ), c );
}

// res = -a*b - c
static inline mm_double mm_fnmsub_double( mm_double a, mm_double b, mm_double c ) {

  mm_double na = mm_sub_double( mm_setzero_double(), a );
  return mm_sub_double( mm_mul_double( na, b ), c );
}
static inline mm_float mm_fnmsub_float( mm_float a, mm_float b, mm_float c ) {

  mm_float na = mm_sub_float( mm_setzero_float(), a );
  return mm_sub_float( mm_mul_float( na, b ), c );
}

#endif

#endif
