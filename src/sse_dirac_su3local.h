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

#ifdef BOUNDARY  
  for ( int i=start; i<end; i+=2 ) {
#else
  for ( int i=start; i<end; i++ ) {
#endif
    
#ifdef MINUSDIR
#define CURRENT_SIGN (+1)
#define CURRENT_UPDATE _mm_add_ps
#define CURRENT_MUL cmul_conj
#define CURRENT_MADD cfmadd_conj
#else
#define CURRENT_SIGN (-1)
#define CURRENT_UPDATE _mm_sub_ps
#define CURRENT_MUL cmul
#define CURRENT_MADD cfmadd
#endif
    
    __m128 res1[2];
    __m128 res2[2];
    
    {
#ifdef BOUNDARY
      float *pt = phi+24*ind[i+1];
#else
#ifdef MINUSDIR
      float *pt = phi+24*ind[i];
#else
      float *pt = phi+24*neighbor[4*ind[i]+MU];
#endif
#endif
      // calc spin0 projection
      res1[0] = _mm_setr_ps( pt[0], pt[2], pt[4], 0 );
      res1[1] = _mm_setr_ps( index_re(pt,MU,0), index_re(pt+2,MU,0), index_re(pt+4,MU,0), 0 );
      res2[0] = _mm_setr_ps( pt[1], pt[3], pt[5], 0 );
      res2[1] = _mm_setr_ps( index_im(pt,MU,0), index_im(pt+2,MU,0), index_im(pt+4,MU,0), 0 );
      __m128 in1_re = CURRENT_UPDATE( res1[0], res1[1] );
      __m128 in1_im = CURRENT_UPDATE( res2[0], res2[1] );
      
      // calc spin1 projection
      res1[0] = _mm_setr_ps( pt[6], pt[8], pt[10], 0 );
      res1[1] = _mm_setr_ps( index_re(pt,MU,1), index_re(pt+2,MU,1), index_re(pt+4,MU,1), 0 );
      res2[0] = _mm_setr_ps( pt[7], pt[9], pt[11], 0 );
      res2[1] = _mm_setr_ps( index_im(pt,MU,1), index_im(pt+2,MU,1), index_im(pt+4,MU,1), 0 );   
      __m128 in2_re = CURRENT_UPDATE( res1[0], res1[1] );
      __m128 in2_im = CURRENT_UPDATE( res2[0], res2[1] );
      
      { // perform su(3) matrix vector multiplication
#ifdef BOUNDARY
#ifdef MINUSDIR
        pt = D + 96*ind[i+1] + 24*MU;
#else
        pt = D + 96*ind[i] + 24*MU;
#endif
#else
        pt = D + 96*ind[i] + 24*MU;
#endif
        // load 1st part of su(3) matrix and multiply
        {          
          __m128 buf1 = _mm_loadu_ps( pt );
          __m128 buf2 = _mm_loadu_ps( pt+SIMD_LENGTH_float );
          {
            __m128 buf3 = _mm_shuffle_ps( in1_re, in1_re, _MM_SHUFFLE(0,0,0,0) );
            __m128 buf4 = _mm_shuffle_ps( in1_im, in1_im, _MM_SHUFFLE(0,0,0,0) );
            CURRENT_MUL( buf1, buf2, buf3, buf4, &res1[0], &res1[1] );
          }
          {
            __m128 buf3 = _mm_shuffle_ps( in2_re, in2_re, _MM_SHUFFLE(0,0,0,0) );
            __m128 buf4 = _mm_shuffle_ps( in2_im, in2_im, _MM_SHUFFLE(0,0,0,0) );
            CURRENT_MUL( buf1, buf2, buf3, buf4, &res2[0], &res2[1] );
          }
        }
        // load 2nd part of su(3) matrix and multiply
        {
          __m128 buf1 = _mm_loadu_ps( pt+2*SIMD_LENGTH_float );
          __m128 buf2 = _mm_loadu_ps( pt+3*SIMD_LENGTH_float );
          {
            __m128 buf3 = _mm_shuffle_ps( in1_re, in1_re, _MM_SHUFFLE(1,1,1,1) );
            __m128 buf4 = _mm_shuffle_ps( in1_im, in1_im, _MM_SHUFFLE(1,1,1,1) );
            CURRENT_MADD( buf1, buf2, buf3, buf4, &res1[0], &res1[1] );
          }
          {
            __m128 buf3 = _mm_shuffle_ps( in2_re, in2_re, _MM_SHUFFLE(1,1,1,1) );
            __m128 buf4 = _mm_shuffle_ps( in2_im, in2_im, _MM_SHUFFLE(1,1,1,1) );
            CURRENT_MADD( buf1, buf2, buf3, buf4, &res2[0], &res2[1] );
          }
        }
        // load 3rd part of su(3) matrix and multiply
        {
          __m128 buf1 = _mm_loadu_ps( pt+4*SIMD_LENGTH_float );
          __m128 buf2 = _mm_loadu_ps( pt+5*SIMD_LENGTH_float );
          {
            __m128 buf3 = _mm_shuffle_ps( in1_re, in1_re, _MM_SHUFFLE(2,2,2,2) );
            __m128 buf4 = _mm_shuffle_ps( in1_im, in1_im, _MM_SHUFFLE(2,2,2,2) );
            CURRENT_MADD( buf1, buf2, buf3, buf4, &res1[0], &res1[1] );
          }
          {
            __m128 buf3 = _mm_shuffle_ps( in2_re, in2_re, _MM_SHUFFLE(2,2,2,2) );
            __m128 buf4 = _mm_shuffle_ps( in2_im, in2_im, _MM_SHUFFLE(2,2,2,2) );
            CURRENT_MADD( buf1, buf2, buf3, buf4, &res2[0], &res2[1] );
          }
        }
      }
    }
    
    { // store result
#ifdef BOUNDARY
      float *pt = eta+24*ind[i];
#else
#ifdef MINUSDIR
      float *pt = eta+24*neighbor[4*ind[i]+MU];
#else
      float *pt = eta+24*ind[i];
#endif
#endif
      {
        // store spin0 contribution
        __m128 buf1 = _mm_unpacklo_ps( res1[0], res1[1] );
        __m128 buf3 = _mm_loadu_ps( pt );
        __m128 buf2 = _mm_unpackhi_ps( res1[0], res1[1] );
        __m128 buf4 = _mm_loadu_ps( pt+4 );
        buf3 = UPD( buf3, buf1 );
        buf4 = UPD( buf4, buf2 );
        _mm_storeu_ps( pt, buf3 );
        _mm_storeu_ps( pt+4, buf4 );
        
        { // store contribution from 1st SU(3) multiplication to either spin2 or spin3
          __m128 *buf[2] = {&buf3,&buf4};
          *buf[gamma_offset[MU][0]] = _mm_set1_ps( CURRENT_SIGN*gamma_re_sign[MU][gamma_co[MU][0]] );
          *buf[1-gamma_offset[MU][0]] = _mm_set1_ps( CURRENT_SIGN*gamma_im_sign[MU][gamma_co[MU][0]] );
          buf1 = _mm_mul_ps( *buf[gamma_offset[MU][0]], res1[gamma_offset[MU][0]] );
          buf2 = _mm_mul_ps( *buf[1-gamma_offset[MU][0]], res1[1-gamma_offset[MU][0]] );
        }
        
        buf3 = _mm_unpacklo_ps( buf1, buf2 );
        buf4 = _mm_unpackhi_ps( buf1, buf2 );
        buf4 = _mm_shuffle_ps( buf3, buf4, _MM_SHUFFLE(1,0,3,2) ); // hi(buf3), lo(buf4)
        buf1 = _mm_loadu_ps( pt + 6*gamma_co[MU][0]   );
        buf2 = _mm_loadu_ps( pt + 6*gamma_co[MU][0]+2 );
        res1[0] = UPD( buf1, buf3 );
        res1[1] = UPD( buf2, buf4 );
        _mm_storeu_ps( pt + 6*gamma_co[MU][0], res1[0] );
        _mm_storeu_ps( pt + 6*gamma_co[MU][0]+2, res1[1] );
        
        // store spin1 contribution
        buf3 = _mm_unpacklo_ps( res2[0], res2[1] );
        buf1 = _mm_loadu_ps( pt+6 );
        buf4 = _mm_unpackhi_ps( res2[0], res2[1] );
        buf2 = _mm_loadu_ps( pt+10 );
        buf1 = UPD( buf1, buf3 );
        buf2 = UPD( buf2, buf4 );
        _mm_storeu_ps( pt+6, buf1 );
        _mm_storeu_ps( pt+10, buf2 );
        
        { // store contribution from 1st SU(3) multiplication to either spin2 or spin3
          __m128 *buf[2] = {&buf3,&buf4};
          *buf[gamma_offset[MU][1]] = _mm_set1_ps( CURRENT_SIGN*gamma_re_sign[MU][gamma_co[MU][1]] );
          *buf[1-gamma_offset[MU][1]] = _mm_set1_ps( CURRENT_SIGN*gamma_im_sign[MU][gamma_co[MU][1]] );
          buf1 = _mm_mul_ps( *buf[gamma_offset[MU][1]], res2[gamma_offset[MU][1]] );
          buf2 = _mm_mul_ps( *buf[1-gamma_offset[MU][1]], res2[1-gamma_offset[MU][1]] );
        }
        
        buf3 = _mm_unpacklo_ps( buf1, buf2 );
        buf4 = _mm_unpackhi_ps( buf1, buf2 );
        buf4 = _mm_shuffle_ps( buf3, buf4, _MM_SHUFFLE(1,0,3,2) ); // hi(buf3), lo(buf4)
        buf1 = _mm_loadu_ps( pt + 6*gamma_co[MU][1]   );
        buf2 = _mm_loadu_ps( pt + 6*gamma_co[MU][1]+2 );
        res2[0] = UPD( buf1, buf3 );
        res2[1] = UPD( buf2, buf4 );
        _mm_storeu_ps( pt + 6*gamma_co[MU][1], res2[0] );
        _mm_storeu_ps( pt + 6*gamma_co[MU][1]+2, res2[1] );
      }
    }
#undef CURRENT_SIGN
#undef CURRENT_UPDATE
#undef CURRENT_MUL
#undef CURRENT_MADD
#undef LINE1
#undef LINE2
#undef LINE3
  }