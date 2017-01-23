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

#define index_re(phi,mu,spin) (gamma_re_sign[mu][spin]) * (phi)[ 6*gamma_co[mu][spin] + gamma_offset[mu][spin] ]
#define index_im(phi,mu,spin) (gamma_im_sign[mu][spin]) * (phi)[ 6*gamma_co[mu][spin] - gamma_offset[mu][spin] +1 ]

#define neighbor_coupling_file "sse_dirac_su3local.h"

void prp_double( complex_double *prn[4], complex_double *phi, int start, int end ) {   
  
  double *phi_pt = (double*)(phi+start);
  double *phi_end = (double*)(phi+end);
  double *pr[4] = {(double*)(prn[0]+start/2),(double*)(prn[1]+start/2),(double*)(prn[2]+start/2),(double*)(prn[3]+start/2)};
  
  while ( phi_pt < phi_end ) {
    
    __m128d phi_pt1_re; __m128d phi_pt1_im;
    
    sse_complex_deinterleaved_load_pd( phi_pt, &phi_pt1_re, &phi_pt1_im );
    for ( int mu=0; mu<4; mu++) {
      __m128d phi_pt2_re = _mm_setr_pd( index_re(phi_pt,mu,0), index_re(phi_pt+2,mu,0) );
      __m128d phi_pt2_im = _mm_setr_pd( index_im(phi_pt,mu,0), index_im(phi_pt+2,mu,0) );
      __m128d res_re = _mm_sub_pd( phi_pt1_re, phi_pt2_re );
      __m128d res_im = _mm_sub_pd( phi_pt1_im, phi_pt2_im );
      sse_complex_interleaved_store_pd( res_re, res_im, pr[mu] );
      pr[mu] += 2*SIMD_LENGTH_double;
    }
    
    sse_complex_deinterleaved_load_pd( phi_pt+4, &phi_pt1_re, &phi_pt1_im );
    for ( int mu=0; mu<4; mu++) {
      __m128d phi_pt2_re = _mm_setr_pd( index_re(phi_pt+4,mu,0), index_re(phi_pt,mu,1) );
      __m128d phi_pt2_im = _mm_setr_pd( index_im(phi_pt+4,mu,0), index_im(phi_pt,mu,1) );
      __m128d res_re = _mm_sub_pd( phi_pt1_re, phi_pt2_re );
      __m128d res_im = _mm_sub_pd( phi_pt1_im, phi_pt2_im );
      sse_complex_interleaved_store_pd( res_re, res_im, pr[mu] );
      pr[mu] += 2*SIMD_LENGTH_double;
    }
    
    sse_complex_deinterleaved_load_pd( phi_pt+8, &phi_pt1_re, &phi_pt1_im );
    for ( int mu=0; mu<4; mu++) {
      __m128d phi_pt2_re = _mm_setr_pd( index_re(phi_pt+2,mu,1), index_re(phi_pt+4,mu,1) );
      __m128d phi_pt2_im = _mm_setr_pd( index_im(phi_pt+2,mu,1), index_im(phi_pt+4,mu,1) );
      __m128d res_re = _mm_sub_pd( phi_pt1_re, phi_pt2_re );
      __m128d res_im = _mm_sub_pd( phi_pt1_im, phi_pt2_im );
      sse_complex_interleaved_store_pd( res_re, res_im, pr[mu] );
      pr[mu] += 2*SIMD_LENGTH_double;
    }
    
    phi_pt += 24;
  }
}


void prp_float( complex_float *prn[4], complex_float *phi, int start, int end ) { 
  
  float *phi_pt = (float*)(phi+start);
  float *phi_end = (float*)(phi+end);
  float *pr[4] = {(float*)(prn[0]+start/2),(float*)(prn[1]+start/2),(float*)(prn[2]+start/2),(float*)(prn[3]+start/2)};
  
  while ( phi_pt < phi_end ) {
    
    __m128 phi_pt1_re = _mm_setr_ps( phi_pt[0], phi_pt[2], phi_pt[4], phi_pt[6] );
    __m128 phi_pt1_im = _mm_setr_ps( phi_pt[1], phi_pt[3], phi_pt[5], phi_pt[7] );
    for ( int mu=0; mu<4; mu++) {
      __m128 phi_pt2_re = _mm_setr_ps( index_re(phi_pt,mu,0), index_re(phi_pt+2,mu,0),
                                       index_re(phi_pt+4,mu,0), index_re(phi_pt,mu,1) );
      __m128 phi_pt2_im = _mm_setr_ps( index_im(phi_pt,mu,0), index_im(phi_pt+2,mu,0),
                                       index_im(phi_pt+4,mu,0), index_im(phi_pt,mu,1) );
      
      __m128 res_re = _mm_sub_ps( phi_pt1_re, phi_pt2_re );
      __m128 res_im = _mm_sub_ps( phi_pt1_im, phi_pt2_im );
      
      sse_complex_interleaved_store( res_re, res_im, pr[mu] );
      pr[mu] += 8;
    }
    
    phi_pt1_re = _mm_setr_ps( phi_pt[8], phi_pt[10], phi_pt[24], phi_pt[26] );
    phi_pt1_im = _mm_setr_ps( phi_pt[9], phi_pt[11], phi_pt[25], phi_pt[27] );
    for ( int mu=0; mu<4; mu++) {
      __m128 phi_pt2_re = _mm_setr_ps( index_re(phi_pt+2,mu,1), index_re(phi_pt+4,mu,1),
                                       index_re(phi_pt+24,mu,0), index_re(phi_pt+26,mu,0) );
      __m128 phi_pt2_im = _mm_setr_ps( index_im(phi_pt+2,mu,1), index_im(phi_pt+4,mu,1),
                                       index_im(phi_pt+24,mu,0), index_im(phi_pt+26,mu,0) );
      __m128 res_re = _mm_sub_ps( phi_pt1_re, phi_pt2_re );
      __m128 res_im = _mm_sub_ps( phi_pt1_im, phi_pt2_im );
      
      sse_complex_interleaved_store( res_re, res_im, pr[mu] );
      pr[mu] += 8;
    }
    
    phi_pt1_re = _mm_setr_ps( phi_pt[28], phi_pt[30], phi_pt[32], phi_pt[34] );
    phi_pt1_im = _mm_setr_ps( phi_pt[29], phi_pt[31], phi_pt[33], phi_pt[35] );
    for ( int mu=0; mu<4; mu++) {
      __m128 phi_pt2_re = _mm_setr_ps( index_re(phi_pt+28,mu,0), index_re(phi_pt+24,mu,1),
                                       index_re(phi_pt+26,mu,1), index_re(phi_pt+28,mu,1) );
      __m128 phi_pt2_im = _mm_setr_ps( index_im(phi_pt+28,mu,0), index_im(phi_pt+24,mu,1),
                                       index_im(phi_pt+26,mu,1), index_im(phi_pt+28,mu,1) );
      __m128 res_re = _mm_sub_ps( phi_pt1_re, phi_pt2_re );
      __m128 res_im = _mm_sub_ps( phi_pt1_im, phi_pt2_im );
      
      sse_complex_interleaved_store( res_re, res_im, pr[mu] );
      pr[mu] += 8;
    }
    
    phi_pt+=48;
  }
}


void prn_su3_double( complex_double *prp[4], complex_double *phi, operator_double_struct *op, int *neighbor, int start, int end ) {
  
  double *phi_pt = (double*)(phi+start);
  double *phi_end_pt = (double*)(phi+end);
  double *pr[4] = {(double*)(prp[0]),(double*)(prp[1]),(double*)(prp[2]),(double*)(prp[3])};
  double *D_pt = ((double*)(op->D))+2*(start*3);
  int *nb_pt = neighbor+((start/12)*4);
  
  while ( phi_pt < phi_end_pt ) {
    
    __m128d in_re[3];
    __m128d in_im[3];
    
    for ( int i=0; i<3; i++ ) {
      in_re[i] = _mm_setr_pd( phi_pt[2*i+0], phi_pt[2*i+6] );
      in_im[i] = _mm_setr_pd( phi_pt[2*i+1], phi_pt[2*i+7] ); 
    }
    
    for ( int mu=0; mu<4; mu++ ) {
      
      __m128d v_re[3];
      __m128d v_im[3];
      
      // calc spin projection
      for ( int i=0; i<3; i++ )  {
        v_re[i] = _mm_setr_pd( index_re(phi_pt+2*i,mu,0), index_re(phi_pt+2*i,mu,1) );
        v_im[i] = _mm_setr_pd( index_im(phi_pt+2*i,mu,0), index_im(phi_pt+2*i,mu,1) );
        v_re[i] = _mm_add_pd( in_re[i], v_re[i] );
        v_im[i] = _mm_add_pd( in_im[i], v_im[i] );
      }
      
      {
        __m128d res_re[3];
        __m128d res_im[3];
        // load su(3) matrix and multiply
        for ( int i=0; i<3; i++ ) {
          __m128d buf_re = _mm_set1_pd( D_pt[0+2*i] );
          __m128d buf_im = _mm_set1_pd( D_pt[1+2*i] );
          cmul_conj_pd( buf_re, buf_im, v_re[0], v_im[0], &res_re[i], &res_im[i] );
          buf_re = _mm_set1_pd( D_pt[6+2*i] );
          buf_im = _mm_set1_pd( D_pt[7+2*i] );
          cfmadd_conj_pd( buf_re, buf_im, v_re[1], v_im[1], &res_re[i], &res_im[i] );
          buf_re = _mm_set1_pd( D_pt[12+2*i] );
          buf_im = _mm_set1_pd( D_pt[13+2*i] );
          cfmadd_conj_pd( buf_re, buf_im, v_re[2], v_im[2], &res_re[i], &res_im[i] );
        }
        
        {
          double *pr_pt = pr[mu]+2*6*(*(nb_pt));
          for ( int i=0; i<3; i++ ) {
            __m128d out1 = _mm_unpacklo_pd( res_re[i], res_im[i] );
            __m128d out2 = _mm_unpackhi_pd( res_re[i], res_im[i] );
            _mm_storeu_pd( pr_pt+0+2*i, out1 );
            _mm_storeu_pd( pr_pt+6+2*i, out2 );
          }
        }
      }
      
      D_pt += 18;
      nb_pt++;
    }
    
    phi_pt += 12*2;
  }
  
}


void prn_su3_float( complex_float *prp[4], complex_float *phi, operator_float_struct *op, int *neighbor, int start, int end ) { 
  
  float *phi_pt = (float*)(phi+start);
  float *phi_end_pt = (float*)(phi+end);
  float *pr[4] = {(float*)(prp[0]),(float*)(prp[1]),(float*)(prp[2]),(float*)(prp[3])};
  float *D_pt = (float*)(op->D_transformed_vectorized+2*(start*4));
  int *nb_pt = neighbor+((start/12)*4);
  
  while ( phi_pt < phi_end_pt ) { 
    
    __m128 in1[2];
    __m128 in2[2];

    in1[0] = _mm_setr_ps( phi_pt[0], phi_pt[2], phi_pt[4], 0 );
    in1[1] = _mm_setr_ps( phi_pt[1], phi_pt[3], phi_pt[5], 0 );
    in2[0] = _mm_setr_ps( phi_pt[6], phi_pt[8], phi_pt[10], 0 );
    in2[1] = _mm_setr_ps( phi_pt[7], phi_pt[9], phi_pt[11], 0 );
    
    for ( int mu=0; mu<4; mu++ ) {
      __m128 res1[2];
      __m128 res2[2];
      
      {
        // calc spin0 projection
        res1[0] = _mm_setr_ps( index_re(phi_pt,mu,0), index_re(phi_pt+2,mu,0), index_re(phi_pt+4,mu,0), 0 );
        res1[1] = _mm_setr_ps( index_im(phi_pt,mu,0), index_im(phi_pt+2,mu,0), index_im(phi_pt+4,mu,0), 0 );
        __m128 in1_re = _mm_add_ps( in1[0], res1[0] );
        __m128 in1_im = _mm_add_ps( in1[1], res1[1] );
        
        // calc spin1 projection
        res1[0] = _mm_setr_ps( index_re(phi_pt,mu,1), index_re(phi_pt+2,mu,1), index_re(phi_pt+4,mu,1), 0 );
        res1[1] = _mm_setr_ps( index_im(phi_pt,mu,1), index_im(phi_pt+2,mu,1), index_im(phi_pt+4,mu,1), 0 );   
        __m128 in2_re = _mm_add_ps( in2[0], res1[0] );
        __m128 in2_im = _mm_add_ps( in2[1], res1[1] );
        
        // load 1st part of su(3) matrix and multiply
        {
          __m128 buf1 = _mm_loadu_ps( D_pt );
          __m128 buf2 = _mm_loadu_ps( D_pt+SIMD_LENGTH_float );
          {
            __m128 buf3 = _mm_shuffle_ps( in1_re, in1_re, _MM_SHUFFLE(0,0,0,0) );
            __m128 buf4 = _mm_shuffle_ps( in1_im, in1_im, _MM_SHUFFLE(0,0,0,0) );
            cmul_conj( buf1, buf2, buf3, buf4, &res1[0], &res1[1] );
          }
          {
            __m128 buf3 = _mm_shuffle_ps( in2_re, in2_re, _MM_SHUFFLE(0,0,0,0) );
            __m128 buf4 = _mm_shuffle_ps( in2_im, in2_im, _MM_SHUFFLE(0,0,0,0) );
            cmul_conj( buf1, buf2, buf3, buf4, &res2[0], &res2[1] );
          }
        }
        // load 2nd part of su(3) matrix and multiply
        {
          __m128 buf1 = _mm_loadu_ps( D_pt+2*SIMD_LENGTH_float );
          __m128 buf2 = _mm_loadu_ps( D_pt+3*SIMD_LENGTH_float );
          {
            __m128 buf3 = _mm_shuffle_ps( in1_re, in1_re, _MM_SHUFFLE(1,1,1,1) );
            __m128 buf4 = _mm_shuffle_ps( in1_im, in1_im, _MM_SHUFFLE(1,1,1,1) );
            cfmadd_conj( buf1, buf2, buf3, buf4, &res1[0], &res1[1] );
          }
          {
            __m128 buf3 = _mm_shuffle_ps( in2_re, in2_re, _MM_SHUFFLE(1,1,1,1) );
            __m128 buf4 = _mm_shuffle_ps( in2_im, in2_im, _MM_SHUFFLE(1,1,1,1) );
            cfmadd_conj( buf1, buf2, buf3, buf4, &res2[0], &res2[1] );
          }
        }
        // load 3rd part of su(3) matrix and multiply
        {
          __m128 buf1 = _mm_loadu_ps( D_pt+4*SIMD_LENGTH_float );
          __m128 buf2 = _mm_loadu_ps( D_pt+5*SIMD_LENGTH_float );
          {
            __m128 buf3 = _mm_shuffle_ps( in1_re, in1_re, _MM_SHUFFLE(2,2,2,2) );
            __m128 buf4 = _mm_shuffle_ps( in1_im, in1_im, _MM_SHUFFLE(2,2,2,2) );
            cfmadd_conj( buf1, buf2, buf3, buf4, &res1[0], &res1[1] );
          }
          {
            __m128 buf3 = _mm_shuffle_ps( in2_re, in2_re, _MM_SHUFFLE(2,2,2,2) );
            __m128 buf4 = _mm_shuffle_ps( in2_im, in2_im, _MM_SHUFFLE(2,2,2,2) );
            cfmadd_conj( buf1, buf2, buf3, buf4, &res2[0], &res2[1] );
          }
        }
      }
      
      {
        __m128 buf1 = _mm_unpacklo_ps( res1[0], res1[1] );
        __m128 buf2 = _mm_unpackhi_ps( res1[0], res1[1] );
        __m128 buf3 = _mm_unpacklo_ps( res2[0], res2[1] );
        
        {
          __m128 buf4 = _mm_unpackhi_ps( res2[0], res2[1] );
          buf2 = _mm_movelh_ps( buf2, buf3 );
          buf3 = _mm_shuffle_ps( buf3, buf4, _MM_SHUFFLE(1,0,3,2) ); // hi(buf3), lo(buf4)
        }
        {
          float *pr_pt = pr[mu]+2*6*(*nb_pt);
          _mm_storeu_ps( pr_pt, buf1 );
          _mm_storeu_ps( pr_pt+4, buf2 );
          _mm_storeu_ps( pr_pt+8, buf3 );
        }
      }
      nb_pt++;
      D_pt += 24;
    }
    
    phi_pt += 24;
  }
}


void pbn_double( complex_double *eta, complex_double *prp[4], int start, int end ) {
  
  double *pr[4] = {(double*)(prp[0]),(double*)(prp[1]),(double*)(prp[2]),(double*)(prp[3])};
  double *eta_pt = (double*)(eta+start);
  
  __m128d gamma0[4];
  __m128d gamma1[4];
  
  for ( int mu=0; mu<4; mu++ ) {
    gamma0[mu] = _mm_setr_pd( gamma_re_sign[mu][gamma_co[mu][0]], gamma_im_sign[mu][gamma_co[mu][0]] );
    gamma1[mu] = _mm_setr_pd( gamma_re_sign[mu][gamma_co[mu][1]], gamma_im_sign[mu][gamma_co[mu][1]] );
  }
  
  for ( int i=start; i<end; i+=12 ) {
    
    __m128d res[12];
    for ( int j=0; j<12; j++ ) {
      res[j] = _mm_loadu_pd( eta_pt + 2*j );
    }
    
    __m128d in[6];
    // mu = T 
    for ( int j=0; j<6; j++ ) {
      in[j] = _mm_loadu_pd( pr[T] + i + 2*j );
      res[j] = _mm_sub_pd( res[j], in[j] );
    }
    for ( int j=0; j<3; j++ ) {
      __m128d buf1 = _mm_mul_pd( gamma0[T], GAMMA_T_SHUFFLE(in[j]) );
      res[3*gamma_co[T][0]+j] = _mm_sub_pd( res[3*gamma_co[T][0]+j], buf1 );
    }
    for ( int j=0; j<3; j++ ) {
      __m128d buf1 = _mm_mul_pd( gamma1[T], GAMMA_T_SHUFFLE(in[3+j]) );
      res[3*gamma_co[T][1]+j] = _mm_sub_pd( res[3*gamma_co[T][1]+j], buf1 );
    }
    // ---------------
    // mu = Z  
    for ( int j=0; j<6; j++ ) {
      in[j] = _mm_loadu_pd( pr[Z] + i + 2*j );
      res[j] = _mm_sub_pd( res[j], in[j] );
    }
    for ( int j=0; j<3; j++ ) {
      __m128d buf1 = _mm_mul_pd( gamma0[Z], GAMMA_Z_SHUFFLE(in[j]) );
      res[3*gamma_co[Z][0]+j] = _mm_sub_pd( res[3*gamma_co[Z][0]+j], buf1 );
    }
    for ( int j=0; j<3; j++ ) {
      __m128d buf1 = _mm_mul_pd( gamma1[Z], GAMMA_Z_SHUFFLE(in[3+j]) );
      res[3*gamma_co[Z][1]+j] = _mm_sub_pd( res[3*gamma_co[Z][1]+j], buf1 );
    }
    // ---------------
    // mu = Y
    for ( int j=0; j<6; j++ ) {
      in[j] = _mm_loadu_pd( pr[Y] + i + 2*j );
      res[j] = _mm_sub_pd( res[j], in[j] );
    }
    for ( int j=0; j<3; j++ ) {
      __m128d buf1 = _mm_mul_pd( gamma0[Y], GAMMA_Y_SHUFFLE(in[j]) );
      res[3*gamma_co[Y][0]+j] = _mm_sub_pd( res[3*gamma_co[Y][0]+j], buf1 );
    }
    for ( int j=0; j<3; j++ ) {
      __m128d buf1 = _mm_mul_pd( gamma1[Y], GAMMA_Y_SHUFFLE(in[3+j]) );
      res[3*gamma_co[Y][1]+j] = _mm_sub_pd( res[3*gamma_co[Y][1]+j], buf1 );
    }
    // ---------------
    // mu = X
    for ( int j=0; j<6; j++ ) {
      in[j] = _mm_loadu_pd( pr[X] + i + 2*j );
      res[j] = _mm_sub_pd( res[j], in[j] );
    }
    for ( int j=0; j<3; j++ ) {
      __m128d buf1 = _mm_mul_pd( gamma0[X], GAMMA_X_SHUFFLE(in[j]) );
      res[3*gamma_co[X][0]+j] = _mm_sub_pd( res[3*gamma_co[X][0]+j], buf1 );
    }
    for ( int j=0; j<3; j++ ) {
      __m128d buf1 = _mm_mul_pd( gamma1[X], GAMMA_X_SHUFFLE(in[3+j]) );
      res[3*gamma_co[X][1]+j] = _mm_sub_pd( res[3*gamma_co[X][1]+j], buf1 );
    }
    // ---------------
    
    for ( int j=0; j<12; j++ ) {
      _mm_storeu_pd( eta_pt + 2*j, res[j] );
    }
    eta_pt += 24;
  }
}


void pbn_float( complex_float *eta, complex_float *prp[4], int start, int end ) {
  
  float *pr[4] = {(float*)(prp[0]),(float*)(prp[1]),(float*)(prp[2]),(float*)(prp[3])};
  float *eta_pt = (float*)(eta+start);
  
  __m128 gamma0[4][2];
  __m128 gamma1[4][2];
  
  for ( int mu=0; mu<4; mu++ ) {
    gamma0[mu][0] = _mm_set1_ps( gamma_re_sign[mu][gamma_co[mu][0]] );
    gamma0[mu][1] = _mm_set1_ps( gamma_im_sign[mu][gamma_co[mu][0]] );
    gamma1[mu][0] = _mm_set1_ps( gamma_re_sign[mu][gamma_co[mu][1]] );
    gamma1[mu][1] = _mm_set1_ps( gamma_im_sign[mu][gamma_co[mu][1]] );
  }
  
  for ( int i=start; i<end; i+=12 ) {
    
    __m128 eta_lo1 = _mm_loadu_ps( eta_pt );
    __m128 eta_lo2 = _mm_loadu_ps( eta_pt + 4 );
    __m128 eta_hi1 = _mm_loadu_ps( eta_pt + 6 );
    __m128 eta_hi2 = _mm_loadu_ps( eta_pt + 10 );
    
    __m128 eta2_lo[2];
    __m128 eta2_hi[2];
    
    eta2_lo[0] = _mm_loadu_ps( eta_pt + 12 );
    eta2_hi[0] = _mm_loadu_ps( eta_pt + 14 );
    eta2_lo[1] = _mm_loadu_ps( eta_pt + 18 );
    eta2_hi[1] = _mm_loadu_ps( eta_pt + 20 );
    
    for ( int mu=0; mu<4; mu++ ) {
      __m128 res1[2];
      __m128 res2[2];
      
      res1[0] = _mm_setr_ps( *(pr[mu]+i), *(pr[mu]+i+2), *(pr[mu]+i+4), 0 );
      res1[1] = _mm_setr_ps( *(pr[mu]+i+1), *(pr[mu]+i+3), *(pr[mu]+i+5), 0 );
      
      res2[0] = _mm_setr_ps( *(pr[mu]+i+6), *(pr[mu]+i+8), *(pr[mu]+i+10), 0 );
      res2[1] = _mm_setr_ps( *(pr[mu]+i+7), *(pr[mu]+i+9), *(pr[mu]+i+11), 0 );
      
      {
        // store spin0 contribution
        {
          __m128 buf1 = _mm_unpacklo_ps( res1[0], res1[1] );
          __m128 buf2 = _mm_unpackhi_ps( res1[0], res1[1] );
          eta_lo1 = _mm_sub_ps( eta_lo1, buf1 );
          eta_lo2 = _mm_sub_ps( eta_lo2, buf2 );
        }
        
        // store contribution from 1st SU(3) multiplication to either spin2 or spin3
        __m128 buf1 = _mm_mul_ps( gamma0[mu][0], res1[gamma_offset[mu][0]] );
        __m128 buf2 = _mm_mul_ps( gamma0[mu][1], res1[1-gamma_offset[mu][0]] );
        __m128 buf3 = _mm_unpacklo_ps( buf1, buf2 );
        __m128 buf4 = _mm_unpackhi_ps( buf1, buf2 );
        buf4 = _mm_shuffle_ps( buf3, buf4, _MM_SHUFFLE(1,0,3,2) ); // hi(buf3), lo(buf4)
        eta2_lo[gamma_co[mu][2]] = _mm_sub_ps( eta2_lo[gamma_co[mu][2]], buf3 );
        eta2_hi[gamma_co[mu][2]] = _mm_sub_ps( eta2_hi[gamma_co[mu][2]], buf4 );
      }
      
      {
        // store spin1 contribution
        {
          __m128 buf1 = _mm_unpacklo_ps( res2[0], res2[1] );
          __m128 buf2 = _mm_unpackhi_ps( res2[0], res2[1] );
          eta_hi1 = _mm_sub_ps( eta_hi1, buf1 );
          eta_hi2 = _mm_sub_ps( eta_hi2, buf2 );
        }
        
        // store contribution from 1st SU(3) multiplication to either spin2 or spin3
        __m128 buf1 = _mm_mul_ps( gamma1[mu][0], res2[gamma_offset[mu][1]] );
        __m128 buf2 = _mm_mul_ps( gamma1[mu][1], res2[1-gamma_offset[mu][1]] );
        __m128 buf3 = _mm_unpacklo_ps( buf1, buf2 );
        __m128 buf4 = _mm_unpackhi_ps( buf1, buf2 );
        buf4 = _mm_shuffle_ps( buf3, buf4, _MM_SHUFFLE(1,0,3,2) ); // hi(buf3), lo(buf4)
        eta2_lo[gamma_co[mu][3]] = _mm_sub_ps( eta2_lo[gamma_co[mu][3]], buf3 );
        eta2_hi[gamma_co[mu][3]] = _mm_sub_ps( eta2_hi[gamma_co[mu][3]], buf4 );
      }
    }
    
    _mm_storeu_ps( eta_pt, eta_lo1 );
    _mm_storeu_ps( eta_pt+4, eta_lo2 );
    _mm_storeu_ps( eta_pt+6, eta_hi1 );
    _mm_storeu_ps( eta_pt+10, eta_hi2 );
    _mm_storeu_ps( eta_pt+12, eta2_lo[0] );
    _mm_storeu_ps( eta_pt+14, eta2_hi[0] );
    _mm_storeu_ps( eta_pt+18, eta2_lo[1] );
    _mm_storeu_ps( eta_pt+20, eta2_hi[1] );
    
    eta_pt += 24;
  }
}



void su3_pbp_double( complex_double* eta, complex_double *prn[4], operator_double_struct *op,
                     int *neighbor, int start, int end ) {  
  
  double *D_pt = ((double*)(op->D))+2*(start*3);
  double *eta_pt = (double*)(eta+start);
  double *eta_end_pt = (double*)(eta+end);
  double *pr[4] = {(double*)(prn[0]),(double*)(prn[1]),(double*)(prn[2]),(double*)(prn[3])};
  int *nb_pt = neighbor+((start/12)*4);
  
  __m128d gamma0[4];
  __m128d gamma1[4];
  
  for ( int mu=0; mu<4; mu++ ) {
    gamma0[mu] = _mm_setr_pd( -gamma_re_sign[mu][gamma_co[mu][0]], -gamma_im_sign[mu][gamma_co[mu][0]] );
    gamma1[mu] = _mm_setr_pd( -gamma_re_sign[mu][gamma_co[mu][1]], -gamma_im_sign[mu][gamma_co[mu][1]] );
  }
  
  while( eta_pt < eta_end_pt ) {
    
    __m128d res[12];
    for ( int i=0; i<12; i++ ) {
      res[i] = _mm_loadu_pd( eta_pt + 2*i );
    }
    
    // ---------------
    // mu = T
    {
      __m128d res_re[3];
      __m128d res_im[3];
      {
        __m128d v_re[3];
        __m128d v_im[3];
        int j = 2*6*(*nb_pt);
        
        for ( int i=0; i<3; i++ )  {
          v_re[i] = _mm_setr_pd( *(pr[T]+j+0+2*i), *(pr[T]+j+6+2*i) );
          v_im[i] = _mm_setr_pd( *(pr[T]+j+1+2*i), *(pr[T]+j+7+2*i) );
        }
        for ( int i=0; i<3; i++ ) {
          __m128d buf_re = _mm_set1_pd( D_pt[0+6*i] );
          __m128d buf_im = _mm_set1_pd( D_pt[1+6*i] );
          cmul_pd( buf_re, buf_im, v_re[0], v_im[0], &res_re[i], &res_im[i] );
          buf_re = _mm_set1_pd( D_pt[2+6*i] );
          buf_im = _mm_set1_pd( D_pt[3+6*i] );
          cfmadd_pd( buf_re, buf_im, v_re[1], v_im[1], &res_re[i], &res_im[i] );
          buf_re = _mm_set1_pd( D_pt[4+6*i] );
          buf_im = _mm_set1_pd( D_pt[5+6*i] );
          cfmadd_pd( buf_re, buf_im, v_re[2], v_im[2], &res_re[i], &res_im[i] );
        }
        D_pt += 18;
        nb_pt++;
      }
      {
        __m128d in[6];
        in[0] = _mm_unpacklo_pd( res_re[0], res_im[0] );
        in[1] = _mm_unpacklo_pd( res_re[1], res_im[1] );
        in[2] = _mm_unpacklo_pd( res_re[2], res_im[2] );
        
        in[3] = _mm_unpackhi_pd( res_re[0], res_im[0] );
        in[4] = _mm_unpackhi_pd( res_re[1], res_im[1] );
        in[5] = _mm_unpackhi_pd( res_re[2], res_im[2] );
        
        for ( int i=0; i<6; i++ ) {
          res[i] = _mm_sub_pd( res[i], in[i] );
        }
        for ( int i=0; i<3; i++ ) {
          __m128d buf1 = _mm_mul_pd( gamma0[T], GAMMA_T_SHUFFLE(in[i]) );
          res[3*gamma_co[T][0]+i] = _mm_sub_pd( res[3*gamma_co[T][0]+i], buf1 );
        }
        for ( int i=0; i<3; i++ ) {
          __m128d buf1 = _mm_mul_pd( gamma1[T], GAMMA_T_SHUFFLE(in[3+i]) );
          res[3*gamma_co[T][1]+i] = _mm_sub_pd( res[3*gamma_co[T][1]+i], buf1 );
        }
      }
    }
    // ---------------
    // mu = Z
    {
      __m128d res_re[3];
      __m128d res_im[3];
      {
        __m128d v_re[3];
        __m128d v_im[3];
        int j = 2*6*(*nb_pt);
        
        for ( int i=0; i<3; i++ )  {
          v_re[i] = _mm_setr_pd( *(pr[Z]+j+0+2*i), *(pr[Z]+j+6+2*i) );
          v_im[i] = _mm_setr_pd( *(pr[Z]+j+1+2*i), *(pr[Z]+j+7+2*i) );
        }
        for ( int i=0; i<3; i++ ) {
          __m128d buf_re = _mm_set1_pd( D_pt[0+6*i] );
          __m128d buf_im = _mm_set1_pd( D_pt[1+6*i] );
          cmul_pd( buf_re, buf_im, v_re[0], v_im[0], &res_re[i], &res_im[i] );
          buf_re = _mm_set1_pd( D_pt[2+6*i] );
          buf_im = _mm_set1_pd( D_pt[3+6*i] );
          cfmadd_pd( buf_re, buf_im, v_re[1], v_im[1], &res_re[i], &res_im[i] );
          buf_re = _mm_set1_pd( D_pt[4+6*i] );
          buf_im = _mm_set1_pd( D_pt[5+6*i] );
          cfmadd_pd( buf_re, buf_im, v_re[2], v_im[2], &res_re[i], &res_im[i] );
        }
        D_pt += 18;
        nb_pt++;
      }
      {
        __m128d in[6];
        in[0] = _mm_unpacklo_pd( res_re[0], res_im[0] );
        in[1] = _mm_unpacklo_pd( res_re[1], res_im[1] );
        in[2] = _mm_unpacklo_pd( res_re[2], res_im[2] );
        
        in[3] = _mm_unpackhi_pd( res_re[0], res_im[0] );
        in[4] = _mm_unpackhi_pd( res_re[1], res_im[1] );
        in[5] = _mm_unpackhi_pd( res_re[2], res_im[2] );
        
        for ( int i=0; i<6; i++ ) {
          res[i] = _mm_sub_pd( res[i], in[i] );
        }
        for ( int i=0; i<3; i++ ) {
          __m128d buf1 = _mm_mul_pd( gamma0[Z], GAMMA_Z_SHUFFLE(in[i]) );
          res[3*gamma_co[Z][0]+i] = _mm_sub_pd( res[3*gamma_co[Z][0]+i], buf1 );
        }
        for ( int i=0; i<3; i++ ) {
          __m128d buf1 = _mm_mul_pd( gamma1[Z], GAMMA_Z_SHUFFLE(in[3+i]) );
          res[3*gamma_co[Z][1]+i] = _mm_sub_pd( res[3*gamma_co[Z][1]+i], buf1 );
        }
      }
    }
    // ---------------
    // mu = Y
    {
      __m128d res_re[3];
      __m128d res_im[3];
      {
        __m128d v_re[3];
        __m128d v_im[3];
        int j = 2*6*(*nb_pt);
        
        for ( int i=0; i<3; i++ )  {
          v_re[i] = _mm_setr_pd( *(pr[Y]+j+0+2*i), *(pr[Y]+j+6+2*i) );
          v_im[i] = _mm_setr_pd( *(pr[Y]+j+1+2*i), *(pr[Y]+j+7+2*i) );
        }
        for ( int i=0; i<3; i++ ) {
          __m128d buf_re = _mm_set1_pd( D_pt[0+6*i] );
          __m128d buf_im = _mm_set1_pd( D_pt[1+6*i] );
          cmul_pd( buf_re, buf_im, v_re[0], v_im[0], &res_re[i], &res_im[i] );
          buf_re = _mm_set1_pd( D_pt[2+6*i] );
          buf_im = _mm_set1_pd( D_pt[3+6*i] );
          cfmadd_pd( buf_re, buf_im, v_re[1], v_im[1], &res_re[i], &res_im[i] );
          buf_re = _mm_set1_pd( D_pt[4+6*i] );
          buf_im = _mm_set1_pd( D_pt[5+6*i] );
          cfmadd_pd( buf_re, buf_im, v_re[2], v_im[2], &res_re[i], &res_im[i] );
        }
        D_pt += 18;
        nb_pt++;
      }
      {
        __m128d in[6];
        in[0] = _mm_unpacklo_pd( res_re[0], res_im[0] );
        in[1] = _mm_unpacklo_pd( res_re[1], res_im[1] );
        in[2] = _mm_unpacklo_pd( res_re[2], res_im[2] );
        
        in[3] = _mm_unpackhi_pd( res_re[0], res_im[0] );
        in[4] = _mm_unpackhi_pd( res_re[1], res_im[1] );
        in[5] = _mm_unpackhi_pd( res_re[2], res_im[2] );
        
        for ( int i=0; i<6; i++ ) {
          res[i] = _mm_sub_pd( res[i], in[i] );
        }
        for ( int i=0; i<3; i++ ) {
          __m128d buf1 = _mm_mul_pd( gamma0[Y], GAMMA_Y_SHUFFLE(in[i]) );
          res[3*gamma_co[Y][0]+i] = _mm_sub_pd( res[3*gamma_co[Y][0]+i], buf1 );
        }
        for ( int i=0; i<3; i++ ) {
          __m128d buf1 = _mm_mul_pd( gamma1[Y], GAMMA_Y_SHUFFLE(in[3+i]) );
          res[3*gamma_co[Y][1]+i] = _mm_sub_pd( res[3*gamma_co[Y][1]+i], buf1 );
        }
      }
    }
    // ---------------  
    // mu = X
    {
      __m128d res_re[3];
      __m128d res_im[3];
      {
        __m128d v_re[3];
        __m128d v_im[3];
        int j = 2*6*(*nb_pt);
        
        for ( int i=0; i<3; i++ )  {
          v_re[i] = _mm_setr_pd( *(pr[X]+j+0+2*i), *(pr[X]+j+6+2*i) );
          v_im[i] = _mm_setr_pd( *(pr[X]+j+1+2*i), *(pr[X]+j+7+2*i) );
        }
        for ( int i=0; i<3; i++ ) {
          __m128d buf_re = _mm_set1_pd( D_pt[0+6*i] );
          __m128d buf_im = _mm_set1_pd( D_pt[1+6*i] );
          cmul_pd( buf_re, buf_im, v_re[0], v_im[0], &res_re[i], &res_im[i] );
          buf_re = _mm_set1_pd( D_pt[2+6*i] );
          buf_im = _mm_set1_pd( D_pt[3+6*i] );
          cfmadd_pd( buf_re, buf_im, v_re[1], v_im[1], &res_re[i], &res_im[i] );
          buf_re = _mm_set1_pd( D_pt[4+6*i] );
          buf_im = _mm_set1_pd( D_pt[5+6*i] );
          cfmadd_pd( buf_re, buf_im, v_re[2], v_im[2], &res_re[i], &res_im[i] );
        }
        D_pt += 18;
        nb_pt++;
      }
      {
        __m128d in[6];
        in[0] = _mm_unpacklo_pd( res_re[0], res_im[0] );
        in[1] = _mm_unpacklo_pd( res_re[1], res_im[1] );
        in[2] = _mm_unpacklo_pd( res_re[2], res_im[2] );
        
        in[3] = _mm_unpackhi_pd( res_re[0], res_im[0] );
        in[4] = _mm_unpackhi_pd( res_re[1], res_im[1] );
        in[5] = _mm_unpackhi_pd( res_re[2], res_im[2] );
        
        for ( int i=0; i<6; i++ ) {
          res[i] = _mm_sub_pd( res[i], in[i] );
        }
        for ( int i=0; i<3; i++ ) {
          __m128d buf1 = _mm_mul_pd( gamma0[X], GAMMA_X_SHUFFLE(in[i]) );
          res[3*gamma_co[X][0]+i] = _mm_sub_pd( res[3*gamma_co[X][0]+i], buf1 );
        }
        for ( int i=0; i<3; i++ ) {
          __m128d buf1 = _mm_mul_pd( gamma1[X], GAMMA_X_SHUFFLE(in[3+i]) );
          res[3*gamma_co[X][1]+i] = _mm_sub_pd( res[3*gamma_co[X][1]+i], buf1 );
        }
      }
    }
    // ---------------
    
    for ( int i=0; i<12; i++ ) {
      _mm_storeu_pd( eta_pt + 2*i, res[i] );
    }
    eta_pt+=24;
  }
  
}


void su3_pbp_float( complex_float* eta, complex_float *prn[4], operator_float_struct *op,
                     int *neighbor, int start, int end ) {
  
  float *D_pt = (float*)(op->D_vectorized+2*(start*4));
  float *eta_pt = (float*)(eta+start);
  float *eta_end_pt = (float*)(eta+end);
  float *pr[4] = {(float*)(prn[0]),(float*)(prn[1]),(float*)(prn[2]),(float*)(prn[3])};
  int *nb_pt = neighbor+((start/12)*4);
  
  __m128 gamma0[4][2];
  __m128 gamma1[4][2];
  
  for ( int mu=0; mu<4; mu++ ) {
    gamma0[mu][0] = _mm_set1_ps( -gamma_re_sign[mu][gamma_co[mu][0]] );
    gamma0[mu][1] = _mm_set1_ps( -gamma_im_sign[mu][gamma_co[mu][0]] );
    gamma1[mu][0] = _mm_set1_ps( -gamma_re_sign[mu][gamma_co[mu][1]] );
    gamma1[mu][1] = _mm_set1_ps( -gamma_im_sign[mu][gamma_co[mu][1]] );
  }
  
  while( eta_pt < eta_end_pt ) {
    
    __m128 eta_lo1 = _mm_loadu_ps( eta_pt );
    __m128 eta_lo2 = _mm_loadu_ps( eta_pt + 4 );
    __m128 eta_hi1 = _mm_loadu_ps( eta_pt + 6 );
    __m128 eta_hi2 = _mm_loadu_ps( eta_pt + 10 );
    
    __m128 eta2_lo[2];
    __m128 eta2_hi[2];
    
    eta2_lo[0] = _mm_loadu_ps( eta_pt + 12 );
    eta2_hi[0] = _mm_loadu_ps( eta_pt + 14 );
    eta2_lo[1] = _mm_loadu_ps( eta_pt + 18 );
    eta2_hi[1] = _mm_loadu_ps( eta_pt + 20 );
    
    for ( int mu=0; mu<4; mu++ ) {
      __m128 res1[2];
      __m128 res2[2];
      
      {
        int j = 2*6*(*nb_pt);
        // load 1st part of su(3) matrix and multiply
        {
          __m128 buf1 = _mm_loadu_ps( D_pt );
          __m128 buf2 = _mm_loadu_ps( D_pt+SIMD_LENGTH_float );
          {
            __m128 buf3 = _mm_set1_ps( *(pr[mu]+j+0) );
            __m128 buf4 = _mm_set1_ps( *(pr[mu]+j+1) );
            cmul( buf1, buf2, buf3, buf4, &res1[0], &res1[1] );
          }
          {
            __m128 buf3 = _mm_set1_ps( *(pr[mu]+j+6) );
            __m128 buf4 = _mm_set1_ps( *(pr[mu]+j+7) );
            cmul( buf1, buf2, buf3, buf4, &res2[0], &res2[1] );
          }
        }
        // load 2nd part of su(3) matrix and multiply
        {
          __m128 buf1 = _mm_loadu_ps( D_pt+2*SIMD_LENGTH_float );
          __m128 buf2 = _mm_loadu_ps( D_pt+3*SIMD_LENGTH_float );
          {
            __m128 buf3 = _mm_set1_ps( *(pr[mu]+j+2) );
            __m128 buf4 = _mm_set1_ps( *(pr[mu]+j+3) );
            cfmadd( buf1, buf2, buf3, buf4, &res1[0], &res1[1] );
          }
          {
            __m128 buf3 = _mm_set1_ps( *(pr[mu]+j+8) );
            __m128 buf4 = _mm_set1_ps( *(pr[mu]+j+9) );
            cfmadd( buf1, buf2, buf3, buf4, &res2[0], &res2[1] );
          }
        }
        // load 3rd part of su(3) matrix and multiply
        {
          __m128 buf1 = _mm_loadu_ps( D_pt+4*SIMD_LENGTH_float );
          __m128 buf2 = _mm_loadu_ps( D_pt+5*SIMD_LENGTH_float );
          {
            __m128 buf3 = _mm_set1_ps( *(pr[mu]+j+4) );
            __m128 buf4 = _mm_set1_ps( *(pr[mu]+j+5) );
            cfmadd( buf1, buf2, buf3, buf4, &res1[0], &res1[1] );
          }
          {
            __m128 buf3 = _mm_set1_ps( *(pr[mu]+j+10) );
            __m128 buf4 = _mm_set1_ps( *(pr[mu]+j+11) );
            cfmadd( buf1, buf2, buf3, buf4, &res2[0], &res2[1] );
          }
        }
      }
            
      {
        // store spin0 contribution
        {
          __m128 buf1 = _mm_unpacklo_ps( res1[0], res1[1] );
          __m128 buf2 = _mm_unpackhi_ps( res1[0], res1[1] );
          eta_lo1 = _mm_sub_ps( eta_lo1, buf1 );
          eta_lo2 = _mm_sub_ps( eta_lo2, buf2 );
        }
        
        // store contribution from 1st SU(3) multiplication to either spin2 or spin3
        __m128 buf1 = _mm_mul_ps( gamma0[mu][0], res1[gamma_offset[mu][0]] );
        __m128 buf2 = _mm_mul_ps( gamma0[mu][1], res1[1-gamma_offset[mu][0]] );
        __m128 buf3 = _mm_unpacklo_ps( buf1, buf2 );
        __m128 buf4 = _mm_unpackhi_ps( buf1, buf2 );
        buf4 = _mm_shuffle_ps( buf3, buf4, _MM_SHUFFLE(1,0,3,2) ); // hi(buf3), lo(buf4)
        eta2_lo[gamma_co[mu][2]] = _mm_sub_ps( eta2_lo[gamma_co[mu][2]], buf3 );
        eta2_hi[gamma_co[mu][2]] = _mm_sub_ps( eta2_hi[gamma_co[mu][2]], buf4 );
      }
      
      {
        // store spin1 contribution
        {
          __m128 buf1 = _mm_unpacklo_ps( res2[0], res2[1] );
          __m128 buf2 = _mm_unpackhi_ps( res2[0], res2[1] );
          eta_hi1 = _mm_sub_ps( eta_hi1, buf1 );
          eta_hi2 = _mm_sub_ps( eta_hi2, buf2 );
        }
        
        // store contribution from 1st SU(3) multiplication to either spin2 or spin3
        __m128 buf1 = _mm_mul_ps( gamma1[mu][0], res2[gamma_offset[mu][1]] );
        __m128 buf2 = _mm_mul_ps( gamma1[mu][1], res2[1-gamma_offset[mu][1]] );
        __m128 buf3 = _mm_unpacklo_ps( buf1, buf2 );
        __m128 buf4 = _mm_unpackhi_ps( buf1, buf2 );
        buf4 = _mm_shuffle_ps( buf3, buf4, _MM_SHUFFLE(1,0,3,2) ); // hi(buf3), lo(buf4)
        eta2_lo[gamma_co[mu][3]] = _mm_sub_ps( eta2_lo[gamma_co[mu][3]], buf3 );
        eta2_hi[gamma_co[mu][3]] = _mm_sub_ps( eta2_hi[gamma_co[mu][3]], buf4 );
      }
      
      nb_pt++;
      D_pt += 24;
    }
    
    _mm_storeu_ps( eta_pt, eta_lo1 );
    _mm_storeu_ps( eta_pt+4, eta_lo2 );
    _mm_storeu_ps( eta_pt+6, eta_hi1 );
    _mm_storeu_ps( eta_pt+10, eta_hi2 );
    _mm_storeu_ps( eta_pt+12, eta2_lo[0] );
    _mm_storeu_ps( eta_pt+14, eta2_hi[0] );
    _mm_storeu_ps( eta_pt+18, eta2_lo[1] );
    _mm_storeu_ps( eta_pt+20, eta2_hi[1] );
  
    eta_pt += 24;
  }
  
}



void block_oddeven_plus_coupling_double( double *eta, double *D, double *phi, int mu, int start, int end, int *ind, int *neighbor ) { }
void block_oddeven_pT_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_sub_ps
#define MU 0
#include neighbor_coupling_file
#undef MU
#undef UPD
}
void block_oddeven_pZ_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_sub_ps
#define MU 1
#include neighbor_coupling_file
#undef MU
#undef UPD
}
void block_oddeven_pY_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_sub_ps
#define MU 2
#include neighbor_coupling_file
#undef MU
#undef UPD
}
void block_oddeven_pX_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_sub_ps
#define MU 3
#include neighbor_coupling_file
#undef MU
#undef UPD
}
void block_oddeven_plus_coupling_float( float *eta, float *D, float *phi, int mu, int start, int end, int *ind, int *neighbor ) {
  switch ( mu ) {
    case T: block_oddeven_pT_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    case Z: block_oddeven_pZ_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    case Y: block_oddeven_pY_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    case X: block_oddeven_pX_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    default: error0("block_oddeven_plus_coupling_float: invalid mu=%d\n", mu );
  }
}


void block_oddeven_nplus_coupling_double( double *eta, double *D, double *phi, int mu, int start, int end, int *ind, int *neighbor ) { }
void block_oddeven_npT_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_add_ps
#define MU 0
#include neighbor_coupling_file
#undef MU
#undef UPD
}
void block_oddeven_npZ_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_add_ps
#define MU 1
#include neighbor_coupling_file
#undef MU
#undef UPD
}
void block_oddeven_npY_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_add_ps
#define MU 2
#include neighbor_coupling_file
#undef MU
#undef UPD
}
void block_oddeven_npX_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_add_ps
#define MU 3
#include neighbor_coupling_file
#undef MU
#undef UPD
}
void block_oddeven_nplus_coupling_float( float *eta, float *D, float *phi, int mu, int start, int end, int *ind, int *neighbor ) {
  switch ( mu ) {
    case T: block_oddeven_npT_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    case Z: block_oddeven_npZ_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    case Y: block_oddeven_npY_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    case X: block_oddeven_npX_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    default: error0("block_oddeven_nplus_coupling_float: invalid mu=%d\n", mu );
  }
}


void block_oddeven_minus_coupling_double( double *eta, double *D, double *phi, int mu, int start, int end, int *ind, int *neighbor ) { }
void block_oddeven_mT_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_sub_ps
#define MU 0
#define MINUSDIR
#include neighbor_coupling_file
#undef MINUSDIR
#undef MU
#undef UPD
}
void block_oddeven_mZ_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_sub_ps
#define MU 1
#define MINUSDIR
#include neighbor_coupling_file
#undef MINUSDIR
#undef MU
#undef UPD
}
void block_oddeven_mY_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_sub_ps
#define MU 2
#define MINUSDIR
#include neighbor_coupling_file
#undef MINUSDIR
#undef MU
#undef UPD
}
void block_oddeven_mX_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_sub_ps
#define MU 3
#define MINUSDIR
#include neighbor_coupling_file
#undef MINUSDIR
#undef MU
#undef UPD
}
void block_oddeven_minus_coupling_float( float *eta, float *D, float *phi, int mu, int start, int end, int *ind, int *neighbor ) {
  switch ( mu ) {
    case T: block_oddeven_mT_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    case Z: block_oddeven_mZ_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    case Y: block_oddeven_mY_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    case X: block_oddeven_mX_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    default: error0("block_oddeven_minus_coupling_float: invalid mu=%d\n", mu );
  }
}


void block_oddeven_nminus_coupling_double( double *eta, double *D, double *phi, int mu, int start, int end, int *ind, int *neighbor ) { }
void block_oddeven_nmT_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_add_ps
#define MU 0
#define MINUSDIR
#include neighbor_coupling_file
#undef MINUSDIR
#undef MU
#undef UPD
}
void block_oddeven_nmZ_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_add_ps
#define MU 1
#define MINUSDIR
#include neighbor_coupling_file
#undef MINUSDIR
#undef MU
#undef UPD
}
void block_oddeven_nmY_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_add_ps
#define MU 2
#define MINUSDIR
#include neighbor_coupling_file
#undef MINUSDIR
#undef MU
#undef UPD
}
void block_oddeven_nmX_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_add_ps
#define MU 3
#define MINUSDIR
#include neighbor_coupling_file
#undef MINUSDIR
#undef MU
#undef UPD
}
void block_oddeven_nminus_coupling_float( float *eta, float *D, float *phi, int mu, int start, int end, int *ind, int *neighbor ) {
  switch ( mu ) {
    case T: block_oddeven_nmT_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    case Z: block_oddeven_nmZ_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    case Y: block_oddeven_nmY_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    case X: block_oddeven_nmX_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    default: error0("block_oddeven_nminus_coupling_float: invalid mu=%d\n", mu );
  }
}


void boundary_nminus_coupling_double( double *eta, double *D, double *phi, int mu, int start, int end, int *ind, int *neighbor ) { }
void boundary_nmT_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_add_ps
#define MU 0
#define BOUNDARY
#define MINUSDIR
#include neighbor_coupling_file
#undef MINUSDIR
#undef BOUNDARY
#undef MU
#undef UPD
}
void boundary_nmZ_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_add_ps
#define MU 1
#define BOUNDARY
#define MINUSDIR
#include neighbor_coupling_file
#undef MINUSDIR
#undef BOUNDARY
#undef MU
#undef UPD
}
void boundary_nmY_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_add_ps
#define MU 2
#define BOUNDARY
#define MINUSDIR
#include neighbor_coupling_file
#undef MINUSDIR
#undef BOUNDARY
#undef MU
#undef UPD
}
void boundary_nmX_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_add_ps
#define MU 3
#define BOUNDARY
#define MINUSDIR
#include neighbor_coupling_file
#undef MINUSDIR
#undef BOUNDARY
#undef MU
#undef UPD
}
void boundary_nminus_coupling_float( float *eta, float *D, float *phi, int mu, int start, int end, int *ind, int *neighbor ) {
  switch ( mu ) {
    case T: boundary_nmT_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    case Z: boundary_nmZ_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    case Y: boundary_nmY_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    case X: boundary_nmX_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    default: error0("boundary_nminus_coupling_float: invalid mu=%d\n", mu );
  }
}


void boundary_nplus_coupling_double( double *eta, double *D, double *phi, int mu, int start, int end, int *ind, int *neighbor ) { }
void boundary_npT_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_add_ps
#define MU 0
#define BOUNDARY
#include neighbor_coupling_file
#undef BOUNDARY
#undef MU
#undef UPD
}
void boundary_npZ_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_add_ps
#define MU 1
#define BOUNDARY
#include neighbor_coupling_file
#undef BOUNDARY
#undef MU
#undef UPD
}
void boundary_npY_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_add_ps
#define MU 2
#define BOUNDARY
#include neighbor_coupling_file
#undef BOUNDARY
#undef MU
#undef UPD
}
void boundary_npX_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_add_ps
#define MU 3
#define BOUNDARY
#include neighbor_coupling_file
#undef BOUNDARY
#undef MU
#undef UPD
}
void boundary_nplus_coupling_float( float *eta, float *D, float *phi, int mu, int start, int end, int *ind, int *neighbor ) {
  switch ( mu ) {
    case T: boundary_npT_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    case Z: boundary_npZ_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    case Y: boundary_npY_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    case X: boundary_npX_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    default: error0("boundary_nplus_coupling_float: invalid mu=%d\n", mu );
  }
}


void boundary_minus_coupling_double( double *eta, double *D, double *phi, int mu, int start, int end, int *ind, int *neighbor ) { }
void boundary_mT_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_sub_ps
#define MU 0
#define BOUNDARY
#define MINUSDIR
#include neighbor_coupling_file
#undef MINUSDIR
#undef BOUNDARY
#undef MU
#undef UPD
}
void boundary_mZ_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_sub_ps
#define MU 1
#define BOUNDARY
#define MINUSDIR
#include neighbor_coupling_file
#undef MINUSDIR
#undef BOUNDARY
#undef MU
#undef UPD
}
void boundary_mY_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_sub_ps
#define MU 2
#define BOUNDARY
#define MINUSDIR
#include neighbor_coupling_file
#undef MINUSDIR
#undef BOUNDARY
#undef MU
#undef UPD
}
void boundary_mX_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_sub_ps
#define MU 3
#define BOUNDARY
#define MINUSDIR
#include neighbor_coupling_file
#undef MINUSDIR
#undef BOUNDARY
#undef MU
#undef UPD
}
void boundary_minus_coupling_float( float *eta, float *D, float *phi, int mu, int start, int end, int *ind, int *neighbor ) {
  switch ( mu ) {
    case T: boundary_mT_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    case Z: boundary_mZ_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    case Y: boundary_mY_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    case X: boundary_mX_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    default: error0("boundary_minus_coupling_float: invalid mu=%d\n", mu );
  }
}


void boundary_plus_coupling_double( double *eta, double *D, double *phi, int mu, int start, int end, int *ind, int *neighbor ) { }
void boundary_pT_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_sub_ps
#define MU 0
#define BOUNDARY
#include neighbor_coupling_file
#undef BOUNDARY
#undef MU
#undef UPD
}
void boundary_pZ_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_sub_ps
#define MU 1
#define BOUNDARY
#include neighbor_coupling_file
#undef BOUNDARY
#undef MU
#undef UPD
}
void boundary_pY_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_sub_ps
#define MU 2
#define BOUNDARY
#include neighbor_coupling_file
#undef BOUNDARY
#undef MU
#undef UPD
}
void boundary_pX_coupling_float( float *eta, float *D, float *phi, int start, int end, int *ind, int *neighbor ) {
#define UPD _mm_sub_ps
#define MU 3
#define BOUNDARY
#include neighbor_coupling_file
#undef BOUNDARY
#undef MU
#undef UPD
}
void boundary_plus_coupling_float( float *eta, float *D, float *phi, int mu, int start, int end, int *ind, int *neighbor ) {
  switch ( mu ) {
    case T: boundary_pT_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    case Z: boundary_pZ_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    case Y: boundary_pY_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    case X: boundary_pX_coupling_float( eta, D, phi, start, end, ind, neighbor ); break;
    default: error0("boundary_plus_coupling_float: invalid mu=%d\n", mu );
  }
}





static inline int sse_clover_real_index( int i, int j ) {
  return (i/SIMD_LENGTH_float)*12*SIMD_LENGTH_float + SIMD_LENGTH_float*j*2 + i%SIMD_LENGTH_float;
}

static inline int sse_clover_imag_index( int i, int j ) {
  return (i/SIMD_LENGTH_float)*12*SIMD_LENGTH_float + SIMD_LENGTH_float*(j*2+1) + i%SIMD_LENGTH_float;
}

void sse_set_clover_double( double *out, complex_double *in ) { }

void sse_set_clover_float( float *out, complex_double *in ) {
    
  int index;
  float sign = 0.0;
  for ( int k=0; k<12; k+=SIMD_LENGTH_float ) {
    for ( int j=0; j<6; j++ ) {
      for ( int i=0; i<SIMD_LENGTH_float; i++ ) {
        if ( i+k == j || i+k-6 == j ) {
          // diagonal entry i+k,i+k
          index = i+k;
          sign = 1.0;
        } else if ( i+k<6 ) {
          // first 6-by-6 matrix
          if ( j > i+k ) {
            // upper triangle
            index = 12 + ( 30 - (5-(k+i))*(6-(k+i)) )/2 + (j-(i+k+1));
            sign = 1.0;
          } else {
            // lower triangle, j < i+k
            index = 12 + ( 30 - (5-(j))*(6-(j)) )/2 + ((i+k)-(j+1));
            sign = -1.0;
          } 
        } else {
          // i+k >= 6
          // second 6-by-6 matrix
          if ( j > i+k-6 ) {
            // upper triangle
            index = 12 + 15 + ( 30 - (5-(k+i-6))*(6-(k+i-6)) )/2 + (j-(i+k-6+1));
            sign = 1.0;
          } else {
            // j < i+k-6
            // lower triangle
            index = 12 + 15 + ( 30 - (5-(j))*(6-(j)) )/2 + ((i+k-6)-(j+1));
            sign = -1.0;
          }
        }
        out[ sse_clover_real_index(i+k,j) ] = creal_float( (complex_float)in[index] );
        out[ sse_clover_imag_index(i+k,j) ] = sign*cimag_float( (complex_float)in[index] );
      }
    }
  }
}


void sse_clover_double( vector_double eta, vector_double phi, operator_double_struct *op, int start, int end,
                        level_struct *l, struct Thread *threading ) {
  
  clover_double( eta+start, phi+start, op->clover+((start/12)*42), end-start, l, threading );
  
}


void sse_clover_float( vector_float eta, vector_float phi, operator_float_struct *op, int start, int end,
                       level_struct *l, struct Thread *threading ) {
  
  float *clov = op->clover_vectorized+start*12;
  for ( int i=start; i<end; i+=12 ) {
    sse_site_clover_float( (float*)(eta+i), (float*)(phi+i), clov+12*i );
  }
}


void sse_site_clover_double( double *eta, const double *phi, const double *clover ) { 

}

void sse_site_clover_float( float *eta, const float *phi, float *clover ) {
  
  __m128 in_re;
  __m128 in_im;
  
  __m128 clov_re;
  __m128 clov_im;
  
  __m128 out_re;
  __m128 out_im;
  
  // lines 1--4
  in_re = _mm_set1_ps( phi[0] );
  in_im = _mm_set1_ps( phi[1] );
  clov_re = _mm_load_ps( clover );
  clov_im = _mm_load_ps( clover+SIMD_LENGTH_float );
  cmul( clov_re, clov_im, in_re, in_im, &out_re, &out_im );
  clover+=2*SIMD_LENGTH_float;
  
  for ( int i=1; i<6; i++ ) {
    in_re = _mm_set1_ps( phi[2*i] );
    in_im = _mm_set1_ps( phi[2*i+1] );
    clov_re = _mm_load_ps( clover );
    clov_im = _mm_load_ps( clover+SIMD_LENGTH_float );
    cfmadd( clov_re, clov_im, in_re, in_im, &out_re, &out_im );
    clover+=2*SIMD_LENGTH_float;
  }
  
  sse_complex_interleaved_store( out_re, out_im, eta );
  
  // lines 5--8 
  in_re = _mm_setr_ps( phi[0], phi[0], phi[12], phi[12] );
  in_im = _mm_setr_ps( phi[1], phi[1], phi[13], phi[13] );
  clov_re = _mm_load_ps( clover );
  clov_im = _mm_load_ps( clover+SIMD_LENGTH_float );
  cmul( clov_re, clov_im, in_re, in_im, &out_re, &out_im );
  clover+=2*SIMD_LENGTH_float;
  
  for ( int i=1; i<6; i++ ) {
    in_re = _mm_setr_ps( phi[2*i], phi[2*i], phi[2*i+12], phi[2*i+12] );
    in_im = _mm_setr_ps( phi[2*i+1], phi[2*i+1], phi[2*i+13], phi[2*i+13] );
    clov_re = _mm_load_ps( clover );
    clov_im = _mm_load_ps( clover+SIMD_LENGTH_float );
    cfmadd( clov_re, clov_im, in_re, in_im, &out_re, &out_im );
    clover+=2*SIMD_LENGTH_float;
  }
  
  sse_complex_interleaved_store( out_re, out_im, eta+8 );
  
  // lines 9--12
  in_re = _mm_set1_ps( phi[12] );
  in_im = _mm_set1_ps( phi[13] );
  clov_re = _mm_load_ps( clover );
  clov_im = _mm_load_ps( clover+SIMD_LENGTH_float );
  cmul( clov_re, clov_im, in_re, in_im, &out_re, &out_im );
  clover+=2*SIMD_LENGTH_float;
  
  for ( int i=1; i<6; i++ ) {
    in_re = _mm_set1_ps( phi[2*i+12] );
    in_im = _mm_set1_ps( phi[2*i+13] );
    clov_re = _mm_load_ps( clover );
    clov_im = _mm_load_ps( clover+SIMD_LENGTH_float );
    cfmadd( clov_re, clov_im, in_re, in_im, &out_re, &out_im );
    clover+=2*SIMD_LENGTH_float;
  }
  
  sse_complex_interleaved_store( out_re, out_im, eta+16 );
}


void sse_site_clover_invert_double( double *clover_in, double *clover_out ) { }

void sse_site_clover_invert_float( float *clover_in, float *clover_out ) {
  
  float M_tmp1[72], M_tmp2[72];
  
  for ( int k=0; k<12; k+=SIMD_LENGTH_float ) {
    for ( int j=0; j<6; j++ ) {
      for ( int i=k; i<k+SIMD_LENGTH_float; i++ ) {
        if ( i<6 ) {
          M_tmp1[12*j+i] = *clover_in;
          M_tmp1[12*j+i+6] = *(clover_in+SIMD_LENGTH_float);
        } else {
          M_tmp2[12*j+i-6] = *clover_in;
          M_tmp2[12*j+i] = *(clover_in+SIMD_LENGTH_float);
        }
        clover_in++;
      }
      clover_in += SIMD_LENGTH_float;
    }
  }
  
  sse_cgem_inverse( 6, M_tmp1, M_tmp1, 6 );  
  sse_cgem_inverse( 6, M_tmp2, M_tmp2, 6 );
  
  for ( int k=0; k<12; k+=SIMD_LENGTH_float ) {
    for ( int j=0; j<6; j++ ) {
      for ( int i=k; i<k+SIMD_LENGTH_float; i++ ) {
        if ( i<6 ) {
          *clover_out = M_tmp1[12*j+i];
          *(clover_out+SIMD_LENGTH_float) = M_tmp1[12*j+i+6];
        } else {
          *clover_out = M_tmp2[12*j+i-6];
          *(clover_out+SIMD_LENGTH_float) = M_tmp2[12*j+i];
        }
        clover_out++;
      }
      clover_out += SIMD_LENGTH_float;
    }
  }
}


#endif
