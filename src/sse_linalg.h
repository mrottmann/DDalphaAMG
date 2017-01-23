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

#ifndef LINALG_SSE_H
#define LINALG_SSE_H
#ifdef SSE


// Standard Gram-Schmidt on aggregates
static inline void sse_aggregate_gram_schmidt_float( complex_float *V, const int num_vec,
    level_struct *l, struct Thread *threading );
// Gram-Schmidt on a block of vectors, used by Block-Gram-Schmidt
static inline void sse_aggregate_gram_schmidt_block_float( float *V,
    int num_vec, int leading_dimension, level_struct *l, struct Thread *threading );
// used by Block-Gram-Schmidt
static inline void sse_aggregate_block_dot_block_float( float *S, float *U, float *B,
    int num_vec, int leading_dimension, level_struct *l, struct Thread *threading );
// used by Block-Gram-Schmidt
static inline void sse_aggregate_block_minus_block_times_dot_float( float *B, float *U, float *S,
    int num_vec, int leading_dimension, level_struct *l, struct Thread *threading );

static inline void sse_aggregate_gram_schmidt_double( complex_double *V, const int num_vec,
    level_struct *l, struct Thread *threading ) {}
static inline void sse_aggregate_gram_schmidt_block_double( double *V,
    int num_vec, int leading_dimension, level_struct *l, struct Thread *threading ) {}
static inline void sse_aggregate_block_dot_block_double( double *S, double *U, double *B,
    int num_vec, int leading_dimension, level_struct *l, struct Thread *threading ) {}
static inline void sse_aggregate_block_minus_block_times_dot_double( double *B, double *U, double *S,
    int num_vec, int leading_dimension, level_struct *l, struct Thread *threading ) {}


static inline void sse_aggregate_gram_schmidt_float( complex_float *V, const int num_vec, level_struct *l, struct Thread *threading ) {

  PROF_float_START( _GRAM_SCHMIDT_ON_AGGREGATES, threading );
  SYNC_CORES(threading)
  SYNC_HYPERTHREADS(threading)
  long int i, j, k, k1, k2, k3, num_aggregates = l->s_float.num_aggregates,
      aggregate_size = l->inner_vector_size / num_aggregates, offset = l->num_lattice_site_var/2;
      
  float *v_pt1;
  float *v_pt2;
  float norm1, norm2;
  float next_norm1;
  float next_norm2;
  int ldv = SIMD_LENGTH_float;
  int V_block_offset = 2*l->vector_size; 
  
  for ( j=threading->n_thread*threading->core+threading->thread; j<num_aggregates; j+=threading->n_thread*threading->n_core ) {

    v_pt1 = (float *)V + 0 + j*aggregate_size*2*ldv;

    next_norm1 = 0.0;
    next_norm2 = 0.0;
    for ( i=0; i<aggregate_size; ) {
      for ( k=0; k<offset; k++, i++ ) {
        float *tmp = v_pt1 + i*2*ldv;
        next_norm1 += tmp[0]*tmp[0] + tmp[ldv]*tmp[ldv];
      }
      for ( k=0; k<offset; k++, i++ ) {
        float *tmp = v_pt1 + i*2*ldv;
        next_norm2 += tmp[0]*tmp[0] + tmp[ldv]*tmp[ldv];
      }
    }
    for ( k1=0; k1<num_vec; k1++ ) {
      v_pt1 = (float *)V + (k1/ldv)*V_block_offset*ldv + k1%ldv + j*aggregate_size*2*ldv;
      v_pt2 = (float *)V + j*aggregate_size*2*ldv;

      norm1 = 1.0/sqrt(next_norm1);
      norm2 = 1.0/sqrt(next_norm2);
      next_norm1 = 0.0;
      next_norm2 = 0.0;

      __m128 alpha1_re[OPERATOR_COMPONENT_OFFSET_float/SIMD_LENGTH_float];
      __m128 alpha1_im[OPERATOR_COMPONENT_OFFSET_float/SIMD_LENGTH_float];
      __m128 alpha2_re[OPERATOR_COMPONENT_OFFSET_float/SIMD_LENGTH_float];
      __m128 alpha2_im[OPERATOR_COMPONENT_OFFSET_float/SIMD_LENGTH_float];
      __m128 v1_re;
      __m128 v1_im;
      __m128 v2_re;
      __m128 v2_im;

      for ( k2=0; k2<num_vec; k2+=SIMD_LENGTH_float ) {
        alpha1_re[k2/SIMD_LENGTH_float] = _mm_setzero_ps();
        alpha1_im[k2/SIMD_LENGTH_float] = _mm_setzero_ps();
        alpha2_re[k2/SIMD_LENGTH_float] = _mm_setzero_ps();
        alpha2_im[k2/SIMD_LENGTH_float] = _mm_setzero_ps();
      }
      for ( i=0; i<aggregate_size; ) {
        // normalize v1 by scaling with previously computed factor
        // this is fused into this dotp loop, to avoid loading everything twice
        for ( k=0; k<offset; k++, i++ ) {
          float *tmp = v_pt1 + i*2*ldv;
          tmp[0]   *= norm1;
          tmp[ldv] *= norm1;
        }
        for ( k=0; k<offset; k++, i++ ) {
          float *tmp = v_pt1 + i*2*ldv;
          tmp[0]   *= norm2;
          tmp[ldv] *= norm2;
        }
        i -= 2*offset;
        // done normalizing

        for ( k=0; k<offset; k++, i++ ) {
          v1_re = _mm_set1_ps(v_pt1[(2*i+0)*ldv]);
          v1_im = _mm_set1_ps(v_pt1[(2*i+1)*ldv]);
          for ( k2=0; k2<OPERATOR_COMPONENT_OFFSET_float; k2+=SIMD_LENGTH_float ) {
            v2_re = _mm_load_ps(v_pt2 + (2*i+0)*ldv + k2*V_block_offset);
            v2_im = _mm_load_ps(v_pt2 + (2*i+1)*ldv + k2*V_block_offset);
            cfmadd_conj(v1_re, v1_im, v2_re, v2_im, &alpha1_re[k2/SIMD_LENGTH_float], &alpha1_im[k2/SIMD_LENGTH_float]);
          }
        }
        for ( k=0; k<offset; k++, i++ ) {
          v1_re = _mm_set1_ps(v_pt1[(2*i+0)*ldv]);
          v1_im = _mm_set1_ps(v_pt1[(2*i+1)*ldv]);
          for ( k2=0; k2<OPERATOR_COMPONENT_OFFSET_float; k2+=SIMD_LENGTH_float ) {
            v2_re = _mm_load_ps(v_pt2 + (2*i+0)*ldv + k2*V_block_offset);
            v2_im = _mm_load_ps(v_pt2 + (2*i+1)*ldv + k2*V_block_offset);
            cfmadd_conj(v1_re, v1_im, v2_re, v2_im, &alpha2_re[k2/SIMD_LENGTH_float], &alpha2_im[k2/SIMD_LENGTH_float]);
          }
        }
      }

      for ( i=0; i<aggregate_size; ) {
        for ( k=0; k<offset; k++, i++ ) {
          
          v1_re = _mm_set1_ps(v_pt1[(2*i+0)*ldv]);
          v1_im = _mm_set1_ps(v_pt1[(2*i+1)*ldv]);

          
          for ( k2=(k1/SIMD_LENGTH_float)*SIMD_LENGTH_float; k2<OPERATOR_COMPONENT_OFFSET_float; k2+=SIMD_LENGTH_float ) {

            v2_re = _mm_load_ps(v_pt2 + (2*i+0)*ldv + k2*V_block_offset);
            v2_im = _mm_load_ps(v_pt2 + (2*i+1)*ldv + k2*V_block_offset);
            
            if(k2 < k1+1) {
              float mask[SIMD_LENGTH_float] __attribute__((aligned(sizeof(float)*SIMD_LENGTH_float)));
              memset( mask, 255, sizeof(float)*SIMD_LENGTH_float );
              
              // emulate storing mask
              for ( k3=k2; k3<MIN(k1+1,k2+SIMD_LENGTH_float); k3++ ) {
                memset( mask+k3-k2, 0, sizeof(float) );
              }
              
              __m128 maskreg = _mm_load_ps(mask);
              
              masked_cfnmadd(alpha1_re[k2/SIMD_LENGTH_float], alpha1_im[k2/SIMD_LENGTH_float], v1_re, v1_im, &v2_re, &v2_im, maskreg );
            } else {
              cfnmadd(alpha1_re[k2/SIMD_LENGTH_float], alpha1_im[k2/SIMD_LENGTH_float], v1_re, v1_im, &v2_re, &v2_im);
            }
             
            _mm_store_ps(v_pt2 + (2*i+0)*ldv + k2*V_block_offset, v2_re);
            _mm_store_ps(v_pt2 + (2*i+1)*ldv + k2*V_block_offset, v2_im);
            
          }
        }
        for ( k=0; k<offset; k++, i++ ) {
          
          v1_re = _mm_set1_ps(v_pt1[(2*i+0)*ldv]);
          v1_im = _mm_set1_ps(v_pt1[(2*i+1)*ldv]);
          
          for ( k2=(k1/SIMD_LENGTH_float)*SIMD_LENGTH_float; k2<OPERATOR_COMPONENT_OFFSET_float; k2+=SIMD_LENGTH_float ) {
            
            v2_re = _mm_load_ps(v_pt2 + (2*i+0)*ldv + k2*V_block_offset);
            v2_im = _mm_load_ps(v_pt2 + (2*i+1)*ldv + k2*V_block_offset);
            
            if(k2 < k1+1) {
              float mask[SIMD_LENGTH_float] __attribute__((aligned(sizeof(float)*SIMD_LENGTH_float)));
              memset( mask, 255, sizeof(float)*SIMD_LENGTH_float );
          
              // emulate storing mask
              for ( k3=k2; k3<MIN(k1+1,k2+SIMD_LENGTH_float); k3++ ) {
                memset( mask+k3-k2, 0, sizeof(float) );
              }         
              
              __m128 maskreg = _mm_load_ps(mask);
              
              masked_cfnmadd(alpha2_re[k2/SIMD_LENGTH_float], alpha2_im[k2/SIMD_LENGTH_float], v1_re, v1_im, &v2_re, &v2_im, maskreg );
            } else {
              cfnmadd(alpha2_re[k2/SIMD_LENGTH_float], alpha2_im[k2/SIMD_LENGTH_float], v1_re, v1_im, &v2_re, &v2_im);
            }
             
            _mm_store_ps(v_pt2 + (2*i+0)*ldv + k2*V_block_offset, v2_re);
            _mm_store_ps(v_pt2 + (2*i+1)*ldv + k2*V_block_offset, v2_im);
            
          }
        }
        // compute norm of v_{k1+1}
        // this is fused into this axpy loop, to avoid loading everything twice
        if ( k1+1<num_vec ) {
          float *v_pt = (float *)V + ((k1+1)/ldv)*V_block_offset*ldv + (k1+1)%ldv + j*aggregate_size*2*ldv;
          i -= 2*offset;
          for ( k=0; k<offset; k++, i++ ) {
            float *tmp = v_pt + i*2*ldv;
            next_norm1 += tmp[0]*tmp[0] + tmp[ldv]*tmp[ldv];
          }
          for ( k=0; k<offset; k++, i++ ) {
            float *tmp = v_pt + i*2*ldv;
            next_norm2 += tmp[0]*tmp[0] + tmp[ldv]*tmp[ldv];
          }
        }
        // end compute norm
      }
    }
  }
  
  SYNC_HYPERTHREADS(threading)
  SYNC_CORES(threading)
  PROF_float_STOP( _GRAM_SCHMIDT_ON_AGGREGATES, 1, threading );
}

static inline void sse_aggregate_gram_schmidt_block_float( float *V,
    int num_vec, int leading_dimension, level_struct *l, struct Thread *threading ) {

  START_NO_HYPERTHREADS(threading)
  SYNC_CORES(threading)
  int i, j, k, k1, k2, num_aggregates = l->s_float.num_aggregates,
      aggregate_size = l->inner_vector_size / num_aggregates, offset = l->num_lattice_site_var/2;

  float *v_pt1;
  float *v_pt2;
  float norm;
  float next_norm;
  int ldv = leading_dimension;
  //offset = 6;


  // current thread chooses an aggregate
  for ( int jp=threading->core; jp<2*num_aggregates; jp+=threading->n_core ) {
    j = jp/2;
    int component = jp%2;


    v_pt1 = V + 2*component*offset*ldv + j*aggregate_size*2*ldv;

    next_norm = 0.0;

    // for the whole aggregate
    for ( i=0; i<aggregate_size; ) {

      // for either the first or the second half of variables
      // (depending on the value of "component")
      for ( k=0; k<offset; k++, i++ ) {
        // data layout contains ldv real parts
        // and thereafter ldv imag parts
        float *tmp = v_pt1 + i*2*ldv;
        // adds square of real part and square of imaginary part to current norm
        next_norm += tmp[0]*tmp[0] + tmp[ldv]*tmp[ldv];
      }
      // skip the other half of variables
      i += offset;
    } // i loop

    // for all test vectors
    for ( k1=0; k1<num_vec; k1++ ) {
      // v_pt1 =  [ component*offset*number of test vectors 
      //          + current test vector index
      //          + current aggregate block of test vectors ] * complex
      //
      // means: current vector
      v_pt1 = V + 2*component*offset*ldv + k1 + j*aggregate_size*2*ldv;
      // v_pt2 =  [ component*offset*number of test vectors
      //          + current aggregate block of test vectors ] * complex
      //
      // means: first vector
      v_pt2 = V + 2*component*offset*ldv + j*aggregate_size*2*ldv;

      norm = 1.0/sqrt(next_norm);
      next_norm = 0.0;

      __m128 alpha_re;
      __m128 alpha_im;
      __m128 v1_re;
      __m128 v1_im;
      __m128 v2_re;
      __m128 v2_im;

      alpha_re = _mm_setzero_ps();
      alpha_im = _mm_setzero_ps();
      for ( i=0; i<aggregate_size; ) {
        // normalize v1 by scaling with previously computed factor
        // this is fused into this dotp loop, to avoid loading everything twice
        {
          for ( k=0; k<offset; k++, i++ ) {
            float *tmp = v_pt1 + i*2*ldv;
            tmp[0]   *= norm;
            tmp[ldv] *= norm;
          }
          i += offset;
          i -= 2*offset;
        }
        // done normalizing current vector

        // calculate inner product of v_pt1 and v_pt2
        for ( k=0; k<offset; k++, i++ ) {
          v1_re = _mm_set1_ps(v_pt1[(2*i+0)*ldv]);
          v1_im = _mm_set1_ps(v_pt1[(2*i+1)*ldv]);
          v2_re = _mm_load_ps(v_pt2 + (2*i+0)*ldv);
          v2_im = _mm_load_ps(v_pt2 + (2*i+1)*ldv);
          cfmadd_conj(v1_re, v1_im, v2_re, v2_im, &alpha_re, &alpha_im);
        }
        i += offset;
      } // i loop

      if(k1 == num_vec-1)
        break; // break k1 loop

      for ( i=0; i<aggregate_size; ) {
        for ( k=0; k<offset; k++, i++ ) {
          v1_re = _mm_set1_ps(v_pt1[(2*i+0)*ldv]);
          v1_im = _mm_set1_ps(v_pt1[(2*i+1)*ldv]);
          
          float buffer[SIMD_LENGTH_float] __attribute__((aligned(sizeof(float)*SIMD_LENGTH_float)));
          
          _mm_store_ps( buffer, alpha_re );
          for ( k2=0; k2<MIN(k1+1,SIMD_LENGTH_float); k2++ )
            buffer[k2]=0.0;
          alpha_re = _mm_load_ps(buffer);
          
          _mm_store_ps( buffer, alpha_im );
          for ( k2=0; k2<MIN(k1+1,SIMD_LENGTH_float); k2++ )
            buffer[k2]=0.0;
          alpha_im = _mm_load_ps(buffer);
          
          v2_re = _mm_load_ps(v_pt2 + (2*i+0)*ldv);
          v2_im = _mm_load_ps(v_pt2 + (2*i+1)*ldv);
          cfnmadd(alpha_re, alpha_im, v1_re, v1_im, &v2_re, &v2_im);
          _mm_store_ps(v_pt2 + (2*i+0)*ldv, v2_re);
          _mm_store_ps(v_pt2 + (2*i+1)*ldv, v2_im);
        }
        i += offset;
        // compute norm of v_{k1+1}
        // this is fused into this axpy loop, to avoid loading everything twice
        {
          float *v_pt = V + 2*component*offset*ldv + k1+1 + j*aggregate_size*2*ldv;
          i -= 2*offset;
          for ( k=0; k<offset; k++, i++ ) {
            float *tmp = v_pt + i*2*ldv;
            next_norm += tmp[0]*tmp[0] + tmp[ldv]*tmp[ldv];
          }
          i += offset;
        }
          // end compute norm
      } // i loop
    } // k1 loop
  } // j loop
  SYNC_CORES(threading)
  END_NO_HYPERTHREADS(threading)
}

static inline void sse_aggregate_block_dot_block_float( float *S, float *U, float *B,
    int num_vec, int leading_dimension, level_struct *l, struct Thread *threading ) {

  START_NO_HYPERTHREADS(threading)

  // we compute S = U^dagger B
  // U has 16 columns, B has num_vec <= 16 columns
  // for vectorization it is more efficient to transpose the MM product:
  // S^T = B^T U^*

  float *Up;
  float *Bp;

  // factor 2 is for counting spin01 and spin23 aggregates separately
  int num_aggregates = 2*l->s_float.num_aggregates;
  int aggregate_size = l->inner_vector_size / num_aggregates;
  int offset = l->num_lattice_site_var/2;

  for ( int jp=threading->core; jp<num_aggregates; jp+=threading->n_core ) {
    int j = jp/2;
    int component = jp%2;
    // factors 2 are for complex and spin01/23 aggregates
    Up = U + 2*component*offset*leading_dimension + 2*2*j*aggregate_size*leading_dimension;
    Bp = B + 2*component*offset*leading_dimension + 2*2*j*aggregate_size*leading_dimension;
    __m128 U_re;
    __m128 U_im;
    __m128 B_re;
    __m128 B_im;
    __m128 S_re[SIMD_LENGTH_float];
    __m128 S_im[SIMD_LENGTH_float];
    for( int i=0; i<SIMD_LENGTH_float; i++) {
      S_re[i] = _mm_setzero_ps();
      S_im[i] = _mm_setzero_ps();
    }
    for ( int i=0; i<aggregate_size; i+=offset ) {
      for ( int k=0; k<offset; k++ ) {
        U_re = _mm_load_ps(Up);
        U_im = _mm_load_ps(Up + leading_dimension);
        for ( int vec=0; vec<num_vec; vec++ ) {
          B_re = _mm_set1_ps(Bp[vec]);
          B_im = _mm_set1_ps(Bp[vec + leading_dimension]);
          cfmadd_conj(U_re, U_im, B_re, B_im, S_re + vec, S_im + vec);
        }
        Bp += 2*leading_dimension;
        Up += 2*leading_dimension;
      }
      Bp += 2*leading_dimension*offset;
      Up += 2*leading_dimension*offset;
    }
    for( int i=0; i<SIMD_LENGTH_float; i++) {
      _mm_store_ps(S+(2*(SIMD_LENGTH_float*jp+i)+0)*SIMD_LENGTH_float, S_re[i]);
      _mm_store_ps(S+(2*(SIMD_LENGTH_float*jp+i)+1)*SIMD_LENGTH_float, S_im[i]);
    }
    // this stored S^T in row-major format == S in column major
  }

  END_NO_HYPERTHREADS(threading)
}

static inline void sse_aggregate_block_minus_block_times_dot_float( float *B, float *U, float *S,
    int num_vec, int leading_dimension, level_struct *l, struct Thread *threading ) {

  START_NO_HYPERTHREADS(threading)

  // we compute B -= U S
  // U has 16 columns, B has num_vec <= 16

  float *Up;
  float *Bp;

  // factor 2 is for counting spin01 and spin23 aggregates separately
  int num_aggregates = 2*l->s_float.num_aggregates;
  int aggregate_size = l->inner_vector_size / num_aggregates;
  int offset = l->num_lattice_site_var/2;

  for ( int jp=threading->core; jp<num_aggregates; jp+=threading->n_core ) {
    int j = jp/2;
    int component = jp%2;
    // factors 2 are for complex and spin01/23 aggregates
    Up = U + 2*component*offset*leading_dimension + 2*2*j*aggregate_size*leading_dimension;
    Bp = B + 2*component*offset*leading_dimension + 2*2*j*aggregate_size*leading_dimension;
    __m128 U_re;
    __m128 U_im;
    __m128 B_re;
    __m128 B_im;
    __m128 S_re[SIMD_LENGTH_float];
    __m128 S_im[SIMD_LENGTH_float];
    for( int i=0; i<SIMD_LENGTH_float; i++) {
      S_re[i] = _mm_load_ps(S+(2*(SIMD_LENGTH_float*jp+i)+0)*SIMD_LENGTH_float);
      S_im[i] = _mm_load_ps(S+(2*(SIMD_LENGTH_float*jp+i)+1)*SIMD_LENGTH_float);
    }
    for ( int i=0; i<aggregate_size; i+=offset ) {
      for ( int k=0; k<offset; k++ ) {
        U_re = _mm_load_ps(Up);
        U_im = _mm_load_ps(Up + leading_dimension);
        for ( int vec=0; vec<num_vec; vec++ ) {
          cmul(U_re, U_im, S_re[vec], S_im[vec], &B_re, &B_im);
                    
          // horizontal add and subtract from Bp
          __m128 tmp1;
          __m128 tmp2;
          
          tmp1 = _mm_add_ps( B_re, _mm_movehl_ps( B_re, B_re ) );
          B_re = _mm_add_ss( tmp1, _mm_shuffle_ps( tmp1, tmp1, 1 ) );
            
          tmp1 = _mm_add_ps( B_im, _mm_movehl_ps( B_im, B_im ) );
          B_im = _mm_add_ss( tmp1, _mm_shuffle_ps( tmp1, tmp1, 1 ) );
          
          tmp1 = _mm_set1_ps(Bp[vec]);
          tmp2 = _mm_set1_ps(Bp[vec + leading_dimension]);
          B_re = _mm_sub_ps(B_re, tmp1);
          B_im = _mm_sub_ps(B_im, tmp2);
          
          _mm_store_ss( Bp+vec, B_re );
          _mm_store_ss( Bp+vec+leading_dimension, B_im );
        }
        Bp += 2*leading_dimension;
        Up += 2*leading_dimension;
      }
      Bp += 2*leading_dimension*offset;
      Up += 2*leading_dimension*offset;
    }
  }

  END_NO_HYPERTHREADS(threading)
}

#endif // SSE
#endif // LINALG_MIC_H
