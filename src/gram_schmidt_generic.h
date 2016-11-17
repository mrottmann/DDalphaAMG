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

#ifndef GRAM_SCHMIDT_PRECISION_HEADER
#define GRAM_SCHMIDT_PRECISION_HEADER

// Gram-Schmidt on full vectors
void gram_schmidt_PRECISION( vector_PRECISION *V, complex_PRECISION *buffer, const int start, const int n,
                             level_struct *l, struct Thread *threading );
// Gram-Schmidt on aggregates
#ifndef OPTIMIZED_INTERPOLATION_SETUP_PRECISION  
void gram_schmidt_on_aggregates_PRECISION( vector_PRECISION *V, const int num_vec,
                                           level_struct *l, struct Thread *threading );
#else // optimized version on the operator layout
void gram_schmidt_on_aggregates_PRECISION( complex_PRECISION *operator, const int num_vec,
                                           level_struct *l, struct Thread *threading );
#endif

#ifdef OPTIMIZED_INTERPOLATION_SETUP_PRECISION  

// SIMD version of gram_schmidt_on_aggregates optimized on a operator layout
// block_gram_schmidt_PRECISION follows, after definition of others inline void functions marked by "used by *IT*"
static inline void block_gram_schmidt_PRECISION( complex_PRECISION *V, int num_vec, level_struct *l, 
                                                 struct Thread *threading );

static inline void aggregate_gram_schmidt_PRECISION( complex_PRECISION *V, const int num_vec, 
                                                     level_struct *l, struct Thread *threading ) {

  SYNC_CORES(threading)
  SYNC_HYPERTHREADS(threading)
  int i, j, k, k1, k2, k3, num_aggregates = l->s_PRECISION.num_aggregates,
      aggregate_size = l->inner_vector_size / num_aggregates, offset = l->num_lattice_site_var/2;
      
  PRECISION *v_pt1;
  PRECISION *v_pt2;
  PRECISION norm1, norm2;
  PRECISION next_norm1;
  PRECISION next_norm2;
  int ldv = SIMD_LENGTH_PRECISION;
  int V_block_offset = 2*l->vector_size; 
  
  for ( j=threading->n_thread*threading->core+threading->thread; j<num_aggregates; j+=threading->n_thread*threading->n_core ) {

    v_pt1 = (PRECISION *)V + 0 + j*aggregate_size*2*ldv;

    next_norm1 = 0.0;
    next_norm2 = 0.0;
    for ( i=0; i<aggregate_size; ) {
      for ( k=0; k<offset; k++, i++ ) {
        PRECISION *tmp = v_pt1 + i*2*ldv;
        next_norm1 += tmp[0]*tmp[0] + tmp[ldv]*tmp[ldv];
      }
      for ( k=0; k<offset; k++, i++ ) {
        PRECISION *tmp = v_pt1 + i*2*ldv;
        next_norm2 += tmp[0]*tmp[0] + tmp[ldv]*tmp[ldv];
      }
    }
    for ( k1=0; k1<num_vec; k1++ ) {
      v_pt1 = (PRECISION *)V + (k1/ldv)*V_block_offset*ldv + k1%ldv + j*aggregate_size*2*ldv;
      v_pt2 = (PRECISION *)V + j*aggregate_size*2*ldv;

      norm1 = 1.0/sqrt(next_norm1);
      norm2 = 1.0/sqrt(next_norm2);
      next_norm1 = 0.0;
      next_norm2 = 0.0;

      mm_PRECISION alpha1_re[OPERATOR_COMPONENT_OFFSET_PRECISION/SIMD_LENGTH_PRECISION];
      mm_PRECISION alpha1_im[OPERATOR_COMPONENT_OFFSET_PRECISION/SIMD_LENGTH_PRECISION];
      mm_PRECISION alpha2_re[OPERATOR_COMPONENT_OFFSET_PRECISION/SIMD_LENGTH_PRECISION];
      mm_PRECISION alpha2_im[OPERATOR_COMPONENT_OFFSET_PRECISION/SIMD_LENGTH_PRECISION];
      mm_PRECISION v1_re;
      mm_PRECISION v1_im;
      mm_PRECISION v2_re;
      mm_PRECISION v2_im;

      for ( k2=0; k2<num_vec; k2+=SIMD_LENGTH_PRECISION ) {
        alpha1_re[k2/SIMD_LENGTH_PRECISION] = mm_setzero_PRECISION();
        alpha1_im[k2/SIMD_LENGTH_PRECISION] = mm_setzero_PRECISION();
        alpha2_re[k2/SIMD_LENGTH_PRECISION] = mm_setzero_PRECISION();
        alpha2_im[k2/SIMD_LENGTH_PRECISION] = mm_setzero_PRECISION();
      }
      for ( i=0; i<aggregate_size; ) {
        // normalize v1 by scaling with previously computed factor
        // this is fused into this dotp loop, to avoid loading everything twice
        for ( k=0; k<offset; k++, i++ ) {
          PRECISION *tmp = v_pt1 + i*2*ldv;
          tmp[0]   *= norm1;
          tmp[ldv] *= norm1;
        }
        for ( k=0; k<offset; k++, i++ ) {
          PRECISION *tmp = v_pt1 + i*2*ldv;
          tmp[0]   *= norm2;
          tmp[ldv] *= norm2;
        }
        i -= 2*offset;
        // done normalizing

        for ( k=0; k<offset; k++, i++ ) {
          v1_re = mm_set1_PRECISION(v_pt1[(2*i+0)*ldv]);
          v1_im = mm_set1_PRECISION(v_pt1[(2*i+1)*ldv]);
          for ( k2=0; k2<OPERATOR_COMPONENT_OFFSET_PRECISION; k2+=SIMD_LENGTH_PRECISION ) {
            v2_re = mm_load_PRECISION(v_pt2 + (2*i+0)*ldv + k2*V_block_offset);
            v2_im = mm_load_PRECISION(v_pt2 + (2*i+1)*ldv + k2*V_block_offset);
            cfmadd_conj_PRECISION(v1_re, v1_im, v2_re, v2_im, &alpha1_re[k2/SIMD_LENGTH_PRECISION], &alpha1_im[k2/SIMD_LENGTH_PRECISION]);
          }
        }
        for ( k=0; k<offset; k++, i++ ) {
          v1_re = mm_set1_PRECISION(v_pt1[(2*i+0)*ldv]);
          v1_im = mm_set1_PRECISION(v_pt1[(2*i+1)*ldv]);
          for ( k2=0; k2<OPERATOR_COMPONENT_OFFSET_PRECISION; k2+=SIMD_LENGTH_PRECISION ) {
            v2_re = mm_load_PRECISION(v_pt2 + (2*i+0)*ldv + k2*V_block_offset);
            v2_im = mm_load_PRECISION(v_pt2 + (2*i+1)*ldv + k2*V_block_offset);
            cfmadd_conj_PRECISION(v1_re, v1_im, v2_re, v2_im, &alpha2_re[k2/SIMD_LENGTH_PRECISION], &alpha2_im[k2/SIMD_LENGTH_PRECISION]);
          }
        }
      }

      for ( i=0; i<aggregate_size; ) {
        for ( k=0; k<offset; k++, i++ ) {
          
          v1_re = mm_set1_PRECISION(v_pt1[(2*i+0)*ldv]);
          v1_im = mm_set1_PRECISION(v_pt1[(2*i+1)*ldv]);

          
          for ( k2=(k1/SIMD_LENGTH_PRECISION)*SIMD_LENGTH_PRECISION; k2<OPERATOR_COMPONENT_OFFSET_PRECISION; k2+=SIMD_LENGTH_PRECISION ) {

            v2_re = mm_load_PRECISION(v_pt2 + (2*i+0)*ldv + k2*V_block_offset);
            v2_im = mm_load_PRECISION(v_pt2 + (2*i+1)*ldv + k2*V_block_offset);
            
            if(k2 < k1+1) {
              PRECISION mask[SIMD_LENGTH_PRECISION] __attribute__((aligned(sizeof(PRECISION)*SIMD_LENGTH_PRECISION)));
              memset( mask, 255, sizeof(PRECISION)*SIMD_LENGTH_PRECISION );
              
              // emulate storing mask
              for ( k3=k2; k3<MIN(k1+1,k2+SIMD_LENGTH_PRECISION); k3++ ) {
                memset( mask+k3-k2, 0, sizeof(PRECISION) );
              }
              
              mm_PRECISION maskreg = mm_load_PRECISION(mask);
              
              masked_cfnmadd_PRECISION(alpha1_re[k2/SIMD_LENGTH_PRECISION], alpha1_im[k2/SIMD_LENGTH_PRECISION], v1_re, v1_im, &v2_re, &v2_im, maskreg );
            } else {
              cfnmadd_PRECISION(alpha1_re[k2/SIMD_LENGTH_PRECISION], alpha1_im[k2/SIMD_LENGTH_PRECISION], v1_re, v1_im, &v2_re, &v2_im);
            }
             
            mm_store_PRECISION(v_pt2 + (2*i+0)*ldv + k2*V_block_offset, v2_re);
            mm_store_PRECISION(v_pt2 + (2*i+1)*ldv + k2*V_block_offset, v2_im);
            
          }
        }
        for ( k=0; k<offset; k++, i++ ) {
          
          v1_re = mm_set1_PRECISION(v_pt1[(2*i+0)*ldv]);
          v1_im = mm_set1_PRECISION(v_pt1[(2*i+1)*ldv]);
          
          for ( k2=(k1/SIMD_LENGTH_PRECISION)*SIMD_LENGTH_PRECISION; k2<OPERATOR_COMPONENT_OFFSET_PRECISION; k2+=SIMD_LENGTH_PRECISION ) {
            
            v2_re = mm_load_PRECISION(v_pt2 + (2*i+0)*ldv + k2*V_block_offset);
            v2_im = mm_load_PRECISION(v_pt2 + (2*i+1)*ldv + k2*V_block_offset);
            
            if(k2 < k1+1) {
              PRECISION mask[SIMD_LENGTH_PRECISION] __attribute__((aligned(sizeof(PRECISION)*SIMD_LENGTH_PRECISION)));
              memset( mask, 255, sizeof(PRECISION)*SIMD_LENGTH_PRECISION );
          
              // emulate storing mask
              for ( k3=k2; k3<MIN(k1+1,k2+SIMD_LENGTH_PRECISION); k3++ ) {
                memset( mask+k3-k2, 0, sizeof(PRECISION) );
              }         
              
              mm_PRECISION maskreg = mm_load_PRECISION(mask);
              
              masked_cfnmadd_PRECISION(alpha2_re[k2/SIMD_LENGTH_PRECISION], alpha2_im[k2/SIMD_LENGTH_PRECISION], v1_re, v1_im, &v2_re, &v2_im, maskreg );
            } else {
              cfnmadd_PRECISION(alpha2_re[k2/SIMD_LENGTH_PRECISION], alpha2_im[k2/SIMD_LENGTH_PRECISION], v1_re, v1_im, &v2_re, &v2_im);
            }
             
            mm_store_PRECISION(v_pt2 + (2*i+0)*ldv + k2*V_block_offset, v2_re);
            mm_store_PRECISION(v_pt2 + (2*i+1)*ldv + k2*V_block_offset, v2_im);
            
          }
        }
        // compute norm of v_{k1+1}
        // this is fused into this axpy loop, to avoid loading everything twice
        if ( k1+1<num_vec ) {
          PRECISION *v_pt = (PRECISION *)V + ((k1+1)/ldv)*V_block_offset*ldv + (k1+1)%ldv + j*aggregate_size*2*ldv;
          i -= 2*offset;
          for ( k=0; k<offset; k++, i++ ) {
            PRECISION *tmp = v_pt + i*2*ldv;
            next_norm1 += tmp[0]*tmp[0] + tmp[ldv]*tmp[ldv];
          }
          for ( k=0; k<offset; k++, i++ ) {
            PRECISION *tmp = v_pt + i*2*ldv;
            next_norm2 += tmp[0]*tmp[0] + tmp[ldv]*tmp[ldv];
          }
        }
        // end compute norm
      }
    }
  }
  
  SYNC_HYPERTHREADS(threading)
  SYNC_CORES(threading)
}

// used by block_gram_schmidt_PRECISION
static inline void aggregate_gram_schmidt_block_PRECISION( PRECISION *V, int num_vec, int leading_dimension, 
                                                           level_struct *l, struct Thread *threading ) {

  START_NO_HYPERTHREADS(threading)
  SYNC_CORES(threading)
  int i, j, k, k1, k2, num_aggregates = l->s_PRECISION.num_aggregates,
      aggregate_size = l->inner_vector_size / num_aggregates, offset = l->num_lattice_site_var/2;

  PRECISION *v_pt1;
  PRECISION *v_pt2;
  PRECISION norm;
  PRECISION next_norm;
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
        PRECISION *tmp = v_pt1 + i*2*ldv;
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

      mm_PRECISION alpha_re;
      mm_PRECISION alpha_im;
      mm_PRECISION v1_re;
      mm_PRECISION v1_im;
      mm_PRECISION v2_re;
      mm_PRECISION v2_im;

      alpha_re = mm_setzero_PRECISION();
      alpha_im = mm_setzero_PRECISION();
      for ( i=0; i<aggregate_size; ) {
        // normalize v1 by scaling with previously computed factor
        // this is fused into this dotp loop, to avoid loading everything twice
        {
          for ( k=0; k<offset; k++, i++ ) {
            PRECISION *tmp = v_pt1 + i*2*ldv;
            tmp[0]   *= norm;
            tmp[ldv] *= norm;
          }
          i += offset;
          i -= 2*offset;
        }
        // done normalizing current vector

        // calculate inner product of v_pt1 and v_pt2
        for ( k=0; k<offset; k++, i++ ) {
          v1_re = mm_set1_PRECISION(v_pt1[(2*i+0)*ldv]);
          v1_im = mm_set1_PRECISION(v_pt1[(2*i+1)*ldv]);
          v2_re = mm_load_PRECISION(v_pt2 + (2*i+0)*ldv);
          v2_im = mm_load_PRECISION(v_pt2 + (2*i+1)*ldv);
          cfmadd_conj_PRECISION(v1_re, v1_im, v2_re, v2_im, &alpha_re, &alpha_im);
        }
        i += offset;
      } // i loop

      if(k1 == num_vec-1)
        break; // break k1 loop

      for ( i=0; i<aggregate_size; ) {
        for ( k=0; k<offset; k++, i++ ) {
          v1_re = mm_set1_PRECISION(v_pt1[(2*i+0)*ldv]);
          v1_im = mm_set1_PRECISION(v_pt1[(2*i+1)*ldv]);
          
          PRECISION buffer[SIMD_LENGTH_PRECISION] __attribute__((aligned(sizeof(PRECISION)*SIMD_LENGTH_PRECISION)));
          
          mm_store_PRECISION( buffer, alpha_re );
          for ( k2=0; k2<MIN(k1+1,SIMD_LENGTH_PRECISION); k2++ )
            buffer[k2]=0.0;
          alpha_re = mm_load_PRECISION(buffer);
          
          mm_store_PRECISION( buffer, alpha_im );
          for ( k2=0; k2<MIN(k1+1,SIMD_LENGTH_PRECISION); k2++ )
            buffer[k2]=0.0;
          alpha_im = mm_load_PRECISION(buffer);
          
          v2_re = mm_load_PRECISION(v_pt2 + (2*i+0)*ldv);
          v2_im = mm_load_PRECISION(v_pt2 + (2*i+1)*ldv);
          cfnmadd_PRECISION(alpha_re, alpha_im, v1_re, v1_im, &v2_re, &v2_im);
          mm_store_PRECISION(v_pt2 + (2*i+0)*ldv, v2_re);
          mm_store_PRECISION(v_pt2 + (2*i+1)*ldv, v2_im);
        }
        i += offset;
        // compute norm of v_{k1+1}
        // this is fused into this axpy loop, to avoid loading everything twice
        {
          PRECISION *v_pt = V + 2*component*offset*ldv + k1+1 + j*aggregate_size*2*ldv;
          i -= 2*offset;
          for ( k=0; k<offset; k++, i++ ) {
            PRECISION *tmp = v_pt + i*2*ldv;
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

// used by block_gram_schmidt_PRECISION
static inline void aggregate_block_dot_block_PRECISION( PRECISION *S, PRECISION *U, PRECISION *B,
    int num_vec, int leading_dimension, level_struct *l, struct Thread *threading ) {

  START_NO_HYPERTHREADS(threading)

  // we compute S = U^dagger B
  // U has 16 columns, B has num_vec <= 16 columns
  // for vectorization it is more efficient to transpose the MM product:
  // S^T = B^T U^*

  PRECISION *Up;
  PRECISION *Bp;

  // factor 2 is for counting spin01 and spin23 aggregates separately
  int num_aggregates = 2*l->s_PRECISION.num_aggregates;
  int aggregate_size = l->inner_vector_size / num_aggregates;
  int offset = l->num_lattice_site_var/2;

  for ( int jp=threading->core; jp<num_aggregates; jp+=threading->n_core ) {
    int j = jp/2;
    int component = jp%2;
    // factors 2 are for complex and spin01/23 aggregates
    Up = U + 2*component*offset*leading_dimension + 2*2*j*aggregate_size*leading_dimension;
    Bp = B + 2*component*offset*leading_dimension + 2*2*j*aggregate_size*leading_dimension;
    mm_PRECISION U_re;
    mm_PRECISION U_im;
    mm_PRECISION B_re;
    mm_PRECISION B_im;
    mm_PRECISION S_re[SIMD_LENGTH_PRECISION];
    mm_PRECISION S_im[SIMD_LENGTH_PRECISION];
    for( int i=0; i<SIMD_LENGTH_PRECISION; i++) {
      S_re[i] = mm_setzero_PRECISION();
      S_im[i] = mm_setzero_PRECISION();
    }
    for ( int i=0; i<aggregate_size; i+=offset ) {
      for ( int k=0; k<offset; k++ ) {
        U_re = mm_load_PRECISION(Up);
        U_im = mm_load_PRECISION(Up + leading_dimension);
        for ( int vec=0; vec<num_vec; vec++ ) {
          B_re = mm_set1_PRECISION(Bp[vec]);
          B_im = mm_set1_PRECISION(Bp[vec + leading_dimension]);
          cfmadd_conj_PRECISION(U_re, U_im, B_re, B_im, S_re + vec, S_im + vec);
        }
        Bp += 2*leading_dimension;
        Up += 2*leading_dimension;
      }
      Bp += 2*leading_dimension*offset;
      Up += 2*leading_dimension*offset;
    }
    for( int i=0; i<SIMD_LENGTH_PRECISION; i++) {
      mm_store_PRECISION(S+(2*(SIMD_LENGTH_PRECISION*jp+i)+0)*SIMD_LENGTH_PRECISION, S_re[i]);
      mm_store_PRECISION(S+(2*(SIMD_LENGTH_PRECISION*jp+i)+1)*SIMD_LENGTH_PRECISION, S_im[i]);
    }
    // this stored S^T in row-major format == S in column major
  }

  END_NO_HYPERTHREADS(threading)
}

// used by block_gram_schmidt_PRECISION
static inline void aggregate_block_minus_block_times_dot_PRECISION( PRECISION *B, PRECISION *U, PRECISION *S,
    int num_vec, int leading_dimension, level_struct *l, struct Thread *threading ) {

  START_NO_HYPERTHREADS(threading)

  // we compute B -= U S
  // U has 16 columns, B has num_vec <= 16

  PRECISION *Up;
  PRECISION *Bp;

  // factor 2 is for counting spin01 and spin23 aggregates separately
  int num_aggregates = 2*l->s_PRECISION.num_aggregates;
  int aggregate_size = l->inner_vector_size / num_aggregates;
  int offset = l->num_lattice_site_var/2;

  for ( int jp=threading->core; jp<num_aggregates; jp+=threading->n_core ) {
    int j = jp/2;
    int component = jp%2;
    // factors 2 are for complex and spin01/23 aggregates
    Up = U + 2*component*offset*leading_dimension + 2*2*j*aggregate_size*leading_dimension;
    Bp = B + 2*component*offset*leading_dimension + 2*2*j*aggregate_size*leading_dimension;
    mm_PRECISION U_re;
    mm_PRECISION U_im;
    mm_PRECISION B_re;
    mm_PRECISION B_im;
    mm_PRECISION S_re[SIMD_LENGTH_PRECISION];
    mm_PRECISION S_im[SIMD_LENGTH_PRECISION];
    for( int i=0; i<SIMD_LENGTH_PRECISION; i++) {
      S_re[i] = mm_load_PRECISION(S+(2*(SIMD_LENGTH_PRECISION*jp+i)+0)*SIMD_LENGTH_PRECISION);
      S_im[i] = mm_load_PRECISION(S+(2*(SIMD_LENGTH_PRECISION*jp+i)+1)*SIMD_LENGTH_PRECISION);
    }
    for ( int i=0; i<aggregate_size; i+=offset ) {
      for ( int k=0; k<offset; k++ ) {
        U_re = mm_load_PRECISION(Up);
        U_im = mm_load_PRECISION(Up + leading_dimension);
        for ( int vec=0; vec<num_vec; vec++ ) {
          cmul(U_re, U_im, S_re[vec], S_im[vec], &B_re, &B_im);
                    
          // horizontal add and subtract from Bp
          PRECISION B_re_hsum = mm_reduce_add_PRECISION(B_re);
          PRECISION B_im_hsum = mm_reduce_add_PRECISION(B_im);

          Bp[vec] = B_re_hsum - Bp[vec];
          Bp[vec + leading_dimension] = B_im_hsum - Bp[vec + leading_dimension];
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

// used by block_gram_schmidt_PRECISION
static inline void aggregate_orthogonalize_block_wrt_orthonormal_block_PRECISION( PRECISION *B, PRECISION *U, int num_vec, level_struct *l, struct Thread *threading ) {
  START_NO_HYPERTHREADS(threading)

  PRECISION *S = NULL;
  START_LOCKED_MASTER(threading)
  // factors 2 are for complex and spin01/23 aggregates
  MALLOC_HUGEPAGES(S, PRECISION, 2*2*l->s_PRECISION.num_aggregates*SIMD_LENGTH_PRECISION*SIMD_LENGTH_PRECISION, 64);
  ((PRECISION **)threading->workspace)[0] = S;
  END_LOCKED_MASTER(threading)
  S = ((PRECISION **)threading->workspace)[0];

  aggregate_block_dot_block_PRECISION(S, U, B, num_vec, SIMD_LENGTH_PRECISION, l , threading);
  aggregate_block_minus_block_times_dot_PRECISION(B, U, S, num_vec, SIMD_LENGTH_PRECISION, l , threading);

  START_LOCKED_MASTER(threading)
  FREE_HUGEPAGES(S, PRECISION, 2*2*l->s_PRECISION.num_aggregates*SIMD_LENGTH_PRECISION*SIMD_LENGTH_PRECISION);
  END_LOCKED_MASTER(threading)

  END_NO_HYPERTHREADS(threading)
}

static inline void block_gram_schmidt_PRECISION( complex_PRECISION *V, int num_vec, level_struct *l, 
                                                 struct Thread *threading ) {
  SYNC_CORES(threading);
  for ( int i=0; i<num_vec; i+=SIMD_LENGTH_PRECISION ) {
    int vecs = SIMD_LENGTH_PRECISION;
    if(num_vec-i < SIMD_LENGTH_PRECISION)
      vecs = num_vec-i;
    
    for ( int j=0; j<i; j+=SIMD_LENGTH_PRECISION )
      aggregate_orthogonalize_block_wrt_orthonormal_block_PRECISION( (PRECISION *)(V + i*l->vector_size),
                                                                     (PRECISION *)(V + j*l->vector_size), vecs,
                                                                     l, threading );
    aggregate_gram_schmidt_block_PRECISION( (PRECISION *)(V + i*l->vector_size), vecs, SIMD_LENGTH_PRECISION, l, threading );
  }
  SYNC_CORES(threading);
}

#endif //OPTIMIZED_INTERPOLATION_SETUP_PRECISION
#endif
