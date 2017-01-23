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

#ifndef SSE_COARSE_OPERATOR_PRECISION_HEADER
  #define SSE_COARSE_OPERATOR_PRECISION_HEADER

  #ifdef SSE
  
  #include "blas_vectorized.h"
  
  void coarse_operator_PRECISION_setup_vectorized( complex_PRECISION *operator, level_struct *l, struct Thread *threading );
  void set_coarse_self_coupling_PRECISION_vectorized( complex_PRECISION *spin_0_1, complex_PRECISION *spin_2_3,
      complex_PRECISION *V, level_struct *l, int site, const int n_rhs, complex_PRECISION *tmp );
  void set_coarse_self_coupling_PRECISION_vectorized_finalize( level_struct *l, int site, const int n_rhs, complex_PRECISION *tmp );
  // here we do not check whether site is really on boundary, caller is responsible for that
  // tmp is used to store coarse operator with padding, until sum over all sites has been done
  void set_coarse_neighbor_coupling_PRECISION_vectorized( complex_PRECISION *buffer1, complex_PRECISION *buffer2,
      complex_PRECISION *V, const int mu, level_struct *l, int site, const int n_rhs, complex_PRECISION *tmp );
  void set_coarse_neighbor_coupling_PRECISION_vectorized_finalize( const int mu, level_struct *l, int site, const int n_rhs, complex_PRECISION *tmp );

  void copy_coarse_operator_to_vectorized_layout_PRECISION(config_PRECISION D,
      OPERATOR_TYPE_PRECISION *D_vectorized, int num_aggregates, int num_eig_vect);
  // fw and bw links have a symmetry that allows constructing one from another, see, e.g., coarse_hopp_PRECISION
  // for vectorization we store the operator for both cases, the "daggered" links need this transformed layout
  void copy_coarse_operator_to_transformed_vectorized_layout_PRECISION(config_PRECISION D,
      OPERATOR_TYPE_PRECISION *D_vectorized, int num_aggregates, int num_eig_vect);
  void copy_coarse_operator_clover_to_vectorized_layout_PRECISION(config_PRECISION clover,
      OPERATOR_TYPE_PRECISION *clover_vectorized, int num_aggregates, int num_eig_vect);
  
  void coarse_spinwise_site_self_couplings_PRECISION_vectorized(
      complex_PRECISION *eta1, complex_PRECISION *eta2,
      complex_PRECISION *phi, config_PRECISION clover, int elements, level_struct *l );
  
  void coarse_aggregate_self_couplings_PRECISION_vectorized( complex_PRECISION *eta1, complex_PRECISION *eta2,
      complex_PRECISION *phi, schwarz_PRECISION_struct *s, level_struct *l,
      int site, int *direction_flags );  
  
  void coarse_aggregate_neighbor_couplings_PRECISION_vectorized( complex_PRECISION *eta1, complex_PRECISION *eta2,
      complex_PRECISION *phi, const int mu, schwarz_PRECISION_struct *s, level_struct *l,
      int site );  
  
  
  static inline void coarse_hopp_PRECISION_vectorized( vector_PRECISION eta, vector_PRECISION phi,
      OPERATOR_TYPE_PRECISION *D, level_struct *l ) {
#ifdef VECTORIZE_COARSE_OPERATOR_PRECISION
    int lda = SIMD_LENGTH_PRECISION*((l->num_lattice_site_var+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
    cgenmv(l->num_lattice_site_var, D, lda, (float *)phi, (float *)eta);
#endif
  }
  static inline void coarse_n_hopp_PRECISION_vectorized( vector_PRECISION eta, vector_PRECISION phi,
      OPERATOR_TYPE_PRECISION *D, level_struct *l ) {
#ifdef VECTORIZE_COARSE_OPERATOR_PRECISION
    int lda = SIMD_LENGTH_PRECISION*((l->num_lattice_site_var+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
    cgemv(l->num_lattice_site_var, D, lda, (float *)phi, (float *)eta);
#endif
  }

  static inline void coarse_self_couplings_PRECISION_vectorized( vector_PRECISION eta, vector_PRECISION phi, 
                                                                 OPERATOR_TYPE_PRECISION *clover, int start, int end, level_struct *l ) {
#ifdef VECTORIZE_COARSE_OPERATOR_PRECISION
    int site_size = l->num_lattice_site_var;
    int lda = SIMD_LENGTH_PRECISION*((site_size+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);

    for(int i=start; i<end; i++) {
      for(int j=0; j<site_size; j++)
        eta[i*site_size+j] = 0.0;
      cgemv(site_size, clover+i*2*site_size*lda, lda, (float *)(phi+i*site_size), (float *)(eta+i*site_size));
    }
#endif
  }
  

  static inline void coarse_spinwise_hopp_PRECISION_vectorized( complex_PRECISION *eta1, complex_PRECISION *eta2,
      complex_PRECISION *phi, config_PRECISION D, int elements, level_struct *l ) {

#ifdef SSE
    int num_eig_vect = l->num_lattice_site_var/2;
    int num_eig_vect2 = num_eig_vect*num_eig_vect;
    complex_PRECISION *eta[2] = {eta1, eta2};
    // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
    //             C D ]                        -B*  D* ]
    // storage order: A, C, B, D

    __m128 D_re;
    __m128 D_im;
    __m128 in_re;
    __m128 in_im;
    __m128 out_re;
    __m128 out_im;
    // s refers to "spin" components 0and1 (->eta1) or 2and3 (->eta2)
    for(int s=0; s<2; s++) {
      // t is the row of the input matrix (in 2x2 block form)
      for(int t=0; t<2; t++) {
        for(int i=0; i<elements; i+=SIMD_LENGTH_PRECISION) {
          for(int column=0; column<num_eig_vect; column++) {
            in_re  = _mm_load_ps((float *)phi + i + (2*column+0)*elements);
            in_im  = _mm_load_ps((float *)phi + i + (2*column+1)*elements);
            for(int row=0; row<num_eig_vect; row++) {
              out_re = _mm_load_ps((float *)eta[s] + i + (2*row+0)*elements);
              out_im = _mm_load_ps((float *)eta[s] + i + (2*row+1)*elements);
              D_re = _mm_set1_ps(creal(D[column*num_eig_vect+row]));
              D_im = _mm_set1_ps(cimag(D[column*num_eig_vect+row]));

              cfmadd(D_re, D_im, in_re, in_im, &out_re, &out_im);

              _mm_store_ps((float *)eta[s] + i + (2*row+0)*elements, out_re);
              _mm_store_ps((float *)eta[s] + i + (2*row+1)*elements, out_im);
            }
          }
        }
        eta[s] += num_eig_vect*elements;
        D += num_eig_vect2;
      }
      phi += num_eig_vect*elements;
    }
#endif
  }  
  
  
  static inline void coarse_spinwise_n_hopp_PRECISION_vectorized( complex_PRECISION *eta1, complex_PRECISION *eta2,
      complex_PRECISION *phi, config_PRECISION D, int elements, level_struct *l ) {

#ifdef SSE
    int num_eig_vect = l->num_lattice_site_var/2;
    int num_eig_vect2 = num_eig_vect*num_eig_vect;
    complex_PRECISION *eta[2] = {eta1, eta2};
    // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
    //             C D ]                        -B*  D* ]
    // storage order: A, C, B, D
    // note: minus sign of D = self_coupling - hopping_term is added here

    __m128 D_re;
    __m128 D_im;
    __m128 in_re;
    __m128 in_im;
    __m128 out_re;
    __m128 out_im;
    // s refers to "spin" components 0and1 (->eta1) or 2and3 (->eta2)
    for(int s=0; s<2; s++) {
      // t is the row of the input matrix (in 2x2 block form)
      for(int t=0; t<2; t++) {
        for(int i=0; i<elements; i+=SIMD_LENGTH_PRECISION) {
          for(int column=0; column<num_eig_vect; column++) {
            in_re  = _mm_load_ps((float *)phi + i + (2*column+0)*elements);
            in_im  = _mm_load_ps((float *)phi + i + (2*column+1)*elements);
            for(int row=0; row<num_eig_vect; row++) {
              out_re = _mm_load_ps((float *)eta[s] + i + (2*row+0)*elements);
              out_im = _mm_load_ps((float *)eta[s] + i + (2*row+1)*elements);
              D_re = _mm_set1_ps(creal(D[column*num_eig_vect+row]));
              D_im = _mm_set1_ps(cimag(D[column*num_eig_vect+row]));

              cfnmadd(D_re, D_im, in_re, in_im, &out_re, &out_im);

              _mm_store_ps((float *)eta[s] + i + (2*row+0)*elements, out_re);
              _mm_store_ps((float *)eta[s] + i + (2*row+1)*elements, out_im);
            }
          }
        }
        eta[s] += num_eig_vect*elements;
        D += num_eig_vect2;
      }
      phi += num_eig_vect*elements;
    }
#endif
  }
  

  static inline void coarse_spinwise_n_daggered_hopp_PRECISION_vectorized( complex_PRECISION *eta1, complex_PRECISION *eta2,
      complex_PRECISION *phi, config_PRECISION D, int elements, level_struct *l ) {

#ifdef SSE
    int num_eig_vect = l->num_lattice_site_var/2;
    int num_eig_vect2 = num_eig_vect*num_eig_vect;
    complex_PRECISION *eta[2] = {eta1, eta2};
    // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
    //             C D ]                        -B*  D* ]
    // storage order: A, C, B, D
    // note: minus sign of D = self_coupling - hopping_term is added here

    __m128 D_re;
    __m128 D_im;
    __m128 in_re;
    __m128 in_im;
    __m128 out_re;
    __m128 out_im;
    // A*
    for(int i=0; i<elements; i+=SIMD_LENGTH_PRECISION) {
      for(int column=0; column<num_eig_vect; column++) {
        in_re  = _mm_load_ps((float *)phi + i + (2*column+0)*elements);
        in_im  = _mm_load_ps((float *)phi + i + (2*column+1)*elements);
        for(int row=0; row<num_eig_vect; row++) {
          out_re = _mm_load_ps((float *)eta[0] + i + (2*row+0)*elements);
          out_im = _mm_load_ps((float *)eta[0] + i + (2*row+1)*elements);
          // load transpose
          D_re = _mm_set1_ps(creal(D[row*num_eig_vect+column]));
          D_im = _mm_set1_ps(cimag(D[row*num_eig_vect+column]));

          cfnmadd_conj(D_re, D_im, in_re, in_im, &out_re, &out_im);

          _mm_store_ps((float *)eta[0] + i + (2*row+0)*elements, out_re);
          _mm_store_ps((float *)eta[0] + i + (2*row+1)*elements, out_im);
        }
      }
    }
    // -C*
    phi += num_eig_vect*elements;
    D += num_eig_vect2;
    for(int i=0; i<elements; i+=SIMD_LENGTH_PRECISION) {
      for(int column=0; column<num_eig_vect; column++) {
        in_re  = _mm_load_ps((float *)phi + i + (2*column+0)*elements);
        in_im  = _mm_load_ps((float *)phi + i + (2*column+1)*elements);
        for(int row=0; row<num_eig_vect; row++) {
          out_re = _mm_load_ps((float *)eta[1] + i + (2*row+0)*elements);
          out_im = _mm_load_ps((float *)eta[1] + i + (2*row+1)*elements);
          // load transpose
          D_re = _mm_set1_ps(creal(D[row*num_eig_vect+column]));
          D_im = _mm_set1_ps(cimag(D[row*num_eig_vect+column]));

          cfmadd_conj(D_re, D_im, in_re, in_im, &out_re, &out_im);

          _mm_store_ps((float *)eta[1] + i + (2*row+0)*elements, out_re);
          _mm_store_ps((float *)eta[1] + i + (2*row+1)*elements, out_im);
        }
      }
    }
    // -B*
    eta[0] += num_eig_vect*elements;
    phi -= num_eig_vect*elements;
    D += num_eig_vect2;
    for(int i=0; i<elements; i+=SIMD_LENGTH_PRECISION) {
      for(int column=0; column<num_eig_vect; column++) {
        in_re  = _mm_load_ps((float *)phi + i + (2*column+0)*elements);
        in_im  = _mm_load_ps((float *)phi + i + (2*column+1)*elements);
        for(int row=0; row<num_eig_vect; row++) {
          out_re = _mm_load_ps((float *)eta[0] + i + (2*row+0)*elements);
          out_im = _mm_load_ps((float *)eta[0] + i + (2*row+1)*elements);
          // load transpose
          D_re = _mm_set1_ps(creal(D[row*num_eig_vect+column]));
          D_im = _mm_set1_ps(cimag(D[row*num_eig_vect+column]));

          cfmadd_conj(D_re, D_im, in_re, in_im, &out_re, &out_im);

          _mm_store_ps((float *)eta[0] + i + (2*row+0)*elements, out_re);
          _mm_store_ps((float *)eta[0] + i + (2*row+1)*elements, out_im);
        }
      }
    }
    // D*
    eta[1] += num_eig_vect*elements;
    phi += num_eig_vect*elements;
    D += num_eig_vect2;
    for(int i=0; i<elements; i+=SIMD_LENGTH_PRECISION) {
      for(int column=0; column<num_eig_vect; column++) {
        in_re  = _mm_load_ps((float *)phi + i + (2*column+0)*elements);
        in_im  = _mm_load_ps((float *)phi + i + (2*column+1)*elements);
        for(int row=0; row<num_eig_vect; row++) {
          out_re = _mm_load_ps((float *)eta[1] + i + (2*row+0)*elements);
          out_im = _mm_load_ps((float *)eta[1] + i + (2*row+1)*elements);
          // load transpose
          D_re = _mm_set1_ps(creal(D[row*num_eig_vect+column]));
          D_im = _mm_set1_ps(cimag(D[row*num_eig_vect+column]));

          cfnmadd_conj(D_re, D_im, in_re, in_im, &out_re, &out_im);

          _mm_store_ps((float *)eta[1] + i + (2*row+0)*elements, out_re);
          _mm_store_ps((float *)eta[1] + i + (2*row+1)*elements, out_im);
        }
      }
    }
#endif
  }
  
#endif
#endif