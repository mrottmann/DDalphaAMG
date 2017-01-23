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

#ifndef COARSE_OPERATOR_SSE_H
#define COARSE_OPERATOR_SSE_H
#ifdef SSE

static inline void sse_set_coarse_self_coupling_float( complex_float *spin_0_1, complex_float *spin_2_3,
    complex_float *V, level_struct *l, int site, const int n_rhs, complex_float *tmp );
static inline void sse_set_coarse_neighbor_coupling_float( complex_float *spin_0_1, complex_float *spin_2_3,
    complex_float *V, const int mu, level_struct *l, int site, const int n_rhs, complex_float *tmp );
static inline void sse_coarse_spinwise_site_self_couplings_float( complex_float *eta1, complex_float *eta2,
    complex_float *phi, config_float clover, int elements, level_struct *l );

// not implemented for double precision
static inline void sse_set_coarse_self_coupling_double( complex_double *spin_0_1, complex_double *spin_2_3,
    complex_double *V, level_struct *l, int site, const int n_rhs, complex_double *tmp ) {}
static inline void sse_set_coarse_neighbor_coupling_double( complex_double *spin_0_1, complex_double *spin_2_3,
    complex_double *V, const int mu, level_struct *l, int site, const int n_rhs, complex_double *tmp ) {}
static inline void sse_coarse_spinwise_site_self_couplings_double( complex_double *eta1, complex_double *eta2,
    complex_double *phi, config_double clover, int elements, level_struct *l ) {}


static inline void sse_set_coarse_self_coupling_float( complex_float *spin_0_1, complex_float *spin_2_3,
    complex_float *V, level_struct *l, int site, const int n_rhs, complex_float *tmp ) {

#ifdef SSE
  int k, m, k1, k2, num_eig_vect = l->next_level->num_lattice_site_var/2,
      offset = l->num_lattice_site_var/2;
  float *spin_0_1_pt;
  float *spin_2_3_pt;
  float *interpolation_data;

  int component_offset = SIMD_LENGTH_float;
  int fine_components = l->num_lattice_site_var;

  // U(x) = [ A B      , A=A*, D=D*, C = -B*
  //          C D ]
  // storage order: upper triangle of A, upper triangle of D, B, columnwise
  // diagonal coupling
  for ( int n=0; n<n_rhs; n++ ) {

    int max = SIMD_LENGTH_float*((n+1+SIMD_LENGTH_float-1)/SIMD_LENGTH_float);
    spin_0_1_pt = (float *)(spin_0_1 + (n/SIMD_LENGTH_float)*2*offset*SIMD_LENGTH_float) + n%SIMD_LENGTH_float;
    spin_2_3_pt = (float *)(spin_2_3 + (n/SIMD_LENGTH_float)*2*offset*SIMD_LENGTH_float) + n%SIMD_LENGTH_float;


    // index k used for vectorization
    // original loop runs to k<=n, we must pad as usual to fill SIMD
    for ( k=0; k<max; k+=SIMD_LENGTH_float ) {
      __m128 buffer_re;
      __m128 buffer_im;

      // this are the packed indices, which we do not use in tmp
      //k1 = (n*(n+1))/2;
      //k2 = (n*(n+1))/2+(num_eig_vect*(num_eig_vect+1))/2;
      k1 = (n+0*num_eig_vect)*OPERATOR_COMPONENT_OFFSET_float;
      k2 = (n+1*num_eig_vect)*OPERATOR_COMPONENT_OFFSET_float;

      interpolation_data = (float *)(V + k*l->vector_size + fine_components*SIMD_LENGTH_float*site);

      // A
      buffer_re = _mm_load_ps((float *)(tmp+k1)+k+0*OPERATOR_COMPONENT_OFFSET_float);
      buffer_im = _mm_load_ps((float *)(tmp+k1)+k+1*OPERATOR_COMPONENT_OFFSET_float);
      for ( m=0; m<offset; m++ ) {
        // spin_0_1 is the same for all k => broadcast
        __m128 spin_0_1_re = _mm_set1_ps(spin_0_1_pt[(2*m+0)*component_offset]);
        __m128 spin_0_1_im = _mm_set1_ps(spin_0_1_pt[(2*m+1)*component_offset]);
        __m128 interpolation_data_re = _mm_load_ps(interpolation_data + (2*m+0)*component_offset);
        __m128 interpolation_data_im = _mm_load_ps(interpolation_data + (2*m+1)*component_offset);

        cfmadd_conj(interpolation_data_re, interpolation_data_im, spin_0_1_re, spin_0_1_im, &buffer_re, &buffer_im);
      }
      _mm_store_ps((float *)(tmp+k1)+k+0*OPERATOR_COMPONENT_OFFSET_float, buffer_re);
      _mm_store_ps((float *)(tmp+k1)+k+1*OPERATOR_COMPONENT_OFFSET_float, buffer_im);

      // D
      buffer_re = _mm_load_ps((float *)(tmp+k2)+k+0*OPERATOR_COMPONENT_OFFSET_float);
      buffer_im = _mm_load_ps((float *)(tmp+k2)+k+1*OPERATOR_COMPONENT_OFFSET_float);
      for ( m=offset; m<2*offset; m++ ) {
        // spin_2_3 is the same for all k => broadcast
        __m128 spin_2_3_re = _mm_set1_ps(spin_2_3_pt[(2*m+0)*component_offset]);
        __m128 spin_2_3_im = _mm_set1_ps(spin_2_3_pt[(2*m+1)*component_offset]);
        __m128 interpolation_data_re = _mm_load_ps(interpolation_data + (2*m+0)*component_offset);
        __m128 interpolation_data_im = _mm_load_ps(interpolation_data + (2*m+1)*component_offset);

        cfmadd_conj(interpolation_data_re, interpolation_data_im, spin_2_3_re, spin_2_3_im, &buffer_re, &buffer_im);
      }
      _mm_store_ps((float *)(tmp+k2)+k+0*OPERATOR_COMPONENT_OFFSET_float, buffer_re);
      _mm_store_ps((float *)(tmp+k2)+k+1*OPERATOR_COMPONENT_OFFSET_float, buffer_im);
    }

    // index k used for vectorization
    for ( k=0; k<OPERATOR_COMPONENT_OFFSET_float; k+=SIMD_LENGTH_float ) {
      __m128 buffer_re;
      __m128 buffer_im;

      // this are the packed indices, which we do not use in tmp
      //k1 = component_offset*(num_eig_vect+1+n);
      k1 = (n+2*num_eig_vect)*OPERATOR_COMPONENT_OFFSET_float;

      interpolation_data = (float *)(V + k*l->vector_size + fine_components*component_offset*site);

      // B
      buffer_re = _mm_load_ps((float *)(tmp+k1)+k+0*OPERATOR_COMPONENT_OFFSET_float);
      buffer_im = _mm_load_ps((float *)(tmp+k1)+k+1*OPERATOR_COMPONENT_OFFSET_float);
      for ( m=0; m<offset; m++ ) {
        // spin_2_3 is the same for all k => broadcast
        __m128 spin_2_3_re = _mm_set1_ps(spin_2_3_pt[(2*m+0)*component_offset]);
        __m128 spin_2_3_im = _mm_set1_ps(spin_2_3_pt[(2*m+1)*component_offset]);
        __m128 interpolation_data_re = _mm_load_ps(interpolation_data + (2*m+0)*component_offset);
        __m128 interpolation_data_im = _mm_load_ps(interpolation_data + (2*m+1)*component_offset);

        cfmadd_conj(interpolation_data_re, interpolation_data_im, spin_2_3_re, spin_2_3_im, &buffer_re, &buffer_im);
      }
      _mm_store_ps((float *)(tmp+k1)+k+0*OPERATOR_COMPONENT_OFFSET_float, buffer_re);
      _mm_store_ps((float *)(tmp+k1)+k+1*OPERATOR_COMPONENT_OFFSET_float, buffer_im);
    }
  }
#endif
}


static inline void sse_set_coarse_neighbor_coupling_float( complex_float *spin_0_1, complex_float *spin_2_3,
    complex_float *V, const int mu, level_struct *l, int site, const int n_rhs, complex_float *tmp ) {

#ifdef SSE
  int k, k1, k2, m, num_eig_vect = l->next_level->num_lattice_site_var/2,
      offset = l->num_lattice_site_var/2;

  float *spin_0_1_pt;
  float *spin_2_3_pt;
  float *interpolation_data;

  int component_offset = SIMD_LENGTH_float;
  int fine_components = l->num_lattice_site_var;

  // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
  //             C D ]                        -B*  D* ]
  // storage order: A, C, B, D, each column wise
  for ( int n=0; n<n_rhs; n++ ) {

    spin_0_1_pt = (float *)(spin_0_1 + (n/SIMD_LENGTH_float)*2*offset*SIMD_LENGTH_float) + n%SIMD_LENGTH_float;
    spin_2_3_pt = (float *)(spin_2_3 + (n/SIMD_LENGTH_float)*2*offset*SIMD_LENGTH_float) + n%SIMD_LENGTH_float;

    // index k used for vectorization
    for ( k=0; k<OPERATOR_COMPONENT_OFFSET_float; k+=SIMD_LENGTH_float ) {
      __m128 buffer_re;
      __m128 buffer_im;

      interpolation_data = (float *)(V + k*l->vector_size + fine_components*component_offset*site);

      k1 = (n+0*num_eig_vect)*OPERATOR_COMPONENT_OFFSET_float;
      k2 = (n+1*num_eig_vect)*OPERATOR_COMPONENT_OFFSET_float;

      // A
      buffer_re = _mm_load_ps((float *)(tmp+k1)+k+0*OPERATOR_COMPONENT_OFFSET_float);
      buffer_im = _mm_load_ps((float *)(tmp+k1)+k+1*OPERATOR_COMPONENT_OFFSET_float);
      for ( m=0; m<offset; m++ ) {
        // spin_0_1 is the same for all k => broadcast
        __m128 spin_0_1_re = _mm_set1_ps(spin_0_1_pt[(2*m+0)*component_offset]);
        __m128 spin_0_1_im = _mm_set1_ps(spin_0_1_pt[(2*m+1)*component_offset]);
        __m128 interpolation_data_re = _mm_load_ps(interpolation_data + (2*m+0)*component_offset);
        __m128 interpolation_data_im = _mm_load_ps(interpolation_data + (2*m+1)*component_offset);

        cfmadd_conj(interpolation_data_re, interpolation_data_im, spin_0_1_re, spin_0_1_im, &buffer_re, &buffer_im);
      }
      _mm_store_ps((float *)(tmp+k1)+k+0*OPERATOR_COMPONENT_OFFSET_float, buffer_re);
      _mm_store_ps((float *)(tmp+k1)+k+1*OPERATOR_COMPONENT_OFFSET_float, buffer_im);

      // C
      buffer_re = _mm_load_ps((float *)(tmp+k2)+k+0*OPERATOR_COMPONENT_OFFSET_float);
      buffer_im = _mm_load_ps((float *)(tmp+k2)+k+1*OPERATOR_COMPONENT_OFFSET_float);
      for ( m=offset; m<2*offset; m++ ) {
        // spin_0_1 is the same for all k => broadcast
        __m128 spin_0_1_re = _mm_set1_ps(spin_0_1_pt[(2*m+0)*component_offset]);
        __m128 spin_0_1_im = _mm_set1_ps(spin_0_1_pt[(2*m+1)*component_offset]);
        __m128 interpolation_data_re = _mm_load_ps(interpolation_data + (2*m+0)*component_offset);
        __m128 interpolation_data_im = _mm_load_ps(interpolation_data + (2*m+1)*component_offset);

        cfmadd_conj(interpolation_data_re, interpolation_data_im, spin_0_1_re, spin_0_1_im, &buffer_re, &buffer_im);
      }
      _mm_store_ps((float *)(tmp+k2)+k+0*OPERATOR_COMPONENT_OFFSET_float, buffer_re);
      _mm_store_ps((float *)(tmp+k2)+k+1*OPERATOR_COMPONENT_OFFSET_float, buffer_im);


      k1 = (n+2*num_eig_vect)*OPERATOR_COMPONENT_OFFSET_float;
      k2 = (n+3*num_eig_vect)*OPERATOR_COMPONENT_OFFSET_float;

      // B
      buffer_re = _mm_load_ps((float *)(tmp+k1)+k+0*OPERATOR_COMPONENT_OFFSET_float);
      buffer_im = _mm_load_ps((float *)(tmp+k1)+k+1*OPERATOR_COMPONENT_OFFSET_float);
      for ( m=0; m<offset; m++ ) {
        // spin_2_3 is the same for all k => broadcast
        __m128 spin_2_3_re = _mm_set1_ps(spin_2_3_pt[(2*m+0)*component_offset]);
        __m128 spin_2_3_im = _mm_set1_ps(spin_2_3_pt[(2*m+1)*component_offset]);
        __m128 interpolation_data_re = _mm_load_ps(interpolation_data + (2*m+0)*component_offset);
        __m128 interpolation_data_im = _mm_load_ps(interpolation_data + (2*m+1)*component_offset);

        cfmadd_conj(interpolation_data_re, interpolation_data_im, spin_2_3_re, spin_2_3_im, &buffer_re, &buffer_im);
      }
      _mm_store_ps((float *)(tmp+k1)+k+0*OPERATOR_COMPONENT_OFFSET_float, buffer_re);
      _mm_store_ps((float *)(tmp+k1)+k+1*OPERATOR_COMPONENT_OFFSET_float, buffer_im);

      // D
      buffer_re = _mm_load_ps((float *)(tmp+k2)+k+0*OPERATOR_COMPONENT_OFFSET_float);
      buffer_im = _mm_load_ps((float *)(tmp+k2)+k+1*OPERATOR_COMPONENT_OFFSET_float);
      for ( m=offset; m<2*offset; m++ ) {
        // spin_2_3 is the same for all k => broadcast
        __m128 spin_2_3_re = _mm_set1_ps(spin_2_3_pt[(2*m+0)*component_offset]);
        __m128 spin_2_3_im = _mm_set1_ps(spin_2_3_pt[(2*m+1)*component_offset]);
        __m128 interpolation_data_re = _mm_load_ps(interpolation_data + (2*m+0)*component_offset);
        __m128 interpolation_data_im = _mm_load_ps(interpolation_data + (2*m+1)*component_offset);

        cfmadd_conj(interpolation_data_re, interpolation_data_im, spin_2_3_re, spin_2_3_im, &buffer_re, &buffer_im);
      }
      _mm_store_ps((float *)(tmp+k2)+k+0*OPERATOR_COMPONENT_OFFSET_float, buffer_re);
      _mm_store_ps((float *)(tmp+k2)+k+1*OPERATOR_COMPONENT_OFFSET_float, buffer_im);
    }
  }
#endif
}


static inline void sse_coarse_spinwise_site_self_couplings_float( complex_float *eta1, complex_float *eta2,
    complex_float *phi, config_float clover, int elements, level_struct *l ) {

#ifdef SSE
  int num_eig_vect = l->num_lattice_site_var/2;
  int clover_step_size1 = (num_eig_vect * (num_eig_vect+1))/2;
  complex_float *eta[2] = {eta1, eta2};
  // U(x) = [ A B      , A=A*, D=D*, C = -B*
  //          C D ]
  // storage order: upper triangle of A, upper triangle of D, B, columnwise

  __m128 clover_re;
  __m128 clover_im;
  __m128 in_re;
  __m128 in_im;
  __m128 out_re;
  __m128 out_im;

  // zero output matrices
  __m128 zero = _mm_setzero_ps();
  for(int s=0; s<2; s++) {
    for(int i=0; i<elements; i+=SIMD_LENGTH_float) {
      for(int row=0; row<2*num_eig_vect; row++) {
        _mm_store_ps((float *)eta[s] + i + (2*row+0)*elements, zero);
        _mm_store_ps((float *)eta[s] + i + (2*row+1)*elements, zero);
      }
    }
  }

  // s refers to "spin" components 0and1 (->eta1) or 2and3 (->eta2)
  eta[1] += num_eig_vect*elements;
  for(int s=0; s<2; s++) {
    // A and D: column major hermitian, stored as upper triangular
    for(int i=0; i<elements; i+=SIMD_LENGTH_float) {
      for(int column=0; column<num_eig_vect; column++) {
        in_re  = _mm_load_ps((float *)phi + i + (2*column+0)*elements);
        in_im  = _mm_load_ps((float *)phi + i + (2*column+1)*elements);
        for(int row=0; row<=column; row++) {
          out_re = _mm_load_ps((float *)eta[s] + i + (2*row+0)*elements);
          out_im = _mm_load_ps((float *)eta[s] + i + (2*row+1)*elements);
          clover_re = _mm_set1_ps(creal(clover[(column*column+column)/2+row]));
          clover_im = _mm_set1_ps(cimag(clover[(column*column+column)/2+row]));

          cfmadd(clover_re, clover_im, in_re, in_im, &out_re, &out_im);

          _mm_store_ps((float *)eta[s] + i + (2*row+0)*elements, out_re);
          _mm_store_ps((float *)eta[s] + i + (2*row+1)*elements, out_im);
        }
        for(int row=column+1; row<num_eig_vect; row++) {
          out_re = _mm_load_ps((float *)eta[s] + i + (2*row+0)*elements);
          out_im = _mm_load_ps((float *)eta[s] + i + (2*row+1)*elements);
          clover_re = _mm_set1_ps(creal(clover[(row*row+row)/2+column]));
          clover_im = _mm_set1_ps(cimag(clover[(row*row+row)/2+column]));

          cfmadd_conj(clover_re, clover_im, in_re, in_im, &out_re, &out_im);

          _mm_store_ps((float *)eta[s] + i + (2*row+0)*elements, out_re);
          _mm_store_ps((float *)eta[s] + i + (2*row+1)*elements, out_im);
        }
      }
    }
    clover += clover_step_size1;
    phi += num_eig_vect*elements;
  }
  // rewind phi back to upper components
  phi -= 2*num_eig_vect*elements;
  eta[0] += num_eig_vect*elements;
  eta[1] -= num_eig_vect*elements;
  // C = -B^{\dagger}
  for(int i=0; i<elements; i+=SIMD_LENGTH_float) {
    for(int column=0; column<num_eig_vect; column++) {
      in_re  = _mm_load_ps((float *)phi + i + (2*column+0)*elements);
      in_im  = _mm_load_ps((float *)phi + i + (2*column+1)*elements);
      for(int row=0; row<num_eig_vect; row++) {
        out_re = _mm_load_ps((float *)eta[0] + i + (2*row+0)*elements);
        out_im = _mm_load_ps((float *)eta[0] + i + (2*row+1)*elements);
        // load transposed B
        clover_re = _mm_set1_ps(creal(clover[row*num_eig_vect+column]));
        clover_im = _mm_set1_ps(cimag(clover[row*num_eig_vect+column]));

        cfnmadd_conj(clover_re, clover_im, in_re, in_im, &out_re, &out_im);

        _mm_store_ps((float *)eta[0] + i + (2*row+0)*elements, out_re);
        _mm_store_ps((float *)eta[0] + i + (2*row+1)*elements, out_im);
      }
    }
  }
  phi += num_eig_vect*elements;
  // B
  for(int i=0; i<elements; i+=SIMD_LENGTH_float) {
    for(int column=0; column<num_eig_vect; column++) {
      in_re  = _mm_load_ps((float *)phi + i + (2*column+0)*elements);
      in_im  = _mm_load_ps((float *)phi + i + (2*column+1)*elements);
      for(int row=0; row<num_eig_vect; row++) {
        out_re = _mm_load_ps((float *)eta[1] + i + (2*row+0)*elements);
        out_im = _mm_load_ps((float *)eta[1] + i + (2*row+1)*elements);
        clover_re = _mm_set1_ps(creal(clover[column*num_eig_vect+row]));
        clover_im = _mm_set1_ps(cimag(clover[column*num_eig_vect+row]));

        cfmadd(clover_re, clover_im, in_re, in_im, &out_re, &out_im);

        _mm_store_ps((float *)eta[1] + i + (2*row+0)*elements, out_re);
        _mm_store_ps((float *)eta[1] + i + (2*row+1)*elements, out_im);
      }
    }
  }
#endif
}

#endif //SSE
#endif