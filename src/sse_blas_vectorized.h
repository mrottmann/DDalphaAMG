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

#ifndef SSE_BLAS_VECTORIZED_H
#define SSE_BLAS_VECTORIZED_H
#ifdef SSE

static inline void sse_cgem_inverse( const int N, float *A_inverse, float *A, int lda ) {
  // generate LU decomp in A
  
  int i, j, k;
  complex_float alpha;
  
  complex_float tmpA[N*N];
  complex_float tmpA_inverse[N*N];
  
  for ( j=0; j<N; j++ ) {
    for ( i=0; i<N; i++ ) {
      tmpA[i+N*j] = A[2*j*lda+i] + _Complex_I * A[(2*j+1)*lda+i];
    }
  }
  
  // LU decomp in A
  for ( k=0; k<N-1; k++ ) {
    for ( i=k+1; i<N; i++ ) {
      // alpha = A_ik/A_kk
      alpha = tmpA[i+k*N]/tmpA[k+k*N];
      tmpA[i+k*N] = alpha;
      for ( j=k+1; j<N; j++ ) {
        // A_ij = A_ij - alpha * A_kj
        tmpA[i+j*N] -= alpha* tmpA[k+j*N];
      }
    }    
  } 
  
  complex_float b[N];
  complex_float *x;
  
  for ( k=0; k<N; k++ ) {
    b[k] = 0;
  }
  
  for ( k=0; k<N; k++ ) {
    x = tmpA_inverse+k*N;
    b[k] = 1;
    if ( k>0 )
      b[k-1] = 0;
    
    for ( i=0; i<N; i++ ) {
      x[i] = b[i];
      for ( j=0; j<i; j++ ) {
        // x_i = x_i - A_ij + x_j
        x[i] = x[i] - tmpA[i+j*N]*x[j];
      }
    } // i
    
    for ( i=N-1; i>=0; i-- ) {
      for ( j=i+1; j<N; j++ ) {
        // x_i = x_i - A_ij * x_j
        x[i] = x[i] - tmpA[i+j*N]*x[j];
      }
      // x_i = x_i / A_ii
      x[i] = x[i]/tmpA[i+i*N];
    } // i
  } // k
  
  for ( j=0; j<N; j++ ) {
    for ( i=0; i<N; i++ ) {
      A_inverse[i+2*j*lda] = creal(tmpA_inverse[i+j*N]);
      A_inverse[i+(2*j+1)*lda] = cimag(tmpA_inverse[i+j*N]);
    }
    for ( i=N; i<lda; i++ ) {
      A_inverse[i+2*j*lda] = 0.0;
      A_inverse[i+(2*j+1)*lda] = 0.0;
    }
  } 
}


static inline void sse_cgemv( const int N, const OPERATOR_TYPE_float *A, int lda, const float *B, float *C ) {
  int i, j;
  
  __m128 A_re;
  __m128 A_im;
  __m128 B_re;
  __m128 B_im;
  __m128 C_re[lda/SIMD_LENGTH_float];
  __m128 C_im[lda/SIMD_LENGTH_float];
  
  // deinterleaved load
  for ( i=0; i<lda; i+= SIMD_LENGTH_float ) {
    C_re[i/SIMD_LENGTH_float] = _mm_setr_ps(C[2*i], C[2*i+2], C[2*i+4], C[2*i+6] );
    C_im[i/SIMD_LENGTH_float] = _mm_setr_ps(C[2*i+1], C[2*i+3], C[2*i+5], C[2*i+7] );
  }
  
  for ( j=0; j<N; j++ ) {
    // load the j-th complex number in B
    B_re = _mm_set1_ps( B[2*j] );
    B_im = _mm_set1_ps( B[2*j+1] );
    
    for ( i=0; i<lda; i+= SIMD_LENGTH_float ) {
       A_re = _mm_load_ps( A + 2*j*lda + i );
       A_im = _mm_load_ps( A + (2*j+1)*lda + i );
       
       // C += A*B
       cfmadd(A_re, A_im, B_re, B_im, &(C_re[i/SIMD_LENGTH_float]), &(C_im[i/SIMD_LENGTH_float]) );
    }
  }  
  
  // interleaves real and imaginary parts and stores the resulting complex numbers in C
  for ( i=0; i<lda; i+= SIMD_LENGTH_float ) {
     A_re = _mm_unpacklo_ps( C_re[i/SIMD_LENGTH_float], C_im[i/SIMD_LENGTH_float] );
     A_im = _mm_unpackhi_ps( C_re[i/SIMD_LENGTH_float], C_im[i/SIMD_LENGTH_float] );
     _mm_store_ps( C+2*i,                   A_re );
     _mm_store_ps( C+2*i+SIMD_LENGTH_float, A_im );
  }
}


static inline void sse_cgenmv( const int N, const OPERATOR_TYPE_float *A,  int lda, const float *B, float *C ) {
  int i, j;
  
  __m128 A_re;
  __m128 A_im;
  __m128 B_re;
  __m128 B_im;
  __m128 C_re[lda/SIMD_LENGTH_float];
  __m128 C_im[lda/SIMD_LENGTH_float];
  
  for ( i=0; i<lda; i+= SIMD_LENGTH_float ) {
    C_re[i/SIMD_LENGTH_float] = _mm_setr_ps(C[2*i], C[2*i+2], C[2*i+4], C[2*i+6] );
    C_im[i/SIMD_LENGTH_float] = _mm_setr_ps(C[2*i+1], C[2*i+3], C[2*i+5], C[2*i+7] );
  }
  
  for ( j=0; j<N; j++ ) {
    
    B_re = _mm_set1_ps( B[2*j] );
    B_im = _mm_set1_ps( B[2*j+1] ); 
    
    for ( i=0; i<lda; i+= SIMD_LENGTH_float ) {
       A_re = _mm_load_ps( A + 2*j*lda + i );
       A_im = _mm_load_ps( A + (2*j+1)*lda + i );
       
       // C -= A*B
       cfnmadd(A_re, A_im, B_re, B_im, &(C_re[i/SIMD_LENGTH_float]), &(C_im[i/SIMD_LENGTH_float]) );
    }
  }  
  
  for ( i=0; i<lda; i+= SIMD_LENGTH_float ) {
     A_re = _mm_unpacklo_ps( C_re[i/SIMD_LENGTH_float], C_im[i/SIMD_LENGTH_float] );
     A_im = _mm_unpackhi_ps( C_re[i/SIMD_LENGTH_float], C_im[i/SIMD_LENGTH_float] );
     _mm_store_ps( C+2*i,                   A_re );
     _mm_store_ps( C+2*i+SIMD_LENGTH_float, A_im );
  }
}

#endif
#endif
