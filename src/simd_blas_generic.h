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

#ifndef SIMD_BLAS_PRECISION_HEADER
#define SIMD_BLAS_PRECISION_HEADER

static inline void cgem_inverse_PRECISION( const int N, PRECISION *A_inverse, PRECISION *A, int lda ) {

  // generate LU decomp in A
  int i, j, k;
  complex_PRECISION alpha;
  
  complex_PRECISION tmpA[N*N];
  complex_PRECISION tmpA_inverse[N*N];
  
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
  
  complex_PRECISION b[N];
  complex_PRECISION *x;
  
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


static inline void cgemv_PRECISION( const int N, const OPERATOR_TYPE_PRECISION *A, int lda,
                                    const PRECISION *B, PRECISION *C ) {
  int i, j;
  
  mm_PRECISION A_re;
  mm_PRECISION A_im;
  mm_PRECISION B_re;
  mm_PRECISION B_im;
  mm_PRECISION C_re[lda/SIMD_LENGTH_PRECISION];
  mm_PRECISION C_im[lda/SIMD_LENGTH_PRECISION];
  
  // deinterleaved load
  for ( i=0; i<lda; i+= SIMD_LENGTH_PRECISION )
    cload_PRECISION( C+2*i, &(C_re[i/SIMD_LENGTH_PRECISION]), &(C_im[i/SIMD_LENGTH_PRECISION]));
  
  for ( j=0; j<N; j++ ) {
    // load the j-th complex number in B
    B_re = mm_set1_PRECISION( B[2*j] );
    B_im = mm_set1_PRECISION( B[2*j+1] );
    
    for ( i=0; i<lda; i+= SIMD_LENGTH_PRECISION ) {
       A_re = mm_load_PRECISION( A + 2*j*lda + i );
       A_im = mm_load_PRECISION( A + (2*j+1)*lda + i );
       
       // C += A*B
       cfmadd_PRECISION( A_re, A_im, B_re, B_im, &(C_re[i/SIMD_LENGTH_PRECISION]), &(C_im[i/SIMD_LENGTH_PRECISION]));
    }
  }  
  
  // interleaves real and imaginary parts and stores the resulting complex numbers in C
  for ( i=0; i<lda; i+= SIMD_LENGTH_PRECISION )
    cstore_PRECISION( C+2*i, C_re[i/SIMD_LENGTH_PRECISION], C_im[i/SIMD_LENGTH_PRECISION] );
}

static inline void cgenmv_PRECISION( const int N, const OPERATOR_TYPE_PRECISION *A, int lda,
                                    const PRECISION *B, PRECISION *C ) {
  int i, j;
  
  mm_PRECISION A_re;
  mm_PRECISION A_im;
  mm_PRECISION B_re;
  mm_PRECISION B_im;
  mm_PRECISION C_re[lda/SIMD_LENGTH_PRECISION];
  mm_PRECISION C_im[lda/SIMD_LENGTH_PRECISION];
  
  // deinterleaved load
  for ( i=0; i<lda; i+= SIMD_LENGTH_PRECISION )
    cload_PRECISION( C+2*i, &(C_re[i/SIMD_LENGTH_PRECISION]), &(C_im[i/SIMD_LENGTH_PRECISION]));
  
  for ( j=0; j<N; j++ ) {
    // load the j-th complex number in B
    B_re = mm_set1_PRECISION( B[2*j] );
    B_im = mm_set1_PRECISION( B[2*j+1] );
    
    for ( i=0; i<lda; i+= SIMD_LENGTH_PRECISION ) {
       A_re = mm_load_PRECISION( A + 2*j*lda + i );
       A_im = mm_load_PRECISION( A + (2*j+1)*lda + i );
       
       // C += A*B
       cfnmadd_PRECISION( A_re, A_im, B_re, B_im, &(C_re[i/SIMD_LENGTH_PRECISION]), &(C_im[i/SIMD_LENGTH_PRECISION]));
    }
  }  
  
  // interleaves real and imaginary parts and stores the resulting complex numbers in C
  for ( i=0; i<lda; i+= SIMD_LENGTH_PRECISION )
    cstore_PRECISION( C+2*i, C_re[i/SIMD_LENGTH_PRECISION], C_im[i/SIMD_LENGTH_PRECISION] );
}

static inline void cgemv_padded_PRECISION( const int N, const OPERATOR_TYPE_PRECISION *A, int lda, int padded,
                                           const PRECISION *B, PRECISION *C ) {
  int i, j, ip;

  int offset = SIMD_LENGTH_PRECISION*((padded+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
#ifdef HAVE_TM1p1
  if( g.n_flavours == 2 ) {
    int jp;
    mm_PRECISION A_re;
    mm_PRECISION A_im;
    mm_PRECISION B1_re;
    mm_PRECISION B1_im;
    mm_PRECISION B2_re;
    mm_PRECISION B2_im;
    mm_PRECISION C1_re[lda/SIMD_LENGTH_PRECISION];
    mm_PRECISION C1_im[lda/SIMD_LENGTH_PRECISION];
    mm_PRECISION C2_re[lda/SIMD_LENGTH_PRECISION];
    mm_PRECISION C2_im[lda/SIMD_LENGTH_PRECISION];
    
    // deinterleaved load
    for ( i=0; i<lda; i+= SIMD_LENGTH_PRECISION ) {
      ip = i%offset + 2*(i/offset)*padded;
      cload_PRECISION( C+2*ip, &(C1_re[i/SIMD_LENGTH_PRECISION]), &(C1_im[i/SIMD_LENGTH_PRECISION]));
      cload_PRECISION( C+2*(ip+padded), &(C2_re[i/SIMD_LENGTH_PRECISION]), &(C2_im[i/SIMD_LENGTH_PRECISION]));
    }
    
    for ( j=0; j<N; j++ ) {
      // load the j-th complex number in B
      jp = j + (j/padded)*padded;
      B1_re = mm_set1_PRECISION( B[2*jp] );
      B1_im = mm_set1_PRECISION( B[2*jp+1] );
      B2_re = mm_set1_PRECISION( B[2*(jp+padded)] );
      B2_im = mm_set1_PRECISION( B[2*(jp+padded)+1] );
      
      for ( i=0; i<lda; i+= SIMD_LENGTH_PRECISION ) {
        A_re = mm_load_PRECISION( A + 2*j*lda + i );
        A_im = mm_load_PRECISION( A + (2*j+1)*lda + i );
        
        // C += A*B
        cfmadd_PRECISION(A_re, A_im, B1_re, B1_im, &(C1_re[i/SIMD_LENGTH_PRECISION]),
                         &(C1_im[i/SIMD_LENGTH_PRECISION]) );
        cfmadd_PRECISION(A_re, A_im, B2_re, B2_im, &(C2_re[i/SIMD_LENGTH_PRECISION]),
                         &(C2_im[i/SIMD_LENGTH_PRECISION]) );
      }
    }  
    
    // interleaves real and imaginary parts and stores the resulting complex numbers in C
    for ( j=0; j<lda/offset; j++ ) {
      // we save it from last to first in order to avoid overriting issues.
      for ( i = (j+1)*offset-SIMD_LENGTH_PRECISION; i >= j*offset; i -= SIMD_LENGTH_PRECISION ) {
        ip = i%offset + 2*(i/offset)*padded;
        cstore_PRECISION( C+2*ip, C1_re[i/SIMD_LENGTH_PRECISION], C1_im[i/SIMD_LENGTH_PRECISION] );
        cstore_PRECISION( C+2*(ip+padded), C2_re[i/SIMD_LENGTH_PRECISION], C2_im[i/SIMD_LENGTH_PRECISION] );
      }
    }
  } else {
#endif
    mm_PRECISION A_re;
    mm_PRECISION A_im;
    mm_PRECISION B_re;
    mm_PRECISION B_im;
    mm_PRECISION C_re[lda/SIMD_LENGTH_PRECISION];
    mm_PRECISION C_im[lda/SIMD_LENGTH_PRECISION];
    
    // deinterleaved load
    for ( i=0; i<lda; i+= SIMD_LENGTH_PRECISION ) {
      ip = i%offset + (i/offset)*padded;
      cload_PRECISION( C+2*ip, &(C_re[i/SIMD_LENGTH_PRECISION]), &(C_im[i/SIMD_LENGTH_PRECISION]));
    }
    
    for ( j=0; j<N; j++ ) {
      // load the j-th complex number in B
      B_re = mm_set1_PRECISION( B[2*j] );
      B_im = mm_set1_PRECISION( B[2*j+1] );
      
      for ( i=0; i<lda; i+= SIMD_LENGTH_PRECISION ) {
        A_re = mm_load_PRECISION( A + 2*j*lda + i );
        A_im = mm_load_PRECISION( A + (2*j+1)*lda + i );
        
        // C += A*B
        cfmadd_PRECISION(A_re, A_im, B_re, B_im, &(C_re[i/SIMD_LENGTH_PRECISION]), &(C_im[i/SIMD_LENGTH_PRECISION]) );
      }
    }  
    
    // interleaves real and imaginary parts and stores the resulting complex numbers in C
    for ( i=0; i<lda; i+= SIMD_LENGTH_PRECISION ) {
      ip = i%offset + (i/offset)*padded;
      cstore_PRECISION( C+2*ip, C_re[i/SIMD_LENGTH_PRECISION], C_im[i/SIMD_LENGTH_PRECISION] );
    }
#ifdef HAVE_TM1p1
  }
#endif
}


static inline void cgenmv_padded_PRECISION( const int N, const OPERATOR_TYPE_PRECISION *A, int lda, int padded,
                                           const PRECISION *B, PRECISION *C ) {
  int i, j, ip;

  int offset = SIMD_LENGTH_PRECISION*((padded+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
#ifdef HAVE_TM1p1
  if( g.n_flavours == 2 ) {
    int jp;
    mm_PRECISION A_re;
    mm_PRECISION A_im;
    mm_PRECISION B1_re;
    mm_PRECISION B1_im;
    mm_PRECISION B2_re;
    mm_PRECISION B2_im;
    mm_PRECISION C1_re[lda/SIMD_LENGTH_PRECISION];
    mm_PRECISION C1_im[lda/SIMD_LENGTH_PRECISION];
    mm_PRECISION C2_re[lda/SIMD_LENGTH_PRECISION];
    mm_PRECISION C2_im[lda/SIMD_LENGTH_PRECISION];
    
    // deinterleaved load
    for ( i=0; i<lda; i+= SIMD_LENGTH_PRECISION ) {
      ip = i%offset + 2*(i/offset)*padded;
      cload_PRECISION( C+2*ip, &(C1_re[i/SIMD_LENGTH_PRECISION]), &(C1_im[i/SIMD_LENGTH_PRECISION]));
      cload_PRECISION( C+2*(ip+padded), &(C2_re[i/SIMD_LENGTH_PRECISION]), &(C2_im[i/SIMD_LENGTH_PRECISION]));
    }
    
    for ( j=0; j<N; j++ ) {
      // load the j-th complex number in B
      jp = j + (j/padded)*padded;
      B1_re = mm_set1_PRECISION( B[2*jp] );
      B1_im = mm_set1_PRECISION( B[2*jp+1] );
      B2_re = mm_set1_PRECISION( B[2*(jp+padded)] );
      B2_im = mm_set1_PRECISION( B[2*(jp+padded)+1] );
      
      for ( i=0; i<lda; i+= SIMD_LENGTH_PRECISION ) {
        A_re = mm_load_PRECISION( A + 2*j*lda + i );
        A_im = mm_load_PRECISION( A + (2*j+1)*lda + i );
        
        // C += A*B
        cfnmadd_PRECISION(A_re, A_im, B1_re, B1_im, &(C1_re[i/SIMD_LENGTH_PRECISION]),
                         &(C1_im[i/SIMD_LENGTH_PRECISION]) );
        cfnmadd_PRECISION(A_re, A_im, B2_re, B2_im, &(C2_re[i/SIMD_LENGTH_PRECISION]),
                         &(C2_im[i/SIMD_LENGTH_PRECISION]) );
      }
    }  
    
    // interleaves real and imaginary parts and stores the resulting complex numbers in C
    for ( j=0; j<lda/offset; j++ ) {
      // we save it from last to first in order to avoid overriting issues.
      for ( i = (j+1)*offset-SIMD_LENGTH_PRECISION; i >= j*offset; i -= SIMD_LENGTH_PRECISION ) {
        ip = i%offset + 2*(i/offset)*padded;
        cstore_PRECISION( C+2*ip, C1_re[i/SIMD_LENGTH_PRECISION], C1_im[i/SIMD_LENGTH_PRECISION] );
        cstore_PRECISION( C+2*(ip+padded), C2_re[i/SIMD_LENGTH_PRECISION], C2_im[i/SIMD_LENGTH_PRECISION] );
      }
    }
  } else {
#endif
    mm_PRECISION A_re;
    mm_PRECISION A_im;
    mm_PRECISION B_re;
    mm_PRECISION B_im;
    mm_PRECISION C_re[lda/SIMD_LENGTH_PRECISION];
    mm_PRECISION C_im[lda/SIMD_LENGTH_PRECISION];
    
    // deinterleaved load
    for ( i=0; i<lda; i+= SIMD_LENGTH_PRECISION ) {
      ip = i%offset + (i/offset)*padded;
      cload_PRECISION( C+2*ip, &(C_re[i/SIMD_LENGTH_PRECISION]), &(C_im[i/SIMD_LENGTH_PRECISION]));
    }
    
    for ( j=0; j<N; j++ ) {
      // load the j-th complex number in B
      B_re = mm_set1_PRECISION( B[2*j] );
      B_im = mm_set1_PRECISION( B[2*j+1] );
      
      for ( i=0; i<lda; i+= SIMD_LENGTH_PRECISION ) {
        A_re = mm_load_PRECISION( A + 2*j*lda + i );
        A_im = mm_load_PRECISION( A + (2*j+1)*lda + i );
        
        // C += A*B
        cfnmadd_PRECISION(A_re, A_im, B_re, B_im, &(C_re[i/SIMD_LENGTH_PRECISION]), &(C_im[i/SIMD_LENGTH_PRECISION]) );
      }
    }  
    
    // interleaves real and imaginary parts and stores the resulting complex numbers in C
    for ( i=0; i<lda; i+= SIMD_LENGTH_PRECISION ) {
      ip = i%offset + (i/offset)*padded;
      cstore_PRECISION( C+2*ip, C_re[i/SIMD_LENGTH_PRECISION], C_im[i/SIMD_LENGTH_PRECISION] );
    }
#ifdef HAVE_TM1p1
  }
#endif
}

#endif
