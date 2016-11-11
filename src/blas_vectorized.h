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

#ifndef BLAS_VECTORIZED_H
#define BLAS_VECTORIZED_H

// BLAS naming convention: LDA = leading dimension of A
#ifdef SSE
#include "sse_blas_vectorized.h"
#endif

// C=A*B+C
static inline void cgemv(const int N, const OPERATOR_TYPE_float *A, int lda, const float *B, float *C)
{
#ifdef SSE
  sse_cgemv( N, A, lda, B, C );
#endif
}

// C=-A*B+C
static inline void cgenmv(const int N, const OPERATOR_TYPE_float *A, int lda, const float *B, float *C)
{
#ifdef SSE
  sse_cgenmv( N, A, lda, B, C );
#endif
}

// C=A*B+C with padded layout
static inline void cgemv_padded(const int N, const OPERATOR_TYPE_float *A, int lda, int padded, const float *B, float *C)
{
#ifdef SSE
  sse_cgemv_padded( N, A, lda, padded, B, C );
#endif
}

// C=-A*B+C with padded layout
static inline void cgenmv_padded(const int N, const OPERATOR_TYPE_float *A, int lda, int padded, const float *B, float *C)
{
#ifdef SSE
  sse_cgenmv_padded( N, A, lda, padded, B, C );
#endif
}


static inline void cgem_inverse(const int N, OPERATOR_TYPE_float *A_inverse, OPERATOR_TYPE_float *A, int lda)
{
#ifdef SSE
  sse_cgem_inverse( N, A_inverse, A, lda );
#endif
}

#endif // BLAS_VECTORIZED_H
