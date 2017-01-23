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

#ifndef VECTORIZATION_CONTROL_H
#define VECTORIZATION_CONTROL_H

#ifdef SSE

#define SIMD_LENGTH_float 4
#define SIMD_LENGTH_double 2

#define INTERPOLATION_OPERATOR_LAYOUT_OPTIMIZED_float
#define INTERPOLATION_SETUP_LAYOUT_OPTIMIZED_float
#define VECTORIZE_COARSE_OPERATOR_float
#define GRAM_SCHMIDT_VECTORIZED_float
#define OPTIMIZED_NEIGHBOR_COUPLING_float
#define OPTIMIZED_SELF_COUPLING_float
#define OPTIMIZED_NEIGHBOR_COUPLING_double
#define OPTIMIZED_LINALG_float
#define OPTIMIZED_LINALG_double

#include "sse_complex_float_intrinsic.h"
#include "sse_complex_double_intrinsic.h"

#endif

#define OPERATOR_COMPONENT_OFFSET_float  (SIMD_LENGTH_float *((l->num_eig_vect+SIMD_LENGTH_float -1)/SIMD_LENGTH_float ))
#define OPERATOR_COMPONENT_OFFSET_double (SIMD_LENGTH_double*((l->num_eig_vect+SIMD_LENGTH_double-1)/SIMD_LENGTH_double))

#define OPERATOR_TYPE_float float
#define OPERATOR_TYPE_double double

#endif // VECTORIZATION_CONTROL_H
