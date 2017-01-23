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

#ifndef SSE_SCHWARZ_PRECISION_H
#define SSE_SCHWARZ_PRECISION_H
#ifdef SSE

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_PRECISION
static inline void set_PRECISION_D_vectorized( PRECISION *out1, PRECISION *out2, complex_PRECISION *in ) {
  // out1: column major, out2: row major
  for ( int i=0; i<3; i++ ) { // column
    for ( int j=0; j<3; j++ ) { // row
      out1[8*i  +j] = creal_PRECISION(in[3*j+i]);
      out1[8*i+4+j] = cimag_PRECISION(in[3*j+i]);
      out2[8*i  +j] = creal_PRECISION(in[j+3*i]);
      out2[8*i+4+j] = cimag_PRECISION(in[j+3*i]);
    }
    out1[8*i+3] = 0.0;
    out1[8*i+7] = 0.0;
    out2[8*i+3] = 0.0;
    out2[8*i+7] = 0.0;
  }
}
#endif

#endif // SSE
#endif // LINALG_MIC_H
