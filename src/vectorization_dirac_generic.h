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

#ifndef VECTORIZATION_DIRAC_PRECISION_HEADER
  #define VECTORIZATION_DIRAC_PRECISION_HEADER

#ifdef SSE
  #include "sse_dirac.h"
#endif

  // caller is responsibel for checking that he needs coupling in this direction for this site
  void d_neighbor_aggregate_PRECISION_vectorized( complex_PRECISION *eta1, complex_PRECISION *eta2,
      complex_PRECISION *phi, const int mu, schwarz_PRECISION_struct *s, level_struct *l,
      int site );
  
  void d_plus_clover_aggregate_PRECISION_vectorized( complex_PRECISION *eta1, complex_PRECISION *eta2,
    complex_PRECISION *phi, schwarz_PRECISION_struct *s, level_struct *l,
    int site, int *direction_flags );

  // spinors are vectorized, gauge is same for all (use for multiple rhs)
  static inline void mvm_PRECISION_vectorized_simd_length(
      const complex_PRECISION *eta, const complex_PRECISION *D, const complex_PRECISION *phi ) {
#ifdef SSE
    sse_mvm_PRECISION_simd_length( eta, D, phi );
#endif
    
  }
  // spinors are vectorized, gauge is same for all (use for multiple rhs)
  static inline void mvm_PRECISION_vectorized(
      const complex_PRECISION *eta, const complex_PRECISION *D, const complex_PRECISION *phi, int elements ) {
#ifdef SSE
    sse_mvm_PRECISION( eta, D, phi, elements );
#endif
  }

  // spinors are vectorized, gauge is same for all (use for multiple rhs)
  static inline void mvmh_PRECISION_vectorized(
      const complex_PRECISION *eta, const complex_PRECISION *D, const complex_PRECISION *phi, int elements ) {
#ifdef SSE
    sse_mvmh_PRECISION( eta, D, phi, elements );
#endif
  }
  
  // mu is according to the enum for T,Z,Y,X defined in clifford.h
  static inline void twospin_PRECISION_vectorized( complex_PRECISION *out_spin0and1, complex_PRECISION *out_spin2and3, const complex_PRECISION *in, int elements, int mu, double sign ) {
#ifdef SSE
    sse_twospin_PRECISION( out_spin0and1, out_spin2and3, in, elements, mu, sign );
#endif
  }
  static inline void twospin_p_T_PRECISION_vectorized( complex_PRECISION *out_spin0and1, complex_PRECISION *out_spin2and3, const complex_PRECISION *in, int elements ) {
    twospin_PRECISION_vectorized( out_spin0and1, out_spin2and3, in, elements, T, 1.0);
  }
  static inline void twospin_n_T_PRECISION_vectorized( complex_PRECISION *out_spin0and1, complex_PRECISION *out_spin2and3, const complex_PRECISION *in, int elements ) {
    twospin_PRECISION_vectorized( out_spin0and1, out_spin2and3, in, elements, T, -1.0);
  }
  static inline void twospin_p_Z_PRECISION_vectorized( complex_PRECISION *out_spin0and1, complex_PRECISION *out_spin2and3, const complex_PRECISION *in, int elements ) {
    twospin_PRECISION_vectorized( out_spin0and1, out_spin2and3, in, elements, Z, 1.0);
  }
  static inline void twospin_n_Z_PRECISION_vectorized( complex_PRECISION *out_spin0and1, complex_PRECISION *out_spin2and3, const complex_PRECISION *in, int elements ) {
    twospin_PRECISION_vectorized( out_spin0and1, out_spin2and3, in, elements, Z, -1.0);
  }
  static inline void twospin_p_Y_PRECISION_vectorized( complex_PRECISION *out_spin0and1, complex_PRECISION *out_spin2and3, const complex_PRECISION *in, int elements ) {
    twospin_PRECISION_vectorized( out_spin0and1, out_spin2and3, in, elements, Y, 1.0);
  }
  static inline void twospin_n_Y_PRECISION_vectorized( complex_PRECISION *out_spin0and1, complex_PRECISION *out_spin2and3, const complex_PRECISION *in, int elements ) {
    twospin_PRECISION_vectorized( out_spin0and1, out_spin2and3, in, elements, Y, -1.0);
  }
  static inline void twospin_p_X_PRECISION_vectorized( complex_PRECISION *out_spin0and1, complex_PRECISION *out_spin2and3, const complex_PRECISION *in, int elements ) {
    twospin_PRECISION_vectorized( out_spin0and1, out_spin2and3, in, elements, X, 1.0);
  }
  static inline void twospin_n_X_PRECISION_vectorized( complex_PRECISION *out_spin0and1, complex_PRECISION *out_spin2and3, const complex_PRECISION *in, int elements ) {
    twospin_PRECISION_vectorized( out_spin0and1, out_spin2and3, in, elements, X, -1.0);
  }

  // mu is according to the enum for T,Z,Y,X defined in clifford.h
  static inline void twospin2_p_PRECISION_vectorized_simd_length( complex_PRECISION *out_spin0and1, complex_PRECISION *out_spin2and3, const complex_PRECISION *in, int mu ) {
#ifdef SSE
    sse_twospin2_p_PRECISION_simd_length( out_spin0and1, out_spin2and3, in, mu );
#endif
  }
  // mu is according to the enum for T,Z,Y,X defined in clifford.h
  static inline void twospin2_p_PRECISION_vectorized( complex_PRECISION *out_spin0and1, complex_PRECISION *out_spin2and3, const complex_PRECISION *in, int elements, int mu ) {
#ifdef SSE
    sse_twospin2_p_PRECISION( out_spin0and1, out_spin2and3, in, elements, mu );
#endif
  }
  static inline void twospin2_p_T_PRECISION_vectorized( complex_PRECISION *out_spin0and1, complex_PRECISION *out_spin2and3, const complex_PRECISION *in, int elements ) {
    twospin2_p_PRECISION_vectorized( out_spin0and1, out_spin2and3, in, elements, T);
  }
  static inline void twospin2_p_Z_PRECISION_vectorized( complex_PRECISION *out_spin0and1, complex_PRECISION *out_spin2and3, const complex_PRECISION *in, int elements ) {
    twospin2_p_PRECISION_vectorized( out_spin0and1, out_spin2and3, in, elements, Z);
  }
  static inline void twospin2_p_Y_PRECISION_vectorized( complex_PRECISION *out_spin0and1, complex_PRECISION *out_spin2and3, const complex_PRECISION *in, int elements ) {
    twospin2_p_PRECISION_vectorized( out_spin0and1, out_spin2and3, in, elements, Y);
  }
  static inline void twospin2_p_X_PRECISION_vectorized( complex_PRECISION *out_spin0and1, complex_PRECISION *out_spin2and3, const complex_PRECISION *in, int elements ) {
    twospin2_p_PRECISION_vectorized( out_spin0and1, out_spin2and3, in, elements, X);
  }
  
  static inline void spin0and1_site_clover_PRECISION_vectorized( const complex_PRECISION *eta, const complex_PRECISION *phi,
                                                                 const config_PRECISION clover, double shift, int elements ) {
#ifdef SSE
    sse_spin0and1_site_clover_PRECISION( eta, phi, clover, shift, elements );
#endif
  }
  
  static inline void spin2and3_site_clover_PRECISION_vectorized( const complex_PRECISION *eta, const complex_PRECISION *phi,
                                                                 const config_PRECISION clover, double shift, int elements ) {
#ifdef SSE
    sse_spin2and3_site_clover_PRECISION( eta, phi, clover, shift, elements );
#endif
  }
  
#endif