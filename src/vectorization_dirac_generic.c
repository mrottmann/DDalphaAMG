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

#include "main.h"

#ifdef SSE
void d_plus_clover_aggregate_PRECISION_vectorized( complex_PRECISION *eta1, complex_PRECISION *eta2,
    complex_PRECISION *phi, schwarz_PRECISION_struct *s, level_struct *l,
    int site, int *direction_flags ) {

  int offset = SIMD_LENGTH_PRECISION;
  int site_offset = 12*offset;
  int index_out;
  int index_bw;
  int index_fw;
  int *neighbor = s->op.neighbor_table;
  int *backward_neighbor = s->op.backward_neighbor_table;
  complex_PRECISION *phi_pt;
  complex_PRECISION buffer1[site_offset] __attribute__((aligned(64)));
  complex_PRECISION buffer2[site_offset] __attribute__((aligned(64)));
  config_PRECISION D_pt;
  config_PRECISION D = s->op.D;

  // add clover term/shift
  spin0and1_site_clover_PRECISION_vectorized( eta1, phi+site_offset*site, s->op.clover+42*site, s->op.shift, offset );
  spin2and3_site_clover_PRECISION_vectorized( eta2, phi+site_offset*site, s->op.clover+42*site, s->op.shift, offset );

  index_out = site;

  for(int mu=0; mu<4; mu++) {
    index_fw  = neighbor[4*index_out + mu];
    index_bw  = backward_neighbor[4*index_out + mu];

    // from backward
    if ( direction_flags[2*mu+0] == 1 ) {
      D_pt = D + 36*index_bw+9*mu;
      phi_pt = phi + site_offset*index_bw;
      mvmh_PRECISION_vectorized( buffer2+0*offset, D_pt, phi_pt+0*offset, offset );
      mvmh_PRECISION_vectorized( buffer2+3*offset, D_pt, phi_pt+3*offset, offset );
      mvmh_PRECISION_vectorized( buffer2+6*offset, D_pt, phi_pt+6*offset, offset );
      mvmh_PRECISION_vectorized( buffer2+9*offset, D_pt, phi_pt+9*offset, offset );
      twospin_PRECISION_vectorized( eta1, eta2, buffer2, offset, mu, -1.0 );
    }

    // from forward
    if ( direction_flags[2*mu+1] == 1 ) {
      D_pt = D + 36*index_out+9*mu;
      phi_pt = phi + site_offset*index_fw;
      mvm_PRECISION_vectorized( buffer1+0*offset, D_pt, phi_pt+0*offset, offset );
      mvm_PRECISION_vectorized( buffer1+3*offset, D_pt, phi_pt+3*offset, offset );
      mvm_PRECISION_vectorized( buffer1+6*offset, D_pt, phi_pt+6*offset, offset );
      mvm_PRECISION_vectorized( buffer1+9*offset, D_pt, phi_pt+9*offset, offset );
      twospin_PRECISION_vectorized( eta1, eta2, buffer1, offset, mu, 1.0 );
    }
  }
}
#endif

#ifdef SSE
void d_neighbor_aggregate_PRECISION_vectorized( complex_PRECISION *eta1, complex_PRECISION *eta2,
      complex_PRECISION *phi, const int mu, schwarz_PRECISION_struct *s, level_struct *l,
      int site ) {

  int offset = SIMD_LENGTH_PRECISION;
  int site_offset = 12*offset;
  int index_out;
  int index_fw;
  int *neighbor = s->op.neighbor_table;
  complex_PRECISION *phi_pt;
  complex_PRECISION buffer[site_offset] __attribute__((aligned(64)));
  config_PRECISION D_pt;
  config_PRECISION D = s->op.D;

  index_out = site;

  // requires the positive boundaries of phi to be communicated befor
  index_fw  = neighbor[4*index_out + mu];
  D_pt = D + 36*index_out+9*mu;
  phi_pt = phi + site_offset*index_fw;
  mvm_PRECISION_vectorized_simd_length( buffer+0*offset, D_pt, phi_pt+0*offset );
  mvm_PRECISION_vectorized_simd_length( buffer+3*offset, D_pt, phi_pt+3*offset );
  mvm_PRECISION_vectorized_simd_length( buffer+6*offset, D_pt, phi_pt+6*offset );
  mvm_PRECISION_vectorized_simd_length( buffer+9*offset, D_pt, phi_pt+9*offset );
  twospin2_p_PRECISION_vectorized_simd_length( eta1, eta2, buffer, mu );
}
#endif
