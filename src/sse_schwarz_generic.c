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

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_PRECISION
void block_PRECISION_boundary_op( vector_PRECISION eta, vector_PRECISION phi, int k,
                                  schwarz_PRECISION_struct *s, level_struct *l ) {
  int *bbl = s->block_boundary_length;
  PRECISION *Dplus = s->op.D_vectorized;
  PRECISION *Dminus = s->op.D_transformed_vectorized;
  
  for ( int mu=0; mu<4; mu++ ) {
    boundary_plus_coupling_PRECISION( (PRECISION*)eta, Dplus, (PRECISION*)phi,
                                              mu, bbl[2*mu], bbl[2*mu+1], s->block[k].bt, NULL );
    boundary_minus_coupling_PRECISION( (PRECISION*)eta, Dminus, (PRECISION*)phi,
                                               mu, bbl[2*mu+1], bbl[2*mu+2], s->block[k].bt, NULL );
  }
}
#endif

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_PRECISION
void n_block_PRECISION_boundary_op( vector_PRECISION eta, vector_PRECISION phi, int k,
                                    schwarz_PRECISION_struct *s, level_struct *l ) {
  int *bbl = s->block_boundary_length;
  PRECISION *Dplus = s->op.D_vectorized;
  PRECISION *Dminus = s->op.D_transformed_vectorized;
  
  for ( int mu=0; mu<4; mu++ ) {
    boundary_nplus_coupling_PRECISION( (PRECISION*)eta, Dplus, (PRECISION*)phi,
                                               mu, bbl[2*mu], bbl[2*mu+1], s->block[k].bt, NULL );
    boundary_nminus_coupling_PRECISION( (PRECISION*)eta, Dminus, (PRECISION*)phi,
                                                mu, bbl[2*mu+1], bbl[2*mu+2], s->block[k].bt, NULL );
  }
}
#endif

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_PRECISION
void coarse_block_PRECISION_boundary_op( vector_PRECISION eta, vector_PRECISION phi,
                                         int k, schwarz_PRECISION_struct *s, level_struct *l ) {
  // k: number of current block
  int *bbl = s->block_boundary_length, n = l->num_lattice_site_var;
  int column_offset = SIMD_LENGTH_PRECISION*((l->num_lattice_site_var+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  int vectorized_link_offset = 2*l->num_lattice_site_var*column_offset;
  
  for ( int mu=0; mu<4; mu++ ) {
    OPERATOR_TYPE_PRECISION *Dplus = s->op.D_vectorized + mu*vectorized_link_offset;
    OPERATOR_TYPE_PRECISION *Dminus = s->op.D_transformed_vectorized + mu*vectorized_link_offset;
    // plus mu direction
    for ( int i=bbl[2*mu]; i<bbl[2*mu+1]; i+=2 ) {
      int index = s->block[k].bt[i];
      int neighbor_index = s->block[k].bt[i+1];
      vector_PRECISION phi_pt = phi + n*neighbor_index;
      vector_PRECISION eta_pt = eta + n*index;
      coarse_hopp_PRECISION_vectorized( eta_pt, phi_pt, Dplus + 4*vectorized_link_offset*index, l );
    }
    // minus mu direction
    for ( int i=bbl[2*mu+1]; i<bbl[2*mu+2]; i+=2 ) {
      int index = s->block[k].bt[i];
      int neighbor_index = s->block[k].bt[i+1];
      vector_PRECISION phi_pt = phi + n*neighbor_index;
      vector_PRECISION eta_pt = eta + n*index;
      coarse_hopp_PRECISION_vectorized( eta_pt, phi_pt, Dminus + 4*vectorized_link_offset*neighbor_index, l );
    }
  }
}
#endif

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_PRECISION
void n_coarse_block_PRECISION_boundary_op( vector_PRECISION eta, vector_PRECISION phi,
                                           int k, schwarz_PRECISION_struct *s, level_struct *l ) {
  // k: number of current block
  int *bbl = s->block_boundary_length, n = l->num_lattice_site_var;
  int column_offset = SIMD_LENGTH_PRECISION*((l->num_lattice_site_var+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  int vectorized_link_offset = 2*l->num_lattice_site_var*column_offset;
  
  for ( int mu=0; mu<4; mu++ ) {
    OPERATOR_TYPE_PRECISION *Dplus = s->op.D_vectorized + mu*vectorized_link_offset;
    OPERATOR_TYPE_PRECISION *Dminus = s->op.D_transformed_vectorized + mu*vectorized_link_offset;
    // plus mu direction
    for ( int i=bbl[2*mu]; i<bbl[2*mu+1]; i+=2 ) {
      int index = s->block[k].bt[i];
      int neighbor_index = s->block[k].bt[i+1];
      vector_PRECISION phi_pt = phi + n*neighbor_index;
      vector_PRECISION eta_pt = eta + n*index;
      coarse_n_hopp_PRECISION_vectorized( eta_pt, phi_pt, Dplus + 4*vectorized_link_offset*index, l );
    }
    // minus mu direction
    for ( int i=bbl[2*mu+1]; i<bbl[2*mu+2]; i+=2 ) {
      int index = s->block[k].bt[i];
      int neighbor_index = s->block[k].bt[i+1];
      vector_PRECISION phi_pt = phi + n*neighbor_index;
      vector_PRECISION eta_pt = eta + n*index;
      coarse_n_hopp_PRECISION_vectorized( eta_pt, phi_pt, Dminus + 4*vectorized_link_offset*neighbor_index, l );
    }
  }
}
#endif

#if defined(OPTIMIZED_NEIGHBOR_COUPLING_PRECISION) || defined(OPTIMIZED_SELF_COUPLING_PRECISION)
void schwarz_PRECISION_setup( schwarz_PRECISION_struct *s, operator_double_struct *op_in, level_struct *l ) {

/*********************************************************************************  
* Copies the Dirac operator and the clover term from op_in into the Schwarz 
* struct (this function is depth 0 only).
* - operator_double_struct *op_in: Input operator.                                  
*********************************************************************************/

  int i, index, n = l->num_inner_lattice_sites, *tt = s->op.translation_table;
  config_PRECISION D_out_pt, clover_out_pt;
  config_double D_in_pt = op_in->D, clover_in_pt = op_in->clover;
  s->op.shift = op_in->shift;
  
  for ( i=0; i<n; i++ ) {
    index = tt[i];
    D_out_pt = s->op.D + 36*index;
    FOR36( *D_out_pt = (complex_PRECISION) *D_in_pt; D_out_pt++; D_in_pt++; )
  }
  
#ifdef OPTIMIZED_NEIGHBOR_COUPLING_PRECISION
  for ( i=0; i<n; i++ ) {
    PRECISION *D_vectorized = s->op.D_vectorized + 96*i;
    PRECISION *D_transformed_vectorized = s->op.D_transformed_vectorized + 96*i;
    complex_PRECISION *D_out_pt = s->op.D + 36*i;
    for ( int mu=0; mu<4; mu++ ) {
      set_PRECISION_D_vectorized( D_vectorized+24*mu, D_transformed_vectorized+24*mu, D_out_pt+9*mu );
    }
  }
#endif
  
  if ( g.csw != 0 ) {
    for ( i=0; i<n; i++ ) {
      index = tt[i];
      clover_out_pt = s->op.clover + 42*index;
#ifdef OPTIMIZED_SELF_COUPLING_PRECISION
      PRECISION *clover_out_vectorized_pt = s->op.clover_vectorized + 144*index;
      sse_set_clover_PRECISION( clover_out_vectorized_pt, clover_in_pt );
#endif
      FOR42( *clover_out_pt = (complex_PRECISION) *clover_in_pt; clover_out_pt++; clover_in_pt++; )
    }
  } else {
    for ( i=0; i<n; i++ ) {
      index = tt[i];
      clover_out_pt = s->op.clover + 12*index;
      FOR12( *clover_out_pt = (complex_PRECISION) *clover_in_pt; clover_out_pt++; clover_in_pt++; )
    }
  }
  
  if ( g.odd_even )
    schwarz_PRECISION_oddeven_setup( &(s->op), l );
  
  schwarz_PRECISION_boundary_update( s, l );

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_PRECISION
  int start = l->num_lattice_sites;
  int end = 2*l->num_lattice_sites - l->num_inner_lattice_sites;
  for ( i=start; i<end; i++ ) {
    PRECISION *D_vectorized = s->op.D_vectorized + 96*i;
    PRECISION *D_transformed_vectorized = s->op.D_transformed_vectorized + 96*i;
    complex_PRECISION *D_out_pt = s->op.D + 36*i;
    for ( int mu=0; mu<4; mu++ ) {
      set_PRECISION_D_vectorized( D_vectorized+24*mu, D_transformed_vectorized+24*mu, D_out_pt+9*mu );
    }
  }
#endif
}
#endif

#endif // SSE
