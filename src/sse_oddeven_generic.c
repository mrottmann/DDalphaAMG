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
void hopping_term_PRECISION( vector_PRECISION eta, vector_PRECISION phi, operator_PRECISION_struct *op,
                             const int amount, level_struct *l, struct Thread *threading ) {
  
  int start_even, end_even, start_odd, end_odd, n = l->num_inner_lattice_sites,
      *neighbor = op->neighbor_table, start=0, plus_dir_param=_FULL_SYSTEM, minus_dir_param=_FULL_SYSTEM;
  
  SYNC_CORES(threading)  
  
  if ( amount == _EVEN_SITES || amount == _ODD_SITES ) {
    compute_core_start_end_custom(0, op->num_even_sites, &start_even, &end_even, l, threading, 1 );
    compute_core_start_end_custom(op->num_even_sites, op->num_even_sites+op->num_odd_sites, &start_odd, &end_odd, l, threading, 1 );
  } else {
    compute_core_start_end_custom(0, l->num_inner_lattice_sites, &start, &n, l, threading, 1 );
  }
  
  if ( amount == _EVEN_SITES ) {
    start = start_odd, n = end_odd;
    minus_dir_param = _ODD_SITES;
    plus_dir_param = _EVEN_SITES;
  } else if ( amount == _ODD_SITES ) {
    start = start_even, n = end_even;
    minus_dir_param = _EVEN_SITES;
    plus_dir_param = _ODD_SITES;
  }
  
  complex_PRECISION *prn[4] = { op->prnT, op->prnZ, op->prnY, op->prnX };
  complex_PRECISION *prp[4] = { op->prpT, op->prpZ, op->prpY, op->prpX };
  
  // project minus dir
  prp_PRECISION( prn, phi, 12*start, 12*n );  
  
  // start communication in negative direction
  START_LOCKED_MASTER(threading)
  ghost_sendrecv_PRECISION( op->prnT, T, -1, &(op->c), minus_dir_param, l );
  ghost_sendrecv_PRECISION( op->prnZ, Z, -1, &(op->c), minus_dir_param, l );
  ghost_sendrecv_PRECISION( op->prnY, Y, -1, &(op->c), minus_dir_param, l );
  ghost_sendrecv_PRECISION( op->prnX, X, -1, &(op->c), minus_dir_param, l );
  END_LOCKED_MASTER(threading)
  
  // project plus dir and multiply with U dagger
  prn_su3_PRECISION( prp, phi, op, neighbor, 12*start, 12*n );

  if ( amount == _EVEN_SITES ) {
    start = start_even, n = end_even;
  } else if ( amount == _ODD_SITES ) {
    start = start_odd, n = end_odd;
  }  
  // start communication in positive direction
  START_LOCKED_MASTER(threading)
  ghost_sendrecv_PRECISION( op->prpT, T, +1, &(op->c), plus_dir_param, l );
  ghost_sendrecv_PRECISION( op->prpZ, Z, +1, &(op->c), plus_dir_param, l );
  ghost_sendrecv_PRECISION( op->prpY, Y, +1, &(op->c), plus_dir_param, l );
  ghost_sendrecv_PRECISION( op->prpX, X, +1, &(op->c), plus_dir_param, l );
  // wait for communication in negative direction
  ghost_wait_PRECISION( op->prnT, T, -1, &(op->c), minus_dir_param, l );
  ghost_wait_PRECISION( op->prnZ, Z, -1, &(op->c), minus_dir_param, l );
  ghost_wait_PRECISION( op->prnY, Y, -1, &(op->c), minus_dir_param, l );
  ghost_wait_PRECISION( op->prnX, X, -1, &(op->c), minus_dir_param, l );
  END_LOCKED_MASTER(threading) 
  
  // multiply with U and lift up minus dir
  su3_pbp_PRECISION( eta, prn, op, neighbor, 12*start, 12*n );
  
  // wait for communication in positive direction
  START_LOCKED_MASTER(threading)
  ghost_wait_PRECISION( op->prpT, T, +1, &(op->c), plus_dir_param, l );
  ghost_wait_PRECISION( op->prpZ, Z, +1, &(op->c), plus_dir_param, l );
  ghost_wait_PRECISION( op->prpY, Y, +1, &(op->c), plus_dir_param, l );
  ghost_wait_PRECISION( op->prpX, X, +1, &(op->c), plus_dir_param, l );
  END_LOCKED_MASTER(threading)
  
  // lift up plus dir
  pbn_PRECISION( eta, prp, 12*start, 12*n );

  SYNC_CORES(threading)
}
#endif

// ---- block odd even ---------------------------------------------------

#ifdef OPTIMIZED_SELF_COUPLING_PRECISION
void schwarz_PRECISION_oddeven_setup( operator_PRECISION_struct *op, level_struct *l ) {
  
  PRECISION *clover_pt = op->clover_vectorized, *oe_clover_pt = op->oe_clover_vectorized;
  int mu, i, d0, c0, b0, a0, d1, c1, b1, a1, t, z, y, x, agg_split[4], block_split[4], block_size[4];
  
  if ( g.csw ) {
    for ( mu=0; mu<4; mu++ ) {
      agg_split[mu] = l->local_lattice[mu]/l->coarsening[mu];
      block_split[mu] = l->coarsening[mu]/l->block_lattice[mu];
      block_size[mu] = l->block_lattice[mu];
    }
    
    for ( d0=0; d0<agg_split[T]; d0++ )
      for ( c0=0; c0<agg_split[Z]; c0++ )
        for ( b0=0; b0<agg_split[Y]; b0++ )
          for ( a0=0; a0<agg_split[X]; a0++ )
            
            for ( d1=d0*block_split[T]; d1<(d0+1)*block_split[T]; d1++ )
              for ( c1=c0*block_split[Z]; c1<(c0+1)*block_split[Z]; c1++ )
                for ( b1=b0*block_split[Y]; b1<(b0+1)*block_split[Y]; b1++ )
                  for ( a1=a0*block_split[X]; a1<(a0+1)*block_split[X]; a1++ ) {
                    
                    for ( t=d1*block_size[T]; t<(d1+1)*block_size[T]; t++ )
                      for ( z=c1*block_size[Z]; z<(c1+1)*block_size[Z]; z++ )
                        for ( y=b1*block_size[Y]; y<(b1+1)*block_size[Y]; y++ )
                          for ( x=a1*block_size[X]; x<(a1+1)*block_size[X]; x++ ) {
                            if (((t-d1*block_size[T])+(z-c1*block_size[Z])+
                                (y-b1*block_size[Y])+(x-a1*block_size[X]))%2 == 0 ) {
                              for ( i=0; i<144; i++ )
                                oe_clover_pt[i] = clover_pt[i];
                              clover_pt += 144;
                              oe_clover_pt += 144;
                            }
                          }
                    for ( t=d1*block_size[T]; t<(d1+1)*block_size[T]; t++ )
                      for ( z=c1*block_size[Z]; z<(c1+1)*block_size[Z]; z++ )
                        for ( y=b1*block_size[Y]; y<(b1+1)*block_size[Y]; y++ )
                          for ( x=a1*block_size[X]; x<(a1+1)*block_size[X]; x++ ) {
                            if (((t-d1*block_size[T])+(z-c1*block_size[Z])+
                                (y-b1*block_size[Y])+(x-a1*block_size[X]))%2 == 1 ) {
                              sse_site_clover_invert_PRECISION( clover_pt, oe_clover_pt );
                              clover_pt += 144;
                              oe_clover_pt += 144;
                            }
                          }
                  }
  } else {
    vector_PRECISION_copy( op->oe_clover, op->clover, 0, l->inner_vector_size, l );
  }
  op->shift = 4+l->dirac_shift;
}
#endif

#ifdef OPTIMIZED_SELF_COUPLING_PRECISION
void block_diag_ee_PRECISION( vector_PRECISION eta, vector_PRECISION phi,
     int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {
  
  START_UNTHREADED_FUNCTION(threading)
  PRECISION *clover_vectorized = s->op.oe_clover_vectorized + (start/12)*144;  
  int i, n1 = s->num_block_even_sites;
  config_PRECISION clover = (g.csw==0.0)?s->op.oe_clover+start:s->op.oe_clover+(start/12)*42;
  vector_PRECISION lphi = phi+start, leta = eta+start;
  // diagonal blocks applied to the even sites of a block
  if ( g.csw ) {
    for ( i=0; i<n1; i++ ) {
      sse_site_clover_PRECISION( (PRECISION*)leta, (PRECISION*)lphi, clover_vectorized );
      leta+=12; lphi+=12; clover_vectorized+=144;
    }
  } else {
    for ( i=0; i<12*n1; i++ )
      leta[i] = lphi[i]*clover[i];
  }

  END_UNTHREADED_FUNCTION(threading)
}
#endif

#ifdef OPTIMIZED_SELF_COUPLING_PRECISION
void block_diag_oo_PRECISION( vector_PRECISION eta, vector_PRECISION phi,
    int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {
  
  START_UNTHREADED_FUNCTION(threading)

  int i, n1 = s->num_block_even_sites, n2 = s->num_block_odd_sites;
  config_PRECISION clover = (g.csw==0.0)?s->op.oe_clover+start:s->op.oe_clover+(start/12)*42;
  vector_PRECISION lphi = phi+start, leta = eta+start;
  // diagonal blocks applied to the odd sites of a block
  if ( g.csw ) {
    error0("block_diag_oo_PRECISION is not available when using SSE\n");
  } else {
    leta += n1*12; lphi += n1*12; clover += n1*12;
    for ( i=0; i<12*n2; i++ )
      leta[i] = lphi[i]*clover[i];
  }

  END_UNTHREADED_FUNCTION(threading)
}
#endif

#ifdef OPTIMIZED_SELF_COUPLING_PRECISION
void block_diag_oo_inv_PRECISION( vector_PRECISION eta, vector_PRECISION phi,
    int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {

  START_UNTHREADED_FUNCTION(threading)
  PRECISION *clover_vectorized = s->op.oe_clover_vectorized + (start/12)*144;
  int i, n1 = s->num_block_even_sites, n2 = s->num_block_odd_sites;
  config_PRECISION clover = (g.csw==0.0)?s->op.oe_clover+start:s->op.oe_clover+(start/12)*42;
  vector_PRECISION lphi = phi+start, leta = eta+start;
  // inverted diagonal blocks applied to the odd sites of a block
  if ( g.csw ) {
    leta += n1*12; lphi += n1*12; clover_vectorized += n1*144;
    for ( i=0; i<n2; i++ ) {
      sse_site_clover_PRECISION( (PRECISION*)leta, (PRECISION*)lphi, clover_vectorized );
      leta+=12; lphi+=12; clover_vectorized+=144;
    }
  } else {
    leta += n1*12; lphi += n1*12; clover += n1*12;
    for ( i=0; i<12*n2; i++ )
      leta[i] = lphi[i]/clover[i];
  }

  END_UNTHREADED_FUNCTION(threading)
}
#endif


#ifdef OPTIMIZED_NEIGHBOR_COUPLING_PRECISION
void block_hopping_term_PRECISION( vector_PRECISION eta, vector_PRECISION phi, 
                                   int start, int amount, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {
  
  START_UNTHREADED_FUNCTION(threading)
  
  int *length_even = s->dir_length_even, *length_odd = s->dir_length_odd,
      **index = s->oe_index, *neighbor = s->op.neighbor_table;
  PRECISION *Dplus = s->op.D_vectorized + (start/12)*96;
  PRECISION *Dminus = s->op.D_transformed_vectorized + (start/12)*96;

  for ( int mu=0; mu<4; mu++ ) {
    int a1, a2, n1, n2;
    if ( amount == _EVEN_SITES ) {
      a1 = 0; n1 = length_even[mu];
      a2 = n1; n2 = a2 + length_odd[mu];
    } else if ( amount == _ODD_SITES ) {
      a1 = length_even[mu]; n1 = a1 + length_odd[mu];
      a2 = 0; n2 = a1;
    } else {
      a1 = 0; n1 = length_even[mu]+length_odd[mu];
      a2 = 0; n2 = n1;
    }
    block_oddeven_plus_coupling_PRECISION( (PRECISION*)(eta+start), Dplus, (PRECISION*)(phi+start),
                                           mu, a1, n1, index[mu], neighbor );
    block_oddeven_minus_coupling_PRECISION( (PRECISION*)(eta+start), Dminus, (PRECISION*)(phi+start),
                                            mu, a2, n2, index[mu], neighbor );
  }
  
  END_UNTHREADED_FUNCTION(threading)
}
#endif

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_PRECISION
void block_n_hopping_term_PRECISION( vector_PRECISION eta, vector_PRECISION phi, 
                                   int start, int amount, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {
  
  START_UNTHREADED_FUNCTION(threading)
  
  int *length_even = s->dir_length_even, *length_odd = s->dir_length_odd,
      **index = s->oe_index, *neighbor = s->op.neighbor_table;
  PRECISION *Dplus = s->op.D_vectorized + (start/12)*96;
  PRECISION *Dminus = s->op.D_transformed_vectorized + (start/12)*96;

  for ( int mu=0; mu<4; mu++ ) {
    int a1, a2, n1, n2;
    if ( amount == _EVEN_SITES ) {
      a1 = 0; n1 = length_even[mu];
      a2 = n1; n2 = a2 + length_odd[mu];
    } else if ( amount == _ODD_SITES ) {
      a1 = length_even[mu]; n1 = a1 + length_odd[mu];
      a2 = 0; n2 = a1;
    } else {
      a1 = 0; n1 = length_even[mu]+length_odd[mu];
      a2 = 0; n2 = n1;
    }
    block_oddeven_nplus_coupling_PRECISION( (PRECISION*)(eta+start), Dplus, (PRECISION*)(phi+start),
                                            mu, a1, n1, index[mu], neighbor );
    block_oddeven_nminus_coupling_PRECISION( (PRECISION*)(eta+start), Dminus, (PRECISION*)(phi+start),
                                             mu, a2, n2, index[mu], neighbor );
  }
  
  END_UNTHREADED_FUNCTION(threading)
}
#endif


#endif // SSE

