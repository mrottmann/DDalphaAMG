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

#include "main.h"

void clover_PRECISION( vector_PRECISION eta, vector_PRECISION phi, operator_PRECISION_struct *op, int start, int end,
                       level_struct *l, struct Thread *threading ) {

  int nv = l->num_lattice_site_var;
  vector_PRECISION lphi = phi+start, leta = eta+start;
  vector_PRECISION leta_end = eta+end;

#ifdef PROFILING
  START_MASTER(threading)
  PROF_PRECISION_START( _SC );
  END_MASTER(threading)
#endif

#ifdef HAVE_TM
  config_PRECISION tm_term = op->tm_term+(start/nv)*12;
#endif

  if ( g.csw == 0.0 ) {

    config_PRECISION clover = op->clover+(start/nv)*12;
#ifdef HAVE_TM1p1
    if( g.n_flavours == 2 ) {
#ifdef HAVE_TM
      if ( g.mu + g.mu_odd_shift != 0.0 || g.mu + g.mu_even_shift != 0.0 ) 
        while ( leta < leta_end ) {
          FOR6( *leta = (*lphi)*((*clover)+(*tm_term)); leta++; lphi++; clover++; tm_term++; );
          clover -= 6;
          tm_term -= 6;
          FOR6( *leta = (*lphi)*((*clover)-(*tm_term)); leta++; lphi++; clover++; tm_term++; );
          FOR6( *leta = (*lphi)*((*clover)+(*tm_term)); leta++; lphi++; clover++; tm_term++; );
          clover -= 6;
          tm_term -= 6;
          FOR6( *leta = (*lphi)*((*clover)-(*tm_term)); leta++; lphi++; clover++; tm_term++; );
        }
      else
#endif
        while ( leta < leta_end ) {
          FOR6( *leta = (*lphi)*(*clover); leta++; lphi++; clover++; );
          clover -= 6;
          FOR12( *leta = (*lphi)*(*clover); leta++; lphi++; clover++; );
          clover -= 6;
          FOR6( *leta = (*lphi)*(*clover); leta++; lphi++; clover++; );
        }
    } else {
#endif
#ifdef HAVE_TM
      if ( g.mu + g.mu_odd_shift != 0.0 || g.mu + g.mu_even_shift != 0.0 ) {
        while ( leta < leta_end )
          FOR12( *leta = (*lphi)*((*clover)+(*tm_term)); leta++; lphi++; clover++; tm_term++; );
      } else
#endif
        while ( leta < leta_end )
          FOR12( *leta = (*lphi)*(*clover); leta++; lphi++; clover++; );
#ifdef HAVE_TM1p1
    }
#endif

  } else {

#ifndef OPTIMIZED_SELF_COUPLING_PRECISION

    config_PRECISION clover = op->clover+(start/nv)*42;
#ifdef HAVE_TM1p1
    if( g.n_flavours == 2 ) {
#ifdef HAVE_TM
      if ( g.mu + g.mu_odd_shift != 0.0 || g.mu + g.mu_even_shift != 0.0 ) 
        while ( leta < leta_end ) {
          doublet_site_clover_PRECISION( leta, lphi, clover );
          clover+=42;
          FOR6( *leta += (*lphi)*(*tm_term); leta++; lphi++; tm_term++; );
          tm_term -= 6;
          FOR6( *leta -=(*lphi)*(*tm_term); leta++; lphi++; tm_term++; );
          FOR6( *leta += (*lphi)*(*tm_term); leta++; lphi++; tm_term++; );
          tm_term -= 6;
          FOR6( *leta -= (*lphi)*(*tm_term); leta++; lphi++; tm_term++; );
        }
      else
#endif
        while ( leta < leta_end ) {
          doublet_site_clover_PRECISION( leta, lphi, clover );
          leta+=24; lphi+=24;
          clover+=42;
        }
    } else {
#endif
#ifdef HAVE_TM
      if ( g.mu + g.mu_odd_shift != 0.0 || g.mu + g.mu_even_shift != 0.0 ) 
        while ( leta < leta_end ) {
          site_clover_PRECISION( leta, lphi, clover );
          FOR12( *leta += (*lphi)*(*tm_term); leta++; lphi++; tm_term++; );
          clover+=42;
        }
      else
#endif
        while ( leta < leta_end ) {
          site_clover_PRECISION( leta, lphi, clover );
          leta+=12; lphi+=12;
          clover+=42;
        }
#ifdef HAVE_TM1p1
    }
#endif

#else

#ifdef HAVE_TM1p1
    PRECISION *clover = ( g.n_flavours == 2 ) ? op->clover_doublet_vectorized : op->clover_vectorized;
#else
    PRECISION *clover = op->clover_vectorized;
#endif
    clover += start*12;
    while ( leta < leta_end ) { // tm_term included in the clover vectorized
      sse_site_clover_PRECISION( (PRECISION*) leta, (PRECISION*) lphi, clover );
      leta += nv; lphi += nv;
      clover += 12*nv;
    }
    
#endif
    
  }

#ifdef HAVE_TM1p1
  config_PRECISION eps_term = op->epsbar_term+(start/nv)*12;  
  lphi = phi+start, leta = eta+start;
  if ( g.n_flavours == 2 &&
       ( g.epsbar != 0 || g.epsbar_ig5_odd_shift != 0 || g.epsbar_ig5_odd_shift != 0 ) )
    while ( leta < leta_end ) { 
      lphi += 6;
      FOR6( *leta += (*lphi)*(*eps_term); leta++; lphi++; eps_term++; )
      lphi -= 12;
      eps_term -= 6;
      FOR6( *leta += (*lphi)*(*eps_term); leta++; lphi++; eps_term++; )
      lphi += 6;
    }
#endif

  
#ifdef PROFILING
  START_MASTER(threading)
  PROF_PRECISION_STOP( _SC, 1 );
  END_MASTER(threading)
#endif
    
}

static void spin0and1_clover_PRECISION( vector_PRECISION eta, vector_PRECISION phi, config_PRECISION clover, level_struct *l ) {
  
  vector_PRECISION eta_end = eta + l->inner_vector_size;
  if ( g.csw == 0.0 ) {
    while ( eta < eta_end ) {
      FOR6( *eta = (*phi)*(*clover); eta++; phi++; clover++; )
      FOR6( *eta = _COMPLEX_PRECISION_ZERO; eta++; )
      phi+=6; clover+=6;
    }
  } else {
    while ( eta < eta_end ) {
      spin0and1_site_clover_PRECISION( eta, phi, clover );
      eta+=12; phi+=12; clover+=42;
    }
  }
}


static void spin2and3_clover_PRECISION( vector_PRECISION eta, vector_PRECISION phi, config_PRECISION clover, level_struct *l ) {
  
  vector_PRECISION eta_end = eta + l->inner_vector_size;
  if ( g.csw == 0.0 ) {
    while ( eta < eta_end ) {
      phi+=6; clover+=6;
      FOR6( *eta = _COMPLEX_PRECISION_ZERO; eta++; )
      FOR6( *eta = (*phi)*(*clover); eta++; phi++; clover++; )
    }
  } else {
    while ( eta < eta_end ) {
      spin2and3_site_clover_PRECISION( eta, phi, clover );
      eta +=12; phi+=12; clover+=42;
    }
  }
}

void block_d_plus_clover_PRECISION( vector_PRECISION eta, vector_PRECISION phi, int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {

  START_UNTHREADED_FUNCTION(threading)
  
  int n = s->num_block_sites, *length = s->dir_length, **index = s->index, *neighbor = s->op.neighbor_table, nv = l->num_lattice_site_var;
  vector_PRECISION lphi = phi+start, leta = eta+start;

  // clover term
  clover_PRECISION(eta, phi, &(s->op), start, start+nv*n, l, no_threading );

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_float // block operator vectorized just in the float environment
  PRECISION *Dplus = s->op.D_vectorized + (start/nv)*96;
  PRECISION *Dminus = s->op.D_transformed_vectorized + (start/nv)*96;
  for ( int mu=0; mu<4; mu++ ) {
    block_oddeven_plus_coupling_PRECISION( (PRECISION*)leta, Dplus, (PRECISION*)lphi, mu, 0, length[mu], index[mu], neighbor );
    block_oddeven_minus_coupling_PRECISION( (PRECISION*)leta, Dminus, (PRECISION*)lphi, mu, 0, length[mu], index[mu], neighbor );
  }
#else
  int i, j, k, *ind;
  config_PRECISION D_pt;
  config_PRECISION D = s->op.D + (start/nv)*36;
#ifdef HAVE_TM1p1
  if( g.n_flavours == 2 ) {
    complex_PRECISION buf1[50]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, *buf2=buf1+12, *buf3=buf2+12, *buf4=buf3+12;
    // inner block couplings
    ind = index[T]; // T direction
    for ( i=0; i<length[T]; i++ ) {
      k = ind[i]; j = neighbor[4*k+T]; D_pt = D + 36*k + 9*T;
      dprn_T_PRECISION( buf1, lphi+24*k ); // (1+gamma_T) phi(x) + projection
      dprp_T_PRECISION( buf2, lphi+24*j ); // (1-gamma_T) phi(x+hat{T}) + projection
      mvmh_PRECISION( buf3, D_pt, buf1 );     // U_T^dagger(x) (1+gamma_T) phi(x) - projected
      mvmh_PRECISION( buf3+3, D_pt, buf1+3 ); // U_T^dagger(x) (1+gamma_T) phi(x) - projected
      mvmh_PRECISION( buf3+6, D_pt, buf1+6 ); // U_T^dagger(x) (1+gamma_T) phi(x) - projected
      mvmh_PRECISION( buf3+9, D_pt, buf1+9 ); // U_T^dagger(x) (1+gamma_T) phi(x) - projected
      mvm_PRECISION( buf4, D_pt, buf2 );      // U_T(x) (1-gamma_T) phi(x+hat{T}) - projected
      mvm_PRECISION( buf4+3, D_pt, buf2+3 );     // U_T(x) (1-gamma_T) phi(x+hat{T}) - projected
      mvm_PRECISION( buf4+6, D_pt, buf2+6 ); // U_T(x) (1-gamma_T) phi(x+hat{T}) - projected
      mvm_PRECISION( buf4+9, D_pt, buf2+9 ); // U_T(x) (1-gamma_T) phi(x+hat{T}) - projected
      dpbn_su3_T_PRECISION( buf3, leta+24*j ); // eta(x+hat{T}) -= U_T(x)^dagger(x) (1+gamma_T) phi(x) + lift back
      dpbp_su3_T_PRECISION( buf4, leta+24*k ); // eta(x) -= U_T(x) (1-gamma_T) phi(x+hat{T}) + lift back
    }
    ind = index[Z]; // Z direction
    for ( i=0; i<length[Z]; i++ ) {
      k = ind[i]; j = neighbor[4*k+Z]; D_pt = D + 36*k + 9*Z;
      dprn_Z_PRECISION( buf1, lphi+24*k ); // (1+gamma_Z) phi(x) + projection
      dprp_Z_PRECISION( buf2, lphi+24*j ); // (1-gamma_Z) phi(x+hat{Z}) + projection
      mvmh_PRECISION( buf3, D_pt, buf1 );     // U_Z^dagger(x) (1+gamma_Z) phi(x) - projected
      mvmh_PRECISION( buf3+3, D_pt, buf1+3 );     // U_Z^dagger(x) (1+gamma_Z) phi(x) - projected
      mvmh_PRECISION( buf3+6, D_pt, buf1+6 ); // U_Z^dagger(x) (1+gamma_Z) phi(x) - projected
      mvmh_PRECISION( buf3+9, D_pt, buf1+9 ); // U_Z^dagger(x) (1+gamma_Z) phi(x) - projected
      mvm_PRECISION( buf4, D_pt, buf2 );     // U_Z(x) (1-gamma_Z) phi(x+hat{Z}) - projected
      mvm_PRECISION( buf4+3, D_pt, buf2+3 );     // U_Z(x) (1-gamma_Z) phi(x+hat{Z}) - projected
      mvm_PRECISION( buf4+6, D_pt, buf2+6 ); // U_Z(x) (1-gamma_Z) phi(x+hat{Z}) - projected
      mvm_PRECISION( buf4+9, D_pt, buf2+9 ); // U_Z(x) (1-gamma_Z) phi(x+hat{Z}) - projected
      dpbn_su3_Z_PRECISION( buf3, leta+24*j ); // eta(x+hat{Z}) -= U_Z(x)^dagger(x) (1+gamma_Z) phi(x) + lift back
      dpbp_su3_Z_PRECISION( buf4, leta+24*k ); // eta(x) -= U_Z(x) (1-gamma_Z) phi(x+hat{Z}) + lift back
    }
    ind = index[Y]; // Y direction
    for ( i=0; i<length[Y]; i++ ) {
      k = ind[i]; j = neighbor[4*k+Y]; D_pt = D + 36*k + 9*Y;
      dprn_Y_PRECISION( buf1, lphi+24*k ); // (1+gamma_Y) phi(x) + projection
      dprp_Y_PRECISION( buf2, lphi+24*j ); // (1-gamma_Y) phi(x+hat{Y}) + projection
      mvmh_PRECISION( buf3, D_pt, buf1 );     // U_Y^dagger(x) (1+gamma_Y) phi(x) - projected
      mvmh_PRECISION( buf3+3, D_pt, buf1+3 );     // U_Y^dagger(x) (1+gamma_Y) phi(x) - projected
      mvmh_PRECISION( buf3+6, D_pt, buf1+6 ); // U_Y^dagger(x) (1+gamma_Y) phi(x) - projected
      mvmh_PRECISION( buf3+9, D_pt, buf1+9 ); // U_Y^dagger(x) (1+gamma_Y) phi(x) - projected
      mvm_PRECISION( buf4, D_pt, buf2 );     // U_Y(x) (1-gamma_Y) phi(x+hat{Y}) - projected
      mvm_PRECISION( buf4+3, D_pt, buf2+3 );     // U_Y(x) (1-gamma_Y) phi(x+hat{Y}) - projected
      mvm_PRECISION( buf4+6, D_pt, buf2+6 ); // U_Y(x) (1-gamma_Y) phi(x+hat{Y}) - projected
      mvm_PRECISION( buf4+9, D_pt, buf2+9 ); // U_Y(x) (1-gamma_Y) phi(x+hat{Y}) - projected
      dpbn_su3_Y_PRECISION( buf3, leta+24*j ); // eta(x+hat{Y}) -= U_Y(x)^dagger(x) (1+gamma_Y) phi(x) + lift back
      dpbp_su3_Y_PRECISION( buf4, leta+24*k ); // eta(x) -= U_Y(x) (1-gamma_Y) phi(x+hat{Y}) + lift back
    }
    ind = index[X]; // X direction
    for ( i=0; i<length[X]; i++ ) {
      k = ind[i]; j = neighbor[4*k+X]; D_pt = D + 36*k + 9*X;
      dprn_X_PRECISION( buf1, lphi+24*k ); // (1+gamma_X) phi(x) + projection
      dprp_X_PRECISION( buf2, lphi+24*j ); // (1-gamma_X) phi(x+hat{X}) + projection
      mvmh_PRECISION( buf3, D_pt, buf1 );     // U_X^dagger(x) (1+gamma_X) phi(x) - projected
      mvmh_PRECISION( buf3+3, D_pt, buf1+3 );     // U_X^dagger(x) (1+gamma_X) phi(x) - projected
      mvmh_PRECISION( buf3+6, D_pt, buf1+6 ); // U_X^dagger(x) (1+gamma_X) phi(x) - projected
      mvmh_PRECISION( buf3+9, D_pt, buf1+9 ); // U_X^dagger(x) (1+gamma_X) phi(x) - projected
      mvm_PRECISION( buf4, D_pt, buf2 );     // U_mu(x) (1-gamma_X) phi(x+hat{X}) - projected
      mvm_PRECISION( buf4+3, D_pt, buf2+3 );     // U_mu(x) (1-gamma_X) phi(x+hat{X}) - projected
      mvm_PRECISION( buf4+6, D_pt, buf2+6 ); // U_mu(x) (1-gamma_X) phi(x+hat{X}) - projected
      mvm_PRECISION( buf4+9, D_pt, buf2+9 ); // U_mu(x) (1-gamma_X) phi(x+hat{X}) - projected
      dpbn_su3_X_PRECISION( buf3, leta+24*j ); // eta(x+hat{X}) -= U_X(x)^dagger(x) (1+gamma_X) phi(x) + lift back
      dpbp_su3_X_PRECISION( buf4, leta+24*k ); // eta(x) -= U_X(x) (1-gamma_X) phi(x+hat{X}) + lift back
    }    
  } else {
#endif   
    complex_PRECISION buf1[25]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, *buf2=buf1+6, *buf3=buf2+6, *buf4=buf3+6;
    // inner block couplings
    ind = index[T]; // T direction
    for ( i=0; i<length[T]; i++ ) {
      k = ind[i]; j = neighbor[4*k+T]; D_pt = D + 36*k + 9*T;
      prn_T_PRECISION( buf1, lphi+12*k ); // (1+gamma_T) phi(x) + projection
      prp_T_PRECISION( buf2, lphi+12*j ); // (1-gamma_T) phi(x+hat{T}) + projection
      mvmh_PRECISION( buf3, D_pt, buf1 );     // U_T^dagger(x) (1+gamma_T) phi(x) - projected
      mvmh_PRECISION( buf3+3, D_pt, buf1+3 ); // U_T^dagger(x) (1+gamma_T) phi(x) - projected
      mvm_PRECISION( buf4, D_pt, buf2 );     // U_T(x) (1-gamma_T) phi(x+hat{T}) - projected
      mvm_PRECISION( buf4+3, D_pt, buf2+3 ); // U_T(x) (1-gamma_T) phi(x+hat{T}) - projected
      pbn_su3_T_PRECISION( buf3, leta+12*j ); // eta(x+hat{T}) -= U_T(x)^dagger(x) (1+gamma_T) phi(x) + lift back
      pbp_su3_T_PRECISION( buf4, leta+12*k ); // eta(x) -= U_T(x) (1-gamma_T) phi(x+hat{T}) + lift back
    }
    ind = index[Z]; // Z direction
    for ( i=0; i<length[Z]; i++ ) {
      k = ind[i]; j = neighbor[4*k+Z]; D_pt = D + 36*k + 9*Z;
      prn_Z_PRECISION( buf1, lphi+12*k ); // (1+gamma_Z) phi(x) + projection
      prp_Z_PRECISION( buf2, lphi+12*j ); // (1-gamma_Z) phi(x+hat{Z}) + projection
      mvmh_PRECISION( buf3, D_pt, buf1 );     // U_Z^dagger(x) (1+gamma_Z) phi(x) - projected
      mvmh_PRECISION( buf3+3, D_pt, buf1+3 ); // U_Z^dagger(x) (1+gamma_Z) phi(x) - projected
      mvm_PRECISION( buf4, D_pt, buf2 );     // U_Z(x) (1-gamma_Z) phi(x+hat{Z}) - projected
      mvm_PRECISION( buf4+3, D_pt, buf2+3 ); // U_Z(x) (1-gamma_Z) phi(x+hat{Z}) - projected
      pbn_su3_Z_PRECISION( buf3, leta+12*j ); // eta(x+hat{Z}) -= U_Z(x)^dagger(x) (1+gamma_Z) phi(x) + lift back
      pbp_su3_Z_PRECISION( buf4, leta+12*k ); // eta(x) -= U_Z(x) (1-gamma_Z) phi(x+hat{Z}) + lift back
    }
    ind = index[Y]; // Y direction
    for ( i=0; i<length[Y]; i++ ) {
      k = ind[i]; j = neighbor[4*k+Y]; D_pt = D + 36*k + 9*Y;
      prn_Y_PRECISION( buf1, lphi+12*k ); // (1+gamma_Y) phi(x) + projection
      prp_Y_PRECISION( buf2, lphi+12*j ); // (1-gamma_Y) phi(x+hat{Y}) + projection
      mvmh_PRECISION( buf3, D_pt, buf1 );     // U_Y^dagger(x) (1+gamma_Y) phi(x) - projected
      mvmh_PRECISION( buf3+3, D_pt, buf1+3 ); // U_Y^dagger(x) (1+gamma_Y) phi(x) - projected
      mvm_PRECISION( buf4, D_pt, buf2 );     // U_Y(x) (1-gamma_Y) phi(x+hat{Y}) - projected
      mvm_PRECISION( buf4+3, D_pt, buf2+3 ); // U_Y(x) (1-gamma_Y) phi(x+hat{Y}) - projected
      pbn_su3_Y_PRECISION( buf3, leta+12*j ); // eta(x+hat{Y}) -= U_Y(x)^dagger(x) (1+gamma_Y) phi(x) + lift back
      pbp_su3_Y_PRECISION( buf4, leta+12*k ); // eta(x) -= U_Y(x) (1-gamma_Y) phi(x+hat{Y}) + lift back
    }
    ind = index[X]; // X direction
    for ( i=0; i<length[X]; i++ ) {
      k = ind[i]; j = neighbor[4*k+X]; D_pt = D + 36*k + 9*X;
      prn_X_PRECISION( buf1, lphi+12*k ); // (1+gamma_X) phi(x) + projection
      prp_X_PRECISION( buf2, lphi+12*j ); // (1-gamma_X) phi(x+hat{X}) + projection
      mvmh_PRECISION( buf3, D_pt, buf1 );     // U_X^dagger(x) (1+gamma_X) phi(x) - projected
      mvmh_PRECISION( buf3+3, D_pt, buf1+3 ); // U_X^dagger(x) (1+gamma_X) phi(x) - projected
      mvm_PRECISION( buf4, D_pt, buf2 );     // U_mu(x) (1-gamma_X) phi(x+hat{X}) - projected
      mvm_PRECISION( buf4+3, D_pt, buf2+3 ); // U_mu(x) (1-gamma_X) phi(x+hat{X}) - projected
      pbn_su3_X_PRECISION( buf3, leta+12*j ); // eta(x+hat{X}) -= U_X(x)^dagger(x) (1+gamma_X) phi(x) + lift back
      pbp_su3_X_PRECISION( buf4, leta+12*k ); // eta(x) -= U_X(x) (1-gamma_X) phi(x+hat{X}) + lift back
    }
#ifdef HAVE_TM1p1
  }
#endif
#endif
  END_UNTHREADED_FUNCTION(threading)
}

void d_plus_clover_PRECISION( vector_PRECISION eta, vector_PRECISION phi, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {
  
  int n = l->num_inner_lattice_sites, *neighbor = op->neighbor_table, start, end, nv = l->num_lattice_site_var;
#ifdef OPTIMIZED_NEIGHBOR_COUPLING_PRECISION
  complex_PRECISION *prn[4] = { op->prnT, op->prnZ, op->prnY, op->prnX };
  complex_PRECISION *prp[4] = { op->prpT, op->prpZ, op->prpY, op->prpX };
#else
  int i, j, *nb_pt;
  vector_PRECISION phi_pt, eta_pt, end_pt;
  config_PRECISION D_pt;
#endif

  compute_core_start_end(0, nv*n, &start, &end, l, threading );

  SYNC_MASTER_TO_ALL(threading)

  clover_PRECISION(eta, phi, op, start, end, l, threading );

  START_MASTER(threading)
  PROF_PRECISION_START( _NC ); 
  END_MASTER(threading)

#ifdef HAVE_TM1p1
  if( g.n_flavours == 2 ) {
#ifdef OPTIMIZED_NEIGHBOR_COUPLING_PRECISION
    dprp_PRECISION( prn, phi, start, end );
#else
    complex_PRECISION pbuf[12];  
    for ( i=start/2, phi_pt=phi+start; i<end/2; i+=12, phi_pt+=24 ) {
      dprp_T_PRECISION( op->prnT+i, phi_pt );
      dprp_Z_PRECISION( op->prnZ+i, phi_pt );
      dprp_Y_PRECISION( op->prnY+i, phi_pt );
      dprp_X_PRECISION( op->prnX+i, phi_pt );
    }
#endif
    // start communication in negative direction
    START_LOCKED_MASTER(threading)
    ghost_sendrecv_PRECISION( op->prnT, T, -1, &(op->c), _FULL_SYSTEM, l );
    ghost_sendrecv_PRECISION( op->prnZ, Z, -1, &(op->c), _FULL_SYSTEM, l );
    ghost_sendrecv_PRECISION( op->prnY, Y, -1, &(op->c), _FULL_SYSTEM, l );
    ghost_sendrecv_PRECISION( op->prnX, X, -1, &(op->c), _FULL_SYSTEM, l );
    END_LOCKED_MASTER(threading) 
  
#ifdef OPTIMIZED_NEIGHBOR_COUPLING_PRECISION
    dprn_su3_PRECISION( prp, phi, op, neighbor, start, end );
#else
    // project plus dir and multiply with U dagger
    for ( phi_pt=phi+start, end_pt=phi+end, D_pt = op->D+((start/nv)*36), nb_pt=neighbor+((start/nv)*4); phi_pt<end_pt; phi_pt+=24 ) {
      // T dir
      j = 12*(*nb_pt); nb_pt++;
      dprn_T_PRECISION( pbuf, phi_pt );
      mvmh_PRECISION( op->prpT+j, D_pt, pbuf );
      mvmh_PRECISION( op->prpT+j+3, D_pt, pbuf+3 );
      mvmh_PRECISION( op->prpT+j+6, D_pt, pbuf+6 );
      mvmh_PRECISION( op->prpT+j+9, D_pt, pbuf+9 ); D_pt += 9;
      // Z dir
      j = 12*(*nb_pt); nb_pt++;
      dprn_Z_PRECISION( pbuf, phi_pt );
      mvmh_PRECISION( op->prpZ+j, D_pt, pbuf );
      mvmh_PRECISION( op->prpZ+j+3, D_pt, pbuf+3 );
      mvmh_PRECISION( op->prpZ+j+6, D_pt, pbuf+6 );
      mvmh_PRECISION( op->prpZ+j+9, D_pt, pbuf+9 ); D_pt += 9;
      // Y dir
      j = 12*(*nb_pt); nb_pt++;
      dprn_Y_PRECISION( pbuf, phi_pt );
      mvmh_PRECISION( op->prpY+j, D_pt, pbuf );
      mvmh_PRECISION( op->prpY+j+3, D_pt, pbuf+3 );
      mvmh_PRECISION( op->prpY+j+6, D_pt, pbuf+6 );
      mvmh_PRECISION( op->prpY+j+9, D_pt, pbuf+9 ); D_pt += 9;
      // X dir
      j = 12*(*nb_pt); nb_pt++;
      dprn_X_PRECISION( pbuf, phi_pt );
      mvmh_PRECISION( op->prpX+j, D_pt, pbuf );
      mvmh_PRECISION( op->prpX+j+3, D_pt, pbuf+3 );
      mvmh_PRECISION( op->prpX+j+6, D_pt, pbuf+6 );
      mvmh_PRECISION( op->prpX+j+9, D_pt, pbuf+9 ); D_pt += 9;
    }
#endif

    // start communication in positive direction
    START_LOCKED_MASTER(threading)
    ghost_sendrecv_PRECISION( op->prpT, T, +1, &(op->c), _FULL_SYSTEM, l );
    ghost_sendrecv_PRECISION( op->prpZ, Z, +1, &(op->c), _FULL_SYSTEM, l );
    ghost_sendrecv_PRECISION( op->prpY, Y, +1, &(op->c), _FULL_SYSTEM, l );
    ghost_sendrecv_PRECISION( op->prpX, X, +1, &(op->c), _FULL_SYSTEM, l );
    // wait for communication in negative direction
    ghost_wait_PRECISION( op->prnT, T, -1, &(op->c), _FULL_SYSTEM, l );
    ghost_wait_PRECISION( op->prnZ, Z, -1, &(op->c), _FULL_SYSTEM, l );
    ghost_wait_PRECISION( op->prnY, Y, -1, &(op->c), _FULL_SYSTEM, l );
    ghost_wait_PRECISION( op->prnX, X, -1, &(op->c), _FULL_SYSTEM, l );
    END_LOCKED_MASTER(threading)
     
#ifdef OPTIMIZED_NEIGHBOR_COUPLING_PRECISION
    su3_dpbp_PRECISION( eta, prn, op, neighbor, start, end );
#else 
    // multiply with U and lift up minus dir
    for ( eta_pt=eta+start, end_pt=eta+end, D_pt = op->D+(start/nv)*36, nb_pt=neighbor+(start/nv)*4; eta_pt<end_pt; eta_pt+=24 ) {
      // T dir
      j = 12*(*nb_pt); nb_pt++;
      mvm_PRECISION( pbuf, D_pt, op->prnT+j );
      mvm_PRECISION( pbuf+3, D_pt, op->prnT+j+3 );
      mvm_PRECISION( pbuf+6, D_pt, op->prnT+j+6 );
      mvm_PRECISION( pbuf+9, D_pt, op->prnT+j+9 );
      dpbp_su3_T_PRECISION( pbuf, eta_pt ); D_pt += 9;
      // Z dir
      j = 12*(*nb_pt); nb_pt++;
      mvm_PRECISION( pbuf, D_pt, op->prnZ+j );
      mvm_PRECISION( pbuf+3, D_pt, op->prnZ+j+3 );
      mvm_PRECISION( pbuf+6, D_pt, op->prnZ+j+6 );
      mvm_PRECISION( pbuf+9, D_pt, op->prnZ+j+9 );
      dpbp_su3_Z_PRECISION( pbuf, eta_pt ); D_pt += 9;
      // Y dir
      j = 12*(*nb_pt); nb_pt++;
      mvm_PRECISION( pbuf, D_pt, op->prnY+j );
      mvm_PRECISION( pbuf+3, D_pt, op->prnY+j+3 );
      mvm_PRECISION( pbuf+6, D_pt, op->prnY+j+6 );
      mvm_PRECISION( pbuf+9, D_pt, op->prnY+j+9 );
      dpbp_su3_Y_PRECISION( pbuf, eta_pt ); D_pt += 9;
      // X dir
      j = 12*(*nb_pt); nb_pt++;
      mvm_PRECISION( pbuf, D_pt, op->prnX+j );
      mvm_PRECISION( pbuf+3, D_pt, op->prnX+j+3 );
      mvm_PRECISION( pbuf+6, D_pt, op->prnX+j+6 );
      mvm_PRECISION( pbuf+9, D_pt, op->prnX+j+9 );
      dpbp_su3_X_PRECISION( pbuf, eta_pt ); D_pt += 9;
    }
#endif

    // wait for communication in positive direction
    START_LOCKED_MASTER(threading)
    ghost_wait_PRECISION( op->prpT, T, +1, &(op->c), _FULL_SYSTEM, l );
    ghost_wait_PRECISION( op->prpZ, Z, +1, &(op->c), _FULL_SYSTEM, l );
    ghost_wait_PRECISION( op->prpY, Y, +1, &(op->c), _FULL_SYSTEM, l );
    ghost_wait_PRECISION( op->prpX, X, +1, &(op->c), _FULL_SYSTEM, l );
    END_LOCKED_MASTER(threading)
      
    // lift up plus dir
#ifdef OPTIMIZED_NEIGHBOR_COUPLING_PRECISION
    dpbn_PRECISION( eta, prp, start, end );
#else
    for ( i=start/2, eta_pt=eta+start; i<end/2; i+=12, eta_pt+=24 ) {
      dpbn_su3_T_PRECISION( op->prpT+i, eta_pt );
      dpbn_su3_Z_PRECISION( op->prpZ+i, eta_pt );
      dpbn_su3_Y_PRECISION( op->prpY+i, eta_pt );
      dpbn_su3_X_PRECISION( op->prpX+i, eta_pt );
    }
#endif
  } else {
#endif

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_PRECISION
  prp_PRECISION( prn, phi, start, end );
#else
  complex_PRECISION pbuf[6];
  for ( i=start/2, phi_pt=phi+start; i<end/2; i+=6, phi_pt+=12 ) {
    prp_T_PRECISION( op->prnT+i, phi_pt );
    prp_Z_PRECISION( op->prnZ+i, phi_pt );
    prp_Y_PRECISION( op->prnY+i, phi_pt );
    prp_X_PRECISION( op->prnX+i, phi_pt );
  }
#endif
  // start communication in negative direction
  START_LOCKED_MASTER(threading)
  ghost_sendrecv_PRECISION( op->prnT, T, -1, &(op->c), _FULL_SYSTEM, l );
  ghost_sendrecv_PRECISION( op->prnZ, Z, -1, &(op->c), _FULL_SYSTEM, l );
  ghost_sendrecv_PRECISION( op->prnY, Y, -1, &(op->c), _FULL_SYSTEM, l );
  ghost_sendrecv_PRECISION( op->prnX, X, -1, &(op->c), _FULL_SYSTEM, l );
  END_LOCKED_MASTER(threading) 
  
  // project plus dir and multiply with U dagger
#ifdef OPTIMIZED_NEIGHBOR_COUPLING_PRECISION
  prn_su3_PRECISION( prp, phi, op, neighbor, start, end );
#else
  for ( phi_pt=phi+start, end_pt=phi+end, D_pt = op->D+(start*3), nb_pt=neighbor+((start/12)*4); phi_pt<end_pt; phi_pt+=12 ) {
    // T dir
    j = 6*(*nb_pt); nb_pt++;
    prn_T_PRECISION( pbuf, phi_pt );
    mvmh_PRECISION( op->prpT+j, D_pt, pbuf );
    mvmh_PRECISION( op->prpT+j+3, D_pt, pbuf+3 ); D_pt += 9;
    // Z dir
    j = 6*(*nb_pt); nb_pt++;
    prn_Z_PRECISION( pbuf, phi_pt );
    mvmh_PRECISION( op->prpZ+j, D_pt, pbuf );
    mvmh_PRECISION( op->prpZ+j+3, D_pt, pbuf+3 ); D_pt += 9;
    // Y dir
    j = 6*(*nb_pt); nb_pt++;
    prn_Y_PRECISION( pbuf, phi_pt );
    mvmh_PRECISION( op->prpY+j, D_pt, pbuf );
    mvmh_PRECISION( op->prpY+j+3, D_pt, pbuf+3 ); D_pt += 9;
    // X dir
    j = 6*(*nb_pt); nb_pt++;
    prn_X_PRECISION( pbuf, phi_pt );
    mvmh_PRECISION( op->prpX+j, D_pt, pbuf );
    mvmh_PRECISION( op->prpX+j+3, D_pt, pbuf+3 ); D_pt += 9;
  }
#endif
  
  // start communication in positive direction
  START_LOCKED_MASTER(threading)
  ghost_sendrecv_PRECISION( op->prpT, T, +1, &(op->c), _FULL_SYSTEM, l );
  ghost_sendrecv_PRECISION( op->prpZ, Z, +1, &(op->c), _FULL_SYSTEM, l );
  ghost_sendrecv_PRECISION( op->prpY, Y, +1, &(op->c), _FULL_SYSTEM, l );
  ghost_sendrecv_PRECISION( op->prpX, X, +1, &(op->c), _FULL_SYSTEM, l );
  // wait for communication in negative direction
  ghost_wait_PRECISION( op->prnT, T, -1, &(op->c), _FULL_SYSTEM, l );
  ghost_wait_PRECISION( op->prnZ, Z, -1, &(op->c), _FULL_SYSTEM, l );
  ghost_wait_PRECISION( op->prnY, Y, -1, &(op->c), _FULL_SYSTEM, l );
  ghost_wait_PRECISION( op->prnX, X, -1, &(op->c), _FULL_SYSTEM, l );
  END_LOCKED_MASTER(threading)
  
  // multiply with U and lift up minus dir
#ifdef OPTIMIZED_NEIGHBOR_COUPLING_PRECISION
  su3_pbp_PRECISION( eta, prn, op, neighbor, start, end );
#else
  for ( eta_pt=eta+start, end_pt=eta+end, D_pt = op->D+start*3, nb_pt=neighbor+(start/12)*4; eta_pt<end_pt; eta_pt+=12 ) {
    // T dir
    j = 6*(*nb_pt); nb_pt++;
    mvm_PRECISION( pbuf, D_pt, op->prnT+j );
    mvm_PRECISION( pbuf+3, D_pt, op->prnT+j+3 );
    pbp_su3_T_PRECISION( pbuf, eta_pt ); D_pt += 9;
    // Z dir
    j = 6*(*nb_pt); nb_pt++;
    mvm_PRECISION( pbuf, D_pt, op->prnZ+j );
    mvm_PRECISION( pbuf+3, D_pt, op->prnZ+j+3 );
    pbp_su3_Z_PRECISION( pbuf, eta_pt ); D_pt += 9;
    // Y dir
    j = 6*(*nb_pt); nb_pt++;
    mvm_PRECISION( pbuf, D_pt, op->prnY+j );
    mvm_PRECISION( pbuf+3, D_pt, op->prnY+j+3 );
    pbp_su3_Y_PRECISION( pbuf, eta_pt ); D_pt += 9;
    // X dir
    j = 6*(*nb_pt); nb_pt++;
    mvm_PRECISION( pbuf, D_pt, op->prnX+j );
    mvm_PRECISION( pbuf+3, D_pt, op->prnX+j+3 );
    pbp_su3_X_PRECISION( pbuf, eta_pt ); D_pt += 9;
  }
#endif
  
  // wait for communication in positive direction
  START_LOCKED_MASTER(threading)
  ghost_wait_PRECISION( op->prpT, T, +1, &(op->c), _FULL_SYSTEM, l );
  ghost_wait_PRECISION( op->prpZ, Z, +1, &(op->c), _FULL_SYSTEM, l );
  ghost_wait_PRECISION( op->prpY, Y, +1, &(op->c), _FULL_SYSTEM, l );
  ghost_wait_PRECISION( op->prpX, X, +1, &(op->c), _FULL_SYSTEM, l );
  END_LOCKED_MASTER(threading)
  
  // lift up plus dir
#ifdef OPTIMIZED_NEIGHBOR_COUPLING_PRECISION
  pbn_PRECISION( eta, prp, start, end );
#else
  for ( i=start/2, eta_pt=eta+start; i<end/2; i+=6, eta_pt+=12 ) {
    pbn_su3_T_PRECISION( op->prpT+i, eta_pt );
    pbn_su3_Z_PRECISION( op->prpZ+i, eta_pt );
    pbn_su3_Y_PRECISION( op->prpY+i, eta_pt );
    pbn_su3_X_PRECISION( op->prpX+i, eta_pt );
  }
#endif
#ifdef HAVE_TM1p1
  }
#endif

  START_MASTER(threading)
  PROF_PRECISION_STOP( _NC, 1 );
  END_MASTER(threading)
  
  SYNC_MASTER_TO_ALL(threading)
}


void gamma5_PRECISION( vector_PRECISION eta, vector_PRECISION phi, level_struct *l, struct Thread *threading ) {
  
  ASSERT(l->depth == 0);

  vector_PRECISION eta_end = eta + threading->end_index[l->depth];
  eta += threading->start_index[l->depth];
  phi += threading->start_index[l->depth];
#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 ) {
    while ( eta < eta_end ) {
      FOR12( *eta = -(*phi); phi++; eta++; )
      FOR12( *eta =  (*phi); phi++; eta++; )
    }
  } else
#endif
  while ( eta < eta_end ) {
    FOR6( *eta = -(*phi); phi++; eta++; )
    FOR6( *eta =  (*phi); phi++; eta++; )
  }
}

void tau1_gamma5_PRECISION( vector_PRECISION eta, vector_PRECISION phi, level_struct *l, struct Thread *threading ) {
  
  ASSERT(l->depth == 0);

#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 ) {
    vector_PRECISION eta_end = eta + threading->end_index[l->depth];
    complex_PRECISION b[6];
    eta += threading->start_index[l->depth];
    phi += threading->start_index[l->depth];
    while ( eta < eta_end ) {
      int i = 0;
      FOR6( b[i] =  (*phi); phi++; i++;   );
      FOR6( *eta = -(*phi); phi++; eta++; );
      i = 0;
      FOR6( *eta = - b[i] ; eta++; i++;   );
      i = 0;
      FOR6( b[i] =  (*phi); phi++; i++;   );
      FOR6( *eta =  (*phi); phi++; eta++; );
      i = 0;
      FOR6( *eta =   b[i] ; eta++; i++;   );
    }
  } else 
#endif
    {
      START_MASTER(threading)
      warning0("tau1_gamma5_PRECISION called with g.n_flavours != 2\n");
      END_MASTER(threading)
      gamma5_PRECISION( eta, phi, l, threading );
    }
}

void set_even_to_zero_PRECISION( vector_PRECISION eta, vector_PRECISION phi, level_struct *l, struct Thread *threading ) {

  ASSERT(l->depth == 0);
  
  int i = threading->start_site[l->depth];
  vector_PRECISION eta_end = eta + threading->end_index[l->depth];
  eta += threading->start_index[l->depth];
  phi += threading->start_index[l->depth];

#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 )
    while ( eta < eta_end ) {
      if(g.odd_even_table[i]==_ODD) {
        FOR24( *eta = (*phi); phi++; eta++; );
      }
      else if(g.odd_even_table[i]==_EVEN) {
        FOR24( *eta = _COMPLEX_PRECISION_ZERO; phi++; eta++; );
      }
      i++;
    }
  else
#endif
    while ( eta < eta_end ) {
      if(g.odd_even_table[i]==_ODD) {
        FOR12( *eta = (*phi); phi++; eta++; );
      }
      else if(g.odd_even_table[i]==_EVEN) {
        FOR12( *eta = _COMPLEX_PRECISION_ZERO; phi++; eta++; );
      }
      i++;
    }
}

void gamma5_set_even_to_zero_PRECISION( vector_PRECISION eta, vector_PRECISION phi, level_struct *l, struct Thread *threading ) {

  ASSERT(l->depth == 0);
  
  int i = threading->start_site[l->depth];
  vector_PRECISION eta_end = eta + threading->end_index[l->depth];
  eta += threading->start_index[l->depth];
  phi += threading->start_index[l->depth];

#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 )
    while ( eta < eta_end ) {
      if(g.odd_even_table[i]==_ODD){
        FOR12( *eta = -(*phi); phi++; eta++; );
        FOR12( *eta = (*phi); phi++; eta++; );
      }
      else if(g.odd_even_table[i]==_EVEN){
        FOR24( *eta = _COMPLEX_PRECISION_ZERO; phi++; eta++; );
      }
      i++;
    }
  else
#endif
    while ( eta < eta_end ) {
      if(g.odd_even_table[i]==_ODD){
        FOR6( *eta = -(*phi); phi++; eta++; );
        FOR6( *eta = (*phi); phi++; eta++; );
      }
      else if(g.odd_even_table[i]==_EVEN){
        FOR12( *eta = _COMPLEX_PRECISION_ZERO; phi++; eta++; );
      }
      i++;
    }
}

void tau1_gamma5_set_even_to_zero_PRECISION( vector_PRECISION eta, vector_PRECISION phi, level_struct *l, struct Thread *threading ) {

  ASSERT(l->depth == 0);
  
#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 ) {
    int i = threading->start_site[l->depth];
    vector_PRECISION eta_end = eta + threading->end_index[l->depth];
    eta += threading->start_index[l->depth];
    phi += threading->start_index[l->depth];

    complex_PRECISION b[6];
    while ( eta < eta_end ) {
      if(g.odd_even_table[i]==_ODD){
        int i = 0;
        FOR6( b[i] =  (*phi); phi++; i++;   );
        FOR6( *eta = -(*phi); phi++; eta++; );
        i = 0;
        FOR6( *eta = - b[i] ; eta++; i++;   );
        i = 0;
        FOR6( b[i] =  (*phi); phi++; i++;   );
        FOR6( *eta =  (*phi); phi++; eta++; );
        i = 0;
        FOR6( *eta =   b[i] ; eta++; i++;   );
      } else if(g.odd_even_table[i]==_EVEN){
        FOR24( *eta = _COMPLEX_PRECISION_ZERO; phi++; eta++; );
      }
      i++;
    }
  } else 
#endif
    {
      START_MASTER(threading)
      warning0("tau1_gamma5_set_even_to_zero_PRECISION called with g.n_flavours != 2\n");
      END_MASTER(threading)
      gamma5_set_even_to_zero_PRECISION( eta, phi, l, threading );
    }
}

void set_odd_to_zero_PRECISION( vector_PRECISION eta, vector_PRECISION phi, level_struct *l, struct Thread *threading ) {
   
  int i = threading->start_site[l->depth];
  vector_PRECISION eta_end = eta + threading->end_index[l->depth];
  eta += threading->start_index[l->depth];
  phi += threading->start_index[l->depth];

#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 )
    while ( eta < eta_end ) {
      if(g.odd_even_table[i]==_EVEN){
        FOR24( *eta = (*phi); phi++; eta++; );
      }
      else if(g.odd_even_table[i]==_ODD){
        FOR24( *eta = 0; phi++; eta++; );
      }
      i++;
    }
  else
#endif
    while ( eta < eta_end ) {
      if(g.odd_even_table[i]==_EVEN) {
        FOR12( *eta = (*phi); phi++; eta++; );
      }
      else if(g.odd_even_table[i]==_ODD) {
        FOR12( *eta = 0; phi++; eta++; );
      }
      i++;
    }
}

void gamma5_set_odd_to_zero_PRECISION( vector_PRECISION eta, vector_PRECISION phi, level_struct *l, struct Thread *threading ) {
  
  int i = threading->start_site[l->depth];
  vector_PRECISION eta_end = eta + threading->end_index[l->depth];
  eta += threading->start_index[l->depth];
  phi += threading->start_index[l->depth];

#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 )
    while ( eta < eta_end ) {
      if(g.odd_even_table[i]==_EVEN){
        FOR12( *eta = -(*phi); phi++; eta++; );
        FOR12( *eta = (*phi); phi++; eta++; );
      }
      else if(g.odd_even_table[i]==_ODD){
        FOR24( *eta = 0; phi++; eta++; );
      }
      i++;
    }
  else
#endif
    while ( eta < eta_end ) {
      if(g.odd_even_table[i]==_EVEN){
        FOR6( *eta = -(*phi); phi++; eta++; );
        FOR6( *eta = (*phi); phi++; eta++; );
      }
      else if(g.odd_even_table[i]==_ODD){
        FOR12( *eta = 0; phi++; eta++; );
      }
      i++;
    }
}

void tau1_gamma5_set_odd_to_zero_PRECISION( vector_PRECISION eta, vector_PRECISION phi, level_struct *l, struct Thread *threading ) {

  ASSERT(l->depth == 0);
  
#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 ) {
    int i = threading->start_site[l->depth];
    vector_PRECISION eta_end = eta + threading->end_index[l->depth];
    eta += threading->start_index[l->depth];
    phi += threading->start_index[l->depth];
    
    complex_PRECISION b[6];
    while ( eta < eta_end ) {
      if(g.odd_even_table[i]==_EVEN){
        int i = 0;
        FOR6( b[i] =  (*phi); phi++; i++;   );
        FOR6( *eta = -(*phi); phi++; eta++; );
        i = 0;
        FOR6( *eta = - b[i] ; eta++; i++;   );
        i = 0;
        FOR6( b[i] =  (*phi); phi++; i++;   );
        FOR6( *eta =  (*phi); phi++; eta++; );
        i = 0;
        FOR6( *eta =   b[i] ; eta++; i++;   );
      } else if(g.odd_even_table[i]==_ODD){
        FOR24( *eta = _COMPLEX_PRECISION_ZERO; phi++; eta++; );
      }
      i++;
    }
  } else 
#endif
    {
      START_MASTER(threading)
      warning0("tau1_gamma5_set_odd_to_zero_PRECISION called with g.n_flavours != 2\n");
      END_MASTER(threading)
      gamma5_set_odd_to_zero_PRECISION( eta, phi, l, threading );
    }
}

void scale_even_odd_PRECISION( vector_PRECISION eta, vector_PRECISION phi, complex_double even, complex_double odd, 
                               level_struct *l, struct Thread *threading ) {
   
  int i = threading->start_site[l->depth];
  vector_PRECISION eta_end = eta + threading->end_index[l->depth];
  eta += threading->start_index[l->depth];
  phi += threading->start_index[l->depth];

#ifdef HAVE_TM1p1
  if ( g.n_flavours == 2 )
    while ( eta < eta_end ) {
      if(g.odd_even_table[i]==_EVEN){
        FOR24( *eta = even*(*phi); phi++; eta++; );
      }
      else if(g.odd_even_table[i]==_ODD){
        FOR24( *eta = odd*(*phi); phi++; eta++; );
      }
      i++;
    }
  else
#endif
    while ( eta < eta_end ) {
      if(g.odd_even_table[i]==_EVEN) {
        FOR12( *eta = even*(*phi); phi++; eta++; );
      }
      else if(g.odd_even_table[i]==_ODD) {
        FOR12( *eta = odd*(*phi); phi++; eta++; );
      }
      i++;
    }
}


void two_flavours_to_serial_PRECISION( vector_PRECISION flav1, vector_PRECISION flav2, vector_PRECISION serial, level_struct *l, struct Thread *threading ) {

#ifdef HAVE_TM1p1

  /*
   * Order: spin0and1 of flav1
   *        spin0and1 of flav2
   *        spin2and3 of flav1
   *        spin2and3 of flav2
   */
  vector_PRECISION serial_end;
  
  if( g.n_flavours == 2 ) {
    serial_end = serial + threading->end_index[l->depth];
    serial += threading->start_index[l->depth];
    flav1 += threading->start_index[l->depth]/2;
    flav2 += threading->start_index[l->depth]/2;
  }
  else {
    serial_end = serial + threading->end_index[l->depth]*2;
    serial += threading->start_index[l->depth]*2;
    flav1 += threading->start_index[l->depth];
    flav2 += threading->start_index[l->depth];
  }

  while ( serial < serial_end ) {
    FOR6( *serial = (*flav1); serial++; flav1++; )
    FOR6( *serial = (*flav2); serial++; flav2++; )
    FOR6( *serial = (*flav1); serial++; flav1++; )
    FOR6( *serial = (*flav2); serial++; flav2++; )
  }
#else
  START_MASTER(threading)
  warning0("two_flavours_to_serial_PRECISION called without HAVE_TM1p1 defined\n");
  END_MASTER(threading)
#endif
    
}

void serial_to_two_flavours_PRECISION( vector_PRECISION flav1, vector_PRECISION flav2, vector_PRECISION serial, level_struct *l, struct Thread *threading ) {

#ifdef HAVE_TM1p1
  vector_PRECISION serial_end;
  
  if( g.n_flavours == 2 ) {
    serial_end = serial + threading->end_index[l->depth];
    serial += threading->start_index[l->depth];
    flav1 += threading->start_index[l->depth]/2;
    flav2 += threading->start_index[l->depth]/2;
  }
  else {
    serial_end = serial + threading->end_index[l->depth]*2;
    serial += threading->start_index[l->depth]*2;
    flav1 += threading->start_index[l->depth];
    flav2 += threading->start_index[l->depth];
  }

  while ( serial < serial_end ) {
    FOR6( *flav1 = (*serial); serial++; flav1++; )
    FOR6( *flav2 = (*serial); serial++; flav2++; )
    FOR6( *flav1 = (*serial); serial++; flav1++; )
    FOR6( *flav2 = (*serial); serial++; flav2++; )
  }
#else
  START_MASTER(threading)
  warning0("two_flavours_to_serial_PRECISION called without HAVE_TM1p1 defined\n");
  END_MASTER(threading)
#endif
    
}

void g5D_plus_clover_PRECISION( vector_PRECISION eta, vector_PRECISION phi, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {
  d_plus_clover_PRECISION( eta, phi, op, l, threading );
  SYNC_CORES(threading)
  gamma5_PRECISION( eta, eta, l, threading );
  SYNC_CORES(threading)
}

void diagonal_aggregate_PRECISION( vector_PRECISION eta1, vector_PRECISION eta2, vector_PRECISION phi, config_PRECISION diag, level_struct *l ) {

  vector_PRECISION eta_end = eta1 + l->inner_vector_size;
  
  while ( eta1 < eta_end ) {
    FOR6( *eta1 = (*phi)*(*diag); *eta2 = _COMPLEX_PRECISION_ZERO; eta1++; eta2++; phi++; diag++; );
    FOR6( *eta2 = (*phi)*(*diag); *eta1 = _COMPLEX_PRECISION_ZERO; eta1++; eta2++; phi++; diag++; );
  }
}


void d_plus_clover_aggregate_PRECISION( vector_PRECISION eta1, vector_PRECISION eta2, vector_PRECISION phi, schwarz_PRECISION_struct *s, level_struct *l ) {
  
  int i, length, index1, index2, *index_dir, *neighbor = s->op.neighbor_table;
  vector_PRECISION eta1_pt, eta2_pt, phi_pt;
  complex_PRECISION buffer1[12], buffer2[12];
  config_PRECISION D_pt, D = s->op.D;
    
  // add clover term/shift
  spin0and1_clover_PRECISION( eta1, phi, s->op.clover, l );
  spin2and3_clover_PRECISION( eta2, phi, s->op.clover, l );
  
  // T dir
  length = l->is_PRECISION.agg_length[T]; index_dir = l->is_PRECISION.agg_index[T];
  for ( i=0; i<length; i++ ) {
    index1 = index_dir[i]; index2 = neighbor[4*index1 + T];
    phi_pt = phi + 12*index2; D_pt = D + 36*index1+9*T;
    mvm_PRECISION( buffer1, D_pt, phi_pt );
    mvm_PRECISION( buffer1+3, D_pt, phi_pt+3 );
    mvm_PRECISION( buffer1+6, D_pt, phi_pt+6 );
    mvm_PRECISION( buffer1+9, D_pt, phi_pt+9 );
    phi_pt = phi + 12*index1;
    mvmh_PRECISION( buffer2, D_pt, phi_pt );
    mvmh_PRECISION( buffer2+3, D_pt, phi_pt+3 );
    mvmh_PRECISION( buffer2+6, D_pt, phi_pt+6 );
    mvmh_PRECISION( buffer2+9, D_pt, phi_pt+9 );
    eta1_pt = eta1 + 12*index1; eta2_pt = eta2 + 12*index1;
    twospin_p_T_PRECISION( eta1_pt, eta2_pt, buffer1 );
    eta1_pt = eta1 + 12*index2; eta2_pt = eta2 + 12*index2;
    twospin_n_T_PRECISION( eta1_pt, eta2_pt, buffer2 );
  }
  // Z dir
  length = l->is_PRECISION.agg_length[Z]; index_dir = l->is_PRECISION.agg_index[Z];
  for ( i=0; i<length; i++ ) {
    index1 = index_dir[i]; index2 = neighbor[4*index1 + Z];
    phi_pt = phi + 12*index2; D_pt = D + 36*index1+9*Z;
    mvm_PRECISION( buffer1, D_pt, phi_pt );
    mvm_PRECISION( buffer1+3, D_pt, phi_pt+3 );
    mvm_PRECISION( buffer1+6, D_pt, phi_pt+6 );
    mvm_PRECISION( buffer1+9, D_pt, phi_pt+9 );
    phi_pt = phi + 12*index1;
    mvmh_PRECISION( buffer2, D_pt, phi_pt );
    mvmh_PRECISION( buffer2+3, D_pt, phi_pt+3 );
    mvmh_PRECISION( buffer2+6, D_pt, phi_pt+6 );
    mvmh_PRECISION( buffer2+9, D_pt, phi_pt+9 );
    eta1_pt = eta1 + 12*index1; eta2_pt = eta2 + 12*index1;
    twospin_p_Z_PRECISION( eta1_pt, eta2_pt, buffer1 );
    eta1_pt = eta1 + 12*index2; eta2_pt = eta2 + 12*index2;
    twospin_n_Z_PRECISION( eta1_pt, eta2_pt, buffer2 );
  }
  // Y dir
  length = l->is_PRECISION.agg_length[Y]; index_dir = l->is_PRECISION.agg_index[Y];
  for ( i=0; i<length; i++ ) {
    index1 = index_dir[i]; index2 = neighbor[4*index1 + Y];
    phi_pt = phi + 12*index2; D_pt = D + 36*index1+9*Y;
    mvm_PRECISION( buffer1, D_pt, phi_pt );
    mvm_PRECISION( buffer1+3, D_pt, phi_pt+3 );
    mvm_PRECISION( buffer1+6, D_pt, phi_pt+6 );
    mvm_PRECISION( buffer1+9, D_pt, phi_pt+9 );
    phi_pt = phi + 12*index1;
    mvmh_PRECISION( buffer2, D_pt, phi_pt );
    mvmh_PRECISION( buffer2+3, D_pt, phi_pt+3 );
    mvmh_PRECISION( buffer2+6, D_pt, phi_pt+6 );
    mvmh_PRECISION( buffer2+9, D_pt, phi_pt+9 );
    eta1_pt = eta1 + 12*index1; eta2_pt = eta2 + 12*index1;
    twospin_p_Y_PRECISION( eta1_pt, eta2_pt, buffer1 );
    eta1_pt = eta1 + 12*index2; eta2_pt = eta2 + 12*index2;
    twospin_n_Y_PRECISION( eta1_pt, eta2_pt, buffer2 );
  }
  // X dir
  length = l->is_PRECISION.agg_length[X]; index_dir = l->is_PRECISION.agg_index[X];
  for ( i=0; i<length; i++ ) {
    index1 = index_dir[i]; index2 = neighbor[4*index1 + X];
    phi_pt = phi + 12*index2; D_pt = D + 36*index1+9*X;
    mvm_PRECISION( buffer1, D_pt, phi_pt );
    mvm_PRECISION( buffer1+3, D_pt, phi_pt+3 );
    mvm_PRECISION( buffer1+6, D_pt, phi_pt+6 );
    mvm_PRECISION( buffer1+9, D_pt, phi_pt+9 );
    phi_pt = phi + 12*index1;
    mvmh_PRECISION( buffer2, D_pt, phi_pt );
    mvmh_PRECISION( buffer2+3, D_pt, phi_pt+3 );
    mvmh_PRECISION( buffer2+6, D_pt, phi_pt+6 );
    mvmh_PRECISION( buffer2+9, D_pt, phi_pt+9 );
    eta1_pt = eta1 + 12*index1; eta2_pt = eta2 + 12*index1;
    twospin_p_X_PRECISION( eta1_pt, eta2_pt, buffer1 );
    eta1_pt = eta1 + 12*index2; eta2_pt = eta2 + 12*index2;
    twospin_n_X_PRECISION( eta1_pt, eta2_pt, buffer2 );
  }
}

void d_neighbor_aggregate_PRECISION( vector_PRECISION eta1, vector_PRECISION eta2, vector_PRECISION phi, const int mu, schwarz_PRECISION_struct *s, level_struct *l ) {
  
  int i, length, index1, index2, *index_dir, *neighbor;
  vector_PRECISION eta1_pt, eta2_pt, phi_pt;
  complex_PRECISION buffer1[12];
  config_PRECISION D_pt, D = s->op.D;
  
  length = l->is_PRECISION.agg_boundary_length[mu];
  index_dir = l->is_PRECISION.agg_boundary_index[mu];
  neighbor = l->is_PRECISION.agg_boundary_neighbor[mu];
  
  // requires the positive boundaries of phi to be communicated befor
  if ( mu == T ) {
    // T dir
    for ( i=0; i<length; i++ ) {
      index1 = index_dir[i]; index2 = neighbor[i];
      phi_pt = phi + 12*index2; D_pt = D + 36*index1 + 9*T;
      mvm_PRECISION( buffer1, D_pt, phi_pt );
      mvm_PRECISION( buffer1+3, D_pt, phi_pt+3 );
      mvm_PRECISION( buffer1+6, D_pt, phi_pt+6 );
      mvm_PRECISION( buffer1+9, D_pt, phi_pt+9 );
      eta1_pt = eta1 + 12*index1; eta2_pt = eta2 + 12*index1;   
      twospin2_p_T_PRECISION( eta1_pt, eta2_pt, buffer1 );
    }
  } else if ( mu == Z ) {
    // Z dir
    for ( i=0; i<length; i++ ) {
      index1 = index_dir[i]; index2 = neighbor[i];
      phi_pt = phi + 12*index2; D_pt = D + 36*index1 + 9*Z;
      mvm_PRECISION( buffer1, D_pt, phi_pt );
      mvm_PRECISION( buffer1+3, D_pt, phi_pt+3 );
      mvm_PRECISION( buffer1+6, D_pt, phi_pt+6 );
      mvm_PRECISION( buffer1+9, D_pt, phi_pt+9 );
      eta1_pt = eta1 + 12*index1; eta2_pt = eta2 + 12*index1;
      twospin2_p_Z_PRECISION( eta1_pt, eta2_pt, buffer1 );
    }
  } else if ( mu == Y ) {
    // Y dir
    for ( i=0; i<length; i++ ) {
      index1 = index_dir[i]; index2 = neighbor[i];
      phi_pt = phi + 12*index2; D_pt = D + 36*index1 + 9*Y;
      mvm_PRECISION( buffer1, D_pt, phi_pt );
      mvm_PRECISION( buffer1+3, D_pt, phi_pt+3 );
      mvm_PRECISION( buffer1+6, D_pt, phi_pt+6 );
      mvm_PRECISION( buffer1+9, D_pt, phi_pt+9 );
      eta1_pt = eta1 + 12*index1; eta2_pt = eta2 + 12*index1;
      twospin2_p_Y_PRECISION( eta1_pt, eta2_pt, buffer1 );
    }
  } else if ( mu == X ) {
    // X dir
    for ( i=0; i<length; i++ ) {
      index1 = index_dir[i]; index2 = neighbor[i];
      phi_pt = phi + 12*index2; D_pt = D + 36*index1 + 9*X;
      mvm_PRECISION( buffer1, D_pt, phi_pt );
      mvm_PRECISION( buffer1+3, D_pt, phi_pt+3 );
      mvm_PRECISION( buffer1+6, D_pt, phi_pt+6 );
      mvm_PRECISION( buffer1+9, D_pt, phi_pt+9 );
      eta1_pt = eta1 + 12*index1; eta2_pt = eta2 + 12*index1;
      twospin2_p_X_PRECISION( eta1_pt, eta2_pt, buffer1 );
    }
  }
}
 
void apply_twisted_bc_to_vector_PRECISION( vector_PRECISION eta, vector_PRECISION phi, double *theta, level_struct *l) {
  int t, z, y, x, i;
  int *gl=l->global_lattice, sl[4];
  double phase[4];
  complex_double twisted_bc;
  for (i=0; i<4; i++)
    sl[i] = l->local_lattice[i]*g.my_coords[i];
  
  for (t=0; t<l->local_lattice[0]; t++) {
    phase[T] = theta[T]*((double)sl[T]+t)/(double)gl[T];
    for (z=0; z<l->local_lattice[1]; z++) {
      phase[Z] = phase[T] + theta[Z]*((double)sl[Z]+z)/(double)gl[Z];
      for (y=0; y<l->local_lattice[2]; y++) {
        phase[Y] = phase[Z] + theta[Y]*((double)sl[Y]+y)/(double)gl[Y];
        for (x=0; x<l->local_lattice[3]; x++) {
          phase[X] = phase[Y] + theta[X]*((double)sl[X]+x)/(double)gl[X];
          twisted_bc = exp(I*phase[X]);
#ifdef HAVE_TM1p1
          if( g.n_flavours == 2 ) {
            FOR24( *eta = (*phi)*twisted_bc; phi++; eta++; );
          } else
#endif
            { FOR12( *eta = (*phi)*twisted_bc; phi++; eta++; ) }
        }
      }
    }
  }
}

void operator_updates_PRECISION( level_struct *l, struct Thread *threading ) {

  if ( l->level > 0 ) {
    if ( !l->idle ) {
#ifdef INTERPOLATION_SETUP_LAYOUT_OPTIMIZED_PRECISION
      coarse_operator_PRECISION_setup_vectorized( l->is_PRECISION.operator, l, threading );
      START_LOCKED_MASTER(threading)
#else
      START_LOCKED_MASTER(threading)
      coarse_operator_PRECISION_setup( l->is_PRECISION.interpolation, l );
#endif
      conf_PRECISION_gather( &(l->next_level->s_PRECISION.op), &(l->next_level->op_PRECISION), l->next_level );
      END_LOCKED_MASTER(threading)
      if ( !l->next_level->idle && l->next_level->level > 0 ) {
        START_LOCKED_MASTER(threading)
        schwarz_PRECISION_boundary_update( &(l->next_level->s_PRECISION), l->next_level );
        END_LOCKED_MASTER(threading)
        if ( g.method >= 4 && g.odd_even ) {
          coarse_oddeven_setup_PRECISION( &(l->next_level->s_PRECISION.op), _REORDER, l->next_level, threading );
        } else {
          coarse_operator_PRECISION_set_couplings( &(l->next_level->s_PRECISION.op), l->next_level, threading );
        }
      }
      if ( !l->next_level->idle && l->next_level->level == 0 && g.odd_even ) {
        coarse_oddeven_setup_PRECISION( &(l->next_level->s_PRECISION.op), _NO_REORDERING, l->next_level, threading );
      } else if ( !l->next_level->idle && l->next_level->level == 0 ) {
        coarse_operator_PRECISION_set_couplings( &(l->next_level->s_PRECISION.op), l->next_level, threading );
      }
      operator_updates_PRECISION( l->next_level, threading );
    }
  }  
}


void m0_update_PRECISION( PRECISION m0, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {
  
  // no hyperthreading in this function
  if(threading->thread != 0)
    return;

  config_PRECISION clover = op->clover;
  
  if ( clover != NULL && op->m0 != m0 ) {
    int i, j;
    complex_PRECISION m0_diff = m0 - op->m0;

    START_MASTER(threading)
    op->m0 = m0;
    END_MASTER(threading)

    if( m0_diff != 0 ) {
      if ( l->depth == 0 ) {
        int start = threading->start_site[l->depth];
        int n     = threading->n_site[l->depth];
        clover += start*(g.csw?42:12);
        for ( i=0; i<n; i++ ) {
          for ( j=0; j<12; j++ ) {
            clover[j] += m0_diff;
          }
          // clover term diag also stored as complex, so size is 2*15+2*6 = 42
          clover += (g.csw?42:12);
        }
      } else {
        int start = threading->start_site[l->depth];
        int n     = threading->n_site[l->depth];
        int k = l->num_parent_eig_vect;
        int sc_size = (l->num_parent_eig_vect)*(l->num_parent_eig_vect*2+1);
        clover += start*sc_size;
        for ( i=0; i<n; i++ ) {
          for ( j=0; j<k; j++ ) {
            if ( j>0 ) clover += j+1; 
            *clover += m0_diff;
          }
          clover ++;
          for ( j=0; j<k; j++ ) {
            if ( j>0 ) clover += j+1;
            *clover += m0_diff;
          }
          clover += 1 + SQUARE(k);
        }
      }
    }
  }
}

void tm_term_PRECISION_setup( PRECISION mu, PRECISION even, PRECISION odd, operator_PRECISION_struct *op,
                              level_struct *l, struct Thread *threading ) {
   
#ifdef HAVE_TM
  if(threading->thread != 0)
    return;

  config_PRECISION tm_term = op->tm_term;
  if ( tm_term != NULL ) {
    config_PRECISION odd_proj = op->odd_proj;
    complex_PRECISION shift = I*mu;
    complex_PRECISION even_shift = I*even;
    complex_PRECISION odd_shift = I*odd;

    START_MASTER(threading)
    op->mu = mu;
    op->mu_even_shift = even;
    op->mu_odd_shift = odd;
    END_MASTER(threading)

    int i, j;
    int start, end;
    compute_core_start_end(0, l->num_inner_lattice_sites, &start, &end, l, threading);
    int n = end-start;
          
    if ( l->depth == 0 ) {
      complex_PRECISION tm_shift;
      tm_term += start*12;
      odd_proj += start*12;
      
      for ( i=0; i<n; i++ ) {
        if( cimag(even_shift) == 0. && cimag(odd_shift) == 0. )
          tm_shift = shift;
        else
          tm_shift = shift + even_shift + odd_proj[0]*(odd_shift - even_shift);
        FOR6( *tm_term = - tm_shift; tm_term++; )
        FOR6( *tm_term = tm_shift; tm_term++; )
        odd_proj += 12;
      }
    } else {
      int k, m  = l->num_parent_eig_vect;
      int tm_size = m*(m+1);
      
      tm_term += start*tm_size;
      odd_proj += start*tm_size;
      
      if( cimag(even_shift) == 0. && cimag(odd_shift) == 0. ) {
        
        for ( i=0; i<n; i++ ) {
          for ( j=0; j<m; j++ ) {
            for ( k=0; k<j; k++ )
              tm_term[k] = _COMPLEX_PRECISION_ZERO;
            tm_term += j;
            *tm_term = -1.* shift;
            tm_term++;
          }
          
          for ( j=0; j<m; j++ ) {
            for ( k=0; k<j; k++ )
              tm_term[k] = _COMPLEX_PRECISION_ZERO;
            tm_term += j;
            *tm_term = shift;
            tm_term++;
          }
        }
      } else {
        complex_PRECISION odd_factor = odd_shift - even_shift;
        
        for ( i=0; i<n; i++ ) {
          for ( j=0; j<m; j++ ) {
            for ( k=0; k<j; k++ ) 
              tm_term[k] = -1. * odd_factor * odd_proj[k] ;
            tm_term += j;
            odd_proj += j;
            *tm_term = -1.* ( shift + even_shift + odd_factor * (*odd_proj));
            tm_term++;
            odd_proj++;
          } 
          
          for ( j=0; j<m; j++ ) {
            for ( k=0; k<j; k++ ) 
              tm_term[k] = odd_factor * odd_proj[k] ;
            tm_term += j;
            odd_proj += j;
            *tm_term = ( shift + even_shift + odd_factor * (*odd_proj));
            tm_term++;
            odd_proj++;
          } 
        }
      }
    }  
  }
#endif
}

void epsbar_term_PRECISION_setup( PRECISION epsbar, PRECISION even, PRECISION odd, operator_PRECISION_struct *op,
                                  level_struct *l, struct Thread *threading ) {
  
#ifdef HAVE_TM1p1
  if(threading->thread != 0)
    return;

  config_PRECISION eps_term = op->epsbar_term;
  if ( eps_term != NULL ) {
    config_PRECISION odd_proj = op->odd_proj;
    complex_PRECISION shift = -epsbar;
    complex_PRECISION even_shift = I*even;
    complex_PRECISION odd_shift = I*odd;

    START_MASTER(threading)
    op->epsbar = epsbar;
    op->epsbar_ig5_even_shift = even;
    op->epsbar_ig5_odd_shift = odd;
    END_MASTER(threading)

    int i, j;
    int start, end;
    compute_core_start_end(0, l->num_inner_lattice_sites, &start, &end, l, threading);
    int n = end-start;
          
    if ( l->depth == 0 ) {
      eps_term += start*12;
      odd_proj += start*12;

      if( cimag(even_shift) == 0. && cimag(odd_shift) == 0. )
        for ( i=0; i<n; i++ ) {
          FOR12( *eps_term = shift; eps_term++; );
        }
      else
        for ( i=0; i<n; i++ ) {
          complex_PRECISION ig5_shift = even_shift + (*odd_proj)*(odd_shift - even_shift);
          FOR6( *eps_term = shift-ig5_shift; eps_term++; );
          FOR6( *eps_term = shift+ig5_shift; eps_term++; );
          odd_proj += 12;
      }
    } else {
      int k, m  = l->num_parent_eig_vect;
      int eps_size = m*(m+1);
      
      eps_term += start*eps_size;
      odd_proj += start*eps_size;
      
      if( cimag(even_shift) == 0. && cimag(odd_shift) == 0. ) {
        for ( i=0; i<2*n; i++ ) {
          for ( j=0; j<m; j++ ) {
            for ( k=0; k<j; k++ )
              eps_term[k] = _COMPLEX_PRECISION_ZERO;
            eps_term += j;
            *eps_term = shift;
            eps_term++;
          }
        } 
      } else {
        complex_PRECISION odd_factor = odd_shift - even_shift;
         
        for ( i=0; i<n; i++ ) {
          for ( j=0; j<m; j++ ) {
            for ( k=0; k<j; k++ ) 
              eps_term[k] = -1.* odd_factor*odd_proj[k] ;
            eps_term += j;
            odd_proj += j;
            *eps_term = shift - (even_shift + odd_factor * (*odd_proj));
            eps_term++;
            odd_proj++;
          } 
          
          for ( j=0; j<m; j++ ) {
            for ( k=0; k<j; k++ ) 
              eps_term[k] = odd_factor*odd_proj[k] ;
            eps_term += j;
            odd_proj += j;
            *eps_term = shift + (even_shift + odd_factor * (*odd_proj));
            eps_term++;
            odd_proj++;
          } 
        }
      }
    }  
  }
#endif
}


void two_flavours_test_PRECISION( operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {

#ifdef HAVE_TM1p1
  double diff;
  
  vector_double vd1=NULL, vd2, vd3, vd4, vdd1, vdd2, vdd3, vdd4;
  vector_PRECISION vpp1=NULL, vpp2;

  ASSERT(g.n_flavours==2);

  data_layout_n_flavours( 1, l, threading );

  int ivs = l->inner_vector_size;
  
  PUBLIC_MALLOC( vd1, complex_double, 4*ivs + 2*4*ivs );
  PUBLIC_MALLOC( vpp1, complex_PRECISION, 2*2*ivs );

  vd2 = vd1 + ivs; vd3 = vd2 + ivs; vd4 = vd3 + ivs;
  vdd1 = vd4 + ivs; vdd2 = vdd1 + 2*ivs; vdd3 = vdd2 + 2*ivs; vdd4 = vdd3 + 2*ivs;
  vpp2 = vpp1 + 2*ivs;
  
  START_LOCKED_MASTER(threading)

  vector_double_define_random( vd1, 0, l->inner_vector_size, l );
  vector_double_define_random( vd2, 0, l->inner_vector_size, l );
  apply_operator_double( vd3, vd1, &(g.p), l, no_threading );
#ifdef HAVE_TM
  vector_double_real_scale( g.op_double.tm_term, g.op_double.tm_term, -1, 0, l->inner_vector_size, l ); 
#endif
  apply_operator_double( vd4, vd2, &(g.p), l, no_threading );
#ifdef HAVE_TM
  vector_double_real_scale( g.op_double.tm_term, g.op_double.tm_term, -1, 0, l->inner_vector_size, l ); 
#endif
  add_diagonal_double( vd3, vd2, g.op_double.epsbar_term, l->inner_vector_size );
  add_diagonal_double( vd4, vd1, g.op_double.epsbar_term, l->inner_vector_size );

  two_flavours_to_serial_double( vd1, vd2, vdd1, l, no_threading );
  two_flavours_to_serial_double( vd3, vd4, vdd2, l, no_threading );

  data_layout_n_flavours( 2, l, threading );

  trans_PRECISION( vpp1, vdd1, op->translation_table, l, no_threading );
  apply_operator_PRECISION( vpp2, vpp1, &(l->p_PRECISION), l, no_threading );
  trans_back_PRECISION( vdd3, vpp2, op->translation_table, l, no_threading );
  
  vector_double_minus( vdd4, vdd3, vdd2, 0, l->inner_vector_size, l );
  diff = global_norm_double( vdd4, 0, l->inner_vector_size, l, no_threading ) /
    global_norm_double( vdd3, 0, l->inner_vector_size, l, no_threading );
  
  test0_PRECISION("depth: %d, correctness of doublet Dirac operator PRECISION: %le\n", l->depth, diff );
  END_LOCKED_MASTER(threading)

  if(threading->n_core > 1) {
    trans_PRECISION( vpp1, vdd1, op->translation_table, l, threading );
    apply_operator_PRECISION( vpp2, vpp1, &(l->p_PRECISION), l, threading );
    trans_back_PRECISION( vdd3, vpp2, op->translation_table, l, threading );
    
    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)

    START_LOCKED_MASTER(threading)
    vector_double_minus( vdd4, vdd3, vdd2, 0, l->inner_vector_size, l );
    diff = global_norm_double( vdd4, 0, l->inner_vector_size, l, no_threading ) /
      global_norm_double( vdd3, 0, l->inner_vector_size, l, no_threading );
    
    test0_PRECISION("depth: %d, correctness of doublet Dirac operator PRECISION with threading: %le\n", l->depth, diff );
    END_LOCKED_MASTER(threading)
  }    
  
  PUBLIC_FREE( vd1, complex_double, 4*ivs + 2*4*ivs );
  PUBLIC_FREE( vpp1, complex_PRECISION, 2*2*ivs );

  START_LOCKED_MASTER(threading)
  if ( g.method >=4 && g.odd_even )
    oddeven_PRECISION_test( l );
  END_LOCKED_MASTER(threading) 
#endif
    
}
