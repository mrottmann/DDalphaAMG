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

#if defined(OPTIMIZED_NEIGHBOR_COUPLING_PRECISION) || defined(OPTIMIZED_SELF_COUPLING_PRECISION)
void block_d_plus_clover_PRECISION( vector_PRECISION eta, vector_PRECISION phi, int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {

  START_UNTHREADED_FUNCTION(threading)
  
  int i, n = s->num_block_sites, *length = s->dir_length, **index = s->index, *neighbor = s->op.neighbor_table;
  vector_PRECISION lphi = phi+start, leta = eta+start;
  config_PRECISION clover = (g.csw==0.0)?s->op.clover+start:s->op.clover+(start/12)*42;
#ifdef OPTIMIZED_NEIGHBOR_COUPLING_PRECISION
  PRECISION *Dplus = s->op.D_vectorized + (start/12)*96;
  PRECISION *Dminus = s->op.D_transformed_vectorized + (start/12)*96;
#else
  int j, k, *ind;
  complex_PRECISION buf1[25]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, *buf2=buf1+6, *buf3=buf2+6, *buf4=buf3+6;
  config_PRECISION D_pt;
  config_PRECISION D = s->op.D + (start/12)*36;
#endif
#ifdef OPTIMIZED_SELF_COUPLING_PRECISION
  PRECISION* clover_vectorized = (PRECISION*) (s->op.clover_vectorized+start*12);
#endif
  
  // clover term
  if ( g.csw == 0.0 ) {
    clover_PRECISION( leta, lphi, clover, 12*n, l, threading ); 
  } else {
    for ( i=0; i<n; i++ ) {
#ifdef OPTIMIZED_SELF_COUPLING_PRECISION
      sse_site_clover_PRECISION( (PRECISION*)(leta+12*i), (PRECISION*)(lphi+12*i), clover_vectorized+144*i );
#else
      site_clover_PRECISION( leta+12*i, lphi+12*i, clover+42*i );
#endif
    }
  }
  
#ifdef OPTIMIZED_NEIGHBOR_COUPLING_PRECISION
  for ( int mu=0; mu<4; mu++ ) {
    block_oddeven_plus_coupling_PRECISION( (PRECISION*)leta, Dplus, (PRECISION*)lphi, mu, 0, length[mu], index[mu], neighbor );
    block_oddeven_minus_coupling_PRECISION( (PRECISION*)leta, Dminus, (PRECISION*)lphi, mu, 0, length[mu], index[mu], neighbor );
  }
#else
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
#endif
  END_UNTHREADED_FUNCTION(threading)
}
#endif


#if defined(OPTIMIZED_NEIGHBOR_COUPLING_PRECISION) || defined(OPTIMIZED_SELF_COUPLING_PRECISION)
void sse_clover_PRECISION( vector_PRECISION eta, vector_PRECISION phi, operator_PRECISION_struct *op,
                           int start, int end, level_struct *l, struct Thread *threading );
void d_plus_clover_PRECISION( vector_PRECISION eta, vector_PRECISION phi, operator_PRECISION_struct *op,
                              level_struct *l, struct Thread *threading ) {
  
  int n = l->num_inner_lattice_sites, *neighbor = op->neighbor_table, start, end;
#ifndef OPTIMIZED_NEIGHBOR_COUPLING_PRECISION
  int i, j, *nb_pt;
  complex_PRECISION pbuf[6];
  vector_PRECISION phi_pt, eta_pt, end_pt;
  config_PRECISION D_pt;
#endif
  
  compute_core_start_end(0, 12*n, &start, &end, l, threading );
  
  SYNC_MASTER_TO_ALL(threading)

  if ( g.csw == 0.0 ) {
    vector_PRECISION_scale( eta, phi, op->shift, start, end, l );
  } else {
    clover_PRECISION( eta+start, phi+start, op->clover+((start/12)*42), end-start, l, threading );
  }
  
  START_MASTER(threading)
  PROF_PRECISION_START( _NC ); 
  END_MASTER(threading)
  
#ifdef OPTIMIZED_NEIGHBOR_COUPLING_PRECISION
  complex_PRECISION *prn[4] = { op->prnT, op->prnZ, op->prnY, op->prnX };
  prp_PRECISION( prn, phi, start, end );
#else
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
  complex_PRECISION *prp[4] = { op->prpT, op->prpZ, op->prpY, op->prpX };
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
  
  START_MASTER(threading)
  PROF_PRECISION_STOP( _NC, 1 );
  END_MASTER(threading)
  
  SYNC_MASTER_TO_ALL(threading)
}
#endif


#endif

