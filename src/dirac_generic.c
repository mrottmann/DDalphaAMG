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

void clover_PRECISION( vector_PRECISION eta, vector_PRECISION phi, config_PRECISION clover, int length,
                       level_struct *l, struct Thread *threading ) {
  
  vector_PRECISION eta_end = eta + length;
  if ( g.csw == 0.0 ) {
    while ( eta < eta_end ) {
      FOR12( *eta = (*phi)*(*clover); eta++; phi++; clover++; )
    }
  } else {
    START_MASTER(threading)
    PROF_PRECISION_START( _SC );
    END_MASTER(threading)
    while ( eta < eta_end ) {
      site_clover_PRECISION( eta, phi, clover );
      eta+=12; phi+=12; clover+=42;
    }
    START_MASTER(threading)
    PROF_PRECISION_STOP( _SC, 1 );
    END_MASTER(threading)
  }
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

#if !defined(OPTIMIZED_NEIGHBOR_COUPLING_PRECISION) && !defined(OPTIMIZED_SELF_COUPLING_PRECISION)
void block_d_plus_clover_PRECISION( vector_PRECISION eta, vector_PRECISION phi, int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {

  START_UNTHREADED_FUNCTION(threading)
  
  int i, n = s->num_block_sites, *length = s->dir_length, **index = s->index, *neighbor = s->op.neighbor_table;
  vector_PRECISION lphi = phi+start, leta = eta+start;
  config_PRECISION clover = (g.csw==0.0)?s->op.clover+start:s->op.clover+(start/12)*42;
  int j, k, *ind;
  complex_PRECISION buf1[25]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}, *buf2=buf1+6, *buf3=buf2+6, *buf4=buf3+6;
  config_PRECISION D_pt;
  config_PRECISION D = s->op.D + (start/12)*36;
  
  // clover term
  clover_PRECISION( leta, lphi, clover, 12*n, l, no_threading ); 
#ifdef HAVE_TM
  config_PRECISION tm_term = s->op.tm_term+start;
  if (g.tm_mu + g.tm_mu_odd_shift != 0.0 || g.tm_mu + g.tm_mu_even_shift != 0.0 )
    add_diagonal_PRECISION( leta, lphi, tm_term, 12*n );
#endif
    
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
  END_UNTHREADED_FUNCTION(threading)
}
#endif


#if !defined(OPTIMIZED_NEIGHBOR_COUPLING_PRECISION) && !defined(OPTIMIZED_SELF_COUPLING_PRECISION)
void d_plus_clover_PRECISION( vector_PRECISION eta, vector_PRECISION phi, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {
  
  int n = l->num_inner_lattice_sites, *neighbor = op->neighbor_table, start, end;
  int i, j, *nb_pt;
  complex_PRECISION pbuf[6];
  vector_PRECISION phi_pt, eta_pt, end_pt;
  config_PRECISION D_pt;
  compute_core_start_end(0, 12*n, &start, &end, l, threading );
  vector_PRECISION lphi = phi+start, leta = eta+start;
  config_PRECISION clover = (g.csw==0.0)?op->clover+start:op->clover+(start/12)*42;
    
  SYNC_MASTER_TO_ALL(threading)

  // clover term
  clover_PRECISION( leta, lphi, clover, end-start, l, threading ); 
#ifdef HAVE_TM
  if (g.tm_mu + g.tm_mu_odd_shift != 0.0 || g.tm_mu + g.tm_mu_even_shift != 0.0 )
    add_diagonal_PRECISION( leta, lphi, op->tm_term+start, end-start );
#endif
    
  START_MASTER(threading)
  PROF_PRECISION_START( _NC ); 
  END_MASTER(threading)
  
  for ( i=start/2, phi_pt=phi+start; i<end/2; i+=6, phi_pt+=12 ) {
    prp_T_PRECISION( op->prnT+i, phi_pt );
    prp_Z_PRECISION( op->prnZ+i, phi_pt );
    prp_Y_PRECISION( op->prnY+i, phi_pt );
    prp_X_PRECISION( op->prnX+i, phi_pt );
  }
  // start communication in negative direction
  START_LOCKED_MASTER(threading)
  ghost_sendrecv_PRECISION( op->prnT, T, -1, &(op->c), _FULL_SYSTEM, l );
  ghost_sendrecv_PRECISION( op->prnZ, Z, -1, &(op->c), _FULL_SYSTEM, l );
  ghost_sendrecv_PRECISION( op->prnY, Y, -1, &(op->c), _FULL_SYSTEM, l );
  ghost_sendrecv_PRECISION( op->prnX, X, -1, &(op->c), _FULL_SYSTEM, l );
  END_LOCKED_MASTER(threading) 
  
  // project plus dir and multiply with U dagger
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
  
  // wait for communication in positive direction
  START_LOCKED_MASTER(threading)
  ghost_wait_PRECISION( op->prpT, T, +1, &(op->c), _FULL_SYSTEM, l );
  ghost_wait_PRECISION( op->prpZ, Z, +1, &(op->c), _FULL_SYSTEM, l );
  ghost_wait_PRECISION( op->prpY, Y, +1, &(op->c), _FULL_SYSTEM, l );
  ghost_wait_PRECISION( op->prpX, X, +1, &(op->c), _FULL_SYSTEM, l );
  END_LOCKED_MASTER(threading)
  
  // lift up plus dir
  for ( i=start/2, eta_pt=eta+start; i<end/2; i+=6, eta_pt+=12 ) {
    pbn_su3_T_PRECISION( op->prpT+i, eta_pt );
    pbn_su3_Z_PRECISION( op->prpZ+i, eta_pt );
    pbn_su3_Y_PRECISION( op->prpY+i, eta_pt );
    pbn_su3_X_PRECISION( op->prpX+i, eta_pt );
  }
  
  START_MASTER(threading)
  PROF_PRECISION_STOP( _NC, 1 );
  END_MASTER(threading)
  
  SYNC_MASTER_TO_ALL(threading)
}
#endif


void d_plus_clover_dagger_PRECISION( vector_PRECISION eta, vector_PRECISION phi, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {  
  gamma5_PRECISION( l->vbuf_PRECISION[6], phi, l, threading );
  d_plus_clover_PRECISION( l->vbuf_PRECISION[7], l->vbuf_PRECISION[6], op, l, threading );
  gamma5_PRECISION( eta, l->vbuf_PRECISION[7], l, threading );
}


void gamma5_PRECISION( vector_PRECISION eta, vector_PRECISION phi, level_struct *l, struct Thread *threading ) {
  
  vector_PRECISION eta_end = eta + threading->end_index[l->depth];
  eta += threading->start_index[l->depth];
  phi += threading->start_index[l->depth];
  while ( eta < eta_end ) {
    FOR6( *eta = -(*phi); phi++; eta++; )
    FOR6( *eta =  (*phi); phi++; eta++; )
  }
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
    FOR6( *eta1 = (*phi)*(*diag); *eta2 = _COMPLEX_PRECISION_ZERO; eta1++; eta2++; phi++; diag++; )
    FOR6( *eta2 = (*phi)*(*diag); *eta1 = _COMPLEX_PRECISION_ZERO; eta1++; eta2++; phi++; diag++; )
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
	  FOR12( *eta = (*phi)*twisted_bc; phi++; eta++; );
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
#ifdef HAVE_TM
      l->next_level->tm_shift = g.tm_mu*g.tm_mu_factor[l->next_level->depth];
      l->next_level->tm_even_shift = g.tm_mu_even_shift*g.tm_mu_factor[l->next_level->depth];
      l->next_level->tm_odd_shift = g.tm_mu_odd_shift*g.tm_mu_factor[l->next_level->depth];
      
      if( g.tm_mu_factor[l->next_level->depth]!=g.tm_mu_factor[l->depth] )
	tm_term_PRECISION_setup( l->next_level->op_PRECISION.tm_term, l->next_level->op_PRECISION.odd_proj, l->next_level, no_threading );
#endif
      conf_PRECISION_gather( &(l->next_level->s_PRECISION.op), &(l->next_level->op_PRECISION), l->next_level );
      END_LOCKED_MASTER(threading)
      if ( !l->next_level->idle && l->next_level->level > 0 ) {
        START_LOCKED_MASTER(threading)
        schwarz_PRECISION_boundary_update( &(l->next_level->s_PRECISION), l->next_level );
        END_LOCKED_MASTER(threading)
        if ( g.method >= 4 && g.odd_even ) {
          coarse_oddeven_re_setup_PRECISION( &(l->next_level->s_PRECISION.op), _REORDER, l->next_level, threading );
        } else {
          coarse_operator_PRECISION_set_couplings( &(l->next_level->s_PRECISION.op), l->next_level, threading );
        }
      }
      if ( !l->next_level->idle && l->next_level->level == 0 && g.odd_even ) {
        coarse_oddeven_re_setup_PRECISION( &(l->next_level->s_PRECISION.op), _NO_REORDERING, l->next_level, threading );
      } else if ( !l->next_level->idle && l->next_level->level == 0 ) {
        coarse_operator_PRECISION_set_couplings( &(l->next_level->s_PRECISION.op), l->next_level, threading );
      }
      operator_updates_PRECISION( l->next_level, threading );
    }
  }  
}


void shift_update_PRECISION( operator_PRECISION_struct *op, complex_PRECISION shift, level_struct *l, struct Thread *threading ) {
  
  // no hyperthreading in this function
  if(threading->thread != 0)
    return;

  config_PRECISION clover = op->clover;
  
  if ( clover != NULL ) {
    int i, j;
    complex_PRECISION old_shift = (complex_PRECISION) l->dirac_shift;
    complex_PRECISION shift_diff = shift - old_shift;
    
    if ( l->depth == 0 ) {
      int start = threading->start_site[l->depth];
      int n     = threading->n_site[l->depth];
      clover += start*(g.csw?42:12);
      for ( i=0; i<n; i++ ) {
        for ( j=0; j<12; j++ ) {
          clover[j] += shift_diff;
        }
        // clover term diag also stored as complex, so size is 2*15+2*6 = 42
        clover += (g.csw?42:12);
      }
    } else {
      int start = threading->start_site[l->depth];
      int n     = threading->n_site[l->depth];
      int k = l->num_lattice_site_var/2;
      int sc_size = (l->num_lattice_site_var/2)*(l->num_lattice_site_var+1);
      clover += start*sc_size;
      for ( i=0; i<n; i++ ) {
        for ( j=0; j<k; j++ ) {
          if ( j>0 ) clover += j+1; 
          *clover += shift_diff;
        }
        clover ++;
        for ( j=0; j<k; j++ ) {
          if ( j>0 ) clover += j+1;
          *clover += shift_diff;
        }
        clover += 1 + SQUARE(k);
      }
    }
  }
}

void tm_term_PRECISION_setup( config_PRECISION tm_term, config_PRECISION odd_proj, level_struct *l, struct Thread *threading ) {
   
#ifdef HAVE_TM
  if(threading->thread != 0)
    return;

  complex_PRECISION shift = I*l->tm_shift;
  complex_PRECISION even_shift = I*l->tm_even_shift;
  complex_PRECISION odd_shift = I*l->tm_odd_shift;

  if ( tm_term != NULL ) {
    int i, j;
    int start, end;
    compute_core_start_end(0, l->num_inner_lattice_sites, &start, &end, l, threading);
    int n = end-start;
    complex_PRECISION tm_shift;
          
    if ( l->depth == 0 ) {
      tm_term += start*12;
      odd_proj += start*12;
      
      for ( i=0; i<n; i++ ) {
	if( cimag(even_shift) == 0. && cimag(odd_shift) == 0. )
	  tm_shift = shift;
	else
	  tm_shift = shift + even_shift + odd_proj[0]*(odd_shift - even_shift);
        for ( j=0; j<6; j++ ) {
          tm_term[j] = - tm_shift;
        }
        for ( j=6; j<12; j++ ) {
          tm_term[j] = tm_shift;
        }
	tm_term += 12;
	odd_proj += 12;
      }
    } else {
      int k, m  = l->num_lattice_site_var/2;
      int tm_size = (l->num_lattice_site_var/2)*(l->num_lattice_site_var/2+1);
            
      tm_term += start*tm_size;
      odd_proj += start*tm_size;

      if( cimag(even_shift) == 0. && cimag(odd_shift) == 0. ) {
	
	tm_shift = shift;
	
	for ( i=0; i<n; i++ ) {
	  for ( j=0; j<m; j++ ) {
	    for ( k=0; k<j; k++ )
	      tm_term[k] = _COMPLEX_PRECISION_ZERO;
	    tm_term += j;
	    *tm_term = -1.* tm_shift;
	    tm_term++;
	  }
	  
	  for ( j=0; j<m; j++ ) {
	    for ( k=0; k<j; k++ )
	      tm_term[k] = _COMPLEX_PRECISION_ZERO;
	    tm_term += j;
	    *tm_term = tm_shift;
	    tm_term++;
	  }
	}
      } else {
	complex_PRECISION odd_factor = odd_shift - even_shift;
	tm_shift = shift + even_shift;
	
	for ( i=0; i<n; i++ ) {
	  for ( j=0; j<m; j++ ) {
	    for ( k=0; k<j; k++ ) 
	      tm_term[k] = -1.* odd_factor*odd_proj[k] ;
	    tm_term += j;
	    odd_proj += j;
	    *tm_term = -1.* (tm_shift + odd_factor * (*odd_proj));
	    tm_term++;
	    odd_proj++;
	  } 
        
	  for ( j=0; j<m; j++ ) {
	    for ( k=0; k<j; k++ ) 
	      tm_term[k] = odd_factor*odd_proj[k] ;
	    tm_term += j;
	    odd_proj += j;
	    *tm_term = (tm_shift + odd_factor * (*odd_proj));
	    tm_term++;
	    odd_proj++;
	  } 
	}
      }
    }
  }
#endif
}

void optimized_shift_update_PRECISION( complex_PRECISION mass_shift, level_struct *l, struct Thread *threading ) {

  if ( !l->idle ) {

    if ( mass_shift !=  l->dirac_shift ) {
      shift_update_PRECISION( &(l->op_PRECISION), mass_shift, l, threading );
      shift_update_PRECISION( &(l->s_PRECISION.op), mass_shift, l, threading );
      START_LOCKED_MASTER(threading)
      l->dirac_shift = mass_shift;
      l->real_shift = creal(mass_shift);
      END_LOCKED_MASTER(threading)
    }
    
#ifdef HAVE_TM
    if ( l->tm_shift != g.tm_mu*g.tm_mu_factor[l->depth] ||
	 l->tm_even_shift != g.tm_mu_even_shift*g.tm_mu_factor[l->depth] ||
	 l->tm_odd_shift != g.tm_mu_odd_shift*g.tm_mu_factor[l->depth] ) {
      START_LOCKED_MASTER(threading)
      if( g.tm_mu_even_shift == g.tm_mu_odd_shift )
	  printf0("depth: %d, updating mu to %f \n", (l->depth), cimag(g.tm_mu+g.tm_mu_even_shift));
      else  
	printf0("depth: %d, updating mu to %f on even sites and %f on odd sites \n", l->depth, cimag(g.tm_mu+g.tm_mu_even_shift), cimag(g.tm_mu+g.tm_mu_even_shift));
      
      l->tm_shift = g.tm_mu*g.tm_mu_factor[l->depth];
      l->tm_even_shift = g.tm_mu_even_shift*g.tm_mu_factor[l->depth];
      l->tm_odd_shift = g.tm_mu_odd_shift*g.tm_mu_factor[l->depth];
      END_LOCKED_MASTER(threading)
	
      tm_term_PRECISION_setup( l->op_PRECISION.tm_term, l->op_PRECISION.odd_proj, l, threading ); 
      tm_term_PRECISION_setup( l->s_PRECISION.op.tm_term, l->s_PRECISION.op.odd_proj, l, threading );
    }
#endif
  
    SYNC_CORES(threading)

    if ( !l->idle && g.method >= 4 && l->level > 0 && g.odd_even ) 
      coarse_oddeven_re_setup_PRECISION( &(l->s_PRECISION.op), _REORDER, l, threading );
    else if ( !l->idle && l->level == 0 && g.odd_even)
      coarse_oddeven_re_setup_PRECISION( &(l->s_PRECISION.op), _NO_REORDERING, l, threading );
    else
      coarse_operator_PRECISION_set_couplings_clover( &(l->s_PRECISION.op), l, threading );
   
    if(l->level > 0)
      optimized_shift_update_PRECISION( mass_shift, l->next_level, threading );
  }
}
