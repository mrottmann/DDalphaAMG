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
      FOR6( *eta = 0; eta++; )
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
      FOR6( *eta = 0; eta++; )
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
  if ( g.csw == 0.0 ) {
    clover_PRECISION( leta, lphi, clover, 12*n, l, threading ); 
  } else {
    for ( i=0; i<n; i++ ) {
      site_clover_PRECISION( leta+12*i, lphi+12*i, clover+42*i );
    }
  }
  
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
  
  SYNC_MASTER_TO_ALL(threading)

  if ( g.csw == 0.0 ) {
    vector_PRECISION_scale( eta, phi, op->shift, start, end, l );
  } else {
    clover_PRECISION( eta+start, phi+start, op->clover+((start/12)*42), end-start, l, threading );
  }
  
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


void d_plus_clover_aggregate_PRECISION( vector_PRECISION eta1, vector_PRECISION eta2, vector_PRECISION phi, schwarz_PRECISION_struct *s, level_struct *l ) {
  
  int i, length, index1, index2, *index_dir, *neighbor = s->op.neighbor_table;
  vector_PRECISION eta1_pt, eta2_pt, phi_pt;
  complex_PRECISION buffer1[12], buffer2[12];
  config_PRECISION D_pt, D = s->op.D;
    
  // add clover term/shift
  if ( g.csw == 0.0 ) {
    spinwise_PRECISION_skalarmultiply( eta1, eta2, phi, s->op.shift, 0, l->inner_vector_size, l );
  } else {
    spin0and1_clover_PRECISION( eta1, phi, s->op.clover, l );
    spin2and3_clover_PRECISION( eta2, phi, s->op.clover, l );
  }
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


void operator_updates_PRECISION( level_struct *l ) {
  
  if ( !l->idle ) {
    if ( l->depth == 0 ) {
      schwarz_PRECISION_setup( &(l->s_PRECISION), &(g.op_double), l );
      if ( g.method >= 4 ) {
        oddeven_free_PRECISION( l );
        oddeven_setup_PRECISION( &(g.op_double), l );
      }
    }
    
    if ( l->level > 0 ) {
#ifndef INTERPOLATION_SETUP_LAYOUT_OPTIMIZED_PRECISION
      coarse_operator_PRECISION_setup( l->is_PRECISION.interpolation, l );
#else
      coarse_operator_PRECISION_setup_vectorized( l->is_PRECISION.operator, l, no_threading );
#endif
      conf_PRECISION_gather( &(l->next_level->s_PRECISION.op), &(l->next_level->op_PRECISION), l->next_level );
    }
    
    if ( l->level > 0 && !l->next_level->idle && l->next_level->level > 0 ) {
      schwarz_PRECISION_boundary_update( &(l->next_level->s_PRECISION), l->next_level );
      if ( g.method >= 4 && g.odd_even ) {
        coarse_oddeven_free_PRECISION( l->next_level );
        coarse_oddeven_setup_PRECISION( &(l->next_level->s_PRECISION.op), _REORDER, l->next_level );
      }
    }
    
    if (  l->level > 0 && !l->next_level->idle && l->next_level->level == 0 && g.odd_even ) {
      coarse_oddeven_free_PRECISION( l->next_level );
      coarse_oddeven_setup_PRECISION( &(l->next_level->s_PRECISION.op), _NO_REORDERING, l->next_level );
    }
    
    if ( l->level > 0 )
      operator_updates_PRECISION( l->next_level );   
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
    START_LOCKED_MASTER(threading)
    op->shift = 4+shift;
    END_LOCKED_MASTER(threading)
  }
}


void g5D_shift_update_PRECISION( operator_PRECISION_struct *op, complex_PRECISION shift, level_struct *l, struct Thread *threading ) {
  
  // no hyperthreading in this function
  if(threading->thread != 0)
    return;

  config_PRECISION clover = op->clover;
  
  if ( clover != NULL ) {
    int i, j;
    complex_PRECISION old_shift = (complex_PRECISION) g.g5D_shift;
    complex_PRECISION shift_diff = shift - old_shift;
    
    if ( l->depth == 0 ) {
      int start = threading->start_site[l->depth];
      int n     = threading->n_site[l->depth];
      clover += start*(g.csw?42:12);
      for ( i=0; i<n; i++ ) {
        for ( j=0; j<6; j++ ) {
          clover[j] -= shift_diff;
        }
        for ( j=6; j<12; j++ ) {
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
          *clover -= shift_diff;
        }
        clover ++;
        for ( j=0; j<k; j++ ) {
          if ( j>0 ) clover += j+1;
          *clover += shift_diff;
        }
        clover += 1 + SQUARE(k);
      }
    }
    
    SYNC_CORES(threading)
  }
}


