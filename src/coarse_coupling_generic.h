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

#ifndef COARSE_COUPLING_PRECISION_HEADER
  #define COARSE_COUPLING_PRECISION_HEADER

  #undef COMM_HIDING_COARSEOP

  struct Thread;

  // eta +/-= D*phi, D stored columnwise
  static inline void pnmv_PRECISION( const vector_PRECISION eta, const complex_PRECISION *D,
                                     const vector_PRECISION phi, const register int n,
                                     const int sign ) {
    register int i, j, k=0;
    
    for ( i=0; i<n; i++ )
      for ( j=0; j<n; j++, k++ )
        eta[j] += sign*D[k]*phi[i];
  }

  // eta +/-= D^Dagger*phi, D stored columnwise
  static inline void pnmvh_PRECISION( const vector_PRECISION eta, const complex_PRECISION *D,
                                      const vector_PRECISION phi, const register int n,
                                      const int sign ) {
    register int i, j, k=0;    

    for ( i=0; i<n; i++ )
      for ( j=0; j<n; j++, k++ )
        eta[i] += sign*conj_PRECISION(D[k])*phi[j];
  }

  static inline void coarse_pn_hopp_PRECISION( vector_PRECISION eta, vector_PRECISION phi,
                                               config_PRECISION D, const int sign,
                                               level_struct *l ) {
  
    int num_eig_vect = l->num_parent_eig_vect,
        num_eig_vect2 = SQUARE(l->num_parent_eig_vect);
    // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
    //             C D ]                        -B*  D* ]
    // storage order: A, C, B, D
    // note: minus sign of D = self_coupling - hopping_term is added here

#ifdef HAVE_TM1p1
    if( g.n_flavours == 2 ) {
      // A  
      pnmv_PRECISION( eta, D, phi, num_eig_vect, -sign );
      eta += num_eig_vect;//1
      phi += num_eig_vect;//1
      pnmv_PRECISION( eta, D, phi, num_eig_vect, -sign );
      // C
      eta += num_eig_vect;//2
      phi -= num_eig_vect;//0
      D += num_eig_vect2;
      pnmv_PRECISION( eta, D, phi, num_eig_vect, -sign );
      eta += num_eig_vect;//3
      phi += num_eig_vect;//1
      pnmv_PRECISION( eta, D, phi, num_eig_vect, -sign );
      // B
      eta -= 3*num_eig_vect;//0
      phi += num_eig_vect;//2
      D += num_eig_vect2;
      pnmv_PRECISION( eta, D, phi, num_eig_vect, -sign );
      eta += num_eig_vect;//1
      phi += num_eig_vect;//3
      pnmv_PRECISION( eta, D, phi, num_eig_vect, -sign );
      // D
      eta += num_eig_vect;//2
      phi -= num_eig_vect;//2
      D += num_eig_vect2;
      pnmv_PRECISION( eta, D, phi, num_eig_vect, -sign );
      eta += num_eig_vect;//3
      phi += num_eig_vect;//3
      pnmv_PRECISION( eta, D, phi, num_eig_vect, -sign );
    } else {
#endif
      // A  
      pnmv_PRECISION( eta, D, phi, num_eig_vect, -sign );
      // C
      eta += num_eig_vect;
      D += num_eig_vect2;
      pnmv_PRECISION( eta, D, phi, num_eig_vect, -sign );
      // B
      phi += num_eig_vect;
      eta -= num_eig_vect;
      D += num_eig_vect2;
      pnmv_PRECISION( eta, D, phi, num_eig_vect, -sign );
      // D
      eta += num_eig_vect;
      D += num_eig_vect2;
      pnmv_PRECISION( eta, D, phi, num_eig_vect, -sign );
#ifdef HAVE_TM1p1
    }
#endif
  }

  static inline void coarse_pn_daggered_hopp_PRECISION( vector_PRECISION eta, vector_PRECISION phi,
                                                        config_PRECISION D, const int sign,
                                                        level_struct *l ) {
    
    int num_eig_vect = l->num_parent_eig_vect,
        num_eig_vect2 = SQUARE(l->num_parent_eig_vect);
    // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
    //             C D ]                        -B*  D* ]
    // storage order: A, C, B, D
    // note: minus sign of D = self_coupling - hopping_term is added here

#ifdef HAVE_TM1p1
    if( g.n_flavours == 2 ) {
      // A* 
      pnmvh_PRECISION( eta, D, phi, num_eig_vect, -sign );
      eta += num_eig_vect;//1
      phi += num_eig_vect;//1
      pnmvh_PRECISION( eta, D, phi, num_eig_vect, -sign );
      // -C*
      eta -= num_eig_vect;//0
      phi += num_eig_vect;//2
      D += num_eig_vect2;
      pnmvh_PRECISION( eta, D, phi, num_eig_vect, sign );
      eta += num_eig_vect;//1
      phi += num_eig_vect;//3
      pnmvh_PRECISION( eta, D, phi, num_eig_vect, sign );
      // -B*
      eta += num_eig_vect;//2
      phi -= 3*num_eig_vect;//0
      D += num_eig_vect2;
      pnmvh_PRECISION( eta, D, phi, num_eig_vect, sign );
      eta += num_eig_vect;//3
      phi += num_eig_vect;//1
      pnmvh_PRECISION( eta, D, phi, num_eig_vect, sign );
      // D*
      eta -= num_eig_vect;//2
      phi += num_eig_vect;//2
      D += num_eig_vect2;
      pnmvh_PRECISION( eta, D, phi, num_eig_vect, -sign );
      eta += num_eig_vect;//3
      phi += num_eig_vect;//3
      pnmvh_PRECISION( eta, D, phi, num_eig_vect, -sign );
    } else {
#endif
      // A* 
      pnmvh_PRECISION( eta, D, phi, num_eig_vect, -sign );
      // -C*
      phi += num_eig_vect;
      D += num_eig_vect2;
      pnmvh_PRECISION( eta, D, phi, num_eig_vect, sign );
      // -B*
      eta += num_eig_vect;
      phi -= num_eig_vect;
      D += num_eig_vect2;
      pnmvh_PRECISION( eta, D, phi, num_eig_vect, sign );
      // D*
      phi += num_eig_vect;
      D += num_eig_vect2;
      pnmvh_PRECISION( eta, D, phi, num_eig_vect, -sign );
#ifdef HAVE_TM1p1
    }
#endif
  }

  static inline void coarse_pn_hopp_PRECISION_vectorized( vector_PRECISION eta, vector_PRECISION phi, OPERATOR_TYPE_PRECISION *D, const int sign, level_struct *l ) {
#ifdef OPTIMIZED_COARSE_NEIGHBOR_COUPLING_PRECISION
    int nv = l->num_parent_eig_vect;
    int lda = 2*SIMD_LENGTH_PRECISION*((nv+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
    switch (sign) {
    case -1:
      cgemv_padded_PRECISION( 2*nv, D, lda, nv, (float *)phi, (float *)eta);
      break;
    case +1:
    default:
      cgenmv_padded_PRECISION( 2*nv, D, lda, nv, (float *)phi, (float *)eta);
      break;
    }
#endif
  }

  static inline void coarse_pn_hopping_term_PRECISION( vector_PRECISION out, vector_PRECISION in, operator_PRECISION_struct *op, const int amount, const int sign, level_struct *l, struct Thread *threading ) {

    START_NO_HYPERTHREADS(threading)
      
    int mu, i, num_site_var=l->num_lattice_site_var,
      num_eig_vect = l->num_parent_eig_vect,
      num_lattice_sites, start, end, core_start, core_end,
      plus_dir_param, minus_dir_param;
    vector_PRECISION in_pt, out_pt;
    
#ifndef OPTIMIZED_COARSE_NEIGHBOR_COUPLING_PRECISION
    int num_link_var = SQUARE(2*num_eig_vect),
      num_4link_var = 4*num_link_var;
    config_PRECISION D_pt;
    // dagger applied by functions daggered_hopp below
    config_PRECISION D = op->D, D_dagger = op->D;
#else
    int column_offset = 2*SIMD_LENGTH_PRECISION*((num_eig_vect+SIMD_LENGTH_PRECISION-1)/
                                                 SIMD_LENGTH_PRECISION),
      num_link_var = 2*2*num_eig_vect*column_offset,
      num_4link_var = 4*num_link_var;
    OPERATOR_TYPE_PRECISION *D_pt;
    // dagger applied in D_dagger
    OPERATOR_TYPE_PRECISION *D = op->D_vectorized, *D_dagger = op->D_transformed_vectorized;
#endif
    
#ifndef COMM_HIDING_COARSEOP
    int communicate = ( l->num_processes > 1 && op->c.comm ) ? 1:0;
    int *neighbor_fw = op->neighbor_table;
    int *neighbor_bw = op->backward_neighbor_table;
#else
    int communicate = ( op->c.comm ) ? 1:0;
    int *neighbor_fw = op->neighbor_table;
#endif
    
    switch (amount) {
    case _EVEN_SITES:
      minus_dir_param = _ODD_SITES;
      plus_dir_param = _EVEN_SITES;
      break;
    case _ODD_SITES:
      minus_dir_param = _EVEN_SITES;
      plus_dir_param = _ODD_SITES;
      break;
    case _FULL_SYSTEM:
    default:
      minus_dir_param = _FULL_SYSTEM;
      plus_dir_param = _FULL_SYSTEM;
      break;
    }
    
    // assumptions (1) self coupling has already been performed
    //          OR (2) "out" is initialized with zeros
    set_boundary_PRECISION( out, 0, l, threading );
    
    // communicate in -mu direction
    MASTER(threading)
      if ( communicate ) 
        for ( mu=0; mu<4; mu++ ) 
          ghost_sendrecv_PRECISION( in, mu, -1, &(op->c), minus_dir_param, l );
    SYNC_CORES(threading);
    
    switch (amount) {
    case _EVEN_SITES:
      start = op->num_even_sites;
      num_lattice_sites = op->num_odd_sites;
      break;
    case _ODD_SITES:
      start = 0;
      num_lattice_sites = op->num_even_sites;
      break;
    case _FULL_SYSTEM:
    default:
      start=0;
      num_lattice_sites=l->num_inner_lattice_sites;
      break;
    }
    end = start + num_lattice_sites;
    compute_core_start_end_custom( start, end, &core_start, &core_end, l, threading, 1 );
    
#ifndef COMM_HIDING_COARSEOP
    // prepare for sending to fw: compute hopping terms into forward boundary buffer
    if ( communicate ) 
      for ( i=core_start; i<core_end; i++ ) {
        in_pt = in + num_site_var*neighbor_fw[5*i];
        D_pt = D_dagger + num_4link_var*neighbor_fw[5*i];
        for ( mu=0; mu<4; mu++ ) {
          if(neighbor_fw[5*i+1+mu] < l->num_inner_lattice_sites) //num_lattice_sites?
            continue;
          out_pt = out + num_site_var*neighbor_fw[5*i+1+mu];
#ifndef OPTIMIZED_COARSE_NEIGHBOR_COUPLING_PRECISION
          coarse_pn_daggered_hopp_PRECISION( out_pt, in_pt, D_pt+mu*num_link_var, sign, l );
#else
          coarse_pn_hopp_PRECISION_vectorized( out_pt, in_pt, D_pt+mu*num_link_var, sign, l );
#endif
        }
      }
#else
    // compute U_mu^dagger coupling
    for ( mu=0; mu<4; mu++ ) {
      for ( i=core_start; i<core_end; i++ ) {
        in_pt = in + num_site_var*neighbor_fw[5*i];
        D_pt = D_dagger + num_4link_var*neighbor_fw[5*i] + mu*num_link_var;
        out_pt = out + num_site_var*neighbor_fw[5*i+1+mu];
#ifndef OPTIMIZED_COARSE_NEIGHBOR_COUPLING_PRECISION
        coarse_pn_daggered_hopp_PRECISION( out_pt, in_pt, D_pt, sign, l );
#else
        coarse_pn_hopp_PRECISION_vectorized( out_pt, in_pt, D_pt, sign, l );
#endif
      }
    SYNC_CORES(threading);
    }
#endif
    
    if ( communicate ) {
      START_LOCKED_MASTER(threading)
        for ( mu=0; mu<4; mu++ ) {
          // communicate in +mu direction
          ghost_sendrecv_PRECISION( out, mu, +1, &(op->c), plus_dir_param, l );
        }
      for ( mu=0; mu<4; mu++ ) {
        // wait for -mu direction
        ghost_wait_PRECISION( in, mu, -1, &(op->c), minus_dir_param, l );
      }
      END_LOCKED_MASTER(threading);
    }
    else
      SYNC_CORES(threading);
    
    switch (amount) {
    case _EVEN_SITES:
      start = 0;
      num_lattice_sites = op->num_even_sites;
      break;
    case _ODD_SITES:
      start = op->num_even_sites;
      num_lattice_sites = op->num_odd_sites;
      break;
    case _FULL_SYSTEM:
    default:
      start=0;
      num_lattice_sites=l->num_inner_lattice_sites;
      break;
    }
    end = start + num_lattice_sites;
    compute_core_start_end_custom( start, end, &core_start, &core_end, l, threading, 1 );
  
#ifndef COMM_HIDING_COARSEOP
    for ( i=core_start; i<core_end; i++ ) {
      out_pt = out + num_site_var*neighbor_fw[5*i];
      
      // compute U_mu^dagger coupling
      for( mu=0; mu<4; mu++ ) {
        // terms coming from backward boundary buffer are done by the ghost_wait_PRECISION below
        if(neighbor_bw[5*i+1+mu] >= l->num_inner_lattice_sites)
          continue;
        in_pt = in + num_site_var*neighbor_bw[5*i+1+mu];
        D_pt = D_dagger + num_4link_var*neighbor_bw[5*i+1+mu] + mu*num_link_var;
#ifndef OPTIMIZED_COARSE_NEIGHBOR_COUPLING_PRECISION
        coarse_pn_daggered_hopp_PRECISION( out_pt, in_pt, D_pt, sign, l );
#else
        coarse_pn_hopp_PRECISION_vectorized( out_pt, in_pt, D_pt, sign, l );
#endif
      }
      
      // compute U_mu couplings
      D_pt = D + num_4link_var*neighbor_fw[5*i];
      for( mu=0; mu<4; mu++ ) {
        in_pt = in + num_site_var*neighbor_fw[5*i+1+mu];
#ifndef OPTIMIZED_COARSE_NEIGHBOR_COUPLING_PRECISION
        coarse_pn_hopp_PRECISION( out_pt, in_pt, D_pt+mu*num_link_var, sign, l );
#else
        coarse_pn_hopp_PRECISION_vectorized( out_pt, in_pt, D_pt+mu*num_link_var, sign, l );
#endif
      }
    }
#else
    // compute U_mu couplings
    for ( i=core_start; i<core_end; i++ ) {
      out_pt = out + num_site_var*neighbor_fw[5*i];
      D_pt = D + num_4link_var*neighbor_fw[5*i];
      for( mu=0; mu<4; mu++ ) {
        in_pt = in + num_site_var*neighbor_fw[5*i+1+mu];    
#ifndef OPTIMIZED_COARSE_NEIGHBOR_COUPLING_PRECISION
        coarse_pn_hopp_PRECISION( out_pt, in_pt, D_pt+mu*num_link_var, sign, l );
#else
        coarse_pn_hopp_PRECISION_vectorized( out_pt, in_pt, D_pt+mu*num_link_var, sign, l );
#endif
      }
    }
#endif
    
    // wait for terms from bw and add them
    if ( communicate ) {
      START_LOCKED_MASTER(threading)
        for ( mu=0; mu<4; mu++ ) {
          // wait for +mu direction
          ghost_wait_PRECISION( out, mu, +1, &(op->c), plus_dir_param, l );
        }
      END_LOCKED_MASTER(threading);
    }
    else
      SYNC_CORES(threading);
    
    END_NO_HYPERTHREADS(threading);
  }

#endif
