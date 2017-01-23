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

void coarse_operator_PRECISION_alloc( level_struct *l ) {
  
  int nd = l->next_level->num_inner_lattice_sites,
      k = l->next_level->num_lattice_site_var;  
  l->next_level->D_size = k*k*4*nd;
  l->next_level->clover_size = ((k*(k+1))/2)*nd;
  
  operator_PRECISION_alloc( &(l->next_level->op_PRECISION), _ORDINARY, l->next_level );
}


void coarse_operator_PRECISION_free( level_struct *l ) {
  
  operator_PRECISION_free( &(l->next_level->op_PRECISION), _ORDINARY, l->next_level );
  
#ifdef VECTORIZE_COARSE_OPERATOR_PRECISION
  operator_PRECISION_struct *op = &(l->next_level->s_PRECISION.op);
  if( op->D_vectorized != NULL ) {
    int n2 = 2*l->next_level->num_lattice_sites-l->next_level->num_inner_lattice_sites, n = l->next_level->num_inner_lattice_sites;
    int column_offset = SIMD_LENGTH_PRECISION*((l->next_level->num_lattice_site_var+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
    // 2 is for complex, 4 is for 4 directions
    FREE_HUGEPAGES( op->D_vectorized, OPERATOR_TYPE_PRECISION, 2*4*l->next_level->num_lattice_site_var*column_offset*n2 );
    FREE_HUGEPAGES( op->D_transformed_vectorized, OPERATOR_TYPE_PRECISION, 2*4*l->next_level->num_lattice_site_var*column_offset*n2 );
    FREE_HUGEPAGES( op->clover_vectorized, OPERATOR_TYPE_PRECISION, 2*l->next_level->num_lattice_site_var*column_offset*n );
  }
#endif
}


void coarse_operator_PRECISION_setup( vector_PRECISION *V, level_struct *l ) {
  
  double t0, t1;
  t0 = MPI_Wtime();
  
  vector_PRECISION buffer1 = l->vbuf_PRECISION[4], buffer2 = l->vbuf_PRECISION[5];
  
  int mu, n = l->num_eig_vect, i, j,
      D_size = l->next_level->D_size,
      clover_size = l->next_level->clover_size;
  void (*aggregate_self_coupling)() = (l->depth==0)?d_plus_clover_aggregate_PRECISION:coarse_aggregate_self_couplings_PRECISION,
       (*aggregate_neighbor_coupling)() = (l->depth==0)?d_neighbor_aggregate_PRECISION:coarse_aggregate_neighbor_couplings_PRECISION;
   
  operator_PRECISION_define( &(l->next_level->op_PRECISION), l->next_level );
    
  for ( j=0; j<D_size; j++ )
    l->next_level->op_PRECISION.D[j] = _COMPLEX_PRECISION_ZERO;
  for ( j=0; j<clover_size; j++ )
    l->next_level->op_PRECISION.clover[j] = _COMPLEX_PRECISION_ZERO;
  
  // for all test vectors V[i]:
  for ( i=0; i<n; i++ ) {
    for ( mu=0; mu<4; mu++ ) {
      // update ghost cells of V[i]
      negative_sendrecv_PRECISION( V[i], mu, &(l->s_PRECISION.op.c), l );
    }
    // apply self coupling of block-and-2spin-restricted dirac operator for each aggregate
    aggregate_self_coupling( buffer1, buffer2, V[i], &(l->s_PRECISION), l );
    // calculate selfcoupling entries of the coarse grid operator
    set_coarse_self_coupling_PRECISION( buffer1, buffer2, V, i, l );
    
    for ( mu=0; mu<4; mu++ ) {
      // finish updating ghostcells of V[i]
      negative_wait_PRECISION( mu, &(l->s_PRECISION.op.c), l );      
      // apply 2spin-restricted dirac operator for direction mu for all aggregates
      aggregate_neighbor_coupling( buffer1, buffer2, V[i], mu, &(l->s_PRECISION), l );      
      set_coarse_neighbor_coupling_PRECISION( buffer1, buffer2, V, mu, i, l );
    }
  }
  
  t1 = MPI_Wtime();
  if ( g.print > 0 ) printf0("depth: %d, time spent for setting up next coarser operator: %lf seconds\n", l->depth, t1-t0 );

}


void set_coarse_self_coupling_PRECISION( vector_PRECISION spin_0_1, vector_PRECISION spin_2_3, 
                                         vector_PRECISION *V, const int n, level_struct *l ) {
  
  int i, j, k, m, k1, k2, num_aggregates = l->is_PRECISION.num_agg,
      num_eig_vect = l->next_level->num_lattice_site_var/2,
      aggregate_size = l->inner_vector_size / num_aggregates, offset = l->num_lattice_site_var/2,
      clover_site_size = (l->next_level->num_lattice_site_var*(l->next_level->num_lattice_site_var+1))/2;
  vector_PRECISION spin_0_1_pt, spin_2_3_pt, interpolation_data;
  config_PRECISION clover_pt, clover = l->next_level->op_PRECISION.clover;  
  
  // U(x) = [ A B      , A=A*, D=D*, C = -B*
  //          C D ]
  // storage order: upper triangle of A, upper triangle of D, B, columnwise
  // diagonal coupling
  for ( k=0; k<=n; k++ ) {
    k1 = (n*(n+1))/2+k; k2 = (n*(n+1))/2+k+(num_eig_vect*(num_eig_vect+1))/2;
  
    for ( j=0; j<num_aggregates; j++ ) {
      spin_0_1_pt = spin_0_1 + j*aggregate_size;
      spin_2_3_pt = spin_2_3 + j*aggregate_size;
      interpolation_data = V[k] + j*aggregate_size;
      clover_pt = clover + j*clover_site_size;
      
      for ( i=0; i<aggregate_size; ) {
        // A
        for ( m=0; m<offset; m++, i++ )
          clover_pt[ k1 ] += conj_PRECISION( interpolation_data[i] ) * spin_0_1_pt[i];
        // D
        for ( m=0; m<offset; m++, i++ )
          clover_pt[ k2 ] += conj_PRECISION( interpolation_data[i] ) * spin_2_3_pt[i];
      }
    }
  }
  
  for ( k=0; k<num_eig_vect; k++ ) {
    k1 = num_eig_vect*(num_eig_vect+1+n) + k;
  
    for ( j=0; j<num_aggregates; j++ ) {
      spin_0_1_pt = spin_0_1 + j*aggregate_size;
      spin_2_3_pt = spin_2_3 + j*aggregate_size;
      interpolation_data = V[k] + j*aggregate_size;
      clover_pt = clover + j*clover_site_size;
      
      for ( i=0; i<aggregate_size; ) {
        // B
        for ( m=0; m<offset; m++, i++ )
          clover_pt[ k1 ] += conj_PRECISION( interpolation_data[i] ) * spin_2_3_pt[i];
        i += offset;
      }
    }
  }
}


void set_coarse_neighbor_coupling_PRECISION( vector_PRECISION spin_0_1, vector_PRECISION spin_2_3, 
                                             vector_PRECISION *V, const int mu, const int n, level_struct *l ) {
  
  int i, i1, j, k, k1, k2, m, num_aggregates = l->is_PRECISION.num_agg,
      num_eig_vect = l->next_level->num_lattice_site_var/2,
      offset = l->num_lattice_site_var/2, nlsv = l->num_lattice_site_var,
      D_link_size = num_eig_vect*num_eig_vect*4, *index_dir = l->is_PRECISION.agg_boundary_index[mu],
      aggregate_boundary_sites = l->is_PRECISION.agg_boundary_length[mu]/num_aggregates;
      
  vector_PRECISION spin_0_1_pt, spin_2_3_pt, interpolation_data;
  config_PRECISION D_pt, D = l->next_level->op_PRECISION.D;
  
  // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
  //             C D ]                        -B*  D* ]
  // storage order: A, C, B, D, each column wise
  for ( k=0; k<num_eig_vect; k++ ) {
    k1 = n*num_eig_vect + k;
    k2 = (n+num_eig_vect)*num_eig_vect + k;
    i1 = 0;
    for ( j=0; j<num_aggregates; j++ ) {
      D_pt = D+(j*4+mu)*D_link_size;
      
      for ( i=0; i<aggregate_boundary_sites; i++ ) {
        spin_0_1_pt = spin_0_1 + nlsv*index_dir[i1];
        interpolation_data = V[k] + nlsv*index_dir[i1]; i1++;
        // A
        for ( m=0; m<offset; m++ )
          D_pt[ k1 ] += conj_PRECISION( interpolation_data[m] ) * spin_0_1_pt[m];
        // C
        for ( ; m<2*offset; m++ )
          D_pt[ k2 ] += conj_PRECISION( interpolation_data[m] ) * spin_0_1_pt[m];
      }
    }
    
    k1 = (n+2*num_eig_vect)*num_eig_vect + k;
    k2 = (n+3*num_eig_vect)*num_eig_vect + k;
    i1 = 0;
    for ( j=0; j<num_aggregates; j++ ) {
      D_pt = D+(j*4+mu)*D_link_size;
      
      for ( i=0; i<aggregate_boundary_sites; i++ ) {
        spin_2_3_pt = spin_2_3 + nlsv*index_dir[i1];
        interpolation_data = V[k] + nlsv*index_dir[i1]; i1++;
        // B
        for ( m=0; m<offset; m++ )
          D_pt[ k1 ] += conj_PRECISION( interpolation_data[m] ) * spin_2_3_pt[m];
        // D
        for ( ; m<2*offset; m++ )
          D_pt[ k2 ] += conj_PRECISION( interpolation_data[m] ) * spin_2_3_pt[m];
      }
    }
  }
}

#ifndef VECTORIZE_COARSE_OPERATOR_PRECISION
void coarse_block_operator_PRECISION( vector_PRECISION eta, vector_PRECISION phi, int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {
  
  START_UNTHREADED_FUNCTION(threading)

  int n = s->num_block_sites, *length = s->dir_length, **index = s->index,
      *ind, *neighbor = s->op.neighbor_table, m = l->num_lattice_site_var;
  vector_PRECISION lphi = phi+start, leta = eta+start;
  int hopp_size = 4 * SQUARE( l->num_lattice_site_var );
  config_PRECISION D_pt, D = s->op.D + (start/m)*hopp_size;
  
  // site-wise self coupling
  int clov_size = ( (l->num_lattice_site_var*(l->num_lattice_site_var+1))/2 );
  config_PRECISION clover = s->op.clover + (start/m)*clov_size;
  coarse_self_couplings_PRECISION( leta, lphi, clover, n*m, l );
  
  // inner block couplings
  for ( int mu=0; mu<4; mu++ ) {
    ind = index[mu]; // mu direction
    for ( int i=0; i<length[mu]; i++ ) {
      int k = ind[i]; int j = neighbor[5*k+mu+1];
      D_pt = D + hopp_size*k + (hopp_size/4)*mu;
      coarse_hopp_PRECISION( leta+m*k, lphi+m*j, D_pt, l );
      coarse_daggered_hopp_PRECISION( leta+m*j, lphi+m*k, D_pt, l );
    }
  }

  END_UNTHREADED_FUNCTION(threading)
}
#endif


void coarse_aggregate_self_couplings_PRECISION( vector_PRECISION eta1, vector_PRECISION eta2, vector_PRECISION phi,
                                                schwarz_PRECISION_struct *s, level_struct *l ) {
  
  int i, mu, index1, index2, length, *index_dir, *neighbor = s->op.neighbor_table,
      n = l->num_lattice_site_var, Dls = n*n, Dss = 4*n*n;
  vector_PRECISION eta1_pt, eta2_pt, phi_pt;
  config_PRECISION D_pt, D = s->op.D;
  
  vector_PRECISION_define( eta1, 0, 0, l->vector_size, l );
  vector_PRECISION_define( eta2, 0, 0, l->vector_size, l );  
  coarse_spinwise_self_couplings_PRECISION( eta1, eta2, phi, s->op.clover, l->inner_vector_size, l );
  
  for ( mu=0; mu<4; mu++ ) { // direction mu
    length = l->is_PRECISION.agg_length[mu]; index_dir = l->is_PRECISION.agg_index[mu];
    for ( i=0; i<length; i++ ) {
      index1 = index_dir[i]; index2 = neighbor[5*index1+mu+1]; D_pt = D + Dss*index1 + Dls*mu;
      phi_pt = phi + n*index2; eta1_pt = eta1 + n*index1; eta2_pt = eta2 + n*index1;
      coarse_spinwise_n_hopp_PRECISION( eta1_pt, eta2_pt, phi_pt, D_pt, l );
      phi_pt = phi + n*index1; eta1_pt = eta1 + n*index2; eta2_pt = eta2 + n*index2;
      coarse_spinwise_n_daggered_hopp_PRECISION( eta1_pt, eta2_pt, phi_pt, D_pt, l );
    }
  }
}


void coarse_aggregate_neighbor_couplings_PRECISION( vector_PRECISION eta1, vector_PRECISION eta2, vector_PRECISION phi,
                                                    const int mu, schwarz_PRECISION_struct *s, level_struct *l ) {
  
  int i, index1, index2, length = l->is_PRECISION.agg_boundary_length[mu],
      *index_dir = l->is_PRECISION.agg_boundary_index[mu],
      *neighbor = l->is_PRECISION.agg_boundary_neighbor[mu],
      n = l->num_lattice_site_var, Dls = n*n, Dss = 4*n*n;
  vector_PRECISION eta1_pt, eta2_pt, phi_pt;
  config_PRECISION D_pt, D = s->op.D;
  
  vector_PRECISION_define( eta1, 0, 0, l->vector_size, l );
  vector_PRECISION_define( eta2, 0, 0, l->vector_size, l ); 
  
  // requires the positive boundaries of phi to be communicated befor
  for ( i=0; i<length; i++ ) {
    index1 = index_dir[i];
    index2 = neighbor[i];
    D_pt = D + Dss*index1 + Dls*mu;
    phi_pt = phi + n*index2; eta1_pt = eta1 + n*index1; eta2_pt = eta2 + n*index1;
    coarse_spinwise_hopp_PRECISION( eta1_pt, eta2_pt, phi_pt, D_pt, l );
  }
}


void coarse_self_couplings_PRECISION( vector_PRECISION eta, vector_PRECISION phi,
                                      config_PRECISION clover, int length, level_struct *l ) {

  int site_var = l->num_lattice_site_var,
      num_eig_vect = l->num_lattice_site_var/2,
      clover_step_size1 = (num_eig_vect * (num_eig_vect+1))/2,
      clover_step_size2 = SQUARE(l->num_lattice_site_var/2);
  config_PRECISION clover_pt = clover;
  vector_PRECISION phi_pt=phi, eta_pt=eta, phi_end_pt=phi+length;
  // U(x) = [ A B      , A=A*, D=D*, C = -B*
  //          C D ]
  // storage order: upper triangle of A, upper triangle of D, B, columnwise
  // diagonal coupling
  while ( phi_pt < phi_end_pt ) {
    // A
    mvp_PRECISION( eta_pt, clover_pt, phi_pt, num_eig_vect );
    clover_pt += clover_step_size1; eta_pt += num_eig_vect; phi_pt += num_eig_vect;
    // D
    mvp_PRECISION( eta_pt, clover_pt, phi_pt, num_eig_vect );
    clover_pt += clover_step_size1; phi_pt -= num_eig_vect;
    // C = -B*
    nmvh_PRECISION( eta_pt, clover_pt, phi_pt, num_eig_vect );
    phi_pt += num_eig_vect; eta_pt -= num_eig_vect;
    // B
    mv_PRECISION( eta_pt, clover_pt, phi_pt, num_eig_vect );
    clover_pt += clover_step_size2; phi_pt += num_eig_vect; eta_pt += site_var;
  }
}


void coarse_spinwise_self_couplings_PRECISION( vector_PRECISION eta1, vector_PRECISION eta2, vector_PRECISION phi,
                                               config_PRECISION clover, int length, level_struct *l ) {
  
  int num_eig_vect = l->num_lattice_site_var/2,
      clover_step_size1 = (num_eig_vect * (num_eig_vect+1))/2,
      clover_step_size2 = SQUARE(l->num_lattice_site_var/2);
  config_PRECISION clover_pt = clover;
  vector_PRECISION phi_pt=phi, eta1_pt=eta1, eta2_pt=eta2+num_eig_vect, phi_end_pt=phi+length;
  // U(x) = [ A B      , A=A*, D=D*, C = -B*
  //          C D ]
  // storage order: upper triangle of A, upper triangle of D, B, columnwise
  // diagonal coupling
  while ( phi_pt < phi_end_pt ) {
    // A
    mvp_PRECISION( eta1_pt, clover_pt, phi_pt, num_eig_vect );
    clover_pt += clover_step_size1; phi_pt += num_eig_vect; eta1_pt += num_eig_vect; 
    // D
    mvp_PRECISION( eta2_pt, clover_pt, phi_pt, num_eig_vect );
    clover_pt += clover_step_size1; phi_pt -= num_eig_vect; eta2_pt -= num_eig_vect; 
    // C = -B*
    nmvh_PRECISION( eta1_pt, clover_pt, phi_pt, num_eig_vect );
    phi_pt += num_eig_vect; eta1_pt += num_eig_vect;
    // B
    mv_PRECISION( eta2_pt, clover_pt, phi_pt, num_eig_vect );
    clover_pt += clover_step_size2; phi_pt += num_eig_vect; eta2_pt += 3*num_eig_vect;
  }
}


#ifndef VECTORIZE_COARSE_OPERATOR_PRECISION
void coarse_operator_PRECISION_set_couplings( operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) { }
#endif

void coarse_gamma5_PRECISION( vector_PRECISION eta, vector_PRECISION phi, int start, int end, level_struct *l ) {
  
  int j, k=l->num_lattice_site_var/2;
  vector_PRECISION eta_end;
  
  eta_end = eta+end;
  phi += start;
  eta += start;
  
  if ( eta != phi ) {
    while ( eta < eta_end ) {
      for ( j=0; j<k; j++ ) {
        *eta = -(*phi);
        eta++; phi++;
      }
      for ( j=0; j<k; j++ ) {
        *eta = *phi;
        eta++; phi++;
      }
    }
  } else {
    while ( eta < eta_end ) {
      for ( j=0; j<k; j++ ) {
        *eta = -(*eta);
        eta++;
      }
      eta+=k;
    }
  }
}

#ifndef VECTORIZE_COARSE_OPERATOR_PRECISION
void apply_coarse_operator_PRECISION( vector_PRECISION eta, vector_PRECISION phi, operator_PRECISION_struct *op,
                                      level_struct *l, struct Thread *threading ) {
  
  PROF_PRECISION_START( _SC, threading );
  START_LOCKED_MASTER(threading)
  coarse_self_couplings_PRECISION( eta, phi, op->clover, l->inner_vector_size, l );
  END_LOCKED_MASTER(threading)
  PROF_PRECISION_STOP( _SC, 1, threading );
  PROF_PRECISION_START( _NC, threading );
  coarse_hopping_term_PRECISION( eta, phi, op, _FULL_SYSTEM, l, threading );
  PROF_PRECISION_STOP( _NC, 1, threading );
}
#endif

void g5D_apply_coarse_operator_PRECISION( vector_PRECISION eta, vector_PRECISION phi, operator_PRECISION_struct *op,
                                          level_struct *l, struct Thread *threading ) {
  int start, end;
  compute_core_start_end_custom(0, l->inner_vector_size, &start, &end, l, threading, l->num_lattice_site_var );
  apply_coarse_operator_PRECISION( eta, phi, op, l, threading );
  SYNC_CORES(threading)
  coarse_gamma5_PRECISION( eta, eta, start, end, l );
  SYNC_CORES(threading)
}


void apply_coarse_operator_dagger_PRECISION( vector_PRECISION eta, vector_PRECISION phi, operator_PRECISION_struct *op,
                                             level_struct *l, struct Thread *threading ) {
  
  coarse_gamma5_PRECISION( l->vbuf_PRECISION[3], phi, threading->start_index[l->depth], threading->end_index[l->depth], l );
  apply_coarse_operator_PRECISION( eta, l->vbuf_PRECISION[3], op, l, threading );
  coarse_gamma5_PRECISION( eta, eta, threading->start_index[l->depth], threading->end_index[l->depth], l );
}


void coarse_operator_PRECISION_test_routine( level_struct *l, struct Thread *threading ) {

  if ( !l->idle ) {
    int vs = l->vector_size, ivs = l->inner_vector_size,
        cvs = l->next_level->vector_size, civs = l->next_level->inner_vector_size;
    double diff;
    vector_PRECISION vp1=NULL, vp2, vp3, vp4, vc1=NULL, vc2, vc3;

    PUBLIC_MALLOC( vp1, complex_PRECISION, 4*vs );
    PUBLIC_MALLOC( vc1, complex_PRECISION, 3*cvs );
    
    SYNC_MASTER_TO_ALL(threading)
    
    vp2 = vp1 + vs; vp3 = vp2 + vs; vp4 = vp3 + vs; vc2 = vc1 + cvs; vc3 = vc2 + cvs; 

    START_LOCKED_MASTER(threading)
#ifdef INTERPOLATION_OPERATOR_LAYOUT_OPTIMIZED_PRECISION
    double norm = 0.0;
    double dot = 0.0;
    float *op = (float *)l->is_PRECISION.operator;
    float *op2 = (float *)(l->is_PRECISION.operator+0*SIMD_LENGTH_PRECISION*l->vector_size)+1;
    for ( int i=0; i<l->inner_vector_size; i++ )
      norm += (op[2*i*SIMD_LENGTH_PRECISION+0] + I*op[2*i*SIMD_LENGTH_PRECISION+SIMD_LENGTH_PRECISION])*conj(op[2*i*SIMD_LENGTH_PRECISION+0] + I*op[2*i*SIMD_LENGTH_PRECISION+SIMD_LENGTH_PRECISION]);
    for ( int i=0; i<l->inner_vector_size; i++ )
      dot += (op[2*i*SIMD_LENGTH_PRECISION+0] + I*op[2*i*SIMD_LENGTH_PRECISION+SIMD_LENGTH_PRECISION])*conj(op2[2*i*SIMD_LENGTH_PRECISION+0] + I*op2[2*i*SIMD_LENGTH_PRECISION+SIMD_LENGTH_PRECISION]);
    diff = dot/norm;
#else
    diff = global_inner_product_PRECISION( l->is_PRECISION.interpolation[0], l->is_PRECISION.interpolation[1], 0, ivs, l, no_threading )
         / global_norm_PRECISION( l->is_PRECISION.interpolation[0], 0, ivs, l, no_threading );
#endif
    printf0("depth: %d, correctness of block_gram_schmidt: %le\n", l->depth, cabs(diff) );
    
    if ( !l->next_level->idle )
      vector_PRECISION_define_random( vc1, 0, civs, l->next_level );
    vector_PRECISION_distribute( vc2, vc1, l->next_level );
    vector_PRECISION_gather( vc3, vc2, l->next_level );
    if ( !l->next_level->idle ) {
      vector_PRECISION_minus( vc2, vc1, vc3, 0, civs, l->next_level );
      diff = global_norm_PRECISION( vc2, 0, civs, l->next_level, no_threading ) / global_norm_PRECISION( vc1, 0, civs, l->next_level, no_threading );
    }
    printf0("depth: %d, correctness of gather( distribute( phi_c ) ) : %le\n", l->depth, diff );
    
    if ( !l->next_level->idle )
      vector_PRECISION_define_random( vc1, 0, civs, l->next_level );
    interpolate3_PRECISION( vp1, vc1, l, no_threading );
    restrict_PRECISION( vc2, vp1, l, no_threading );
    if ( !l->next_level->idle ) {
      vector_PRECISION_minus( vc3, vc1, vc2, 0, civs, l->next_level );
      diff = global_norm_PRECISION( vc3, 0, civs, l->next_level, no_threading ) / global_norm_PRECISION( vc1, 0, civs, l->next_level, no_threading );
      printf0("depth: %d, correctness of ( P* P - 1 ) phi_c: %le\n", l->depth, diff );
    }    

    END_LOCKED_MASTER(threading)
    if(threading->n_core>1) {
      interpolate3_PRECISION( vp1, vc1, l, threading );
      restrict_PRECISION( vc2, vp1, l, threading );
      START_LOCKED_MASTER(threading)
      if ( !l->next_level->idle ) {
        vector_PRECISION_minus( vc3, vc1, vc2, 0, civs, l->next_level );
        diff = global_norm_PRECISION( vc3, 0, civs, l->next_level, no_threading ) / global_norm_PRECISION( vc1, 0, civs, l->next_level, no_threading );
        printf0("depth: %d, correctness of ( P* P - 1 ) phi_c with threading: %le\n", l->depth, diff );
      }
      END_LOCKED_MASTER(threading)
    }
    
    if ( l->level > 0 ) {
      START_LOCKED_MASTER(threading)
      interpolate3_PRECISION( vp1, vc1, l, no_threading );
      apply_operator_PRECISION( vp2, vp1, &(l->p_PRECISION), l, no_threading );      
      restrict_PRECISION( vc2, vp2, l, no_threading );
      if ( !l->next_level->idle ) {
        if ( l->level==1 && g.odd_even )
          coarse_odd_even_PRECISION_test( vc3, vc1, l->next_level, no_threading );
        else
          apply_operator_PRECISION( vc3, vc1, &(l->next_level->p_PRECISION), l->next_level, no_threading );

        vector_PRECISION_minus( vc3, vc2, vc3, 0, civs, l->next_level );
        diff = global_norm_PRECISION( vc3, 0, civs, l->next_level, no_threading ) / global_norm_PRECISION( vc2, 0, civs, l->next_level, no_threading );
        if ( l->level==1 && g.odd_even )
          printf0("depth: %d, correctness of odd even preconditioned ( P* D P - D_c ) phi_c: %le\n", l->depth, diff );
        else
          printf0("depth: %d, correctness of ( P* D P - D_c ) phi_c: %le\n", l->depth, diff );
      }      
      END_LOCKED_MASTER(threading)

      if(threading->n_core>1) {
        if ( !l->next_level->idle ) {
          if ( l->level==1 && g.odd_even )
            coarse_odd_even_PRECISION_test( vc3, vc1, l->next_level, threading );
          else
            apply_operator_PRECISION( vc3, vc1, &(l->next_level->p_PRECISION), l->next_level, threading );
        }
        START_LOCKED_MASTER(threading)
        if ( !l->next_level->idle ) {
          vector_PRECISION_minus( vc3, vc2, vc3, 0, civs, l->next_level );
          diff = global_norm_PRECISION( vc3, 0, civs, l->next_level, no_threading ) / global_norm_PRECISION( vc2, 0, civs, l->next_level, no_threading );
          if ( l->level==1 && g.odd_even )
            printf0("depth: %d, correctness of odd even preconditioned ( P* D P - D_c ) phi_c with D_c threaded: %le\n", l->depth, diff );
          else
            printf0("depth: %d, correctness of ( P* D P - D_c ) phi_c with D_c threaded: %le\n", l->depth, diff );
        }
        END_LOCKED_MASTER(threading)
      }
    }
    START_LOCKED_MASTER(threading)


    if ( l->level > 0 && l->depth > 0 && g.method == 3 && g.odd_even ) {
      vector_PRECISION_define_random( vp1, 0, ivs, l );
      block_to_oddeven_PRECISION( vp4, vp1, l, no_threading );
      coarse_diag_ee_PRECISION( vp3, vp4, &(l->oe_op_PRECISION), l, no_threading );
      coarse_diag_oo_PRECISION( vp3, vp4, &(l->oe_op_PRECISION), l, no_threading );
      coarse_hopping_term_PRECISION( vp3, vp4, &(l->oe_op_PRECISION), _FULL_SYSTEM, l, no_threading );
      oddeven_to_block_PRECISION( vp4, vp3, l, no_threading );
      apply_operator_PRECISION( vp2, vp1, &(l->p_PRECISION), l, no_threading );
      vector_PRECISION_minus( vp4, vp4, vp2, 0, ivs, l );
      diff = global_norm_PRECISION( vp4, 0, ivs, l, no_threading ) / global_norm_PRECISION( vp2, 0, ivs, l, no_threading );
      printf0("depth: %d, correctness of odd even layout (smoother): %le\n", l->depth, diff );
      
      block_to_oddeven_PRECISION( vp4, vp1, l, no_threading );
      coarse_odd_even_PRECISION_test( vp3, vp4, l, no_threading );
      oddeven_to_block_PRECISION( vp4, vp3, l, no_threading );
      apply_operator_PRECISION( vp2, vp1, &(l->p_PRECISION), l, no_threading );
      vector_PRECISION_minus( vp4, vp4, vp2, 0, ivs, l );
      diff = global_norm_PRECISION( vp4, 0, ivs, l, no_threading ) / global_norm_PRECISION( vp2, 0, ivs, l, no_threading );
      printf0("depth: %d, correctness of odd even preconditioned operator (smoother): %le\n", l->depth, diff );
    }
    
    FREE( vp1, complex_PRECISION, 4*vs );
    FREE( vc1, complex_PRECISION, 3*cvs );
    END_LOCKED_MASTER(threading)
    
    if ( g.method != 6 && l->next_level->level > 0  && !l->next_level->idle ) {
      schwarz_PRECISION_mvm_testfun( &(l->next_level->s_PRECISION), l->next_level, threading );
    }

    if ( l->next_level->level > 0 && !l->next_level->idle )
      coarse_operator_PRECISION_test_routine( l->next_level, threading );
    
    SYNC_CORES(threading)
    SYNC_MASTER_TO_ALL(threading)
  }
}
