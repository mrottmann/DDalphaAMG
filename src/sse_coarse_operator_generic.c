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

#ifdef SSE

#include "sse_coarse_operator.h"

#ifdef OPTIMIZED_INTERPOLATION_SETUP_PRECISION
void coarse_operator_PRECISION_setup_vectorized( complex_PRECISION *operator, level_struct *l, struct Thread *threading ) {
  
  SYNC_HYPERTHREADS(threading)
  SYNC_CORES(threading)

  double t0, t1;
  t0 = MPI_Wtime();

  int mu, n = l->num_eig_vect, j, num_aggregates = l->is_PRECISION.num_agg,
      aggregate_sites = l->num_inner_lattice_sites / num_aggregates,
      clover_site_size = (l->next_level->num_lattice_site_var*(l->next_level->num_lattice_site_var+1))/2,
      block_site_size = (l->next_level->num_parent_eig_vect*(l->next_level->num_parent_eig_vect+1)),
      D_link_size = 4*l->num_eig_vect*l->num_eig_vect*4,  // size of links in all 4 directions
      fine_components = l->num_lattice_site_var;



  START_LOCKED_MASTER(threading)
  operator_PRECISION_define( &(l->next_level->op_PRECISION), l->next_level );
  END_LOCKED_MASTER(threading)
  SYNC_HYPERTHREADS(threading)

  // each thread loops overs its aggregates and then over internal d.o.f.
  for ( int a=threading->n_core*threading->thread+threading->core; a<num_aggregates; a+=threading->n_core*threading->n_thread ) {
    for ( j=0; j<D_link_size; j++ )
      l->next_level->op_PRECISION.D[j+a*D_link_size] = _COMPLEX_PRECISION_ZERO;
    for ( j=0; j<clover_site_size; j++ )
      l->next_level->op_PRECISION.clover[j+a*clover_site_size] = _COMPLEX_PRECISION_ZERO;
    for ( j=0; j<block_site_size; j++ )
      l->next_level->op_PRECISION.odd_proj[j+a*block_site_size] = _COMPLEX_PRECISION_ZERO;
  }

  complex_PRECISION *mpi_buffer = NULL;
  START_MASTER(threading)
  MALLOC_HUGEPAGES( mpi_buffer, complex_PRECISION, OPERATOR_COMPONENT_OFFSET_PRECISION*(l->vector_size-l->inner_vector_size), 64 );
  END_MASTER(threading)

  int direction_flags[8*l->block_lattice[T]*l->block_lattice[Z]*l->block_lattice[Y]*l->block_lattice[X]];

  // set up table for direction flags
  int *flags = direction_flags;
  if(l->depth == 0) {
    // even sites
    for(int t=0; t < l->block_lattice[T]; t++) {
      for(int z=0; z < l->block_lattice[Z]; z++) {
        for(int y=0; y < l->block_lattice[Y]; y++) {
          for(int x=0; x < l->block_lattice[X]/2; x++) {
            flags[2*X+0] = 1;
            flags[2*X+1] = 1;
            if((y+z+t)%2 == 0) {
              if(x == 0)
                flags[2*X+0] = 0;
            } else {
              if(x == l->block_lattice[X]/2-1)
                flags[2*X+1] = 0;
            }
            flags[2*Y+0] = (y ==                     0)?0:1;
            flags[2*Y+1] = (y == l->block_lattice[Y]-1)?0:1;
            flags[2*Z+0] = (z ==                     0)?0:1;
            flags[2*Z+1] = (z == l->block_lattice[Z]-1)?0:1;
            flags[2*T+0] = (t ==                     0)?0:1;
            flags[2*T+1] = (t == l->block_lattice[T]-1)?0:1;
            flags += 8;
          }
        }
      }
    }
    // odd sites
    for(int t=0; t < l->block_lattice[T]; t++) {
      for(int z=0; z < l->block_lattice[Z]; z++) {
        for(int y=0; y < l->block_lattice[Y]; y++) {
          for(int x=0; x < l->block_lattice[X]/2; x++) {
            flags[2*X+0] = 1;
            flags[2*X+1] = 1;
            if((y+z+t)%2 == 1) {
              if(x == 0)
                flags[2*X+0] = 0;
            } else {
              if(x == l->block_lattice[X]/2-1)
                flags[2*X+1] = 0;
            }
            flags[2*Y+0] = (y ==                     0)?0:1;
            flags[2*Y+1] = (y == l->block_lattice[Y]-1)?0:1;
            flags[2*Z+0] = (z ==                     0)?0:1;
            flags[2*Z+1] = (z == l->block_lattice[Z]-1)?0:1;
            flags[2*T+0] = (t ==                     0)?0:1;
            flags[2*T+1] = (t == l->block_lattice[T]-1)?0:1;
            flags += 8;
          }
        }
      }
    }
  } else {
    for(int t=0; t < l->block_lattice[T]; t++) {
      for(int z=0; z < l->block_lattice[Z]; z++) {
        for(int y=0; y < l->block_lattice[Y]; y++) {
          for(int x=0; x < l->block_lattice[X]; x++) {
            flags[2*X+0] = (x ==                     0)?0:1;
            flags[2*X+1] = (x == l->block_lattice[X]-1)?0:1;
            flags[2*Y+0] = (y ==                     0)?0:1;
            flags[2*Y+1] = (y == l->block_lattice[Y]-1)?0:1;
            flags[2*Z+0] = (z ==                     0)?0:1;
            flags[2*Z+1] = (z == l->block_lattice[Z]-1)?0:1;
            flags[2*T+0] = (t ==                     0)?0:1;
            flags[2*T+1] = (t == l->block_lattice[T]-1)?0:1;
            flags += 8;
          }
        }
      }
    }
  }

  complex_PRECISION eta1[fine_components*OPERATOR_COMPONENT_OFFSET_PRECISION] __attribute__((aligned(64)));
  complex_PRECISION eta2[fine_components*OPERATOR_COMPONENT_OFFSET_PRECISION] __attribute__((aligned(64)));
  complex_PRECISION tmp[4*4*n*OPERATOR_COMPONENT_OFFSET_PRECISION] __attribute__((aligned(64)));

  for ( int a=threading->n_core*threading->thread+threading->core; a<num_aggregates; a+=threading->n_core*threading->n_thread ) {

    // new aggregate is starting, zero out tmp
    for(int i=0; i<4*4*n*OPERATOR_COMPONENT_OFFSET_PRECISION; i++)
      tmp[i] = 0.0;

    for ( int site=a*aggregate_sites; site<(a+1)*aggregate_sites; site++ ) {
      if(l->depth == 0) {
        for ( int c=0; c<l->num_eig_vect; c+=SIMD_LENGTH_PRECISION )
          d_plus_clover_aggregate_PRECISION_vectorized( eta1+c*fine_components, eta2+c*fine_components,
              operator+c*l->vector_size, &(l->s_PRECISION), l, site,
              direction_flags+8*(site%(l->block_lattice[T]*l->block_lattice[Z]*l->block_lattice[Y]*l->block_lattice[X])) );
      } else {
        for ( int c=0; c<l->num_eig_vect; c+=SIMD_LENGTH_PRECISION )
          coarse_aggregate_self_couplings_PRECISION_vectorized( eta1+c*fine_components, eta2+c*fine_components,
              operator+c*l->vector_size, &(l->s_PRECISION), l, site,
              direction_flags+8*(site%(l->block_lattice[T]*l->block_lattice[Z]*l->block_lattice[Y]*l->block_lattice[X])) );
      }
      set_coarse_self_coupling_PRECISION_vectorized( eta1, eta2, operator, l, site, n, tmp );
    }

    // aggregate is done, finalize
    set_coarse_self_coupling_PRECISION_vectorized_finalize( l, a*aggregate_sites, n, tmp );

  }


  SYNC_HYPERTHREADS(threading)
  START_LOCKED_MASTER(threading)
  // neighbors
  for ( int c=0; c<l->num_eig_vect; c+=SIMD_LENGTH_PRECISION ) {
    for ( mu=0; mu<4; mu++ ) {
      // determine start of buffer for this mu
      int start = 0;
      for ( int j=0; j<mu; j++ )
        start += l->s_PRECISION.op.c.num_boundary_sites[2*j];

      // update ghost cells of V[i]
      negative_sendrecv_PRECISION_vectorized( operator+c*l->vector_size, mu, &(l->s_PRECISION.op.c), l,
          SIMD_LENGTH_PRECISION, mpi_buffer+c*(l->vector_size-l->inner_vector_size)+fine_components*start*SIMD_LENGTH_PRECISION );
    }
    for ( mu=0; mu<4; mu++ ) {
      // finish updating ghostcells of V[i]
      negative_wait_PRECISION( mu, &(l->s_PRECISION.op.c), l );
    }
  }
  END_LOCKED_MASTER(threading)
  SYNC_HYPERTHREADS(threading)


  for ( int a=threading->n_core*threading->thread+threading->core; a<num_aggregates; a+=threading->n_core*threading->n_thread ) {

    // new aggregate is starting, zero out tmp
    for(int i=0; i<4*4*n*OPERATOR_COMPONENT_OFFSET_PRECISION; i++)
      tmp[i] = 0.0;

    for ( int site=a*aggregate_sites; site<(a+1)*aggregate_sites; site++ ) {
      for ( mu=0; mu<4; mu++ ) {
        if( (direction_flags+8*(site%(l->block_lattice[T]*l->block_lattice[Z]*l->block_lattice[Y]*l->block_lattice[X])))[2*mu+1] != 0)
          continue;

        if(l->depth == 0)
          for ( int c=0; c<l->num_eig_vect; c+=SIMD_LENGTH_PRECISION )
            d_neighbor_aggregate_PRECISION_vectorized( eta1+c*fine_components, eta2+c*fine_components,
                operator+c*l->vector_size, mu, &(l->s_PRECISION), l, site );
        else
          for ( int c=0; c<l->num_eig_vect; c+=SIMD_LENGTH_PRECISION )
            coarse_aggregate_neighbor_couplings_PRECISION_vectorized( eta1+c*fine_components, eta2+c*fine_components,
                operator+c*l->vector_size, mu, &(l->s_PRECISION), l, site );
        set_coarse_neighbor_coupling_PRECISION_vectorized( eta1, eta2, operator, mu, l, site, n, tmp+mu*4*n*OPERATOR_COMPONENT_OFFSET_PRECISION );
      }
    }

    // aggregate is done, finalize
    for ( mu=0; mu<4; mu++ )
      set_coarse_neighbor_coupling_PRECISION_vectorized_finalize( mu, l, a*aggregate_sites, n, tmp+mu*4*n*OPERATOR_COMPONENT_OFFSET_PRECISION );

    // new aggregate is starting, zero out tmp
    for(int i=0; i<4*4*n*OPERATOR_COMPONENT_OFFSET_PRECISION; i++)
      tmp[i] = 0.0;

    for ( int site=a*aggregate_sites; site<(a+1)*aggregate_sites; site++ ) {
      if(l->depth == 0) {
        for ( int c=0; c<l->num_eig_vect; c+=SIMD_LENGTH_PRECISION )
          diagonal_aggregate_PRECISION_vectorized( eta1+c*fine_components, eta2+c*fine_components,
                                                   operator+c*l->vector_size, &(l->s_PRECISION), l, site );
      } else {
        for ( int c=0; c<l->num_eig_vect; c+=SIMD_LENGTH_PRECISION )
          coarse_aggregate_block_diagonal_PRECISION_vectorized( eta1+c*fine_components, eta2+c*fine_components,
              operator+c*l->vector_size, &(l->s_PRECISION), l, site );
      }
      set_coarse_block_diagonal_PRECISION_vectorized( eta1, eta2, operator, l, site, n, tmp );
    }

    // aggregate is done, finalize
    set_coarse_block_diagonal_PRECISION_vectorized_finalize( l, a*aggregate_sites, n, tmp );

  }

  SYNC_HYPERTHREADS(threading)
  SYNC_CORES(threading)

  coarse_operator_PRECISION_setup_finalize( l, threading );

  START_MASTER(threading)
  FREE_HUGEPAGES( mpi_buffer, complex_PRECISION, OPERATOR_COMPONENT_OFFSET_PRECISION*(l->vector_size-l->inner_vector_size) );

  t1 = MPI_Wtime();
  if ( g.print > 0 ) printf0("depth: %d, time spent for setting up next coarser operator: %lf seconds\n", l->depth, t1-t0 );
  END_MASTER(threading)

  SYNC_HYPERTHREADS(threading)
  SYNC_CORES(threading)
}
#endif
 
void set_coarse_self_coupling_PRECISION_vectorized( complex_PRECISION *spin_0_1, complex_PRECISION *spin_2_3,
    complex_PRECISION *V, level_struct *l, int site, const int n_rhs, complex_PRECISION *tmp ) {

  sse_set_coarse_self_coupling_PRECISION( spin_0_1, spin_2_3, V, l, site, n_rhs, tmp );
}

void set_coarse_block_diagonal_PRECISION_vectorized( complex_PRECISION *spin_0_1, complex_PRECISION *spin_2_3,
    complex_PRECISION *V, level_struct *l, int site, const int n_rhs, complex_PRECISION *tmp ) {

  sse_set_coarse_block_diagonal_PRECISION( spin_0_1, spin_2_3, V, l, site, n_rhs, tmp );
}

void set_coarse_self_coupling_PRECISION_vectorized_finalize( level_struct *l, int site, const int n_rhs, complex_PRECISION *tmp ) {

  int k, k1, k2, num_aggregates = l->is_PRECISION.num_agg,
      num_eig_vect = l->next_level->num_lattice_site_var/2,
      aggregate_size = l->inner_vector_size / num_aggregates,
      clover_site_size = (l->next_level->num_lattice_site_var*(l->next_level->num_lattice_site_var+1))/2;
  int t1, t2;

  config_PRECISION clover_pt, clover = l->next_level->op_PRECISION.clover;

  // just an abbreviation
  int component_offset = OPERATOR_COMPONENT_OFFSET_PRECISION;
  int fine_components = l->num_lattice_site_var;

  int aggregate = (fine_components*site)/aggregate_size;
  clover_pt = clover + aggregate*clover_site_size;

  // U(x) = [ A B      , A=A*, D=D*, C = -B*
  //          C D ]
  // storage order: upper triangle of A, upper triangle of D, B, columnwise
  // diagonal coupling
  for ( int n=0; n<n_rhs; n++ ) {

    // index k used for vectorization
    for ( k=0; k<=n; k++ ) {

      k1 = (n*(n+1))/2;
      k2 = (n*(n+1))/2+(num_eig_vect*(num_eig_vect+1))/2;
      t1 = (n+0*num_eig_vect)*component_offset;
      t2 = (n+1*num_eig_vect)*component_offset;

      // A
      clover_pt[ k1+k ] += ((float *)(tmp+t1))[k] + I * ((float *)(tmp+t1)+component_offset)[k];

      // D
      clover_pt[ k2+k ] += ((float *)(tmp+t2))[k] + I * ((float *)(tmp+t2)+component_offset)[k];
    }

    // index k used for vectorization
    for ( k=0; k<num_eig_vect; k++ ) {

      k1 = num_eig_vect*(num_eig_vect+1+n);
      t1 = (n+2*num_eig_vect)*component_offset;

      // B
      clover_pt[ k1+k ] += ((float *)(tmp+t1))[k] + I * ((float *)(tmp+t1)+component_offset)[k];
    }
  }
}


void set_coarse_neighbor_coupling_PRECISION_vectorized( complex_PRECISION *spin_0_1, complex_PRECISION *spin_2_3, 
                                                        complex_PRECISION *V, const int mu, level_struct *l, int site,
                                                        const int n_rhs, complex_PRECISION *tmp ) {

  sse_set_coarse_neighbor_coupling_PRECISION( spin_0_1, spin_2_3, V, mu, l, site, n_rhs, tmp );
}


void set_coarse_neighbor_coupling_PRECISION_vectorized_finalize( const int mu, level_struct *l, int site,
                                                                 const int n_rhs, complex_PRECISION *tmp ) {

  int k, k1, k2, num_eig_vect = l->next_level->num_lattice_site_var/2,
      D_link_size = num_eig_vect*num_eig_vect*4;
  int t1, t2;

  config_PRECISION D_pt, D = l->next_level->op_PRECISION.D;

  // just an abbreviation
  int component_offset = OPERATOR_COMPONENT_OFFSET_PRECISION;
  int fine_components = l->num_lattice_site_var;

  int aggregate = (fine_components*site)/(l->inner_vector_size / l->is_PRECISION.num_agg);
  D_pt = D + (4*aggregate+mu)*D_link_size;

  // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
  //             C D ]                        -B*  D* ]
  // storage order: A, C, B, D, each column wise
  for ( int n=0; n<n_rhs; n++ ) {

    // index k used for vectorization
    for ( k=0; k<num_eig_vect; k++ ) {

      k1 = (n+0*num_eig_vect)*num_eig_vect;
      k2 = (n+1*num_eig_vect)*num_eig_vect;
      t1 = (n+0*num_eig_vect)*component_offset;
      t2 = (n+1*num_eig_vect)*component_offset;

      
      // A
      D_pt[ k1+k ] += ((float *)(tmp+t1))[k] + I * ((float *)(tmp+t1)+component_offset)[k];

      // C
      D_pt[ k2+k ] += ((float *)(tmp+t2))[k] + I * ((float *)(tmp+t2)+component_offset)[k];


      k1 = (n+2*num_eig_vect)*num_eig_vect;
      k2 = (n+3*num_eig_vect)*num_eig_vect;
      t1 = (n+2*num_eig_vect)*component_offset;
      t2 = (n+3*num_eig_vect)*component_offset;

      // B
      D_pt[ k1+k ] += ((float *)(tmp+t1))[k] + I * ((float *)(tmp+t1)+component_offset)[k];

      // D
      D_pt[ k2+k ] += ((float *)(tmp+t2))[k] + I * ((float *)(tmp+t2)+component_offset)[k];
    }
  }
}

void set_coarse_block_diagonal_PRECISION_vectorized_finalize( level_struct *l, int site, const int n_rhs, complex_PRECISION *tmp ) {

  int k, k1, k2, num_aggregates = l->is_PRECISION.num_agg,
      num_eig_vect = l->next_level->num_parent_eig_vect,
      aggregate_size = l->inner_vector_size / num_aggregates,
      block_site_size = (l->next_level->num_parent_eig_vect*(l->next_level->num_parent_eig_vect+1));
  int t1, t2;

  config_PRECISION block_pt, block = l->next_level->op_PRECISION.odd_proj;

  // just an abbreviation
  int component_offset = OPERATOR_COMPONENT_OFFSET_PRECISION;
  int fine_components = l->num_lattice_site_var;

  int aggregate = (fine_components*site)/aggregate_size;
  block_pt = block + aggregate*block_site_size;

  // U(x) = [ A 0      , A=A*, D=D*
  //          0 D ]
  // storage order: upper triangle of A, upper triangle of D, B, columnwise
  // diagonal coupling
  for ( int n=0; n<n_rhs; n++ ) {

    // index k used for vectorization
    for ( k=0; k<=n; k++ ) {

      k1 = (n*(n+1))/2;
      k2 = (n*(n+1))/2+(num_eig_vect*(num_eig_vect+1))/2;
      t1 = (n+0*num_eig_vect)*component_offset;
      t2 = (n+1*num_eig_vect)*component_offset;

      // A
      block_pt[ k1+k ] += ((float *)(tmp+t1))[k] + I * ((float *)(tmp+t1)+component_offset)[k];

      // D
      block_pt[ k2+k ] += ((float *)(tmp+t2))[k] + I * ((float *)(tmp+t2)+component_offset)[k];
    }
  }
}

void copy_coarse_operator_to_vectorized_layout_PRECISION( config_PRECISION D,
                                                          OPERATOR_TYPE_PRECISION *D_vectorized, 
                                                          int num_aggregates, int num_eig_vect) {

  int vecs = num_eig_vect;
  // in vectorized layout D is stored column wise, but not split into ABCD
  // each column is padded, such that next column can also start at 64B boundary
  int half_offset = SIMD_LENGTH_PRECISION*((vecs+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  int column_offset = 2*half_offset;
  // offset between blocks in D
  int block_offset = vecs*vecs;

  PRECISION *out_tmp = D_vectorized;

  // we zero out the padded area (o) to avoid potential floating-point errors
  // D_vectorized is
  // AB
  // oo
  // CD
  // oo
  // (column wise, size of zeros such that columns length is multiple of 64B)

  // 4 directions
  for ( int a=0; a<4*num_aggregates; a++ ) {
    for(int i=0; i<vecs; i++) {
      for(int j=0; j<vecs; j++) {
        // A
        out_tmp[(2*i+0)*column_offset + j] = creal(D[0*block_offset + i*vecs+j]);
        out_tmp[(2*i+1)*column_offset + j] = cimag(D[0*block_offset + i*vecs+j]);
        // C
        out_tmp[(2*i+0)*column_offset + j + half_offset] = creal(D[1*block_offset + i*vecs+j]);
        out_tmp[(2*i+1)*column_offset + j + half_offset] = cimag(D[1*block_offset + i*vecs+j]);
      }
      // zero
      for(int j=vecs; j<half_offset; j++) {
        out_tmp[(2*i+0)*column_offset + j] = 0.0;
        out_tmp[(2*i+1)*column_offset + j] = 0.0;
        out_tmp[(2*i+0)*column_offset + j + half_offset] = 0.0;
        out_tmp[(2*i+1)*column_offset + j + half_offset] = 0.0;
      }
    }

    for(int i=0; i<vecs; i++) {
      for(int j=0; j<vecs; j++) {
        // B
        out_tmp[(2*(i+vecs)+0)*column_offset + j] = creal(D[2*block_offset + i*vecs+j]);
        out_tmp[(2*(i+vecs)+1)*column_offset + j] = cimag(D[2*block_offset + i*vecs+j]);
        // D
        out_tmp[(2*(i+vecs)+0)*column_offset + j + half_offset] = creal(D[3*block_offset + i*vecs+j]);
        out_tmp[(2*(i+vecs)+1)*column_offset + j + half_offset] = cimag(D[3*block_offset + i*vecs+j]);
      }
      // zero
      for(int j=vecs; j<half_offset; j++) {
        out_tmp[(2*(i+vecs)+0)*column_offset + j] = 0.0;
        out_tmp[(2*(i+vecs)+1)*column_offset + j] = 0.0;
        out_tmp[(2*(i+vecs)+0)*column_offset + j + half_offset] = 0.0;
        out_tmp[(2*(i+vecs)+1)*column_offset + j + half_offset] = 0.0;
      }
    }
    D += 2*vecs*2*vecs;
    // out_tmp is an alias for the actual output
    out_tmp += 2*column_offset*2*vecs;
  }
}


void copy_coarse_operator_to_transformed_vectorized_layout_PRECISION(config_PRECISION D,
                                                                     OPERATOR_TYPE_PRECISION *D_vectorized,
                                                                     int num_aggregates, int num_eig_vect) {

  int vecs = num_eig_vect;
  // in vectorized layout D is stored column wise, but not split into ABCD
  // output is transposed
  // each column is padded, such that the next column can also start at 64bit boundary
  int half_offset = SIMD_LENGTH_PRECISION*((vecs+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  int column_offset = 2*half_offset;
  // offset between blocks in D
  int block_offset = vecs*vecs;

  PRECISION *out_tmp = D_vectorized;

  // we zero out the padded area to avoid potential floating-point errors
  // D_vectorized is
  // A^T C^T
  //  o   o
  // B^T D^T
  //  o   o
  // (column wise, size of zeros such that columns length is multiple of 64B)

  // 4 directions
  for ( int a=0; a<4*num_aggregates; a++ ) {
    for(int i=0; i<vecs; i++) {
      for(int j=0; j<vecs; j++) {
        // A
        out_tmp[(2*i+0)*column_offset + j] = creal(D[0*block_offset + j*vecs+i]);
        out_tmp[(2*i+1)*column_offset + j] = -cimag(D[0*block_offset + j*vecs+i]);
        // B
        out_tmp[(2*i+0)*column_offset + j + half_offset] = -creal(D[2*block_offset + j*vecs+i]);
        out_tmp[(2*i+1)*column_offset + j + half_offset] = cimag(D[2*block_offset + j*vecs+i]);
      }
      // zero
      for(int j=vecs; j<half_offset; j++) {
        out_tmp[(2*i+0)*column_offset + j] = 0.0;
        out_tmp[(2*i+1)*column_offset + j] = 0.0;
        out_tmp[(2*i+0)*column_offset + j + half_offset] = 0.0;
        out_tmp[(2*i+1)*column_offset + j + half_offset] = 0.0;
      }
    }

    for(int i=0; i<vecs; i++) {
      for(int j=0; j<vecs; j++) {
        // C
        out_tmp[(2*(i+vecs)+0)*column_offset + j] = -creal(D[1*block_offset + j*vecs+i]);
        out_tmp[(2*(i+vecs)+1)*column_offset + j] = cimag(D[1*block_offset + j*vecs+i]);
        // D
        out_tmp[(2*(i+vecs)+0)*column_offset + j + half_offset] = creal(D[3*block_offset + j*vecs+i]);
        out_tmp[(2*(i+vecs)+1)*column_offset + j + half_offset] = -cimag(D[3*block_offset + j*vecs+i]);
      }
      // zero
      for(int j=vecs; j<half_offset; j++) {
        out_tmp[(2*(i+vecs)+0)*column_offset + j] = 0.0;
        out_tmp[(2*(i+vecs)+1)*column_offset + j] = 0.0;
        out_tmp[(2*(i+vecs)+0)*column_offset + j + half_offset] = 0.0;
        out_tmp[(2*(i+vecs)+1)*column_offset + j + half_offset] = 0.0;
      }
    }
    D += 2*vecs*2*vecs;
    // out_tmp is an alias for the actual output
    out_tmp += 2*column_offset*2*vecs;
  }
}


void copy_coarse_operator_clover_to_vectorized_layout_PRECISION(config_PRECISION clover,
                                                                OPERATOR_TYPE_PRECISION *clover_vectorized,
                                                                int num_aggregates, int num_eig_vect) {

  int vecs = num_eig_vect;
  // in vectorized layout clover is stored column wise, but not split into ABCD
  // each column is padded, such that next column can also start at 64B boundary
  int column_offset = SIMD_LENGTH_PRECISION*((2*num_eig_vect+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  // offset between blocks in clover
  int offset_to_D = (vecs*vecs+vecs)/2; // upper triangle of A including diagonal
  int offset_to_B = 2*offset_to_D; // B comes after A and D

  PRECISION *out_tmp = clover_vectorized;

  // we zero out the padded area to avoid potential floating-point errors
  // cloverD_vectorized is
  // AB
  // CD
  // 00
  // (column wise, size of zeros such that columns length is multiple of 64B)

  // 4 directions
  for ( int a=0; a<num_aggregates; a++ ) {
    for(int i=0; i<vecs; i++) {
      for(int j=0; j<vecs; j++) {
        // A
        // primed indices to transpose input when necessary to get lower triangle of output
        int ip = i;
        int jp = j;
        PRECISION sign = 1.0;
        if(j > i) {
          ip = j;
          jp = i;
          sign = -1.0;
        }
        int offset_to_column = (ip*ip+ip)/2; // upper triangle including diagonal
        out_tmp[(2*i+0)*column_offset + j] = creal(clover[offset_to_column+jp]);
        out_tmp[(2*i+1)*column_offset + j] = sign*cimag(clover[offset_to_column+jp]);
        // C = -B^dagger
        out_tmp[(2*i+0)*column_offset + j + vecs] = -creal(clover[offset_to_B + j*vecs+i]);
        out_tmp[(2*i+1)*column_offset + j + vecs] =  cimag(clover[offset_to_B + j*vecs+i]);
      }
      // zero
      for(int j=2*vecs; j<column_offset; j++) {
        out_tmp[(2*i+0)*column_offset + j] = 0.0;
        out_tmp[(2*i+1)*column_offset + j] = 0.0;
      }
    }

    for(int i=0; i<vecs; i++) {
      for(int j=0; j<vecs; j++) {
        // B
        out_tmp[(2*(i+vecs)+0)*column_offset + j] = creal(clover[offset_to_B + i*vecs+j]);
        out_tmp[(2*(i+vecs)+1)*column_offset + j] = cimag(clover[offset_to_B + i*vecs+j]);
        // D
        // primed indices to transpose input when necessary to get lower triangle of output
        int ip = i;
        int jp = j;
        PRECISION sign = 1.0;
        if(j > i) {
          ip = j;
          jp = i;
          sign = -1.0;
        }
        int offset_to_column = (ip*ip+ip)/2; // upper triangle including diagonal
        out_tmp[(2*(i+vecs)+0)*column_offset + j + vecs] = creal(clover[offset_to_D + offset_to_column+jp]);
        out_tmp[(2*(i+vecs)+1)*column_offset + j + vecs] = sign*cimag(clover[offset_to_D + offset_to_column+jp]);
      }
      // zero
      for(int j=2*vecs; j<column_offset; j++) {
        out_tmp[(2*(i+vecs)+0)*column_offset + j] = 0.0;
        out_tmp[(2*(i+vecs)+1)*column_offset + j] = 0.0;
      }
    }
    clover += offset_to_B + vecs*vecs;
    // out_tmp is an alias for the actual output
    out_tmp += 2*column_offset*2*vecs;
  }
}

void copy_coarse_operator_clover_to_doublet_vectorized_layout_PRECISION(config_PRECISION clover,
                                                                        OPERATOR_TYPE_PRECISION *clover_vectorized,
                                                                        int num_aggregates, int num_eig_vect) {

  int vecs = num_eig_vect;
  // in vectorized layout clover is stored column wise, but not split into ABCD
  // each column is padded, such that next column can also start at 64B boundary
  int column_offset = SIMD_LENGTH_PRECISION*((4*vecs+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  // offset between blocks in clover
  int offset_to_D = (vecs*vecs+vecs)/2; // upper triangle of A including diagonal
  int offset_to_B = 2*offset_to_D; // B comes after A and D

  PRECISION *out_tmp = clover_vectorized;

  // cloverD_vectorized is
  // A0B0
  // 0A0B
  // C0D0
  // 0C0D
  // 0000  we zero out the padded area to avoid potential floating-point errors
  // (column wise, size of zeros such that columns length is multiple of 64B)

  // 4 directions
  for ( int a=0; a<num_aggregates; a++ ) {
    for(int i=0; i<vecs; i++) {
      for(int j=0; j<vecs; j++) {
        // primed indices to transpose input when necessary to get lower triangle of output
        int ip = i;
        int jp = j;
        PRECISION sign = 1.0;
        if(j > i) {
          ip = j;
          jp = i;
          sign = -1.0;
        }
        int offset_to_column = (ip*ip+ip)/2; // upper triangle including diagonal
        // A
        out_tmp[(2*(i+0*vecs)+0)*column_offset + j + 0*vecs] =
        out_tmp[(2*(i+1*vecs)+0)*column_offset + j + 1*vecs] = creal(clover[offset_to_column+jp]);
        out_tmp[(2*(i+0*vecs)+1)*column_offset + j + 0*vecs] =
        out_tmp[(2*(i+1*vecs)+1)*column_offset + j + 1*vecs] = sign*cimag(clover[offset_to_column+jp]);
        // B
        out_tmp[(2*(i+2*vecs)+0)*column_offset + j + 0*vecs] =
        out_tmp[(2*(i+3*vecs)+0)*column_offset + j + 1*vecs] = creal(clover[offset_to_B + i*vecs+j]);
        out_tmp[(2*(i+2*vecs)+1)*column_offset + j + 0*vecs] = 
        out_tmp[(2*(i+3*vecs)+1)*column_offset + j + 1*vecs] = cimag(clover[offset_to_B + i*vecs+j]);
        // C = -B^dagger
        out_tmp[(2*(i+0*vecs)+0)*column_offset + j + 2*vecs] =
        out_tmp[(2*(i+1*vecs)+0)*column_offset + j + 3*vecs] = -creal(clover[offset_to_B + j*vecs+i]);
        out_tmp[(2*(i+0*vecs)+1)*column_offset + j + 2*vecs] =
        out_tmp[(2*(i+1*vecs)+1)*column_offset + j + 3*vecs] = cimag(clover[offset_to_B + j*vecs+i]);
        // D
        out_tmp[(2*(i+2*vecs)+0)*column_offset + j + 2*vecs] = 
        out_tmp[(2*(i+3*vecs)+0)*column_offset + j + 3*vecs] = creal(clover[offset_to_D + offset_to_column+jp]);
        out_tmp[(2*(i+2*vecs)+1)*column_offset + j + 2*vecs] = 
        out_tmp[(2*(i+3*vecs)+1)*column_offset + j + 3*vecs] = sign*cimag(clover[offset_to_D + offset_to_column+jp]);
        // 0
        out_tmp[(2*(i+0*vecs)+0)*column_offset + j + 1*vecs] =
        out_tmp[(2*(i+0*vecs)+1)*column_offset + j + 1*vecs] = 
        out_tmp[(2*(i+0*vecs)+0)*column_offset + j + 3*vecs] =
        out_tmp[(2*(i+0*vecs)+1)*column_offset + j + 3*vecs] = 
        out_tmp[(2*(i+1*vecs)+0)*column_offset + j + 0*vecs] =
        out_tmp[(2*(i+1*vecs)+1)*column_offset + j + 0*vecs] =
        out_tmp[(2*(i+1*vecs)+0)*column_offset + j + 2*vecs] =
        out_tmp[(2*(i+1*vecs)+1)*column_offset + j + 2*vecs] =
        out_tmp[(2*(i+2*vecs)+0)*column_offset + j + 1*vecs] =
        out_tmp[(2*(i+2*vecs)+1)*column_offset + j + 1*vecs] =
        out_tmp[(2*(i+2*vecs)+0)*column_offset + j + 3*vecs] =
        out_tmp[(2*(i+2*vecs)+1)*column_offset + j + 3*vecs] =
        out_tmp[(2*(i+3*vecs)+0)*column_offset + j + 0*vecs] =
        out_tmp[(2*(i+3*vecs)+1)*column_offset + j + 0*vecs] =
        out_tmp[(2*(i+3*vecs)+0)*column_offset + j + 2*vecs] =
        out_tmp[(2*(i+3*vecs)+1)*column_offset + j + 2*vecs] = 0.0;
      }
      // zero
      for(int j=4*vecs; j<column_offset; j++) {
        out_tmp[(2*(i+0*vecs)+0)*column_offset + j] = 
        out_tmp[(2*(i+0*vecs)+1)*column_offset + j] = 
        out_tmp[(2*(i+1*vecs)+0)*column_offset + j] = 
        out_tmp[(2*(i+1*vecs)+1)*column_offset + j] = 
        out_tmp[(2*(i+2*vecs)+0)*column_offset + j] = 
        out_tmp[(2*(i+2*vecs)+1)*column_offset + j] = 
        out_tmp[(2*(i+3*vecs)+0)*column_offset + j] = 
        out_tmp[(2*(i+3*vecs)+1)*column_offset + j] = 0.0;
      }
    }
    clover += offset_to_B + vecs*vecs;
    // out_tmp is an alias for the actual output
    out_tmp += 2*4*vecs*column_offset;
  }
}


void add_tm_term_to_vectorized_layout_PRECISION(config_PRECISION tm_term, OPERATOR_TYPE_PRECISION *clover_vectorized,
                                                int num_aggregates, int num_eig_vect) {
#ifdef HAVE_TM
  int vecs = num_eig_vect;
  // in vectorized layout clover is stored column wise, but not split into ABCD
  // each column is padded, such that next column can also start at 64B boundary
  int column_offset = SIMD_LENGTH_PRECISION*((2*num_eig_vect+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  // offset between blocks in clover
  int offset_to_F = (vecs*vecs+vecs)/2; // upper triangle of A including diagonal

  PRECISION *out_tmp = clover_vectorized;

  // we add the tm term to cloverD_vectorized
  // AB   E0
  // CD + 0F
  // 00   00
  // (column wise, size of zeros such that columns length is multiple of 64B)

  // 4 directions
  for ( int a=0; a<num_aggregates; a++ ) {
    for(int i=0; i<vecs; i++) {
      for(int j=0; j<vecs; j++) {
        // primed indices to transpose input when necessary to get lower triangle of output
        int ip = i;
        int jp = j;
        PRECISION sign = 1.0;
        if(j > i) {
          ip = j;
          jp = i;
          sign = -1.0;
        }
        int offset_to_column = (ip*ip+ip)/2; // upper triangle including diagonal
        // E
        out_tmp[(2*i+0)*column_offset + j] += sign*creal(tm_term[offset_to_column+jp]);
        out_tmp[(2*i+1)*column_offset + j] += cimag(tm_term[offset_to_column+jp]);
        // F
        out_tmp[(2*(i+vecs)+0)*column_offset + j + vecs] += sign*creal(tm_term[offset_to_F + offset_to_column+jp]);
        out_tmp[(2*(i+vecs)+1)*column_offset + j + vecs] += cimag(tm_term[offset_to_F + offset_to_column+jp]);
      }
    }
    tm_term += 2*offset_to_F;
    // out_tmp is an alias for the actual output
    out_tmp += 2*column_offset*2*vecs;
  }
#endif
}

void add_tm_term_to_doublet_vectorized_layout_PRECISION(config_PRECISION tm_term, OPERATOR_TYPE_PRECISION *clover_vectorized,
                                                int num_aggregates, int num_eig_vect) {
#ifdef HAVE_TM
  int vecs = num_eig_vect;
  // in vectorized layout clover is stored column wise, but not split into ABCD
  // each column is padded, such that next column can also start at 64B boundary
  int column_offset = SIMD_LENGTH_PRECISION*((4*vecs+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  // offset between blocks in clover
  int offset_to_F = (vecs*vecs+vecs)/2; // upper triangle of A including diagonal

  PRECISION *out_tmp = clover_vectorized;

  // we add/sub the tm term to cloverD_vectorized
  // A0B0   E000   0000
  // 0A0B + 0000 - 0E00
  // C0D0   00F0   0000
  // 0C0D   0000   000F
  // 0000   0000   0000
  // (column wise, size of zeros such that columns length is multiple of 64B)

  // 4 directions
  for ( int a=0; a<num_aggregates; a++ ) {
    for(int i=0; i<vecs; i++) {
      for(int j=0; j<vecs; j++) {
        // primed indices to transpose input when necessary to get lower triangle of output
        int ip = i;
        int jp = j;
        PRECISION sign = 1.0;
        if(j > i) {
          ip = j;
          jp = i;
          sign = -1.0;
        }
        int offset_to_column = (ip*ip+ip)/2; // upper triangle including diagonal
        // E
        out_tmp[(2*(i+0*vecs)+0)*column_offset + j + 0*vecs] += sign*creal(tm_term[offset_to_column+jp]);
        out_tmp[(2*(i+0*vecs)+1)*column_offset + j + 0*vecs] += cimag(tm_term[offset_to_column+jp]); 
        out_tmp[(2*(i+1*vecs)+0)*column_offset + j + 1*vecs] -= sign*creal(tm_term[offset_to_column+jp]);
        out_tmp[(2*(i+1*vecs)+1)*column_offset + j + 1*vecs] -= cimag(tm_term[offset_to_column+jp]);
        // F
        out_tmp[(2*(i+2*vecs)+0)*column_offset + j + 2*vecs] += sign*creal(tm_term[offset_to_F+offset_to_column+jp]);
        out_tmp[(2*(i+2*vecs)+1)*column_offset + j + 2*vecs] += cimag(tm_term[offset_to_F+offset_to_column+jp]);
        out_tmp[(2*(i+3*vecs)+0)*column_offset + j + 3*vecs] -= sign*creal(tm_term[offset_to_F+offset_to_column+jp]);
        out_tmp[(2*(i+3*vecs)+1)*column_offset + j + 3*vecs] -= cimag(tm_term[offset_to_F+offset_to_column+jp]);
     }
    }
    tm_term += 2*offset_to_F;
    // out_tmp is an alias for the actual output
    out_tmp += 2*4*vecs*column_offset;
  }
#endif
}

void add_epsbar_term_to_doublet_vectorized_layout_PRECISION(config_PRECISION eps_term, OPERATOR_TYPE_PRECISION *clover_vectorized,
                                                            int num_aggregates, int num_eig_vect) {
#ifdef HAVE_TM1p1
  int vecs = num_eig_vect;
  // in vectorized layout clover is stored column wise, but not split into ABCD
  // each column is padded, such that next column can also start at 64B boundary
  int column_offset = SIMD_LENGTH_PRECISION*((4*vecs+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  // offset between blocks in clover
  int offset_to_F = (vecs*vecs+vecs)/2; // upper triangle of A including diagonal

  PRECISION *out_tmp = clover_vectorized;

  // we add the eps term to cloverD_vectorized
  // A0B0   0E00
  // 0A0B + E000
  // C0D0   000F
  // 0C0D   00F0
  // 0000   0000
  // (column wise, size of zeros such that columns length is multiple of 64B)

  // 4 directions
  for ( int a=0; a<num_aggregates; a++ ) {
    for(int i=0; i<vecs; i++) {
      for(int j=0; j<vecs; j++) {
        // primed indices to transpose input when necessary to get lower triangle of output
        int ip = i;
        int jp = j;
        PRECISION sign = 1.0;
        if(j > i) {
          ip = j;
          jp = i;
          sign = -1.0;
        }
        int offset_to_column = (ip*ip+ip)/2; // upper triangle including diagonal
        // E
        out_tmp[(2*(i+0*vecs)+0)*column_offset + j + 1*vecs] += sign*creal(eps_term[offset_to_column+jp]);
        out_tmp[(2*(i+0*vecs)+1)*column_offset + j + 1*vecs] += cimag(eps_term[offset_to_column+jp]); 
        out_tmp[(2*(i+1*vecs)+0)*column_offset + j + 0*vecs] += sign*creal(eps_term[offset_to_column+jp]);
        out_tmp[(2*(i+1*vecs)+1)*column_offset + j + 0*vecs] += cimag(eps_term[offset_to_column+jp]);
        // F
        out_tmp[(2*(i+2*vecs)+0)*column_offset + j + 3*vecs] += sign*creal(eps_term[offset_to_F+offset_to_column+jp]);
        out_tmp[(2*(i+2*vecs)+1)*column_offset + j + 3*vecs] += cimag(eps_term[offset_to_F+offset_to_column+jp]);
        out_tmp[(2*(i+3*vecs)+0)*column_offset + j + 2*vecs] += sign*creal(eps_term[offset_to_F+offset_to_column+jp]);
        out_tmp[(2*(i+3*vecs)+1)*column_offset + j + 2*vecs] += cimag(eps_term[offset_to_F+offset_to_column+jp]);
     }
    }
    eps_term += 2*offset_to_F;
    // out_tmp is an alias for the actual output
    out_tmp += 2*4*vecs*column_offset;
  }
#endif
}

void coarse_aggregate_self_couplings_PRECISION_vectorized( complex_PRECISION *eta1, complex_PRECISION *eta2, 
                                                           complex_PRECISION *phi, schwarz_PRECISION_struct *s,
                                                           level_struct *l, int site, int *direction_flags ) {

  int offset = SIMD_LENGTH_PRECISION;
  int site_offset = l->num_lattice_site_var*offset;
  int index_bw;
  int index_fw;
  int *neighbor = s->op.neighbor_table;
  int *backward_neighbor = s->op.backward_neighbor_table;
  complex_PRECISION *phi_pt;
  config_PRECISION D_pt;
  config_PRECISION D = s->op.D;
  int n = l->num_lattice_site_var;
  int D_site_offset = 4*n*n;
  int D_link_offset = n*n;
  int clover_offset = (n*(n+1))/2*site;

  coarse_spinwise_site_self_couplings_PRECISION_vectorized( eta1, eta2, phi+site_offset*site, s->op.clover+clover_offset, offset, l );

  for(int mu=0; mu<4; mu++) {
    index_fw  = neighbor[5*site+1 + mu];
    index_bw  = backward_neighbor[5*site+1 + mu];

    // from backward
    if ( direction_flags[2*mu+0] == 1 ) {
      D_pt = D + D_site_offset*index_bw + D_link_offset*mu;
      phi_pt = phi + site_offset*index_bw;
      coarse_spinwise_n_daggered_hopp_PRECISION_vectorized( eta1, eta2, phi_pt, D_pt, offset, l );
    }

    // from forward
    if ( direction_flags[2*mu+1] == 1 ) {
      D_pt = D + D_site_offset*site + D_link_offset*mu;
      phi_pt = phi + site_offset*index_fw;
      coarse_spinwise_n_hopp_PRECISION_vectorized( eta1, eta2, phi_pt, D_pt, offset, l );
    }
  }
}

void coarse_aggregate_block_diagonal_PRECISION_vectorized( complex_PRECISION *eta1, complex_PRECISION *eta2, 
                                                           complex_PRECISION *phi, schwarz_PRECISION_struct *s,
                                                           level_struct *l, int site ) {

  int offset = SIMD_LENGTH_PRECISION;
  int site_offset = l->num_lattice_site_var*offset;
  int n = l->num_parent_eig_vect;
  int block_offset = (n*(n+1))*site;
 
  sse_coarse_aggregate_block_diagonal_PRECISION( eta1, eta2, phi+site_offset*site, s->op.odd_proj+block_offset, offset, l );
}

void coarse_aggregate_neighbor_couplings_PRECISION_vectorized( complex_PRECISION *eta1, complex_PRECISION *eta2, 
                                                               complex_PRECISION *phi, const int mu,
                                                               schwarz_PRECISION_struct *s, level_struct *l, int site ) {

  int offset = SIMD_LENGTH_PRECISION;
  int site_offset = l->num_lattice_site_var*offset;
  int index_fw;
  int *neighbor = s->op.neighbor_table;
  complex_PRECISION *phi_pt;
  config_PRECISION D_pt;
  config_PRECISION D = s->op.D;
  int n = l->num_lattice_site_var;
  int D_site_offset = 4*n*n;
  int D_link_offset = n*n;

  vector_PRECISION_define_zero( eta1, 0, n*offset, l, no_threading );
  vector_PRECISION_define_zero( eta2, 0, n*offset, l, no_threading );

  // requires the positive boundaries of phi to be communicated before
  index_fw  = neighbor[5*site+1 + mu];
  D_pt = D + D_site_offset*site + D_link_offset*mu;
  phi_pt = phi + site_offset*index_fw;
  coarse_spinwise_hopp_PRECISION_vectorized( eta1, eta2, phi_pt, D_pt, offset, l );
}


void coarse_spinwise_site_self_couplings_PRECISION_vectorized(
    complex_PRECISION *eta1, complex_PRECISION *eta2,
    complex_PRECISION *phi, config_PRECISION clover, int elements, level_struct *l ) {
  
  sse_coarse_spinwise_site_self_couplings_PRECISION( eta1, eta2, phi, clover, elements, l );
}

#endif
