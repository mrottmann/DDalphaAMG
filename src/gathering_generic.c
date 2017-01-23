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

void gathering_PRECISION_next_level_init( gathering_PRECISION_struct *gs, level_struct *l ) {
  
  int mu;
  
  gs->permutation = NULL;
  gs->gather_list = NULL;
  gs->reqs = NULL;
  gs->buffer = NULL;
  gs->transfer_buffer = NULL;
  
  gs->dist_inner_lattice_sites = 1;
  gs->gather_list_length = 1;
  for ( mu=0; mu<4; mu++ ) {
    gs->dist_local_lattice[mu] = l->local_lattice[mu] / l->coarsening[mu];
    gs->dist_inner_lattice_sites *= gs->dist_local_lattice[mu];
    gs->gather_list_length *= l->next_level->local_lattice[mu] / gs->dist_local_lattice[mu];
  }
}


void gathering_PRECISION_setup( gathering_PRECISION_struct *gs, level_struct *l ) {
  
  // define data merging
  // define data gathering permutation
  int i, mu, current_rank, offset, offset_sum,
      process_coords[4] = {0,0,0,0}, parent_coords[4] = {0,0,0,0}, *process_list = NULL;
  MALLOC( process_list, int, l->num_processes );
  MALLOC( gs->transfer_buffer, complex_PRECISION, gs->dist_inner_lattice_sites * l->num_lattice_site_var );
  
  l->idle = 0;
  i = 0;
  for ( process_coords[T]=0; process_coords[T]<g.process_grid[T]; process_coords[T]++ )
    for ( process_coords[Z]=0; process_coords[Z]<g.process_grid[Z]; process_coords[Z]++ )
      for ( process_coords[Y]=0; process_coords[Y]<g.process_grid[Y]; process_coords[Y]++ )
        for ( process_coords[X]=0; process_coords[X]<g.process_grid[X]; process_coords[X]++ ) {
          MPI_Cart_rank( g.comm_cart, process_coords, &current_rank );
          offset_sum = 0;
          for ( mu=0; mu<4; mu++ ) {
            offset = process_coords[mu] % l->comm_offset[mu];
            parent_coords[mu] = process_coords[mu] - offset;
            offset_sum+=offset;
          }
          // get parent rank
          if ( current_rank == g.my_rank ) {
            MPI_Cart_rank( g.comm_cart, parent_coords, &(l->parent_rank) );
            // find out if current process is supposed to idle
            if ( offset_sum > 0 )
              l->idle = 1;
          }
          // store current rank in case the respective process is not supposed to idle
          if ( offset_sum == 0 ) {
            process_list[i] = current_rank;
            i++;
          }
        }
        
  // creating sub communicator for inner products
  MPI_Comm_group( g.comm_cart, &(gs->level_comm_group) );
  MPI_Group_incl( g.global_comm_group, l->num_processes, process_list, &(gs->level_comm_group) );
  MPI_Comm_create( g.comm_cart, gs->level_comm_group, &(gs->level_comm) );
  FREE( process_list, int, l->num_processes );
  
  if ( !l->idle ) {
    int d0, c0, b0, a0, d1, c1, b1, a1, t, z, y, x, k, j, *field1=NULL, *field2=NULL, *count[4],
    merge[4], block_size[4], block_split[4], agg_split[4];
    MALLOC( gs->gather_list, int, gs->gather_list_length );
    MALLOC( gs->permutation, int, l->num_inner_lattice_sites );
    MALLOC( gs->reqs, MPI_Request, gs->gather_list_length );
    MALLOC( gs->buffer, complex_PRECISION, l->inner_vector_size );
    MALLOC( field1, int, l->num_inner_lattice_sites );
    MALLOC( field2, int, l->num_inner_lattice_sites );
    
    for ( mu=0; mu<4; mu++ ) {
      block_size[mu] = gs->dist_local_lattice[mu];
      merge[mu] = l->local_lattice[mu] / block_size[mu];
    }
    
    // process wise linearly ordered, process blocks linearly ordered
    // ---> completely linearly ordered
    count[0]=&d0; count[1]=&c0; count[2]=&b0; count[3]=&a0;
    i=0; j=0;
    for ( d0=0; d0<merge[T]; d0++)
      for ( c0=0; c0<merge[Z]; c0++)
        for ( b0=0; b0<merge[Y]; b0++ )
          for ( a0=0; a0<merge[X]; a0++ ) {
            // determine list off processes for gathering and distribution data
            for ( mu=0; mu<4; mu++ )
              process_coords[mu] = g.my_coords[mu] + *(count[mu]) * (l->comm_offset[mu]/merge[mu]);
            MPI_Cart_rank( g.comm_cart, process_coords, gs->gather_list + j );
            j++;
            
            for ( t=d0*block_size[T]; t<(d0+1)*block_size[T]; t++ )
              for ( z=c0*block_size[Z]; z<(c0+1)*block_size[Z]; z++ )
                for ( y=b0*block_size[Y]; y<(b0+1)*block_size[Y]; y++ )
                  for ( x=a0*block_size[X]; x<(a0+1)*block_size[X]; x++ ) {
                    k = lex_index( t, z, y, x, l->local_lattice );
                    field1[k] = i;
                    i++;
                  }
          }
    
    // completely linearly ordered ----> aggregate wise, Schwarz block wise and
    // within Schwarz blocks again linearly
    if ( l->level > 0 ) {
      
      for ( mu=0; mu<4; mu++ ) {
        agg_split[mu] = l->local_lattice[mu]/l->coarsening[mu];
        block_split[mu] = l->coarsening[mu]/l->block_lattice[mu];
        block_size[mu] = l->block_lattice[mu];
      }
      
      i = 0;
      for ( d0=0; d0<agg_split[T]; d0++ )
        for ( c0=0; c0<agg_split[Z]; c0++ )
          for ( b0=0; b0<agg_split[Y]; b0++ )
            for ( a0=0; a0<agg_split[X]; a0++ )
              
              for ( d1=d0*block_split[T]; d1<(d0+1)*block_split[T]; d1++ )
                for ( c1=c0*block_split[Z]; c1<(c0+1)*block_split[Z]; c1++ )
                  for ( b1=b0*block_split[Y]; b1<(b0+1)*block_split[Y]; b1++ )
                    for ( a1=a0*block_split[X]; a1<(a0+1)*block_split[X]; a1++ )
                      
                      for ( t=d1*block_size[T]; t<(d1+1)*block_size[T]; t++ )
                        for ( z=c1*block_size[Z]; z<(c1+1)*block_size[Z]; z++ )
                          for ( y=b1*block_size[Y]; y<(b1+1)*block_size[Y]; y++ )
                            for ( x=a1*block_size[X]; x<(a1+1)*block_size[X]; x++ ) {
                              k = lex_index( t, z, y, x, l->local_lattice );
                              field2[k] = i;
                              i++;
                            }
    } else {
      if ( g.odd_even ) {
        int oe_offset=0, *le = l->local_lattice;
        
        for ( mu=0; mu<4; mu++ )
          oe_offset += (le[mu]*(g.my_coords[mu]/l->comm_offset[mu]))%2;
        oe_offset = oe_offset%2;
        
        i=0;
        for ( t=0; t<le[T]; t++ )
          for ( z=0; z<le[Z]; z++ )
            for ( y=0; y<le[Y]; y++ )
              for ( x=0; x<le[X]; x++ )
                if ( (t+z+y+x+oe_offset)%2 == 0 ) {
                  field2[ lex_index( t, z, y, x, le ) ] = i;
                  i++;
                }
                
        for ( t=0; t<le[T]; t++ )
          for ( z=0; z<le[Z]; z++ )
            for ( y=0; y<le[Y]; y++ )
              for ( x=0; x<le[X]; x++ )
                if ( (t+z+y+x+oe_offset)%2 == 1 ) {
                  field2[ lex_index( t, z, y, x, le ) ] = i;
                  i++;
                }
                        
      } else {
        for ( i=0; i<l->num_inner_lattice_sites; i++ )
          field2[i] = i;
      }
    }      
    
    // build the composition
    for ( i=0; i<l->num_inner_lattice_sites; i++ )
      gs->permutation[ field1[i] ] = field2[i];
    
    FREE( field1, int, l->num_inner_lattice_sites );
    FREE( field2, int, l->num_inner_lattice_sites );
  }  
}


void gathering_PRECISION_free( gathering_PRECISION_struct *gs, level_struct *l ) {
  
  if ( !l->idle ) {
    FREE( gs->gather_list, int, gs->gather_list_length );
    FREE( gs->permutation, int, l->num_inner_lattice_sites );
    FREE( gs->reqs, MPI_Request, gs->gather_list_length );
    FREE( gs->buffer, complex_PRECISION, l->inner_vector_size );
    MPI_Comm_free( &(gs->level_comm) );
    MPI_Group_free( &(gs->level_comm_group) );
  }
  
  FREE( gs->transfer_buffer, complex_PRECISION, gs->dist_inner_lattice_sites * l->num_lattice_site_var );
}


void conf_PRECISION_gather( operator_PRECISION_struct *out, operator_PRECISION_struct *in, level_struct *l ) {
  
  int send_size_hopp = l->gs_PRECISION.dist_inner_lattice_sites * 4 * SQUARE( l->num_lattice_site_var ),
      send_size_clov = l->gs_PRECISION.dist_inner_lattice_sites * ( (l->num_lattice_site_var*(l->num_lattice_site_var+1))/2 );
      
  if ( g.my_rank != l->parent_rank ) {
    MPI_Request req;
    MPI_Isend( in->D, send_size_hopp, MPI_COMPLEX_PRECISION, l->parent_rank, 0, g.comm_cart, &req );
    MPI_Send( in->clover, send_size_clov, MPI_COMPLEX_PRECISION, l->parent_rank, 1, g.comm_cart );
    MPI_Wait( &req, MPI_STATUS_IGNORE );
  } else {
    int i, j, n=l->gs_PRECISION.gather_list_length, s=l->num_inner_lattice_sites,
        t, *pi = l->gs_PRECISION.permutation;
    vector_PRECISION buffer_hopp = NULL, buffer_clov = NULL;
    MPI_Request *hopp_reqs = NULL, *clov_reqs = NULL;
    
    MALLOC( buffer_hopp, complex_PRECISION, n*send_size_hopp );
    MALLOC( buffer_clov, complex_PRECISION, n*send_size_clov );
    MALLOC( hopp_reqs, MPI_Request, n );
    MALLOC( clov_reqs, MPI_Request, n );
    
    PROF_PRECISION_START( _GD_COMM );
    for ( i=1; i<n; i++ ) {
      MPI_Irecv( buffer_hopp+i*send_size_hopp, send_size_hopp, MPI_COMPLEX_PRECISION,
                 l->gs_PRECISION.gather_list[i], 0, g.comm_cart, &(hopp_reqs[i]) );
      MPI_Irecv( buffer_clov+i*send_size_clov, send_size_clov, MPI_COMPLEX_PRECISION,
                 l->gs_PRECISION.gather_list[i], 1, g.comm_cart, &(clov_reqs[i]) );
    }
    PROF_PRECISION_STOP( _GD_COMM, 2*n-2 );
    
    for ( i=0; i<send_size_hopp; i++ )
      buffer_hopp[i] = in->D[i];
    
    for ( i=0; i<send_size_clov; i++ )
      buffer_clov[i] = in->clover[i];
    
    PROF_PRECISION_START( _GD_IDLE );
    for ( i=1; i<n; i++ )
      MPI_Wait( &(hopp_reqs[i]), MPI_STATUS_IGNORE );
    PROF_PRECISION_STOP( _GD_IDLE, n-1 );
    
    t = (send_size_hopp*n)/s;
    for ( i=0; i<s; i++ )
      for ( j=0; j<t; j++ )
        out->D[ t*pi[i] + j ] = buffer_hopp[ t*i + j ];
    
    PROF_PRECISION_START( _GD_IDLE );
    for ( i=1; i<n; i++ )
      MPI_Wait( &(clov_reqs[i]), MPI_STATUS_IGNORE );
    PROF_PRECISION_STOP( _GD_IDLE, n-1 );
    
    t = (send_size_clov*n)/s;
    for ( i=0; i<s; i++ )
      for ( j=0; j<t; j++ )
        out->clover[ t*pi[i] + j ] = buffer_clov[ t*i + j ];
      
    FREE( buffer_hopp, complex_PRECISION, n*send_size_hopp );
    FREE( buffer_clov, complex_PRECISION, n*send_size_clov );
    FREE( hopp_reqs, MPI_Request, n );
    FREE( clov_reqs, MPI_Request, n );
  }
  l->dummy_p_PRECISION.op = out;
  l->dummy_p_PRECISION.shift = 0;
  l->dummy_p_PRECISION.v_start = 0;
  l->dummy_p_PRECISION.v_end = l->inner_vector_size;
  
  if ( g.method == 6 )
    l->dummy_p_PRECISION.eval_operator = g5D_apply_coarse_operator_PRECISION;
  else
    l->dummy_p_PRECISION.eval_operator = apply_coarse_operator_PRECISION;
}


void vector_PRECISION_gather( vector_PRECISION gath, vector_PRECISION dist, level_struct *l ) {
  
  int send_size = l->gs_PRECISION.dist_inner_lattice_sites * l->num_lattice_site_var;
  
  if ( g.my_rank != l->parent_rank ) {
    MPI_Send( dist, send_size, MPI_COMPLEX_PRECISION, l->parent_rank, g.my_rank, g.comm_cart );
  } else {
    int i, j, n=l->gs_PRECISION.gather_list_length, s=l->num_inner_lattice_sites,
        t=l->num_lattice_site_var, *pi = l->gs_PRECISION.permutation;
    vector_PRECISION buffer = l->gs_PRECISION.buffer;

    PROF_PRECISION_START( _GD_COMM );
    for ( i=1; i<n; i++ )
      MPI_Irecv( buffer+i*send_size, send_size, MPI_COMPLEX_PRECISION, l->gs_PRECISION.gather_list[i],
                 l->gs_PRECISION.gather_list[i], g.comm_cart, &(l->gs_PRECISION.reqs[i]) );
    PROF_PRECISION_STOP( _GD_COMM, n-1 );

    for ( i=0; i<send_size; i++ )
      buffer[i] = dist[i];
    
    PROF_PRECISION_START( _GD_IDLE );
    for ( i=1; i<n; i++ )
      MPI_Wait( &(l->gs_PRECISION.reqs[i]), MPI_STATUS_IGNORE );
    PROF_PRECISION_STOP( _GD_IDLE, n-1 );
    // permute data according to desired data layout for parent process
    for ( i=0; i<s; i++ )
      for ( j=0; j<t; j++ )
        gath[ t*pi[i] + j ] = buffer[ t*i + j ];
  }  
}


void vector_PRECISION_distribute( vector_PRECISION dist, vector_PRECISION gath, level_struct *l ) {
  
  int send_size = l->gs_PRECISION.dist_inner_lattice_sites * l->num_lattice_site_var;
  
  if ( g.my_rank != l->parent_rank ) {
    MPI_Recv( dist, send_size, MPI_COMPLEX_PRECISION, l->parent_rank, g.my_rank, g.comm_cart, MPI_STATUS_IGNORE );
  } else {
    int i, j, n=l->gs_PRECISION.gather_list_length, s=l->num_inner_lattice_sites,
        t=l->num_lattice_site_var, *pi = l->gs_PRECISION.permutation;
    vector_PRECISION buffer = l->gs_PRECISION.buffer;
    // permute data according to desired distributed data layout
    for ( i=0; i<s; i++ )
      for ( j=0; j<t; j++ )
        buffer[ t*i+j ] = gath[ t*pi[i]+j ];
    
    PROF_PRECISION_START( _GD_COMM );
    for ( i=1; i<n; i++ )
      MPI_Isend( buffer+i*send_size, send_size, MPI_COMPLEX_PRECISION, l->gs_PRECISION.gather_list[i], 
                 l->gs_PRECISION.gather_list[i], g.comm_cart, &(l->gs_PRECISION.reqs[i]) );
    PROF_PRECISION_STOP( _GD_COMM, n-1 );
      
    for ( i=0; i<send_size; i++ )
      dist[i] = buffer[i];
    
    PROF_PRECISION_START( _GD_IDLE );
    for ( i=1; i<n; i++ )
      MPI_Wait( &(l->gs_PRECISION.reqs[i]), MPI_STATUS_IGNORE );
    PROF_PRECISION_STOP( _GD_IDLE, n-1 );
  }  
}
