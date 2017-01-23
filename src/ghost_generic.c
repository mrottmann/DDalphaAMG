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

void negative_sendrecv_PRECISION( vector_PRECISION phi, const int mu, comm_PRECISION_struct *c, level_struct *l ) {
  // send dir = -1
  if( l->global_splitting[mu] > 1 ) {    
    
    int i, j, num_boundary_sites = c->num_boundary_sites[2*mu+1], boundary_start,
        *boundary_table = c->boundary_table[2*mu+1], n = l->num_lattice_site_var;
    
    vector_PRECISION buffer, tmp_pt, buffer_pt;
    
    boundary_start = l->num_inner_lattice_sites;
    for ( i=0; i<mu; i++ )
      boundary_start += c->num_boundary_sites[2*i];

    buffer = l->vbuf_PRECISION[8]+n*(boundary_start-l->num_inner_lattice_sites);
    buffer_pt = buffer;
    
    for ( i=0; i<num_boundary_sites; i++ ) {
      tmp_pt = phi + n*boundary_table[i];
      for ( j=0; j<n; j++, buffer_pt++, tmp_pt++ )
        *buffer_pt = *tmp_pt;
    }
    
    MPI_Irecv( phi+n*boundary_start, n*num_boundary_sites, MPI_COMPLEX_PRECISION,
               l->neighbor_rank[2*mu], 2*mu+1, g.comm_cart, &(c->rreqs[2*mu+1]) );
    MPI_Isend( buffer, n*num_boundary_sites, MPI_COMPLEX_PRECISION,
               l->neighbor_rank[2*mu+1], 2*mu+1, g.comm_cart, &(c->sreqs[2*mu+1]) );
  }
}


void negative_sendrecv_PRECISION_vectorized( complex_PRECISION *phi, const int mu, comm_PRECISION_struct *c,
                                             level_struct *l, int count, complex_PRECISION *buffer ) {
  // send dir = -1
  if( l->global_splitting[mu] > 1 ) {

    int i, j, num_boundary_sites = c->num_boundary_sites[2*mu+1], boundary_start,
        *boundary_table = c->boundary_table[2*mu+1], n = l->num_lattice_site_var;

    complex_PRECISION *tmp_pt;
    complex_PRECISION *buffer_pt;

    boundary_start = l->num_inner_lattice_sites;
    for ( i=0; i<mu; i++ )
      boundary_start += c->num_boundary_sites[2*i];

    buffer_pt = buffer;

    for ( i=0; i<num_boundary_sites; i++ ) {
      tmp_pt = phi + count*n*boundary_table[i];
      for ( j=0; j<count*n; j++, buffer_pt++, tmp_pt++ )
        *buffer_pt = *tmp_pt;
    }

    MPI_Irecv( phi+count*n*boundary_start, count*n*num_boundary_sites, MPI_COMPLEX_PRECISION,
               l->neighbor_rank[2*mu], 2*mu+1, g.comm_cart, &(c->rreqs[2*mu+1]) );
    MPI_Isend( buffer, count*n*num_boundary_sites, MPI_COMPLEX_PRECISION,
               l->neighbor_rank[2*mu+1], 2*mu+1, g.comm_cart, &(c->sreqs[2*mu+1]) );
  }
}


void negative_wait_PRECISION( const int mu, comm_PRECISION_struct *c, level_struct *l ) {
 
  if( l->global_splitting[mu] > 1 ) {
    MPI_Wait( &(c->sreqs[2*mu+1]), MPI_STATUS_IGNORE );
    MPI_Wait( &(c->rreqs[2*mu+1]), MPI_STATUS_IGNORE );
  }
}


void ghost_alloc_PRECISION( int buffer_size, comm_PRECISION_struct *c, level_struct *l ) {
  
  int mu, nu, factor=1;
  
  if ( l->depth > 0 ) {
    c->offset = l->num_lattice_site_var;
  } else {
    c->offset = l->num_lattice_site_var/2;
    if ( g.method < 5 )
      factor = 2;
  }
  
  if ( buffer_size <= 0 ) {
    c->comm_start[0] = c->offset*l->num_inner_lattice_sites;
    c->comm_start[1] = c->offset*l->num_inner_lattice_sites;
    for ( mu=0; mu<4; mu++ ) {
      if ( mu > 0 ) {
        c->comm_start[2*mu] = c->comm_start[2*(mu-1)] + buffer_size;
        c->comm_start[2*mu+1] = c->comm_start[2*(mu-1)+1] + buffer_size;
      }
      buffer_size = c->offset;
      for ( nu=0; nu<4; nu++ ) {
        if ( nu != mu ) {
          buffer_size *= l->local_lattice[nu];
        }
      }
      c->length[2*mu] = buffer_size;
      c->length[2*mu+1] = buffer_size;
      c->max_length[mu] = factor*buffer_size;
      MALLOC( c->buffer[2*mu], complex_PRECISION, factor*buffer_size );
      MALLOC( c->buffer[2*mu+1], complex_PRECISION, factor*buffer_size );
      c->in_use[2*mu] = 0;
      c->in_use[2*mu+1] = 0;
    }
  } else {
    for ( mu=0; mu<4; mu++ ) {
      c->max_length[mu] = buffer_size;
      MALLOC( c->buffer[2*mu], complex_PRECISION, buffer_size );
      MALLOC( c->buffer[2*mu+1], complex_PRECISION, buffer_size );
    }
  }
  
  if ( l->vbuf_PRECISION[8] == NULL ) {
    MALLOC( l->vbuf_PRECISION[8], complex_PRECISION, l->vector_size );
  }
}


void ghost_free_PRECISION( comm_PRECISION_struct *c, level_struct *l ) {
  
  int mu;
  
  for ( mu=0; mu<4; mu++ ) {    
    FREE( c->buffer[2*mu], complex_PRECISION, c->max_length[mu] );
    FREE( c->buffer[2*mu+1], complex_PRECISION, c->max_length[mu] );
  }
  
  if ( l->vbuf_PRECISION[8] != NULL ) {
    FREE( l->vbuf_PRECISION[8], complex_PRECISION, l->vector_size );
  }
}


void ghost_sendrecv_init_PRECISION( const int type, comm_PRECISION_struct *c, level_struct *l ) {
  
  int mu; 
  
  if ( type == _COARSE_GLOBAL ) {
    c->comm = 1;
    for ( mu=0; mu<4; mu++ ) {
      ASSERT( c->in_use[2*mu] == 0 );
      ASSERT( c->in_use[2*mu+1] == 0 );
    }
  }
}


void ghost_sendrecv_PRECISION( vector_PRECISION phi, const int mu, const int dir,
                               comm_PRECISION_struct *c, const int amount, level_struct *l ) {
  // does not allow sending in both directions at the same time
  if( l->global_splitting[mu] > 1 ) {
    
    int i, j, *table=NULL, mu_dir = 2*mu-MIN(dir,0), offset = c->offset,
        length[2] = {0,0}, comm_start = 0, table_start = 0;
    vector_PRECISION buffer, phi_pt;
    
    if ( amount == _FULL_SYSTEM ) {
      length[0] = (c->num_boundary_sites[2*mu])*offset;
      length[1] = (c->num_boundary_sites[2*mu+1])*offset;
      comm_start = c->comm_start[mu_dir];
      table_start = 0;
    } else if ( amount == _EVEN_SITES ) {
      length[0] = c->num_even_boundary_sites[2*mu]*offset;
      length[1] = c->num_even_boundary_sites[2*mu+1]*offset;
      comm_start = c->comm_start[mu_dir];
      table_start = 0;
    } else if ( amount == _ODD_SITES ) {
      length[0] = c->num_odd_boundary_sites[2*mu]*offset;
      length[1] = c->num_odd_boundary_sites[2*mu+1]*offset;
      comm_start = c->comm_start[mu_dir]+c->num_even_boundary_sites[mu_dir]*offset;
      table_start = c->num_even_boundary_sites[mu_dir];
    }
    
    ASSERT( c->in_use[mu_dir] == 0 );
    c->in_use[mu_dir] = 1;
    
    if ( MAX(length[0],length[1]) > c->max_length[mu] ) {
      printf("CAUTION: my_rank: %d, not enough comm buffer\n", g.my_rank ); fflush(0);
      ghost_free_PRECISION( c, l );
      ghost_alloc_PRECISION( MAX(length[0],length[1]), c, l );
    }
    
    buffer = (vector_PRECISION)c->buffer[mu_dir];
    
    // dir = senddir
    if ( dir == 1 ) {
      // data to be communicated is stored serially in the vector phi
      // recv target is a buffer
      // afterwards (in ghost_wait) the data has to be distributed onto the correct sites
      // touching the respective boundary in -mu direction
      
      phi_pt = phi + comm_start;
      if ( length[1] > 0 ) {
        PROF_PRECISION_START( _OP_COMM );
        MPI_Irecv( buffer, length[1], MPI_COMPLEX_PRECISION,
                   l->neighbor_rank[2*mu+1], 2*mu, g.comm_cart, &(c->rreqs[2*mu]) );
        PROF_PRECISION_STOP( _OP_COMM, 1 );
      }
      if ( length[0] > 0 ) {
        PROF_PRECISION_START( _OP_COMM );
        MPI_Isend( phi_pt, length[0], MPI_COMPLEX_PRECISION,
                   l->neighbor_rank[2*mu], 2*mu, g.comm_cart, &(c->sreqs[2*mu]) );
        PROF_PRECISION_STOP( _OP_COMM, 0 );
      }
      
      
    } else if ( dir == -1 ) {
      // data to be communicated is stored on the sites touching the boundary in -mu direction
      // this data is gathered in a buffer in the correct ordering
      // which is required on the boundary of the vector phi
      int num_boundary_sites = length[1]/offset;
      
      table = c->boundary_table[2*mu+1]+table_start;
      for ( j=0; j<num_boundary_sites; j++ ) {
        phi_pt = phi + table[j]*offset;
        for ( i=0; i<offset; i++ ) {
          buffer[i] = phi_pt[i];
        }
        buffer += offset;
      }
      
      buffer = (vector_PRECISION)c->buffer[mu_dir];      
      phi_pt = phi + comm_start;
      
      if ( length[0] > 0 ) {
        PROF_PRECISION_START( _OP_COMM );
        MPI_Irecv( phi_pt, length[0], MPI_COMPLEX_PRECISION,
                   l->neighbor_rank[2*mu], 2*mu+1, g.comm_cart, &(c->rreqs[2*mu+1]) );
        PROF_PRECISION_STOP( _OP_COMM, 1 );
      }
      if ( length[1] > 0 ) {
        PROF_PRECISION_START( _OP_COMM );
        MPI_Isend( buffer, length[1], MPI_COMPLEX_PRECISION,
                   l->neighbor_rank[2*mu+1], 2*mu+1, g.comm_cart, &(c->sreqs[2*mu+1]) );
        PROF_PRECISION_STOP( _OP_COMM, 0 );
      }
      
    } else ASSERT( dir == 1 || dir == -1 );
  }
}


void ghost_wait_PRECISION( vector_PRECISION phi, const int mu, const int dir,
                           comm_PRECISION_struct *c, const int amount, level_struct *l ) {
  
  if( l->global_splitting[mu] > 1 ) {    
    int mu_dir = 2*mu-MIN(dir,0);
    int i, j, *table, offset = c->offset, length[2]={0,0}, table_start = 0;
    vector_PRECISION buffer, phi_pt;
      
    if ( amount == _FULL_SYSTEM ) {
      length[0] = (c->num_boundary_sites[2*mu])*offset;
      length[1] = (c->num_boundary_sites[2*mu+1])*offset;
      table_start = 0;
    } else if ( amount == _EVEN_SITES ) {
      length[0] = c->num_even_boundary_sites[2*mu]*offset;
      length[1] = c->num_even_boundary_sites[2*mu+1]*offset;
      table_start = 0;
    } else if ( amount == _ODD_SITES ) {
      length[0] = c->num_odd_boundary_sites[2*mu]*offset;
      length[1] = c->num_odd_boundary_sites[2*mu+1]*offset;
      table_start = c->num_even_boundary_sites[mu_dir];
    }
    
    ASSERT( c->in_use[mu_dir] == 1 );
    
    if ( dir == 1 ) {
      
      int num_boundary_sites = length[0]/offset;
      
      buffer = (vector_PRECISION)c->buffer[mu_dir];      
      table = c->boundary_table[2*mu+1] + table_start;
      
      if ( length[0] > 0 ) {
        PROF_PRECISION_START( _OP_IDLE );
        MPI_Wait( &(c->sreqs[2*mu]), MPI_STATUS_IGNORE );
        PROF_PRECISION_STOP( _OP_IDLE, 0 );
      }
      if ( length[1] > 0 ) {
        PROF_PRECISION_START( _OP_IDLE );
        MPI_Wait( &(c->rreqs[2*mu]), MPI_STATUS_IGNORE );
        PROF_PRECISION_STOP( _OP_IDLE, 1 );
      }
      
      if ( l->depth == 0 ) {
        for ( j=0; j<num_boundary_sites; j++ ) {
          phi_pt = phi + table[j]*offset;
          
          for ( i=0; i<offset; i++ )
            phi_pt[i] = buffer[i];
          
          buffer += offset;
        }
      } else {
        for ( j=0; j<num_boundary_sites; j++ ) {
          phi_pt = phi + table[j]*offset;
          
          for ( i=0; i<offset; i++ )
            phi_pt[i] += buffer[i];
          
          buffer += offset;
        }
      }
    } else if ( dir == -1 ) {
      
      if ( length[1] > 0 ) {
        PROF_PRECISION_START( _OP_IDLE );
        MPI_Wait( &(c->sreqs[2*mu+1]), MPI_STATUS_IGNORE );
        PROF_PRECISION_STOP( _OP_IDLE, 0 );
      }
      if ( length[0] > 0 ) {
        PROF_PRECISION_START( _OP_IDLE );
        MPI_Wait( &(c->rreqs[2*mu+1]), MPI_STATUS_IGNORE );
        PROF_PRECISION_STOP( _OP_IDLE, 1 );
      }
      
    } else ASSERT( dir == 1 || dir == -1 );
    
    c->in_use[mu_dir] = 0;
  }
}


void ghost_update_PRECISION( vector_PRECISION phi, const int mu, const int dir, comm_PRECISION_struct *c, level_struct *l ) {
  
  if( l->global_splitting[mu] > 1 ) {
    int i, j, mu_dir = 2*mu-MIN(dir,0), nu, inv_mu_dir = 2*mu+1+MIN(dir,0), length, *table=NULL,
        comm_start, num_boundary_sites, site_var;
    vector_PRECISION buffer, recv_pt, phi_pt;
    
    site_var = l->num_lattice_site_var;
    length = c->num_boundary_sites[mu_dir]*l->num_lattice_site_var;
    num_boundary_sites = c->num_boundary_sites[mu_dir];
    buffer = (vector_PRECISION)c->buffer[mu_dir];
    
    if ( dir == -1 )
      comm_start = l->vector_size;
    else
      comm_start = l->inner_vector_size;
    for ( nu=0; nu<mu; nu++ ) {
      comm_start += c->num_boundary_sites[2*nu]*l->num_lattice_site_var;
    }
    
    ASSERT( c->in_use[mu_dir] == 0 );
    c->in_use[mu_dir] = 1;
    
    recv_pt = phi + comm_start;
    if ( length > 0 ) {
      PROF_PRECISION_START( _OP_COMM );
      MPI_Irecv( recv_pt, length, MPI_COMPLEX_PRECISION,
                 l->neighbor_rank[mu_dir], mu_dir, g.comm_cart, &(c->rreqs[mu_dir]) );
      PROF_PRECISION_STOP( _OP_COMM, 1 );
    }
    
    table = c->boundary_table[inv_mu_dir];
    for ( j=0; j<num_boundary_sites; j++ ) {
      phi_pt = phi + table[j]*site_var;
      
      for ( i=0; i<site_var; i++ ) {
        buffer[i] = phi_pt[i];
      }
      buffer += site_var;
    }
    buffer = (vector_PRECISION)c->buffer[mu_dir];
    
    if ( length > 0 ) {
      PROF_PRECISION_START( _OP_COMM );
      MPI_Isend( buffer, length, MPI_COMPLEX_PRECISION,
                 l->neighbor_rank[inv_mu_dir], mu_dir, g.comm_cart, &(c->sreqs[mu_dir]) );
      PROF_PRECISION_STOP( _OP_COMM, 0 );
    }
  }
}


void ghost_update_wait_PRECISION( vector_PRECISION phi, const int mu, const int dir, comm_PRECISION_struct *c, level_struct *l ) {
  
  if( l->global_splitting[mu] > 1 ) {
    int mu_dir = 2*mu-MIN(dir,0), length = c->num_boundary_sites[mu_dir]*l->num_lattice_site_var;
    
    ASSERT( c->in_use[mu_dir] == 1 );
      
    if ( length > 0 ) {
      PROF_PRECISION_START( _OP_IDLE );
      MPI_Wait( &(c->sreqs[mu_dir]), MPI_STATUS_IGNORE );
      MPI_Wait( &(c->rreqs[mu_dir]), MPI_STATUS_IGNORE );
      PROF_PRECISION_STOP( _OP_IDLE, 1 );
    }
    c->in_use[mu_dir] = 0;
  }
}
