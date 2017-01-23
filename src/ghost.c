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
 
void neighbor_define( level_struct *l ) {
  
  int mu, neighbor_coords[4];
  
  for ( mu=0; mu<4; mu++ ) {
    neighbor_coords[mu] = g.my_coords[mu];
  }
  
  for ( mu=0; mu<4; mu++ ) {    
    neighbor_coords[mu]+=l->comm_offset[mu];
    MPI_Cart_rank( g.comm_cart, neighbor_coords, &(l->neighbor_rank[2*mu]) );
    neighbor_coords[mu]-=2*l->comm_offset[mu];
    MPI_Cart_rank( g.comm_cart, neighbor_coords, &(l->neighbor_rank[2*mu+1]) );
    neighbor_coords[mu]+=l->comm_offset[mu];
  }
}


void predefine_rank( void ) {
  MPI_Comm_rank( MPI_COMM_WORLD, &(g.my_rank) );
}


void cart_define( level_struct *l ) {
  
  int mu, num_processes;
  
  MPI_Comm_size( MPI_COMM_WORLD, &num_processes ); 
  if (num_processes != g.num_processes) {
    error0("Error: Number of processes has to be %d\n", g.num_processes);
  }
  MPI_Cart_create( MPI_COMM_WORLD, 4, l->global_splitting, l->periodic_bc, 1, &(g.comm_cart) );
  MPI_Comm_rank( g.comm_cart, &(g.my_rank) );
  MPI_Cart_coords( g.comm_cart, g.my_rank, 4, g.my_coords );
  
  for ( mu=0; mu<4; mu++ ) {
    l->num_processes_dir[mu] = l->global_lattice[mu]/l->local_lattice[mu];
    l->comm_offset[mu] = 1;
  }

  neighbor_define( l );
  MPI_Comm_group( g.comm_cart, &(g.global_comm_group) );
}


void cart_free( level_struct *l ) {
  MPI_Group_free( &(g.global_comm_group) );
  MPI_Comm_free( &(g.comm_cart) );  
}
