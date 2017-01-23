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

#ifndef GHOST_PRECISION_HEADER
  #define GHOST_PRECISION_HEADER
    
  void negative_sendrecv_PRECISION( vector_PRECISION phi, const int mu, comm_PRECISION_struct *c, level_struct *l );
  
  // as negative_sendrecv_PRECISION, but for count vectors stored in phi in vector-fused data layout
  // buffer must be big enough to hold the surface data for count vectors (in one direction)
  void negative_sendrecv_PRECISION_vectorized( complex_PRECISION *phi, const int mu, comm_PRECISION_struct *c, level_struct *l, int count, complex_PRECISION *buffer );
  void negative_wait_PRECISION( const int mu, comm_PRECISION_struct *c, level_struct *l );
  
  void ghost_alloc_PRECISION( int buffer_size, comm_PRECISION_struct *c, level_struct *l );
  void ghost_free_PRECISION( comm_PRECISION_struct *c, level_struct *l );
  void ghost_sendrecv_init_PRECISION( const int type, comm_PRECISION_struct *c, level_struct *l );
  void ghost_sendrecv_PRECISION( vector_PRECISION phi, const int mu, const int dir,
                                 comm_PRECISION_struct *c, const int amount, level_struct *l );
  void ghost_wait_PRECISION( vector_PRECISION phi, const int mu, const int dir,
                             comm_PRECISION_struct *c, const int amount, level_struct *l );
  
  void ghost_update_PRECISION( vector_PRECISION phi, const int mu, const int dir, comm_PRECISION_struct *c, level_struct *l );
  void ghost_update_wait_PRECISION( vector_PRECISION phi, const int mu, const int dir, comm_PRECISION_struct *c, level_struct *l );

#endif
