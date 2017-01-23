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

#ifndef DIRAC_HEADER
  #define DIRAC_HEADER

struct Thread;

  typedef complex_double ******SU3_storage;

  void dirac_setup( config_double hopp, config_double clover, level_struct *l );
  double *dirac_setup_get_gauge_pointer();
  double *dirac_setup_get_clover_pointer();
  
  void SU3_storage_alloc( SU3_storage *U, level_struct *l );
  void SU3_storage_free( SU3_storage *U, level_struct *l );
  void SU3_ghost_update( SU3_storage *U, level_struct *l );
  
  void spin_alloc( int num_spin, int n );
  void spin_free( int num_spin, int n );
  void spin_define( void );
  
  void mat_alloc( complex_double **A, int n );
  void mat_free( complex_double **A, int n );
  void zeroMat( complex_double *A, int n );
  
  void apply_anti_pbc( SU3_storage U, level_struct *l );
  void calc_plaq( SU3_storage U, level_struct *l );
  
  void Qdiff( complex_double *q_store, int mu, int nu, int t, int z, int y, int x, SU3_storage U );
  void set_clover( complex_double *q_store, int mu, int nu, int index, config_double clover );
  
  void define_odd_even_table( level_struct *l );
  void scale_clover( operator_double_struct *op, double scale_even, double scale_odd, level_struct *l );
  void shift_update( complex_double shift, level_struct *l, struct Thread *threading );
  void g5D_shift_update( complex_double shift, level_struct *l, struct Thread *threading );

#endif
