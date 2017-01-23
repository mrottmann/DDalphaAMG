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

#ifndef DATA_LAYOUT_HEADER
  #define DATA_LAYOUT_HEADER
  
  void data_layout_init( level_struct *l );
  void define_eot( int *eot, int *N, level_struct *l );
  void define_eo_bt( int **bt, int *eot, int *n_ebs, int *n_obs, int *n_bs, int *N, level_struct *l );
  void define_nt_bt_tt( int *nt, int *backward_nt, int **bt, int *tt, int *it, int *dt, level_struct *l );
  
  static inline int lex_index( int t, int z, int y, int x, int N[3] ) {  
    return x + N[X]*( y + N[Y]*(z + N[Z]*t ) );
  }
  
  static inline int lex_mod_index( int t, int z, int y, int x, int N[4] ) {
    return (x+N[X])%N[X] + N[X]*( (y+N[Y])%N[Y] + N[Y]*((z+N[Z])%N[Z] + N[Z]*((t+N[T])%N[T]) ) );    
  }
  
  static inline int site_index( int t, int z, int y, int x, int N[3], int *index_table ) {
    return index_table[ lex_index( t, z, y, x, N ) ];
  }
  
  static inline int site_mod_index( int t, int z, int y, int x, int N[4], int *index_table ) {
    return index_table[ lex_mod_index( t, z, y, x, N ) ];
  }

#endif
