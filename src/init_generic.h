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

#ifndef INIT_PRECISION_HEADER
  #define INIT_PRECISION_HEADER

struct Thread;
  
  void prof_PRECISION_init( level_struct *l );
  double prof_PRECISION_print( level_struct *l );
  
  void fine_level_PRECISION_alloc( level_struct *l );
  void fine_level_PRECISION_free( level_struct *l );
  void next_level_PRECISION_setup( level_struct *l );
  void next_level_PRECISION_free( level_struct *l );
  void level_PRECISION_init( level_struct *l );
  
  void vcycle_timing_PRECISION( int n, level_struct *l, struct Thread *threading );

#endif
