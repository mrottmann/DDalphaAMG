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

#ifndef SSE_INTERPOLATION_PRECISION_HEADER
  #define SSE_INTERPOLATION_PRECISION_HEADER  
  
  #ifdef SSE
  void interpolation_PRECISION_alloc( level_struct *l );
  void interpolation_PRECISION_free( level_struct *l );
  void interpolation_PRECISION_dummy_alloc( level_struct *l );
  void interpolation_PRECISION_dummy_free( level_struct *l );
  
  void interpolate_PRECISION( vector_PRECISION phi, vector_PRECISION phi_c, level_struct *l, Thread *threading );
  void interpolate3_PRECISION( vector_PRECISION phi, vector_PRECISION phi_c, level_struct *l, Thread *threading );
  void restrict_PRECISION( vector_PRECISION phi_c, vector_PRECISION phi, level_struct *l, Thread *threading );
#endif
  
#endif