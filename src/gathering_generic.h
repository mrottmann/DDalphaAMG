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

#ifndef GATHERING_PRECISION_HEADER
  #define GATHERING_PRECISION_HEADER
  
  void gathering_PRECISION_next_level_init( gathering_PRECISION_struct *gs, level_struct *l );
  void gathering_PRECISION_setup( gathering_PRECISION_struct *gs, level_struct *l );
  void gathering_PRECISION_free( gathering_PRECISION_struct *gs, level_struct *l );
  
  void conf_PRECISION_gather( operator_PRECISION_struct *out, operator_PRECISION_struct *in, level_struct *l );
  void vector_PRECISION_gather( vector_PRECISION gath, vector_PRECISION dist, level_struct *l );
  void vector_PRECISION_distribute( vector_PRECISION dist, vector_PRECISION gath, level_struct *l );
  
  void distribution_PRECISION_next_level_test( level_struct *l );
  
#endif 
