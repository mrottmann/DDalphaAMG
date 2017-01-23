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

#ifndef LINALG_HEADER
  #define LINALG_HEADER

  struct Thread;
  
  void vector_double_multi_saxpy( vector_double z, vector_double *V, complex_double *alpha,
                               int sign, int count, int start, int end, level_struct *l );
  
  void vector_float_multi_saxpy( vector_float z, vector_float *V, complex_float *alpha,
                                 int sign, int count, int start, int end, level_struct *l );
  
  void process_multi_inner_product_MP( int count, complex_double *results, vector_float *phi,
                                       vector_float psi, int start, int end, level_struct *l,
                                       struct Thread *threading );
                                       
  double global_norm_MP( vector_float x, int start, int end, level_struct *l, struct Thread *threading );
  
#endif
