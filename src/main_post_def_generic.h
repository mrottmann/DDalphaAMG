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

#ifndef MAIN_POST_DEF_PRECISION_HEADER
  #define MAIN_POST_DEF_PRECISION_HEADER
  
  #include "coarse_oddeven_PRECISION.h"
  #include "dirac_PRECISION.h"
  #include "coarse_operator_PRECISION.h"


  static inline void apply_operator_PRECISION( vector_PRECISION output, vector_PRECISION input, gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading ) {
    p->eval_operator( output, input, p->op, l, threading );
    if ( p->shift ) {
      int start, end;
      compute_core_start_end_custom(p->v_start, p->v_end, &start, &end, l, threading, l->num_lattice_site_var );
      vector_PRECISION_saxpy( output, output, input, -p->shift, start, end, l );
    }
  }
  
  static inline void apply_operator_dagger_PRECISION( vector_PRECISION output, vector_PRECISION input, gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading ) {
    if ( l->depth > 0 ) apply_coarse_operator_dagger_PRECISION( output, input, &(l->s_PRECISION.op), l, threading );
    else d_plus_clover_dagger_PRECISION( output, input, p->op, l, threading );
  }
  
#endif
