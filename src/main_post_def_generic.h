/*
 * Copyright (C) 2016, Matthias Rottmann, Artur Strebel, Simon Heybrock, Simone Bacchio, Bjoern Leder.
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

  }
  
  static inline void apply_operator_dagger_PRECISION( vector_PRECISION output, vector_PRECISION input, gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading ) {

#ifdef HAVE_TM1p1
    if( g.n_flavours == 2 ) {
      tau1_gamma5_PRECISION( l->vbuf_PRECISION[6], input, l, threading );
    } else
#endif
      {
        gamma5_PRECISION( l->vbuf_PRECISION[6], input, l, threading );
#ifdef HAVE_TM
        //TODO: change_mu_sign_PRECISION( p->op, l, threading );
#endif
      }

    apply_operator_PRECISION( l->vbuf_PRECISION[7], l->vbuf_PRECISION[6], p, l, threading );

#ifdef HAVE_TM1p1
    if( g.n_flavours == 2 ) {
      tau1_gamma5_PRECISION( output, l->vbuf_PRECISION[7], l, threading );
    } else
#endif
      {
        gamma5_PRECISION( output, l->vbuf_PRECISION[7], l, threading );
#ifdef HAVE_TM
        //TODO: change_mu_sign_PRECISION( p->op, l, threading );
#endif
      }
    
  }

  static inline void test0_PRECISION( char* format, int depth, PRECISION test ) {
    if ( g.my_rank == 0 && g.print >= 0 ) {
      if ( test > EPS_PRECISION )
        printf("\x1b[31m");
      printf(format, depth, test);
      if ( test > EPS_PRECISION )
        printf("\x1b[0m");
      if ( test > g.test )
        g.test = test;
      fflush(0);
    }
  }
  
#endif
