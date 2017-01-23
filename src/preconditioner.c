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
#include "preconditioner.h"

void preconditioner( vector_double phi, vector_double Dphi, vector_double eta,
                      const int res, level_struct *l, struct Thread *threading ) {
  if ( g.method == 0 )
    vector_double_copy( phi, eta, threading->start_index[l->depth], threading->end_index[l->depth], l );
  else if ( g.method < 5 || g.method == 6 || !g.odd_even ) {
    if ( g.mixed_precision ) {
      trans_float( l->sbuf_float[0], eta, l->s_float.op.translation_table, l, threading );
      vcycle_float( l->sbuf_float[1], NULL, l->sbuf_float[0], res, l, threading );
      trans_back_float( phi, l->sbuf_float[1], l->s_float.op.translation_table, l, threading );
    } else {
      trans_double( l->sbuf_double[0], eta, l->s_double.op.translation_table, l, threading );
      vcycle_double( l->sbuf_double[1], NULL, l->sbuf_double[0], res, l, threading );
      trans_back_double( phi, l->sbuf_double[1], l->s_double.op.translation_table, l, threading );
    }
  } else {
    if ( g.mixed_precision ) {
      START_LOCKED_MASTER(threading)
      l->sp_float.num_restart = l->n_cy;
      l->sp_float.initial_guess_zero = res;
      END_LOCKED_MASTER(threading)
      serial_to_oddeven_float( l->sp_float.b, eta, l, threading );
      if ( g.method == 6 ) {
        g5D_solve_oddeven_float( &(l->sp_float), &(l->oe_op_float), l, threading );
      } else {
        solve_oddeven_float( &(l->sp_float), &(l->oe_op_float), l, threading );
      }
      oddeven_to_serial_float( phi, l->sp_float.x, l, threading );
    } else {
      START_LOCKED_MASTER(threading)
      l->sp_double.num_restart = l->n_cy;
      l->sp_double.initial_guess_zero = res;
      END_LOCKED_MASTER(threading)
      serial_to_oddeven_double( l->sp_double.b, eta, l, threading );
      if ( g.method == 6 ) {
        g5D_solve_oddeven_double( &(l->sp_double), &(l->oe_op_double), l, threading );
      } else {
        solve_oddeven_double( &(l->sp_double), &(l->oe_op_double), l, threading );
      }
      oddeven_to_serial_double( phi, l->sp_double.x, l, threading );
    }
    
  }
  ASSERT( g.mixed_precision != 2 );
  ASSERT( Dphi == NULL );
}

