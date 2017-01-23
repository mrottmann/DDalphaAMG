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

#ifndef VAR_TABLE_HEADER
  #define VAR_TABLE_HEADER

  void var_table_init( var_table *t );
  void var_table_insert( var_table *t, var_table_entry e );
  void var_table_free( var_table *t );
  void scan_var( var_table *t, level_struct *l );
  void new_plot_table_line( var_table *t );
  void plot_table( var_table *t );
  
  #define SCAN_VAR( var_pt, kind, start_val, end_val, step_size, mult, name, l ) do{ \
    warning0("SCAN_VAR does not support threading, yet.\n"); \
    kind *tmp_var = (kind*)(var_pt); \
    kind signum = (start_val<end_val)?1:-1; \
    vector_double v = NULL; \
    double norm_v = 0.0, tt0, tt1; \
    vector_double x = (g.mixed_precision==2)?g.p_MP.dp.x:g.p.x; \
    vector_double b = (g.mixed_precision==2)?g.p_MP.dp.b:g.p.b; \
    tt0 = MPI_Wtime(); \
    \
    if ( g.vt.track_error ) { \
      MALLOC( v, complex_double, l->inner_vector_size ); \
      if (g.mixed_precision==2) fgmres_MP( &(g.p_MP), l, no_threading ); \
      else fgmres_double( &(g.p), l, no_threading ); \
      vector_double_copy( v, x, 0, l->inner_vector_size, l ); \
      norm_v = global_norm_double( v, 0, l->inner_vector_size, l, no_threading ); \
    } \
    \
    for ( *tmp_var = (kind)start_val; signum*(*tmp_var) <= signum*((kind)end_val) + EPS_double; \
          *tmp_var = *tmp_var * (kind)(mult?step_size:1)  + (kind)((mult?0:signum)*step_size) ) { \
      prof_init( l ); \
      new_plot_table_line( &(g.vt) ); \
      for ( int i=0; i<g.vt.average_over; i++ ) { \
        g.vt.p_end->values[_TRCKD_VAL] = *tmp_var; \
        parameter_update( l ); \
        if ( g.vt.shift_update ) \
          shift_update( *tmp_var, l, no_threading ); \
        if ( g.vt.re_setup ) { \
          double t0, t1; \
          t0 = MPI_Wtime(); \
          method_re_setup( l, no_threading ); \
          method_update( g.setup_iter[0], l, no_threading ); \
          t1 = MPI_Wtime();\
          if ( g.vt.p_end != NULL ) g.vt.p_end->values[_STP_TIME] += (t1-t0) / ((double)g.vt.average_over); \
        } \
        printf0("scanning variable \"%s\", value: %lf, run %d of %d\n", name, (double)(*tmp_var), i+1, g.vt.average_over ); \
        if ( g.vt.track_error ) { \
          apply_operator_double( b, v, &(g.p), l, no_threading ); \
          vector_double_define( x, 0, 0, l->inner_vector_size, l ); \
          if ( g.vt.track_cgn_error ) { \
            ASSERT( g.method >=0 && g.p.restart_length >= 4 ); \
            vector_double_define( x, 0, 0, l->inner_vector_size, l ); \
            cgn_double( &(g.p), l, no_threading ); \
            vector_double_minus( x, x, v, 0, l->inner_vector_size, l ); \
            g.vt.p_end->values[_CGNR_ERR] += ( global_norm_double( x, 0, l->inner_vector_size, l, no_threading ) / norm_v ) / ((double)g.vt.average_over); \
            printf0("CGN: error norm: %le\n", g.vt.p_end->values[_CGNR_ERR] ); \
            vector_double_define( x, 0, 0, l->inner_vector_size, l ); \
            } \
        } else {\
          rhs_define( b, l, no_threading );\
        } \
        vector_double_define( x, 0, 0, l->inner_vector_size, l ); \
        if (g.mixed_precision==2) fgmres_MP( &(g.p_MP), l, no_threading ); \
        else fgmres_double( &(g.p), l, no_threading ); \
        if ( i == g.vt.average_over-1 ) prof_print( l ); \
        if ( g.vt.track_error ) { \
          vector_double_minus( x, x, v, 0, l->inner_vector_size, l ); \
          g.vt.p_end->values[_SLV_ERR] += ( global_norm_double( x, 0, l->inner_vector_size, l, no_threading ) / norm_v ) / ((double)g.vt.average_over); \
        } \
      } \
    } \
    if ( g.vt.track_error ) { \
      FREE( v, complex_double, l->inner_vector_size ); \
    } \
    tt1 = MPI_Wtime(); \
    printf0("\n\ntotal time for parameter scan: %d minutes and %d seconds\n", \
    (int)((tt1-tt0)/60), ((int)round(tt1-tt0))%60 + 1 ); \
  }while(0)
  
#endif
