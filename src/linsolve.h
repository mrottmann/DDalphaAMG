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

#ifndef LINSOLVE_HEADER
  #define LINSOLVE_HEADER

  #include "linsolve_float.h"
  #include "linsolve_double.h"

  struct Thread;
                                     
  void arnoldi_step_MP( vector_float *V, vector_float *Z, vector_float w,
                        complex_double **H, complex_double* buffer, int j, void (*prec)(),
                        complex_float shift, gmres_float_struct *p, level_struct *l,
                        struct Thread *threading );
                      
  void compute_solution_MP( vector_float x, vector_float *V, complex_double *y,
                            complex_double *gamma, complex_double **H, int j,
                            gmres_float_struct *p, level_struct *l, struct Thread *threading );
  
  void fgmres_MP_struct_init( gmres_MP_struct *p );
  void fgmres_MP_struct_alloc( int m, int n, int vl, double tol, const int prec_kind,
                               void (*precond)(), gmres_MP_struct *p, level_struct* l );
  void fgmres_MP_struct_free( gmres_MP_struct *p );
  
  int fgmres_MP( gmres_MP_struct *p, level_struct *l, struct Thread *threading );
  
#endif
