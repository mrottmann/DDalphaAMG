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

#ifndef LINSOLVE_PRECISION_HEADER
  #define LINSOLVE_PRECISION_HEADER

  struct Thread;
  
  void fgmres_PRECISION_struct_init( gmres_PRECISION_struct *p );
  void fgmres_PRECISION_struct_alloc( int m, int n, int vl, PRECISION tol, const int type, const int prec_kind,
                                      void (*precond)(), void (*eval_op)(), gmres_PRECISION_struct *p, level_struct* l );
  void fgmres_PRECISION_struct_free( gmres_PRECISION_struct *p, level_struct *l );
  
  int fgmres_PRECISION( gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading );
  void fgcr_PRECISION( gmres_PRECISION_struct *p, level_struct *l );
  void cgn_PRECISION( gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading );
  void bicgstab_PRECISION( gmres_PRECISION_struct *ps, level_struct *l, struct Thread *threading );
  void local_minres_PRECISION( vector_PRECISION phi, vector_PRECISION eta, vector_PRECISION latest_iter,
                               int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading );
  int arnoldi_step_PRECISION( vector_PRECISION *V, vector_PRECISION *Z, vector_PRECISION w,
                               complex_PRECISION **H, complex_PRECISION* buffer, int j, void (*prec)(),
                               complex_PRECISION shift, gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading );
  void qr_update_PRECISION( complex_PRECISION **H, complex_PRECISION *s,
                            complex_PRECISION *c, complex_PRECISION *gamma, int j,
                            level_struct *l, struct Thread *threading );
  void compute_solution_PRECISION( vector_PRECISION x, vector_PRECISION *V, complex_PRECISION *y, complex_PRECISION *gamma,
                                   complex_PRECISION **H, int j, int ol, gmres_PRECISION_struct *p, level_struct *l, struct Thread *threading );
  
#endif
