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

#ifndef COARSE_ODDEVEN_PRECISION_HEADER
  #define COARSE_ODDEVEN_PRECISION_HEADER

  struct Thread;
  
  void coarse_selfcoupling_LU_decomposition_PRECISION( config_PRECISION output, config_PRECISION input, level_struct *l );
  void coarse_perform_fwd_bwd_subs_PRECISION( vector_PRECISION x, vector_PRECISION b, config_PRECISION A, level_struct *l );
  
  void coarse_oddeven_setup_PRECISION( operator_PRECISION_struct *in, int reorder, level_struct *l );
  void coarse_oddeven_re_setup_PRECISION( operator_PRECISION_struct *in, int reorder, level_struct *l, struct Thread *threading );
  void coarse_oddeven_free_PRECISION( level_struct *l );
  
  void coarse_solve_odd_even_PRECISION( gmres_PRECISION_struct *p, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );
  void coarse_apply_schur_complement_PRECISION( vector_PRECISION out, vector_PRECISION in, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );
  void g5D_coarse_solve_odd_even_PRECISION( gmres_PRECISION_struct *p, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );
  void g5D_coarse_apply_schur_complement_PRECISION( vector_PRECISION out, vector_PRECISION in, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );
  
  void coarse_diag_ee_PRECISION( vector_PRECISION y, vector_PRECISION x, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );
  void coarse_diag_oo_PRECISION( vector_PRECISION y, vector_PRECISION x, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );
  void coarse_diag_oo_inv_PRECISION( vector_PRECISION y, vector_PRECISION x, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );
  void coarse_hopping_term_PRECISION( vector_PRECISION out, vector_PRECISION in, operator_PRECISION_struct *op,
                                      const int amount, level_struct *l, struct Thread *threading );
  void coarse_n_hopping_term_PRECISION( vector_PRECISION out, vector_PRECISION in, operator_PRECISION_struct *op,
                                        const int amount, level_struct *l, struct Thread *threading );
  void coarse_hopping_term_PRECISION_vectorized( vector_PRECISION out, vector_PRECISION in, operator_PRECISION_struct *op,
                                      const int amount, level_struct *l, struct Thread *threading );
  void coarse_pn_hopping_term_PRECISION_vectorized( vector_PRECISION out, vector_PRECISION in, operator_PRECISION_struct *op,
                                      const int amount, level_struct *l, int sign, struct Thread *threading );
  void coarse_n_hopping_term_PRECISION_vectorized( vector_PRECISION out, vector_PRECISION in, operator_PRECISION_struct *op,
                                        const int amount, level_struct *l, struct Thread *threading );
  
  void coarse_odd_even_PRECISION_test( vector_PRECISION c4, vector_PRECISION c1, level_struct *l, struct Thread *threading );
  
#endif
