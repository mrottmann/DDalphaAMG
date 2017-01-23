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

#ifndef DIRAC_PRECISION_HEADER
  #define DIRAC_PRECISION_HEADER

  struct Thread;
  
  void gamma5_PRECISION( vector_PRECISION eta, vector_PRECISION phi, level_struct *l, struct Thread *threading );
  
  void clover_PRECISION( vector_PRECISION eta, vector_PRECISION phi, config_PRECISION clover, int length,
                         level_struct *l, struct Thread *threading );
  
  void d_plus_clover_PRECISION( vector_PRECISION eta, vector_PRECISION phi, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );
  void d_plus_clover_dagger_PRECISION( vector_PRECISION eta, vector_PRECISION phi, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );
  void g5D_plus_clover_PRECISION( vector_PRECISION eta, vector_PRECISION phi, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading );
  void block_d_plus_clover_PRECISION( vector_PRECISION eta, vector_PRECISION phi, int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading );
  void d_plus_clover_aggregate_PRECISION( vector_PRECISION eta1, vector_PRECISION eta2, vector_PRECISION phi, schwarz_PRECISION_struct *s, level_struct *l );
  void d_neighbor_aggregate_PRECISION( vector_PRECISION eta1, vector_PRECISION eta2, vector_PRECISION phi, const int mu, schwarz_PRECISION_struct *s, level_struct *l );
  void operator_updates_PRECISION( level_struct *l );
  void shift_update_PRECISION( operator_PRECISION_struct *op, complex_PRECISION shift, level_struct *l, struct Thread *threading );
  void g5D_shift_update_PRECISION( operator_PRECISION_struct *op, complex_PRECISION shift, level_struct *l, struct Thread *threading );
  
  static inline void zero12_PRECISION( const vector_PRECISION phi ) {
    phi[ 0] = _COMPLEX_PRECISION_ZERO;
    phi[ 1] = _COMPLEX_PRECISION_ZERO;
    phi[ 2] = _COMPLEX_PRECISION_ZERO;
    phi[ 3] = _COMPLEX_PRECISION_ZERO;
    phi[ 4] = _COMPLEX_PRECISION_ZERO;
    phi[ 5] = _COMPLEX_PRECISION_ZERO;
    phi[ 6] = _COMPLEX_PRECISION_ZERO;
    phi[ 7] = _COMPLEX_PRECISION_ZERO;
    phi[ 8] = _COMPLEX_PRECISION_ZERO;
    phi[ 9] = _COMPLEX_PRECISION_ZERO;
    phi[10] = _COMPLEX_PRECISION_ZERO;
    phi[11] = _COMPLEX_PRECISION_ZERO;
  }

  // eta = D*phi
  static inline void mvm_PRECISION( const vector_PRECISION eta, const complex_PRECISION *D, const vector_PRECISION phi ) {
    eta[0]  = D[0]*phi[0];
    eta[0] += D[1]*phi[1];
    eta[0] += D[2]*phi[2];
    eta[1]  = D[3]*phi[0];
    eta[1] += D[4]*phi[1];
    eta[1] += D[5]*phi[2];
    eta[2]  = D[6]*phi[0];
    eta[2] += D[7]*phi[1];
    eta[2] += D[8]*phi[2];
  }
  
  // eta = D**H*phi
  static inline void mvmh_PRECISION( const vector_PRECISION eta, const complex_PRECISION *D, const vector_PRECISION phi ) {
    eta[0]  = conj_PRECISION(D[0])*phi[0];
    eta[1]  = conj_PRECISION(D[1])*phi[0];
    eta[2]  = conj_PRECISION(D[2])*phi[0];
    eta[0] += conj_PRECISION(D[3])*phi[1];
    eta[1] += conj_PRECISION(D[4])*phi[1];
    eta[2] += conj_PRECISION(D[5])*phi[1];
    eta[0] += conj_PRECISION(D[6])*phi[2];
    eta[1] += conj_PRECISION(D[7])*phi[2];
    eta[2] += conj_PRECISION(D[8])*phi[2];
  }
  
  // eta = -D*phi
  static inline void nmvm_PRECISION( const vector_PRECISION eta, const complex_PRECISION *D, const vector_PRECISION phi ) {
    eta[0]  = - D[0]*phi[0];
    eta[0] -=   D[1]*phi[1];
    eta[0] -=   D[2]*phi[2];
    eta[1]  = - D[3]*phi[0];
    eta[1] -=   D[4]*phi[1];
    eta[1] -=   D[5]*phi[2];
    eta[2]  = - D[6]*phi[0];
    eta[2] -=   D[7]*phi[1];
    eta[2] -=   D[8]*phi[2];
  }

  // eta = -D**H*phi
  static inline void nmvmh_PRECISION( const vector_PRECISION eta, const complex_PRECISION *D, const vector_PRECISION phi ) {
    eta[0]  = - conj_PRECISION(D[0])*phi[0];
    eta[1]  = - conj_PRECISION(D[1])*phi[0];
    eta[2]  = - conj_PRECISION(D[2])*phi[0];
    eta[0] -=   conj_PRECISION(D[3])*phi[1];
    eta[1] -=   conj_PRECISION(D[4])*phi[1];
    eta[2] -=   conj_PRECISION(D[5])*phi[1];
    eta[0] -=   conj_PRECISION(D[6])*phi[2];
    eta[1] -=   conj_PRECISION(D[7])*phi[2];
    eta[2] -=   conj_PRECISION(D[8])*phi[2];
  }

  // 1 - gamma_T
  static inline void prp_T_PRECISION( const vector_PRECISION prp_pt, const vector_PRECISION l_pt ) {
    prp_pt[0] = l_pt[0] -GAMMA_T_SPIN0_VAL*l_pt[3*GAMMA_T_SPIN0_CO];
    prp_pt[1] = l_pt[1] -GAMMA_T_SPIN0_VAL*l_pt[3*GAMMA_T_SPIN0_CO+1];
    prp_pt[2] = l_pt[2] -GAMMA_T_SPIN0_VAL*l_pt[3*GAMMA_T_SPIN0_CO+2];
    prp_pt[3] = l_pt[3] -GAMMA_T_SPIN1_VAL*l_pt[3*GAMMA_T_SPIN1_CO];
    prp_pt[4] = l_pt[4] -GAMMA_T_SPIN1_VAL*l_pt[3*GAMMA_T_SPIN1_CO+1];
    prp_pt[5] = l_pt[5] -GAMMA_T_SPIN1_VAL*l_pt[3*GAMMA_T_SPIN1_CO+2];
  }

  // 1 + gamma_T
  static inline void prn_T_PRECISION( const vector_PRECISION prn_pt, const vector_PRECISION l_pt ) {
    prn_pt[0] = l_pt[0] +GAMMA_T_SPIN0_VAL*l_pt[3*GAMMA_T_SPIN0_CO];
    prn_pt[1] = l_pt[1] +GAMMA_T_SPIN0_VAL*l_pt[3*GAMMA_T_SPIN0_CO+1];
    prn_pt[2] = l_pt[2] +GAMMA_T_SPIN0_VAL*l_pt[3*GAMMA_T_SPIN0_CO+2];
    prn_pt[3] = l_pt[3] +GAMMA_T_SPIN1_VAL*l_pt[3*GAMMA_T_SPIN1_CO];
    prn_pt[4] = l_pt[4] +GAMMA_T_SPIN1_VAL*l_pt[3*GAMMA_T_SPIN1_CO+1];
    prn_pt[5] = l_pt[5] +GAMMA_T_SPIN1_VAL*l_pt[3*GAMMA_T_SPIN1_CO+2];
  }

  // - (1 - gamma_T)
  static inline void pbp_su3_T_PRECISION( const vector_PRECISION prp_su3_pt, const vector_PRECISION l_pt ) {
    l_pt[ 0] -= prp_su3_pt[0];
    l_pt[ 1] -= prp_su3_pt[1];
    l_pt[ 2] -= prp_su3_pt[2];
    l_pt[ 3] -= prp_su3_pt[3];
    l_pt[ 4] -= prp_su3_pt[4];
    l_pt[ 5] -= prp_su3_pt[5];
    l_pt[ 6] += GAMMA_T_SPIN2_VAL*prp_su3_pt[3*GAMMA_T_SPIN2_CO];
    l_pt[ 7] += GAMMA_T_SPIN2_VAL*prp_su3_pt[3*GAMMA_T_SPIN2_CO+1];
    l_pt[ 8] += GAMMA_T_SPIN2_VAL*prp_su3_pt[3*GAMMA_T_SPIN2_CO+2];
    l_pt[ 9] += GAMMA_T_SPIN3_VAL*prp_su3_pt[3*GAMMA_T_SPIN3_CO];
    l_pt[10] += GAMMA_T_SPIN3_VAL*prp_su3_pt[3*GAMMA_T_SPIN3_CO+1];
    l_pt[11] += GAMMA_T_SPIN3_VAL*prp_su3_pt[3*GAMMA_T_SPIN3_CO+2];
  }

  // -(1 + gamma_T)
  static inline void pbn_su3_T_PRECISION( const vector_PRECISION prn_su3_pt, const vector_PRECISION l_pt ) {
    l_pt[ 0] -= prn_su3_pt[0];
    l_pt[ 1] -= prn_su3_pt[1];
    l_pt[ 2] -= prn_su3_pt[2];
    l_pt[ 3] -= prn_su3_pt[3];
    l_pt[ 4] -= prn_su3_pt[4];
    l_pt[ 5] -= prn_su3_pt[5];
    l_pt[ 6] -= GAMMA_T_SPIN2_VAL*prn_su3_pt[3*GAMMA_T_SPIN2_CO];
    l_pt[ 7] -= GAMMA_T_SPIN2_VAL*prn_su3_pt[3*GAMMA_T_SPIN2_CO+1];
    l_pt[ 8] -= GAMMA_T_SPIN2_VAL*prn_su3_pt[3*GAMMA_T_SPIN2_CO+2];
    l_pt[ 9] -= GAMMA_T_SPIN3_VAL*prn_su3_pt[3*GAMMA_T_SPIN3_CO];
    l_pt[10] -= GAMMA_T_SPIN3_VAL*prn_su3_pt[3*GAMMA_T_SPIN3_CO+1];
    l_pt[11] -= GAMMA_T_SPIN3_VAL*prn_su3_pt[3*GAMMA_T_SPIN3_CO+2];
  }

  static inline void prp_Z_PRECISION( const vector_PRECISION prp_pt, const vector_PRECISION l_pt ) {
    prp_pt[0] = l_pt[0] -GAMMA_Z_SPIN0_VAL*l_pt[3*GAMMA_Z_SPIN0_CO];
    prp_pt[1] = l_pt[1] -GAMMA_Z_SPIN0_VAL*l_pt[3*GAMMA_Z_SPIN0_CO+1];
    prp_pt[2] = l_pt[2] -GAMMA_Z_SPIN0_VAL*l_pt[3*GAMMA_Z_SPIN0_CO+2];
    prp_pt[3] = l_pt[3] -GAMMA_Z_SPIN1_VAL*l_pt[3*GAMMA_Z_SPIN1_CO];
    prp_pt[4] = l_pt[4] -GAMMA_Z_SPIN1_VAL*l_pt[3*GAMMA_Z_SPIN1_CO+1];
    prp_pt[5] = l_pt[5] -GAMMA_Z_SPIN1_VAL*l_pt[3*GAMMA_Z_SPIN1_CO+2];
  }

  static inline void prn_Z_PRECISION( const vector_PRECISION prn_pt, const vector_PRECISION l_pt ) {
    prn_pt[0] = l_pt[0] +GAMMA_Z_SPIN0_VAL*l_pt[3*GAMMA_Z_SPIN0_CO];
    prn_pt[1] = l_pt[1] +GAMMA_Z_SPIN0_VAL*l_pt[3*GAMMA_Z_SPIN0_CO+1];
    prn_pt[2] = l_pt[2] +GAMMA_Z_SPIN0_VAL*l_pt[3*GAMMA_Z_SPIN0_CO+2];
    prn_pt[3] = l_pt[3] +GAMMA_Z_SPIN1_VAL*l_pt[3*GAMMA_Z_SPIN1_CO];
    prn_pt[4] = l_pt[4] +GAMMA_Z_SPIN1_VAL*l_pt[3*GAMMA_Z_SPIN1_CO+1];
    prn_pt[5] = l_pt[5] +GAMMA_Z_SPIN1_VAL*l_pt[3*GAMMA_Z_SPIN1_CO+2];
  }

  static inline void pbp_su3_Z_PRECISION( const vector_PRECISION prp_su3_pt, const vector_PRECISION l_pt ) {
    l_pt[ 0] -= prp_su3_pt[0];
    l_pt[ 1] -= prp_su3_pt[1];
    l_pt[ 2] -= prp_su3_pt[2];
    l_pt[ 3] -= prp_su3_pt[3];
    l_pt[ 4] -= prp_su3_pt[4];
    l_pt[ 5] -= prp_su3_pt[5];
    l_pt[ 6] += GAMMA_Z_SPIN2_VAL*prp_su3_pt[3*GAMMA_Z_SPIN2_CO];
    l_pt[ 7] += GAMMA_Z_SPIN2_VAL*prp_su3_pt[3*GAMMA_Z_SPIN2_CO+1];
    l_pt[ 8] += GAMMA_Z_SPIN2_VAL*prp_su3_pt[3*GAMMA_Z_SPIN2_CO+2];
    l_pt[ 9] += GAMMA_Z_SPIN3_VAL*prp_su3_pt[3*GAMMA_Z_SPIN3_CO];
    l_pt[10] += GAMMA_Z_SPIN3_VAL*prp_su3_pt[3*GAMMA_Z_SPIN3_CO+1];
    l_pt[11] += GAMMA_Z_SPIN3_VAL*prp_su3_pt[3*GAMMA_Z_SPIN3_CO+2];
  }

  static inline void pbn_su3_Z_PRECISION( const vector_PRECISION prn_su3_pt, const vector_PRECISION l_pt ) {
    l_pt[ 0] -= prn_su3_pt[0];
    l_pt[ 1] -= prn_su3_pt[1];
    l_pt[ 2] -= prn_su3_pt[2];
    l_pt[ 3] -= prn_su3_pt[3];
    l_pt[ 4] -= prn_su3_pt[4];
    l_pt[ 5] -= prn_su3_pt[5];
    l_pt[ 6] -= GAMMA_Z_SPIN2_VAL*prn_su3_pt[3*GAMMA_Z_SPIN2_CO];
    l_pt[ 7] -= GAMMA_Z_SPIN2_VAL*prn_su3_pt[3*GAMMA_Z_SPIN2_CO+1];
    l_pt[ 8] -= GAMMA_Z_SPIN2_VAL*prn_su3_pt[3*GAMMA_Z_SPIN2_CO+2];
    l_pt[ 9] -= GAMMA_Z_SPIN3_VAL*prn_su3_pt[3*GAMMA_Z_SPIN3_CO];
    l_pt[10] -= GAMMA_Z_SPIN3_VAL*prn_su3_pt[3*GAMMA_Z_SPIN3_CO+1];
    l_pt[11] -= GAMMA_Z_SPIN3_VAL*prn_su3_pt[3*GAMMA_Z_SPIN3_CO+2];
  }

  static inline void prp_Y_PRECISION( const vector_PRECISION prp_pt, const vector_PRECISION l_pt ) {
    prp_pt[0] = l_pt[0] -GAMMA_Y_SPIN0_VAL*l_pt[3*GAMMA_Y_SPIN0_CO];
    prp_pt[1] = l_pt[1] -GAMMA_Y_SPIN0_VAL*l_pt[3*GAMMA_Y_SPIN0_CO+1];
    prp_pt[2] = l_pt[2] -GAMMA_Y_SPIN0_VAL*l_pt[3*GAMMA_Y_SPIN0_CO+2];
    prp_pt[3] = l_pt[3] -GAMMA_Y_SPIN1_VAL*l_pt[3*GAMMA_Y_SPIN1_CO];
    prp_pt[4] = l_pt[4] -GAMMA_Y_SPIN1_VAL*l_pt[3*GAMMA_Y_SPIN1_CO+1];
    prp_pt[5] = l_pt[5] -GAMMA_Y_SPIN1_VAL*l_pt[3*GAMMA_Y_SPIN1_CO+2];
  }

  static inline void prn_Y_PRECISION( const vector_PRECISION prn_pt, const vector_PRECISION l_pt ) {
    prn_pt[0] = l_pt[0] +GAMMA_Y_SPIN0_VAL*l_pt[3*GAMMA_Y_SPIN0_CO];
    prn_pt[1] = l_pt[1] +GAMMA_Y_SPIN0_VAL*l_pt[3*GAMMA_Y_SPIN0_CO+1];
    prn_pt[2] = l_pt[2] +GAMMA_Y_SPIN0_VAL*l_pt[3*GAMMA_Y_SPIN0_CO+2];
    prn_pt[3] = l_pt[3] +GAMMA_Y_SPIN1_VAL*l_pt[3*GAMMA_Y_SPIN1_CO];
    prn_pt[4] = l_pt[4] +GAMMA_Y_SPIN1_VAL*l_pt[3*GAMMA_Y_SPIN1_CO+1];
    prn_pt[5] = l_pt[5] +GAMMA_Y_SPIN1_VAL*l_pt[3*GAMMA_Y_SPIN1_CO+2];
  }

  static inline void pbp_su3_Y_PRECISION( const vector_PRECISION prp_su3_pt, const vector_PRECISION l_pt ) {
    l_pt[ 0] -= prp_su3_pt[0];
    l_pt[ 1] -= prp_su3_pt[1];
    l_pt[ 2] -= prp_su3_pt[2];
    l_pt[ 3] -= prp_su3_pt[3];
    l_pt[ 4] -= prp_su3_pt[4];
    l_pt[ 5] -= prp_su3_pt[5];
    l_pt[ 6] += GAMMA_Y_SPIN2_VAL*prp_su3_pt[3*GAMMA_Y_SPIN2_CO];
    l_pt[ 7] += GAMMA_Y_SPIN2_VAL*prp_su3_pt[3*GAMMA_Y_SPIN2_CO+1];
    l_pt[ 8] += GAMMA_Y_SPIN2_VAL*prp_su3_pt[3*GAMMA_Y_SPIN2_CO+2];
    l_pt[ 9] += GAMMA_Y_SPIN3_VAL*prp_su3_pt[3*GAMMA_Y_SPIN3_CO];
    l_pt[10] += GAMMA_Y_SPIN3_VAL*prp_su3_pt[3*GAMMA_Y_SPIN3_CO+1];
    l_pt[11] += GAMMA_Y_SPIN3_VAL*prp_su3_pt[3*GAMMA_Y_SPIN3_CO+2];
  }

  static inline void pbn_su3_Y_PRECISION( const vector_PRECISION prn_su3_pt, const vector_PRECISION l_pt ) {
    l_pt[ 0] -= prn_su3_pt[0];
    l_pt[ 1] -= prn_su3_pt[1];
    l_pt[ 2] -= prn_su3_pt[2];
    l_pt[ 3] -= prn_su3_pt[3];
    l_pt[ 4] -= prn_su3_pt[4];
    l_pt[ 5] -= prn_su3_pt[5];
    l_pt[ 6] -= GAMMA_Y_SPIN2_VAL*prn_su3_pt[3*GAMMA_Y_SPIN2_CO];
    l_pt[ 7] -= GAMMA_Y_SPIN2_VAL*prn_su3_pt[3*GAMMA_Y_SPIN2_CO+1];
    l_pt[ 8] -= GAMMA_Y_SPIN2_VAL*prn_su3_pt[3*GAMMA_Y_SPIN2_CO+2];
    l_pt[ 9] -= GAMMA_Y_SPIN3_VAL*prn_su3_pt[3*GAMMA_Y_SPIN3_CO];
    l_pt[10] -= GAMMA_Y_SPIN3_VAL*prn_su3_pt[3*GAMMA_Y_SPIN3_CO+1];
    l_pt[11] -= GAMMA_Y_SPIN3_VAL*prn_su3_pt[3*GAMMA_Y_SPIN3_CO+2];
  }

  static inline void prp_X_PRECISION( const vector_PRECISION prp_pt, const vector_PRECISION l_pt ) {
    prp_pt[0] = l_pt[0] -GAMMA_X_SPIN0_VAL*l_pt[3*GAMMA_X_SPIN0_CO];
    prp_pt[1] = l_pt[1] -GAMMA_X_SPIN0_VAL*l_pt[3*GAMMA_X_SPIN0_CO+1];
    prp_pt[2] = l_pt[2] -GAMMA_X_SPIN0_VAL*l_pt[3*GAMMA_X_SPIN0_CO+2];
    prp_pt[3] = l_pt[3] -GAMMA_X_SPIN1_VAL*l_pt[3*GAMMA_X_SPIN1_CO];
    prp_pt[4] = l_pt[4] -GAMMA_X_SPIN1_VAL*l_pt[3*GAMMA_X_SPIN1_CO+1];
    prp_pt[5] = l_pt[5] -GAMMA_X_SPIN1_VAL*l_pt[3*GAMMA_X_SPIN1_CO+2];
  }

  static inline void prn_X_PRECISION( const vector_PRECISION prn_pt, const vector_PRECISION l_pt ) {
    prn_pt[0] = l_pt[0] +GAMMA_X_SPIN0_VAL*l_pt[3*GAMMA_X_SPIN0_CO];
    prn_pt[1] = l_pt[1] +GAMMA_X_SPIN0_VAL*l_pt[3*GAMMA_X_SPIN0_CO+1];
    prn_pt[2] = l_pt[2] +GAMMA_X_SPIN0_VAL*l_pt[3*GAMMA_X_SPIN0_CO+2];
    prn_pt[3] = l_pt[3] +GAMMA_X_SPIN1_VAL*l_pt[3*GAMMA_X_SPIN1_CO];
    prn_pt[4] = l_pt[4] +GAMMA_X_SPIN1_VAL*l_pt[3*GAMMA_X_SPIN1_CO+1];
    prn_pt[5] = l_pt[5] +GAMMA_X_SPIN1_VAL*l_pt[3*GAMMA_X_SPIN1_CO+2];
  }

  static inline void pbp_su3_X_PRECISION( const vector_PRECISION prp_su3_pt, const vector_PRECISION l_pt ) {
    l_pt[ 0] -= prp_su3_pt[0];
    l_pt[ 1] -= prp_su3_pt[1];
    l_pt[ 2] -= prp_su3_pt[2];
    l_pt[ 3] -= prp_su3_pt[3];
    l_pt[ 4] -= prp_su3_pt[4];
    l_pt[ 5] -= prp_su3_pt[5];
    l_pt[ 6] += GAMMA_X_SPIN2_VAL*prp_su3_pt[3*GAMMA_X_SPIN2_CO];
    l_pt[ 7] += GAMMA_X_SPIN2_VAL*prp_su3_pt[3*GAMMA_X_SPIN2_CO+1];
    l_pt[ 8] += GAMMA_X_SPIN2_VAL*prp_su3_pt[3*GAMMA_X_SPIN2_CO+2];
    l_pt[ 9] += GAMMA_X_SPIN3_VAL*prp_su3_pt[3*GAMMA_X_SPIN3_CO];
    l_pt[10] += GAMMA_X_SPIN3_VAL*prp_su3_pt[3*GAMMA_X_SPIN3_CO+1];
    l_pt[11] += GAMMA_X_SPIN3_VAL*prp_su3_pt[3*GAMMA_X_SPIN3_CO+2];
  }

  static inline void pbn_su3_X_PRECISION( const vector_PRECISION prn_su3_pt, const vector_PRECISION l_pt ) {
    l_pt[ 0] -= prn_su3_pt[0];
    l_pt[ 1] -= prn_su3_pt[1];
    l_pt[ 2] -= prn_su3_pt[2];
    l_pt[ 3] -= prn_su3_pt[3];
    l_pt[ 4] -= prn_su3_pt[4];
    l_pt[ 5] -= prn_su3_pt[5];
    l_pt[ 6] -= GAMMA_X_SPIN2_VAL*prn_su3_pt[3*GAMMA_X_SPIN2_CO];
    l_pt[ 7] -= GAMMA_X_SPIN2_VAL*prn_su3_pt[3*GAMMA_X_SPIN2_CO+1];
    l_pt[ 8] -= GAMMA_X_SPIN2_VAL*prn_su3_pt[3*GAMMA_X_SPIN2_CO+2];
    l_pt[ 9] -= GAMMA_X_SPIN3_VAL*prn_su3_pt[3*GAMMA_X_SPIN3_CO];
    l_pt[10] -= GAMMA_X_SPIN3_VAL*prn_su3_pt[3*GAMMA_X_SPIN3_CO+1];
    l_pt[11] -= GAMMA_X_SPIN3_VAL*prn_su3_pt[3*GAMMA_X_SPIN3_CO+2];
  }

  static inline void twospin_p_T_PRECISION( const vector_PRECISION out_spin0and1, const vector_PRECISION out_spin2and3, const vector_PRECISION in ) {
    out_spin0and1[ 0] -= in[ 0];
    out_spin0and1[ 1] -= in[ 1];
    out_spin0and1[ 2] -= in[ 2];
    out_spin0and1[ 3] -= in[ 3];
    out_spin0and1[ 4] -= in[ 4];
    out_spin0and1[ 5] -= in[ 5];
    out_spin0and1[ 6] += GAMMA_T_SPIN2_VAL * in[3*GAMMA_T_SPIN2_CO  ];
    out_spin0and1[ 7] += GAMMA_T_SPIN2_VAL * in[3*GAMMA_T_SPIN2_CO+1];
    out_spin0and1[ 8] += GAMMA_T_SPIN2_VAL * in[3*GAMMA_T_SPIN2_CO+2];
    out_spin0and1[ 9] += GAMMA_T_SPIN3_VAL * in[3*GAMMA_T_SPIN3_CO  ];
    out_spin0and1[10] += GAMMA_T_SPIN3_VAL * in[3*GAMMA_T_SPIN3_CO+1];
    out_spin0and1[11] += GAMMA_T_SPIN3_VAL * in[3*GAMMA_T_SPIN3_CO+2];
    out_spin2and3[ 0] += GAMMA_T_SPIN0_VAL * in[3*GAMMA_T_SPIN0_CO  ];
    out_spin2and3[ 1] += GAMMA_T_SPIN0_VAL * in[3*GAMMA_T_SPIN0_CO+1];
    out_spin2and3[ 2] += GAMMA_T_SPIN0_VAL * in[3*GAMMA_T_SPIN0_CO+2];
    out_spin2and3[ 3] += GAMMA_T_SPIN1_VAL * in[3*GAMMA_T_SPIN1_CO  ];
    out_spin2and3[ 4] += GAMMA_T_SPIN1_VAL * in[3*GAMMA_T_SPIN1_CO+1];
    out_spin2and3[ 5] += GAMMA_T_SPIN1_VAL * in[3*GAMMA_T_SPIN1_CO+2];
    out_spin2and3[ 6] -= in[ 6];
    out_spin2and3[ 7] -= in[ 7];
    out_spin2and3[ 8] -= in[ 8];
    out_spin2and3[ 9] -= in[ 9];
    out_spin2and3[10] -= in[10];
    out_spin2and3[11] -= in[11];
  }

  static inline void twospin2_p_T_PRECISION( const vector_PRECISION out_spin0and1, const vector_PRECISION out_spin2and3, const vector_PRECISION in ) {
    out_spin0and1[ 0] = in[ 0];
    out_spin0and1[ 1] = in[ 1];
    out_spin0and1[ 2] = in[ 2];
    out_spin0and1[ 3] = in[ 3];
    out_spin0and1[ 4] = in[ 4];
    out_spin0and1[ 5] = in[ 5];
    out_spin0and1[ 6] = -GAMMA_T_SPIN2_VAL * in[3*GAMMA_T_SPIN2_CO  ];
    out_spin0and1[ 7] = -GAMMA_T_SPIN2_VAL * in[3*GAMMA_T_SPIN2_CO+1];
    out_spin0and1[ 8] = -GAMMA_T_SPIN2_VAL * in[3*GAMMA_T_SPIN2_CO+2];
    out_spin0and1[ 9] = -GAMMA_T_SPIN3_VAL * in[3*GAMMA_T_SPIN3_CO  ];
    out_spin0and1[10] = -GAMMA_T_SPIN3_VAL * in[3*GAMMA_T_SPIN3_CO+1];
    out_spin0and1[11] = -GAMMA_T_SPIN3_VAL * in[3*GAMMA_T_SPIN3_CO+2];
    out_spin2and3[ 0] = -GAMMA_T_SPIN0_VAL * in[3*GAMMA_T_SPIN0_CO  ];
    out_spin2and3[ 1] = -GAMMA_T_SPIN0_VAL * in[3*GAMMA_T_SPIN0_CO+1];
    out_spin2and3[ 2] = -GAMMA_T_SPIN0_VAL * in[3*GAMMA_T_SPIN0_CO+2];
    out_spin2and3[ 3] = -GAMMA_T_SPIN1_VAL * in[3*GAMMA_T_SPIN1_CO  ];
    out_spin2and3[ 4] = -GAMMA_T_SPIN1_VAL * in[3*GAMMA_T_SPIN1_CO+1];
    out_spin2and3[ 5] = -GAMMA_T_SPIN1_VAL * in[3*GAMMA_T_SPIN1_CO+2];
    out_spin2and3[ 6] = in[ 6];
    out_spin2and3[ 7] = in[ 7];
    out_spin2and3[ 8] = in[ 8];
    out_spin2and3[ 9] = in[ 9];
    out_spin2and3[10] = in[10];
    out_spin2and3[11] = in[11];
  }

  static inline void twospin_n_T_PRECISION( const vector_PRECISION out_spin0and1, const vector_PRECISION out_spin2and3, const vector_PRECISION in ) {
    out_spin0and1[ 0] -= in[ 0];
    out_spin0and1[ 1] -= in[ 1];
    out_spin0and1[ 2] -= in[ 2];
    out_spin0and1[ 3] -= in[ 3];
    out_spin0and1[ 4] -= in[ 4];
    out_spin0and1[ 5] -= in[ 5];
    out_spin0and1[ 6] -= GAMMA_T_SPIN2_VAL * in[3*GAMMA_T_SPIN2_CO  ];
    out_spin0and1[ 7] -= GAMMA_T_SPIN2_VAL * in[3*GAMMA_T_SPIN2_CO+1];
    out_spin0and1[ 8] -= GAMMA_T_SPIN2_VAL * in[3*GAMMA_T_SPIN2_CO+2];
    out_spin0and1[ 9] -= GAMMA_T_SPIN3_VAL * in[3*GAMMA_T_SPIN3_CO  ];
    out_spin0and1[10] -= GAMMA_T_SPIN3_VAL * in[3*GAMMA_T_SPIN3_CO+1];
    out_spin0and1[11] -= GAMMA_T_SPIN3_VAL * in[3*GAMMA_T_SPIN3_CO+2];
    out_spin2and3[ 0] -= GAMMA_T_SPIN0_VAL * in[3*GAMMA_T_SPIN0_CO  ];
    out_spin2and3[ 1] -= GAMMA_T_SPIN0_VAL * in[3*GAMMA_T_SPIN0_CO+1];
    out_spin2and3[ 2] -= GAMMA_T_SPIN0_VAL * in[3*GAMMA_T_SPIN0_CO+2];
    out_spin2and3[ 3] -= GAMMA_T_SPIN1_VAL * in[3*GAMMA_T_SPIN1_CO  ];
    out_spin2and3[ 4] -= GAMMA_T_SPIN1_VAL * in[3*GAMMA_T_SPIN1_CO+1];
    out_spin2and3[ 5] -= GAMMA_T_SPIN1_VAL * in[3*GAMMA_T_SPIN1_CO+2];
    out_spin2and3[ 6] -= in[ 6];
    out_spin2and3[ 7] -= in[ 7];
    out_spin2and3[ 8] -= in[ 8];
    out_spin2and3[ 9] -= in[ 9];
    out_spin2and3[10] -= in[10];
    out_spin2and3[11] -= in[11];
  }

  static inline void twospin_p_Z_PRECISION( const vector_PRECISION out_spin0and1, const vector_PRECISION out_spin2and3, const vector_PRECISION in ) {
    out_spin0and1[ 0] -= in[ 0];
    out_spin0and1[ 1] -= in[ 1];
    out_spin0and1[ 2] -= in[ 2];
    out_spin0and1[ 3] -= in[ 3];
    out_spin0and1[ 4] -= in[ 4];
    out_spin0and1[ 5] -= in[ 5];
    out_spin0and1[ 6] += GAMMA_Z_SPIN2_VAL * in[3*GAMMA_Z_SPIN2_CO  ];
    out_spin0and1[ 7] += GAMMA_Z_SPIN2_VAL * in[3*GAMMA_Z_SPIN2_CO+1];
    out_spin0and1[ 8] += GAMMA_Z_SPIN2_VAL * in[3*GAMMA_Z_SPIN2_CO+2];
    out_spin0and1[ 9] += GAMMA_Z_SPIN3_VAL * in[3*GAMMA_Z_SPIN3_CO  ];
    out_spin0and1[10] += GAMMA_Z_SPIN3_VAL * in[3*GAMMA_Z_SPIN3_CO+1];
    out_spin0and1[11] += GAMMA_Z_SPIN3_VAL * in[3*GAMMA_Z_SPIN3_CO+2];
    out_spin2and3[ 0] += GAMMA_Z_SPIN0_VAL * in[3*GAMMA_Z_SPIN0_CO  ];
    out_spin2and3[ 1] += GAMMA_Z_SPIN0_VAL * in[3*GAMMA_Z_SPIN0_CO+1];
    out_spin2and3[ 2] += GAMMA_Z_SPIN0_VAL * in[3*GAMMA_Z_SPIN0_CO+2];
    out_spin2and3[ 3] += GAMMA_Z_SPIN1_VAL * in[3*GAMMA_Z_SPIN1_CO  ];
    out_spin2and3[ 4] += GAMMA_Z_SPIN1_VAL * in[3*GAMMA_Z_SPIN1_CO+1];
    out_spin2and3[ 5] += GAMMA_Z_SPIN1_VAL * in[3*GAMMA_Z_SPIN1_CO+2];
    out_spin2and3[ 6] -= in[ 6];
    out_spin2and3[ 7] -= in[ 7];
    out_spin2and3[ 8] -= in[ 8];
    out_spin2and3[ 9] -= in[ 9];
    out_spin2and3[10] -= in[10];
    out_spin2and3[11] -= in[11];
  }

  static inline void twospin2_p_Z_PRECISION( const vector_PRECISION out_spin0and1, const vector_PRECISION out_spin2and3, const vector_PRECISION in ) {
    out_spin0and1[ 0] = in[ 0];
    out_spin0and1[ 1] = in[ 1];
    out_spin0and1[ 2] = in[ 2];
    out_spin0and1[ 3] = in[ 3];
    out_spin0and1[ 4] = in[ 4];
    out_spin0and1[ 5] = in[ 5];
    out_spin0and1[ 6] = -GAMMA_Z_SPIN2_VAL * in[3*GAMMA_Z_SPIN2_CO  ];
    out_spin0and1[ 7] = -GAMMA_Z_SPIN2_VAL * in[3*GAMMA_Z_SPIN2_CO+1];
    out_spin0and1[ 8] = -GAMMA_Z_SPIN2_VAL * in[3*GAMMA_Z_SPIN2_CO+2];
    out_spin0and1[ 9] = -GAMMA_Z_SPIN3_VAL * in[3*GAMMA_Z_SPIN3_CO  ];
    out_spin0and1[10] = -GAMMA_Z_SPIN3_VAL * in[3*GAMMA_Z_SPIN3_CO+1];
    out_spin0and1[11] = -GAMMA_Z_SPIN3_VAL * in[3*GAMMA_Z_SPIN3_CO+2];
    out_spin2and3[ 0] = -GAMMA_Z_SPIN0_VAL * in[3*GAMMA_Z_SPIN0_CO  ];
    out_spin2and3[ 1] = -GAMMA_Z_SPIN0_VAL * in[3*GAMMA_Z_SPIN0_CO+1];
    out_spin2and3[ 2] = -GAMMA_Z_SPIN0_VAL * in[3*GAMMA_Z_SPIN0_CO+2];
    out_spin2and3[ 3] = -GAMMA_Z_SPIN1_VAL * in[3*GAMMA_Z_SPIN1_CO  ];
    out_spin2and3[ 4] = -GAMMA_Z_SPIN1_VAL * in[3*GAMMA_Z_SPIN1_CO+1];
    out_spin2and3[ 5] = -GAMMA_Z_SPIN1_VAL * in[3*GAMMA_Z_SPIN1_CO+2];
    out_spin2and3[ 6] = in[ 6];
    out_spin2and3[ 7] = in[ 7];
    out_spin2and3[ 8] = in[ 8];
    out_spin2and3[ 9] = in[ 9];
    out_spin2and3[10] = in[10];
    out_spin2and3[11] = in[11];
  }

  static inline void twospin_n_Z_PRECISION( const vector_PRECISION out_spin0and1, const vector_PRECISION out_spin2and3, const vector_PRECISION in ) {
    out_spin0and1[ 0] -= in[ 0];
    out_spin0and1[ 1] -= in[ 1];
    out_spin0and1[ 2] -= in[ 2];
    out_spin0and1[ 3] -= in[ 3];
    out_spin0and1[ 4] -= in[ 4];
    out_spin0and1[ 5] -= in[ 5];
    out_spin0and1[ 6] -= GAMMA_Z_SPIN2_VAL * in[3*GAMMA_Z_SPIN2_CO  ];
    out_spin0and1[ 7] -= GAMMA_Z_SPIN2_VAL * in[3*GAMMA_Z_SPIN2_CO+1];
    out_spin0and1[ 8] -= GAMMA_Z_SPIN2_VAL * in[3*GAMMA_Z_SPIN2_CO+2];
    out_spin0and1[ 9] -= GAMMA_Z_SPIN3_VAL * in[3*GAMMA_Z_SPIN3_CO  ];
    out_spin0and1[10] -= GAMMA_Z_SPIN3_VAL * in[3*GAMMA_Z_SPIN3_CO+1];
    out_spin0and1[11] -= GAMMA_Z_SPIN3_VAL * in[3*GAMMA_Z_SPIN3_CO+2];
    out_spin2and3[ 0] -= GAMMA_Z_SPIN0_VAL * in[3*GAMMA_Z_SPIN0_CO  ];
    out_spin2and3[ 1] -= GAMMA_Z_SPIN0_VAL * in[3*GAMMA_Z_SPIN0_CO+1];
    out_spin2and3[ 2] -= GAMMA_Z_SPIN0_VAL * in[3*GAMMA_Z_SPIN0_CO+2];
    out_spin2and3[ 3] -= GAMMA_Z_SPIN1_VAL * in[3*GAMMA_Z_SPIN1_CO  ];
    out_spin2and3[ 4] -= GAMMA_Z_SPIN1_VAL * in[3*GAMMA_Z_SPIN1_CO+1];
    out_spin2and3[ 5] -= GAMMA_Z_SPIN1_VAL * in[3*GAMMA_Z_SPIN1_CO+2];
    out_spin2and3[ 6] -= in[ 6];
    out_spin2and3[ 7] -= in[ 7];
    out_spin2and3[ 8] -= in[ 8];
    out_spin2and3[ 9] -= in[ 9];
    out_spin2and3[10] -= in[10];
    out_spin2and3[11] -= in[11];
  }

  static inline void twospin_p_Y_PRECISION( const vector_PRECISION out_spin0and1, const vector_PRECISION out_spin2and3, const vector_PRECISION in ) {
    out_spin0and1[ 0] -= in[ 0];
    out_spin0and1[ 1] -= in[ 1];
    out_spin0and1[ 2] -= in[ 2];
    out_spin0and1[ 3] -= in[ 3];
    out_spin0and1[ 4] -= in[ 4];
    out_spin0and1[ 5] -= in[ 5];
    out_spin0and1[ 6] += GAMMA_Y_SPIN2_VAL * in[3*GAMMA_Y_SPIN2_CO  ];
    out_spin0and1[ 7] += GAMMA_Y_SPIN2_VAL * in[3*GAMMA_Y_SPIN2_CO+1];
    out_spin0and1[ 8] += GAMMA_Y_SPIN2_VAL * in[3*GAMMA_Y_SPIN2_CO+2];
    out_spin0and1[ 9] += GAMMA_Y_SPIN3_VAL * in[3*GAMMA_Y_SPIN3_CO  ];
    out_spin0and1[10] += GAMMA_Y_SPIN3_VAL * in[3*GAMMA_Y_SPIN3_CO+1];
    out_spin0and1[11] += GAMMA_Y_SPIN3_VAL * in[3*GAMMA_Y_SPIN3_CO+2];
    out_spin2and3[ 0] += GAMMA_Y_SPIN0_VAL * in[3*GAMMA_Y_SPIN0_CO  ];
    out_spin2and3[ 1] += GAMMA_Y_SPIN0_VAL * in[3*GAMMA_Y_SPIN0_CO+1];
    out_spin2and3[ 2] += GAMMA_Y_SPIN0_VAL * in[3*GAMMA_Y_SPIN0_CO+2];
    out_spin2and3[ 3] += GAMMA_Y_SPIN1_VAL * in[3*GAMMA_Y_SPIN1_CO  ];
    out_spin2and3[ 4] += GAMMA_Y_SPIN1_VAL * in[3*GAMMA_Y_SPIN1_CO+1];
    out_spin2and3[ 5] += GAMMA_Y_SPIN1_VAL * in[3*GAMMA_Y_SPIN1_CO+2];
    out_spin2and3[ 6] -= in[ 6];
    out_spin2and3[ 7] -= in[ 7];
    out_spin2and3[ 8] -= in[ 8];
    out_spin2and3[ 9] -= in[ 9];
    out_spin2and3[10] -= in[10];
    out_spin2and3[11] -= in[11];
  }

  static inline void twospin2_p_Y_PRECISION( const vector_PRECISION out_spin0and1, const vector_PRECISION out_spin2and3, const vector_PRECISION in ) {
    out_spin0and1[ 0] = in[ 0];
    out_spin0and1[ 1] = in[ 1];
    out_spin0and1[ 2] = in[ 2];
    out_spin0and1[ 3] = in[ 3];
    out_spin0and1[ 4] = in[ 4];
    out_spin0and1[ 5] = in[ 5];
    out_spin0and1[ 6] = -GAMMA_Y_SPIN2_VAL * in[3*GAMMA_Y_SPIN2_CO  ];
    out_spin0and1[ 7] = -GAMMA_Y_SPIN2_VAL * in[3*GAMMA_Y_SPIN2_CO+1];
    out_spin0and1[ 8] = -GAMMA_Y_SPIN2_VAL * in[3*GAMMA_Y_SPIN2_CO+2];
    out_spin0and1[ 9] = -GAMMA_Y_SPIN3_VAL * in[3*GAMMA_Y_SPIN3_CO  ];
    out_spin0and1[10] = -GAMMA_Y_SPIN3_VAL * in[3*GAMMA_Y_SPIN3_CO+1];
    out_spin0and1[11] = -GAMMA_Y_SPIN3_VAL * in[3*GAMMA_Y_SPIN3_CO+2];
    out_spin2and3[ 0] = -GAMMA_Y_SPIN0_VAL * in[3*GAMMA_Y_SPIN0_CO  ];
    out_spin2and3[ 1] = -GAMMA_Y_SPIN0_VAL * in[3*GAMMA_Y_SPIN0_CO+1];
    out_spin2and3[ 2] = -GAMMA_Y_SPIN0_VAL * in[3*GAMMA_Y_SPIN0_CO+2];
    out_spin2and3[ 3] = -GAMMA_Y_SPIN1_VAL * in[3*GAMMA_Y_SPIN1_CO  ];
    out_spin2and3[ 4] = -GAMMA_Y_SPIN1_VAL * in[3*GAMMA_Y_SPIN1_CO+1];
    out_spin2and3[ 5] = -GAMMA_Y_SPIN1_VAL * in[3*GAMMA_Y_SPIN1_CO+2];
    out_spin2and3[ 6] = in[ 6];
    out_spin2and3[ 7] = in[ 7];
    out_spin2and3[ 8] = in[ 8];
    out_spin2and3[ 9] = in[ 9];
    out_spin2and3[10] = in[10];
    out_spin2and3[11] = in[11];
  }

  static inline void twospin_n_Y_PRECISION( const vector_PRECISION out_spin0and1, const vector_PRECISION out_spin2and3, const vector_PRECISION in ) {
    out_spin0and1[ 0] -= in[ 0];
    out_spin0and1[ 1] -= in[ 1];
    out_spin0and1[ 2] -= in[ 2];
    out_spin0and1[ 3] -= in[ 3];
    out_spin0and1[ 4] -= in[ 4];
    out_spin0and1[ 5] -= in[ 5];
    out_spin0and1[ 6] -= GAMMA_Y_SPIN2_VAL * in[3*GAMMA_Y_SPIN2_CO  ];
    out_spin0and1[ 7] -= GAMMA_Y_SPIN2_VAL * in[3*GAMMA_Y_SPIN2_CO+1];
    out_spin0and1[ 8] -= GAMMA_Y_SPIN2_VAL * in[3*GAMMA_Y_SPIN2_CO+2];
    out_spin0and1[ 9] -= GAMMA_Y_SPIN3_VAL * in[3*GAMMA_Y_SPIN3_CO  ];
    out_spin0and1[10] -= GAMMA_Y_SPIN3_VAL * in[3*GAMMA_Y_SPIN3_CO+1];
    out_spin0and1[11] -= GAMMA_Y_SPIN3_VAL * in[3*GAMMA_Y_SPIN3_CO+2];
    out_spin2and3[ 0] -= GAMMA_Y_SPIN0_VAL * in[3*GAMMA_Y_SPIN0_CO  ];
    out_spin2and3[ 1] -= GAMMA_Y_SPIN0_VAL * in[3*GAMMA_Y_SPIN0_CO+1];
    out_spin2and3[ 2] -= GAMMA_Y_SPIN0_VAL * in[3*GAMMA_Y_SPIN0_CO+2];
    out_spin2and3[ 3] -= GAMMA_Y_SPIN1_VAL * in[3*GAMMA_Y_SPIN1_CO  ];
    out_spin2and3[ 4] -= GAMMA_Y_SPIN1_VAL * in[3*GAMMA_Y_SPIN1_CO+1];
    out_spin2and3[ 5] -= GAMMA_Y_SPIN1_VAL * in[3*GAMMA_Y_SPIN1_CO+2];
    out_spin2and3[ 6] -= in[ 6];
    out_spin2and3[ 7] -= in[ 7];
    out_spin2and3[ 8] -= in[ 8];
    out_spin2and3[ 9] -= in[ 9];
    out_spin2and3[10] -= in[10];
    out_spin2and3[11] -= in[11];
  }

  static inline void twospin_p_X_PRECISION( const vector_PRECISION out_spin0and1, const vector_PRECISION out_spin2and3, const vector_PRECISION in ) {
    out_spin0and1[ 0] -= in[ 0];
    out_spin0and1[ 1] -= in[ 1];
    out_spin0and1[ 2] -= in[ 2];
    out_spin0and1[ 3] -= in[ 3];
    out_spin0and1[ 4] -= in[ 4];
    out_spin0and1[ 5] -= in[ 5];
    out_spin0and1[ 6] += GAMMA_X_SPIN2_VAL * in[3*GAMMA_X_SPIN2_CO  ];
    out_spin0and1[ 7] += GAMMA_X_SPIN2_VAL * in[3*GAMMA_X_SPIN2_CO+1];
    out_spin0and1[ 8] += GAMMA_X_SPIN2_VAL * in[3*GAMMA_X_SPIN2_CO+2];
    out_spin0and1[ 9] += GAMMA_X_SPIN3_VAL * in[3*GAMMA_X_SPIN3_CO  ];
    out_spin0and1[10] += GAMMA_X_SPIN3_VAL * in[3*GAMMA_X_SPIN3_CO+1];
    out_spin0and1[11] += GAMMA_X_SPIN3_VAL * in[3*GAMMA_X_SPIN3_CO+2];
    out_spin2and3[ 0] += GAMMA_X_SPIN0_VAL * in[3*GAMMA_X_SPIN0_CO  ];
    out_spin2and3[ 1] += GAMMA_X_SPIN0_VAL * in[3*GAMMA_X_SPIN0_CO+1];
    out_spin2and3[ 2] += GAMMA_X_SPIN0_VAL * in[3*GAMMA_X_SPIN0_CO+2];
    out_spin2and3[ 3] += GAMMA_X_SPIN1_VAL * in[3*GAMMA_X_SPIN1_CO  ];
    out_spin2and3[ 4] += GAMMA_X_SPIN1_VAL * in[3*GAMMA_X_SPIN1_CO+1];
    out_spin2and3[ 5] += GAMMA_X_SPIN1_VAL * in[3*GAMMA_X_SPIN1_CO+2];
    out_spin2and3[ 6] -= in[ 6];
    out_spin2and3[ 7] -= in[ 7];
    out_spin2and3[ 8] -= in[ 8];
    out_spin2and3[ 9] -= in[ 9];
    out_spin2and3[10] -= in[10];
    out_spin2and3[11] -= in[11];
  }

  static inline void twospin2_p_X_PRECISION( const vector_PRECISION out_spin0and1, const vector_PRECISION out_spin2and3, const vector_PRECISION in ) {
    out_spin0and1[ 0] = in[ 0];
    out_spin0and1[ 1] = in[ 1];
    out_spin0and1[ 2] = in[ 2];
    out_spin0and1[ 3] = in[ 3];
    out_spin0and1[ 4] = in[ 4];
    out_spin0and1[ 5] = in[ 5];
    out_spin0and1[ 6] = -GAMMA_X_SPIN2_VAL * in[3*GAMMA_X_SPIN2_CO  ];
    out_spin0and1[ 7] = -GAMMA_X_SPIN2_VAL * in[3*GAMMA_X_SPIN2_CO+1];
    out_spin0and1[ 8] = -GAMMA_X_SPIN2_VAL * in[3*GAMMA_X_SPIN2_CO+2];
    out_spin0and1[ 9] = -GAMMA_X_SPIN3_VAL * in[3*GAMMA_X_SPIN3_CO  ];
    out_spin0and1[10] = -GAMMA_X_SPIN3_VAL * in[3*GAMMA_X_SPIN3_CO+1];
    out_spin0and1[11] = -GAMMA_X_SPIN3_VAL * in[3*GAMMA_X_SPIN3_CO+2];
    out_spin2and3[ 0] = -GAMMA_X_SPIN0_VAL * in[3*GAMMA_X_SPIN0_CO  ];
    out_spin2and3[ 1] = -GAMMA_X_SPIN0_VAL * in[3*GAMMA_X_SPIN0_CO+1];
    out_spin2and3[ 2] = -GAMMA_X_SPIN0_VAL * in[3*GAMMA_X_SPIN0_CO+2];
    out_spin2and3[ 3] = -GAMMA_X_SPIN1_VAL * in[3*GAMMA_X_SPIN1_CO  ];
    out_spin2and3[ 4] = -GAMMA_X_SPIN1_VAL * in[3*GAMMA_X_SPIN1_CO+1];
    out_spin2and3[ 5] = -GAMMA_X_SPIN1_VAL * in[3*GAMMA_X_SPIN1_CO+2];
    out_spin2and3[ 6] = in[ 6];
    out_spin2and3[ 7] = in[ 7];
    out_spin2and3[ 8] = in[ 8];
    out_spin2and3[ 9] = in[ 9];
    out_spin2and3[10] = in[10];
    out_spin2and3[11] = in[11];
  }

  static inline void twospin_n_X_PRECISION( const vector_PRECISION out_spin0and1, const vector_PRECISION out_spin2and3, const vector_PRECISION in ) {
    out_spin0and1[ 0] -= in[ 0];
    out_spin0and1[ 1] -= in[ 1];
    out_spin0and1[ 2] -= in[ 2];
    out_spin0and1[ 3] -= in[ 3];
    out_spin0and1[ 4] -= in[ 4];
    out_spin0and1[ 5] -= in[ 5];
    out_spin0and1[ 6] -= GAMMA_X_SPIN2_VAL * in[3*GAMMA_X_SPIN2_CO  ];
    out_spin0and1[ 7] -= GAMMA_X_SPIN2_VAL * in[3*GAMMA_X_SPIN2_CO+1];
    out_spin0and1[ 8] -= GAMMA_X_SPIN2_VAL * in[3*GAMMA_X_SPIN2_CO+2];
    out_spin0and1[ 9] -= GAMMA_X_SPIN3_VAL * in[3*GAMMA_X_SPIN3_CO  ];
    out_spin0and1[10] -= GAMMA_X_SPIN3_VAL * in[3*GAMMA_X_SPIN3_CO+1];
    out_spin0and1[11] -= GAMMA_X_SPIN3_VAL * in[3*GAMMA_X_SPIN3_CO+2];
    out_spin2and3[ 0] -= GAMMA_X_SPIN0_VAL * in[3*GAMMA_X_SPIN0_CO  ];
    out_spin2and3[ 1] -= GAMMA_X_SPIN0_VAL * in[3*GAMMA_X_SPIN0_CO+1];
    out_spin2and3[ 2] -= GAMMA_X_SPIN0_VAL * in[3*GAMMA_X_SPIN0_CO+2];
    out_spin2and3[ 3] -= GAMMA_X_SPIN1_VAL * in[3*GAMMA_X_SPIN1_CO  ];
    out_spin2and3[ 4] -= GAMMA_X_SPIN1_VAL * in[3*GAMMA_X_SPIN1_CO+1];
    out_spin2and3[ 5] -= GAMMA_X_SPIN1_VAL * in[3*GAMMA_X_SPIN1_CO+2];
    out_spin2and3[ 6] -= in[ 6];
    out_spin2and3[ 7] -= in[ 7];
    out_spin2and3[ 8] -= in[ 8];
    out_spin2and3[ 9] -= in[ 9];
    out_spin2and3[10] -= in[10];
    out_spin2and3[11] -= in[11];
  }

  static inline void spin0and1_site_clover_PRECISION( const vector_PRECISION eta, const vector_PRECISION phi, const config_PRECISION clover ) {
    // diagonal
    eta[ 0] = clover[ 0]*phi[ 0];
    eta[ 1] = clover[ 1]*phi[ 1];
    eta[ 2] = clover[ 2]*phi[ 2];
    eta[ 3] = clover[ 3]*phi[ 3];
    eta[ 4] = clover[ 4]*phi[ 4];
    eta[ 5] = clover[ 5]*phi[ 5];
    eta[ 6] = _COMPLEX_PRECISION_ZERO;
    eta[ 7] = _COMPLEX_PRECISION_ZERO;
    eta[ 8] = _COMPLEX_PRECISION_ZERO;
    eta[ 9] = _COMPLEX_PRECISION_ZERO;
    eta[10] = _COMPLEX_PRECISION_ZERO;
    eta[11] = _COMPLEX_PRECISION_ZERO;
    // spin 0 and 1
    eta[0] += clover[12]*phi[1];
    eta[0] += clover[13]*phi[2];
    eta[0] += clover[14]*phi[3];
    eta[0] += clover[15]*phi[4];
    eta[0] += clover[16]*phi[5];
    eta[1] += clover[17]*phi[2];
    eta[1] += clover[18]*phi[3];
    eta[1] += clover[19]*phi[4];
    eta[1] += clover[20]*phi[5];
    eta[2] += clover[21]*phi[3];
    eta[2] += clover[22]*phi[4];
    eta[2] += clover[23]*phi[5];
    eta[3] += clover[24]*phi[4];
    eta[3] += clover[25]*phi[5];
    eta[4] += clover[26]*phi[5];
    eta[1] += conj_PRECISION(clover[12])*phi[0];
    eta[2] += conj_PRECISION(clover[13])*phi[0];
    eta[3] += conj_PRECISION(clover[14])*phi[0];
    eta[4] += conj_PRECISION(clover[15])*phi[0];
    eta[5] += conj_PRECISION(clover[16])*phi[0];
    eta[2] += conj_PRECISION(clover[17])*phi[1];
    eta[3] += conj_PRECISION(clover[18])*phi[1];
    eta[4] += conj_PRECISION(clover[19])*phi[1];
    eta[5] += conj_PRECISION(clover[20])*phi[1];
    eta[3] += conj_PRECISION(clover[21])*phi[2];
    eta[4] += conj_PRECISION(clover[22])*phi[2];
    eta[5] += conj_PRECISION(clover[23])*phi[2];
    eta[4] += conj_PRECISION(clover[24])*phi[3];
    eta[5] += conj_PRECISION(clover[25])*phi[3];
    eta[5] += conj_PRECISION(clover[26])*phi[4];
  }

  static inline void spin2and3_site_clover_PRECISION( const vector_PRECISION eta, const vector_PRECISION phi, const config_PRECISION clover ) {
    // diagonal
    eta[ 0] = _COMPLEX_PRECISION_ZERO;
    eta[ 1] = _COMPLEX_PRECISION_ZERO;
    eta[ 2] = _COMPLEX_PRECISION_ZERO;
    eta[ 3] = _COMPLEX_PRECISION_ZERO;
    eta[ 4] = _COMPLEX_PRECISION_ZERO;
    eta[ 5] = _COMPLEX_PRECISION_ZERO;
    eta[ 6] = clover[ 6]*phi[ 6];
    eta[ 7] = clover[ 7]*phi[ 7];
    eta[ 8] = clover[ 8]*phi[ 8];
    eta[ 9] = clover[ 9]*phi[ 9];
    eta[10] = clover[10]*phi[10];
    eta[11] = clover[11]*phi[11];
    // spin 2 and 3
    eta[ 6] += clover[28]*phi[ 8];
    eta[ 6] += clover[27]*phi[ 7];
    eta[ 6] += clover[29]*phi[ 9];
    eta[ 6] += clover[30]*phi[10];
    eta[ 6] += clover[31]*phi[11];
    eta[ 7] += clover[32]*phi[ 8];
    eta[ 7] += clover[33]*phi[ 9];
    eta[ 7] += clover[34]*phi[10];
    eta[ 7] += clover[35]*phi[11];
    eta[ 8] += clover[36]*phi[ 9];
    eta[ 8] += clover[37]*phi[10];
    eta[ 8] += clover[38]*phi[11];
    eta[ 9] += clover[39]*phi[10];
    eta[ 9] += clover[40]*phi[11];
    eta[10] += clover[41]*phi[11];
    eta[ 7] += conj_PRECISION(clover[27])*phi[ 6];
    eta[ 8] += conj_PRECISION(clover[28])*phi[ 6];
    eta[ 9] += conj_PRECISION(clover[29])*phi[ 6];
    eta[10] += conj_PRECISION(clover[30])*phi[ 6];
    eta[11] += conj_PRECISION(clover[31])*phi[ 6];
    eta[ 8] += conj_PRECISION(clover[32])*phi[ 7];
    eta[ 9] += conj_PRECISION(clover[33])*phi[ 7];
    eta[10] += conj_PRECISION(clover[34])*phi[ 7];
    eta[11] += conj_PRECISION(clover[35])*phi[ 7];
    eta[ 9] += conj_PRECISION(clover[36])*phi[ 8];
    eta[10] += conj_PRECISION(clover[37])*phi[ 8];
    eta[11] += conj_PRECISION(clover[38])*phi[ 8];
    eta[10] += conj_PRECISION(clover[39])*phi[ 9];
    eta[11] += conj_PRECISION(clover[40])*phi[ 9];
    eta[11] += conj_PRECISION(clover[41])*phi[10];
  }
  
  static inline void site_clover_PRECISION( const vector_PRECISION eta, const vector_PRECISION phi, const config_PRECISION clover ) {
    // diagonal
    eta[ 0] = clover[ 0]*phi[ 0];
    eta[ 1] = clover[ 1]*phi[ 1];
    eta[ 2] = clover[ 2]*phi[ 2];
    eta[ 3] = clover[ 3]*phi[ 3];
    eta[ 4] = clover[ 4]*phi[ 4];
    eta[ 5] = clover[ 5]*phi[ 5];
    eta[ 6] = clover[ 6]*phi[ 6];
    eta[ 7] = clover[ 7]*phi[ 7];
    eta[ 8] = clover[ 8]*phi[ 8];
    eta[ 9] = clover[ 9]*phi[ 9];
    eta[10] = clover[10]*phi[10];
    eta[11] = clover[11]*phi[11];
    // spin 0 and 1, row major
    eta[0] += clover[12]*phi[1];
    eta[0] += clover[13]*phi[2];
    eta[0] += clover[14]*phi[3];
    eta[0] += clover[15]*phi[4];
    eta[0] += clover[16]*phi[5];
    eta[1] += clover[17]*phi[2];
    eta[1] += clover[18]*phi[3];
    eta[1] += clover[19]*phi[4];
    eta[1] += clover[20]*phi[5];
    eta[2] += clover[21]*phi[3];
    eta[2] += clover[22]*phi[4];
    eta[2] += clover[23]*phi[5];
    eta[3] += clover[24]*phi[4];
    eta[3] += clover[25]*phi[5];
    eta[4] += clover[26]*phi[5];
    eta[1] += conj_PRECISION(clover[12])*phi[0];
    eta[2] += conj_PRECISION(clover[13])*phi[0];
    eta[3] += conj_PRECISION(clover[14])*phi[0];
    eta[4] += conj_PRECISION(clover[15])*phi[0];
    eta[5] += conj_PRECISION(clover[16])*phi[0];
    eta[2] += conj_PRECISION(clover[17])*phi[1];
    eta[3] += conj_PRECISION(clover[18])*phi[1];
    eta[4] += conj_PRECISION(clover[19])*phi[1];
    eta[5] += conj_PRECISION(clover[20])*phi[1];
    eta[3] += conj_PRECISION(clover[21])*phi[2];
    eta[4] += conj_PRECISION(clover[22])*phi[2];
    eta[5] += conj_PRECISION(clover[23])*phi[2];
    eta[4] += conj_PRECISION(clover[24])*phi[3];
    eta[5] += conj_PRECISION(clover[25])*phi[3];
    eta[5] += conj_PRECISION(clover[26])*phi[4];
    // spin 2 and 3, row major
    eta[ 6] += clover[27]*phi[ 7];
    eta[ 6] += clover[28]*phi[ 8];
    eta[ 6] += clover[29]*phi[ 9];
    eta[ 6] += clover[30]*phi[10];
    eta[ 6] += clover[31]*phi[11];
    eta[ 7] += clover[32]*phi[ 8];
    eta[ 7] += clover[33]*phi[ 9];
    eta[ 7] += clover[34]*phi[10];
    eta[ 7] += clover[35]*phi[11];
    eta[ 8] += clover[36]*phi[ 9];
    eta[ 8] += clover[37]*phi[10];
    eta[ 8] += clover[38]*phi[11];
    eta[ 9] += clover[39]*phi[10];
    eta[ 9] += clover[40]*phi[11];
    eta[10] += clover[41]*phi[11];
    eta[ 7] += conj_PRECISION(clover[27])*phi[ 6];
    eta[ 8] += conj_PRECISION(clover[28])*phi[ 6];
    eta[ 9] += conj_PRECISION(clover[29])*phi[ 6];
    eta[10] += conj_PRECISION(clover[30])*phi[ 6];
    eta[11] += conj_PRECISION(clover[31])*phi[ 6];
    eta[ 8] += conj_PRECISION(clover[32])*phi[ 7];
    eta[ 9] += conj_PRECISION(clover[33])*phi[ 7];
    eta[10] += conj_PRECISION(clover[34])*phi[ 7];
    eta[11] += conj_PRECISION(clover[35])*phi[ 7];
    eta[ 9] += conj_PRECISION(clover[36])*phi[ 8];
    eta[10] += conj_PRECISION(clover[37])*phi[ 8];
    eta[11] += conj_PRECISION(clover[38])*phi[ 8];
    eta[10] += conj_PRECISION(clover[39])*phi[ 9];
    eta[11] += conj_PRECISION(clover[40])*phi[ 9];
    eta[11] += conj_PRECISION(clover[41])*phi[10];
  }
  
#endif 
