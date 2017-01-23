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
 
#ifndef CLIFFORD_HEADER
  #define CLIFFORD_HEADER
  // assertion: gamma5 = (+/-) diag( 1, 1, -1, -1 )
  
  // choose basis:
  // BASIS0: Basis used within the OPENQCD/DD-HMC code
  // BASIS1: Basis used for BMW-c lattices
  // BASIS2: Basis used for QCDSF lattices
  // BASIS3: Basis used in the QOPQDP Code
  #define BASIS0 // change here
  
  enum { T, Z, Y, X };
  
  #ifndef I
    #define I _Complex_I
  #endif
  
  #ifdef BASIS0
    // Basis used within the OPENQCD/DD-HMC code
    #define CLIFFORD_BASIS "BASIS0:OPENQCD/DD-HMC BASIS"
    /* gamma_T = 
     *  0  0 -1  0
     *  0  0  0 -1
     * -1  0  0  0
     *  0 -1  0  0
     */
    #define GAMMA_T_SPIN0_CO   2
    #define GAMMA_T_SPIN0_VAL -1
    #define GAMMA_T_SPIN1_CO   3
    #define GAMMA_T_SPIN1_VAL -1
    #define GAMMA_T_SPIN2_CO   0
    #define GAMMA_T_SPIN2_VAL -1
    #define GAMMA_T_SPIN3_CO   1
    #define GAMMA_T_SPIN3_VAL -1
    
    /* gamma_Z =
     *  0  0  0 -I
     *  0  0 -I  0
     *  0  I  0  0
     *  I  0  0  0
     */
    #define GAMMA_Z_SPIN0_CO   3
    #define GAMMA_Z_SPIN0_VAL -I
    #define GAMMA_Z_SPIN1_CO   2
    #define GAMMA_Z_SPIN1_VAL -I
    #define GAMMA_Z_SPIN2_CO   1
    #define GAMMA_Z_SPIN2_VAL  I
    #define GAMMA_Z_SPIN3_CO   0
    #define GAMMA_Z_SPIN3_VAL  I
    
    /* gamma_Y =
     *  0  0  0 -1
     *  0  0  1  0
     *  0  1  0  0
     * -1  0  0  0
     */
    #define GAMMA_Y_SPIN0_CO   3
    #define GAMMA_Y_SPIN0_VAL -1
    #define GAMMA_Y_SPIN1_CO   2
    #define GAMMA_Y_SPIN1_VAL  1
    #define GAMMA_Y_SPIN2_CO   1
    #define GAMMA_Y_SPIN2_VAL  1
    #define GAMMA_Y_SPIN3_CO   0
    #define GAMMA_Y_SPIN3_VAL -1   
    
    /* gamma_X = 
     *  0  0 -I  0
     *  0  0  0  I
     *  I  0  0  0
     *  0 -I  0  0
     */
    #define GAMMA_X_SPIN0_CO   2
    #define GAMMA_X_SPIN0_VAL -I
    #define GAMMA_X_SPIN1_CO   3
    #define GAMMA_X_SPIN1_VAL  I
    #define GAMMA_X_SPIN2_CO   0
    #define GAMMA_X_SPIN2_VAL  I
    #define GAMMA_X_SPIN3_CO   1
    #define GAMMA_X_SPIN3_VAL -I
#ifdef SSE
    #define GAMMA_T_SPIN0_RE_SIGN -1
    #define GAMMA_T_SPIN1_RE_SIGN -1
    #define GAMMA_T_SPIN2_RE_SIGN -1
    #define GAMMA_T_SPIN3_RE_SIGN -1
    #define GAMMA_T_SPIN0_IM_SIGN -1
    #define GAMMA_T_SPIN1_IM_SIGN -1
    #define GAMMA_T_SPIN2_IM_SIGN -1
    #define GAMMA_T_SPIN3_IM_SIGN -1
    #define GAMMA_T_SPIN0_OFFSET 0
    #define GAMMA_T_SPIN1_OFFSET 0
    #define GAMMA_T_SPIN2_OFFSET 0
    #define GAMMA_T_SPIN3_OFFSET 0
    
    #define GAMMA_Z_SPIN0_RE_SIGN +1
    #define GAMMA_Z_SPIN1_RE_SIGN +1
    #define GAMMA_Z_SPIN2_RE_SIGN -1
    #define GAMMA_Z_SPIN3_RE_SIGN -1
    #define GAMMA_Z_SPIN0_IM_SIGN -1
    #define GAMMA_Z_SPIN1_IM_SIGN -1
    #define GAMMA_Z_SPIN2_IM_SIGN +1
    #define GAMMA_Z_SPIN3_IM_SIGN +1
    #define GAMMA_Z_SPIN0_OFFSET 1
    #define GAMMA_Z_SPIN1_OFFSET 1
    #define GAMMA_Z_SPIN2_OFFSET 1
    #define GAMMA_Z_SPIN3_OFFSET 1
    
    #define GAMMA_Y_SPIN0_RE_SIGN -1
    #define GAMMA_Y_SPIN1_RE_SIGN +1
    #define GAMMA_Y_SPIN2_RE_SIGN +1
    #define GAMMA_Y_SPIN3_RE_SIGN -1
    #define GAMMA_Y_SPIN0_IM_SIGN -1
    #define GAMMA_Y_SPIN1_IM_SIGN +1
    #define GAMMA_Y_SPIN2_IM_SIGN +1
    #define GAMMA_Y_SPIN3_IM_SIGN -1
    #define GAMMA_Y_SPIN0_OFFSET 0
    #define GAMMA_Y_SPIN1_OFFSET 0
    #define GAMMA_Y_SPIN2_OFFSET 0
    #define GAMMA_Y_SPIN3_OFFSET 0
    
    #define GAMMA_X_SPIN0_RE_SIGN +1
    #define GAMMA_X_SPIN1_RE_SIGN -1
    #define GAMMA_X_SPIN2_RE_SIGN -1
    #define GAMMA_X_SPIN3_RE_SIGN +1
    #define GAMMA_X_SPIN0_IM_SIGN -1
    #define GAMMA_X_SPIN1_IM_SIGN +1
    #define GAMMA_X_SPIN2_IM_SIGN +1
    #define GAMMA_X_SPIN3_IM_SIGN -1
    #define GAMMA_X_SPIN0_OFFSET 1
    #define GAMMA_X_SPIN1_OFFSET 1
    #define GAMMA_X_SPIN2_OFFSET 1
    #define GAMMA_X_SPIN3_OFFSET 1
    
    #define GAMMA_T_SHUFFLE(A) A
    #define GAMMA_Z_SHUFFLE(A) _mm_shuffle_pd(A,A,_MM_SHUFFLE2(0,1))
    #define GAMMA_Y_SHUFFLE(A) A
    #define GAMMA_X_SHUFFLE(A) _mm_shuffle_pd(A,A,_MM_SHUFFLE2(0,1))
#endif
  
  /* ------------------------------------------------- */
  #else
    #ifdef BASIS1
      // Basis used for BMW-c lattices
      #define CLIFFORD_BASIS "BASIS1:BMW-c BASIS"
      // Stored negatively due to different interpretations
      // of the Wilson Dirac operator
      /* gamma_T =
      *  0  0  1  0 
      *  0  0  0  1
      *  1  0  0  0
      *  0  1  0  0  
      */
      #define GAMMA_T_SPIN0_CO   2
      #define GAMMA_T_SPIN0_VAL -1
      #define GAMMA_T_SPIN1_CO   3
      #define GAMMA_T_SPIN1_VAL -1
      #define GAMMA_T_SPIN2_CO   0
      #define GAMMA_T_SPIN2_VAL -1
      #define GAMMA_T_SPIN3_CO   1
      #define GAMMA_T_SPIN3_VAL -1

      /* gamma_Z =
      *  0  0  I  0
      *  0  0  0 -I
      * -I  0  0  0
      *  0  I  0  0  
      */
      #define GAMMA_Z_SPIN0_CO   2
      #define GAMMA_Z_SPIN0_VAL -I
      #define GAMMA_Z_SPIN1_CO   3
      #define GAMMA_Z_SPIN1_VAL  I
      #define GAMMA_Z_SPIN2_CO   0
      #define GAMMA_Z_SPIN2_VAL  I
      #define GAMMA_Z_SPIN3_CO   1
      #define GAMMA_Z_SPIN3_VAL -I

      /* gamma_Y =
      *  0  0  0 -1
      *  0  0  1  0
      *  0  1  0  0
      * -1  0  0  0  
      */
      #define GAMMA_Y_SPIN0_CO   3
      #define GAMMA_Y_SPIN0_VAL  1
      #define GAMMA_Y_SPIN1_CO   2
      #define GAMMA_Y_SPIN1_VAL -1
      #define GAMMA_Y_SPIN2_CO   1
      #define GAMMA_Y_SPIN2_VAL -1
      #define GAMMA_Y_SPIN3_CO   0
      #define GAMMA_Y_SPIN3_VAL  1

      /* gamma_X =
      *  0  0  0  I
      *  0  0  I  0
      *  0 -I  0  0
      * -I  0  0  0  
      */
      #define GAMMA_X_SPIN0_CO   3
      #define GAMMA_X_SPIN0_VAL -I
      #define GAMMA_X_SPIN1_CO   2
      #define GAMMA_X_SPIN1_VAL -I
      #define GAMMA_X_SPIN2_CO   1
      #define GAMMA_X_SPIN2_VAL  I
      #define GAMMA_X_SPIN3_CO   0
      #define GAMMA_X_SPIN3_VAL  I
#ifdef SSE
      #define GAMMA_T_SPIN0_RE_SIGN -1
      #define GAMMA_T_SPIN1_RE_SIGN -1
      #define GAMMA_T_SPIN2_RE_SIGN -1
      #define GAMMA_T_SPIN3_RE_SIGN -1
      #define GAMMA_T_SPIN0_IM_SIGN -1
      #define GAMMA_T_SPIN1_IM_SIGN -1
      #define GAMMA_T_SPIN2_IM_SIGN -1
      #define GAMMA_T_SPIN3_IM_SIGN -1
      #define GAMMA_T_SPIN0_OFFSET 0
      #define GAMMA_T_SPIN1_OFFSET 0
      #define GAMMA_T_SPIN2_OFFSET 0
      #define GAMMA_T_SPIN3_OFFSET 0
      
      #define GAMMA_Z_SPIN0_RE_SIGN +1
      #define GAMMA_Z_SPIN1_RE_SIGN -1
      #define GAMMA_Z_SPIN2_RE_SIGN -1
      #define GAMMA_Z_SPIN3_RE_SIGN +1
      #define GAMMA_Z_SPIN0_IM_SIGN -1
      #define GAMMA_Z_SPIN1_IM_SIGN +1
      #define GAMMA_Z_SPIN2_IM_SIGN +1
      #define GAMMA_Z_SPIN3_IM_SIGN -1
      #define GAMMA_Z_SPIN0_OFFSET 1
      #define GAMMA_Z_SPIN1_OFFSET 1
      #define GAMMA_Z_SPIN2_OFFSET 1
      #define GAMMA_Z_SPIN3_OFFSET 1
      
      #define GAMMA_Y_SPIN0_RE_SIGN +1
      #define GAMMA_Y_SPIN1_RE_SIGN -1
      #define GAMMA_Y_SPIN2_RE_SIGN -1
      #define GAMMA_Y_SPIN3_RE_SIGN +1
      #define GAMMA_Y_SPIN0_IM_SIGN +1
      #define GAMMA_Y_SPIN1_IM_SIGN -1
      #define GAMMA_Y_SPIN2_IM_SIGN -1
      #define GAMMA_Y_SPIN3_IM_SIGN +1
      #define GAMMA_Y_SPIN0_OFFSET 0
      #define GAMMA_Y_SPIN1_OFFSET 0
      #define GAMMA_Y_SPIN2_OFFSET 0
      #define GAMMA_Y_SPIN3_OFFSET 0
      
      #define GAMMA_X_SPIN0_RE_SIGN +1
      #define GAMMA_X_SPIN1_RE_SIGN +1
      #define GAMMA_X_SPIN2_RE_SIGN -1
      #define GAMMA_X_SPIN3_RE_SIGN -1
      #define GAMMA_X_SPIN0_IM_SIGN -1
      #define GAMMA_X_SPIN1_IM_SIGN -1
      #define GAMMA_X_SPIN2_IM_SIGN +1
      #define GAMMA_X_SPIN3_IM_SIGN +1
      #define GAMMA_X_SPIN0_OFFSET 1
      #define GAMMA_X_SPIN1_OFFSET 1
      #define GAMMA_X_SPIN2_OFFSET 1
      #define GAMMA_X_SPIN3_OFFSET 1
      
      #define GAMMA_T_SHUFFLE(A) A
      #define GAMMA_Z_SHUFFLE(A) _mm_shuffle_pd(A,A,_MM_SHUFFLE2(0,1))
      #define GAMMA_Y_SHUFFLE(A) A
      #define GAMMA_X_SHUFFLE(A) _mm_shuffle_pd(A,A,_MM_SHUFFLE2(0,1))
#endif    
  /* ------------------------------------------------- */
    #else
      #ifdef BASIS2
        // Bais used for QCDSF lattices
        #define CLIFFORD_BASIS "BASIS2:QCDSF BASIS"
        /* gamma_T =
        *  0  0  1  0
        *  0  0  0  1
        *  1  0  0  0
        *  0  1  0  0  
        */
        #define GAMMA_T_SPIN0_CO   2
        #define GAMMA_T_SPIN0_VAL  1
        #define GAMMA_T_SPIN1_CO   3
        #define GAMMA_T_SPIN1_VAL  1
        #define GAMMA_T_SPIN2_CO   0
        #define GAMMA_T_SPIN2_VAL  1
        #define GAMMA_T_SPIN3_CO   1
        #define GAMMA_T_SPIN3_VAL  1

        /* gamma_Z =
        *  0  0  I  0
        *  0  0  0 -I
        * -I  0  0  0
        *  0  I  0  0  
        */
        #define GAMMA_Z_SPIN0_CO   2
        #define GAMMA_Z_SPIN0_VAL  I
        #define GAMMA_Z_SPIN1_CO   3
        #define GAMMA_Z_SPIN1_VAL -I
        #define GAMMA_Z_SPIN2_CO   0
        #define GAMMA_Z_SPIN2_VAL -I
        #define GAMMA_Z_SPIN3_CO   1
        #define GAMMA_Z_SPIN3_VAL  I

        /* gamma_Y =
        *  0  0  0 -1
        *  0  0  1  0
        *  0  1  0  0
        * -1  0  0  0  
        */
        #define GAMMA_Y_SPIN0_CO   3
        #define GAMMA_Y_SPIN0_VAL -1
        #define GAMMA_Y_SPIN1_CO   2
        #define GAMMA_Y_SPIN1_VAL  1
        #define GAMMA_Y_SPIN2_CO   1
        #define GAMMA_Y_SPIN2_VAL  1
        #define GAMMA_Y_SPIN3_CO   0
        #define GAMMA_Y_SPIN3_VAL -1

        /* gamma_X =
        *  0  0  0  I
        *  0  0  I  0
        *  0 -I  0  0
        * -I  0  0  0  
        */
        #define GAMMA_X_SPIN0_CO   3
        #define GAMMA_X_SPIN0_VAL  I
        #define GAMMA_X_SPIN1_CO   2
        #define GAMMA_X_SPIN1_VAL  I
        #define GAMMA_X_SPIN2_CO   1
        #define GAMMA_X_SPIN2_VAL -I
        #define GAMMA_X_SPIN3_CO   0
        #define GAMMA_X_SPIN3_VAL -I
#ifdef SSE
        #define GAMMA_T_SPIN0_RE_SIGN +1
        #define GAMMA_T_SPIN1_RE_SIGN +1
        #define GAMMA_T_SPIN2_RE_SIGN +1
        #define GAMMA_T_SPIN3_RE_SIGN +1
        #define GAMMA_T_SPIN0_IM_SIGN +1
        #define GAMMA_T_SPIN1_IM_SIGN +1
        #define GAMMA_T_SPIN2_IM_SIGN +1
        #define GAMMA_T_SPIN3_IM_SIGN +1
        #define GAMMA_T_SPIN0_OFFSET 0
        #define GAMMA_T_SPIN1_OFFSET 0
        #define GAMMA_T_SPIN2_OFFSET 0
        #define GAMMA_T_SPIN3_OFFSET 0
        
        #define GAMMA_Z_SPIN0_RE_SIGN -1
        #define GAMMA_Z_SPIN1_RE_SIGN +1
        #define GAMMA_Z_SPIN2_RE_SIGN +1
        #define GAMMA_Z_SPIN3_RE_SIGN -1
        #define GAMMA_Z_SPIN0_IM_SIGN +1
        #define GAMMA_Z_SPIN1_IM_SIGN -1
        #define GAMMA_Z_SPIN2_IM_SIGN -1
        #define GAMMA_Z_SPIN3_IM_SIGN +1
        #define GAMMA_Z_SPIN0_OFFSET 1
        #define GAMMA_Z_SPIN1_OFFSET 1
        #define GAMMA_Z_SPIN2_OFFSET 1
        #define GAMMA_Z_SPIN3_OFFSET 1
        
        #define GAMMA_Y_SPIN0_RE_SIGN -1
        #define GAMMA_Y_SPIN1_RE_SIGN +1
        #define GAMMA_Y_SPIN2_RE_SIGN +1
        #define GAMMA_Y_SPIN3_RE_SIGN -1
        #define GAMMA_Y_SPIN0_IM_SIGN -1
        #define GAMMA_Y_SPIN1_IM_SIGN +1
        #define GAMMA_Y_SPIN2_IM_SIGN +1
        #define GAMMA_Y_SPIN3_IM_SIGN -1
        #define GAMMA_Y_SPIN0_OFFSET 0
        #define GAMMA_Y_SPIN1_OFFSET 0
        #define GAMMA_Y_SPIN2_OFFSET 0
        #define GAMMA_Y_SPIN3_OFFSET 0
        
        #define GAMMA_X_SPIN0_RE_SIGN -1
        #define GAMMA_X_SPIN1_RE_SIGN -1
        #define GAMMA_X_SPIN2_RE_SIGN +1
        #define GAMMA_X_SPIN3_RE_SIGN +1
        #define GAMMA_X_SPIN0_IM_SIGN +1
        #define GAMMA_X_SPIN1_IM_SIGN +1
        #define GAMMA_X_SPIN2_IM_SIGN -1
        #define GAMMA_X_SPIN3_IM_SIGN -1
        #define GAMMA_X_SPIN0_OFFSET 1
        #define GAMMA_X_SPIN1_OFFSET 1
        #define GAMMA_X_SPIN2_OFFSET 1
        #define GAMMA_X_SPIN3_OFFSET 1
 
        #define GAMMA_T_SHUFFLE(A) A
        #define GAMMA_Z_SHUFFLE(A) _mm_shuffle_pd(A,A,_MM_SHUFFLE2(0,1))
        #define GAMMA_Y_SHUFFLE(A) A
        #define GAMMA_X_SHUFFLE(A) _mm_shuffle_pd(A,A,_MM_SHUFFLE2(0,1))
#endif
      #else
        #ifdef BASIS3
          // Basis used in the QOPQDP Code (by James Osborn/USQCD)
          #define CLIFFORD_BASIS "BASIS3:QOPQDP BASIS"
          /* gamma_T =
          *  0  0  1  0 
          *  0  0  0  1
          *  1  0  0  0
          *  0  1  0  0  
          */
          #define GAMMA_T_SPIN0_CO   2
          #define GAMMA_T_SPIN0_VAL  1
          #define GAMMA_T_SPIN1_CO   3
          #define GAMMA_T_SPIN1_VAL  1
          #define GAMMA_T_SPIN2_CO   0
          #define GAMMA_T_SPIN2_VAL  1
          #define GAMMA_T_SPIN3_CO   1
          #define GAMMA_T_SPIN3_VAL  1

          /* gamma_Z =
          *  0  0  0  I
          *  0  0  I  0
          *  0 -I  0  0
          * -I  0  0  0  
          */
          #define GAMMA_Z_SPIN0_CO   3
          #define GAMMA_Z_SPIN0_VAL  I
          #define GAMMA_Z_SPIN1_CO   2
          #define GAMMA_Z_SPIN1_VAL  I
          #define GAMMA_Z_SPIN2_CO   1
          #define GAMMA_Z_SPIN2_VAL -I
          #define GAMMA_Z_SPIN3_CO   0
          #define GAMMA_Z_SPIN3_VAL -I
          
          /* gamma_Y =
          *  0  0  0 -1
          *  0  0  1  0
          *  0  1  0  0
          * -1  0  0  0  
          */
          #define GAMMA_Y_SPIN0_CO   3
          #define GAMMA_Y_SPIN0_VAL -1
          #define GAMMA_Y_SPIN1_CO   2
          #define GAMMA_Y_SPIN1_VAL  1
          #define GAMMA_Y_SPIN2_CO   1
          #define GAMMA_Y_SPIN2_VAL  1
          #define GAMMA_Y_SPIN3_CO   0
          #define GAMMA_Y_SPIN3_VAL -1
          
          /* gamma_X =
          *  0  0  I  0
          *  0  0  0 -I
          * -I  0  0  0
          *  0  I  0  0  
          */
          #define GAMMA_X_SPIN0_CO   2
          #define GAMMA_X_SPIN0_VAL  I
          #define GAMMA_X_SPIN1_CO   3
          #define GAMMA_X_SPIN1_VAL -I
          #define GAMMA_X_SPIN2_CO   0
          #define GAMMA_X_SPIN2_VAL -I
          #define GAMMA_X_SPIN3_CO   1
          #define GAMMA_X_SPIN3_VAL  I
#ifdef SSE
          #define GAMMA_T_SPIN0_RE_SIGN +1
          #define GAMMA_T_SPIN1_RE_SIGN +1
          #define GAMMA_T_SPIN2_RE_SIGN +1
          #define GAMMA_T_SPIN3_RE_SIGN +1
          #define GAMMA_T_SPIN0_IM_SIGN +1
          #define GAMMA_T_SPIN1_IM_SIGN +1
          #define GAMMA_T_SPIN2_IM_SIGN +1
          #define GAMMA_T_SPIN3_IM_SIGN +1
          #define GAMMA_T_SPIN0_OFFSET 0
          #define GAMMA_T_SPIN1_OFFSET 0
          #define GAMMA_T_SPIN2_OFFSET 0
          #define GAMMA_T_SPIN3_OFFSET 0
          
          #define GAMMA_Z_SPIN0_RE_SIGN -1
          #define GAMMA_Z_SPIN1_RE_SIGN -1
          #define GAMMA_Z_SPIN2_RE_SIGN +1
          #define GAMMA_Z_SPIN3_RE_SIGN +1
          #define GAMMA_Z_SPIN0_IM_SIGN +1
          #define GAMMA_Z_SPIN1_IM_SIGN +1
          #define GAMMA_Z_SPIN2_IM_SIGN -1
          #define GAMMA_Z_SPIN3_IM_SIGN -1
          #define GAMMA_Z_SPIN0_OFFSET 1
          #define GAMMA_Z_SPIN1_OFFSET 1
          #define GAMMA_Z_SPIN2_OFFSET 1
          #define GAMMA_Z_SPIN3_OFFSET 1
          
          #define GAMMA_Y_SPIN0_RE_SIGN -1
          #define GAMMA_Y_SPIN1_RE_SIGN +1
          #define GAMMA_Y_SPIN2_RE_SIGN +1
          #define GAMMA_Y_SPIN3_RE_SIGN -1
          #define GAMMA_Y_SPIN0_IM_SIGN -1
          #define GAMMA_Y_SPIN1_IM_SIGN +1
          #define GAMMA_Y_SPIN2_IM_SIGN +1
          #define GAMMA_Y_SPIN3_IM_SIGN -1
          #define GAMMA_Y_SPIN0_OFFSET 0
          #define GAMMA_Y_SPIN1_OFFSET 0
          #define GAMMA_Y_SPIN2_OFFSET 0
          #define GAMMA_Y_SPIN3_OFFSET 0
          
          #define GAMMA_X_SPIN0_RE_SIGN -1
          #define GAMMA_X_SPIN1_RE_SIGN +1
          #define GAMMA_X_SPIN2_RE_SIGN +1
          #define GAMMA_X_SPIN3_RE_SIGN -1
          #define GAMMA_X_SPIN0_IM_SIGN +1
          #define GAMMA_X_SPIN1_IM_SIGN -1
          #define GAMMA_X_SPIN2_IM_SIGN -1
          #define GAMMA_X_SPIN3_IM_SIGN +1
          #define GAMMA_X_SPIN0_OFFSET 1
          #define GAMMA_X_SPIN1_OFFSET 1
          #define GAMMA_X_SPIN2_OFFSET 1
          #define GAMMA_X_SPIN3_OFFSET 1
          
          #define GAMMA_T_SHUFFLE(A) A
          #define GAMMA_Z_SHUFFLE(A) _mm_shuffle_pd(A,A,_MM_SHUFFLE2(0,1))
          #define GAMMA_Y_SHUFFLE(A) A
          #define GAMMA_X_SHUFFLE(A) _mm_shuffle_pd(A,A,_MM_SHUFFLE2(0,1))
#endif
        #endif
      #endif
    #endif
  #endif

#ifdef SSE
static const int gamma_co[4][4] = {
  {GAMMA_T_SPIN0_CO, GAMMA_T_SPIN1_CO, GAMMA_T_SPIN2_CO, GAMMA_T_SPIN3_CO},
  {GAMMA_Z_SPIN0_CO, GAMMA_Z_SPIN1_CO, GAMMA_Z_SPIN2_CO, GAMMA_Z_SPIN3_CO},
  {GAMMA_Y_SPIN0_CO, GAMMA_Y_SPIN1_CO, GAMMA_Y_SPIN2_CO, GAMMA_Y_SPIN3_CO},
  {GAMMA_X_SPIN0_CO, GAMMA_X_SPIN1_CO, GAMMA_X_SPIN2_CO, GAMMA_X_SPIN3_CO}};

static const double complex gamma_val[4][4] = {
  {GAMMA_T_SPIN0_VAL, GAMMA_T_SPIN1_VAL, GAMMA_T_SPIN2_VAL, GAMMA_T_SPIN3_VAL},
  {GAMMA_Z_SPIN0_VAL, GAMMA_Z_SPIN1_VAL, GAMMA_Z_SPIN2_VAL, GAMMA_Z_SPIN3_VAL},
  {GAMMA_Y_SPIN0_VAL, GAMMA_Y_SPIN1_VAL, GAMMA_Y_SPIN2_VAL, GAMMA_Y_SPIN3_VAL},
  {GAMMA_X_SPIN0_VAL, GAMMA_X_SPIN1_VAL, GAMMA_X_SPIN2_VAL, GAMMA_X_SPIN3_VAL}};
  
static const int gamma_offset[4][4] = { 
  {GAMMA_T_SPIN0_OFFSET,GAMMA_T_SPIN1_OFFSET,GAMMA_T_SPIN2_OFFSET,GAMMA_T_SPIN3_OFFSET},
  {GAMMA_Z_SPIN0_OFFSET,GAMMA_Z_SPIN1_OFFSET,GAMMA_Z_SPIN2_OFFSET,GAMMA_Z_SPIN3_OFFSET},
  {GAMMA_Y_SPIN0_OFFSET,GAMMA_Y_SPIN1_OFFSET,GAMMA_Y_SPIN2_OFFSET,GAMMA_Y_SPIN3_OFFSET},
  {GAMMA_X_SPIN0_OFFSET,GAMMA_X_SPIN1_OFFSET,GAMMA_X_SPIN2_OFFSET,GAMMA_X_SPIN3_OFFSET}};
  
static const int gamma_re_sign[4][4] = { 
  {GAMMA_T_SPIN0_RE_SIGN,GAMMA_T_SPIN1_RE_SIGN,GAMMA_T_SPIN2_RE_SIGN,GAMMA_T_SPIN3_RE_SIGN},
  {GAMMA_Z_SPIN0_RE_SIGN,GAMMA_Z_SPIN1_RE_SIGN,GAMMA_Z_SPIN2_RE_SIGN,GAMMA_Z_SPIN3_RE_SIGN},
  {GAMMA_Y_SPIN0_RE_SIGN,GAMMA_Y_SPIN1_RE_SIGN,GAMMA_Y_SPIN2_RE_SIGN,GAMMA_Y_SPIN3_RE_SIGN},
  {GAMMA_X_SPIN0_RE_SIGN,GAMMA_X_SPIN1_RE_SIGN,GAMMA_X_SPIN2_RE_SIGN,GAMMA_X_SPIN3_RE_SIGN}};
  
static const int gamma_im_sign[4][4] = { 
  {GAMMA_T_SPIN0_IM_SIGN,GAMMA_T_SPIN1_IM_SIGN,GAMMA_T_SPIN2_IM_SIGN,GAMMA_T_SPIN3_IM_SIGN},
  {GAMMA_Z_SPIN0_IM_SIGN,GAMMA_Z_SPIN1_IM_SIGN,GAMMA_Z_SPIN2_IM_SIGN,GAMMA_Z_SPIN3_IM_SIGN},
  {GAMMA_Y_SPIN0_IM_SIGN,GAMMA_Y_SPIN1_IM_SIGN,GAMMA_Y_SPIN2_IM_SIGN,GAMMA_Y_SPIN3_IM_SIGN},
  {GAMMA_X_SPIN0_IM_SIGN,GAMMA_X_SPIN1_IM_SIGN,GAMMA_X_SPIN2_IM_SIGN,GAMMA_X_SPIN3_IM_SIGN}};
#endif
  
#endif
