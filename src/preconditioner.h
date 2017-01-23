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

#ifndef PRECONDITIONER_HEADER
  #define PRECONDITIONER_HEADER
  
  #include "linsolve_double.h"
  #include "linsolve_float.h"
  #include "vcycle_double.h"
  #include "vcycle_float.h"
  #include "schwarz_float.h"
  #include "schwarz_double.h"

  void preconditioner( vector_double phi, vector_double Dphi, vector_double eta,
                       const int res, level_struct *l, struct Thread *threading );
#endif
