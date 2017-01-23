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

#ifndef DDaplhaAMG_PARAMETERS_H
#define DDaplhaAMG_PARAMETERS_H

#define MAX_MG_LEVELS 4
typedef struct dd_alpha_amg_parameters {

    int number_of_levels;

    int global_lattice[MAX_MG_LEVELS][4];
    int local_lattice[MAX_MG_LEVELS][4];
    int block_lattice[MAX_MG_LEVELS][4];

    int mg_basis_vectors[MAX_MG_LEVELS];
    int setup_iterations[MAX_MG_LEVELS];
    int discard_setup_after;
    int update_setup_iterations[MAX_MG_LEVELS];
    int update_setup_after;

    int post_smooth_iterations[MAX_MG_LEVELS];
    int post_smooth_block_iterations[MAX_MG_LEVELS];

    int coarse_grid_iterations;
    int coarse_grid_maximum_number_of_restarts;
    double coarse_grid_tolerance;

    double solver_mass;
    double setup_mass;
    double c_sw;

} dd_alpha_amg_parameters;

#endif
