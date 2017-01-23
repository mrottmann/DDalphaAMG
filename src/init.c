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

complex_double _COMPLEX_double_ONE = (complex_double)1.0;
complex_double _COMPLEX_double_MINUS_ONE = (complex_double)(-1.0);
complex_double _COMPLEX_double_ZERO = (complex_double)0.0;
complex_float  _COMPLEX_float_ONE = (complex_float)1.0;
complex_float  _COMPLEX_float_MINUS_ONE = (complex_float)(-1.0);
complex_float  _COMPLEX_float_ZERO = (complex_float)0.0;


void next_level_setup( vector_double *V, level_struct *l, struct Thread *threading ) {
  
  if ( l->level > 0 ) {
    int mu;

    START_LOCKED_MASTER(threading)

    // allocate storage for next level parameters and initialize them
    MALLOC( l->next_level, level_struct, 1 );
    l_init( l->next_level );
    
    // define next level parameters
    l->next_level->level = l->level-1;
    l->next_level->depth = l->depth+1;
    l->next_level->real_shift = l->real_shift;
    l->next_level->dirac_shift = l->dirac_shift;
    l->next_level->tol = l->tol;
    l->next_level->post_smooth_iter = g.post_smooth_iter[l->depth+1];
    l->next_level->relax_fac = g.relax_fac[l->depth+1];
    l->next_level->block_iter = g.block_iter[l->depth+1];
    l->next_level->setup_iter = g.setup_iter[l->depth+1];
    l->next_level->num_eig_vect = l->level==1?l->num_eig_vect:g.num_eig_vect[l->depth+1];
    l->next_level->num_lattice_site_var = 2 * l->num_eig_vect;
    l->next_level->n_cy = g.ncycle[l->depth+1];
    l->next_level->global_lattice = g.global_lattice[l->depth+1];
    l->next_level->local_lattice = g.local_lattice[l->depth+1];
    l->next_level->block_lattice = g.block_lattice[l->depth+1];
    l->next_level->num_processes = 1;

    for (mu=0; mu<4; mu++) {
      if ( l->depth+2 < g.num_levels )
        l->next_level->coarsening[mu] = g.global_lattice[l->depth+1][mu]/g.global_lattice[l->depth+2][mu];
      else
        l->next_level->coarsening[mu] = g.local_lattice[l->depth+1][mu];
      
      l->next_level->num_processes_dir[mu] = l->next_level->global_lattice[mu]/l->next_level->local_lattice[mu];
      l->next_level->comm_offset[mu] = (l->num_processes_dir[mu]/l->next_level->num_processes_dir[mu])*l->comm_offset[mu];
      l->next_level->global_splitting[mu] = l->next_level->global_lattice[mu] / l->next_level->local_lattice[mu];
      l->next_level->periodic_bc[mu] = l->periodic_bc[mu];
      l->next_level->num_processes *= l->next_level->num_processes_dir[mu];
    }    

    data_layout_init( l->next_level );
    neighbor_define( l->next_level );

    // update threading struct with size info of level
    update_threading(no_threading, l);

    END_LOCKED_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)

    // update threading struct with size info of level
    update_threading(threading, l);
    
    if ( g.mixed_precision ) {
      START_LOCKED_MASTER(threading)
      if ( l->depth == 0 ) fine_level_float_alloc( l );
      level_float_init( l->next_level );
      next_level_float_setup( l );
      END_LOCKED_MASTER(threading)
      if ( l->depth == 0 ) {
        START_LOCKED_MASTER(threading)
        interpolation_float_alloc( l );
        END_LOCKED_MASTER(threading)
        interpolation_float_define( V, l, threading );
        coarse_grid_correction_float_setup( l, threading );
      }
    } else {
      START_LOCKED_MASTER(threading)
      if ( l->depth == 0 ) fine_level_double_alloc( l );
      level_double_init( l->next_level );
      next_level_double_setup( l );
      END_LOCKED_MASTER(threading)
      if ( l->depth == 0 ) {
        START_LOCKED_MASTER(threading)
        interpolation_double_alloc( l );
        END_LOCKED_MASTER(threading)
        interpolation_double_define( V, l, threading );
        coarse_grid_correction_double_setup( l, threading );
      }
    }
  }
  
  if ( l->depth == 0 ) printf0("\ninitial coarse grid correction is defined\n");
}


void next_level_free( level_struct *l ) {
  
  if ( l->level > 0 ) {
    if ( g.mixed_precision ) {
      if ( l->depth == 0 ) fine_level_float_free( l );
      next_level_float_free( l );
    } else {
      if ( l->depth == 0 ) fine_level_double_free( l );
      next_level_double_free( l );
    }
    FREE( l->next_level, level_struct, 1 );
  }
}


void method_setup( vector_double *V, level_struct *l, struct Thread *threading ) {
  
  double t0=0, t1=0;
  
  START_LOCKED_MASTER(threading)
  if ( g.vt.evaluation ) {
    l->dirac_shift = l->real_shift;
    l->level = g.num_levels-1;
  }
  
  prof_float_init( l );
  prof_double_init( l );
  if ( l->depth==0 )
    prof_init( l );
  
  if ( g.method > 0 ) {
    if ( g.mixed_precision == 2 ) {
      fgmres_MP_struct_alloc( g.restart, g.max_restart, l->inner_vector_size,
                              g.tol, _RIGHT, vcycle_float, &(g.p_MP), l );
      g.p.op = &(g.op_double);
#if defined (DEBUG) || defined (TEST_VECTOR_ANALYSIS)
      MALLOC( g.p.b, complex_double, l->inner_vector_size );
      MALLOC( g.p.x, complex_double, l->inner_vector_size );
#endif
    } else {
      fgmres_double_struct_alloc( g.restart, g.max_restart, l->inner_vector_size, g.tol,
                                  _GLOBAL_FGMRES, _RIGHT, preconditioner,
                                  g.method==6?g5D_plus_clover_double:d_plus_clover_double, &(g.p), l );
    }
  }
  else if ( g.method == 0 ) {
    if ( g.mixed_precision == 2 ) {
      fgmres_MP_struct_alloc( g.restart, g.max_restart, l->inner_vector_size,
                              g.tol, _NOTHING, NULL, &(g.p_MP), l );
      g.p.op = &(g.op_double);
#if defined (DEBUG) || defined (TEST_VECTOR_ANALYSIS)
      MALLOC( g.p.b, complex_double, l->inner_vector_size );
      MALLOC( g.p.x, complex_double, l->inner_vector_size );
#endif
    } else {
      fgmres_double_struct_alloc( g.restart, g.max_restart, l->inner_vector_size, g.tol,
                                  _GLOBAL_FGMRES, _NOTHING, NULL, d_plus_clover_double,
                                  &(g.p), l );
    }
  } else if ( g.method == -1 ) {
    fgmres_double_struct_alloc( 4, g.restart*g.max_restart, l->inner_vector_size, g.tol,
                                _GLOBAL_FGMRES, _NOTHING, NULL, d_plus_clover_double, &(g.p), l );
    fine_level_double_alloc( l );
  }
  END_LOCKED_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  
  if ( g.method >= 0 ) {
    START_LOCKED_MASTER(threading)
    t0 = MPI_Wtime();
    if ( g.mixed_precision ) {
      smoother_float_def( l );
      if ( g.method >= 4 && g.odd_even )
        oddeven_setup_float( &(g.op_double), l );
    } else {
      smoother_double_def( l );
      if ( g.method >= 4 && g.odd_even )
        oddeven_setup_double( &(g.op_double), l );
    }
    END_LOCKED_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)
    if ( g.method > 0 )
      if ( g.interpolation && g.num_levels > 1 )
        next_level_setup( V, l, threading );
    START_LOCKED_MASTER(threading)
    t1 = MPI_Wtime();
    g.total_time = t1-t0;
    printf0("elapsed time: %lf seconds\n", t1-t0 );
    END_LOCKED_MASTER(threading)
  }
  START_LOCKED_MASTER(threading)
#ifdef PARAMOUTPUT  
  if ( g.method >= -1 && g.print > 0 && !( g.vt.evaluation && g.vt.re_setup ) ) {
    if(threading->n_core>1) {
      printf0("\nrunning with %d openmp threads per core", threading->n_core );
    }
    printf0("\n+----------------------------------------------------------+\n");
    if ( g.method > 0 ) {
      printf0("| %d-level method                                           |\n", l->level+1);
      if ( l->level > 0 )
        printf0("| postsmoothing %s-cycle                                    |\n", g.kcycle?"K":(l->n_cy>1?"W":"V") );
    }
    switch( g.method ) {
      case -1: printf0("| pure CGN                                                 |\n"); break;
      case  0: printf0("| pure GMRES                                               |\n"); break;
      case  1: printf0("| FGMRES + additive Schwarz                                |\n"); break;
      case  2: printf0("| FGMRES + red-black multiplicative Schwarz                |\n"); break;
      case  3: printf0("| FGMRES + sixteen color multiplicative Schwarz            |\n"); break;
      default: printf0("| FGMRES + GMRES                                           |\n"); break;
    }
    if ( g.method >=0  )
      printf0("|          restart length: %-3d                             |\n", g.restart );
    printf0("|                      m0: %+9.6lf                       |\n", creal(l->dirac_shift) );
    printf0("|                     csw: %+9.6lf                       |\n", g.csw );
    if ( g.method > 0 ) {
      printf0("+----------------------------------------------------------+\n");
      printf0("|%17s cycles: %-6d                          |\n", "preconditioner", l->n_cy );
      printf0("|            inner solver: %-26s      |\n", g.method==4?"GMRES":"minimal residual iteration" );
      printf0("|               precision: %6s                          |\n", g.mixed_precision?"single":"double" );
    }
    for ( int i=0; i<g.num_levels; i++ ) {
      int *gl = g.global_lattice[i], *ll = g.local_lattice[i], *bl = g.block_lattice[i];
      printf0("+---------------------- depth %2d --------------------------+\n", i );
      printf0("|          global lattice: %-3d %-3d %-3d %-3d                 |\n", gl[0], gl[1], gl[2], gl[3] );
      printf0("|           local lattice: %-3d %-3d %-3d %-3d                 |\n", ll[0], ll[1], ll[2], ll[3] );
      if ( g.method > 0 ) {
        printf0("|           block lattice: %-3d %-3d %-3d %-3d                 |\n", bl[0], bl[1], bl[2], bl[3] );
        if ( i+1 < g.num_levels ) {
          printf0("|        post smooth iter: %-3d                             |\n", g.post_smooth_iter[i] );
          printf0("|     smoother inner iter: %-3d                             |\n", g.block_iter[i] );
          printf0("|              setup iter: %-3d                             |\n", g.setup_iter[i] );
          printf0("|            test vectors: %-3d                             |\n", g.num_eig_vect[i] );
        } else {
          printf0("|      coarge grid solver: %-30s  |\n", g.odd_even?"odd even GMRES":"GMRES" );
          printf0("|              iterations: %-6d                          |\n", g.coarse_iter );
          printf0("|                  cycles: %-6d                          |\n", g.coarse_restart );
          printf0("|               tolerance: %-5.0le                           |\n", g.coarse_tol );
        }
      }
    }
    if ( g.method > 0 && g.kcycle > 0 ) {
      printf0("+----------------------------------------------------------+\n");
      printf0("|          K-cycle length: %-6d                          |\n", g.kcycle_restart );
      printf0("|        K-cycle restarts: %-6d                          |\n", g.kcycle_max_restart );
      printf0("|       K-cycle tolerance: %-5.0le                           |\n", g.kcycle_tol );
    }
    printf0("+----------------------------------------------------------+\n");
    printf0("\n");
  }
#endif
  END_LOCKED_MASTER(threading)

  
  
#ifdef DEBUG
  test_routine( l, threading );
#endif
  START_LOCKED_MASTER(threading)
  g.mg_setup_status.gauge_updates_since_last_setup = 0;
  g.mg_setup_status.gauge_updates_since_last_setup_update = 0;
    if ( l->depth==0 )
      prof_print( l );
  END_LOCKED_MASTER(threading)
}


void method_free( level_struct *l ) {
  
  if ( g.method>=0 ) {
    if ( g.mixed_precision ) {
      if ( g.method >= 4 && g.odd_even )
        oddeven_free_float( l );
      smoother_float_free( l );
    } else {
      if ( g.method >= 4 && g.odd_even )
        oddeven_free_double( l );
      smoother_double_free( l );
    }
    if ( g.method > 0 )
      if ( g.interpolation )
        next_level_free( l );
  } else if ( g.method == -1 ) {
    fine_level_double_free( l );
  }

  if ( g.mixed_precision == 2 && g.method >= 0 ) {
    fgmres_MP_struct_free( &(g.p_MP) );
#if defined (DEBUG) || defined (TEST_VECTOR_ANALYSIS)
    FREE( g.p.b, complex_double, l->inner_vector_size );
    FREE( g.p.x, complex_double, l->inner_vector_size );
#endif
  } else {
    fgmres_double_struct_free( &(g.p), l );
  }

}


void method_re_setup( level_struct *l, struct Thread *threading ) {
  START_LOCKED_MASTER(threading)
  method_free( l );
  END_LOCKED_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  method_setup( NULL, l, threading );
}


void method_update( int setup_iter, level_struct *l, struct Thread *threading ) {
  
  if ( g.method > 0 && g.interpolation && g.num_levels > 1 && setup_iter > 0 ) {
    
    double t0=0, t1=0, shift = creal(l->dirac_shift);
    
    START_LOCKED_MASTER(threading)
    g.in_setup = 1;
      if ( l->depth==0 )
        prof_init( l );
    END_LOCKED_MASTER(threading)

    START_MASTER(threading)
    t0 = MPI_Wtime();
    END_MASTER(threading)
    
    if ( g.setup_m0 != shift )
      shift_update( (complex_double)g.setup_m0, l, threading );
    
    if ( g.mixed_precision )
      iterative_float_setup( setup_iter, l, threading );
    else
      iterative_double_setup( setup_iter, l, threading );
        
    if ( g.setup_m0 != shift )
      shift_update( (complex_double)shift, l, threading );
    
    START_MASTER(threading)
    t1 = MPI_Wtime();
    g.total_time = t1-t0;
    printf0("\nperformed %d iterative setup steps\n", setup_iter );
    printf0("elapsed time: %lf seconds (%lf seconds on coarse grid)\n\n", t1-t0, g.coarse_time );
    END_MASTER(threading)
    
#ifdef DEBUG
    test_routine( l, threading );
#endif

    START_LOCKED_MASTER(threading)
    g.in_setup = 0;
    g.mg_setup_status.gauge_updates_since_last_setup_update = 0;
      if ( l->depth==0 )
        prof_print( l );
    END_LOCKED_MASTER(threading)
    
    
  }
}


void method_init( int *argc, char ***argv, level_struct *l ) {
  
/********************************************************************************* 
* Sets up the global and level struct for the method and assignes values 
* according to the inputfile. 
* - int *argc: Argument count of main function. Determines if inputfile is 
*   provided.
* - char ***argv: In case inputfile is provided, contains name of this file.
* 
* CAUTION: changes in this function have no influence on the interface         
* library, cf. dd_alpha_amg_init.                                                    
*********************************************************************************/

  char inputfile[STRINGLENGTH];
  
  if ( *argc > 1 ) {
    strcpy( inputfile, (*argv)[1] );
  } else {
    strcpy( inputfile, "sample.ini" );
  }
  
#ifdef WRITE_LOGFILE
  g.logfile = fopen( "output.log", "w" );
  fprintf(g.logfile,"---------- log file -----------\n\n");
  fflush(g.logfile);
#endif
  
  predefine_rank();
  l_init( l );
  g_init( l );
  lg_in( inputfile, l );
  cart_define( l );
  data_layout_init( l );
  operator_double_alloc( &(g.op_double), _ORDINARY, l ); 
  operator_double_define( &(g.op_double), l );
  MALLOC( g.odd_even_table, int, l->num_inner_lattice_sites );
  define_odd_even_table( l );
}


void method_finalize( level_struct *l ) {
  
  int ls = MAX(g.num_desired_levels,2);
  
  operator_double_free( &(g.op_double), _ORDINARY, l );
  FREE( g.odd_even_table, int, l->num_inner_lattice_sites );
  FREE( g.global_lattice[0], int, 4*ls );
  FREE( g.local_lattice[0], int, 4*ls );
  FREE( g.block_lattice[0], int, 4*ls );
  FREE( g.global_lattice, int*, ls );
  FREE( g.local_lattice, int*, ls );
  FREE( g.block_lattice, int*, ls );
  FREE( g.post_smooth_iter, int, ls );
  FREE( g.ncycle, int, ls );
  FREE( g.relax_fac, double, ls );
  FREE( g.block_iter, int, ls );
  FREE( g.setup_iter, int, ls );
  FREE( g.num_eig_vect, int, ls );
  cart_free( l );
  var_table_free( &(g.vt) );
  
  if ( g.cur_storage )
    warning0("amount of not freed memory/MPIproc: %lf MB\n", g.cur_storage );
  
#ifdef WRITE_LOGFILE
  fprintf(g.logfile,"---------- end of log file -----------\n\n");
  fflush(g.logfile);
  fclose(g.logfile);
#endif
}


int read_parameter( void **save_at, char *search_pattern, char *read_format, int number, FILE *read_from, int set_default ) {

/********************************************************************************* 
* Reads input parameters from provided inputfile.                                                                                                         
* - void **save_at: Points to variable, where parameters are stored.               
* - char *search_pattern: Gives the name of parameter to be read.                      
* - char *read_format: Specifies datatype of parameter.                                
* - int number: Specifies how many numbers need to be read.  
* - FILE *read_from: The inputfile.                         
* - int set_default: Either _NO_DEFAULT_SET for no default value, or DEFAULT_SET,  
*   if a default value is assigned to this parameter.                            
*   If set_default = _NO_DEFAULT_SET, then this parameter MUST be set in the     
*   inputfile.
*                                                                   
* For examples, see lg_in(...) further below.                                  
*********************************************************************************/ 
  
  int i=0, j, k, n=strlen( search_pattern ), match=0;
  char read_pattern[100000], *read_pattern_pt, buffer[50];
  var_table_entry e;
  
  fseek( read_from, 0L, SEEK_SET );
  
  while ( !match && fgets( read_pattern, 100000, read_from ) ) {
  
    k = strlen( read_pattern );
    j = 0;
    for ( i=0; i<k && !match; i++ ) {
      if ( search_pattern[j] == read_pattern[i] )
        j++;
      else
        j = 0;
      if ( j == n ) {
        match = 1;
      }
    }
  }
  
  read_pattern_pt = read_pattern+i;
  while ( *read_pattern_pt == ' ' )
    read_pattern_pt++;
  
  if ( match ) {
    if ( strcmp(read_format,"%s") != 0 ) {
      e.pt = *save_at;
      for ( i=0; i<n-1; i++ ) {
        e.name[i] = search_pattern[i];
      }
      e.name[n-1]='\0';
      if ( strcmp(read_format,"%d") == 0 ) {
        // int
        for ( j=0; j<number; j++ ) {
          sscanf( read_pattern_pt, read_format, &(((int*)*save_at)[j]) );
          sscanf( read_pattern_pt, "%s", buffer );
          read_pattern_pt += strlen( buffer );
          while ( *read_pattern_pt == ' ' )
            read_pattern_pt++;
        }
        sprintf( e.datatype, "int" );
      } else {
        // double
        for ( j=0; j<number; j++ ) {
          sscanf( read_pattern_pt, read_format, &(((double*)*save_at)[j]) );
          sscanf( read_pattern_pt, "%s", buffer );
          read_pattern_pt += strlen( buffer );
          while ( *read_pattern_pt == ' ' )
            read_pattern_pt++;
        }
        sprintf( e.datatype, "double" );
      }
      if ( number == 1 ) {
        var_table_insert( &(g.vt), e );
      }
    } else {
      // string
      sprintf( ((char*)*save_at), "%s", read_pattern_pt );
      ((char*)*save_at)[strlen(read_pattern_pt)-1] = '\0';
    }
  } else {
    if ( !set_default )
      error0("unable to find string \"%s\" --- fatal error\n", search_pattern);
  }
  return match;
}


void l_init( level_struct *l ) {

  level_double_init( l );
  level_float_init( l );
  
  l->x = NULL;
  l->next_level = NULL;
  l->reqs = NULL;
}


void g_init( level_struct *l ) {

  var_table_init( &(g.vt) );
  operator_double_init( &(g.op_double) );
  operator_float_init( &(g.op_float) );
  fgmres_double_struct_init( &(g.p) );
  fgmres_MP_struct_init( &(g.p_MP) );
  g.global_lattice = NULL;
  g.local_lattice = NULL;
  g.block_lattice = NULL;
  g.post_smooth_iter = NULL;
  g.block_iter = NULL;
  g.setup_iter = NULL;
  g.relax_fac = NULL;
  g.gamma = NULL;
  g.odd_even_table = NULL;
  g.two_cnfgs = 0;
  g.cur_storage = 0;
  g.max_storage = 0;
  g.in_setup = 0;
}


void parameter_update( level_struct *l ) {
  
  l->level = g.num_levels-1-l->depth;
  l->post_smooth_iter = g.post_smooth_iter[l->depth];
  l->block_iter = g.block_iter[l->depth];
  l->setup_iter = g.setup_iter[l->depth];
  l->num_eig_vect = g.num_eig_vect[l->depth];
  
  if ( l->level > 0 ) 
    parameter_update( l->next_level );
}


void set_default_global_info() {

  g.in[0] = '\0';
  g.rhs = 1;
  g.propagator_coords[0] = 0;
  g.propagator_coords[1] = 0;
  g.propagator_coords[2] = 0;
  g.propagator_coords[3] = 0;
  g.anti_pbc = 1;
}

void read_global_info( FILE *in ) {

  void *save_pt;
  int match;

  // Note: There is actually no default set for the three following values
  // Though, when using the code as a library, no configuration paths are required.
  save_pt = &(g.in);
  match = read_parameter( &save_pt, "configuration:", "%s", 1, in, _DEFAULT_SET );
  if ( !match ) {
    save_pt = &(g.in);
    match = read_parameter( &save_pt, "hopp cnfg:", "%s", 1, in, _DEFAULT_SET );
    save_pt = &(g.in_clov);
    match = match && read_parameter( &save_pt, "clov cnfg:", "%s", 1, in, _DEFAULT_SET ); 
    g.two_cnfgs = 1;
    if ( !match ) {
      printf00("Caution: No configuration read from disk!\n");
    }  
  }
  save_pt = &(g.in_format); g.in_format = _STANDARD;
  read_parameter( &save_pt, "format:", "%d", 1, in, _DEFAULT_SET );
  
  save_pt = &(g.rhs);
  read_parameter( &save_pt, "right hand side:", "%d", 1, in, _DEFAULT_SET );
  save_pt = g.propagator_coords;
  read_parameter( &save_pt, "propagator coordinates:", "%d", 4, in, _DEFAULT_SET );
  
  if ( g.rhs == 4 ) {
    save_pt = &(g.source_list);
    read_parameter( &save_pt, "source list:", "%s", 1, in, _NO_DEFAULT_SET );
  }
  
  save_pt = &(g.num_levels); g.num_levels = 2;
  read_parameter( &save_pt, "number of levels:", "%d", 1, in, _DEFAULT_SET );
  g.num_desired_levels = g.num_levels;
  save_pt = &(g.anti_pbc); g.anti_pbc = 0;
  read_parameter( &save_pt, "antiperiodic boundary conditions:", "%d", 1, in, _DEFAULT_SET );
  
  save_pt = &(g.num_openmp_processes); g.num_openmp_processes = 1;
  read_parameter( &save_pt, "number of openmp threads:", "%d", 1, in, _DEFAULT_SET );
}

void set_global_info( struct dd_alpha_amg_parameters *amg_params ) {
  g.num_levels = amg_params->number_of_levels;
  g.anti_pbc = 0;
}


int shortest_dir( int* data ) {
  int min=0, mu;
  for ( mu=1; mu<4; mu++ )
    if ( data[mu] < data[min] )
      min = mu;
  return min;
}


int gcd( int a, int b ) {
  if ( b==0 )
    return a;
  return gcd( b, a%b );
}
int lcm( int a, int b ) {
  return ( a*b / gcd( a, b ) );
}


void read_geometry_data( FILE *in, int ls ) {

  void *save_pt;
  char inputstr[STRINGLENGTH];
  int i, mu, nb, nls, nlls, flag;
  
  for ( i=0; i<ls; i++ ) {
    
    // global lattice
    sprintf( inputstr, "d%d global lattice:", i );
    save_pt = g.global_lattice[i];
    if ( i > 0 ) {
      if ( ! read_parameter( &save_pt, inputstr, "%d", 4, in, _DEFAULT_SET ) ) {
        nls = 1;
        for ( mu=0; mu<4; mu++ ) {
          g.global_lattice[i][mu] = g.global_lattice[i-1][mu]/g.block_lattice[i-1][mu];
          nls *= g.global_lattice[i][mu];
        }
        if ( g.odd_even && nls < 2 ) {
          warning0("lattice dimensions not valid for a %d-level method, choosing a %d-level method\n", g.num_levels, i );
          g.num_levels = i; ls = i;
          break;
        }
      }
    } else {
      read_parameter( &save_pt, inputstr, "%d", 4, in, _NO_DEFAULT_SET );
    }
    
    // local lattice
    sprintf( inputstr, "d%d local lattice:", i );
    save_pt = g.local_lattice[i];
    if ( i > 0 ) {
      if ( ! read_parameter( &save_pt, inputstr, "%d", 4, in, _DEFAULT_SET ) ) {
        nls = 1;
        nlls = 1;
        for ( mu=0; mu<4; mu++ ) {
          g.local_lattice[i][mu] = g.local_lattice[i-1][mu]/g.block_lattice[i-1][mu];
          nlls *= g.local_lattice[i][mu];
          nls *= g.global_lattice[i][mu];
        }
        if ( g.odd_even && nlls < 2 ) {
          if ( nls/nlls > 1 ) {
            mu = shortest_dir( g.local_lattice[i] );
            if ( g.global_lattice[i][mu] > g.local_lattice[i][mu] ) {
              g.local_lattice[i][mu] *= lcm( g.local_lattice[i][mu],
                                             g.global_lattice[i][mu]/g.local_lattice[i][mu] );
            }
          }
        }
      }
    } else {
      read_parameter( &save_pt, inputstr, "%d", 4, in, _NO_DEFAULT_SET );
    }
    
    // block lattice
    for ( mu=0; mu<4; mu++ )
      g.block_lattice[i][mu] = 1;
    if ( i<ls-1 ) {
      sprintf( inputstr, "d%d block lattice:", i );
      save_pt = g.block_lattice[i];
      if ( i > 0 ) {
        if ( ! read_parameter( &save_pt, inputstr, "%d", 4, in, _DEFAULT_SET ) ) {
          nls = 1;
          nb = 1;
          flag = 1;
          for ( mu=0; mu<4; mu++ )  {
            if ( DIVIDES( 2, g.global_lattice[i][mu] ) ) {
              g.block_lattice[i][mu] = 2;
            } else if ( DIVIDES( 3, g.global_lattice[i][mu] ) ) {
              g.block_lattice[i][mu] = 3;
            } else {
              warning0("lattice dimensions not valid for a %d-level method, choosing a %d-level method\n", g.num_levels, i+1 );
              g.num_levels = i+1; ls=i+1;
              g.block_lattice[i][mu] = 1;
              flag = 0;
              break;
            }
            nb *= g.local_lattice[i][mu]/g.block_lattice[i][mu];
                        
            if ( g.local_lattice[i][mu] < g.block_lattice[i][mu] ) {
              g.local_lattice[i][mu] *= g.block_lattice[i][mu];
              if ( ! DIVIDES( g.local_lattice[i][mu], g.global_lattice[i][mu] ) ) {
                g.local_lattice[i][mu] /= g.block_lattice[i][mu];
              }
              warning0("lattice dimensions not valid for a %d-level method, choosing a %d-level method\n", g.num_levels, i+1 );
              g.num_levels = i+1; ls=i+1;
              g.block_lattice[i][mu] = 1;
              flag = 0;
              break;
            }
          }
          
          if ( flag == 1 && g.method == 2 && nb == 1 ) {
            mu = shortest_dir( g.local_lattice[i] );
            if ( g.global_lattice[i][mu] > g.local_lattice[i][mu] ) {
              g.local_lattice[i][mu] *= lcm( g.local_lattice[i][mu],
                                            g.global_lattice[i][mu]/g.local_lattice[i][mu] );
            }
          }
          
        }
      } else {
        read_parameter( &save_pt, inputstr, "%d", 4, in, _NO_DEFAULT_SET );
      }
    }
    
#ifdef DEBUG
    printf00("level: %d, gl: %3d %3d %3d %3d\n", i, g.global_lattice[i][0],
            g.global_lattice[i][1],g.global_lattice[i][2],g.global_lattice[i][3] );
    
    printf00("level: %d, ll: %3d %3d %3d %3d\n", i, g.local_lattice[i][0],
            g.local_lattice[i][1],g.local_lattice[i][2],g.local_lattice[i][3] );
            
    printf00("level: %d, bl: %3d %3d %3d %3d\n\n", i, g.block_lattice[i][0],
            g.block_lattice[i][1],g.block_lattice[i][2],g.block_lattice[i][3] );
#endif
            
            
    sprintf( inputstr, "d%d post smooth iter:", i );
    save_pt = &(g.post_smooth_iter[i]); g.post_smooth_iter[i] = 2;
    read_parameter( &save_pt, inputstr, "%d", 1, in, _DEFAULT_SET );
    
    sprintf( inputstr, "d%d preconditioner cycles:", i );
    save_pt = &(g.ncycle[i]); g.ncycle[i] = 1;
    read_parameter( &save_pt, inputstr, "%d", 1, in, _DEFAULT_SET );
    
    
    sprintf( inputstr, "d%d relaxation factor:", i );
    save_pt = &(g.relax_fac[i]); g.relax_fac[i] = 1.0;
    read_parameter( &save_pt, inputstr, "%lf", 1, in, _DEFAULT_SET );
    
    sprintf( inputstr, "d%d block iter:", i );
    save_pt = &(g.block_iter[i]);
    g.block_iter[i] = 4;
    read_parameter( &save_pt, inputstr, "%d", 1, in, _DEFAULT_SET );
    
    sprintf( inputstr, "d%d setup iter:", i );
    save_pt = &(g.setup_iter[i]);
    if ( i==0 ) g.setup_iter[i] = 6;
    else if ( i==1 ) g.setup_iter[i] = 3;
    else if ( i>1 ) g.setup_iter[i] = 2;
    read_parameter( &save_pt, inputstr, "%d", 1, in, _DEFAULT_SET );
    
    sprintf( inputstr, "d%d test vectors:", i );
    
    save_pt = &(g.num_eig_vect[i]);
    if ( i==0 ) g.num_eig_vect[i] = 20;
    else
#ifdef SSE
      g.num_eig_vect[i] = MAX(4,((int)(1.5*g.num_eig_vect[0])-(((int)(1.5*g.num_eig_vect[0]))%SIMD_LENGTH_float)) );
#else
      g.num_eig_vect[i] = 1.5*g.num_eig_vect[0];
#endif
    read_parameter( &save_pt, inputstr, "%d", 1, in, _DEFAULT_SET );
        
  }
}

void set_geometry_data( struct dd_alpha_amg_parameters *amg_params ) {
  int ls = MAX(g.num_levels, 2);
  for ( int i=0; i<ls; i++ ) {
    for ( int mu=0; mu<4; mu++ ) {
      g.global_lattice[i][mu] = amg_params->global_lattice[i][3-mu];
      g.local_lattice[i][mu]  = amg_params->local_lattice[i][3-mu];
      g.block_lattice[i][mu]  = amg_params->block_lattice[i][3-mu];
    }
    g.num_eig_vect[i]     = amg_params->mg_basis_vectors[i];
  }
}

void read_solver_parameters( FILE *in, level_struct *l ) {

  void *save_pt;

  save_pt = &(g.mixed_precision); g.mixed_precision = 2;
  read_parameter( &save_pt, "mixed precision:", "%d", 1, in, _DEFAULT_SET );
  if ( g.num_levels == 1 ) g.interpolation = 0; else {
    save_pt = &(g.interpolation); g.interpolation = 2;
    read_parameter( &save_pt, "interpolation:", "%d", 1, in, _DEFAULT_SET );
  }
  
  save_pt = &(g.randomize); g.randomize = 0;
  read_parameter( &save_pt, "randomize test vectors:", "%d", 1, in, _DEFAULT_SET );
  save_pt = &(g.coarse_iter); g.coarse_iter = 25;
  read_parameter( &save_pt, "coarse grid iterations:", "%d", 1, in, _DEFAULT_SET );
  save_pt = &(g.coarse_restart); g.coarse_restart = 40;
  read_parameter( &save_pt, "coarse grid restarts:", "%d", 1, in, _DEFAULT_SET );
  save_pt = &(g.coarse_tol); g.coarse_tol = 5E-2;
  read_parameter( &save_pt, "coarse grid tolerance:", "%le", 1, in, _DEFAULT_SET );
  save_pt = &(g.odd_even); g.odd_even = 1;
  read_parameter( &save_pt, "odd even preconditioning:", "%d", 1, in, _DEFAULT_SET );

  save_pt = &(l->real_shift);
  read_parameter( &save_pt, "m0:", "%lf", 1, in, _NO_DEFAULT_SET ); // ensuring downward compatibility
  read_parameter( &save_pt, "solver m0:", "%lf", 1, in, _DEFAULT_SET );
  save_pt = &(g.csw);
  read_parameter( &save_pt, "csw:", "%lf", 1, in, _NO_DEFAULT_SET );
  save_pt = &(g.setup_m0);
  read_parameter( &save_pt, "setup m0:", "%lf", 1, in, _DEFAULT_SET );
  
  save_pt = &(g.method); g.method = 2;
  read_parameter( &save_pt, "method:", "%d", 1, in, _DEFAULT_SET );
  save_pt = &(g.restart); g.restart = 10;
  read_parameter( &save_pt, "iterations between restarts:", "%d", 1, in, _DEFAULT_SET );
  save_pt = &(g.max_restart); g.max_restart = 100;
  read_parameter( &save_pt, "maximum of restarts:", "%d", 1, in, _DEFAULT_SET );
  save_pt = &(g.tol); g.tol = 1E-10;
  read_parameter( &save_pt, "tolerance for relative residual:", "%le", 1, in, _DEFAULT_SET );
  save_pt = &(g.print); g.print = 0;
  read_parameter( &save_pt, "print mode:", "%d", 1, in, _DEFAULT_SET );
  
  if ( g.randomize ) {
    srand( time( 0 ) + 1000*g.my_rank );
  } else 
    srand( 1000*g.my_rank );
}

void set_solver_parameters( struct dd_alpha_amg_parameters *amg_params, level_struct *l ) {

  g.mixed_precision = 1;
  g.interpolation = 2;
  g.randomize = 0;
  g.coarse_iter    = amg_params->coarse_grid_iterations;
  g.coarse_restart = amg_params->coarse_grid_maximum_number_of_restarts;
  g.coarse_tol     = amg_params->coarse_grid_tolerance;
  g.odd_even = 1;

  l->real_shift    = amg_params->solver_mass;
  g.setup_m0       = amg_params->setup_mass;
  g.csw            = amg_params->c_sw;

  g.method = 2;
  l->n_cy = 1;
  // outer solver in MG not used not used, -1 should disable allocations
  g.restart = -1;
  g.max_restart = 100;
  g.tol = 1e-10;
  g.print = 1;

  // setting status to maximum to force initial setup
  g.mg_setup_status.gauge_updates_since_last_setup = amg_params->discard_setup_after;
  g.mg_setup_status.gauge_updates_since_last_setup_update = amg_params->update_setup_after;
}


void read_testvector_io_data_if_necessary( FILE *in ) {

  void *save_pt;
  if ( g.interpolation == 4 ) {
    save_pt = &(g.tv_io_single_file); g.tv_io_single_file = 1;
    read_parameter( &save_pt, "test vector io from single file:", "%d", 1, in, _DEFAULT_SET );
    save_pt = &(g.tv_io_file_name);
    read_parameter( &save_pt, "test vector io file name:", "%s", 1, in, _NO_DEFAULT_SET );
  }
}
void read_evaluation_parameters_if_necessary( FILE *in ) {

  void *save_pt;
  save_pt = &(g.vt.evaluation); g.vt.evaluation = 0;
  read_parameter( &save_pt, "evaluation:", "%d", 1, in, _DEFAULT_SET );
  if ( g.vt.evaluation ) {
    save_pt = &(g.vt.scan_var);
    read_parameter( &save_pt, "scan variable:", "%s", 1, in, _NO_DEFAULT_SET );  
    save_pt = &(g.vt.start_val);
    read_parameter( &save_pt, "start value:", "%lf", 1, in, _NO_DEFAULT_SET );
    save_pt = &(g.vt.end_val);
    read_parameter( &save_pt, "end value:", "%lf", 1, in, _NO_DEFAULT_SET );
    save_pt = &(g.vt.step_size);
    read_parameter( &save_pt, "step size:", "%lf", 1, in, _NO_DEFAULT_SET );
    save_pt = &(g.vt.multiplicative);
    read_parameter( &save_pt, "multiplicative:", "%d", 1, in, _NO_DEFAULT_SET );
    save_pt = &(g.vt.shift_update);
    read_parameter( &save_pt, "shift update:", "%d", 1, in, _NO_DEFAULT_SET );
    save_pt = &(g.vt.re_setup);
    read_parameter( &save_pt, "setup update:", "%d", 1, in, _NO_DEFAULT_SET );
    save_pt = &(g.vt.track_error); g.vt.track_error = 0;
    read_parameter( &save_pt, "track error:", "%d", 1, in, _NO_DEFAULT_SET );
    save_pt = &(g.vt.track_cgn_error); g.vt.track_cgn_error = 0;
    read_parameter( &save_pt, "compare with CGN error:", "%d", 1, in, _DEFAULT_SET );
    save_pt = &(g.vt.average_over); g.vt.average_over = 1;
    read_parameter( &save_pt, "average over:", "%d", 1, in, _DEFAULT_SET );
  }
}

void read_kcycle_data( FILE *in ) {

  void *save_pt;
  save_pt = &(g.kcycle); g.kcycle = 1;
  read_parameter( &save_pt, "kcycle:", "%d", 1, in, _DEFAULT_SET );
  save_pt = &(g.kcycle_restart); g.kcycle_restart = 5;
  read_parameter( &save_pt, "kcycle length:", "%d", 1, in, _DEFAULT_SET );
  save_pt = &(g.kcycle_max_restart); g.kcycle_max_restart = 2;
  read_parameter( &save_pt, "kcycle restarts:", "%d", 1, in, _DEFAULT_SET );
  save_pt = &(g.kcycle_tol); g.kcycle_tol = 1E-1;
  read_parameter( &save_pt, "kcycle tolerance:", "%le", 1, in, _DEFAULT_SET );
}

void set_kcycle_data() {

  g.kcycle = 1;
  g.kcycle_restart = 5;
  g.kcycle_max_restart = 2;
  g.kcycle_tol = 1e-1;
}

void validate_parameters( int ls, level_struct *l ) {

  int i;
  int mu;
  
#ifdef SSE
  if ( !g.odd_even )
    warning0("The SSE implementation is based on the odd-even preconditioned code.\
    \n         Switch on odd-even preconditioning in the input file.\n");
  ASSERT( g.odd_even );
#endif
  
  if ( g.method == 5 && g.interpolation != 0 ) {
    warning0("Multigrid with BiCGstab smoothing is not supported.\n         Switching to FGMRES preconditioned with BiCGstab (g.interpolation=0).\n");
    g.interpolation = 0;
  }

  ASSERT( ASCENDING( 0, g.rhs, 2 ) );
  ASSERT( ASCENDING( -1, g.method, 5 ) );
  
  ASSERT( IMPLIES( g.vt.evaluation, g.rhs <= 2 ) );
#ifdef _20TV
  ASSERT( g.num_eig_vect == 20 );
#endif
  
  for ( i=0; i<g.num_levels; i++ )
    for ( mu=0; mu<4; mu++)
      ASSERT( DIVIDES( g.local_lattice[i][mu], g.global_lattice[i][mu] ) );
  
  for ( i=0; i<g.num_levels-1; i++ )
    for ( mu=0; mu<4; mu++) {
      ASSERT( DIVIDES( g.global_lattice[i+1][mu], g.global_lattice[i][mu] ) );
      ASSERT( DIVIDES( g.block_lattice[i][mu], g.local_lattice[i][mu] ) );
      ASSERT( DIVIDES( g.global_lattice[i][mu]/g.global_lattice[i+1][mu], g.local_lattice[i][mu] ) ); 
      ASSERT( DIVIDES( g.block_lattice[i][mu], g.global_lattice[i][mu]/g.global_lattice[i+1][mu] ) );
#ifdef SSE
      if ( ! g.block_lattice[i][mu] == g.global_lattice[i][mu]/g.global_lattice[i+1][mu] )
        warning0("when using SSE, Schwarz block size and aggregate size have to match.\n");
      ASSERT( g.block_lattice[i][mu] == g.global_lattice[i][mu]/g.global_lattice[i+1][mu] );
#endif
    }
    
  for ( i=0; i<g.num_levels-2; i++ )
    ASSERT( g.num_eig_vect[i] <= g.num_eig_vect[i+1] );
  
  if ( g.odd_even ) {
    int coarse_sites_per_core = 1;
    for ( mu=0; mu<4; mu++ ) {
      ASSERT( DIVIDES( 2, g.global_lattice[g.num_levels-1][mu] ) );
      coarse_sites_per_core *= g.local_lattice[g.num_levels-1][mu];
    }
    ASSERT( DIVIDES( 2, coarse_sites_per_core ) );
  }

  if ( g.method == 2 ) {
    for ( i=0; i<ls-1; i++ ) {
      int num_blocks = 1;
      for ( mu=0; mu<4; mu++) {
        num_blocks *= ( g.local_lattice[i][mu]/g.block_lattice[i][mu] );
      }
      ASSERT( num_blocks >= 2 );
    }
  }
  
  for ( mu=0; mu<4; mu++ )
    ASSERT( IMPLIES( g.rhs == 3, ASCENDING( 0, g.propagator_coords[mu], l->global_lattice[mu]-1 ) ) );
  
  ASSERT( IMPLIES( g.method > 0 && g.interpolation > 0, g.coarse_iter > 0 ) );
  ASSERT( IMPLIES( g.method > 0 && g.interpolation > 0, g.coarse_restart > 0 ) );
  ASSERT( IMPLIES( g.method > 0 && g.interpolation > 0, 0 < g.coarse_tol && g.coarse_tol < 1 ) );
  ASSERT( IMPLIES( g.method > 0, l->n_cy > 0 ) );
  ASSERT( g.max_restart > 0 );
  ASSERT( 0 < g.tol && g.tol < 1 );
  ASSERT( ASCENDING( 0, g.kcycle, 1 ) );
  ASSERT( IMPLIES( g.kcycle && g.method > 0, g.kcycle_restart > 0 ) );
  ASSERT( IMPLIES( g.kcycle && g.method > 0, g.kcycle_max_restart > 0 ) );
  ASSERT( IMPLIES( g.kcycle && g.method > 0, 0 < g.kcycle_tol && g.kcycle_tol < 1 ) );
  
#ifdef SSE
  ASSERT( g.mixed_precision );
//   ASSERT( DIVIDES( 4, g.num_eig_vect[0] ) );
#endif
}

void allocate_for_global_struct_after_read_global_info( int ls ) {

  int i;
  MALLOC( g.global_lattice, int*, ls );
  MALLOC( g.local_lattice, int*, ls );
  MALLOC( g.block_lattice, int*, ls );
  g.global_lattice[0] = NULL;
  g.local_lattice[0] = NULL;
  g.block_lattice[0] = NULL;
  MALLOC( g.global_lattice[0], int, 4*ls );
  MALLOC( g.local_lattice[0], int, 4*ls );
  MALLOC( g.block_lattice[0], int, 4*ls );
  MALLOC( g.post_smooth_iter, int, ls );
  MALLOC( g.ncycle, int, ls );
  MALLOC( g.relax_fac, double, ls );
  MALLOC( g.block_iter, int, ls );
  MALLOC( g.setup_iter, int, ls );
  MALLOC( g.num_eig_vect, int, ls );
  
  for ( i=1; i<ls; i++ ) {
    g.global_lattice[i] = g.global_lattice[0] + i*4;
    g.local_lattice[i] = g.local_lattice[0] + i*4;
    g.block_lattice[i] = g.block_lattice[0] + i*4;
  }
}

void set_level_and_global_structs_according_to_global_struct( level_struct *l ) {

  int mu;
  
  l->level = g.num_levels-1; l->depth = 0; l->idle = 0;
  l->global_lattice = g.global_lattice[0];
  l->local_lattice = g.local_lattice[0];
  l->block_lattice = g.block_lattice[0];
  l->post_smooth_iter = g.post_smooth_iter[0];
  l->n_cy = g.ncycle[0];
  l->relax_fac = g.relax_fac[0];
  l->block_iter = g.block_iter[0];
  l->setup_iter = g.setup_iter[0];
  l->num_eig_vect = g.num_eig_vect[0];
  
  // compute some additional values
  l->num_lattice_site_var = 12;
  g.num_processes = 1;
  for ( mu=0; mu<4; mu++ ) {
    l->comm_offset[mu] = 1;
    l->coarsening[mu] = g.global_lattice[0][mu]/MAX(1,g.global_lattice[1][mu]);
    l->global_splitting[mu] = l->global_lattice[mu]/l->local_lattice[mu];
    g.process_grid[mu] = l->global_splitting[mu];
    l->periodic_bc[mu] = 1;
    g.num_processes *= l->global_splitting[mu];
  }
  
  l->dirac_shift = l->real_shift;
  l->even_shift = l->dirac_shift;
  l->odd_shift = l->dirac_shift;
  g.solve_m0 = l->dirac_shift;
  g.setup_m0 = l->dirac_shift;
}

void lg_in( char *inputfile, level_struct *l ) {

  FILE *in;

  ASSERT( (in = fopen( inputfile, "r" )) != NULL );

  set_default_global_info();
  read_global_info( in );

  int ls = MAX(g.num_levels,2);
  allocate_for_global_struct_after_read_global_info( ls );

  read_solver_parameters( in, l );
  read_geometry_data( in, ls );
  
  ls = MAX(g.num_levels,2);
  
  set_level_and_global_structs_according_to_global_struct( l );

  read_testvector_io_data_if_necessary( in );
  read_evaluation_parameters_if_necessary( in );
  read_kcycle_data( in );
  
  validate_parameters( ls, l );

  printf00("configuration: %s\n", g.in );
  if( g.rhs == 4 )
    printf00("source list: %s\n", g.source_list );
  fclose(in);
}

void update_parameters_in_global_struct( const struct dd_alpha_amg_parameters *amg_params ) {
  // changes only parameters that can also be changed after initial setup
  int ls = MAX(g.num_levels, 2);
  for ( int i=0; i<ls; i++ ) {
    g.post_smooth_iter[i] = amg_params->post_smooth_iterations[i];
    g.block_iter[i]       = amg_params->post_smooth_block_iterations[i];
    g.setup_iter[i]       = amg_params->setup_iterations[i];
  }
  g.mass_for_next_solve = amg_params->solver_mass;
}
void update_parameters_in_level_structs( level_struct *l ) {

  l->post_smooth_iter = g.post_smooth_iter[l->depth];
  l->block_iter = g.block_iter[l->depth];
  l->setup_iter = g.setup_iter[l->depth];

  // for some reason coarser levels are initialized late, so we cannot always recurse here
  if ( l->level > 0 && l->next_level != NULL )
    update_parameters_in_level_structs( l->next_level );
}

void set_dd_alpha_amg_parameters( struct dd_alpha_amg_parameters *amg_params, level_struct *l ) {

  g.amg_params = *amg_params;
  set_default_global_info();
  set_global_info( amg_params );
  int ls = MAX(g.num_levels,2);
  allocate_for_global_struct_after_read_global_info( ls );

  set_geometry_data( amg_params );
  update_parameters_in_global_struct( amg_params );
  set_solver_parameters( amg_params, l );

  set_level_and_global_structs_according_to_global_struct( l );

  set_kcycle_data();

  validate_parameters( ls, l );
}

void update_dd_alpha_amg_parameters( const struct dd_alpha_amg_parameters *amg_params, level_struct *l ) {
  update_parameters_in_global_struct( amg_params );
  update_parameters_in_level_structs( l );
}
