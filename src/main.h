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

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include <stdarg.h>

#ifdef JUROPA
#include <mkl.h>
#endif

#include "dd_alpha_amg_parameters.h"
#include "dd_alpha_amg_setup_status.h"

#ifndef MAIN_HEADER
  #define MAIN_HEADER
  
  #define STRINGLENGTH 500
  
  #define _FILE_OFFSET_BITS 64
  #define EPS_float 1E-6
  #define EPS_double 1E-14
  
  #define FOR2( e ) { e e }
  #define FOR3( e ) { e e e }
  #define FOR4( e ) { e e e e }
  #define FOR10( e ) { e e e e e  e e e e e }
  #define FOR20( e ) { e e e e e  e e e e e  e e e e e  e e e e e  }
  #define FOR40( e ) { e e e e e  e e e e e  e e e e e  e e e e e  e e e e e  e e e e e  e e e e e  e e e e e }
  #define FOR6( e )  { e e e  e e e }
  #define FOR12( e ) { e e e  e e e  e e e  e e e }
  #define FOR24( e ) { e e e  e e e  e e e  e e e  e e e  e e e  e e e  e e e }
  #define FOR36( e ) { FOR12( e ) FOR12( e ) FOR12( e ) }
  #define FOR42( e ) { FOR36( e ) FOR6( e ) }
  
  #define SQUARE( e ) (e)*(e)
  #define NORM_SQUARE_float( e ) SQUARE(crealf( e ))+SQUARE(cimagf( e ))
  #define NORM_SQUARE_double( e ) SQUARE(creal( e ))+SQUARE(cimag( e ))
  #define CSPLIT( e ) creal(e), cimag(e)
  
  #define MPI_double MPI_DOUBLE
  #define MPI_float MPI_FLOAT
  #define MPI_COMPLEX_double MPI_DOUBLE_COMPLEX
  #define MPI_COMPLEX_float MPI_COMPLEX
  #define I _Complex_I
  #define conj_double conj
  #define conj_float conjf
  #define cabs_double cabs
  #define cabs_float cabsf
  #define creal_double creal
  #define creal_float crealf
  #define cimag_double cimag
  #define cimag_float cimagf
  #define csqrt_double csqrt
  #define csqrt_float csqrtf
  #define cpow_double cpow
  #define cpow_float cpowf
  #define pow_double pow
  #define pow_float powf
  #define abs_float fabs
  #define abs_double abs
  
#ifdef SSE
  #define MALLOC( variable, kind, length ) do{ if ( variable != NULL ) { \
  printf0("malloc of \"%s\" failed: pointer is not NULL (%s:%d).\n", #variable, __FILE__, __LINE__ ); } \
  if ( (length) > 0 ) { variable = (kind*) memalign( 64, sizeof(kind) * (length) ); } \
  if ( variable == NULL && (length) > 0 ) { \
  error0("malloc of \"%s\" failed: no memory allocated (%s:%d), current memory used: %lf GB.\n", \
  #variable, __FILE__, __LINE__, g.cur_storage/1024.0 ); } \
  g.cur_storage += (sizeof(kind) * (length))/(1024.0*1024.0); \
  if ( g.cur_storage > g.max_storage ) g.max_storage = g.cur_storage; }while(0)
#else
  #define MALLOC( variable, kind, length ) do{ if ( variable != NULL ) { \
  printf0("malloc of \"%s\" failed: pointer is not NULL (%s:%d).\n", #variable, __FILE__, __LINE__ ); } \
  if ( (length) > 0 ) { variable = (kind*) malloc( sizeof(kind) * (length) ); } \
  if ( variable == NULL && (length) > 0 ) { \
  error0("malloc of \"%s\" failed: no memory allocated (%s:%d), current memory used: %lf GB.\n", \
  #variable, __FILE__, __LINE__, g.cur_storage/1024.0 ); } \
  g.cur_storage += (sizeof(kind) * (length))/(1024.0*1024.0); \
  if ( g.cur_storage > g.max_storage ) g.max_storage = g.cur_storage; }while(0)
#endif

  // allocate and deallocate macros (hugepages, aligned)
  #include <fcntl.h>
  #include <sys/mman.h>
  #define HUGE_PAGE_SIZE (2 * 1024 * 1024)
  #define ROUND_UP_TO_FULL_PAGE(x) \
    (((x) + HUGE_PAGE_SIZE - 1) / HUGE_PAGE_SIZE * HUGE_PAGE_SIZE)
  
  
//   void *tmppointer = (void*)(variable);
//   posix_memalign( &tmppointer, alignment, sizeof(kind)*(length));
//   variable = (kind*)tmppointer; }
  
  #define MALLOC_HUGEPAGES( variable, kind, length, alignment ) do { if ( variable != NULL ) { \
  printf0("malloc of \"%s\" failed: pointer is not NULL (%s:%d).\n", #variable, __FILE__, __LINE__ ); } \
  if ( (length) > 0 ) { \
  variable = (kind*)memalign( alignment, sizeof(kind)*((size_t)length)); } \
  if ( variable == NULL && (length) > 0 ) { \
  error0("malloc of \"%s\" failed: no memory allocated (%s:%d), current memory used: %lf GB.\n", \
  #variable, __FILE__, __LINE__, g.cur_storage/1024.0 ); } \
  g.cur_storage += (sizeof(kind) * (length))/(1024.0*1024.0); \
  if ( g.cur_storage > g.max_storage ) g.max_storage = g.cur_storage; }while(0)
  
  #define FREE_HUGEPAGES( addr, kind, length ) do { free( addr ); addr = NULL; \
  g.cur_storage -= (sizeof(kind) * (length))/(1024.0*1024.0); }while(0)
  
  #ifdef DEBUG
    #define DPRINTF0 printf0
  #else
    #define DPRINTF0( ARGS, ... )
  #endif
  
  #define FREE( variable, kind, length ) do{ if ( variable != NULL ) { \
  free( variable ); variable = NULL; g.cur_storage -= (sizeof(kind) * (length))/(1024.0*1024.0); } else { \
  printf0("multiple free of \"%s\"? pointer is already NULL (%s:%d).\n", #variable, __FILE__, __LINE__ ); } }while(0)
  
  #define PUBLIC_MALLOC( variable, kind, size ) do{ START_MASTER(threading) MALLOC( variable, kind, size ); \
  ((kind**)threading->workspace)[0] = variable; END_MASTER(threading) SYNC_MASTER_TO_ALL(threading) \
  variable = ((kind**)threading->workspace)[0]; SYNC_MASTER_TO_ALL(threading) }while(0)
  
  #define PUBLIC_FREE( variable, kind, size ) do{ SYNC_MASTER_TO_ALL(threading) \
  START_MASTER(threading) FREE( variable, kind, size ); END_MASTER(threading) SYNC_MASTER_TO_ALL(threading) variable = NULL; }while(0)
  
  #define PUBLIC_MALLOC2( variable, kind, size, thrdng ) do{ START_MASTER(thrdng) MALLOC( variable, kind, size ); \
  ((kind**)thrdng->workspace)[0] = variable; END_MASTER(thrdng) SYNC_MASTER_TO_ALL(thrdng) \
  variable = ((kind**)thrdng->workspace)[0]; SYNC_MASTER_TO_ALL(thrdng) }while(0)
  
  #define PUBLIC_FREE2( variable, kind, size, thrdng ) do{ SYNC_MASTER_TO_ALL(thrdng) \
  START_MASTER(thrdng) FREE( variable, kind, size ); END_MASTER(thrdng) SYNC_MASTER_TO_ALL(thrdng) variable = NULL; }while(0)
  
  
  #define ASSERT( expression ) do{ if ( !(expression) ) { \
  error0("assertion \"%s\" failed (%s:%d)\n       bad choice of input parameters (please read the user manual in /doc).\n", \
  #expression, __FILE__, __LINE__ ); } }while(0)
  
  #define IMPLIES( A, B ) !( A ) || ( B )
  #define XOR( A, B ) (( A ) && !( B )) || (!( A ) && ( B ))
  #define NAND( A, B ) !( (A) && (B) )
  #define DIVIDES( A, B ) A == 0 || ((double)(B)/(double)(A) - (double)((int)(B)/(int)(A))) == 0 
  #define ASCENDING( A, B, C ) ( (A)<=(B) ) && ( (B)<=(C) )
  #define MAX( A, B ) ( (A > B) ? A : B )
  #define MIN( A, B ) ( (A < B) ? A : B )
  
  #ifdef DEBUG
  #define DEBUGOUTPUT_ARRAY( A, FORMAT, INDEX ) do{ \
  char TMPSTRING[100]; sprintf( TMPSTRING, "%s[%d] = %s\n", #A, INDEX, FORMAT ); \
  printf0( TMPSTRING, A[INDEX] ); }while(0)
  #else
  #define DEBUGOUTPUT_ARRAY( A, FORMAT, INDEX )
  #endif
  
  #ifdef DEBUG
  #define DEBUGOUTPUT( A, FORMAT ) do{ \
  char TMPSTRING[100]; sprintf( TMPSTRING, "%s = %s\n", #A, FORMAT ); \
  printf0( TMPSTRING, A ); }while(0)
  #else
  #define DEBUGOUTPUT( A, FORMAT )
  #endif

  #include "vectorization_control.h"
  #include "threading.h"

  // enumerations
  enum { _EVEN, _ODD };
  enum { _NO_DEFAULT_SET, _DEFAULT_SET };
  enum { _NO_REORDERING, _REORDER };
  enum { _ADD, _COPY };
  enum { _ORDINARY, _SCHWARZ };
  enum { _RES, _NO_RES };
  enum { _STANDARD, _LIME, _MULTI}; //formats
  enum { _READ, _WRITE };
  enum { _NO_SHIFT };
  enum { _BTWN_ORTH = 20 };
  enum { _GLOBAL_FGMRES, _K_CYCLE, _COARSE_GMRES, _SMOOTHER };
  enum { _COARSE_GLOBAL };
  enum { _FULL_SYSTEM, _EVEN_SITES, _ODD_SITES };
  enum { _LEFT, _RIGHT, _NOTHING };
  enum { _GIP, _PIP, _LA2, _LA6, _LA8, _LA, _CPY, _SET, _PR, _SC, _NC, _SM, _OP_COMM, _OP_IDLE, _ALLR, _GD_COMM, _GD_IDLE, _GRAM_SCHMIDT, _GRAM_SCHMIDT_ON_AGGREGATES,
      _SM1, _SM2, _SM3, _SM4, _SMALL1, _SMALL2, _NUM_PROF }; // _NUM_PROF has always to be the last constant!
  enum { _VTS = 20 };
  enum { _TRCKD_VAL, _STP_TIME, _SLV_ITER, _SLV_TIME, _CRS_ITER, _CRS_TIME, _SLV_ERR, _CGNR_ERR, _NUM_OPTB };
  
  typedef struct block_struct {
    int start, color, no_comm, *bt;
  } block_struct;
  
  #include "main_pre_def_float.h"
  #include "main_pre_def_double.h"
  
  extern complex_double _COMPLEX_double_ONE;
  extern complex_double _COMPLEX_double_ZERO;
  extern complex_double _COMPLEX_double_MINUS_ONE;
  extern complex_float  _COMPLEX_float_ONE;
  extern complex_float  _COMPLEX_float_ZERO;
  extern complex_float  _COMPLEX_float_MINUS_ONE;
  
  typedef struct plot_table_line {
    
    double values[_NUM_OPTB];
    struct plot_table_line *next;
    
  } plot_table_line;
  
  typedef struct var_table_entry {
    
    void *pt;
    char name[STRINGLENGTH];
    char datatype[20];
    struct var_table_entry *next;
    
  } var_table_entry;
  
  typedef struct var_table {
    
    int evaluation, multiplicative, shift_update, re_setup,
        track_error, track_cgn_error, average_over;
    char scan_var[STRINGLENGTH];
    double start_val, end_val, step_size, *output_table[6];
    var_table_entry *entry, *iterator;
    plot_table_line *p, *p_end;
    
  } var_table;
  
  typedef struct confbuffer_struct {
    
    double *data;
    struct confbuffer_struct *next;
    
  } confbuffer_struct;
  
  typedef struct {
    
    gmres_float_struct sp;
    gmres_double_struct dp;
    
  } gmres_MP_struct;

  typedef struct level_struct {    
    
    // distributed: non-idling processes of previos level
    // gathered: non-idling processes of current level
    
    // distributed
    operator_double_struct op_double;
    operator_float_struct op_float;
    // odd_even
    operator_double_struct oe_op_double;
    operator_float_struct oe_op_float;
    // gathered / schwarz
    schwarz_double_struct s_double;
    schwarz_float_struct s_float;
    // interpolation / aggregation
    interpolation_double_struct is_double;
    interpolation_float_struct is_float;
    // gathering parameters and buffers
    gathering_double_struct gs_double;
    gathering_float_struct gs_float;
    // k cycle
    gmres_float_struct p_float;
    gmres_double_struct p_double;
    // gmres as a smoother
    gmres_float_struct sp_float;
    gmres_double_struct sp_double;
    // dummy gmres struct
    gmres_float_struct dummy_p_float;
    gmres_double_struct dummy_p_double;
    //profiling
    profiling_float_struct prof_float;
    profiling_double_struct prof_double;
    
    // communication
    MPI_Request *reqs;
    int parent_rank, idle, neighbor_rank[8], num_processes, num_processes_dir[4];
    // lattice
    int *global_lattice;
    int *local_lattice;
    int *block_lattice;
    int num_eig_vect;
    int coarsening[4];
    int global_splitting[4];
    int periodic_bc[4];
    int comm_offset[4];
    // degrees of freedom on a site
    // 12 on fine lattice (i.e., complex d.o.f.)
    // 2*num_eig_vect on coarser lattices
    int num_lattice_site_var;
    int level;
    int depth;
    // number of sites in local volume + ghost shell (either fw or bw)
    int num_lattice_sites;
    // number of sites in local volume
    int num_inner_lattice_sites;
    int num_boundary_sites[4];
    // complex d.o.f. in local volume + ghost shell = num_lattice_sites * num_lattice_site_var
    int vector_size;
    // complex d.o.f. in local volume = num_inner_lattice_sites * num_lattice_site_var
    int inner_vector_size;
    int schwarz_vector_size;
    int D_size;
    int clover_size;
    // operator
    double real_shift;
    complex_double dirac_shift, even_shift, odd_shift;
    // buffer vectors
    vector_float vbuf_float[9], sbuf_float[2];
    vector_double vbuf_double[9], sbuf_double[2];
    // storage + daggered-operator bufferes
    vector_double x;
    // local solver parameters
    double tol, relax_fac;
    int n_cy, post_smooth_iter, block_iter, setup_iter;
    
    // next coarser level
    struct level_struct *next_level;
    
  } level_struct;


  typedef struct global_struct {
    
    FILE *logfile;
    
    gmres_double_struct p;
    gmres_MP_struct p_MP;
    operator_double_struct op_double;
    operator_float_struct op_float;

    // communication
    MPI_Comm comm_cart;
    MPI_Group global_comm_group;
    MPI_Request sreqs[8], rreqs[8];
    int num_processes, my_rank, my_coords[4], two_cnfgs, tv_io_single_file, num_openmp_processes;
    // string buffers
    char in[STRINGLENGTH], in_clov[STRINGLENGTH], source_list[STRINGLENGTH], tv_io_file_name[STRINGLENGTH];
    // geometry, method parameters
    int num_levels, num_desired_levels, process_grid[4], in_format,
        **global_lattice, **local_lattice, **block_lattice, 
        *post_smooth_iter, *block_iter, *setup_iter, *ncycle,
        method, odd_even, anti_pbc, rhs, propagator_coords[4],
        interpolation, randomize, *num_eig_vect, num_coarse_eig_vect, kcycle, mixed_precision,
        restart, max_restart, kcycle_restart, kcycle_max_restart, coarse_iter, coarse_restart;
    double tol, coarse_tol, kcycle_tol, csw, rho, *relax_fac;

    // profiling, analysis, output
    int coarse_iter_count, iter_count, iterator, print, conf_flag, setup_flag, in_setup;
    double coarse_time, prec_time, *output_table[8], cur_storage, max_storage, total_time,
           plaq_hopp, plaq_clov, norm_res, plaq, setup_m0, solve_m0, bicgstab_tol;
           
    // index functions for external usage
    int (*conf_index_fct)(), (*vector_index_fct)();
    int *odd_even_table;
    
    // bc: 0 dirichlet, 1 periodic, 2 anti-periodic
    int bc; 
    
    complex_double **gamma, g5D_shift;
    var_table vt;

    struct dd_alpha_amg_parameters amg_params;
    struct dd_alpha_amg_setup_status mg_setup_status;
    double mass_for_next_solve;
    
  } global_struct;

  extern global_struct g;
  
  static inline void printf0( char* format, ... ) {
    START_MASTER(no_threading)
    if ( g.my_rank == 0 && g.print >= 0 ) {
      va_list argpt;
      va_start(argpt,format);
      vprintf(format,argpt);
#ifdef WRITE_LOGFILE
      vfprintf(g.logfile,format,argpt);
      fflush(g.logfile);
#endif
      va_end(argpt);
      fflush(0);
    }
    END_MASTER(no_threading)
  }
  
  static inline void warning0( char* format, ... ) {
    if ( g.my_rank == 0 && g.print >= 0 ) {
      printf("\x1b[31mwarning: ");
      va_list argpt;
      va_start(argpt,format);
      vprintf(format,argpt);
#ifdef WRITE_LOGFILE
      vfprintf(g.logfile,format,argpt);
      fflush(g.logfile);
#endif
      va_end(argpt);
      printf("\x1b[0m");
      fflush(0);
    }
  }
  
  static inline void error0( char* format, ... ) {
    if ( g.my_rank == 0 ) {
      printf("\x1b[31merror: ");
      va_list argpt;
      va_start(argpt,format);
      vprintf(format,argpt);
#ifdef WRITE_LOGFILE
      vfprintf(g.logfile,format,argpt);
      fflush(g.logfile);
#endif
      va_end(argpt);
      printf("\x1b[0m");
      fflush(0);
      MPI_Abort( MPI_COMM_WORLD, 0 );
    }
  }
  
  static inline void printf00( char* format, ... ) {
    if ( g.my_rank == 0 && g.print >= 0 ) {
      va_list argpt;
      va_start(argpt,format);
      vprintf(format,argpt);
#ifdef WRITE_LOGFILE
      vfprintf(g.logfile,format,argpt);
      fflush(g.logfile);
#endif
      va_end(argpt);
      fflush(0);
    }
  }
  
#endif

// functions
#include "clifford.h"

#ifdef SSE
#include "vectorization_dirac_float.h"
#include "vectorization_dirac_double.h"
#include "blas_vectorized.h"
#include "sse_blas_vectorized.h"
#include "sse_complex_float_intrinsic.h"
#include "sse_complex_double_intrinsic.h"
#include "sse_coarse_operator_float.h"
#include "sse_coarse_operator_double.h"
#include "sse_linalg_float.h"
#include "sse_linalg_double.h"
#include "sse_interpolation_float.h"
#include "sse_interpolation_double.h"
#include "sse_schwarz_float.h"
#include "sse_schwarz_double.h"
#else
//no intrinsics
#include "interpolation_float.h"
#include "interpolation_double.h"
#endif

#include "data_float.h"
#include "data_double.h"
#include "data_layout.h"
#include "io.h"
#include "init.h"
#include "operator_float.h"
#include "operator_double.h"
#include "dirac.h"
#include "dirac_float.h"
#include "dirac_double.h"
#include "oddeven_float.h"
#include "oddeven_double.h"
#include "linalg.h"
#include "linalg_float.h"
#include "linalg_double.h"
#include "ghost_float.h"
#include "ghost_double.h"
#include "linsolve_float.h"
#include "linsolve_double.h"
#include "linsolve.h"
#include "preconditioner.h"
#include "vcycle_float.h"
#include "vcycle_double.h"
#include "solver_analysis.h"
#include "top_level.h"
#include "ghost.h"
#include "init_float.h"
#include "init_double.h"
#include "schwarz_double.h"
#include "schwarz_float.h"
#include "setup_float.h"
#include "setup_double.h"
#include "coarsening_float.h"
#include "coarsening_double.h"
#include "gathering_float.h"
#include "gathering_double.h"
#include "coarse_operator_float.h"
#include "coarse_operator_double.h"
#include "coarse_oddeven_float.h"
#include "coarse_oddeven_double.h"
#include "var_table.h"
#include "main_post_def_float.h"
#include "main_post_def_double.h"
#ifdef HAVE_LIME
#include <lime.h>
#include <lime_config.h>
#include <dcap-overload.h>
#include <lime_defs.h>
#include <lime_header.h>
#include <lime_writer.h>
#include <lime_reader.h>
#endif
#include "lime_io.h"
