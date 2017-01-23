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

void coarse_selfcoupling_LU_decomposition_PRECISION( const config_PRECISION output, config_PRECISION input, level_struct *l ) {
  // input = [ A B      , A=A*, D=D*, C = -B*
  //           C D ]
  //
  // order: upper triangle of A, upper triangle of D, B, each column major

  register int i, j, k, n = l->num_lattice_site_var/2, n2 = l->num_lattice_site_var;
  
  // set the matrix up
  // A
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<j; i++ ) {
      output[n2*i+j] = *input;
      output[i+n2*j] = conj_PRECISION(*input);
      input++;      
    }
    output[(n2+1)*j] = *input;
    input++; // diagonal entry
  }
  // D
  for ( j=n; j<n2; j++ ) {
    for ( i=n; i<j; i++ ) {
      output[n2*i+j] = *input;
      output[i+n2*j] = conj_PRECISION(*input);
      input++;      
    }
    output[(n2+1)*j] = *input;
    input++; // diagonal entry
  }
  // B
  for ( j=n; j<n2; j++ ) {
    for ( i=0; i<n; i++ ) {
      output[n2*i+j] = *input;
      output[i+n2*j] = -conj_PRECISION(*input);
      input++;      
    }
  }
  
  // compute LU decomposition
  // output = triu(L,1) + tril(U,0)
  // i.e., output contains L and U without the diagonal of L which is equal to 1
  // order: row major
  for ( k=0; k<n2; k++ ) {
    for ( i=k+1; i<n2; i++ ) {
      output[n2*i+k] = output[n2*i+k]/output[(n2+1)*k]; // output(i,k) = output(i,k)/output(k,k)
      for ( j=k+1; j<n2; j++ )
        output[n2*i+j] = output[n2*i+j]-output[n2*i+k]*output[n2*k+j]; // output(i,j) = output(i,j)-output(i,k)*output(k,j)
    }
  }
}


void coarse_perform_fwd_bwd_subs_PRECISION( vector_PRECISION x, vector_PRECISION b, config_PRECISION A, level_struct *l ) {
  
  register int i, j, n2 = l->num_lattice_site_var;
  
  // solve x = U^(-1) L^(-1) b
  // forward substitution with L
  for ( i=0; i<n2; i++ ) {
    x[i] = b[i];
    for ( j=0; j<i; j++ ) {
      x[i] = x[i] - A[i*n2+j]*x[j];
    }
  }
  // backward substitution with U
  for ( i=n2-1; i>=0; i-- ) {
    for ( j=i+1; j<n2; j++ ) {
      x[i] = x[i] - A[i*n2+j]*x[j];
    }
    x[i] = x[i]/A[i*(n2+1)];
  }
}


void coarse_LU_multiply_PRECISION( vector_PRECISION y, vector_PRECISION x, config_PRECISION A, level_struct *l ) {
  
  register int i, j, n2 = l->num_lattice_site_var;
  
  // y = Ax
  // multiplication with U
  for ( i=0; i<n2; i++ ) {
    y[i] = A[i*(n2+1)]*x[i];
    for ( j=i+1; j<n2; j++ )
      y[i] += A[i*n2+j]*x[j];
  }
  // multiplication with L
  for ( i=n2-1; i>0; i-- )
    for ( j=0; j<i; j++ )
      y[i] += A[i*n2+j]*y[j];
}


void coarse_diag_PRECISION( vector_PRECISION y, vector_PRECISION x, operator_PRECISION_struct *op, level_struct *l ) {
  
  coarse_diag_ee_PRECISION( y, x, op, l, no_threading );
  coarse_diag_oo_PRECISION( y, x, op, l, no_threading );
}


void coarse_diag_ee_PRECISION( vector_PRECISION y, vector_PRECISION x, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {
  
  int n1 = op->num_even_sites;
  int start;
  int end;
  compute_core_start_end_custom(0, n1, &start, &end, l, threading, 1);
  // even sites
#ifndef VECTORIZE_COARSE_OPERATOR_PRECISION
  int offset = l->num_lattice_site_var;
  coarse_self_couplings_PRECISION( y+start*offset, x+start*offset, op->clover+start*(offset*offset+offset)/2, (end-start)*offset, l );
#else
  coarse_self_couplings_PRECISION_vectorized( y, x, op->clover_vectorized, start, end, l );
#endif
}


void coarse_diag_oo_PRECISION( vector_PRECISION y, vector_PRECISION x, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {
  
  int n1 = op->num_even_sites, n2 = op->num_odd_sites,
      offset = l->num_lattice_site_var, ess = (l->num_lattice_site_var/2)*(l->num_lattice_site_var+1);
  config_PRECISION sc = op->clover;
  int start;
  int end;
  compute_core_start_end_custom(n1, n1+n2, &start, &end, l, threading, 1);
  
  x += start*offset;
  y += start*offset;
  sc += start*ess;
  
  // odd sites
#ifndef VECTORIZE_COARSE_OPERATOR_PRECISION
  int oss = l->num_lattice_site_var*l->num_lattice_site_var;
  for ( int i=start; i<end; i++ ) {
    coarse_LU_multiply_PRECISION( y, x, sc, l );
    x += offset;
    y += offset;
    sc += oss;
  }
#else
  // take care on last level:
  // - vectorized, but we have stored oo^{-1}, so we cannot use it
  // - when vectorization is used LU decomposition is not computed, so we also cannot use coarse_LU_multiply_PRECISION
  // => use standard non-vectorized multiplication
  if ( l->level == 0 )
    coarse_self_couplings_PRECISION( y, x, sc, (end-start)*offset, l );
  else
    coarse_self_couplings_PRECISION_vectorized( y-start*offset, x-start*offset, op->clover_vectorized, start, end, l );
#endif
}


void coarse_diag_oo_inv_PRECISION( vector_PRECISION y, vector_PRECISION x, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {
  
  int n1 = op->num_even_sites, n2 = op->num_odd_sites, start, end;

  compute_core_start_end_custom(n1, n1+n2, &start, &end, l, threading, 1);
  
  // odd sites
#ifndef VECTORIZE_COARSE_OPERATOR_PRECISION
  int offset = l->num_lattice_site_var, ess = (l->num_lattice_site_var/2)*(l->num_lattice_site_var+1),
      oss = l->num_lattice_site_var*l->num_lattice_site_var;
  config_PRECISION sc = op->clover;
  x += start*offset;
  y += start*offset;
  sc += n1*ess + (start-n1)*oss;
  
  for ( int i=start; i<end; i++ ) {
    coarse_perform_fwd_bwd_subs_PRECISION( y, x, sc, l );
    x += offset;
    y += offset;
    sc += oss;
  }
#else
  coarse_self_couplings_PRECISION_vectorized( y, x, op->clover_vectorized, start, end, l );
#endif
}

void coarse_oddeven_setup_PRECISION_set_couplings( operator_PRECISION_struct *in, int reorder, level_struct *l, struct Thread *threading ) {

  int i, j, n=l->num_inner_lattice_sites, sc_size = (l->num_lattice_site_var/2)*(l->num_lattice_site_var+1),
      nc_size = SQUARE(l->num_lattice_site_var),
      t, z, y, x;
  operator_PRECISION_struct *op = &(l->oe_op_PRECISION);
  config_PRECISION sc_in = in->clover, nc_in = in->D, Aee = NULL, Aoo = NULL;
  int *le = l->local_lattice;
  int oe_offset = op->oe_offset;
  
#ifndef VECTORIZE_COARSE_OPERATOR_PRECISION
  int lu_dec_size = SQUARE(l->num_lattice_site_var);
#endif  
  
  Aee = op->clover;
  Aoo = op->clover + op->num_even_sites*sc_size;
  
  START_LOCKED_MASTER(threading)
  // self coupling  
  if ( reorder ) {
    int k=0, index, *it = in->index_table, *dt = in->table_dim;
    j=0;
    for ( t=0; t<le[T]; t++ )
      for ( z=0; z<le[Z]; z++ )
        for ( y=0; y<le[Y]; y++ )
          for ( x=0; x<le[X]; x++ ) {
            index = site_index( t, z, y, x, dt, it );
            if ( (t+z+y+x+oe_offset)%2 == 1 ) {
#ifndef VECTORIZE_COARSE_OPERATOR_PRECISION
              coarse_selfcoupling_LU_decomposition_PRECISION( Aoo+j, sc_in+sc_size*index, l );
              j+=lu_dec_size;
#else
              for ( i=0; i<sc_size; i++, j++ )
                Aoo[j] = sc_in[ sc_size*index+i ];
#endif
            } else {
              for ( i=0; i<sc_size; i++, k++ )
                Aee[k] = sc_in[ sc_size*index+i ];
            }
          }
          
  } else {
    j = op->num_even_sites*sc_size;
    for ( i=0; i<j; i++ )
      Aee[i] = sc_in[i]; // even sites
    
#ifndef VECTORIZE_COARSE_OPERATOR_PRECISION
    sc_in += j;
    j = op->num_odd_sites;
    for ( i=0; i<j; i++ ) {
      coarse_selfcoupling_LU_decomposition_PRECISION( Aoo, sc_in, l ); // odd sites, compute LU decomposition
      sc_in += sc_size; Aoo += lu_dec_size;
    }
#else
    for ( i=op->num_even_sites*sc_size; i<n*sc_size; i++ )
      Aee[i] = sc_in[i]; // even sites
#endif
  }
  
  // neighbor couplings
  if ( reorder ) {
    int k=0, index, *it = in->index_table, *dt = in->table_dim, site_size=4*nc_size;
    config_PRECISION oAe=op->D, eAo=(op->D)+site_size*op->num_even_sites;
    j=0;
    for ( t=0; t<le[T]; t++ )
      for ( z=0; z<le[Z]; z++ )
        for ( y=0; y<le[Y]; y++ )
          for ( x=0; x<le[X]; x++ ) {
            index = site_index( t, z, y, x, dt, it );
            if ( (t+z+y+x+oe_offset)%2 == 1 ) {
              for ( i=0; i<site_size; i++, j++ ) {
                eAo[j] = nc_in[ site_size*index+i ];
              }
            } else {
              for ( i=0; i<site_size; i++, k++ ) {
                oAe[k] = nc_in[ site_size*index+i ];
              }
            }
          }
          
  } else {
    j = n*4*nc_size;
    for ( i=0; i<j; i++ )
      op->D[i] = nc_in[i];
  }
  END_LOCKED_MASTER(threading)

#ifdef VECTORIZE_COARSE_OPERATOR_PRECISION
  int start;
  int end;
  compute_core_start_end_custom(0, n, &start, &end, l, threading, 1);
  int n_per_core = end-start;
  int column_offset = SIMD_LENGTH_PRECISION*((l->num_lattice_site_var+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  int offset_v = 2*l->num_lattice_site_var*column_offset;
  copy_coarse_operator_to_vectorized_layout_PRECISION(
      op->D + 4*start*nc_size,
      op->D_vectorized + 4*start*offset_v,
      n_per_core, l->num_lattice_site_var/2);
  copy_coarse_operator_to_transformed_vectorized_layout_PRECISION(
      op->D + 4*start*nc_size,
      op->D_transformed_vectorized + 4*start*offset_v,
      n_per_core, l->num_lattice_site_var/2);
  copy_coarse_operator_clover_to_vectorized_layout_PRECISION(
      op->clover + start*sc_size,
      op->clover_vectorized + start*offset_v,
      n_per_core, l->num_lattice_site_var/2);
  SYNC_CORES(threading)

  compute_core_start_end_custom(op->num_even_sites, n, &start, &end, l, threading, 1);
  OPERATOR_TYPE_PRECISION tmp[offset_v] __attribute__((aligned(64)));
  for(int a=start; a<end; a++) {
    for(int i=0; i<offset_v; i++)
      tmp[i] = (op->clover_vectorized + a*offset_v)[i];
    cgem_inverse(l->num_lattice_site_var, op->clover_vectorized + a*offset_v, tmp, column_offset);
  }

  SYNC_CORES(threading)
#endif
}

void coarse_oddeven_setup_PRECISION( operator_PRECISION_struct *in, int reorder, level_struct *l ) {

  int n=l->num_inner_lattice_sites, oe_offset=0, mu, nu,
      lu_dec_size = SQUARE(l->num_lattice_site_var),
      nc_size = SQUARE(l->num_lattice_site_var), bs, **bt = NULL,
      *eot = NULL, *nt = NULL, *tt = NULL, t, z, y, x, le[4], N[4];
  operator_PRECISION_struct *op = &(l->oe_op_PRECISION);

  for ( mu=0; mu<4; mu++ ) {
    le[mu] = l->local_lattice[mu];
    N[mu] = le[mu]+1;
    op->table_dim[mu] = N[mu];
  }

  for ( mu=0; mu<4; mu++ )
    oe_offset += (l->local_lattice[mu]*(g.my_coords[mu]/l->comm_offset[mu]))%2;
  oe_offset = oe_offset%2;

  // estimate site numbers
  op->num_even_sites = 0;
  op->num_odd_sites = 0;
  op->oe_offset = oe_offset;
  for ( t=0; t<le[T]; t++ )
    for ( z=0; z<le[Z]; z++ )
      for ( y=0; y<le[Y]; y++ )
        for ( x=0; x<le[X]; x++ ) {
          if ( (t+z+y+x+oe_offset)%2 == 1 ) {
            op->num_odd_sites++;
          } else {
            op->num_even_sites++;
          }
        }

  MALLOC( op->D, complex_PRECISION, 4*nc_size*n );
  MALLOC( op->clover, complex_PRECISION, lu_dec_size*n );
#ifdef VECTORIZE_COARSE_OPERATOR_PRECISION
  int column_offset = SIMD_LENGTH_PRECISION*((l->num_lattice_site_var+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  // 2 is for complex, 4 is for 4 directions
  MALLOC_HUGEPAGES( op->D_vectorized, OPERATOR_TYPE_PRECISION, 2*4*l->num_lattice_site_var*column_offset*n, 64 );
  MALLOC_HUGEPAGES( op->D_transformed_vectorized, OPERATOR_TYPE_PRECISION, 2*4*l->num_lattice_site_var*column_offset*n, 64 );
  MALLOC_HUGEPAGES( op->clover_vectorized, OPERATOR_TYPE_PRECISION, 2*l->num_lattice_site_var*column_offset*n, 64 );
#endif

  coarse_oddeven_setup_PRECISION_set_couplings( in, reorder, l, no_threading );
    
  // define data layout
  MALLOC( op->index_table, int, N[T]*N[Z]*N[Y]*N[X] );
  eot = op->index_table;
  
  define_eot( eot, N, l );
    
  // neighbor table, translation table
  MALLOC( op->neighbor_table, int, 5*N[T]*N[Z]*N[Y]*N[X] );
  MALLOC( op->backward_neighbor_table, int, 5*N[T]*N[Z]*N[Y]*N[X] );
  MALLOC( op->translation_table, int, le[T]*le[Z]*le[Y]*le[X] );
  nt = op->neighbor_table;
  tt = op->translation_table;
  
  define_nt_bt_tt( nt, op->backward_neighbor_table, NULL, tt, eot, N, l );
  
  // boundary table
  for ( mu=0; mu<4; mu++ ) {
    bs = 1;
    le[mu] = 1;
    for ( nu=0; nu<4; nu++ )
      bs *= le[nu];
    
    MALLOC( op->c.boundary_table[2*mu], int, bs );
    op->c.boundary_table[2*mu+1] = op->c.boundary_table[2*mu];
    
    le[mu] = l->local_lattice[mu];
  }
  
  bt = op->c.boundary_table;
  define_eo_bt( bt, eot, op->c.num_even_boundary_sites, op->c.num_odd_boundary_sites, op->c.num_boundary_sites, N, l );

  MALLOC( op->buffer, complex_PRECISION*, 2 );
  op->buffer[0] = NULL;
  MALLOC( op->buffer[0], complex_PRECISION, 2*l->vector_size );
  op->buffer[1] = op->buffer[0] + l->vector_size;  
  ghost_alloc_PRECISION( 0, &(op->c), l );
  ghost_sendrecv_init_PRECISION( _COARSE_GLOBAL, &(op->c), l ) ;
  if ( l->level == 0 )
    l->p_PRECISION.v_end = op->num_even_sites*l->num_lattice_site_var;
  else
    l->sp_PRECISION.v_end = op->num_even_sites*l->num_lattice_site_var;
}

void coarse_oddeven_re_setup_PRECISION( operator_PRECISION_struct *in, int reorder, level_struct *l, struct Thread *threading ) {
  coarse_oddeven_setup_PRECISION_set_couplings( in, reorder, l, threading );
}


void coarse_oddeven_free_PRECISION( level_struct *l ) {
  
  int mu, nu, nc_size = SQUARE(l->num_lattice_site_var),
      *ll = l->local_lattice, n = l->num_inner_lattice_sites, bs;
  
  ghost_free_PRECISION( &(l->oe_op_PRECISION.c), l );
  FREE( l->oe_op_PRECISION.D, complex_PRECISION, 4*nc_size*n );
#ifdef VECTORIZE_COARSE_OPERATOR_PRECISION
  int column_offset = SIMD_LENGTH_PRECISION*((l->num_lattice_site_var+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  FREE_HUGEPAGES( l->oe_op_PRECISION.D_vectorized, OPERATOR_TYPE_PRECISION, 2*4*l->num_lattice_site_var*column_offset*n );
  FREE_HUGEPAGES( l->oe_op_PRECISION.D_transformed_vectorized, OPERATOR_TYPE_PRECISION, 2*4*l->num_lattice_site_var*column_offset*n );
  FREE_HUGEPAGES( l->oe_op_PRECISION.clover_vectorized, OPERATOR_TYPE_PRECISION, 2*l->num_lattice_site_var*column_offset*n );
#endif
  FREE( l->oe_op_PRECISION.clover, complex_PRECISION, nc_size*n );
  FREE( l->oe_op_PRECISION.index_table, int, (ll[T]+1)*(ll[Z]+1)*(ll[Y]+1)*(ll[X]+1) );
  FREE( l->oe_op_PRECISION.neighbor_table, int, 5*(ll[T]+1)*(ll[Z]+1)*(ll[Y]+1)*(ll[X]+1) );
  FREE( l->oe_op_PRECISION.backward_neighbor_table, int, 5*(ll[T]+1)*(ll[Z]+1)*(ll[Y]+1)*(ll[X]+1) );
  FREE( l->oe_op_PRECISION.translation_table, int, ll[T]*ll[Z]*ll[Y]*ll[X] );
  
  for ( mu=0; mu<4; mu++ ) {
    bs = 1;
    for ( nu=0; nu<4; nu++ )
      if ( mu != nu )
        bs *= ll[nu];
      
      FREE( l->oe_op_PRECISION.c.boundary_table[2*mu], int, bs );
    l->oe_op_PRECISION.c.boundary_table[2*mu+1] = NULL;
  }
  
  FREE( l->oe_op_PRECISION.buffer[0], complex_PRECISION, 2*l->vector_size );
  FREE( l->oe_op_PRECISION.buffer, complex_PRECISION*, 2 );
}


void coarse_hopping_term_PRECISION( vector_PRECISION out, vector_PRECISION in, operator_PRECISION_struct *op,
                                    const int amount, level_struct *l, struct Thread *threading ) {

  START_NO_HYPERTHREADS(threading)

  int mu, i, index, num_site_var=l->num_lattice_site_var,
      num_4link_var=4*l->num_lattice_site_var*l->num_lattice_site_var,
      num_link_var=l->num_lattice_site_var*l->num_lattice_site_var,
      start=0, num_lattice_sites=l->num_inner_lattice_sites,
      plus_dir_param=_FULL_SYSTEM, minus_dir_param=_FULL_SYSTEM;
  vector_PRECISION in_pt, out_pt;
  config_PRECISION D_pt;

  int core_start;
  int core_end;
  
  // assumptions (1) self coupling has already been performed
  //          OR (2) "out" is initialized with zeros
  set_boundary_PRECISION( out, 0, l, threading );
  
  if ( amount == _EVEN_SITES ) {
    minus_dir_param = _ODD_SITES;
    plus_dir_param = _EVEN_SITES;
  } else if ( amount == _ODD_SITES ) {
    minus_dir_param = _EVEN_SITES;
    plus_dir_param = _ODD_SITES;
  }
  
  START_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // communicate in -mu direction
      ghost_sendrecv_PRECISION( in, mu, -1, &(op->c), minus_dir_param, l );
    }
  }
  END_MASTER(threading)
  SYNC_CORES(threading)
  
  if ( amount == _EVEN_SITES ) {
    start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
  } else if ( amount == _ODD_SITES ) {
    start = 0; num_lattice_sites = op->num_even_sites;
  }
  compute_core_start_end_custom(start, start+num_lattice_sites, &core_start, &core_end, l, threading, 1);
  
  // compute U_mu^dagger coupling
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 0*num_link_var;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+T];
    coarse_daggered_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 1*num_link_var;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+Z];
    coarse_daggered_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 2*num_link_var;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+Y];
    coarse_daggered_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 3*num_link_var;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+X];
    coarse_daggered_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }
  
  START_LOCKED_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // communicate in +mu direction
      ghost_sendrecv_PRECISION( out, mu, +1, &(op->c), plus_dir_param, l );    
    }
    for ( mu=0; mu<4; mu++ ) {
      // wait for -mu direction
      ghost_wait_PRECISION( in, mu, -1, &(op->c), minus_dir_param, l );    
    }
  }
  END_LOCKED_MASTER(threading)
  
  if ( amount == _EVEN_SITES ) {
    start = 0; num_lattice_sites = op->num_even_sites;
  } else if ( amount == _ODD_SITES ) {
    start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
  }
  compute_core_start_end_custom(start, start+num_lattice_sites, &core_start, &core_end, l, threading, 1);
  
  // compute U_mu couplings
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    out_pt = out + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index];
    index++;
    in_pt = in + num_site_var*op->neighbor_table[index+T];
    coarse_hopp_PRECISION( out_pt, in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt = in + num_site_var*op->neighbor_table[index+Z];
    coarse_hopp_PRECISION( out_pt, in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt = in + num_site_var*op->neighbor_table[index+Y];
    coarse_hopp_PRECISION( out_pt, in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt = in + num_site_var*op->neighbor_table[index+X];
    coarse_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }
  
  START_LOCKED_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // wait for +mu direction
      ghost_wait_PRECISION( out, mu, +1, &(op->c), plus_dir_param, l );
    }
  }
  END_LOCKED_MASTER(threading)

  END_NO_HYPERTHREADS(threading)
}


void coarse_n_hopping_term_PRECISION( vector_PRECISION out, vector_PRECISION in, operator_PRECISION_struct *op,
                                      const int amount, level_struct *l, struct Thread *threading ) {

#ifdef VECTORIZE_COARSE_OPERATOR_PRECISION
#ifndef COMM_HIDING_COARSEOP
  int sign = -1;
  coarse_pn_hopping_term_PRECISION_vectorized( out, in, op, amount, l, sign, threading);
#else
  coarse_n_hopping_term_PRECISION_vectorized( out, in, op, amount, l, threading );
#endif
  return;
#else
  START_NO_HYPERTHREADS(threading)

  int mu, i, index, num_site_var=l->num_lattice_site_var,
      num_4link_var=4*l->num_lattice_site_var*l->num_lattice_site_var,
      num_link_var=l->num_lattice_site_var*l->num_lattice_site_var,
      start=0, num_lattice_sites=l->num_inner_lattice_sites,
      plus_dir_param=_FULL_SYSTEM, minus_dir_param=_FULL_SYSTEM;
  vector_PRECISION in_pt, out_pt;
  config_PRECISION D_pt;

  int core_start;
  int core_end;
  
  // assumptions (1) self coupling has already been performed
  //          OR (2) "out" is initialized with zeros
  set_boundary_PRECISION( out, 0, l, threading );
  
  if ( amount == _EVEN_SITES ) {
    minus_dir_param = _ODD_SITES;
    plus_dir_param = _EVEN_SITES;
  } else if ( amount == _ODD_SITES ) {
    minus_dir_param = _EVEN_SITES;
    plus_dir_param = _ODD_SITES;
  }
  
  START_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // communicate in -mu direction
      ghost_sendrecv_PRECISION( in, mu, -1, &(op->c), minus_dir_param, l );
    }
  }
  END_MASTER(threading)
  SYNC_CORES(threading)
  
  if ( amount == _EVEN_SITES ) {
    start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
  } else if ( amount == _ODD_SITES ) {
    start = 0; num_lattice_sites = op->num_even_sites;
  }
  compute_core_start_end_custom(start, start+num_lattice_sites, &core_start, &core_end, l, threading, 1);
  
  // compute U_mu^dagger coupling
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 0*num_link_var;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+T];
    coarse_n_daggered_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 1*num_link_var;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+Z];
    coarse_n_daggered_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 2*num_link_var;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+Y];
    coarse_n_daggered_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index] + 3*num_link_var;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+X];
    coarse_n_daggered_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }
  
  START_LOCKED_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // communicate in +mu direction
      ghost_sendrecv_PRECISION( out, mu, +1, &(op->c), plus_dir_param, l );    
    }
    for ( mu=0; mu<4; mu++ ) {
      // wait for -mu direction
      ghost_wait_PRECISION( in, mu, -1, &(op->c), minus_dir_param, l );    
    }
  }
  END_LOCKED_MASTER(threading)
  
  if ( amount == _EVEN_SITES ) {
    start = 0; num_lattice_sites = op->num_even_sites;
  } else if ( amount == _ODD_SITES ) {
    start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
  }
  compute_core_start_end_custom(start, start+num_lattice_sites, &core_start, &core_end, l, threading, 1);
  
  // compute U_mu couplings
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    out_pt = out + num_site_var*op->neighbor_table[index];
    D_pt = op->D + num_4link_var*op->neighbor_table[index];
    index++;
    in_pt = in + num_site_var*op->neighbor_table[index+T];
    coarse_n_hopp_PRECISION( out_pt, in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt = in + num_site_var*op->neighbor_table[index+Z];
    coarse_n_hopp_PRECISION( out_pt, in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt = in + num_site_var*op->neighbor_table[index+Y];
    coarse_n_hopp_PRECISION( out_pt, in_pt, D_pt, l );
    
    D_pt += num_link_var;
    in_pt = in + num_site_var*op->neighbor_table[index+X];
    coarse_n_hopp_PRECISION( out_pt, in_pt, D_pt, l );
  }
  
  START_LOCKED_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // wait for +mu direction
      ghost_wait_PRECISION( out, mu, +1, &(op->c), plus_dir_param, l );
    }
  }
  END_LOCKED_MASTER(threading)

  END_NO_HYPERTHREADS(threading)
#endif
}


void coarse_hopping_term_PRECISION_vectorized( vector_PRECISION out, vector_PRECISION in, operator_PRECISION_struct *op,
                                    const int amount, level_struct *l, struct Thread *threading ) {

#ifdef VECTORIZE_COARSE_OPERATOR_PRECISION
  START_NO_HYPERTHREADS(threading)

  int mu, i, index, num_site_var=l->num_lattice_site_var,
      start=0, num_lattice_sites=l->num_inner_lattice_sites,
      plus_dir_param=_FULL_SYSTEM, minus_dir_param=_FULL_SYSTEM;
  vector_PRECISION in_pt, out_pt;

  OPERATOR_TYPE_PRECISION *D_vectorized;
  int column_offset = SIMD_LENGTH_PRECISION*((l->num_lattice_site_var+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  int vectorized_link_offset = 2*l->num_lattice_site_var*column_offset;

  int core_start;
  int core_end;

  // assumptions (1) self coupling has already been performed
  //          OR (2) "out" is initialized with zeros
  set_boundary_PRECISION( out, 0, l, threading );

  if ( amount == _EVEN_SITES ) {
    minus_dir_param = _ODD_SITES;
    plus_dir_param = _EVEN_SITES;
  } else if ( amount == _ODD_SITES ) {
    minus_dir_param = _EVEN_SITES;
    plus_dir_param = _ODD_SITES;
  }

  START_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // communicate in -mu direction
      ghost_sendrecv_PRECISION( in, mu, -1, &(op->c), minus_dir_param, l );
    }
  }
  END_MASTER(threading)
  SYNC_CORES(threading)

  if ( amount == _EVEN_SITES ) {
    start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
  } else if ( amount == _ODD_SITES ) {
    start = 0; num_lattice_sites = op->num_even_sites;
  }
  compute_core_start_end_custom(start, start+num_lattice_sites, &core_start, &core_end, l, threading, 1);

  // compute U_mu^dagger coupling
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_vectorized = op->D_transformed_vectorized + 4*vectorized_link_offset*op->neighbor_table[index] + 0*vectorized_link_offset;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+T];
    coarse_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_vectorized = op->D_transformed_vectorized + 4*vectorized_link_offset*op->neighbor_table[index] + 1*vectorized_link_offset;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+Z];
    coarse_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_vectorized = op->D_transformed_vectorized + 4*vectorized_link_offset*op->neighbor_table[index] + 2*vectorized_link_offset;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+Y];
    coarse_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_vectorized = op->D_transformed_vectorized + 4*vectorized_link_offset*op->neighbor_table[index] + 3*vectorized_link_offset;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+X];
    coarse_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );
  }

  START_LOCKED_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // communicate in +mu direction
      ghost_sendrecv_PRECISION( out, mu, +1, &(op->c), plus_dir_param, l );
    }
    for ( mu=0; mu<4; mu++ ) {
      // wait for -mu direction
      ghost_wait_PRECISION( in, mu, -1, &(op->c), minus_dir_param, l );
    }
  }
  END_LOCKED_MASTER(threading)

  if ( amount == _EVEN_SITES ) {
    start = 0; num_lattice_sites = op->num_even_sites;
  } else if ( amount == _ODD_SITES ) {
    start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
  }
  compute_core_start_end_custom(start, start+num_lattice_sites, &core_start, &core_end, l, threading, 1);

  // compute U_mu couplings
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    out_pt = out + num_site_var*op->neighbor_table[index];
    D_vectorized = op->D_vectorized + 4*vectorized_link_offset*op->neighbor_table[index];
    index++;
    in_pt = in + num_site_var*op->neighbor_table[index+T];
    coarse_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );

    D_vectorized += vectorized_link_offset;
    in_pt = in + num_site_var*op->neighbor_table[index+Z];
    coarse_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );

    in_pt = in + num_site_var*op->neighbor_table[index+Y];
    D_vectorized += vectorized_link_offset;
    coarse_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );

    in_pt = in + num_site_var*op->neighbor_table[index+X];
    D_vectorized += vectorized_link_offset;
    coarse_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );
  }

  START_LOCKED_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // wait for +mu direction
      ghost_wait_PRECISION( out, mu, +1, &(op->c), plus_dir_param, l );
    }
  }
  END_LOCKED_MASTER(threading)

  END_NO_HYPERTHREADS(threading)
#endif
}


void coarse_pn_hopping_term_PRECISION_vectorized( vector_PRECISION out, vector_PRECISION in, operator_PRECISION_struct *op,
                                    const int amount, level_struct *l, int sign, struct Thread *threading ) {

#ifdef VECTORIZE_COARSE_OPERATOR_PRECISION
  START_NO_HYPERTHREADS(threading)

  int mu, i, num_site_var=l->num_lattice_site_var,
      start=0, num_lattice_sites=l->num_inner_lattice_sites,
      plus_dir_param=_FULL_SYSTEM, minus_dir_param=_FULL_SYSTEM;
  vector_PRECISION in_pt, out_pt;

  OPERATOR_TYPE_PRECISION *D_vectorized;
  int column_offset = SIMD_LENGTH_PRECISION*((l->num_lattice_site_var+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  int link_offset = 2*l->num_lattice_site_var*column_offset;
  int *neighbor_fw = op->neighbor_table;
  int *neighbor_bw = op->backward_neighbor_table;

  int core_start;
  int core_end;

  void (*coarse_hopp)(vector_PRECISION eta, vector_PRECISION phi, OPERATOR_TYPE_PRECISION *D, level_struct *l);
  if(sign == +1)
    coarse_hopp = coarse_hopp_PRECISION_vectorized;
  else
    coarse_hopp = coarse_n_hopp_PRECISION_vectorized;


  if ( l->num_processes > 1 && op->c.comm ) {
    set_boundary_PRECISION( out, 0, l, threading );

    if ( amount == _EVEN_SITES ) {
      minus_dir_param = _ODD_SITES;
      plus_dir_param = _EVEN_SITES;
    } else if ( amount == _ODD_SITES ) {
      minus_dir_param = _EVEN_SITES;
      plus_dir_param = _ODD_SITES;
    }

    START_MASTER(threading)
    for ( mu=0; mu<4; mu++ ) {
      // send in -mu direction
      ghost_sendrecv_PRECISION( in, mu, -1, &(op->c), minus_dir_param, l );
    }
    END_MASTER(threading)

    if ( amount == _EVEN_SITES ) {
      start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
    } else if ( amount == _ODD_SITES ) {
      start = 0; num_lattice_sites = op->num_even_sites;
    }
    compute_core_start_end_custom(start, start+num_lattice_sites, &core_start, &core_end, l, threading, 1);

    // prepare for sending to fw: compute hopping terms into forward boundary buffer
    for ( i=core_start; i<core_end; i++ ) {
      for(int mu=0; mu<4; mu++) {
        if(neighbor_fw[5*i+1+mu] < l->num_inner_lattice_sites)
          continue;
        out_pt = out + num_site_var*neighbor_fw[5*i+1+mu];
        in_pt = in + num_site_var*neighbor_fw[5*i];
        D_vectorized = op->D_transformed_vectorized + 4*link_offset*neighbor_fw[5*i] + mu*link_offset;
        coarse_hopp( out_pt, in_pt, D_vectorized, l );
      }
    }
    START_LOCKED_MASTER(threading)
    for ( mu=0; mu<4; mu++ ) {
      // send in +mu direction
      ghost_sendrecv_PRECISION( out, mu, +1, &(op->c), plus_dir_param, l );
    }
    for ( mu=0; mu<4; mu++ ) {
      // wait for -mu direction
      ghost_wait_PRECISION( in, mu, -1, &(op->c), minus_dir_param, l );
    }
    END_LOCKED_MASTER(threading)
  }
  else
    SYNC_CORES(threading)


  if ( amount == _EVEN_SITES ) {
    start = 0; num_lattice_sites = op->num_even_sites;
  } else if ( amount == _ODD_SITES ) {
    start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
  }
  compute_core_start_end_custom(start, start+num_lattice_sites, &core_start, &core_end, l, threading, 1);

  // assumptions (1) self coupling has already been performed
  //          OR (2) "out" is initialized with zeros
  for ( i=core_start; i<core_end; i++ ) {
    out_pt = out + num_site_var*neighbor_fw[5*i];

    // U_mu^dagger coupling
    for(int mu=0; mu<4; mu++) {
      // terms coming from backward boundary buffer are done by the ghost_wait_PRECISION call below
      if(neighbor_bw[5*i+1+mu] >= l->num_inner_lattice_sites)
        continue;
      D_vectorized = op->D_transformed_vectorized + 4*link_offset*neighbor_bw[5*i+1+mu] + mu*link_offset;
      in_pt = in + num_site_var*neighbor_bw[5*i+1+mu];
      coarse_hopp( out_pt, in_pt, D_vectorized, l );
    }

    // compute U_mu couplings
    for(int mu=0; mu<4; mu++) {
      D_vectorized = op->D_vectorized + 4*link_offset*neighbor_fw[5*i] + mu*link_offset;
      in_pt = in + num_site_var*neighbor_fw[5*i+1+mu];
      coarse_hopp( out_pt, in_pt, D_vectorized, l );
    }
  }


  // wait for terms from bw and add them
  if ( l->num_processes > 1 && op->c.comm ) {
    START_LOCKED_MASTER(threading)
    for ( mu=0; mu<4; mu++ ) {
      // wait for +mu direction
      ghost_wait_PRECISION( out, mu, +1, &(op->c), plus_dir_param, l );
    }
    END_LOCKED_MASTER(threading)
  }
  else
    SYNC_CORES(threading)

  END_NO_HYPERTHREADS(threading)
#endif
}


void coarse_n_hopping_term_PRECISION_vectorized( vector_PRECISION out, vector_PRECISION in, operator_PRECISION_struct *op,
                                      const int amount, level_struct *l, struct Thread *threading ) {

#ifdef VECTORIZE_COARSE_OPERATOR_PRECISION
  START_NO_HYPERTHREADS(threading)

  int mu, i, index, num_site_var=l->num_lattice_site_var,
      start=0, num_lattice_sites=l->num_inner_lattice_sites,
      plus_dir_param=_FULL_SYSTEM, minus_dir_param=_FULL_SYSTEM;
  vector_PRECISION in_pt, out_pt;

  OPERATOR_TYPE_PRECISION *D_vectorized;
  int column_offset = SIMD_LENGTH_PRECISION*((l->num_lattice_site_var+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  int vectorized_link_offset = 2*l->num_lattice_site_var*column_offset;

  int core_start;
  int core_end;

  // assumptions (1) self coupling has already been performed
  //          OR (2) "out" is initialized with zeros
  set_boundary_PRECISION( out, 0, l, threading );

  if ( amount == _EVEN_SITES ) {
    minus_dir_param = _ODD_SITES;
    plus_dir_param = _EVEN_SITES;
  } else if ( amount == _ODD_SITES ) {
    minus_dir_param = _EVEN_SITES;
    plus_dir_param = _ODD_SITES;
  }

  START_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // communicate in -mu direction
      ghost_sendrecv_PRECISION( in, mu, -1, &(op->c), minus_dir_param, l );
    }
  }
  END_MASTER(threading)
  SYNC_CORES(threading)

  if ( amount == _EVEN_SITES ) {
    start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
  } else if ( amount == _ODD_SITES ) {
    start = 0; num_lattice_sites = op->num_even_sites;
  }
  compute_core_start_end_custom(start, start+num_lattice_sites, &core_start, &core_end, l, threading, 1);

  // D is applied in an input-centric way
  // this makes threading a bit ugly, is there a better way?
  // compute U_mu^dagger coupling
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_vectorized = op->D_transformed_vectorized + 4*vectorized_link_offset*op->neighbor_table[index] + 0*vectorized_link_offset;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+T];
    coarse_n_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_vectorized = op->D_transformed_vectorized + 4*vectorized_link_offset*op->neighbor_table[index] + 1*vectorized_link_offset;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+Z];
    coarse_n_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_vectorized = op->D_transformed_vectorized + 4*vectorized_link_offset*op->neighbor_table[index] + 2*vectorized_link_offset;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+Y];
    coarse_n_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    D_vectorized = op->D_transformed_vectorized + 4*vectorized_link_offset*op->neighbor_table[index] + 3*vectorized_link_offset;
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+X];
    coarse_n_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );
  }

  START_LOCKED_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // communicate in +mu direction
      ghost_sendrecv_PRECISION( out, mu, +1, &(op->c), plus_dir_param, l );
    }
    for ( mu=0; mu<4; mu++ ) {
      // wait for -mu direction
      ghost_wait_PRECISION( in, mu, -1, &(op->c), minus_dir_param, l );
    }
  }
  END_LOCKED_MASTER(threading)

  if ( amount == _EVEN_SITES ) {
    start = 0; num_lattice_sites = op->num_even_sites;
  } else if ( amount == _ODD_SITES ) {
    start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
  }
  compute_core_start_end_custom(start, start+num_lattice_sites, &core_start, &core_end, l, threading, 1);

  // compute U_mu couplings
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    out_pt = out + num_site_var*op->neighbor_table[index];
    D_vectorized = op->D_vectorized + 4*vectorized_link_offset*op->neighbor_table[index];
    index++;
    in_pt = in + num_site_var*op->neighbor_table[index+T];
    coarse_n_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );

    D_vectorized += vectorized_link_offset;
    in_pt = in + num_site_var*op->neighbor_table[index+Z];
    coarse_n_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );

    D_vectorized += vectorized_link_offset;
    in_pt = in + num_site_var*op->neighbor_table[index+Y];
    coarse_n_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );

    D_vectorized += vectorized_link_offset;
    in_pt = in + num_site_var*op->neighbor_table[index+X];
    coarse_n_hopp_PRECISION_vectorized( out_pt, in_pt, D_vectorized, l );
  }

  START_LOCKED_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // wait for +mu direction
      ghost_wait_PRECISION( out, mu, +1, &(op->c), plus_dir_param, l );
    }
  }
  END_LOCKED_MASTER(threading)

  END_NO_HYPERTHREADS(threading)
#endif
}


void coarse_solve_odd_even_PRECISION( gmres_PRECISION_struct *p, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {
  
  SYNC_CORES(threading)
  PROF_PRECISION_START( _SC, threading );
  coarse_diag_oo_inv_PRECISION( p->x, p->b, op, l, threading );
  PROF_PRECISION_STOP( _SC, 0, threading );
  PROF_PRECISION_START( _NC, threading );
  coarse_n_hopping_term_PRECISION( p->b, p->x, op, _EVEN_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 0, threading );
  
  fgmres_PRECISION( p, l, threading );
  
  // even to odd
  PROF_PRECISION_START( _NC, threading );
  coarse_n_hopping_term_PRECISION( p->b, p->x, op, _ODD_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 1, threading );
  PROF_PRECISION_START( _SC, threading );
  coarse_diag_oo_inv_PRECISION( p->x, p->b, op, l, threading );
  PROF_PRECISION_STOP( _SC, 1, threading );
  SYNC_CORES(threading)
}


void coarse_apply_schur_complement_PRECISION( vector_PRECISION out, vector_PRECISION in, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {
    
  // start and end indices for vector functions depending on thread
  int start;
  int end;
  // compute start and end indices for core
  // this puts zero for all other hyperthreads, so we can call functions below with all hyperthreads
  compute_core_start_end(op->num_even_sites*l->num_lattice_site_var, l->inner_vector_size, &start, &end, l, threading);

  vector_PRECISION *tmp = op->buffer;
  
  SYNC_CORES(threading)
  PROF_PRECISION_START( _SC, threading );
  coarse_diag_ee_PRECISION( out, in, op, l, threading );
  PROF_PRECISION_STOP( _SC, 0, threading );
  SYNC_CORES(threading)
  vector_PRECISION_define( tmp[0], 0, start, end, l );
  SYNC_CORES(threading)
  PROF_PRECISION_START( _NC, threading );
  coarse_hopping_term_PRECISION( tmp[0], in, op, _ODD_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 0, threading );
  PROF_PRECISION_START( _SC, threading );
  coarse_diag_oo_inv_PRECISION( tmp[1], tmp[0], op, l, threading );
  PROF_PRECISION_STOP( _SC, 1, threading );
  PROF_PRECISION_START( _NC, threading );
  coarse_n_hopping_term_PRECISION( out, tmp[1], op, _EVEN_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 1, threading );
}


void g5D_coarse_solve_odd_even_PRECISION( gmres_PRECISION_struct *p, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {
  
  int start_even, end_even, start_odd, end_odd;
  compute_core_start_end_custom(0, op->num_even_sites*l->num_lattice_site_var, &start_even, &end_even, l, threading, l->num_lattice_site_var );
  compute_core_start_end_custom(op->num_even_sites*l->num_lattice_site_var, l->inner_vector_size, &start_odd, &end_odd, l, threading, l->num_lattice_site_var );
  
  vector_PRECISION tmp = op->buffer[0];
  
  SYNC_CORES(threading)
  vector_PRECISION_define( tmp, 0, start_even, end_even, l );
  
  SYNC_CORES(threading)
  PROF_PRECISION_START( _SC, threading );
  coarse_gamma5_PRECISION( p->b, p->b, start_odd, end_odd, l );
  SYNC_CORES(threading)
  coarse_diag_oo_inv_PRECISION( p->x, p->b, op, l, threading );
  SYNC_CORES(threading)
  coarse_gamma5_PRECISION( p->b, p->b, start_odd, end_odd, l );
  PROF_PRECISION_STOP( _SC, 0, threading );
  PROF_PRECISION_START( _NC, threading );
  coarse_n_hopping_term_PRECISION( tmp, p->x, op, _EVEN_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 0, threading );
  coarse_gamma5_PRECISION( tmp, tmp, start_even, end_even, l );
  SYNC_CORES(threading)
  vector_PRECISION_plus( p->b, p->b, tmp, start_even, end_even, l );
  
  fgmres_PRECISION( p, l, threading );
  SYNC_CORES(threading)
  coarse_gamma5_PRECISION( p->b, p->b, start_odd, end_odd, l );
  SYNC_CORES(threading)
  coarse_diag_oo_inv_PRECISION( p->x, p->b, op, l, threading );
  SYNC_CORES(threading)
  
  // even to odd
  PROF_PRECISION_START( _NC, threading );
  vector_PRECISION_define( tmp, 0, start_odd, end_odd, l );
  SYNC_CORES(threading)
  coarse_n_hopping_term_PRECISION( tmp, p->x, op, _ODD_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 1, threading );
  SYNC_CORES(threading)
  PROF_PRECISION_START( _SC, threading );
  coarse_diag_oo_inv_PRECISION( p->b, tmp, op, l, threading );
  vector_PRECISION_plus( p->x, p->x, p->b, start_odd, end_odd, l );
  
  PROF_PRECISION_STOP( _SC, 1, threading );
  SYNC_CORES(threading)
}


void g5D_coarse_apply_schur_complement_PRECISION( vector_PRECISION out, vector_PRECISION in, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {
    
  int start_even, end_even, start_odd, end_odd;
  compute_core_start_end_custom(0, op->num_even_sites*l->num_lattice_site_var, &start_even, &end_even, l, threading, l->num_lattice_site_var );
  compute_core_start_end_custom(op->num_even_sites*l->num_lattice_site_var, l->inner_vector_size, &start_odd, &end_odd, l, threading, l->num_lattice_site_var );

  vector_PRECISION *tmp = op->buffer;
  
  SYNC_CORES(threading)
  PROF_PRECISION_START( _SC, threading );
  coarse_diag_ee_PRECISION( out, in, op, l, threading );
  PROF_PRECISION_STOP( _SC, 0, threading );
  SYNC_CORES(threading)
  vector_PRECISION_define( tmp[0], 0, start_odd, end_odd, l );
  SYNC_CORES(threading)
  PROF_PRECISION_START( _NC, threading );
  coarse_hopping_term_PRECISION( tmp[0], in, op, _ODD_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 0, threading );
  PROF_PRECISION_START( _SC, threading );
  coarse_diag_oo_inv_PRECISION( tmp[1], tmp[0], op, l, threading );
  PROF_PRECISION_STOP( _SC, 1, threading );
  PROF_PRECISION_START( _NC, threading );
  coarse_n_hopping_term_PRECISION( out, tmp[1], op, _EVEN_SITES, l, threading );
  PROF_PRECISION_STOP( _NC, 1, threading );
  SYNC_CORES(threading)
  coarse_gamma5_PRECISION( out, out, start_even, end_even, l );
  SYNC_CORES(threading)
}


void coarse_odd_even_PRECISION_test( vector_PRECISION out, vector_PRECISION in, level_struct *l, struct Thread *threading ) {
  
  if ( g.odd_even ) {
    vector_PRECISION buf1 = NULL, buf2 = NULL;
    
    PUBLIC_MALLOC( buf1, complex_PRECISION, 2*l->vector_size );
    buf2 = buf1 + l->vector_size;

    START_LOCKED_MASTER(threading)
    // transformation part
    vector_PRECISION_copy( buf1, in, 0, l->inner_vector_size, l );
    // even to odd
    vector_PRECISION_define( out, 0, l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l->inner_vector_size, l );
    END_LOCKED_MASTER(threading)

    coarse_hopping_term_PRECISION( out, buf1, &(l->oe_op_PRECISION), _ODD_SITES, l, threading );
    coarse_diag_oo_inv_PRECISION( buf2, out, &(l->oe_op_PRECISION), l, threading );

    START_LOCKED_MASTER(threading)
    vector_PRECISION_plus( buf1, buf1, buf2, l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l->inner_vector_size, l );
    END_LOCKED_MASTER(threading)
    
    // block diagonal part
    if ( g.method == 6 ) {
      g5D_coarse_apply_schur_complement_PRECISION( out, buf1, &(l->oe_op_PRECISION), l, threading );
    } else {
      coarse_apply_schur_complement_PRECISION( out, buf1, &(l->oe_op_PRECISION), l, threading );
    }
    
    coarse_diag_oo_PRECISION( out, buf1, &(l->oe_op_PRECISION), l, threading );
    
    // back transformation part
    coarse_diag_oo_inv_PRECISION( buf2, out, &(l->oe_op_PRECISION), l, threading );
    
    if ( g.method == 6 ) {
      START_LOCKED_MASTER(threading)
      coarse_gamma5_PRECISION( out, out, l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l->inner_vector_size, l );
      vector_PRECISION_define( buf1, 0, 0, l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l );
      coarse_hopping_term_PRECISION( buf1, buf2, &(l->oe_op_PRECISION), _EVEN_SITES, l, no_threading );
      coarse_gamma5_PRECISION( buf1, buf1, 0, l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l );
      vector_PRECISION_plus( out, out, buf1, 0, l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l );
      END_LOCKED_MASTER(threading)
    } else {
      coarse_hopping_term_PRECISION( out, buf2, &(l->oe_op_PRECISION), _EVEN_SITES, l, threading );
    }

    PUBLIC_FREE( buf1, complex_PRECISION, 2*l->vector_size );
  }
}
