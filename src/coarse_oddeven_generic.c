/*
 * Copyright (C) 2016, Matthias Rottmann, Artur Strebel, Simon Heybrock, Simone Bacchio, Bjoern Leder.
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

void coarse_selfcoupling_LU_decomposition_PRECISION( config_PRECISION output, operator_PRECISION_struct *op, int index, level_struct *l ) {

  // clover = [ A B      , A=A*, D=D*, C = -B*
  //            C D ]
  //
  // order: upper triangle of A, upper triangle of D, B, each column major
  //
  // tm_term = [ E 0      , E=-E*, F=-F* diag. excluded
  //             0 F ]
  //
  // order: upper triangle of E, upper triangle of F
  //
  // output = [ A+E  B   
  //             C  D+F ] LU decomposed

  register int i, j, k, n = l->num_parent_eig_vect, n2 = 2*n;
  config_PRECISION clover = op->clover + n*(n2+1)*index;
  // A
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<j; i++ ) {
      output[n2*i+j] = *clover;
      output[i+n2*j] = conj_PRECISION(*clover);
      clover++;      
    }
    output[(n2+1)*j] = *clover;
    clover++; // diagonal entry
  }
  // D
  for ( j=n; j<n2; j++ ) {
    for ( i=n; i<j; i++ ) {
      output[n2*i+j] = *clover;
      output[i+n2*j] = conj_PRECISION(*clover);
      clover++;      
    }
    output[(n2+1)*j] = *clover;
    clover++; // diagonal entry
  }
  // B and C
  for ( j=n; j<n2; j++ ) {
    for ( i=0; i<n; i++ ) {
      output[n2*i+j] = *clover;
      output[i+n2*j] = -conj_PRECISION(*clover);
      clover++;      
    }
  }

#ifdef HAVE_TM
  config_PRECISION tm_term = op->tm_term + n*(n+1)*index;
  if (op->mu + op->mu_odd_shift != 0.0 || op->mu + op->mu_even_shift != 0.0 ) {
    // E
    for ( j=0; j<n; j++ ) {
      for ( i=0; i<j; i++ ) {
        output[n2*i+j] += *tm_term;
        output[i+n2*j] += -conj_PRECISION(*tm_term);
        tm_term++;      
      }
      output[(n2+1)*j] += *tm_term;
      tm_term++; // diagonal entry
    }
    // F
    for ( j=n; j<n2; j++ ) {
      for ( i=n; i<j; i++ ) {
        output[n2*i+j] += *tm_term;
        output[i+n2*j] += -conj_PRECISION(*tm_term);
        tm_term++;      
      }
      output[(n2+1)*j] += *tm_term;
      tm_term++; // diagonal entry
    }
  }
#endif
    
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

#ifdef HAVE_TM1p1
void coarse_selfcoupling_LU_doublet_decomposition_PRECISION( config_PRECISION output, operator_PRECISION_struct *op, int index, level_struct *l ) {

  // clover = [ A B      , A=A*, D=D*, C = -B*
  //            C D ]
  //
  // order: upper triangle of A, upper triangle of D, B, each column major
  //
  // tm_term = [ E 0      , E=-E*, F=-F* diag. excluded
  //             0 F ]
  //
  // order: upper triangle of E, upper triangle of F
  //
  // epsbar_term = [ G 0      , G=-G*, H=-H* diag. excluded
  //                 0 H ]
  //
  // order: upper triangle of G, upper triangle of H
  //
  // output = [ A+E  G   B   0
  //             G  A-E  0   B
  //             C   0  D+F  H
  //             0   C   H  D-F ]  LU decomposed

  register int i, j, k, n = l->num_parent_eig_vect, n2 = 2*n, n3 = 3*n, n4 = 4*n;
  // set the matrix up
  // 0
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<n; i++ ) {
      output[n4*(i+0 )+(j+n3)] = _COMPLEX_PRECISION_ZERO;
      output[n4*(i+n )+(j+n2)] = _COMPLEX_PRECISION_ZERO;
      output[n4*(i+n2)+(j+n )] = _COMPLEX_PRECISION_ZERO;
      output[n4*(i+n3)+(j+0 )] = _COMPLEX_PRECISION_ZERO;
    }
  }

  config_PRECISION clover = op->clover + n*(n2+1)*index;
  // A
  for ( j=0; j<n; j++ ) {
    for ( i=0; i<j; i++ ) {
      output[n4*i+j] = *clover;
      output[i+n4*j] = conj_PRECISION(*clover);
      output[n4*(i+n)+(j+n)] = *clover;
      output[(i+n)+n4*(j+n)] = conj_PRECISION(*clover);
      clover++;      
    }
    output[(n4+1)*j] = *clover;
    output[(n4+1)*(j+n)] = *clover;
    clover++; // diagonal entry
  }
  // D
  for ( j=n2; j<n3; j++ ) {
    for ( i=n2; i<j; i++ ) {
      output[n4*i+j] = *clover;
      output[i+n4*j] = conj_PRECISION(*clover);
      output[n4*(i+n)+(j+n)] = *clover;
      output[(i+n)+n4*(j+n)] = conj_PRECISION(*clover);
      clover++;      
    }
    output[(n4+1)*j] = *clover;
    output[(n4+1)*(j+n)] = *clover;
    clover++; // diagonal entry
  }
  // B and C
  for ( j=n2; j<n3; j++ ) {
    for ( i=0; i<n; i++ ) {
      output[n4*i+j] = *clover;
      output[i+n4*j] = -conj_PRECISION(*clover);
      output[n4*(i+n)+(j+n)] = *clover;
      output[(i+n)+n4*(j+n)] = -conj_PRECISION(*clover);
      clover++;      
    }
  }

#ifdef HAVE_TM
  config_PRECISION tm_term = op->tm_term + n*(n+1)*index;
  if (op->mu + op->mu_odd_shift != 0.0 || op->mu + op->mu_even_shift != 0.0 ) {
    // E
    for ( j=0; j<n; j++ ) {
      for ( i=0; i<j; i++ ) {
        output[n4*i+j] += *tm_term;
        output[i+n4*j] += -conj_PRECISION(*tm_term);
        output[n4*(i+n)+(j+n)] -= *tm_term;
        output[(i+n)+n4*(j+n)] -= -conj_PRECISION(*tm_term);
        tm_term++;      
      }
      output[(n4+1)*j] += *tm_term;
      output[(n4+1)*(j+n)] -= *tm_term;
      tm_term++; // diagonal entry
    }
    // F
    for ( j=n2; j<n3; j++ ) {
      for ( i=n2; i<j; i++ ) {
        output[n4*i+j] += *tm_term;
        output[i+n4*j] += -conj_PRECISION(*tm_term);
        output[n4*(i+n)+(j+n)] -= *tm_term;
        output[(i+n)+n4*(j+n)] -= -conj_PRECISION(*tm_term);
        tm_term++;      
      }
      output[(n4+1)*j] += *tm_term;
      output[(n4+1)*(j+n)] -= *tm_term;
      tm_term++; // diagonal entry
    }
  }
#endif

  config_PRECISION epsbar_term = op->epsbar_term + n*(n+1)*index;
  // G
  for ( j=n; j<n2; j++ ) {
    for ( i=0; i<(j-n); i++ ) {
      output[n4*i+j] = (*epsbar_term);
      output[(i+n)+n4*(j-n)] = -conj_PRECISION(*epsbar_term);
      output[n4*(i+n)+(j-n)] = (*epsbar_term);
      output[i+n4*j] = -conj_PRECISION(*epsbar_term);
      epsbar_term++;      
    }
    output[(n4+1)*(j-n)+n] = (*epsbar_term);
    output[(n4+1)*j-n] = (*epsbar_term);
    epsbar_term++; // diagonal entry
  }
  // H
  for ( j=n3; j<n4; j++ ) {
    for ( i=n2; i<(j-n); i++ ) {
      output[n4*i+j] = (*epsbar_term);
      output[(i+n)+n4*(j-n)] = -conj_PRECISION(*epsbar_term);
      output[n4*(i+n)+(j-n)] = (*epsbar_term);
      output[i+n4*j] = -conj_PRECISION(*epsbar_term);
      epsbar_term++;      
    }
    output[(n4+1)*(j-n)+n] = (*epsbar_term);
    output[(n4+1)*j-n] = (*epsbar_term);
    epsbar_term++; // diagonal entry
  }
    
  // compute LU decomposition
  // output = triu(L,1) + tril(U,0)
  // i.e., output contains L and U without the diagonal of L which is equal to 1
  // order: row major
  for ( k=0; k<n4; k++ ) {
    for ( i=k+1; i<n4; i++ ) {
      output[n4*i+k] = output[n4*i+k]/output[(n4+1)*k]; // output(i,k) = output(i,k)/output(k,k)
      for ( j=k+1; j<n4; j++ )
        output[n4*i+j] = output[n4*i+j]-output[n4*i+k]*output[n4*k+j]; // output(i,j) = output(i,j)-output(i,k)*output(k,j)
    }
  }
}
#endif


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


void coarse_diag_ee_PRECISION( vector_PRECISION y, vector_PRECISION x, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {
  
  int start, end;
  compute_core_start_end_custom( 0, op->num_even_sites, &start, &end, l, threading, 1 );
  // even sites
#ifndef OPTIMIZED_COARSE_SELF_COUPLING_PRECISION
  coarse_self_couplings_PRECISION( y, x, op, start, end, l );
#else
  coarse_self_couplings_PRECISION_vectorized( y, x, op, start, end, l );
#endif
}

void coarse_diag_oo_PRECISION( vector_PRECISION y, vector_PRECISION x, operator_PRECISION_struct *op, 
                               level_struct *l, struct Thread *threading ) {
  
  int start, end;
#ifndef OPTIMIZED_COARSE_SELF_COUPLING_PRECISION 
  int num_site_var=l->num_lattice_site_var,
    oo_inv_size = SQUARE(num_site_var);
#ifdef HAVE_TM1p1
  config_PRECISION sc = (g.n_flavours==2) ? op->clover_doublet_oo_inv:op->clover_oo_inv;
#else
  config_PRECISION sc = op->clover_oo_inv;
#endif

  compute_core_start_end_custom( 0, op->num_odd_sites, &start, &end, l, threading, 1 );

  x += num_site_var*(op->num_even_sites+start);
  y += num_site_var*(op->num_even_sites+start);  
  sc += oo_inv_size*start;

  for ( int i=start; i<end; i++ ) {
    coarse_LU_multiply_PRECISION( y, x, sc, l );
    x += num_site_var;
    y += num_site_var;
    sc += oo_inv_size;
  }
  
#else
  compute_core_start_end_custom( op->num_even_sites, l->num_inner_lattice_sites, &start, &end, l, threading, 1 );
  coarse_self_couplings_PRECISION_vectorized( y, x, op, start, end, l );
#endif
}

void coarse_diag_PRECISION( vector_PRECISION y, vector_PRECISION x, operator_PRECISION_struct *op, level_struct *l ) {
  
  coarse_diag_ee_PRECISION( y, x, op, l, no_threading );
  coarse_diag_oo_PRECISION( y, x, op, l, no_threading );
}

void coarse_diag_oo_inv_PRECISION( vector_PRECISION y, vector_PRECISION x, operator_PRECISION_struct *op, 
                               level_struct *l, struct Thread *threading ) {
  
  int start, end;
  compute_core_start_end_custom( 0, op->num_odd_sites, &start, &end, l, threading, 1 );
  
  // odd sites
  int num_site_var = l->num_lattice_site_var,
    oo_inv_size = SQUARE(num_site_var);

#ifndef OPTIMIZED_COARSE_SELF_COUPLING_PRECISION
#ifdef HAVE_TM1p1
  config_PRECISION sc = (g.n_flavours==2) ? op->clover_doublet_oo_inv:op->clover_oo_inv;
#else
  config_PRECISION sc = op->clover_oo_inv;
#endif
#else
  int lda = SIMD_LENGTH_PRECISION*((num_site_var+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  oo_inv_size = 2*num_site_var*lda;
#ifdef HAVE_TM1p1
  OPERATOR_TYPE_PRECISION *sc = (g.n_flavours==2) ? op->clover_doublet_oo_inv_vectorized:op->clover_oo_inv_vectorized;
#else
  OPERATOR_TYPE_PRECISION *sc = op->clover_oo_inv_vectorized;
#endif
#endif

  x += num_site_var*(op->num_even_sites+start);
  y += num_site_var*(op->num_even_sites+start);  
  sc += oo_inv_size*start;

  for ( int i=start; i<end; i++ ) {
#ifndef OPTIMIZED_COARSE_SELF_COUPLING_PRECISION
    coarse_perform_fwd_bwd_subs_PRECISION( y, x, sc, l );
#else
    for(int j=0; j<num_site_var; j++)
      y[j] = _COMPLEX_PRECISION_ZERO;
    cgemv_PRECISION( num_site_var, sc, lda, (PRECISION *)x, (PRECISION *)y);
#endif
    x += num_site_var;
    y += num_site_var;
    sc += oo_inv_size;
  }
}


void coarse_oddeven_PRECISION_set_self_couplings( level_struct *l, struct Thread *threading ) {

  operator_PRECISION_struct *op = &(l->oe_op_PRECISION);
  int nv = l->num_parent_eig_vect, start, end;

  coarse_operator_PRECISION_set_self_couplings( op, l, threading );
  compute_core_start_end_custom( 0, op->num_odd_sites, &start, &end, l, threading, 1);

#ifndef OPTIMIZED_COARSE_SELF_COUPLING_PRECISION

  int size = SQUARE(2*nv);
  for( int i=start; i<end; i++ )
    coarse_selfcoupling_LU_decomposition_PRECISION( op->clover_oo_inv+i*size, op, op->num_even_sites+i, l );

#ifdef HAVE_TM1p1
  int size_doublet = SQUARE(4*nv);
  for( int i=start; i<end; i++ )
    coarse_selfcoupling_LU_doublet_decomposition_PRECISION( op->clover_doublet_oo_inv+i*size_doublet, op, 
                                                            op->num_even_sites+i, l );
#endif

#else

  int column_offset = SIMD_LENGTH_PRECISION*((2*nv+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  int size_v = 2*2*nv*column_offset;
  for( int i=start; i<end; i++ )
    cgem_inverse_PRECISION( 2*nv, op->clover_oo_inv_vectorized + i*size_v, 
                  op->clover_vectorized + (op->num_even_sites+i)*size_v, column_offset );

#ifdef HAVE_TM1p1
  int column_doublet_offset = SIMD_LENGTH_PRECISION*((4*nv+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  int size_doublet_v = 2*4*nv*column_doublet_offset;
  for( int i=start; i<end; i++ )
    cgem_inverse_PRECISION( 4*nv, op->clover_doublet_oo_inv_vectorized + i*size_doublet_v, 
                  op->clover_doublet_vectorized + (op->num_even_sites+i)*size_doublet_v, column_doublet_offset );
#endif

#endif
}

void coarse_oddeven_PRECISION_set_couplings( level_struct *l, struct Thread *threading ) {

  coarse_oddeven_PRECISION_set_self_couplings( l, threading );
  coarse_operator_PRECISION_set_neighbor_couplings( &(l->oe_op_PRECISION), l, threading );

}

void coarse_oddeven_alloc_PRECISION( level_struct *l ) {

  int nv = l->num_parent_eig_vect,
    oe_offset=0, mu, **bt = NULL,
    *eot = NULL, *nt = NULL, *tt = NULL, t, z, y, x, le[4], N[4];
  operator_PRECISION_struct *op = &(l->oe_op_PRECISION);

  operator_PRECISION_alloc( op, _ODDEVEN, l );

  // buffers
  MALLOC( op->buffer, complex_PRECISION*, 2 );
  op->buffer[0] = NULL;
#ifdef HAVE_TM1p1
  MALLOC( op->buffer[0], complex_PRECISION, 4*l->vector_size );
  op->buffer[1] = op->buffer[0] + 2*l->vector_size;  
#else
  MALLOC( op->buffer[0], complex_PRECISION, 2*l->vector_size );
  op->buffer[1] = op->buffer[0] + l->vector_size;  
#endif

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
  
#ifndef OPTIMIZED_COARSE_SELF_COUPLING_PRECISION

  MALLOC( op->clover_oo_inv, complex_PRECISION, SQUARE(2*nv)*op->num_odd_sites );
#ifdef HAVE_TM1p1
  MALLOC( op->clover_doublet_oo_inv, complex_PRECISION, SQUARE(4*nv)*op->num_odd_sites );
#endif

#else
  int column_offset = SIMD_LENGTH_PRECISION*((2*nv+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  MALLOC_HUGEPAGES( op->clover_oo_inv_vectorized, PRECISION, 2*2*nv*column_offset*op->num_odd_sites, 4*SIMD_LENGTH_PRECISION );
#ifdef HAVE_TM1p1
  int column_doublet_offset = SIMD_LENGTH_PRECISION*((4*nv+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  MALLOC_HUGEPAGES( op->clover_doublet_oo_inv_vectorized, PRECISION, 2*4*nv*column_doublet_offset*op->num_odd_sites, 4*SIMD_LENGTH_PRECISION );
#endif

#endif

  // define data layout
  eot = op->index_table;
  define_eot( eot, N, l );

  // neighbor table, translation table
  nt = op->neighbor_table;
  tt = op->translation_table;
  define_nt_bt_tt( nt, op->backward_neighbor_table, NULL, tt, eot, N, l );

  // boundary table
  bt = op->c.boundary_table;
  define_eo_bt( bt, eot, op->c.num_even_boundary_sites, op->c.num_odd_boundary_sites, op->c.num_boundary_sites, N, l );

  // ghost
  ghost_sendrecv_init_PRECISION( _COARSE_GLOBAL, &(op->c), l ) ;

  // solver
  if ( l->level == 0 )
    l->p_PRECISION.v_end = op->num_even_sites*l->num_lattice_site_var;
  else
    l->sp_PRECISION.v_end = op->num_even_sites*l->num_lattice_site_var;

}

void coarse_oddeven_setup_PRECISION( operator_PRECISION_struct *in, int reorder, level_struct *l, 
                                     struct Thread *threading ) {

  operator_PRECISION_struct *op = &(l->oe_op_PRECISION);

  START_LOCKED_MASTER(threading)
    int ns=l->num_inner_lattice_sites, nv = l->num_parent_eig_vect, i,
    D_size = 4*SQUARE(2*nv),
    clover_size = (nv)*(nv*2+1),
    block_size = (nv)*(nv+1);
  config_PRECISION D_in = in->D,
    clover_in = in->clover,
    odd_proj_in = in->odd_proj;

  // neighbor couplings
  if ( reorder ) {
    int t, z, y, x, index, *le = l->local_lattice, oe_offset = op->oe_offset,
      *it = in->index_table, *dt = in->table_dim;
    config_PRECISION D_oe = op->D, 
      D_eo = (op->D)+D_size*op->num_even_sites,
      clover_ee = op->clover,
      clover_oo = (op->clover)+clover_size*op->num_even_sites,
      odd_proj_ee = op->odd_proj,
      odd_proj_oo = op->odd_proj+block_size*op->num_even_sites;

    for ( t=0; t<le[T]; t++ )
      for ( z=0; z<le[Z]; z++ )
        for ( y=0; y<le[Y]; y++ )
          for ( x=0; x<le[X]; x++ ) {
            index = site_index( t, z, y, x, dt, it );
            if ( (t+z+y+x+oe_offset)%2 == 1 ) {
              for ( i=0; i<D_size; i++ ) 
                D_eo[i] = D_in[ index*D_size+i ];
              for ( i=0; i<clover_size; i++ )
                clover_oo[i] = clover_in[ index*clover_size+i ];
              for ( i=0; i<block_size; i++ )
                odd_proj_oo[i] = odd_proj_in[ index*block_size+i ];
              D_eo += D_size;
              clover_oo += clover_size;
              odd_proj_oo += block_size;
            } else {
              for ( i=0; i<D_size; i++ )
                D_oe[i] = D_in[ index*D_size+i ];
              for ( i=0; i<clover_size; i++ )
                clover_ee[i] = clover_in[ index*clover_size+i ];
              for ( i=0; i<block_size; i++ )
                odd_proj_ee[i] = odd_proj_in[ index*block_size+i ];
              D_oe += D_size;
              clover_ee += clover_size;
              odd_proj_ee += block_size;
            }
          }
    
  } else {
    for ( i=0; i<D_size*ns; i++ )
      op->D[i] = D_in[i];
    for ( i=0; i<clover_size*ns; i++ )
      op->clover[i] = clover_in[i];
    for ( i=0; i<block_size*ns; i++ ) {
      op->odd_proj[i] = odd_proj_in[i];
    }
    
  }
  END_LOCKED_MASTER(threading)
  
  op->m0 = in->m0;

#ifdef HAVE_TM
  tm_term_PRECISION_setup( in->mu, in->mu_even_shift, in->mu_odd_shift, op, l, threading );
#endif  
#ifdef HAVE_TM1p1
  epsbar_term_PRECISION_setup( in->epsbar, in->epsbar_ig5_even_shift, in->epsbar_ig5_odd_shift, op, l, threading );
#endif
  
  coarse_oddeven_PRECISION_set_couplings( l, threading );
  
}


void coarse_oddeven_free_PRECISION( level_struct *l ) {
  
  int nv = l->num_parent_eig_vect, vs = l->vector_size;
  operator_PRECISION_struct *op = &(l->oe_op_PRECISION);

  operator_PRECISION_free( op, _ODDEVEN, l );
  coarse_operator_PRECISION_free_vectorized( op, l );

#ifndef OPTIMIZED_COARSE_SELF_COUPLING_PRECISION

  FREE( op->clover_oo_inv, complex_PRECISION, SQUARE(2*nv)*op->num_odd_sites );
#ifdef HAVE_TM1p1
  FREE( op->clover_doublet_oo_inv, complex_PRECISION, SQUARE(4*nv)*op->num_odd_sites );
#endif

#else
  int column_offset = SIMD_LENGTH_PRECISION*((2*nv+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  FREE_HUGEPAGES( op->clover_oo_inv_vectorized, PRECISION, 2*2*nv*column_offset*op->num_odd_sites );
#ifdef HAVE_TM1p1
  int column_doublet_offset = SIMD_LENGTH_PRECISION*((4*nv+SIMD_LENGTH_PRECISION-1)/SIMD_LENGTH_PRECISION);
  FREE_HUGEPAGES( op->clover_doublet_oo_inv_vectorized, PRECISION, 2*4*nv*column_doublet_offset*op->num_odd_sites );
#endif

#endif
  
#ifdef HAVE_TM1p1
  FREE( op->buffer[0], complex_PRECISION, 4*vs );
#else
  FREE( op->buffer[0], complex_PRECISION, 2*vs );
#endif
  FREE( op->buffer, complex_PRECISION*, 2 );
}

void coarse_solve_odd_even_PRECISION( gmres_PRECISION_struct *p, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {
  
  SYNC_CORES(threading)
  PROF_PRECISION_START( _SC, threading );
  coarse_diag_oo_inv_PRECISION( p->x, p->b, op, l, threading );
  PROF_PRECISION_STOP( _SC, 0, threading );
  PROF_PRECISION_START( _NC, threading );
  coarse_pn_hopping_term_PRECISION( p->b, p->x, op, _EVEN_SITES, -1, l, threading );
  PROF_PRECISION_STOP( _NC, 0, threading );
  
  fgmres_PRECISION( p, l, threading );
  
  // even to odd
  PROF_PRECISION_START( _NC, threading );
  coarse_pn_hopping_term_PRECISION( p->b, p->x, op, _ODD_SITES, -1, l, threading );
  PROF_PRECISION_STOP( _NC, 1, threading );
  PROF_PRECISION_START( _SC, threading );
  coarse_diag_oo_inv_PRECISION( p->x, p->b, op, l, threading );
  PROF_PRECISION_STOP( _SC, 1, threading );
  SYNC_CORES(threading)
}

void coarse_apply_schur_complement_PRECISION( vector_PRECISION out, vector_PRECISION in, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {
    
  int start = op->num_even_sites*l->num_lattice_site_var;
  int end = l->inner_vector_size;
  vector_PRECISION *tmp = op->buffer;
  
  SYNC_CORES(threading)
  PROF_PRECISION_START( _SC, threading );
  coarse_diag_ee_PRECISION( out, in, op, l, threading );
  PROF_PRECISION_STOP( _SC, 0, threading );
  SYNC_CORES(threading);
  vector_PRECISION_define_zero( tmp[0], start, end, l, threading );
  SYNC_CORES(threading)
  PROF_PRECISION_START( _NC, threading );
  coarse_pn_hopping_term_PRECISION( tmp[0], in, op, _ODD_SITES, +1, l, threading );
  PROF_PRECISION_STOP( _NC, 0, threading );
  PROF_PRECISION_START( _SC, threading );
  coarse_diag_oo_inv_PRECISION( tmp[1], tmp[0], op, l, threading );
  PROF_PRECISION_STOP( _SC, 1, threading );
  PROF_PRECISION_START( _NC, threading );
  coarse_pn_hopping_term_PRECISION( out, tmp[1], op, _EVEN_SITES, -1, l, threading );
  PROF_PRECISION_STOP( _NC, 1, threading );
}


void g5D_coarse_solve_odd_even_PRECISION( gmres_PRECISION_struct *p, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {
  
  int start_even = 0, end_even = op->num_even_sites*l->num_lattice_site_var, 
    start_odd = end_even, end_odd = l->inner_vector_size;
  int thread_start_even, thread_end_even, thread_start_odd, thread_end_odd;
  compute_core_start_end_custom( start_even, end_even, &thread_start_even, &thread_end_even, l, threading, l->num_lattice_site_var );
  compute_core_start_end_custom( start_odd, end_odd, &thread_start_odd, &thread_end_odd, l, threading, l->num_lattice_site_var );
  
  vector_PRECISION tmp = op->buffer[0];
  
  SYNC_CORES(threading);
  vector_PRECISION_define_zero( tmp, start_even, end_even, l, threading );
  
  SYNC_CORES(threading)
  PROF_PRECISION_START( _SC, threading );
  coarse_gamma5_PRECISION( p->b, p->b, thread_start_odd, thread_end_odd, l );
  SYNC_CORES(threading)
  coarse_diag_oo_inv_PRECISION( p->x, p->b, op, l, threading );
  SYNC_CORES(threading)
  coarse_gamma5_PRECISION( p->b, p->b, thread_start_odd, thread_end_odd, l );
  PROF_PRECISION_STOP( _SC, 0, threading );
  PROF_PRECISION_START( _NC, threading );
  coarse_pn_hopping_term_PRECISION( tmp, p->x, op, _EVEN_SITES, -1, l, threading );
  PROF_PRECISION_STOP( _NC, 0, threading );
  coarse_gamma5_PRECISION( tmp, tmp, thread_start_even, thread_end_even, l );
  SYNC_CORES(threading)
  vector_PRECISION_plus( p->b, p->b, tmp, thread_start_even, thread_end_even, l );
  
  fgmres_PRECISION( p, l, threading );
  SYNC_CORES(threading)
  coarse_gamma5_PRECISION( p->b, p->b, thread_start_odd, thread_end_odd, l );
  SYNC_CORES(threading)
  coarse_diag_oo_inv_PRECISION( p->x, p->b, op, l, threading );
  SYNC_CORES(threading)
  
  // even to odd
  PROF_PRECISION_START( _NC, threading );
  vector_PRECISION_define_zero( tmp, start_odd, end_odd, l, threading );
  SYNC_CORES(threading);
  coarse_pn_hopping_term_PRECISION( tmp, p->x, op, _ODD_SITES, -1, l, threading );
  PROF_PRECISION_STOP( _NC, 1, threading );
  SYNC_CORES(threading)
  PROF_PRECISION_START( _SC, threading );
  coarse_diag_oo_inv_PRECISION( p->b, tmp, op, l, threading );
  vector_PRECISION_plus( p->x, p->x, p->b, thread_start_odd, thread_end_odd, l );
  
  PROF_PRECISION_STOP( _SC, 1, threading );
  SYNC_CORES(threading)
}


void g5D_coarse_apply_schur_complement_PRECISION( vector_PRECISION out, vector_PRECISION in, operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {
    
  int start_even = 0, end_even = op->num_even_sites*l->num_lattice_site_var, 
    start_odd = end_even, end_odd = l->inner_vector_size;
  int thread_start_even, thread_end_even;
  compute_core_start_end_custom( start_even, end_even, &thread_start_even, &thread_end_even, l, threading, l->num_lattice_site_var );

  vector_PRECISION *tmp = op->buffer;
  
  SYNC_CORES(threading)
  PROF_PRECISION_START( _SC, threading );
  coarse_diag_ee_PRECISION( out, in, op, l, threading );
  PROF_PRECISION_STOP( _SC, 0, threading );
  SYNC_CORES(threading);
  vector_PRECISION_define_zero( tmp[0], start_odd, end_odd, l, threading );
  SYNC_CORES(threading)
  PROF_PRECISION_START( _NC, threading );
  coarse_pn_hopping_term_PRECISION( tmp[0], in, op, _ODD_SITES, +1, l, threading );
  PROF_PRECISION_STOP( _NC, 0, threading );
  PROF_PRECISION_START( _SC, threading );
  coarse_diag_oo_inv_PRECISION( tmp[1], tmp[0], op, l, threading );
  PROF_PRECISION_STOP( _SC, 1, threading );
  PROF_PRECISION_START( _NC, threading );
  coarse_pn_hopping_term_PRECISION( out, tmp[1], op, _EVEN_SITES, -1, l, threading );
  PROF_PRECISION_STOP( _NC, 1, threading );
  SYNC_CORES(threading)
  coarse_gamma5_PRECISION( out, out, thread_start_even, thread_end_even, l );
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
    vector_PRECISION_define_zero( out, l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l->inner_vector_size, l, no_threading );
    END_LOCKED_MASTER(threading);

    coarse_pn_hopping_term_PRECISION( out, buf1, &(l->oe_op_PRECISION), _ODD_SITES, +1, l, threading );
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
      vector_PRECISION_define_zero( buf1, 0, l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l, no_threading );
      coarse_pn_hopping_term_PRECISION( buf1, buf2, &(l->oe_op_PRECISION), _EVEN_SITES, +1, l, no_threading );
      coarse_gamma5_PRECISION( buf1, buf1, 0, l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l );
      vector_PRECISION_plus( out, out, buf1, 0, l->oe_op_PRECISION.num_even_sites*l->num_lattice_site_var, l );
      END_LOCKED_MASTER(threading)
    } else {
      coarse_pn_hopping_term_PRECISION( out, buf2, &(l->oe_op_PRECISION), _EVEN_SITES, +1, l, threading );
    }

    PUBLIC_FREE( buf1, complex_PRECISION, 2*l->vector_size );
  }
}
