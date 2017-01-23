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

void operator_PRECISION_init( operator_PRECISION_struct *op ) {
  
  op->prnT = NULL;
  op->index_table = NULL;
  op->neighbor_table = NULL;
  op->backward_neighbor_table = NULL;
  op->translation_table = NULL;
  op->D = NULL;
  op->clover = NULL;
  op->oe_clover = NULL;
  op->oe_clover_vectorized = NULL;
  
  for ( int mu=0; mu<4; mu++ )
    op->config_boundary_table[mu] = NULL;
  
  for ( int i=0; i<8; i++ ) {
    op->c.boundary_table[i] = NULL;
    op->c.buffer[i] = NULL;
    op->c.in_use[i] = 0;
  }
  op->c.comm = 1;
  op->buffer = NULL;
#ifdef VECTORIZE_COARSE_OPERATOR_PRECISION
  op->D_vectorized = NULL;
  op->D_transformed_vectorized = NULL;
  op->clover_vectorized = NULL;
#endif
}


void operator_PRECISION_alloc_projection_buffers( operator_PRECISION_struct *op, level_struct *l ) {

  // when used as preconditioner we usually do not need the projection buffers, unless
  // g.method >= 4: then oddeven_setup_float() is called in init.c, method_setup().
  if ( l->depth == 0 ) {
    int its = (l->num_lattice_site_var/2)*l->num_lattice_sites;
    MALLOC( op->prnT, complex_PRECISION, its*8 );
    op->prnZ = op->prnT + its; op->prnY = op->prnZ + its; op->prnX = op->prnY + its;
    op->prpT = op->prnX + its; op->prpZ = op->prpT + its; op->prpY = op->prpZ + its; op->prpX = op->prpY + its;
  }
}
void operator_PRECISION_free_projection_buffers( operator_PRECISION_struct *op, level_struct *l ) {

  if ( l->depth == 0 ) {
    int its = (l->num_lattice_site_var/2)*l->num_lattice_sites;
    FREE( op->prnT, complex_PRECISION, its*8 );
  }
}

void operator_PRECISION_alloc( operator_PRECISION_struct *op, const int type, level_struct *l ) {
  
/*********************************************************************************
* Allocates space for setting up an operator.
* - operator_PRECISION_struct *op: operator struct for which space is allocated.
* - const int type: Defines the data layout type of the operator.
* Possible values are: { _ORDINARY, _SCHWARZ }
*********************************************************************************/

  int mu, nu, its = 1, its_boundary, nls, clover_site_size, coupling_site_size;
  
  if ( l->depth == 0 ) {
    clover_site_size = 42;
    coupling_site_size = 4*9;
  } else {
    clover_site_size = (l->num_lattice_site_var*(l->num_lattice_site_var+1))/2;
    coupling_site_size = 4*l->num_lattice_site_var*l->num_lattice_site_var;
  }
  
  if ( type ==_SCHWARZ ) {
    its_boundary = 2;
  } else {
    its_boundary = 1;
  }
  for ( mu=0; mu<4; mu++ ) {
    its *= (l->local_lattice[mu]+its_boundary);
  }
  
  nls = (type==_ORDINARY)?l->num_inner_lattice_sites:2*l->num_lattice_sites-l->num_inner_lattice_sites;
  MALLOC( op->D, complex_PRECISION, coupling_site_size*nls );
  MALLOC( op->clover, complex_PRECISION, clover_site_size*l->num_inner_lattice_sites );
  if ( type == _SCHWARZ && l->depth == 0 && g.odd_even )
    MALLOC( op->oe_clover, complex_PRECISION, clover_site_size*l->num_inner_lattice_sites );
  MALLOC( op->index_table, int, its );
  MALLOC( op->neighbor_table, int, (l->depth==0?4:5)*l->num_inner_lattice_sites );
  MALLOC( op->backward_neighbor_table, int, (l->depth==0?4:5)*l->num_inner_lattice_sites );
  MALLOC( op->translation_table, int, l->num_inner_lattice_sites );
#ifdef SSE
  if ( l->depth == 0 ) {
    MALLOC( op->oe_clover_vectorized, PRECISION, 144*l->num_inner_lattice_sites );    
  }
#endif
  
  operator_PRECISION_alloc_projection_buffers( op, l );
  
  ghost_alloc_PRECISION( 0, &(op->c), l );
  
  for ( mu=0; mu<4; mu++ ) {
    its = 1;
    for ( nu=0; nu<4; nu++ ) {
      if ( mu != nu ) {
        its *= l->local_lattice[nu];
      }
    }
    op->c.num_boundary_sites[2*mu] = its;
    op->c.num_boundary_sites[2*mu+1] = its;
    MALLOC( op->c.boundary_table[2*mu], int, its );
    if ( type == _SCHWARZ ) {
      MALLOC( op->c.boundary_table[2*mu+1], int, its );
      MALLOC( op->config_boundary_table[mu], int, its );
    } else {
      op->c.boundary_table[2*mu+1] = op->c.boundary_table[2*mu];
    }
  }
}


void operator_PRECISION_free( operator_PRECISION_struct *op, const int type, level_struct *l ) {
  
  int mu, nu, its = 1, clover_site_size, coupling_site_size;
  
  if ( l->depth == 0 ) {
    clover_site_size = 42;
    coupling_site_size = 4*9;
  } else {
    clover_site_size = (l->num_lattice_site_var*(l->num_lattice_site_var+1))/2;
    coupling_site_size = 4*l->num_lattice_site_var*l->num_lattice_site_var;
  }
  
  int its_boundary;
  if ( type ==_SCHWARZ ) {
    its_boundary = 2;
  } else {
    its_boundary = 1;
  }
  for ( mu=0; mu<4; mu++ ) {
    its *= (l->local_lattice[mu]+its_boundary);
  }
  
  int nls = (type==_ORDINARY)?l->num_inner_lattice_sites:2*l->num_lattice_sites-l->num_inner_lattice_sites;
  FREE( op->D, complex_PRECISION, coupling_site_size*nls );
  FREE( op->clover, complex_PRECISION, clover_site_size*l->num_inner_lattice_sites );
  if ( type == _SCHWARZ && l->depth == 0 && g.odd_even )
    FREE( op->oe_clover, complex_PRECISION, clover_site_size*l->num_inner_lattice_sites );
  FREE( op->index_table, int, its );
  FREE( op->neighbor_table, int, (l->depth==0?4:5)*l->num_inner_lattice_sites );
  FREE( op->backward_neighbor_table, int, (l->depth==0?4:5)*l->num_inner_lattice_sites );
  FREE( op->translation_table, int, l->num_inner_lattice_sites );
#ifdef SSE
  if ( l->depth == 0 ) {
    FREE( op->oe_clover_vectorized, PRECISION, 144*l->num_inner_lattice_sites );    
  }
#endif
  
  operator_PRECISION_free_projection_buffers( op, l );
  
  ghost_free_PRECISION( &(op->c), l );
  
  for ( mu=0; mu<4; mu++ ) {
    its = 1;
    for ( nu=0; nu<4; nu++ ) {
      if ( mu != nu ) {
        its *= l->local_lattice[nu];
      }
    }
    
    FREE( op->c.boundary_table[2*mu], int, its );
    if ( type == _SCHWARZ ) {
      FREE( op->c.boundary_table[2*mu+1], int, its );
      FREE( op->config_boundary_table[mu], int, its );
    } else {
      op->c.boundary_table[2*mu+1] = NULL;
    }
  }
}


void operator_PRECISION_define( operator_PRECISION_struct *op, level_struct *l ) {
  
   int i, mu, t, z, y, x, *it = op->index_table,
      ls[4], le[4], l_st[4], l_en[4], *dt = op->table_dim;
  
  for ( mu=0; mu<4; mu++ ) {
    dt[mu] = l->local_lattice[mu]+1;
    ls[mu] = 0;
    le[mu] = ls[mu] + l->local_lattice[mu];
    l_st[mu] = ls[mu];
    l_en[mu] = le[mu];
  }
  
  op->shift = 4+l->dirac_shift;
  
  // define index table
  // lexicographic inner cuboid and
  // lexicographic +T,+Z,+Y,+X boundaries
  i=0;
  // inner hyper cuboid
  for ( t=ls[T]; t<le[T]; t++ )
    for ( z=ls[Z]; z<le[Z]; z++ )
      for ( y=ls[Y]; y<le[Y]; y++ )
        for ( x=ls[X]; x<le[X]; x++ ) {
          it[ lex_index( t, z, y, x, dt ) ] = i; i++;
        }
  // boundaries (buffers)
  for ( mu=0; mu<4; mu++ ) {
    l_st[mu] = le[mu];
    l_en[mu] = le[mu]+1;
    
    for ( t=l_st[T]; t<l_en[T]; t++ )
      for ( z=l_st[Z]; z<l_en[Z]; z++ )
        for ( y=l_st[Y]; y<l_en[Y]; y++ )
          for ( x=l_st[X]; x<l_en[X]; x++ ) {
            it[ lex_index( t, z, y, x, dt ) ] = i; i++;
          }
          
    l_st[mu] = ls[mu];
    l_en[mu] = le[mu];
  }
  
  // define neighbor table (for the application of the entire operator),
  // negative inner boundary table (for communication),
  // translation table (for translation to lexicographical site ordnering)
  define_nt_bt_tt( op->neighbor_table, op->backward_neighbor_table, op->c.boundary_table, op->translation_table, it, dt, l );
}


void operator_PRECISION_test_routine( operator_PRECISION_struct *op, level_struct *l, struct Thread *threading ) {

/*********************************************************************************
* Checks for correctness of operator data layout by doing:
* - Applying D_W in double precision to a double vector.
* - Translates the same vector into PRECISION and apply D_W in PRECISION to this
*   vector and translate it back to double
* - Compare solutions ( Difference should be close to 0 ).
* If enabled, also tests odd even preconditioning.
*********************************************************************************/ 

  int ivs = l->inner_vector_size;
  double diff;
  
  vector_double vd1=NULL, vd2, vd3, vd4;
  vector_PRECISION vp1=NULL, vp2;

  PUBLIC_MALLOC( vd1, complex_double, 4*ivs );
  PUBLIC_MALLOC( vp1, complex_PRECISION, 2*ivs );

  vd2 = vd1+ivs; vd3 = vd2+ivs; vd4 = vd3 + ivs; vp2 = vp1 + ivs;

  START_LOCKED_MASTER(threading)
  
  vector_double_define_random( vd1, 0, l->inner_vector_size, l );
  apply_operator_double( vd2, vd1, &(g.p), l, no_threading );
  
  trans_PRECISION( vp1, vd1, op->translation_table, l, no_threading );
  apply_operator_PRECISION( vp2, vp1, &(l->p_PRECISION), l, no_threading );
  trans_back_PRECISION( vd3, vp2, op->translation_table, l, no_threading );
  
  vector_double_minus( vd4, vd3, vd2, 0, l->inner_vector_size, l );
  diff = global_norm_double( vd4, 0, ivs, l, no_threading )/global_norm_double( vd3, 0, ivs, l, no_threading );
  printf0("depth: 0, correctness of schwarz PRECISION Dirac operator: %le\n", diff );
  END_LOCKED_MASTER(threading)

  if(threading->n_core > 1) {
    apply_operator_PRECISION( vp2, vp1, &(l->p_PRECISION), l, threading );

    SYNC_MASTER_TO_ALL(threading)
    SYNC_CORES(threading)

    START_LOCKED_MASTER(threading)
    trans_back_PRECISION( vd3, vp2, op->translation_table, l, no_threading );
    
    vector_double_minus( vd4, vd3, vd2, 0, l->inner_vector_size, l );
    diff = global_norm_double( vd4, 0, ivs, l, no_threading )/global_norm_double( vd3, 0, ivs, l, no_threading );

    printf0("depth: 0, correctness of schwarz PRECISION Dirac operator with threading: %le\n", diff );
    END_LOCKED_MASTER(threading) 
  }    

  PUBLIC_FREE( vd1, complex_double, 4*ivs );
  PUBLIC_FREE( vp1, complex_PRECISION, 2*ivs );

  START_LOCKED_MASTER(threading)
  if ( g.method >=4 && g.odd_even )
    oddeven_PRECISION_test( l );
  END_LOCKED_MASTER(threading) 
}
