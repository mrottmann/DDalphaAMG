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

#if defined( SSE ) && defined( INTERPOLATION_OPERATOR_LAYOUT_OPTIMIZED_PRECISION )

void interpolation_PRECISION_alloc( level_struct *l ) {
  
  int k, n = l->num_eig_vect;
  
  MALLOC( l->is_PRECISION.eigenvalues, complex_PRECISION, n );
  MALLOC( l->is_PRECISION.test_vector, complex_PRECISION*, n );
  
#ifndef INTERPOLATION_SETUP_LAYOUT_OPTIMIZED_PRECISION  
  MALLOC( l->is_PRECISION.interpolation, complex_PRECISION*, n );
  l->is_PRECISION.interpolation[0] = NULL;
  MALLOC_HUGEPAGES( l->is_PRECISION.interpolation[0], complex_PRECISION, n*l->vector_size, 128 );
  for ( k=1; k<n; k++ )
    l->is_PRECISION.interpolation[k] = l->is_PRECISION.interpolation[0] + k*l->vector_size;
#endif
  // ghost shell is communicated in coarse_operator_setup, so we need size=vector_size, not inner_vector_size
  MALLOC_HUGEPAGES( l->is_PRECISION.operator, complex_PRECISION,
                    ((size_t)OPERATOR_COMPONENT_OFFSET_PRECISION)*((size_t)l->vector_size), 128 );

  l->is_PRECISION.test_vector[0] = NULL;
  MALLOC_HUGEPAGES( l->is_PRECISION.test_vector[0], complex_PRECISION, n*l->inner_vector_size, 128 );
  for ( k=1; k<n; k++ ) {
    l->is_PRECISION.test_vector[k] = l->is_PRECISION.test_vector[0] + k*l->inner_vector_size;
  }    
}


void interpolation_PRECISION_dummy_alloc( level_struct *l ) {
  
  MALLOC( l->is_PRECISION.test_vector, complex_PRECISION*, l->num_eig_vect );
  MALLOC( l->is_PRECISION.interpolation, complex_PRECISION*, l->num_eig_vect );  
}


void interpolation_PRECISION_dummy_free( level_struct *l ) {
  
  FREE( l->is_PRECISION.test_vector, complex_PRECISION*, l->num_eig_vect );
  FREE( l->is_PRECISION.interpolation, complex_PRECISION*, l->num_eig_vect );  
}


void interpolation_PRECISION_free( level_struct *l ) {
  
  int n = l->num_eig_vect;
  
  FREE_HUGEPAGES( l->is_PRECISION.test_vector[0], complex_PRECISION, n*l->inner_vector_size );
  FREE( l->is_PRECISION.eigenvalues, complex_PRECISION, n );
  FREE( l->is_PRECISION.test_vector, complex_PRECISION*, n );
#ifndef INTERPOLATION_SETUP_LAYOUT_OPTIMIZED_PRECISION  
  FREE_HUGEPAGES( l->is_PRECISION.interpolation[0], complex_PRECISION, n*l->vector_size );
  FREE( l->is_PRECISION.interpolation, complex_PRECISION*, n );
#endif
  FREE_HUGEPAGES( l->is_PRECISION.operator, complex_PRECISION, OPERATOR_COMPONENT_OFFSET_PRECISION*l->vector_size );
}


void swap8_PRECISION( PRECISION* data ) {
  
  int i;
  PRECISION tmp[8];
  
  for ( i=0; i<4; i++ ) {
    tmp[i] = data[2*i];
    tmp[i+4] = data[2*i+1];
  }
  
  for ( i=0; i<8; i++ ) {
    data[i] = tmp[i];
  }
}


void define_interpolation_PRECISION_operator( complex_PRECISION **interpolation, level_struct *l, struct Thread *threading ) {
  
  int j, num_eig_vect = l->num_eig_vect;
  complex_PRECISION *operator = l->is_PRECISION.operator;

  int start = threading->start_index[l->depth];
  int end = threading->end_index[l->depth];
      
  SYNC_CORES(threading)
  int offset = SIMD_LENGTH_PRECISION;
  for ( j=0; j<num_eig_vect; j+=offset ) {
    int j_end = j+offset;
    if(j_end > num_eig_vect)
      j_end = num_eig_vect;
    
    operator = l->is_PRECISION.operator + j*l->vector_size + start*offset;
    
    for ( int i=start; i<end; i+=offset/2 ) {
      __m128 data[offset];
      for ( int j2=j; j2<j_end; j2++ )
        data[j2-j] = _mm_load_ps((float *)(interpolation[j2]+i));
      for ( int j2=j_end; j2<j+offset; j2++ )
        data[j2-j] = _mm_setzero_ps();
      
      transpose_4_registers(data);

      for ( int k=0; k<offset; k++) {
        _mm_store_ps((float *)operator, data[k]);
        // operator type is complex, so offset only by SIMD_LENGTH_PRECISION over *two*
        operator += offset/2;
      }
    }
  }
  SYNC_CORES(threading)
}


void interpolate_PRECISION( vector_PRECISION phi, vector_PRECISION phi_c, level_struct *l, struct Thread *threading ) {
  
  PROF_PRECISION_START( _PR, threading );
  int i, j, k, k1, k2, num_aggregates = l->is_PRECISION.num_agg, num_eig_vect = l->num_eig_vect,
      num_parent_eig_vect = l->num_lattice_site_var/2, aggregate_sites = l->num_inner_lattice_sites / num_aggregates;
  complex_PRECISION *operator = l->is_PRECISION.operator, *phi_pt = phi,
                    *phi_c_pt = l->next_level->gs_PRECISION.transfer_buffer;
                    
  START_LOCKED_MASTER(threading)
  vector_PRECISION_distribute( phi_c_pt, phi_c, l->next_level );
  END_LOCKED_MASTER(threading)
  SYNC_HYPERTHREADS(threading)
  
  for ( i=threading->n_thread*threading->core + threading->thread; i<num_aggregates; i+=threading->n_core*threading->n_thread ) {
    phi_pt   = phi + i*2*num_parent_eig_vect*aggregate_sites;
    phi_c_pt = l->next_level->gs_PRECISION.transfer_buffer + i*2*num_eig_vect;
    float tmp_phi_c_re[2*OPERATOR_COMPONENT_OFFSET_float];
    float tmp_phi_c_im[2*OPERATOR_COMPONENT_OFFSET_float];
    __m128 zero =  _mm_setzero_ps();
    for ( j=0; j<2*OPERATOR_COMPONENT_OFFSET_PRECISION; j+=SIMD_LENGTH_PRECISION ) {
      _mm_store_ps(tmp_phi_c_re+j, zero);
      _mm_store_ps(tmp_phi_c_im+j, zero);
    }
    // copy phi_c into temporary
    for ( j=0; j<num_eig_vect; j++ ) {
      tmp_phi_c_re[j] = creal(phi_c_pt[j]);
      tmp_phi_c_im[j] = cimag(phi_c_pt[j]);
    }
    for ( j=0; j<num_eig_vect; j++ ) {
      tmp_phi_c_re[j+OPERATOR_COMPONENT_OFFSET_PRECISION] = creal(phi_c_pt[j+num_eig_vect]);
      tmp_phi_c_im[j+OPERATOR_COMPONENT_OFFSET_PRECISION] = cimag(phi_c_pt[j+num_eig_vect]);
    }

    int offset = SIMD_LENGTH_PRECISION;
    // loop over blocks of SIMD_LENGTH_PRECISION vectors
    for ( j=0; j<num_eig_vect; j+=offset ) {
      phi_pt   = phi + i*2*num_parent_eig_vect*aggregate_sites;
      operator = l->is_PRECISION.operator + j*l->vector_size + i*2*offset*num_parent_eig_vect*aggregate_sites;

      for ( k=0; k<aggregate_sites; k++ ) {
        // offset used for 2 components of gamma5-symmetry stuff
        int low_high_offset = 0;
        for ( k1=0; k1<2; k1++ ) {
          for ( k2=0; k2<num_parent_eig_vect; k2++ ) {
            __m128 phi_re = _mm_setzero_ps();
            __m128 phi_im = _mm_setzero_ps();
            
            __m128 operator_re = _mm_load_ps((float *)operator);
            __m128 operator_im = _mm_load_ps((float *)operator+offset);
            __m128 phi_c_re = _mm_load_ps(tmp_phi_c_re+j+low_high_offset);
            __m128 phi_c_im = _mm_load_ps(tmp_phi_c_im+j+low_high_offset);

            cfmadd(operator_re, operator_im, phi_c_re, phi_c_im, &phi_re, &phi_im);

            // skip to next real line of matrix
            operator += offset;
            // horizontal sum for phi
            __m128 tmp;
            tmp = _mm_add_ps( phi_re, _mm_movehl_ps( phi_re, phi_re ) );
            phi_re = _mm_add_ss( tmp, _mm_shuffle_ps( tmp, tmp, 1 ) );
            
            tmp = _mm_add_ps( phi_im, _mm_movehl_ps( phi_im, phi_im ) );
            phi_im = _mm_add_ss( tmp, _mm_shuffle_ps( tmp, tmp, 1 ) );
            
            __m128 tmp1; tmp1 = _mm_set1_ps(((float *)phi_pt)[0]);
            __m128 tmp2; tmp2 = _mm_set1_ps(((float *)phi_pt)[1]);
            phi_re = _mm_add_ps(phi_re, tmp1);
            phi_im = _mm_add_ps(phi_im, tmp2);
            _mm_store_ss( (float*)phi_pt, phi_re );
            _mm_store_ss( ((float*)phi_pt)+1, phi_im );
            phi_pt++;
          }
          low_high_offset = OPERATOR_COMPONENT_OFFSET_PRECISION;
        }
      }
    }
  }
    
  PROF_PRECISION_STOP( _PR, 1, threading );

  SYNC_HYPERTHREADS(threading)
}



void interpolate3_PRECISION( vector_PRECISION phi, vector_PRECISION phi_c, level_struct *l, struct Thread *threading ) {
  
  PROF_PRECISION_START( _PR, threading );
  int i, j, k, k1, k2, num_aggregates = l->is_PRECISION.num_agg, num_eig_vect = l->num_eig_vect,
      num_parent_eig_vect = l->num_lattice_site_var/2, aggregate_sites = l->num_inner_lattice_sites / num_aggregates;
  complex_PRECISION *operator = l->is_PRECISION.operator, *phi_pt = phi,
                    *phi_c_pt = l->next_level->gs_PRECISION.transfer_buffer;
  
  START_LOCKED_MASTER(threading)
  vector_PRECISION_distribute( phi_c_pt, phi_c, l->next_level );
  END_LOCKED_MASTER(threading)
  SYNC_HYPERTHREADS(threading)
  
  for ( i=threading->n_thread*threading->core + threading->thread; i<num_aggregates; i+=threading->n_core*threading->n_thread ) {
    phi_pt   = phi + i*2*num_parent_eig_vect*aggregate_sites;
    phi_c_pt = l->next_level->gs_PRECISION.transfer_buffer + i*2*num_eig_vect;
    float tmp_phi_c_re[2*OPERATOR_COMPONENT_OFFSET_float];
    float tmp_phi_c_im[2*OPERATOR_COMPONENT_OFFSET_float];
    __m128 zero =  _mm_setzero_ps();
    for ( j=0; j<2*OPERATOR_COMPONENT_OFFSET_PRECISION; j+=SIMD_LENGTH_PRECISION ) {
      _mm_store_ps(tmp_phi_c_re+j, zero);
      _mm_store_ps(tmp_phi_c_im+j, zero);
    }
    // copy phi_c into temporary
    for ( j=0; j<num_eig_vect; j++ ) {
      tmp_phi_c_re[j] = creal(phi_c_pt[j]);
      tmp_phi_c_im[j] = cimag(phi_c_pt[j]);
    }
    for ( j=0; j<num_eig_vect; j++ ) {
      tmp_phi_c_re[j+OPERATOR_COMPONENT_OFFSET_PRECISION] = creal(phi_c_pt[j+num_eig_vect]);
      tmp_phi_c_im[j+OPERATOR_COMPONENT_OFFSET_PRECISION] = cimag(phi_c_pt[j+num_eig_vect]);
    }

    int offset = SIMD_LENGTH_PRECISION;
    // loop over blocks of SIMD_LENGTH_PRECISION vectors
    for ( j=0; j<num_eig_vect; j+=offset ) {
      phi_pt   = phi + i*2*num_parent_eig_vect*aggregate_sites;
      operator = l->is_PRECISION.operator + j*l->vector_size + i*2*offset*num_parent_eig_vect*aggregate_sites;

      for ( k=0; k<aggregate_sites; k++ ) {
        // offset used for 2 components of gamma5-symmetry stuff
        int low_high_offset = 0;
        for ( k1=0; k1<2; k1++ ) {
          for ( k2=0; k2<num_parent_eig_vect; k2++ ) {
            __m128 phi_re = _mm_setzero_ps();
            __m128 phi_im = _mm_setzero_ps();

            __m128 operator_re = _mm_load_ps((float *)operator);
            __m128 operator_im = _mm_load_ps((float *)operator+offset);
            __m128 phi_c_re = _mm_load_ps(tmp_phi_c_re+j+low_high_offset);
            __m128 phi_c_im = _mm_load_ps(tmp_phi_c_im+j+low_high_offset);

            cfmadd(operator_re, operator_im, phi_c_re, phi_c_im, &phi_re, &phi_im);

            // skip to next real line of matrix
            operator += offset;
            // horizontal sum for phi            
            __m128 tmp;
            tmp = _mm_add_ps( phi_re, _mm_movehl_ps( phi_re, phi_re ) );
            phi_re = _mm_add_ss( tmp, _mm_shuffle_ps( tmp, tmp, 1 ) );
            
            tmp = _mm_add_ps( phi_im, _mm_movehl_ps( phi_im, phi_im ) );
            phi_im = _mm_add_ss( tmp, _mm_shuffle_ps( tmp, tmp, 1 ) );
            
            if ( j!= 0 ) {
              __m128 tmp1; tmp1 = _mm_set1_ps(((float *)phi_pt)[0]);
              __m128 tmp2; tmp2 = _mm_set1_ps(((float *)phi_pt)[1]);
              phi_re = _mm_add_ps(phi_re, tmp1);
              phi_im = _mm_add_ps(phi_im, tmp2);
            }
            _mm_store_ss( (float*)phi_pt, phi_re );
            _mm_store_ss( ((float*)phi_pt)+1, phi_im );
            phi_pt++;
          }
          low_high_offset = OPERATOR_COMPONENT_OFFSET_PRECISION;
        }
      }
    }
  }
    
  PROF_PRECISION_STOP( _PR, 1, threading );

  SYNC_HYPERTHREADS(threading)
}


void restrict_PRECISION( vector_PRECISION phi_c, vector_PRECISION phi, level_struct *l, struct Thread *threading ) {
  
  SYNC_CORES(threading)
  SYNC_HYPERTHREADS(threading)

  PROF_PRECISION_START( _PR, threading );
  int i, j, k, k1, k2, num_aggregates = l->is_PRECISION.num_agg, num_eig_vect = l->num_eig_vect,
      num_parent_eig_vect = l->num_lattice_site_var/2, aggregate_sites = l->num_inner_lattice_sites / num_aggregates;
  complex_PRECISION *operator = l->is_PRECISION.operator, *phi_pt = phi,
                    *phi_c_pt = l->next_level->gs_PRECISION.transfer_buffer;

  for ( i=threading->n_thread*threading->core + threading->thread; i<num_aggregates; i+=threading->n_core*threading->n_thread ) {
    
    phi_pt   = phi + i*2*num_parent_eig_vect*aggregate_sites;
    phi_c_pt = l->next_level->gs_PRECISION.transfer_buffer + i*2*num_eig_vect;

    int offset = SIMD_LENGTH_PRECISION;
    // loop over blocks of SIMD_LENGTH_PRECISION vectors
    for ( j=0; j<num_eig_vect; j+=offset ) {
      phi_pt   = phi + i*2*num_parent_eig_vect*aggregate_sites;
      phi_c_pt = l->next_level->gs_PRECISION.transfer_buffer + i*2*num_eig_vect;
      operator = l->is_PRECISION.operator + j*l->vector_size + i*2*offset*num_parent_eig_vect*aggregate_sites;

      // temporary, so we can used aligned load/store, and don't have to mess around with deinterleaving
      // complex components and masking
      // factor 2 is for low/high (refers to spin components being split up for preserving gamma5-symmetry of coarse operator)
      float tmp_phi_c_re[2*offset];
      float tmp_phi_c_im[2*offset];
      __m128 zero =  _mm_setzero_ps();
      for ( k1=0; k1<2*offset; k1+=offset ) {
        _mm_store_ps(tmp_phi_c_re+k1, zero);
        _mm_store_ps(tmp_phi_c_im+k1, zero);
      }

      for ( k=0; k<aggregate_sites; k++ ) {
        // offset used for 2 components of gamma5-symmetry stuff
        int low_high_offset = 0;
        for ( k1=0; k1<2; k1++ ) {
          for ( k2=0; k2<num_parent_eig_vect; k2++ ) {
            // phi is the same for all eigenvectors -> broadcast
            __m128 phi_re = _mm_set1_ps(((float *)phi_pt)[0]);
            __m128 phi_im = _mm_set1_ps(((float *)phi_pt)[1]);

            __m128 operator_re = _mm_load_ps((float *)operator);
            __m128 operator_im = _mm_load_ps((float *)operator+offset);
            __m128 phi_c_re = _mm_load_ps(tmp_phi_c_re+low_high_offset);
            __m128 phi_c_im = _mm_load_ps(tmp_phi_c_im+low_high_offset);

            cfmadd_conj(operator_re, operator_im, phi_re, phi_im, &phi_c_re, &phi_c_im);

            _mm_store_ps(tmp_phi_c_re+low_high_offset, phi_c_re);
            _mm_store_ps(tmp_phi_c_im+low_high_offset, phi_c_im);
            // skip to next real line of matrix
            operator += offset;
            phi_pt++;
          }
          low_high_offset = offset;
        }
      }

      for ( int m=0; m<offset; m++ ) {
        if ( m+j >= num_eig_vect ) break;
        ((float*)(phi_c_pt+j+m))[0] = tmp_phi_c_re[m];
        ((float*)(phi_c_pt+j+m))[1] = tmp_phi_c_im[m];
      }
      
      for ( int m=0; m<offset; m++ ) {
        if ( m+j >= num_eig_vect ) break;
        ((float*)(phi_c_pt+num_eig_vect+j+m))[0] = tmp_phi_c_re[m+offset];
        ((float*)(phi_c_pt+num_eig_vect+j+m))[1] = tmp_phi_c_im[m+offset];
      }
    }
  }
  
  SYNC_HYPERTHREADS(threading)
  START_LOCKED_MASTER(threading)
  vector_PRECISION_gather( phi_c, l->next_level->gs_PRECISION.transfer_buffer, l->next_level );
  END_LOCKED_MASTER(threading)
  PROF_PRECISION_STOP( _PR, 1, threading );
}

#endif // defined( SSE ) && defined( INTERPOLATION_OPERATOR_LAYOUT_OPTIMIZED_PRECISION )