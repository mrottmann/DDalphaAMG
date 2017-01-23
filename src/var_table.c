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

void var_table_init( var_table *t ) {
  
  t->entry = NULL;
  t->iterator = NULL;
  t->p = NULL;
  t->p_end = NULL;
}


void var_table_insert( var_table *t, var_table_entry e ) {
  
  if ( t->entry == NULL ) {
    MALLOC( t->entry, var_table_entry, 1 );
    *(t->entry) = e;
    t->entry->next = NULL;
    
  } else {
    t->iterator = t->entry;
    while ( t->iterator->next != NULL )
      t->iterator = t->iterator->next;
    
    MALLOC( t->iterator->next, var_table_entry, 1 );
    *(t->iterator->next) = e;
    t->iterator->next->next = NULL;
  }
}


void var_table_free( var_table *t ) {
  
  if ( t->entry != NULL ) {
    t->iterator = t->entry;
    
    while( t->entry->next != NULL ) {
      t->iterator = t->entry->next;
      FREE( t->entry, var_table_entry, 1 );
      t->entry = t->iterator;
    }
    FREE( t->entry, var_table_entry, 1 );
    t->iterator = NULL;
  }
}


void scan_var( var_table *t, level_struct *l ) {
  
  t->iterator = t->entry;
  while( t->iterator != NULL  ) {
    if ( strcmp( t->iterator->name, t->scan_var ) != 0 )
      t->iterator = t->iterator->next;
    else
      break;
  }
  
  if ( t->iterator == NULL ) {
    error0("unable to scan variable \"%s\"\n", t->scan_var );
  } else {
    if ( strcmp( t->iterator->datatype,"int" ) == 0 )
      SCAN_VAR( t->iterator->pt, int, t->start_val, t->end_val, t->step_size, t->multiplicative, t->iterator->name, l );
    else
      SCAN_VAR( t->iterator->pt, double, t->start_val, t->end_val, t->step_size, t->multiplicative, t->iterator->name, l );
  }
  
  plot_table( t );
}


void new_plot_table_line( var_table *t ) {
  
  if ( t->p == NULL ) {
    MALLOC( t->p, plot_table_line, 1 );
    t->p_end = t->p;
    for ( int i=0; i<_NUM_OPTB; i++ )
      t->p->values[i] = 0.0;
    t->p->next = NULL;
  } else {
    ASSERT( t->p_end->next == NULL );
    MALLOC( t->p_end->next, plot_table_line, 1 );
    t->p_end = t->p_end->next;
    for ( int i=0; i<_NUM_OPTB; i++ )
      t->p_end->values[i] = 0.0;
    t->p_end->next = NULL;
  }
}


void plot_table( var_table *t ) {
  
  printf0( "\ntrckd val: \"%s\"\n", t->scan_var );
  printf0( "----------+----------+----------+----------+----------+----------+----------+----------+\n" );
  printf0( "trckd val | stp time | slv iter | slv time | crs iter | crs time | slv err  | CGNR err |\n" );
  printf0( "----------+----------+----------+----------+----------+----------+----------+----------+\n" );
  
  plot_table_line *p_cur;
  while ( t->p != NULL ) {
    for ( int i=0; i<_NUM_OPTB; i++ )
      printf0("%10.4le ", t->p->values[i] );
    printf0("\n");
    p_cur = t->p->next;
    FREE( t->p, plot_table_line, 1 );
    t->p = p_cur;
  }
  
  printf0( "----------------------------------------------------------------------------------------\n\n" );
}
