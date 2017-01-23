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
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#define MALLOC( variable, kind, length ) do{ if ( variable != NULL ) { \
printf("malloc of \"%s\" failed: pointer is not NULL (%s:%d).\n", #variable, __FILE__, __LINE__ ); \
assert( 0 ); } \
variable = (kind*) malloc( sizeof(kind) * (length) ); \
if ( variable == NULL ) { \
printf("malloc of \"%s\" failed: no memory allocated (%s:%d).\n", #variable, __FILE__, __LINE__ ); \
assert( 0 ); } }while(0)

#define FREE( variable ) do{ if ( variable == NULL ) { \
printf("free of \"%s\" failed: pointer is already NULL (%s:%d).\n", #variable, __FILE__, __LINE__ ); \
assert( 0 ); } \
free( variable ); variable = NULL; }while(0)

#define NORM(x) sqrt( creal(x[0])*creal(x[0])+cimag(x[0])*cimag(x[0]) + \
		      creal(x[1])*creal(x[1])+cimag(x[1])*cimag(x[1]) + \
		      creal(x[2])*creal(x[2])+cimag(x[2])*cimag(x[2]) )
		      
#define DET(x) x[0][0]*x[1][1]*x[2][2] + \
	       x[0][1]*x[1][2]*x[2][0] + \
	       x[0][2]*x[1][0]*x[2][1] - \
	       x[0][0]*x[1][2]*x[2][1] - \
	       x[0][1]*x[1][0]*x[2][2] - \
	       x[0][2]*x[1][1]*x[2][0]
	       

enum{ T, Z, Y, X };

typedef double complex complex_double;


double rand_signum( void ) {
  double s = ((double)rand())/((double)RAND_MAX)*2.0 - 1.0;
  if ( s >= 0 ) s = 1;
  else s = -1;
  return s;
}



void rand_su3( double *U ) {
  
  complex_double SU3[3][3], beta;
  double alpha;
  
  SU3[0][0] = ( ((double)rand())/((double)RAND_MAX)* + 10000.0) //* rand_signum()  
	    + ( ((double)rand())/((double)RAND_MAX)*2.0 - 1.0 )*_Complex_I;
	    
  SU3[0][1] = ((double)rand())/((double)RAND_MAX)*2.0 - 1.0
            + ( ((double)rand())/((double)RAND_MAX)*2.0 - 1.0 )*_Complex_I;
  
  SU3[0][2] = ((double)rand())/((double)RAND_MAX)*2.0 - 1.0
            + ( ((double)rand())/((double)RAND_MAX)*2.0 - 1.0 )*_Complex_I;
	    
  SU3[1][0] = ((double)rand())/((double)RAND_MAX)*2.0 - 1.0
            + ( ((double)rand())/((double)RAND_MAX)*2.0 - 1.0 )*_Complex_I;
	    
  SU3[1][1] = ( ((double)rand())/((double)RAND_MAX)* + 10000.0) //* rand_signum()  
	    + ( ((double)rand())/((double)RAND_MAX)*2.0 - 1.0 )*_Complex_I;
	    
  SU3[1][2] = ((double)rand())/((double)RAND_MAX)*2.0 - 1.0
            + ( ((double)rand())/((double)RAND_MAX)*2.0 - 1.0 )*_Complex_I;
   
  alpha = NORM( SU3[0] );
  
  SU3[0][0] /= alpha;
  SU3[0][1] /= alpha;
  SU3[0][2] /= alpha;
  
  beta = conj( SU3[0][0] )*SU3[1][0]
	+ conj( SU3[0][1] )*SU3[1][1]
	+ conj( SU3[0][2] )*SU3[1][2];
	
  SU3[1][0] -= beta*SU3[0][0];
  SU3[1][1] -= beta*SU3[0][1];
  SU3[1][2] -= beta*SU3[0][2];
  
  alpha = NORM( SU3[1] );
  
  SU3[1][0] /= alpha;
  SU3[1][1] /= alpha;
  SU3[1][2] /= alpha;
	       
  SU3[2][0] = conj( SU3[0][1]*SU3[1][2] - SU3[0][2]*SU3[1][1] );
  SU3[2][1] = conj( SU3[0][2]*SU3[1][0] - SU3[0][0]*SU3[1][2] );
  SU3[2][2] = conj( SU3[0][0]*SU3[1][1] - SU3[0][1]*SU3[1][0] );
  
  alpha = NORM( SU3[2] );
  
  SU3[2][0] /= alpha;
  SU3[2][1] /= alpha;
  SU3[2][2] /= alpha;
  
  U[ 0] = 1;     U[ 1] = 0;
  U[ 2] = 0;     U[ 3] = 0;
  U[ 4] = 0;     U[ 5] = 0;
  
  U[ 6] = 0;     U[ 7] = 0;
  U[ 8] = 1;     U[ 9] = 0;
  U[10] = 0;     U[11] = 0;
  
  U[12] = 0;     U[13] = 0;
  U[14] = 0;     U[15] = 0;
  U[16] = 1;     U[17] = 0;
                               
}                              

                               
int main ( int argc, char **argv ) {
                               
  int mu, N[4], k = 5;
  long long int i = 0, n = 0;         
  double plaq = 0.0, U[18];
  char s[100];                 
  FILE *fout = NULL;
  
  srand( time ( 0 ) );
  
  for ( mu=0; mu<4; mu++ )
    N[mu] = atoi(argv[mu+1]);
  
  n = ((long long int)N[T])*((long long int)N[Z])*((long long int)N[Y])*((long long int)N[X]);
  n *= 72;
  
  sprintf( s, "%dx%dx%dx%d_unit", (int)N[T], (int)N[Z], (int)N[Y], (int)N[X] );
  
  assert( ( fout = fopen( s, "wb" ) ) != NULL );
  
  fwrite( N, sizeof(int), 4, fout );
  fwrite( &plaq, sizeof(double), 1, fout );
  
  while ( i < n ) {
    
    if ( ((double)i/(double)n)*100.0 >= k ) {
      printf("%d%%...\n", k );
      k += 5;
    }
    
    rand_su3( U );
    
    fwrite( U, sizeof(double), 18, fout ); fflush(0);
    
    i+=18;
  }
  printf("100%%...\n");
  
  printf("configuration %s written!\n", s );
  
  fclose( fout );
  
  return 0;
  
}