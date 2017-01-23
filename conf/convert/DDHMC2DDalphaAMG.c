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

enum{ T, Z, Y, X };


static inline int grid2index( int t, int z, int y, int x, int ls[4] ) {
  return x + ls[X]*( y + ls[Y]*( z + ls[Z]*t ) ) ;
}


int main ( int argc, char **argv ) {
  
  double *cnfg, plaq;
  int lsize[4], t, z, y, x, i;
  long long cnfg_size;
  FILE *fin, *fout;
  
  printf("converting %s to %s\n\n",argv[1],argv[2]);
  
  fin = fopen(argv[1],"rb");
  
  fread( &(lsize[T]), sizeof(int)   , 1, fin );
  fread( &(lsize[Z]), sizeof(int)   , 1, fin );
  fread( &(lsize[Y]), sizeof(int)   , 1, fin );
  fread( &(lsize[X]), sizeof(int)   , 1, fin );
  fread( &plaq      , sizeof(double), 1, fin );
  
  printf("%d %d %d %d lattice\n", lsize[T], lsize[Z], lsize[Y], lsize[X] );
  cnfg_size = ((((((long long)lsize[T])*lsize[Z])*lsize[Y])*lsize[X])*4*18);
  
  printf("cnfg_size: %ld bytes\n", sizeof(double)*cnfg_size );
    
  cnfg = (double*) malloc( sizeof(double)*cnfg_size );

  printf("reading...\n");
  
  
  
  for ( t=0; t<lsize[T]; t++ ) {
    printf("t=%d\n", t );
    for ( z=0; z<lsize[Z]; z++ )
      for ( y=0; y<lsize[Y]; y++ )
	for ( x=0; x<lsize[X]; x++ )
	  if ( (t+z+y+x)%2 == 1 ) {
	    
	    // +T
	    fread( &(cnfg[ 72*grid2index( t, z, y, x, lsize )                       + 18*T ]) , sizeof(double), 18, fin );
	    // -T
	    fread( &(cnfg[ 72*grid2index( (t+lsize[T]-1)%lsize[T], z, y, x, lsize ) + 18*T ]) , sizeof(double), 18, fin );
	    // +Z
	    fread( &(cnfg[ 72*grid2index( t, z, y, x, lsize )                       + 18*Z ]) , sizeof(double), 18, fin );
	    // -Z
	    fread( &(cnfg[ 72*grid2index( t, (z+lsize[Z]-1)%lsize[Z], y, x, lsize ) + 18*Z ]) , sizeof(double), 18, fin );
	    // +Y
	    fread( &(cnfg[ 72*grid2index( t, z, y, x, lsize )                       + 18*Y ]) , sizeof(double), 18, fin );
	    // -Y
	    fread( &(cnfg[ 72*grid2index( t, z, (y+lsize[Y]-1)%lsize[Y], x, lsize ) + 18*Y ]) , sizeof(double), 18, fin );
	    // +X
	    fread( &(cnfg[ 72*grid2index( t, z, y, x, lsize )                       + 18*X ]) , sizeof(double), 18, fin );
	    // -X
	    fread( &(cnfg[ 72*grid2index( t, z, y, (x+lsize[X]-1)%lsize[X], lsize ) + 18*X ]) , sizeof(double), 18, fin );
	    
	  }
  }
  
  fclose( fin );
  
  fout = fopen( argv[2], "w" );
  
  printf("writing to %s...\n", argv[2] );
  
  fwrite( lsize, sizeof(int)   , 4        , fout );
  fwrite( &plaq, sizeof(double), 1        , fout );
  fwrite(  cnfg, sizeof(double), cnfg_size, fout );
  
  fclose( fout );
  
  if ( argc > 3 ) {
    i=0;
    for ( t=0; t<lsize[T]; t++ )
      for ( z=0; z<lsize[Z]; z++ )
	for ( y=0; y<lsize[Y]; y++ )
	  for ( x=0; x<lsize[X]; x++ ) {
	      printf("site_number=%d, t=%d, z=%d, y=%d, x=%d.\n", i/72, t, z, y, x );
	      printf("dir: T\n");
	      printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", 
		    cnfg[i+0], cnfg[i+1], cnfg[i+2], cnfg[i+3], cnfg[i+4], cnfg[i+5]);
	      printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", 
		    cnfg[i+6], cnfg[i+7], cnfg[i+8], cnfg[i+9], cnfg[i+10], cnfg[i+11]);
	      printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", 
		    cnfg[i+12], cnfg[i+13], cnfg[i+14], cnfg[i+15], cnfg[i+16], cnfg[i+17]);
	      i+=18;
	      printf("dir: Z\n");
	      printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", 
		    cnfg[i+0], cnfg[i+1], cnfg[i+2], cnfg[i+3], cnfg[i+4], cnfg[i+5]);
	      printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", 
		    cnfg[i+6], cnfg[i+7], cnfg[i+8], cnfg[i+9], cnfg[i+10], cnfg[i+11]);
	      printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", 
		    cnfg[i+12], cnfg[i+13], cnfg[i+14], cnfg[i+15], cnfg[i+16], cnfg[i+17]);
	      i+=18;
	      printf("dir: Y\n");
	      printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", 
		    cnfg[i+0], cnfg[i+1], cnfg[i+2], cnfg[i+3], cnfg[i+4], cnfg[i+5]);
	      printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", 
		    cnfg[i+6], cnfg[i+7], cnfg[i+8], cnfg[i+9], cnfg[i+10], cnfg[i+11]);
	      printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", 
		    cnfg[i+12], cnfg[i+13], cnfg[i+14], cnfg[i+15], cnfg[i+16], cnfg[i+17]);
	      i+=18;
	      printf("dir: X\n");
	      printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", 
		    cnfg[i+0], cnfg[i+1], cnfg[i+2], cnfg[i+3], cnfg[i+4], cnfg[i+5]);
	      printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", 
		    cnfg[i+6], cnfg[i+7], cnfg[i+8], cnfg[i+9], cnfg[i+10], cnfg[i+11]);
	      printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", 
		    cnfg[i+12], cnfg[i+13], cnfg[i+14], cnfg[i+15], cnfg[i+16], cnfg[i+17]);
	      i+=18;
	      printf("------------------------------------------------------------\n");
	  }
  }
  
  free( cnfg );
  
  printf("done!\n");
  
  return 0;
}
