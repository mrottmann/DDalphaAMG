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


/*
 * A converter from single file configuration to multifiles for DDalphaAMG
 *   2016 Oct. 3, Issaku Kanamori
 *
 * compile with
 *   -std=c99
 * option. 
 *
 * if you want run it on a machine with opposite endian to that of the
 * configuration, add
 *   -DMACHINE_ENDIAN_IS_OPPOSITE
 * as well.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <string.h>


//#define MACHINE_ENDIAN_IS_OPPOSITE
#define Ndim 4
#define site_dof 72  // 2 * Nc*Nc *Ndim  ( 2 for complex)

#ifdef MACHINE_ENDIAN_IS_OPPOSITE
const bool machine_endian_is_opposite=1;
#else
const bool machine_endian_is_opposite=0;
#endif

void swap_endian4(char *c){
  char buffer[4];
  buffer[0]=c[0];
  buffer[1]=c[1];
  buffer[2]=c[2];
  buffer[3]=c[3];

  c[0]=buffer[3];
  c[1]=buffer[2];
  c[2]=buffer[1];
  c[3]=buffer[0];
  return;
}


void swap_endian8(char *c){
  char buffer[8];
  buffer[0]=c[0];
  buffer[1]=c[1];
  buffer[2]=c[2];
  buffer[3]=c[3];
  buffer[4]=c[4];
  buffer[5]=c[5];
  buffer[6]=c[6];
  buffer[7]=c[7];

  c[0]=buffer[7];
  c[1]=buffer[6];
  c[2]=buffer[5];
  c[3]=buffer[4];
  c[4]=buffer[3];
  c[5]=buffer[2];
  c[6]=buffer[1];
  c[7]=buffer[0];
  return;
}

struct header_type{
  int lsize[Ndim];
  double plaq;
};


//
// read header from the file
//
void read_header(struct header_type *h, FILE *f, const int *L){

  if(machine_endian_is_opposite)
    printf("assumes the machine endian is opposite to the configuration\n");

  assert(sizeof(int)==4);
 
  if( fread(h->lsize, sizeof(int), Ndim, f) < Ndim){
    printf("bad header (lsize), or file is too small\n");
    fclose(f);
    exit(EXIT_FAILURE);
  }
  for(int i=0; i<Ndim; i++){
    int ltmp;
    memcpy(&ltmp, h->lsize+i, sizeof(int));
    //int ltmp=(h->lsize)[i];
    if(machine_endian_is_opposite)
      swap_endian4((char*)&ltmp);

    if(ltmp!=L[i]){
      printf("lattice size does not match: i=%d, given=%d, header=%d\n",i,L[i], ltmp);
      exit(EXIT_FAILURE);
    }
  }
  if( fread(&(h->plaq), sizeof(double), 1, f) <1 ){
    printf("bad header (plaq), or file is too small\n");
    exit(EXIT_FAILURE);
  }
  double plaq_tmp=h->plaq;
  if(machine_endian_is_opposite)
    swap_endian8((char*) &plaq_tmp);

  printf("plaquette [0,3]: %f\n", plaq_tmp);
  printf("plaquette [0,1]: %f\n", plaq_tmp/3.0);
}


//
// write header to file
//
void write_header(const struct header_type *h, FILE *f){
  fwrite(h->lsize, sizeof(int), 4, f);
  fwrite(&(h->plaq), sizeof(double), 1, f);  
}


//
// read configuraion
//
void read_conf(const int *L, const char *infile, double *u, struct header_type *h){

  // open
  FILE *fin=0;
  fin=fopen(infile, "rb");
  if(fin ==0){
    printf("error in openning file: %s\n", infile);
    exit(EXIT_FAILURE);
  }

  // header
  read_header(h,fin,L);

  // body
  long dof=site_dof;
  dof*=L[0];
  dof*=L[1];
  dof*=L[2];
  dof*=L[3];
  if( fread(u,sizeof(double), dof, fin) < dof ){
    printf("seems input file is too small\n");
    exit(EXIT_FAILURE);
  };

  // done
  fclose(fin);
  printf("reading, done\n");
}


//
// save configuration for each node
//
void  save_conf(const int *L, const int *P, const char *ofile, const double *u, struct header_type *h){

  const int localsize_t=L[0]/P[0];
  const int localsize_z=L[1]/P[1];
  const int localsize_y=L[2]/P[2];
  const int localsize_x=L[3]/P[3];

  const long Lt=L[0];
  const long Lz=L[1];
  const long Ly=L[2];
  const long Lx=L[3];

  printf("local lattice size= %d %d %d %d\n", 
	 localsize_t, localsize_z, localsize_y, localsize_x);
  long this_dof=site_dof;
  this_dof*=localsize_t*localsize_z;
  this_dof*=localsize_y*localsize_x;
  double *this_u=(double *)malloc(sizeof(double)*this_dof);
  assert(this_u!=0);
  if(this_u == 0){
    printf("alocation error for this_u\n");
    exit(EXIT_FAILURE);
  }

  for(int pt=0; pt<P[0]; pt++){
    for(int pz=0; pz<P[1]; pz++){
      for(int py=0; py<P[2]; py++){
	for(int px=0; px<P[3]; px++){

	  // make a copy for one node
	  int local_index=0;
	  for(int t=0; t<localsize_t; t++){
	    for(int z=0; z<localsize_z; z++){
	      for(int y=0; y<localsize_y; y++){
		for(int x=0; x<localsize_x; x++){

		  long int gt=localsize_t*pt+t;
		  long int gz=localsize_z*pz+z;
		  long int gy=localsize_y*py+y;
		  long int gx=localsize_x*px+x;
		  long int global_index=(gt*Lz*Ly*Lx + gz*Ly*Lx + gy*Lx +gx);
		  global_index*=site_dof;
		  memcpy(this_u+local_index, u+global_index, sizeof(double)*site_dof);
		  local_index+=site_dof;
	  }}}}

	  // write to file
	  char postfix[128];
	  char filename[1000];
	  sprintf(postfix, ".pt%dpz%dpy%dpx%d",pt,pz,py,px);
	  sprintf(filename, "%s%s",ofile, postfix);
	  printf("writing to %s\n", filename);
	  fflush(0);
	  FILE *fout=fopen(filename, "wb");
	  if(fout == 0){
	    printf("error in opening output file: %s\n", filename);
	    exit(EXIT_FAILURE);
	  }
	  write_header(h,fout);
	  fwrite(this_u, sizeof(double), this_dof, fout);
	  fclose(fout);
	  
	}}}}

  // done
  free(this_u);
  this_u=0;
}
 

//
// main
//
int main(int argc, char **argv){

  if(argc < 11){
    printf("Usage: %s Lt Lz Ly Lx Pt Pz Py Px inputfile outputfile\n",argv[0]);
    printf("  L? : global lattice size in each direction\n");
    printf("  P? : number of processors in each direction\n");
    printf("   total number of processrors is Pt x Pz x Py x Px\n");
    return 0;
  }

  assert(Ndim==4);
  int L[Ndim], P[Ndim];
  const char *infile;
  const char *outfile;

  int count=1;
  for(int i=0; i<Ndim; i++){
    L[i] = atoi(argv[count]);
    count++;
    assert(L[i]>0);
  }
  printf("(Lt,Lz,Ly,Lx)=(%d,%d,%d,%d)\n",L[0],L[1],L[2],L[3]);

  for(int i=0; i<Ndim; i++){
    P[i] = atoi(argv[count]);
    count++;
    assert(P[i]>0);
  }
  printf("(Pt,Pz,Py,Px)=(%d,%d,%d,%d)\n",P[0],P[1],P[2],P[3]);
  printf(" number of total porcessrors in taget: %d\n", P[0]*P[1]*P[2]*P[3]);

  infile=argv[count];
  count++;
  outfile=argv[count];

  printf(" input file: %s\n", infile);
  printf("output file: %s\n", outfile);

  // sanity check
  for(int i=0; i<Ndim; i++){
    if(L[i] % P[i] != 0){
      printf("bad partitioning: i=%d, L[i]=%d, P[i]=%d\n", i, L[i], P[i]);
      exit(EXIT_FAILURE);
    }
  }

  // read configuration
  long dof=site_dof;
  dof*=L[0];
  dof*=L[1];
  dof*=L[2];
  dof*=L[3];
  double *u =(double *) malloc(sizeof(double)*dof);
  assert(u!=0);
  if(u==0){
    printf("allocation error for u\n");
    exit(EXIT_FAILURE);
  }
  printf("dof   = %ld\n",dof);
  printf("u     = %p\n",u);
  printf("u+dof = %p\n",u+dof);

  struct header_type header;
  read_conf(L,infile, u, &header);

  // write configurations for each processror
  save_conf(L,P, outfile, u, &header);
  
  // done
  free(u);
  u=0;
  printf("done!\n");

}
