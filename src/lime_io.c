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

#ifdef HAVE_LIME

typedef struct lime_fileinfo {
  int precision, flavours, l[4], spin, colour;
  double plaquette;
  n_uint64_t data_size, offset;
} lime_fileinfo;

/*
 *  LIME functions
 *                                                                                   
 * In DDalphaAMG format:     
 ** t slowest running index	         
 ** x fastest running index            
 ** all positive directions                
 ** ordering: +T,+Z,+Y,+X         
 ** SU3 matrices stored in row major format             
 **        
 * In LIME format:      
 ** t slowest running index                
 ** x fastest running index                
 ** all positive directions                
 ** ordering: +X,+Y,+Z,+T ---> "[3-mu]"  
 ** SU3 matrices stored in row major format             
 **/

void swap_order(double* data_pt, int tot_size, int fixed_size, int n_ind_to_swap ) {
  
  ASSERT((tot_size%(n_ind_to_swap*fixed_size))==0);
  double * tmp;
  int i, j, k, bin_size=(n_ind_to_swap*fixed_size), bins=tot_size/bin_size;
  tmp=(double*)malloc(fixed_size*(n_ind_to_swap/2)*sizeof(double));
  
  for(i=0; i<bins; i++) {
    for(j=0; j< n_ind_to_swap/2; j++)
      for(k=0; k<fixed_size; k++) {
        tmp[j*fixed_size + k]=data_pt[i*bin_size + j*fixed_size + k];
        data_pt[i*bin_size + j*fixed_size + k] = data_pt[(i+1)*bin_size-(j+1)*fixed_size+k];
      }
      for(j=0; j< n_ind_to_swap/2; j++)
        for(k=0; k<fixed_size; k++) 
          data_pt[(i+1)*bin_size - (j+1)*fixed_size + k] = tmp[j*fixed_size + k]; 
  }
  free(tmp);
}

void swap_spin_in_conf(double* data_pt, int size) {
  swap_order(data_pt, size, 18, 4);
}

void swap_spin_in_vector(double* data_pt, int size) {
  swap_order(data_pt, size, 6, 4);
}

char* lime_getParam(char token[],char* params,int len, int lastchar)
{
  int i,token_len=strlen(token);
  
  for(i=0;i<len-token_len;i++)
  {
    if(memcmp(token,params+i,token_len)==0)
    {
      i+=token_len;
      *(strchr(params+i,lastchar))='\0';
      break;
    }
  }
  return params+i;
}

void lime_read_info(FILE **fin, lime_fileinfo *lime, char* binary_data) {
  LimeReader *limereader;
  char *lime_type,*lime_data;
  
  ASSERT((limereader = limeCreateReader(*fin)));
  
  lime->spin=4;
  lime->colour=3;
  
  while(limeReaderNextRecord(limereader) != LIME_EOF)
  {
    lime_type = limeReaderType(limereader);
    if(strcmp(lime_type, binary_data)==0)
    {
      break;
    }
    if(strcmp(lime_type,"ildg-format")==0)
    {
      lime->data_size = limeReaderBytes(limereader);
      lime_data = (char * )malloc(lime->data_size);
      limeReaderReadData((void *)lime_data,&(lime->data_size), limereader);
      sscanf(lime_getParam("<precision>",lime_data, lime->data_size, '<'),"%i",&(lime->precision));
      sscanf(lime_getParam("<lx>",lime_data, lime->data_size, '<'),"%i",&(lime->l[X]));
      sscanf(lime_getParam("<ly>",lime_data, lime->data_size, '<'),"%i",&(lime->l[Y]));
      sscanf(lime_getParam("<lz>",lime_data, lime->data_size, '<'),"%i",&(lime->l[Z]));
      sscanf(lime_getParam("<lt>",lime_data, lime->data_size, '<'),"%i",&(lime->l[T]));
      free(lime_data);
    }
    if(strcmp(lime_type,"etmc-source-format")==0 || strcmp(lime_type,"etmc-propagator-format")==0)
    {
      lime->data_size = limeReaderBytes(limereader);
      lime_data = (char * )malloc(lime->data_size);
      limeReaderReadData((void *)lime_data,&(lime->data_size), limereader);
      sscanf(lime_getParam("<precision>",lime_data, lime->data_size, '<'),"%i",&(lime->precision));
      sscanf(lime_getParam("<flavours>",lime_data, lime->data_size, '<'),"%i",&(lime->flavours));
      sscanf(lime_getParam("<lx>",lime_data, lime->data_size, '<'),"%i",&(lime->l[X]));
      sscanf(lime_getParam("<ly>",lime_data, lime->data_size, '<'),"%i",&(lime->l[Y]));
      sscanf(lime_getParam("<lz>",lime_data, lime->data_size, '<'),"%i",&(lime->l[Z]));
      sscanf(lime_getParam("<lt>",lime_data, lime->data_size, '<'),"%i",&(lime->l[T]));
      sscanf(lime_getParam("<spin>",lime_data, lime->data_size, '<'),"%i",&(lime->spin));
      sscanf(lime_getParam("<colour>",lime_data, lime->data_size, '<'),"%i",&(lime->colour));
      free(lime_data);
    }
    if(strcmp(lime_type,"xlf-info")==0)
    {
      lime->data_size = limeReaderBytes(limereader);
      lime_data = (char * )malloc(lime->data_size);
      limeReaderReadData((void *)lime_data,&(lime->data_size), limereader);
      sscanf(lime_getParam("plaquette =",lime_data, lime->data_size, '\n'),"%lf",&(lime->plaquette));
      free(lime_data);
    }
  }
  
  lime->data_size = limeReaderBytes(limereader);
  lime->offset = ftell(*fin);
  limeDestroyReader(limereader);
}

void lime_check_info(lime_fileinfo lime, n_uint64_t data_size_per_site) {
  int vol=1, mu;
  for ( mu=0; mu<4; mu++ ){
    ASSERT( lime.l[mu] == g.global_lattice[0][mu] );
    vol*=lime.l[mu];
  }
  ASSERT(lime.precision==64 || lime.precision==32);
  ASSERT(lime.data_size==(n_uint64_t) vol*data_size_per_site
  *((n_uint64_t)(lime.precision==64?sizeof(complex_double):sizeof(complex_float))));
}

void lime_write_info(FILE **fout, char* binary_data, n_uint64_t data_size_per_site, n_uint64_t *offset) {
  
  LimeWriter *limewriter;
  LimeRecordHeader *limeheader = NULL;
  n_uint64_t message_length;
  char tmp_string[STRINGLENGTH];
  
  ASSERT((limewriter = limeCreateWriter(*fout))!=(LimeWriter*)NULL);
  
  sprintf(tmp_string, "Vector_by_DDalphaAMG");
  message_length=(n_uint64_t) strlen(tmp_string);
  limeheader = limeCreateHeader(1, 1, "vector-type", message_length);
  
  ASSERT(limeheader != (LimeRecordHeader*)NULL);
  ASSERT(limeWriteRecordHeader(limeheader, limewriter)>=0);
  ASSERT(limeWriteRecordData(tmp_string, &message_length, limewriter)>=0);
  limeDestroyHeader(limeheader);
  
  sprintf(tmp_string, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<etmcFormat>\n\t<field>diracFermion</field>\n\t<precision>64</precision>\n\t<flavours>1</flavours>\n\t<lx>%d</lx>\n\t<ly>%d</ly>\n\t<lz>%d</lz>\n\t<lt>%d</lt>\n\t<spin>4</spin>\n\t<colour>3</colour>\n</etmcFormat>", g.global_lattice[0][X], g.global_lattice[0][Y], g.global_lattice[0][Z], g.global_lattice[0][T]);
  message_length=(n_uint64_t) strlen(tmp_string);
  limeheader = limeCreateHeader(1, 1, "etmc-propagator-format", message_length);
  
  ASSERT(limeheader != (LimeRecordHeader*)NULL);
  ASSERT(limeWriteRecordHeader(limeheader, limewriter)>=0);
  ASSERT(limeWriteRecordData(tmp_string, &message_length, limewriter)>=0);
  limeDestroyHeader(limeheader);
  
  #ifndef HAVE_TM
  sprintf(tmp_string, "<header>\n\tclifford basis: %s\n\tm0: %.14lf\n\tcsw: %.14lf\n\tclov plaq: %.14lf\n\thopp plaq: %.14lf\n\tcomputed plaq: %.14lf\n\tclov conf name: %s\n\thopp conf name: %s\n\tX: %d\n\tY: %d\n\tZ: %d\n\tT: %d\n\tX local: %d\n\tY local: %d\n\tZ local: %d\n\tT local: %d\n\tsetup iter: %d\n\tpost smoothing iter: %d\n\tblock iter: %d\n</header>", CLIFFORD_BASIS, g.solve_m0, g.csw, g.plaq_clov, g.plaq_hopp, g.plaq, g.in_clov, g.in, g.global_lattice[0][X], g.global_lattice[0][Y], g.global_lattice[0][Z], g.global_lattice[0][T], g.local_lattice[0][X], g.local_lattice[0][Y], g.local_lattice[0][Z], g.local_lattice[0][T], g.setup_iter[0], g.post_smooth_iter[0], g.block_iter[0] );
  #else
  sprintf(tmp_string, "<header>\n\tclifford basis: %s\n\tm0: %.14lf\n\tcsw: %.14lf\n\tmu: %.14lf\n\tclov plaq: %.14lf\n\thopp plaq: %.14lf\n\tcomputed plaq: %.14lf\n\tclov conf name: %s\n\thopp conf name: %s\n\tX: %d\n\tY: %d\n\tZ: %d\n\tT: %d\n\tX local: %d\n\tY local: %d\n\tZ local: %d\n\tT local: %d\n\tsetup iter: %d\n\tpost smoothing iter: %d\n\tblock iter: %d\n</header>", CLIFFORD_BASIS, g.solve_m0, g.csw, g.tm_mu, g.plaq_clov, g.plaq_hopp, g.plaq, g.in_clov, g.in, g.global_lattice[0][X], g.global_lattice[0][Y], g.global_lattice[0][Z], g.global_lattice[0][T], g.local_lattice[0][X], g.local_lattice[0][Y], g.local_lattice[0][Z], g.local_lattice[0][T], g.setup_iter[0], g.post_smooth_iter[0], g.block_iter[0] );
  #endif
  message_length=(n_uint64_t) strlen(tmp_string);
  limeheader = limeCreateHeader(1, 1, "dd_alpha_amg-header", message_length);
  
  ASSERT(limeheader != (LimeRecordHeader*)NULL);
  ASSERT(limeWriteRecordHeader(limeheader, limewriter)>=0);
  ASSERT(limeWriteRecordData(tmp_string, &message_length, limewriter)>=0);
  limeDestroyHeader(limeheader);
  
  int vol=1, mu;
  for ( mu=0; mu<4; mu++ ){
    vol*=g.global_lattice[0][mu];
  }
  message_length=data_size_per_site*vol;
  limeheader = limeCreateHeader(1, 1, "scidac-binary-data", message_length);
  ASSERT(limeWriteRecordHeader(limeheader, limewriter)>=0);
  limeDestroyHeader(limeheader);
  message_length=1;
  //make one fake record-write to set offset
  ASSERT(limeWriteRecordData(tmp_string, &message_length, limewriter)>=0);
  *offset = (n_uint64_t) ftell(*fout)-1;
  
  limeDestroyWriter(limewriter);
}

#endif


void lime_read_conf( double *input_data, char *input_name, double *conf_plaq ) {
  
  /*********************************************************************************
   * Reads in the configuration.
   * - double *input_data: Variable where conf data is stored.
   * - char *input_name: Name of the input file.
   * - double *conf_plaq: Holds the plaquette of given configuration.                            
   *********************************************************************************/
  #ifdef HAVE_LIME
  int t, z, y, x, i, k, desired_rank,
  *gl=g.global_lattice[0], *ll=g.local_lattice[0], read_size = 4*18*ll[X], precision=64;
  double *input_data_pt, plaq;
  float *float_buffer = NULL;
  FILE* fin = NULL;
  MPI_Request sreq;
  confbuffer_struct buffer[2]; // Having two buffers allows communication hiding.
  confbuffer_struct *buffer_pt = NULL;
  
  
  buffer[0].next = &(buffer[1]);
  buffer[1].next = &(buffer[0]);  
  buffer[0].data = NULL;
  buffer[1].data = NULL;
  
  if ( g.my_rank == 0 ) {
    MALLOC( buffer[0].data, double, read_size );
    MALLOC( buffer[1].data, double, read_size );
    buffer_pt = &buffer[0];
  }
  ASSERT( (fin = fopen( input_name, "rb" )) != NULL );
  
  if ( g.my_rank == 0 ) {
    lime_fileinfo lime;
    lime_read_info(&fin, &lime, "ildg-binary-data");
    lime_check_info(lime, 4*9); // n link * matrix size
    precision = lime.precision;
    if ( precision==32 ) {
      MALLOC( float_buffer, float, read_size );
    }
    if(g.print>0) printf0("Lattice size nx=%d ny=%d nz=%d nt=%d\n", lime.l[X], lime.l[Y], lime.l[Z], lime.l[T]);
    if(g.print>0) printf0("Reading configuration in %s precision.\n", precision==64?"double":"single");
    fseek(fin, lime.offset, SEEK_SET);
    plaq=lime.plaquette;
    *conf_plaq = plaq;
    if(g.print>0)  printf0("\nDesired average plaquette: %.13lf in [0,1]\n", plaq );
  }
  input_data_pt = input_data;
  
  // Distribute data to according processes
  if ( g.my_rank == 0 ) {
    if (precision==64) {
      ASSERT( fread( buffer_pt->data, sizeof(double), read_size, fin ) > 0 );
      for ( i=0; i<read_size; i++ )
        byteswap8( (char *) ( buffer_pt->data + i ) );
    }
    else if (precision==32) {
      ASSERT( fread( float_buffer, sizeof(float), read_size, fin ) > 0 );
      for ( int i=0; i<read_size; i++ ) {
        byteswap( (char *) (float_buffer+i) );
        buffer_pt->data[i] = float_buffer[i];
      }
    }
  }
  
  k = 0;
  for ( t=0; t<gl[T]; t++ )
    for ( z=0; z<gl[Z]; z++ )
      for ( y=0; y<gl[Y]; y++ )
        for ( x=0; x<gl[X]; x+=ll[X] ) {
          desired_rank = process_index( t, z, y, x, ll );
          
          if ( g.my_rank == 0 ) {
            MPI_Isend( buffer_pt->data, read_size, MPI_DOUBLE, desired_rank, k, g.comm_cart, &sreq );
            if ( ! ( t == gl[T]-1 && z == gl[Z]-1 && y == gl[Y]-1 && x == gl[X]-ll[X]  ) ) {
              if (precision==64) {
                ASSERT( fread( buffer_pt->next->data, sizeof(double), read_size, fin ) > 0 );
                for ( i=0; i<read_size; i++ )
                  byteswap8( (char *) ( buffer_pt->next->data + i ) );
              }
              else if (precision==32) {
                ASSERT( fread( float_buffer, sizeof(float), read_size, fin ) > 0 );
                for ( i=0; i<read_size; i++ ) {
                  byteswap( (char *) (float_buffer+i) );
                  buffer_pt->next->data[i] =  float_buffer[i];
                }
              }
            }
          }
          
          if ( g.my_rank == desired_rank ) {
            MPI_Recv( input_data_pt, read_size, MPI_DOUBLE, 0, k, g.comm_cart, MPI_STATUS_IGNORE );
            
            swap_spin_in_conf(input_data_pt, read_size);
            
            input_data_pt += read_size;
          }
          
          if ( g.my_rank == 0 ) {
            MPI_Wait( &sreq, MPI_STATUS_IGNORE );
            buffer_pt = buffer_pt->next;
          }
          
          k = (k+1)%10000;
        }
        fclose( fin );
        if ( g.my_rank == 0 ) {
          FREE( buffer[0].data, double, read_size );
          FREE( buffer[1].data, double, read_size );
          if (precision==32)
            FREE( float_buffer, float, read_size );
        }
        #else
        error0("Lime not enabled. Use -DHAVE_LIME during compilation. Aborting...");
        #endif
        
}

void lime_read_vector( double *phi, char *filename ) {
  
  #ifdef HAVE_LIME
  
  int t, z, y, x, *gl=g.global_lattice[0], *ll=g.local_lattice[0], bar_size = 24*ll[X], desired_rank, precision=64;
  double *phi_pt = phi, t0, t1;
  float *float_buffer = NULL;
  FILE* file = NULL;
  MPI_Request sreq;
  confbuffer_struct buffer[2];
  confbuffer_struct *buffer_pt = NULL;
  
  buffer[0].next = &(buffer[1]);
  buffer[1].next = &(buffer[0]);  
  buffer[0].data = NULL;
  buffer[1].data = NULL;
  
  int i;
  
  if ( g.my_rank == 0 ) {
    MALLOC( buffer[0].data, double, bar_size );
    MALLOC( buffer[1].data, double, bar_size );
    buffer_pt = &buffer[0];
  }
  
  t0 = MPI_Wtime();
  
  if ( g.my_rank == 0 ) {
    ASSERT( (file = fopen( filename, "rb" )) != NULL );
    lime_fileinfo lime;
    lime_read_info(&file, &lime, "scidac-binary-data");
    precision=lime.precision;
    printf0("Reading vector in %s precision.\n", precision==64?"double":"single");
    lime_check_info(lime, lime.spin*lime.colour);
    if ( precision==32 )
      MALLOC( float_buffer, float, bar_size );
    fseek(file, lime.offset, SEEK_SET);
  }
  
  printf0("reading from file \"%s\" ...\n", filename );
  
  if ( g.my_rank == 0 ) {
    if (precision==64) {
      ASSERT( fread( buffer_pt->data, sizeof(double), bar_size, file ) > 0 );
      for ( i=0; i<bar_size; i++ )
        byteswap8( (char *) ( buffer_pt->data + i ) );
    }
    else if (precision==32) {
      ASSERT( fread( float_buffer, sizeof(float), bar_size, file ) > 0 );
      for ( int i=0; i<bar_size; i++ ) {
        byteswap( (char *) (float_buffer+i) );
        buffer_pt->data[i] = float_buffer[i];
      }
    }
  }
  
  for ( t=0; t<gl[T]; t++ )
    for ( z=0; z<gl[Z]; z++ )
      for ( y=0; y<gl[Y]; y++ )
        for ( x=0; x<gl[X]; x+=ll[X] ) {
          desired_rank = process_index( t, z, y, x, ll );
          
          if ( g.my_rank == 0 ) {
            MPI_Isend( buffer_pt->data, bar_size, MPI_DOUBLE, desired_rank, 0, g.comm_cart, &sreq );
            if ( ! ( t == gl[T]-1 && z == gl[Z]-1 && y == gl[Y]-1 && x == gl[X]-ll[X]  ) ) {
              if (precision==64) {
                ASSERT( fread( buffer_pt->next->data, sizeof(double), bar_size, file ) > 0 );
                for ( i=0; i<bar_size; i++ )
                  byteswap8( (char *) ( buffer_pt->next->data + i ) );
              }
              else if (precision==32) {
                ASSERT( fread( float_buffer, sizeof(float), bar_size, file ) > 0 );
                for ( int i=0; i<bar_size; i++ ) {
                  byteswap( (char *) (float_buffer+i) );
                  buffer_pt->next->data[i] = float_buffer[i];
                }
              }        
            }
          }
          
          if ( g.my_rank == desired_rank ) {
            MPI_Recv( phi_pt, bar_size, MPI_DOUBLE, 0, 0, g.comm_cart, MPI_STATUS_IGNORE );
            phi_pt += bar_size;
          }
          
          if ( g.my_rank == 0 ) {
            MPI_Wait( &sreq, MPI_STATUS_IGNORE );
            buffer_pt = buffer_pt->next;
          }
        }
        
        if ( g.my_rank == 0 ){
          fclose( file );
        }
        
        t1 = MPI_Wtime();
        
        if ( g.my_rank == 0 ) {
          FREE( buffer[0].data, double, bar_size );
          FREE( buffer[1].data, double, bar_size );
          if (precision==32)
            FREE( float_buffer, float, bar_size );
        }
        
        printf0("...done (%lf seconds)\n\n", t1-t0 ); 
        #else
        error0("Lime not enabled. Use -DHAVE_LIME during compilation. Aborting...");
        #endif
}

void lime_write_vector( double *phi, char *filename ) {
  
  #ifdef HAVE_LIME
  
  int t, z, y, x, *gl=g.global_lattice[0], *ll=g.local_lattice[0], bar_size = 24*ll[X], desired_rank, precision=64;
  double *phi_pt = phi, t0, t1;
  float *float_buffer;
  FILE* file = NULL;
  MPI_Request sreq, rreq;
  confbuffer_struct buffer[2];
  confbuffer_struct *buffer_pt = NULL;
  
  buffer[0].next = &(buffer[1]);
  buffer[1].next = &(buffer[0]);  
  buffer[0].data = NULL;
  buffer[1].data = NULL;
  
  int i;
  
  if ( g.my_rank == 0 ) {
    MALLOC( buffer[0].data, double, bar_size );
    MALLOC( buffer[1].data, double, bar_size );
    buffer_pt = &buffer[0];
  }
  
  t0 = MPI_Wtime();
  if ( g.my_rank == 0 ) {
    printf0("writing file \"%s\" ...\n", filename );
    ASSERT( (file = fopen( filename, "wb" )) != NULL );
    n_uint64_t offset;
    lime_write_info(&file, "scidac-binary-data", (n_uint64_t) 4*3*sizeof(complex_double), &offset);
    fseek(file, offset, SEEK_SET);
  }
  
  for ( t=0; t<gl[T]; t++ ) {
    for ( z=0; z<gl[Z]; z++ )
      for ( y=0; y<gl[Y]; y++ )
        for ( x=0; x<gl[X]; x+=ll[X] ) {
          desired_rank = process_index( t, z, y, x, ll );
          
          if ( g.my_rank == 0 ) {
            MPI_Irecv( buffer_pt->next->data, bar_size, MPI_DOUBLE, desired_rank, 0, g.comm_cart, &rreq );
            if ( ! ( t == 0 && z == 0 && y == 0 && x == 0  ) ) {
              for ( i=0; i<bar_size; i++ ) {
                byteswap8( (char *) ( buffer_pt->data + i ) );
              }
              fwrite( buffer_pt->data, sizeof(double), bar_size, file );
            }
            buffer_pt = buffer_pt->next;
          }
          if ( g.my_rank == desired_rank ) {
            MPI_Isend( phi_pt, bar_size, MPI_DOUBLE, 0, 0, g.comm_cart, &sreq);
            phi_pt += bar_size;
            MPI_Wait( &sreq, MPI_STATUS_IGNORE );
          }
          
          if ( g.my_rank == 0 ) {
            MPI_Wait( &rreq, MPI_STATUS_IGNORE );
          }
        }
  }
  
  if ( g.my_rank == 0 ) {	
    for ( i=0; i<bar_size; i++ ) {
      byteswap8( (char *) ( buffer_pt->data + i ) );
    }
    fwrite( buffer_pt->data, sizeof(double), bar_size, file );
  }
  
  if ( g.my_rank == 0 ){
    fclose( file );
  }
  
  t1 = MPI_Wtime();
  
  if ( g.my_rank == 0 ) {
    FREE( buffer[0].data, double, bar_size );
    FREE( buffer[1].data, double, bar_size );
    if (precision==32)
      FREE( float_buffer, float, bar_size );
  }
  
  printf0("...done (%lf seconds)\n\n", t1-t0 ); 
  #else
  warning0("Lime not enabled. Use -DHAVE_LIME during compilation. Vector not saved...\n");
  #endif
}
