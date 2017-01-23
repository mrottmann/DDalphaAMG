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
 * 2016 Aug. 25: Issaku Kanamori
 *   modified read_conf(), to avoid reading from non-master nodes
 *
 * 2016 Oct. 3:  Issaku Kanamori
 *   added read_conf_multi()
 */

#include "main.h"

#ifdef HAVE_HDF5

hid_t openCreateGroup( hid_t parent, char* groupname ) {
  htri_t group_ex = H5Lexists(parent, groupname, H5P_DEFAULT);
  ASSERT( group_ex >= 0 );
  hid_t group_id, status;
  if(group_ex == 0)
  {
    printf0("group \"%s\" does not exist.\n", groupname);
    hid_t gcpl_id = H5Pcreate(H5P_GROUP_CREATE);
    group_id = H5Gcreate(parent, groupname, H5P_DEFAULT, gcpl_id, H5P_DEFAULT);
    status = H5Pclose(gcpl_id);
    ASSERT( status == 0 );
  }
  else
  {
    group_id = H5Gopen(parent, groupname, H5P_DEFAULT);
  }
  ASSERT( group_id >= 0 );
  return group_id;
}


unsigned int stepIntoEigenmode( int index ) {
  ASSERT( h5info.isOpen == 1 );
  ASSERT( h5info.configgroup_id > 0 );
  double fnEntryTime;
  char eigenmodename[STRINGLENGTH];
  fnEntryTime = MPI_Wtime();
  sprintf(eigenmodename, "eigenmode%i", index);
  if(h5info.thiseigenmodegroup_id > 0)
    H5Gclose(h5info.thiseigenmodegroup_id);
  hid_t thiseigenmodegroup_id = openCreateGroup( h5info.configgroup_id, eigenmodename );
  if(thiseigenmodegroup_id > 0)
    h5info.thiseigenmodegroup_id = thiseigenmodegroup_id;
  else
  {
    printf0("io.c: stepIntoEigenmode: Could not open eigenmode group for eigenmode %i.", index);
    return -1;
  }
  h5info.ioTime += MPI_Wtime()-fnEntryTime;
  return thiseigenmodegroup_id;
}


void writeSmallDataset_double( hid_t parent, char* name, double value) {
  hid_t data_id, space_id, status;
    /*
     * data set already exists. overwrite.
     */
  htri_t link_ex = H5Lexists(parent, name, H5P_DEFAULT);
  ASSERT( link_ex >= 0 );
  if(link_ex == 0)
  {
    space_id = H5Screate(H5S_SCALAR);
    ASSERT ( space_id >= 0 );
    data_id  = H5Dcreate(parent, name, H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Sclose(space_id);
    ASSERT ( status == 0 );
  }
  else 
  {
    data_id = H5Dopen(parent, name, H5P_DEFAULT);
  }
  ASSERT( data_id > 0 );

  status = H5Dwrite(data_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (char*)&value);
  ASSERT ( status == 0 );

  status = H5Dclose(data_id);
  ASSERT ( status == 0 );
}


void writeSmallDataset_doublearr( hid_t parent, char* name, double *value, int len) {
  hid_t data_id, space_id, status;
  htri_t link_ex = H5Lexists(parent, name, H5P_DEFAULT);
  ASSERT( link_ex >= 0 );
  if(link_ex == 0)
  {
    space_id = H5Screate(H5S_SIMPLE);
    hsize_t size[1];
    size[0] = len;
    H5Sset_extent_simple(space_id, 1, size, size);
    ASSERT ( space_id >= 0 );
    data_id  = H5Dcreate(parent, name, H5T_NATIVE_DOUBLE, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space_id);
  }
  else 
  {
    /*
     * data set already exists. overwrite.
     */
    data_id = H5Dopen(parent, name, H5P_DEFAULT);
  }
  ASSERT( data_id > 0 );

  status = H5Dwrite(data_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (char*)value);
  ASSERT(status == 0);

  status = H5Dclose(data_id);
  ASSERT(status == 0);
}


void writeSmallDataset_int( hid_t parent, char* name, int value) {
  hid_t data_id, space_id, status;
    /*
     * data set already exists. overwrite.
     */
  htri_t link_ex = H5Lexists(parent, name, H5P_DEFAULT);
  ASSERT( link_ex >= 0 );
  if(link_ex == 0)
  {
    space_id = H5Screate(H5S_SCALAR);
    ASSERT ( space_id >= 0 );
    data_id  = H5Dcreate(parent, name, H5T_NATIVE_INT, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space_id);
  }
  else 
  {
    data_id = H5Dopen(parent, name, H5P_DEFAULT);
  }
  ASSERT( data_id > 0 );

  status = H5Dwrite(data_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (char*)&value);
  ASSERT ( status == 0 );

  status = H5Dclose(data_id);
  ASSERT ( status == 0 );
}


void writeSmallDataset_intarr( hid_t parent, char* name, int *value, int len) {
  hid_t data_id, space_id, status;
  htri_t link_ex = H5Lexists(parent, name, H5P_DEFAULT);
  ASSERT( link_ex >= 0 );
  if(link_ex == 0)
  {
    space_id = H5Screate(H5S_SIMPLE);
    hsize_t size[1];
    size[0] = len;
    H5Sset_extent_simple(space_id, 1, size, size);
    ASSERT ( space_id >= 0 );
    data_id  = H5Dcreate(parent, name, H5T_NATIVE_INT, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space_id);
  }
  else 
  {
    /*
     * data set already exists. overwrite.
     */
    data_id = H5Dopen(parent, name, H5P_DEFAULT);
  }
  ASSERT( data_id > 0 );

  status = H5Dwrite(data_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (char*)value);
  ASSERT(status == 0);

  status = H5Dclose(data_id);
  ASSERT(status == 0);
}


void writeSmallDataset_str( hid_t parent, char* name, char* value, unsigned int len) {
  hid_t data_id, space_id, status;
    /*
     * data set already exists. overwrite.
     */
  htri_t link_ex = H5Lexists(parent, name, H5P_DEFAULT);
  ASSERT( link_ex >= 0 );
  if(link_ex == 0)
  {
    space_id = H5Screate(H5S_SIMPLE);
    hsize_t size[1];
    size[0] = len;
    H5Sset_extent_simple(space_id, 1, size, size);
    ASSERT ( space_id >= 0 );
    data_id  = H5Dcreate(parent, name, H5T_C_S1, space_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(space_id);
  }
  else 
  {
    data_id = H5Dopen(parent, name, H5P_DEFAULT);
  }
  ASSERT( data_id > 0 );

  status = H5Dwrite(data_id, H5T_C_S1, H5S_ALL, H5S_ALL, H5P_DEFAULT, (char*)value);
  ASSERT(status == 0);

  status = H5Dclose(data_id);
  ASSERT(status == 0);
}


unsigned int initFile( char *filename, const int mode, level_struct *l ) {
  /*
   * This routine opens a hdf5 file and returns the file_id.
   * It checks if the file exists and if it does, it checks that it is an hdf5 file.
   * Parts of this code is taken from qdp++: lib/qdp_hdf5.cc
   */
  hid_t fapl_id = -1, fcpl_id = -1, file_id = -1;
  MPI_Info info = MPI_INFO_NULL;
  FILE* file = NULL;
  /* stripesize of 4M, 4 blocks aligned. 
   * This seems to be a good default for most lustre installations and this filesize
   */
  const unsigned int stripesize = 4*1024*1024;
  const unsigned int align = 4;
  double fnEntryTime;

  //file should not be already open:
  ASSERT(h5info.isOpen == 0);

  fnEntryTime = MPI_Wtime();

  h5info.mode = mode;
  h5info.filename = filename;
  /*
   * the following optimizations are reported to make things faster (file open / close)
   */
  MPI_Info_create(&info);
  MPI_Info_set(info, "romio_ds_read", "disable");
  MPI_Info_set(info, "romio_ds_write", "disable");
  MPI_Info_set(info, "romio_cb_read", "enable");
  MPI_Info_set(info, "romio_cb_write", "enable");
  /*
   * set up file access property lists, one for access, another one for creation:
   */
  fapl_id = H5Pcreate(H5P_FILE_ACCESS);
  fcpl_id = H5Pcreate(H5P_FILE_CREATE);
  if(stripesize > 0)
  {
    H5Pset_alignment(fapl_id, align, stripesize);
    int btree_ik = ceil((stripesize-4096) / 96);
    H5Pset_istore_k(fcpl_id, btree_ik);
  }
  H5Pset_fapl_mpio(fapl_id, g.comm_cart, info);

  /* 
   * check if the file already exists:
   */
  /*
   * exists
   * exists = 0: the file exists but is not a hdf5 file
   * exists > 0: the file exists and is a proper hdf5 file
   * exists < 0: the file does not exists or other error.
   */
  char exists;
  if(g.my_rank == 0)
  {
      file = fopen( filename, "rb" );
      if(file == NULL)
        exists = -1;
      else
      {
        fclose(file);
        // the file exists
        exists = H5Fis_hdf5(filename);
      }
      ASSERT( exists != 0 );
  }
  MPI_Bcast(&exists, 1, MPI_CHAR, 0, g.comm_cart);
  if(exists < 0 && mode == _READ)
  {
    printf0("io.c: initFile: file \"%s\" does not exist.", filename);
    return -1;
  }
  else if( exists < 0 )
  {
    //the file does not exist and needs to be created.
    file_id = H5Fcreate(filename, H5F_ACC_EXCL, fcpl_id, fapl_id);
  }
  else if(exists > 0)
  {
    /*
     * open hdf5 file:
     */
    if(mode == _READ)
      file_id = H5Fopen(filename, H5F_ACC_RDONLY, fapl_id);
    else
      file_id = H5Fopen(filename, H5F_ACC_RDWR, fapl_id);
  }
  if( file_id < 0 )
  {
    printf0("io.c: initFile: file \"%s\" could not be opened for unknown reason.", filename);
    return -1;
  }
  /*
   * close the file access / creation property lists
   */
  H5Pclose(fapl_id);
  H5Pclose(fcpl_id);

  h5info.file_id = file_id;
  /* 
   * writes the header information as attributes / datasets into the hdf5 file
   */
  char strbuffer[256];
  /*
   * create the groups:
   */
  hid_t rootgroup_id = openCreateGroup(file_id, ROOTGROUPNAME);
  if( rootgroup_id < 0 ) {
    printf0("io.c: initFile: rootgroup \"%s\" could not be opened.", ROOTGROUPNAME);
    return -2;
  }
  hid_t eigenmodegroup_id = openCreateGroup(rootgroup_id, EIGENMODEGROUPNAME);
  if( eigenmodegroup_id < 0 ) {
    printf0("io.c: initFile: eigenmodegroup \"%s\" could not be opened.", EIGENMODEGROUPNAME);
    return -2;
  }
  /*
   * build the config name (will be used as group name)
   * configname = g.in + g.in_clov (after removal of the slashes)
   */
  char *pchr = strrchr(g.in, '/');
  if(pchr != NULL)
    strcpy(strbuffer, pchr+1);
  else
    strcpy(strbuffer, g.in);
  pchr = strrchr(g.in_clov, '/');
  if(pchr != NULL)
    strcat(strbuffer, pchr+1);
  else
    strcat(strbuffer, g.in_clov);
  hid_t configgroup_id = openCreateGroup(eigenmodegroup_id, strbuffer);
  if( configgroup_id < 0 ) {
    printf0("io.c: initFile: configgroup \"%s\" could not be opened.", strbuffer);
    return -2;
  }

  h5info.rootgroup_id = rootgroup_id;
  h5info.configgroup_id = configgroup_id;
  h5info.eigenmodegroup_id = eigenmodegroup_id;
  //thiseigenmodegroup is not open:
  h5info.thiseigenmodegroup_id = -1;
  if(mode == _WRITE)
  {
    /*
     * write datasets:
     */
    //reverse the ordering of the lattice sizes for convenience (tzyx -> xyzt)
    int ls[4], lls[4];
    ls[0] = l->global_lattice[3];
    ls[1] = l->global_lattice[1];
    ls[2] = l->global_lattice[2];
    ls[3] = l->global_lattice[0];
    lls[0] = l->local_lattice[3];
    lls[1] = l->local_lattice[1];
    lls[2] = l->local_lattice[2];
    lls[3] = l->local_lattice[0];
    writeSmallDataset_double(configgroup_id, "m0", l->real_shift);
    writeSmallDataset_double(configgroup_id, "csw", g.csw);
    writeSmallDataset_double(configgroup_id, "plaquette_clov", g.plaq_clov);
    writeSmallDataset_double(configgroup_id, "plaquette_hopp", g.plaq_hopp);
    writeSmallDataset_int(configgroup_id, "iter_block", l->block_iter);
    writeSmallDataset_int(configgroup_id, "iter_setup", l->setup_iter);
    writeSmallDataset_int(configgroup_id, "iter_post_smoothing", l->post_smooth_iter);
    writeSmallDataset_int(configgroup_id, "number_vectors", l->num_eig_vect);
    writeSmallDataset_intarr(configgroup_id, "lattice_size", ls, 4);
    writeSmallDataset_intarr(configgroup_id, "local_lattice_size", lls, 4);
    writeSmallDataset_str(configgroup_id, "clifford_basis", CLIFFORD_BASIS, 256);
    writeSmallDataset_str(configgroup_id, "config_name_clover", g.in_clov, 256);
    writeSmallDataset_str(configgroup_id, "config_name_hopp", g.in, 256);
  }

  h5info.isOpen = 1;
  h5info.ioTime = MPI_Wtime() - fnEntryTime;
  return h5info.file_id;
}


void closeFile() {
  double fnEntryTime;
  fnEntryTime = MPI_Wtime();
  if(! h5info.isOpen)
    return;
  if(h5info.thiseigenmodegroup_id > 0)
    H5Gclose(h5info.thiseigenmodegroup_id);
  H5Gclose(h5info.configgroup_id);
  H5Gclose(h5info.eigenmodegroup_id);
  H5Gclose(h5info.rootgroup_id);
  H5Fclose(h5info.file_id);
  h5info.isOpen = 0;
  h5info.ioTime += fnEntryTime - MPI_Wtime();
  printf0("io.c: closed file \"%s\". Time spent in io routines: %lf seconds wct (%lf node hours)\n", h5info.filename, h5info.ioTime, h5info.ioTime/3600*g.num_processes);
  h5info.ioTime = 0.0;
}
#endif


void byteswap( char *in ) {

  char tmp[4];
  tmp[0] = in[3];
  tmp[1] = in[2];
  tmp[2] = in[1];
  tmp[3] = in[0];

  in[0] = tmp[0];
  in[1] = tmp[1];
  in[2] = tmp[2];
  in[3] = tmp[3];
}
void byteswap8(char *in ) {
  
  char tmp[8];
  tmp[0] = in[7];
  tmp[1] = in[6];
  tmp[2] = in[5];
  tmp[3] = in[4];
  tmp[4] = in[3];
  tmp[5] = in[2];
  tmp[6] = in[1];
  tmp[7] = in[0];

  in[0] = tmp[0];
  in[1] = tmp[1];
  in[2] = tmp[2];
  in[3] = tmp[3];
  in[4] = tmp[4];
  in[5] = tmp[5];
  in[6] = tmp[6];
  in[7] = tmp[7];
}

void read_conf( double *input_data, char *input_name, double *conf_plaq, level_struct *l ) {
  
/*********************************************************************************
* Reads in the configuration.
* - double *input_data: Variable where conf data is stored.
* - char *input_name: Name of the input file.
* - double *conf_plaq: Holds the plaquette of given configuration.                                                
*********************************************************************************/

  int t, z, y, x, i, j, k, mu, lsize[4], desired_rank,
      *gl=l->global_lattice, *ll=l->local_lattice, read_size = 4*18*ll[X];
  double *input_data_pt, plaq;
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
  
    ASSERT( (fin = fopen( input_name, "rb" )) != NULL );
  
    // Read in lattice size in each dimension
    ASSERT( fread( lsize, sizeof(int), 4, fin ) > 0 );
#ifdef BIG_ENDIAN_CNFG
    byteswap( (char *) &lsize[0] );  
    byteswap( (char *) &lsize[1] );
    byteswap( (char *) &lsize[2] );
    byteswap( (char *) &lsize[3] );
#endif
    for ( mu=0; mu<4; mu++ )
      ASSERT( lsize[mu] == gl[mu] );
  
  // Read in plaquette
    ASSERT( fread( &plaq, sizeof(double), 1, fin ) > 0 );
#ifdef BIG_ENDIAN_CNFG
    byteswap8( (char *) &plaq );
#endif
    printf0("\nDesired average plaquette: %.13lf in [0,3]\n", plaq );
    printf0("\nDesired average plaquette: %.13lf in [0,1]\n", plaq/3.0 );
    *conf_plaq = plaq; 
  }
  MPI_Bcast(conf_plaq, 1, MPI_DOUBLE, 0, g.comm_cart);

  input_data_pt = input_data;
  
  // Distribute data to according processes
  if ( g.my_rank == 0 ) 
    ASSERT( fread( buffer_pt->data, sizeof(double), read_size, fin ) > 0 );
  
  k = 0;
  for ( t=0; t<gl[T]; t++ )
    for ( z=0; z<gl[Z]; z++ )
      for ( y=0; y<gl[Y]; y++ )
        for ( x=0; x<gl[X]; x+=ll[X] ) {
          desired_rank = process_index( t, z, y, x, ll );
          
          if ( g.my_rank == 0 ) {
            MPI_Isend( buffer_pt->data, read_size, MPI_DOUBLE, desired_rank, k, g.comm_cart, &sreq );
            if ( ! ( t == gl[T]-1 && z == gl[Z]-1 && y == gl[Y]-1 && x == gl[X]-ll[X]  ) )
              ASSERT( fread( buffer_pt->next->data, sizeof(double), read_size, fin ) > 0 );
          }
          
          if ( g.my_rank == desired_rank ) {
            MPI_Recv( input_data_pt, read_size, MPI_DOUBLE, 0, k, g.comm_cart, MPI_STATUS_IGNORE );
#ifdef BIG_ENDIAN_CNFG
            for ( i=0; i<read_size; i++ ) {
              byteswap8( (char *) ( input_data_pt + i ) );
            }
#endif
            if ( t == gl[T]-1 ) {
              if ( g.anti_pbc ) {
                for ( j=0; j<read_size; j+=4*18 )
                  for ( i=0; i<18; i++ )
                    (input_data_pt+j+T*18)[i] = -(input_data_pt+j+T*18)[i];
              }
            }
            input_data_pt += read_size;
          }
          
          if ( g.my_rank == 0 ) {
            MPI_Wait( &sreq, MPI_STATUS_IGNORE );
            buffer_pt = buffer_pt->next;
          }
          
          k = (k+1)%10000;
        }
  
  if( g.my_rank == 0)
    fclose( fin );
  
  if ( g.my_rank == 0 ) {
    FREE( buffer[0].data, double, read_size );
    FREE( buffer[1].data, double, read_size );
  }

}



void read_conf_multi( double *input_data, char *input_name, double *conf_plaq, level_struct *l ) {
  
/*********************************************************************************
* Reads in the configuration from multi files.
* - double *input_data: Variable where conf data is stored.
* - char *input_name: Name of the input file.
* - double *conf_plaq: Holds the plaquette of given configuration.                                                
*********************************************************************************/

  int t, z, y, i, j,  mu, lsize[4],
      *gl=l->global_lattice, *ll=l->local_lattice, read_size = 4*18*ll[X];
  double *input_data_pt, plaq;
  FILE* fin = NULL;

  char postfix[100];
  char full_input_name[1000];

  /****************************************************************
   * 2016 Oct. 3:  Issaku Kanamori
   *   it seems MPI rank = g.my_coords[T] * Pz * Py * Px
   *                     g.my_coords[Z] * Py * Px
   *                   g.my_coords[Y] * Px
   *                    g.my_coords[X]
   *   for (Pt,Pz,Py,Px) processor lattice.
   *  ( It might depend on the implementation of MPI.  At least,
   *    I could not find any specification in the manural of MPI_Cart_* )
   *
   *  If you have a problem in opening the file due to misdistribution
   *  of configuration files to computing nodes, try 
   *  #define READ_CONF_MULTI_CHECKFILE
   *  and find out which node is trying to read which file.
   ******************************************************************/
  printf0("reading from multi file configuration\n");
  sprintf(postfix, ".pt%dpz%dpy%dpx%d", 
          g.my_coords[T], g.my_coords[Z], g.my_coords[Y], g.my_coords[X]);
  sprintf(full_input_name,"%s%s", input_name, postfix);
#ifdef READ_CONF_MULTI_CHECKFILE
  fin = fopen( full_input_name, "rb" );
  int open_flag=0;
  int open_flag_sum=0;
  if(fin == NULL){
    printf("ERROR!  node=%d,  filename=%s\n", g.my_rank, full_input_name);
    fflush(0);
    open_flag=1;
  } else {
    printf("opened: node=%d,  filename=%s\n", g.my_rank, full_input_name);
    fflush(0);
  }
  MPI_Allreduce( &open_flag, &open_flag_sum, 1, MPI_INT, MPI_SUM, g.comm_cart );
  ASSERT(open_flag_sum==0);
#else
  ASSERT( (fin = fopen( full_input_name, "rb" )) != NULL );
#endif

  // Read in lattice size in each dimension
  ASSERT( fread( lsize, sizeof(int), 4, fin ) > 0 );
#ifdef BIG_ENDIAN_CNFG
  byteswap( (char *) &lsize[0] );  
  byteswap( (char *) &lsize[1] );
  byteswap( (char *) &lsize[2] );
  byteswap( (char *) &lsize[3] );
#endif

  for ( mu=0; mu<4; mu++ )
    ASSERT( lsize[mu] == gl[mu] );

  // Read in plaquette
  ASSERT( fread( &plaq, sizeof(double), 1, fin ) > 0 );
#ifdef BIG_ENDIAN_CNFG
  byteswap8( (char *) &plaq );
#endif
  printf0("\nDesired average plaquette: %.13lf in [0,3]\n", plaq );
  printf0("\nDesired average plaquette: %.13lf in [0,1]\n", plaq/3.0 );
  *conf_plaq = plaq; 

  // check if value of the plaauette is the same 
  MPI_Bcast(conf_plaq, 1, MPI_DOUBLE, 0, g.comm_cart);
  ASSERT( plaq == *conf_plaq );

  // read the configuration
  input_data_pt = input_data;
  for ( t=0; t<ll[T]; t++ )
    for ( z=0; z<ll[Z]; z++ )
      for ( y=0; y<ll[Y]; y++ ) {
        // read ll[X] data at once  (see def. of read_size)
        ASSERT( fread( input_data_pt, sizeof(double), read_size, fin ) > 0 );
        #ifdef BIG_ENDIAN_CNFG
        for ( i=0; i<read_size; i++ ) {
          byteswap8( (char *) ( input_data_pt + i ) );
        }
        #endif
        if ( (g.my_coords[T]+1 == gl[T]/ll[T]) && (t == ll[T]-1) ) {
          if ( g.anti_pbc ) {
            for ( j=0; j<read_size; j+=4*18 )
              for ( i=0; i<18; i++ )
                (input_data_pt+j+T*18)[i] = -(input_data_pt+j+T*18)[i];
          }
        }
        input_data_pt += read_size;
      }

  fclose( fin );
}


void write_header( FILE **file, double *lambda, char* vector_type, int n, level_struct *l ) {
  
  fprintf( *file, "<header>\n" );
  fprintf( *file, "%s\n", vector_type );
  fprintf( *file, "clifford basis: %s\n", CLIFFORD_BASIS );
  fprintf( *file, "m0: %.14lf\n", l->real_shift );
  fprintf( *file, "csw: %.14lf\n", g.csw );
  fprintf( *file, "clov plaq: %.14lf\n", g.plaq_clov );
  fprintf( *file, "hopp plaq: %.14lf\n", g.plaq_hopp );
  fprintf( *file, "clov conf name: %s\n", g.in_clov );
  fprintf( *file, "hopp conf name: %s\n", g.in );
  fprintf( *file, "X: %d\n", l->global_lattice[X] );
  fprintf( *file, "Y: %d\n", l->global_lattice[Y] );
  fprintf( *file, "Z: %d\n", l->global_lattice[Z] );
  fprintf( *file, "T: %d\n", l->global_lattice[T] );
  fprintf( *file, "X local: %d\n", l->local_lattice[X] );
  fprintf( *file, "Y local: %d\n", l->local_lattice[Y] );
  fprintf( *file, "Z local: %d\n", l->local_lattice[Z] );
  fprintf( *file, "T local: %d\n", l->local_lattice[T] );
  fprintf( *file, "number of vectors: %d\n", n );
  fprintf( *file, "krylov subspace size: %d\n", 100 );
  fprintf( *file, "clifford basis: %s\n", CLIFFORD_BASIS );
  if ( lambda != NULL ) {
    fprintf( *file, "eigenvalues: " );
    for ( int i=0; i<2*n; i++ ) {
      fprintf( *file, "%.16lf ", lambda[i] );
    }
    fprintf( *file, "\n");
  }
  fprintf( *file, "</header>\n" );
}

#ifndef HAVE_HDF5
void vector_io( double *phi, char *filename, const int mode, level_struct *l ) {
  
  int t, z, y, x, *gl=l->global_lattice, *ll=l->local_lattice, bar_size = 24*ll[X], desired_rank;
  double *phi_pt = phi, t0, t1, norm;
  FILE* file = NULL;
  MPI_Request sreq, rreq;
  confbuffer_struct buffer[2];
  confbuffer_struct *buffer_pt = NULL;
  
  buffer[0].next = &(buffer[1]);
  buffer[1].next = &(buffer[0]);  
  buffer[0].data = NULL;
  buffer[1].data = NULL;
  
#ifdef BIG_ENDIAN_TV
  int i;
#endif
  
  if ( g.my_rank == 0 ) {
    MALLOC( buffer[0].data, double, bar_size );
    MALLOC( buffer[1].data, double, bar_size );
    buffer_pt = &buffer[0];
  }
  
  t0 = MPI_Wtime();
  
  if ( mode == _READ ) {
    
    if ( g.my_rank == 0 ) {
      ASSERT( (file = fopen( filename, "rb" )) != NULL );
      char *cur_line = NULL;
      MALLOC( cur_line, char, STRINGLENGTH );
      fgets( cur_line, STRINGLENGTH-1, file );
      if ( strcmp( cur_line, "<header>\n" ) ) {
        fseek( file, 0L, SEEK_SET );
      } else {
        do {
          ASSERT( fgets( cur_line, STRINGLENGTH-1, file ) );
        } while ( strcmp( cur_line, "</header>\n" ) );
      }
      FREE( cur_line, char, STRINGLENGTH );
    }
    
    printf0("reading from file \"%s\" ...\n", filename );
  
    if ( g.my_rank == 0 ) {
      ASSERT( fread( buffer_pt->data, sizeof(double), bar_size, file ) > 0 );
    }
    
    for ( t=0; t<gl[T]; t++ )
      for ( z=0; z<gl[Z]; z++ )
        for ( y=0; y<gl[Y]; y++ )
          for ( x=0; x<gl[X]; x+=ll[X] ) {
            desired_rank = process_index( t, z, y, x, ll );
            
            if ( g.my_rank == 0 ) {
              MPI_Isend( buffer_pt->data, bar_size, MPI_DOUBLE, desired_rank, 0, g.comm_cart, &sreq );
              if ( ! ( t == gl[T]-1 && z == gl[Z]-1 && y == gl[Y]-1 && x == gl[X]-ll[X]  ) )
                ASSERT( fread( buffer_pt->next->data, sizeof(double), bar_size, file ) > 0 );
            }
            
            if ( g.my_rank == desired_rank ) {
              MPI_Recv( phi_pt, bar_size, MPI_DOUBLE, 0, 0, g.comm_cart, MPI_STATUS_IGNORE );
#ifdef BIG_ENDIAN_TV
              for ( i=0; i<bar_size; i++ ) {
                byteswap8( (char *) ( phi_pt + i ) );
              }
#endif
              phi_pt += bar_size;
            }
            
            if ( g.my_rank == 0 ) {
              MPI_Wait( &sreq, MPI_STATUS_IGNORE );
              buffer_pt = buffer_pt->next;
            }
          }
          
  } else if ( mode == _WRITE ) {
    if ( g.my_rank == 0 ) {
      ASSERT( (file = fopen( filename, "wb" )) != NULL );
      write_header( &file, NULL, filename, 1, l );
    }
    printf0("writing file \"%s\" ...\n", filename );
    
    for ( t=0; t<gl[T]; t++ ) {
      for ( z=0; z<gl[Z]; z++ )
        for ( y=0; y<gl[Y]; y++ )
          for ( x=0; x<gl[X]; x+=ll[X] ) {
            desired_rank = process_index( t, z, y, x, ll );
            
            if ( g.my_rank == desired_rank ) {
              MPI_Isend( phi_pt, bar_size, MPI_DOUBLE, 0, 0, g.comm_cart, &sreq );
              phi_pt += bar_size;
            }
            
            if ( g.my_rank == 0 ) {
              MPI_Irecv( buffer_pt->next->data, bar_size, MPI_DOUBLE, desired_rank, 0, g.comm_cart, &rreq );
              
              if ( ! ( t == 0 && z == 0 && y == 0 && x == 0  ) ) {
#ifdef BIG_ENDIAN_TV
                for ( i=0; i<bar_size; i++ ) {
                  byteswap8( (char *) ( buffer_pt->data + i ) );
                }
#endif
                fwrite( buffer_pt->data, sizeof(double), bar_size, file );
              }
              
              MPI_Wait( &rreq, MPI_STATUS_IGNORE );
              buffer_pt = buffer_pt->next;
            }
            
            if ( g.my_rank == desired_rank ) {
              MPI_Wait( &sreq, MPI_STATUS_IGNORE );
            }
          }
    }
          
    if ( g.my_rank == 0 ) {
#ifdef BIG_ENDIAN_TV
      for ( i=0; i<bar_size; i++ ) {
        byteswap8( (char *) ( buffer_pt->data + i ) );
      }
#endif
      fwrite( buffer_pt->data, sizeof(double), bar_size, file );
    }
  
  } else
    ASSERT( mode == _READ || mode == _WRITE );
  
  if ( g.my_rank == 0 ){
    fclose( file );
  }

  t1 = MPI_Wtime();
  
  if ( g.my_rank == 0 ) {
    FREE( buffer[0].data, double, bar_size );
    FREE( buffer[1].data, double, bar_size );
  }
  
  norm = global_norm_double( (vector_double)phi, 0, l->inner_vector_size, l, no_threading );
  printf0("norm: %e\n", norm );
  printf0("...done (%lf seconds)\n\n", t1-t0 ); 
}
#else
void vector_io( double *phi, char *filename, const int mode, level_struct *l ) {
  double norm;
  // pointers to inner file structures.
  hid_t dset_id, xfer_id;
  // spaces, status
  hid_t memspace, filespace, status;
  // global (g_*) and local (l_*) sizes and start coordinates
  int *gl=l->global_lattice, *ll=l->local_lattice;
  hsize_t g_start[5], g_size[5], l_start[5], l_size[5];
  double fnEntryTime;
  fnEntryTime = MPI_Wtime();

  status = -1;

  // need to be already in thiseigenmodegroup (use stepIntoEigenmode)
  ASSERT(h5info.thiseigenmodegroup_id > 0);

  /*
   * set the parameters:
   */
  for ( hsize_t i = 0; i < 4; i++)
  {
    g_start[i] = g.my_coords[i] * ll[i];
    l_start[i] = 0;
    g_size[i] = gl[i];
    l_size[i] = ll[i];
  }
  g_start[4] = l_start[4] = 0; 
  g_size[4] = l_size[4] = 24; // spin * color * re/im

  // create a memspace:
  memspace = H5Screate_simple(5, l_size, NULL);
  ASSERT( memspace >= 0 );
  // create a filespace:
  filespace = H5Screate_simple(5, g_size, NULL);
  status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, g_start, (hsize_t*)NULL, l_size, (hsize_t*)NULL);
  ASSERT( status >= 0 );
  // create a data transfer property:
  xfer_id = H5Pcreate(H5P_DATASET_XFER);
  status = H5Pset_dxpl_mpio(xfer_id, H5FD_MPIO_COLLECTIVE);
  ASSERT( status >= 0);
  // check if dataset exists:
  htri_t dset_ex = H5Lexists(h5info.thiseigenmodegroup_id, "eigenvector", H5P_DEFAULT);
  ASSERT( dset_ex >= 0 );
  dset_id = -1;

  if ( mode ==  _READ ) {
    printf0("reading hdf5 file \"%s\"... ", filename);

    dset_id = H5Dopen(h5info.thiseigenmodegroup_id, "eigenvector", H5P_DEFAULT);
    ASSERT( dset_id >= 0 );
    //read data:
    status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, xfer_id, phi);
    ASSERT( status >= 0 );
  } else if ( mode == _WRITE ) {
    printf0("writing hdf5 file \"%s\"... ", filename);

    // write_header( file_id, NULL, filename, 1, l );
    printf0("write_header for g.method != 6 not implemented in HDF5!");
    
    if(dset_ex == 0)
    {
      /* dataset does not exist...
       * create it:
       */
      dset_id = H5Dcreate(h5info.thiseigenmodegroup_id, "eigenvector", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    else
    {
      /* dataset already exists...
       * open it for overwriting it...
       */
      dset_id = H5Dopen(h5info.thiseigenmodegroup_id, "eigenvector", H5P_DEFAULT);
    }
    ASSERT( dset_id >= 0 );
    
    // write data:
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, xfer_id, phi);
    ASSERT( status >= 0 );
  }
  /*
   * close all open handles:
   */
  status = H5Dclose(dset_id);
  ASSERT( status >= 0);
  status = H5Sclose(filespace);
  ASSERT( status >= 0);
  status = H5Sclose(memspace);
  ASSERT( status >= 0);
  status = H5Pclose(xfer_id);
  ASSERT( status >= 0);

  h5info.ioTime += MPI_Wtime()-fnEntryTime;

  norm = global_norm_double( (vector_double)phi, 0, l->inner_vector_size, l, no_threading );
  printf0("norm: %e\n", norm );
  printf0("...done\n");
}
#endif


#ifndef HAVE_HDF5
void vector_io_single_file( double *psi, double *lambda, char *filename, const int mode, int n, char *vector_type, level_struct *l ) {
  
  int t, z, y, x, *gl=l->global_lattice, *ll=l->local_lattice, bar_size = 24*ll[X], desired_rank, j;
  double t0, t1;
  double *phi_pt = NULL;
  double *phi = NULL;
  FILE* file = NULL;
  MPI_Request sreq, rreq;
  confbuffer_struct buffer[2];
  confbuffer_struct *buffer_pt = NULL;
    
  buffer[0].next = &(buffer[1]);
  buffer[1].next = &(buffer[0]);  
  buffer[0].data = NULL;
  buffer[1].data = NULL;

#ifdef BIG_ENDIAN_TV
  int i;
#endif
  
  if ( g.my_rank == 0 ) {
    MALLOC( buffer[0].data, double, bar_size );
    MALLOC( buffer[1].data, double, bar_size );
    buffer_pt = &buffer[0];
  }
  
  t0 = MPI_Wtime();
  
  if ( mode == _READ ) {

    if ( g.my_rank == 0 ) {
      ASSERT( (file = fopen( filename, "rb" )) != NULL );
      printf0("reading from file \"%s\" ...\n", filename );

      char *cur_line = NULL;
      MALLOC( cur_line, char, STRINGLENGTH );
      do {
        ASSERT( fgets( cur_line, STRINGLENGTH-1, file ) );
      } while ( strcmp( cur_line, "</header>\n" ) );
      FREE( cur_line, char, STRINGLENGTH );
    }
    
    for ( j=0; j<n; j++ ) {
      if ( g.my_rank == 0 ) {
        ASSERT( fread( buffer_pt->data, sizeof(double), bar_size, file ) );
      }

      phi=(double *) (l->x);
      phi_pt=phi;
      for ( t=0; t<gl[T]; t++ )
        for ( z=0; z<gl[Z]; z++ )
          for ( y=0; y<gl[Y]; y++ )
            for ( x=0; x<gl[X]; x+=ll[X] ) {
              
              desired_rank = process_index( t, z, y, x, ll );
              
              if ( g.my_rank == desired_rank ) {
                MPI_Irecv( phi_pt, bar_size, MPI_DOUBLE, 0, 0, g.comm_cart, &rreq );
              }
              if ( g.my_rank == 0 ) {
                MPI_Isend( buffer_pt->data, bar_size, MPI_DOUBLE, desired_rank, 0, g.comm_cart, &sreq );
              }
              
              if ( g.my_rank == 0 ) {
                if ( ! ( t == gl[T]-1 && z == gl[Z]-1 && y == gl[Y]-1 && x == gl[X]-ll[X] ) )
                  ASSERT( fread( buffer_pt->next->data, sizeof(double), bar_size, file ) );
              }
              
              if ( g.my_rank == desired_rank ) {
                MPI_Wait( &rreq, MPI_STATUS_IGNORE );
#ifdef BIG_ENDIAN_TV
                for ( i=0; i<bar_size; i++ ) {
                  byteswap8( (char *) ( phi_pt + i ) );
                }
#endif
                phi_pt += bar_size;
              }
              if ( g.my_rank == 0 ) {
                MPI_Wait( &sreq, MPI_STATUS_IGNORE );
                buffer_pt = buffer_pt->next;
              }
            }
      if ( psi == NULL ) {
        if ( g.mixed_precision )
          trans_float(l->is_float.test_vector[j], l->x, l->s_float.op.translation_table, l, no_threading);
        else
          trans_double(l->is_double.test_vector[j], l->x, l->s_double.op.translation_table, l, no_threading);
      } else {
        vector_double_copy( ((vector_double)psi)+j*l->inner_vector_size, l->x, 0, l->inner_vector_size, l );
      }
    }
  } else if ( mode == _WRITE ) {
    
    if ( g.my_rank == 0 ) {
      ASSERT( (file = fopen( filename, "wb" )) != NULL );
      write_header( &file, lambda, vector_type, n, l );
    }

    printf0("writing file \"%s\" ...\n", filename );
    
    for ( j=0; j<n; j++ ){
      if ( psi == NULL ) {
        if ( g.mixed_precision )
          trans_back_float( l->x, l->is_float.test_vector[j], l->s_float.op.translation_table, l, no_threading );
        else
          trans_back_double( l->x, l->is_double.test_vector[j], l->s_double.op.translation_table, l, no_threading );
      } else {
        vector_double_copy( l->x, ((complex_double*)psi)+j*l->inner_vector_size, 0, l->inner_vector_size, l );
      }
      phi=(double *)(l->x);
      phi_pt=phi;
      for ( t=0; t<gl[T]; t++ )
        for ( z=0; z<gl[Z]; z++ )
          for ( y=0; y<gl[Y]; y++ )
            for ( x=0; x<gl[X]; x+=ll[X] ) {
              
              desired_rank = process_index( t, z, y, x, ll );
              
              if ( g.my_rank == 0 ) {
                MPI_Irecv( buffer_pt->next->data, bar_size, MPI_DOUBLE, desired_rank, 0, g.comm_cart, &rreq );
              }
              if ( g.my_rank == desired_rank ) {
                MPI_Isend( phi_pt, bar_size, MPI_DOUBLE, 0, 0, g.comm_cart, &sreq );
              }
              
              if ( g.my_rank == 0 ) {
                if ( ! ( t == 0 && z == 0 && y == 0 && x == 0  ) ) {
#ifdef BIG_ENDIAN_TV
                  for ( i=0; i<bar_size; i++ ) {
                    byteswap8( (char *) ( buffer_pt->data + i ) );
                  }
#endif
                  fwrite( buffer_pt->data, sizeof(double), bar_size, file );
                }
              }
              
              if ( g.my_rank == 0 ) {
                MPI_Wait( &rreq, MPI_STATUS_IGNORE );
                buffer_pt = buffer_pt->next;
              }
              if ( g.my_rank == desired_rank ) {
                MPI_Wait( &sreq, MPI_STATUS_IGNORE );
                phi_pt += bar_size;
              }
              
            }
            
      if ( g.my_rank == 0 ) {
        #ifdef BIG_ENDIAN_TV
        for ( i=0; i<bar_size; i++ ) {
          byteswap8( (char *) ( buffer_pt->data + i ) );
        }
        #endif
        fwrite( buffer_pt->data, sizeof(double), bar_size, file );
      }
    }
  } else
    ASSERT( mode == _READ || mode == _WRITE );

  if ( g.my_rank == 0 ) {
    fclose( file );
  }
  
  t1 = MPI_Wtime();
  
  if ( g.my_rank == 0 ) {
    FREE( buffer[0].data, double, bar_size );
    FREE( buffer[1].data, double, bar_size );
  }
  
  printf0("...done (%lf seconds)\n\n", t1-t0 ); 
}
#else
void vector_io_single_file( double *psi, double *lambda, char *filename, const int mode, int n, char *vector_type, level_struct *l ) {
  printf0("io.c: vector_io_single_file: method not implemented using hdf5.");
  return;
}
#endif


void d_dump( config_double D, level_struct *l ) {
  
  if ( g.num_processes == 1 ) {
    
    int i, x, y, z, t;
    
    i=0;
    for ( t=0; t<l->local_lattice[T]; t++ )
      for ( z=0; z<l->local_lattice[Z]; z++ )
        for ( y=0; y<l->local_lattice[Y]; y++ )
          for ( x=0; x<l->local_lattice[X]; x++ ) {
            printf("site_number=%d, t=%d, z=%d, y=%d, x=%d.\n", i/36, t, z, y, x );
            printf("dir: T\n");
            printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", creal(D[i+0]), cimag(D[i+0]), creal(D[i+1]), cimag(D[i+1]), creal(D[i+2]), cimag(D[i+2]));
            printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", creal(D[i+3]), cimag(D[i+3]), creal(D[i+4]), cimag(D[i+4]), creal(D[i+5]), cimag(D[i+5]));
            printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", creal(D[i+6]), cimag(D[i+6]), creal(D[i+7]), cimag(D[i+7]), creal(D[i+8]), cimag(D[i+8]));
            i+=9;
            printf("dir: Z\n");
            printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", creal(D[i+0]), cimag(D[i+0]), creal(D[i+1]), cimag(D[i+1]), creal(D[i+2]), cimag(D[i+2]));
            printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", creal(D[i+3]), cimag(D[i+3]), creal(D[i+4]), cimag(D[i+4]), creal(D[i+5]), cimag(D[i+5]));
            printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", creal(D[i+6]), cimag(D[i+6]), creal(D[i+7]), cimag(D[i+7]), creal(D[i+8]), cimag(D[i+8]));
            i+=9;
            printf("dir: Y\n");
            printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", creal(D[i+0]), cimag(D[i+0]), creal(D[i+1]), cimag(D[i+1]), creal(D[i+2]), cimag(D[i+2]));
            printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", creal(D[i+3]), cimag(D[i+3]), creal(D[i+4]), cimag(D[i+4]), creal(D[i+5]), cimag(D[i+5]));
            printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", creal(D[i+6]), cimag(D[i+6]), creal(D[i+7]), cimag(D[i+7]), creal(D[i+8]), cimag(D[i+8]));
            i+=9;
            printf("dir: X\n");
            printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", creal(D[i+0]), cimag(D[i+0]), creal(D[i+1]), cimag(D[i+1]), creal(D[i+2]), cimag(D[i+2]));
            printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", creal(D[i+3]), cimag(D[i+3]), creal(D[i+4]), cimag(D[i+4]), creal(D[i+5]), cimag(D[i+5]));
            printf("%+lf%+lfi %+lf%+lfi %+lf%+lfi\n", creal(D[i+6]), cimag(D[i+6]), creal(D[i+7]), cimag(D[i+7]), creal(D[i+8]), cimag(D[i+8]));
            i+=9;
            printf("------------------------------------------------------------\n");
          }
  }
}

