
The DDalphaAMG solver library is an inverter for Wilson-Clover fermions from lattice QCD.

## INSTALL:

  The main directory contains a makefile example. The standard makefile
  should work with a few adjustments (compiler and MPI library) on the 
  machine of your choice. Once the necessary   adjustments of the
  makefile are done, the command
  
    "make -f yourmakefile -j numberofthreads wilson"
    
  should compile the whole entire Wilson solver code, such that it is
  ready to be run. The makefile contains additional rules "library" and
  "documentation" for compiling the code as a library and for compiling
  the user documentation. the library interface can be found in the
  "include" folder.

## HOWTO:

  After having compiled the user documentation via
  "make documentation" please consult the compiled PDF in /doc for
  further information.

## NOTE:
 
 This repository has been forked multiple times for different purposes and uses. In particular, find:
 
 -- GPU improvements (Wilson) here : https://github.com/Gustavroot/DDalphaAMG
 
 -- coarsest-level improvements (Twisted Mass) here : https://github.com/JesusEV/DDalphaAMG_ci
