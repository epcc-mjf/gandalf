#!/bin/bash

# numpy is not available in the system Python
# MKL is set up in the intel/compilers module

module load gcc/6.4 python/intel/2.7.14 intel/compilers/18.0.3 fftw/intel/3.3.7 gsl/intel/2.4 intel/mpi/18.0.3

# Is this needed?  A hack to the Intel compilers script.
source ${MKLROOT%/mkl}/bin/compilervars.sh intel64

# This is genrally needed, but is actually OK for the Intel MPI 18 on
# DiAL because MT MPI is the default (lib/libmpi.so is a symbolic link
# to lib/release_mt/libmpi.so).
source $I_MPI_ROOT/intel64/bin/mpivars.sh release_mt
