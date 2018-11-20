#!/bin/bash

# numpy is not available in the system Python
# MKL is set up in the intel/compilers module

if [[ $MEASURE = true ]]; then

    intel_root=/cm/shared/apps/intel

    # Remove gcc so that libstdc++ gets loaded again later.
    module unload gcc

    # Load Advisor before gcc to not pick up Advisor's old libstdc++.
    # Probably Advisor's library directory should not be in
    # LD_LIBRARY_PATH?
    module load intel/advisor/2018
    source $ADVISOR_XE_2018_DIR/advixe-vars.sh
    export ADVISOR_DIR=$ADVISOR_2018_DIR
    
    # Probably don't need Inspector (more a debugging tool)
    #module load intel/inspector/2018
    
    # Load Amplifier before gcc to not pick up Amplifier's old
    # libstdc++ (just in case).
    module load intel/vtune/2018
    source $VTUNE_AMPLIFIER_XE_2018_DIR/amplxe-vars.sh
    export VTUNE_AMPLIFIER_DIR=$VTUNE_AMPLIFIER_2018_DIR
    
    if [[ $ITAC = true ]]; then
	# This sets VT_ROOT
	source $intel_root/itac_2018/bin/itacvars.sh
    fi
fi

module load python/intel/2.7.14 fftw/intel/3.3.7 gsl/intel/2.4

# Load the compilers and MPI last to ensure that LD_LIBRARY_PATH is
# consistent between the plain, measure, and measure+itac versions of
# gandalf.  Not completely:  without measure or measure+itac, the imf,
# intlc, irng and svml libraries (for FFTW and GSL) are found in
# /cm/shared/apps/intel/compilers_and_libraries_2018; with measure or
# measure+itac, they are found in
# /cm/shared/apps/intel/compilers_and_libraries_2018.3.222.  These are
# more-or-less the same, at least for the libraries.
module load gcc/6.4 intel/compilers/18.0.3 intel/mpi/18.0.3

# Is this needed?  A hack to the Intel compilers script.
source ${MKLROOT%/mkl}/bin/compilervars.sh intel64

# This is generally needed, but is actually OK for the Intel MPI 18 on
# DiAL because MT MPI is the default (lib/libmpi.so is a symbolic link
# to lib/release_mt/libmpi.so).
source $I_MPI_ROOT/intel64/bin/mpivars.sh release_mt
