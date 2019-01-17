#!/bin/bash

#PBS -l nodes=1:ppn=36
#PBS -l walltime=01:00:00
#PBS -A dp047
#PBS -q devel

# Export everything.
set -a

# Measuring, not ITAC
unset ITAC
MEASURE=true

if [[ $MEASURE = true ]]; then
    if [[ $ITAC = true ]]; then
	gandalf_program=gandalf_measure_itac
    else
	gandalf_program=gandalf_measure
    fi
else
    gandalf_program=gandalf_plain
fi

cd ~/gandalf/tests

# This isn't the problem with MPI/
input=disc_short_10
sequence=112
rm -rf $input.$sequence
mkdir $input.$sequence
(
    cd $input.$sequence
    ln -s ../../eos.bell.cc.dat .
    ln -s ../../stellar.dat .
    
    . ../../modules.bash
    
    # Intel's KMP_AFFINITY is not needed, used the OpenMP variables.
    OMP_NUM_THREADS=1
    OMP_PROC_BIND=true
    # No need for OMP_PLACES=cores on DIaL - Hyper-Threading is
    # switched off.
    #OMP_PLACES=cores
    
    OMP_NESTED=true
    I_MPI_PIN_DOMAIN=omp
    #OMP_DISPLAY_ENV=true
    #KMP_SETTINGS=true
    
    # Using Intel compilers with OpenMP and MKL will ensure that MKL
    # does not oversubscribe in OpenMP parallel regions, so no need
    # for this.
    #MKL_NUM_THREADS=1
    
    # Note that Gandalf uses a recursive routine for dividing cells,
    # with 2 OpenMP threads per call (like a tree) until the depth
    # corresponds to 2^depth > OMP_NUM_THREADS (Nthreads in the
    # program).  That is, it over-subscribes threads, which means
    # there are 'extra' threads in the Amplifier display - they are
    # correct, and don't do much work.  I'm not sure if these threads
    # will be on any spare cores (e.g. when OMP_NUM_THREADS=18 on a
    # 36-core node).
    
    # Use $input.$sequence as the results dir so that it is easy to
    # distinguish in the Recent Results list in the Amplifier GUI.
    mpirun -n 1 -gtool "amplxe-cl -collect hotspots -knob enable-user-tasks=true -data-limit=0 -result-dir /home/dc-fili1/gandalf/intel/amplxe/projects/gandalf/$input.$sequence:0" /home/dc-fili1/gandalf/bin/gandalf_measure ../DiRAC/$input.dat &> output
)
