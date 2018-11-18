#!/bin/bash

# Export everything.
set -a

# Source the file modules.bash that is in the same directory as this
# file.  Only works if this file is run, not sourced.  Does not follow
# symlinks.
script_dir=$(dirname $0)
source $script_dir/modules.bash

## This probably always works, whether run or sourced, and also if run
## in a bash function, see the bash manpage..  This version also
## follows all symlinks, which may or may not be what you want.
#last=$((${#BASH_SOURCE[*]} - 1))
#script_name=${BASH_SOURCE[$last]}
#script_dir=$(dirname $(readlink -e $script_name))

# How do I use -gtool?  Is that any better/different from calling
# Advisor or VTune directly?
mpirun -np 1 advixe-cl -collect survey --project-dir ./disc_short.advisor ../bin/gandalf DiRAC/disc_short.dat; exit

mpirun -np 1 advixe-cl -collect dependencies -mark-up-list=451 -project-dir disc_short.advisor ../bin/gandalf DiRAC/disc_short.dat; exit

mpirun -np 1 advixe-cl -collect survey --project-dir ./GI_disc_short.advisor ../bin/gandalf DiRAC/GI_disc_short.dat; exit

mpirun -n 1 -gtool "advixe-cl -collect tripcounts -flop -stacks -project-dir /home/dc-fili1/gandalf/tests/GI_disc_short.advisor:0" ../bin/gandalf DiRAC/GI_disc_short.dat; exit


. ../modules.bash
input=disc_short
sequence=001
mpirun -n 1 -gtool "amplxe-cl -collect hpc-performance -knob sampling-interval=10 -knob enable-stack-collection=true -knob collect-memory-bandwidth=false -data-limit=0 -result-dir vtune/$input.$sequence:0" /home/dc-fili1/gandalf/bin/gandalf DiRAC/$input.dat
exit

. ../modules.bash
input=disc_short
sequence=002
mpirun -n 1 -gtool "amplxe-cl -collect hotspots -data-limit=0 -result-dir vtune/$input.$sequence:0" /home/dc-fili1/gandalf/bin/gandalf DiRAC/$input.dat
exit

. ../modules.bash
input=disc_short
sequence=004
# Intel's KMP_AFFINITY is not needed, used the OpenMP variables.
OMP_NUM_THREADS=18
OMP_PROC_BIND=true
# No need for OMP_PLACE=core on DIaL - Hyper-Threading is switched
# off.
OMP_PLACES=cores
# No need to do
#source mpivars.sh release_mt
# because MT MPI is the default (lib/libmpi.so is a symbolic link to
# lib/release_mt/libmpi.so)
I_MPI_PIN_DOMAIN=omp
mpirun -n 1 -gtool "amplxe-cl -collect hotspots -knob analyze-openmp=true -data-limit=0 -result-dir vtune/$input.$sequence:0" /home/dc-fili1/gandalf/bin/gandalf DiRAC/$input.dat
exit

. ../modules.bash
input=disc_short
sequence=005
# Intel's KMP_AFFINITY is not needed, used the OpenMP variables.
OMP_NUM_THREADS=18
OMP_PROC_BIND=true
# No need for OMP_PLACE=core on DIaL - Hyper-Threading is switched
# off.
OMP_PLACES=cores
# No need to do
#source mpivars.sh release_mt
# because MT MPI is the default (lib/libmpi.so is a symbolic link to
# lib/release_mt/libmpi.so)
I_MPI_PIN_DOMAIN=omp
mpirun -n 1 -gtool "amplxe-cl -collect hpc-performance -knob sampling-interval=10 -knob enable-stack-collection=true -knob collect-memory-bandwidth=false -data-limit=0 -result-dir vtune/$input.$sequence:0" /home/dc-fili1/gandalf/bin/gandalf DiRAC/$input.dat
exit
# This gave no OpenMP information.  Re-finalizing on the login node
# gives differnt results (!) but still no OpenMP information.

. ../modules.bash
input=disc_short
sequence=006
# Intel's KMP_AFFINITY is not needed, used the OpenMP variables.
OMP_NUM_THREADS=18
OMP_PROC_BIND=true
# No need for OMP_PLACE=core on DIaL - Hyper-Threading is switched
# off.
OMP_PLACES=cores
# No need to do
#source mpivars.sh release_mt
# because MT MPI is the default (lib/libmpi.so is a symbolic link to
# lib/release_mt/libmpi.so)
I_MPI_PIN_DOMAIN=omp
mpirun -n 1 -gtool "amplxe-cl -collect hpc-performance -knob enable-stack-collection=true -knob collect-memory-bandwidth=false -data-limit=0 -result-dir vtune/$input.$sequence:0" /home/dc-fili1/gandalf/bin/gandalf DiRAC/$input.dat
exit
# Worse - lots of paused regions.

. ../modules.bash
input=disc_short
sequence=007
# Intel's KMP_AFFINITY is not needed, used the OpenMP variables.
OMP_NUM_THREADS=18
OMP_PROC_BIND=true
# No need for OMP_PLACE=core on DIaL - Hyper-Threading is switched
# off.
OMP_PLACES=cores
# No need to do
#source mpivars.sh release_mt
# because MT MPI is the default (lib/libmpi.so is a symbolic link to
# lib/release_mt/libmpi.so)
I_MPI_PIN_DOMAIN=omp
mpirun -n 1 -gtool "amplxe-cl -collect hpc-performance -knob sampling-interval=50 -knob enable-stack-collection=true -knob collect-memory-bandwidth=false -data-limit=0 -result-dir vtune/$input.$sequence:0" /home/dc-fili1/gandalf/bin/gandalf DiRAC/$input.dat
exit
# Worse - only 3 OpenMP regions shown - should be 10.

input=disc_short
sequence=008
mkdir $input.$sequence
cd $input.$sequence
. ../../modules.bash
# Intel's KMP_AFFINITY is not needed, used the OpenMP variables.
OMP_NUM_THREADS=18
OMP_PROC_BIND=true
# No need for OMP_PLACE=core on DIaL - Hyper-Threading is switched
# off.
OMP_PLACES=cores
# No need to do
#source mpivars.sh release_mt
# because MT MPI is the default (lib/libmpi.so is a symbolic link to
# lib/release_mt/libmpi.so)
I_MPI_PIN_DOMAIN=omp
mpirun -n 1 -gtool "amplxe-cl -collect hpc-performance -knob sampling-interval=10 -knob collect-memory-bandwidth=false -data-limit=0 -result-dir ../vtune/$input.$sequence:0" /home/dc-fili1/gandalf/bin/gandalf ../DiRAC/$input.dat
exit

input=disc_short
sequence=009
mkdir $input.$sequence
cd $input.$sequence
. ../../modules.bash
# Intel's KMP_AFFINITY is not needed, used the OpenMP variables.
OMP_NUM_THREADS=18
OMP_PROC_BIND=true
# No need for OMP_PLACES=core on DIaL - Hyper-Threading is switched
# off.
OMP_PLACES=cores
# No need to do
#source mpivars.sh release_mt
# because MT MPI is the default (lib/libmpi.so is a symbolic link to
# lib/release_mt/libmpi.so)
I_MPI_PIN_DOMAIN=omp
mpirun -n 1 -gtool "amplxe-cl -collect hpc-performance -knob sampling-interval=10 -knob collect-memory-bandwidth=false -data-limit=0 -result-dir vtune:0" /home/dc-fili1/gandalf/bin/gandalf ../DiRAC/$input.dat
exit


input=disc_short
sequence=010
mkdir $input.$sequence
cd $input.$sequence
. ../../modules.bash
# Intel's KMP_AFFINITY is not needed, used the OpenMP variables.
OMP_NUM_THREADS=18
OMP_PROC_BIND=true
# No need for OMP_PLACES=core on DIaL - Hyper-Threading is switched
# off.
OMP_PLACES=cores
# No need to do
#source mpivars.sh release_mt
# because MT MPI is the default (lib/libmpi.so is a symbolic link to
# lib/release_mt/libmpi.so)
I_MPI_PIN_DOMAIN=omp
mpirun -n 1 -gtool "amplxe-cl -collect hotspots -data-limit=0 -result-dir vtune:0" /home/dc-fili1/gandalf/bin/gandalf ../DiRAC/$input.dat
exit

input=disc_short
sequence=011
mkdir $input.$sequence
cd $input.$sequence
. ../../modules.bash
# Intel's KMP_AFFINITY is not needed, used the OpenMP variables.
OMP_NUM_THREADS=18
OMP_PROC_BIND=true
# No need for OMP_PLACES=core on DIaL - Hyper-Threading is switched
# off.
OMP_PLACES=cores
# No need to do
#source mpivars.sh release_mt
# because MT MPI is the default (lib/libmpi.so is a symbolic link to
# lib/release_mt/libmpi.so)
OMP_NESTED=true
OMP_DISPLAY_ENV=true
KMP_SETTINGS=true
I_MPI_PIN_DOMAIN=omp
mpirun -n 1 -gtool "amplxe-cl -collect hotspots -data-limit=0 -result-dir vtune:0" /home/dc-fili1/gandalf/bin/gandalf ../DiRAC/$input.dat
exit

input=disc_short
sequence=012
mkdir $input.$sequence
cd $input.$sequence
. ../../modules.bash
# Intel's KMP_AFFINITY is not needed, used the OpenMP variables.
OMP_NUM_THREADS=18
OMP_PROC_BIND=true
# No need for OMP_PLACES=core on DIaL - Hyper-Threading is switched
# off.
OMP_PLACES=cores
# No need to do
#source mpivars.sh release_mt
# because MT MPI is the default (lib/libmpi.so is a symbolic link to
# lib/release_mt/libmpi.so)
OMP_NESTED=true
OMP_DISPLAY_ENV=true
KMP_SETTINGS=true
I_MPI_PIN_DOMAIN=omp
# Use $input.$sequence as the results dir so that it is easy to
# distinguish in the Recent Results list in the Amplifier GUI.
mpirun -n 1 -gtool "amplxe-cl -collect hotspots -knob analyze-openmp=true -data-limit=0 -result-dir $input.$sequence:0" /home/dc-fili1/gandalf/bin/gandalf ../DiRAC/$input.dat &> output
exit

input=disc_short
sequence=013
mkdir $input.$sequence
cd $input.$sequence
. ../../modules.bash
# Intel's KMP_AFFINITY is not needed, used the OpenMP variables.
OMP_NUM_THREADS=1
OMP_PROC_BIND=true
# No need for OMP_PLACES=core on DIaL - Hyper-Threading is switched
# off.
OMP_PLACES=cores
# No need to do
#source mpivars.sh release_mt
# because MT MPI is the default (lib/libmpi.so is a symbolic link to
# lib/release_mt/libmpi.so)
OMP_NESTED=true
OMP_DISPLAY_ENV=true
KMP_SETTINGS=true
I_MPI_PIN_DOMAIN=omp
# Use $input.$sequence as the results dir so that it is easy to
# distinguish in the Recent Results list in the Amplifier GUI.
mpirun -n 1 -gtool "amplxe-cl -collect hotspots -data-limit=0 -result-dir $input.$sequence:0" /home/dc-fili1/gandalf/bin/gandalf ../DiRAC/$input.dat &> output
exit

input=disc_short
sequence=014
mkdir $input.$sequence
cd $input.$sequence
. ../../modules.bash

# Intel's KMP_AFFINITY is not needed, used the OpenMP variables.
OMP_NUM_THREADS=18
OMP_PROC_BIND=true
# No need for OMP_PLACES=core on DIaL - Hyper-Threading is switched
# off.
OMP_PLACES=cores
# No need to do
#source mpivars.sh release_mt
# because MT MPI is the default (lib/libmpi.so is a symbolic link to
# lib/release_mt/libmpi.so)
OMP_NESTED=true
OMP_DISPLAY_ENV=true
KMP_SETTINGS=true
I_MPI_PIN_DOMAIN=omp
# Using Intel compilers with OpenMP and MKL will ensure that MKL does
# not oversubscribe in OpenMP paralle regions.  Check if this variable
# will remove the extra threads seen in Amplifier. No it doesn't.
# MKL_NUM_THREADS=1

# Use $input.$sequence as the results dir so that it is easy to
# distinguish in the Recent Results list in the Amplifier GUI.
mpirun -n 1 -gtool "amplxe-cl -collect hotspots -knob analyze-openmp=true -data-limit=0 -result-dir $input.$sequence:0" /home/dc-fili1/gandalf/bin/gandalf ../DiRAC/$input.dat &> output
exit

input=disc_short
sequence=015
mkdir $input.$sequence
(
cd $input.$sequence
. ../../modules.bash

# Intel's KMP_AFFINITY is not needed, used the OpenMP variables.
OMP_NUM_THREADS=18
OMP_PROC_BIND=true
# No need for OMP_PLACES=core on DIaL - Hyper-Threading is switched
# off.
OMP_PLACES=cores
# No need to do
#source mpivars.sh release_mt
# because MT MPI is the default (lib/libmpi.so is a symbolic link to
# lib/release_mt/libmpi.so)
OMP_NESTED=true
OMP_DISPLAY_ENV=true
KMP_SETTINGS=true
I_MPI_PIN_DOMAIN=omp
# Using Intel compilers with OpenMP and MKL will ensure that MKL does
# not oversubscribe in OpenMP paralle regions.  Check if this variable
# will remove the extra threads seen in Amplifier. No it doesn't.
# MKL_NUM_THREADS=1

# Use $input.$sequence as the results dir so that it is easy to
# distinguish in the Recent Results list in the Amplifier GUI.
mpirun -n 1 -gtool "amplxe-cl -collect hpc-performance -knob sampling-interval=10 -knob collect-memory-bandwidth=false -data-limit=0 -result-dir $input.$sequence:0" /home/dc-fili1/gandalf/bin/gandalf ../DiRAC/$input.dat &> output
)
exit

input=disc_short_100
sequence=016
mkdir $input.$sequence
(
cd $input.$sequence
. ../../modules.bash

# Intel's KMP_AFFINITY is not needed, used the OpenMP variables.
OMP_NUM_THREADS=18
OMP_PROC_BIND=true
# No need for OMP_PLACES=core on DIaL - Hyper-Threading is switched
# off.
OMP_PLACES=cores
# No need to do
#source mpivars.sh release_mt
# because MT MPI is the default (lib/libmpi.so is a symbolic link to
# lib/release_mt/libmpi.so)
OMP_NESTED=true
OMP_DISPLAY_ENV=true
KMP_SETTINGS=true
I_MPI_PIN_DOMAIN=omp
# Using Intel compilers with OpenMP and MKL will ensure that MKL does
# not oversubscribe in OpenMP paralle regions.  Check if this variable
# will remove the extra threads seen in Amplifier. No it doesn't.
# MKL_NUM_THREADS=1

# Use $input.$sequence as the results dir so that it is easy to
# distinguish in the Recent Results list in the Amplifier GUI.
mpirun -n 1 -gtool "amplxe-cl -collect hpc-performance -knob sampling-interval=10 -knob collect-memory-bandwidth=false -data-limit=0 -result-dir $input.$sequence:0" /home/dc-fili1/gandalf/bin/gandalf ../DiRAC/$input.dat &> output
)
exit

input=disc_short_10
sequence=017
mkdir $input.$sequence
(
cd $input.$sequence
. ../../modules.bash

# Intel's KMP_AFFINITY is not needed, used the OpenMP variables.
OMP_NUM_THREADS=18
OMP_PROC_BIND=true
# No need for OMP_PLACES=core on DIaL - Hyper-Threading is switched
# off.
OMP_PLACES=cores
# No need to do
#source mpivars.sh release_mt
# because MT MPI is the default (lib/libmpi.so is a symbolic link to
# lib/release_mt/libmpi.so)
OMP_NESTED=true
OMP_DISPLAY_ENV=true
KMP_SETTINGS=true
I_MPI_PIN_DOMAIN=omp
# Using Intel compilers with OpenMP and MKL will ensure that MKL does
# not oversubscribe in OpenMP paralle regions.  Check if this variable
# will remove the extra threads seen in Amplifier. No it doesn't.
# MKL_NUM_THREADS=1

# Use $input.$sequence as the results dir so that it is easy to
# distinguish in the Recent Results list in the Amplifier GUI.
mpirun -n 1 -gtool "amplxe-cl -collect hotspots -data-limit=0 -result-dir $input.$sequence:0" /home/dc-fili1/gandalf/bin/gandalf ../DiRAC/$input.dat &> output
)

input=disc_short_10
sequence=018
mkdir $input.$sequence
(
cd $input.$sequence
. ../../modules.bash

# Intel's KMP_AFFINITY is not needed, used the OpenMP variables.
OMP_NUM_THREADS=18
OMP_PROC_BIND=true
# No need for OMP_PLACES=core on DIaL - Hyper-Threading is switched
# off.
OMP_PLACES=cores
# No need to do
#source mpivars.sh release_mt
# because MT MPI is the default (lib/libmpi.so is a symbolic link to
# lib/release_mt/libmpi.so)
OMP_NESTED=true
OMP_DISPLAY_ENV=true
KMP_SETTINGS=true
I_MPI_PIN_DOMAIN=omp
# Using Intel compilers with OpenMP and MKL will ensure that MKL does
# not oversubscribe in OpenMP paralle regions.  Check if this variable
# will remove the extra threads seen in Amplifier. No it doesn't.
# MKL_NUM_THREADS=1

# Use $input.$sequence as the results dir so that it is easy to
# distinguish in the Recent Results list in the Amplifier GUI.
mpirun -n 1 -gtool "amplxe-cl -collect hotspots -knob analyze-openmp=true -data-limit=0 -result-dir $input.$sequence:0" /home/dc-fili1/gandalf/bin/gandalf ../DiRAC/$input.dat &> output
)

input=disc_short_10
sequence=019
mkdir $input.$sequence
(
cd $input.$sequence
. ../../modules.bash

# Intel's KMP_AFFINITY is not needed, used the OpenMP variables.
OMP_NUM_THREADS=18
OMP_PROC_BIND=true
# No need for OMP_PLACES=core on DIaL - Hyper-Threading is switched
# off.
OMP_PLACES=cores
# No need to do
#source mpivars.sh release_mt
# because MT MPI is the default (lib/libmpi.so is a symbolic link to
# lib/release_mt/libmpi.so)
OMP_NESTED=true
OMP_DISPLAY_ENV=true
KMP_SETTINGS=true
I_MPI_PIN_DOMAIN=omp
# Using Intel compilers with OpenMP and MKL will ensure that MKL does
# not oversubscribe in OpenMP paralle regions.  Check if this variable
# will remove the extra threads seen in Amplifier. No it doesn't.
# MKL_NUM_THREADS=1

# Use $input.$sequence as the results dir so that it is easy to
# distinguish in the Recent Results list in the Amplifier GUI.
mpirun -n 1 -gtool "amplxe-cl -collect hpc-performance -knob sampling-interval=10 -knob enable-stack-collection=true -knob collect-memory-bandwidth=false -knob analyze-openmp=false -data-limit=0 -result-dir $input.$sequence:0" /home/dc-fili1/gandalf/bin/gandalf ../DiRAC/$input.dat &> output
)

input=disc_short_10
sequence=020
mkdir $input.$sequence
(
cd $input.$sequence
. ../../modules.bash

# Intel's KMP_AFFINITY is not needed, used the OpenMP variables.
OMP_NUM_THREADS=18
OMP_PROC_BIND=true
# No need for OMP_PLACES=core on DIaL - Hyper-Threading is switched
# off.
OMP_PLACES=cores
# No need to do
#source mpivars.sh release_mt
# because MT MPI is the default (lib/libmpi.so is a symbolic link to
# lib/release_mt/libmpi.so)
OMP_NESTED=true
OMP_DISPLAY_ENV=true
KMP_SETTINGS=true
I_MPI_PIN_DOMAIN=omp
# Using Intel compilers with OpenMP and MKL will ensure that MKL does
# not oversubscribe in OpenMP paralle regions.  Check if this variable
# will remove the extra threads seen in Amplifier. No it doesn't.
# MKL_NUM_THREADS=1

# Use $input.$sequence as the results dir so that it is easy to
# distinguish in the Recent Results list in the Amplifier GUI.
mpirun -n 1 -gtool "amplxe-cl -collect hpc-performance -knob sampling-interval=10 -knob enable-stack-collection=true -knob collect-memory-bandwidth=false -data-limit=0 -result-dir $input.$sequence:0" /home/dc-fili1/gandalf/bin/gandalf ../DiRAC/$input.dat &> output
)

input=disc_short_100
sequence=021
mkdir $input.$sequence
(
cd $input.$sequence
. ../../modules.bash

# Intel's KMP_AFFINITY is not needed, used the OpenMP variables.
OMP_NUM_THREADS=18
OMP_PROC_BIND=true
# No need for OMP_PLACES=core on DIaL - Hyper-Threading is switched
# off.
OMP_PLACES=cores
# No need to do
#source mpivars.sh release_mt
# because MT MPI is the default (lib/libmpi.so is a symbolic link to
# lib/release_mt/libmpi.so)
OMP_NESTED=true
OMP_DISPLAY_ENV=true
KMP_SETTINGS=true
I_MPI_PIN_DOMAIN=omp
# Using Intel compilers with OpenMP and MKL will ensure that MKL does
# not oversubscribe in OpenMP paralle regions.  Check if this variable
# will remove the extra threads seen in Amplifier. No it doesn't.
# MKL_NUM_THREADS=1

# Use $input.$sequence as the results dir so that it is easy to
# distinguish in the Recent Results list in the Amplifier GUI.
mpirun -n 1 -gtool "amplxe-cl -collect hotspots -knob analyze-openmp=true -data-limit=0 -result-dir $input.$sequence:0" /home/dc-fili1/gandalf/bin/gandalf ../DiRAC/$input.dat &> output
)

input=disc_short_100
sequence=022
mkdir $input.$sequence
(
cd $input.$sequence
. ../../modules.bash

# Intel's KMP_AFFINITY is not needed, used the OpenMP variables.
OMP_NUM_THREADS=18
OMP_PROC_BIND=true
# No need for OMP_PLACES=core on DIaL - Hyper-Threading is switched
# off.
OMP_PLACES=cores
# No need to do
#source mpivars.sh release_mt
# because MT MPI is the default (lib/libmpi.so is a symbolic link to
# lib/release_mt/libmpi.so)
OMP_NESTED=true
OMP_DISPLAY_ENV=true
KMP_SETTINGS=true
I_MPI_PIN_DOMAIN=omp
# Using Intel compilers with OpenMP and MKL will ensure that MKL does
# not oversubscribe in OpenMP paralle regions.  Check if this variable
# will remove the extra threads seen in Amplifier. No it doesn't.
# MKL_NUM_THREADS=1

# Use $input.$sequence as the results dir so that it is easy to
# distinguish in the Recent Results list in the Amplifier GUI.
mpirun -n 1 -gtool "amplxe-cl -collect hpc-performance -knob sampling-interval=10 -knob enable-stack-collection=true -knob collect-memory-bandwidth=false -data-limit=0 -result-dir $input.$sequence:0" /home/dc-fili1/gandalf/bin/gandalf ../DiRAC/$input.dat &> output
)

input=disc_short_100
sequence=023
mkdir $input.$sequence
(
cd $input.$sequence
. ../../modules.bash

# Intel's KMP_AFFINITY is not needed, used the OpenMP variables.
OMP_NUM_THREADS=36
OMP_PROC_BIND=true
# No need for OMP_PLACES=cores on DIaL - Hyper-Threading is switched
# off.
OMP_PLACES=cores
# No need to do
#source mpivars.sh release_mt
# because MT MPI is the default (lib/libmpi.so is a symbolic link to
# lib/release_mt/libmpi.so)
OMP_NESTED=true
OMP_DISPLAY_ENV=true
KMP_SETTINGS=true
I_MPI_PIN_DOMAIN=omp
# Using Intel compilers with OpenMP and MKL will ensure that MKL does
# not oversubscribe in OpenMP parallel regions, so no need for this.
# MKL_NUM_THREADS=1

# Note that Gandalf uses a recursive routine for dividing cells, with
# 2 OpenMP threads per call (like a tree) until the depth corresponds
# to 2^depth > OMP_NUM_THREADS (Nthreads in the program).  That is, it
# over-subscribes threads, which means there are 'extra' threads in
# the Amplifier display - they are correct, and don't do much work.
# I'm nor sure if these threads will be on any spare cores (e.g. if
# OMP_NUM_THREADS=18 on a 36-core node).

# Use $input.$sequence as the results dir so that it is easy to
# distinguish in the Recent Results list in the Amplifier GUI.
mpirun -n 1 -gtool "amplxe-cl -collect hotspots -knob analyze-openmp=true -data-limit=0 -result-dir $input.$sequence:0" /home/dc-fili1/gandalf/bin/gandalf ../DiRAC/$input.dat &> output
)
