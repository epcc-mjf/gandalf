#==================================================================================================
#  GANDALF v0.4.0 Makefile
#
#  This file is part of GANDALF :
#  Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids
#  https://github.com/gandalfcode/gandalf
#  Contact : gandalfcode@gmail.com
#
#  Copyright (C) 2013  D. A. Hubber, G. Rosotti
#
#  GANDALF is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#
#  GANDALF is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License (http://www.gnu.org/licenses) for more details.
#==================================================================================================


ifneq ($(PYTHON),)
ifeq ($(NUMPY),)
NUMPY = $(shell $(PYTHON) -c "import numpy; print numpy.get_include()")
endif
ifeq ($(PYLIB),)
PYLIB = $(shell $(PYTHON) -c "import distutils.sysconfig; print distutils.sysconfig.get_python_inc()")
endif
endif


#location of the python files produced by swig
PYFOLDER = ../analysis/swig_generated
#GTEST = $(GTEST_DIR)
#GTEST = /Users#/david/astro/code/gtest-1.6.0/include


VPATH = $(PWD)/src:$(PWD)/src/Common:$(PWD)/src/Feedback
VPATH +=$(PWD)/src/GradhSph:$(PWD)/src/Headers:$(PWD)/src/Hydrodynamics
VPATH +=$(PWD)/src/Ic:$(PWD)/src/MeshlessFV:$(PWD)/src/Mpi:$(PWD)/src/Nbody
VPATH +=$(PWD)/src/Radiation:$(PWD)/src/SM2013:$(PWD)/src/Thermal:$(PWD)/src/Tree
VPATH +=$(PWD)/src/UnitTesting
# Use GANDALF_INCLUDE instead of INCLUDE to avoid namespace clash on DIaL
GANDALF_INCLUDE += -I$(PWD)/src/Headers

ifeq ($(strip $(CPP)),)
  CPP=$(CXX)
endif

ifneq ($(strip $(MPICPP)),)
  CPP=$(MPICPP)
endif

# Compiler mode flags
# -------------------------------------------------------------------------------------------------
# For Intel, look at the assembly code with
# -Fa -fsource-asm -fverbose-asm -qopt-report-embed -masm=intel

ifeq ($(COMPILER_MODE),FAST)
GCCOPT   += -O3 -ffast-math -fPIC -fno-exceptions -fno-rtti
# -fast implies -static, which won't find, e.g., libm on DIaL.  And
# -static will conflict with the default -qopenmp-link dynamic for
# -OpenMP.
# The DIaL login nodes are the same as the compute nodes, so -xHost is
# OK when compiling on the login nodes.
# Switch on symbols so that gdb backtrace will work.  Switch on traps
# until seg faults are found.
# for testing the MPI problem ICCOPT   += -O3 -g -ipo -no-prec-div -fp-model fast=2 -xHost -fPIC -fno-exceptions -fp-trap=common
# Try without ipo, and show the assembly code; then it is clear what
# the vectorisation is.
#ICCOPT   += -O3 -ipo -no-prec-div -fp-model fast=2 -xHost -fPIC -fno-exceptions -fno-rtti
ICCOPT   += -O3 -no-prec-div -fp-model fast=2 -xHost -fPIC -fno-exceptions -fno-rtti -Fa -masm=intel -fsource-asm -qopt-report-embed
CLANGOPT += -O3 -ffast-math -fPIC
else ifeq ($(COMPILER_MODE),DEBUG)
GCCOPT   += -Og -g -pg -Wall -Wno-unknown-pragmas -Wno-reorder -fbounds-check -fPIC  -fno-rtti #-fprofile-arcs -ftest-coverage
ICCOPT   += -O0 -g -pg -Wall -Wno-unknown-pragmas -Wno-reorder -xHost -fPIC -fno-rtti
CLANGOPT += -O3 -g -pg -Wall
else ifeq ($(COMPILER_MODE),STANDARD)
GCCOPT   += -O3 -fPIC -fno-exceptions -fno-rtti
ICCOPT   += -O3 -fPIC -fno-exceptions -fno-rtti
else ifeq ($(COMPILER_MODE),PROFILE)
GCCOPT   += -O3 -fPIC -fno-exceptions -fno-rtti -pg
ICCOPT   += -O3 -fPIC -fno-exceptions -fno-rtti -pg
endif

ifeq ($(OPENMP),1)
GCCOPT   += -fopenmp
ICCOPT   += -qopenmp
CLANGOPT += -fopenmp -DOPENMP
else
ICCOPT   += -Wno-unknown-pragmas
endif

# Need to know if MPI is used before possibly adding ITAC
ifeq ($(findstring mpi,$(CPP)),mpi)
MPI = 1
endif

ifeq ($(MEASURE),true)
ICCOPT += -D MEASURE

# Reports
ICCOPT += -qopt-report=5 -qopt-report-phase=vec # -qopt-report-names=unmangled

# Advisor
# -vec is set by -O3
# -simd is set by default
# -qopenmp is set already if OpenMP is used
ICCOPT += -g -debug inline-debug-info -I$(ADVISOR_DIR)/include
# This option is needed for VTune (and maybe Advisor?), but is
# included with -g:  -fno-omit-frame-pointer

# VTune
ICCOPT += -I$(VTUNE_AMPLIFIER_DIR)/include
ICCOPT += -shared-intel
# TBB is not used
#ICCOPT +=-D TBB_USE_THREADING_TOOLS
ICCOPT += -parallel-source-info=2
# See if no inlining helps decipher Advisor output (at a cost!).  None
# of these prevent inlining - why?
# ICCOPT += -fno-inline -fno-inline-functions -ip-no-inlining
# .a files should be found (after .so files but there is only a
# libittnotify.a) but the VTune library directory is only in
# LD_LIBRARY_PATH - would it be seen at the link step?
LIB += -L$(VTUNE_AMPLIFIER_DIR)/lib64 -littnotify

# ITAC
ifeq ($(ITAC),true)
ifeq ($(MPI),1)
ifeq ($(CPP),mpiicpc)
# ITAC
ICCOPT += -D ITAC
ICCOPT += -I$(VT_ROOT)/include
ICCOPT += -tcollect -trace
# No need to use mpirun -trace to trace, since the ITAC libraries are
# statically linked.
endif
endif
else
ICCOPT += -U ITAC
endif
else
# switch off everything
ICCOPT += -U MEASURE -U ITAC
# Or could do:
# Switch off Advisor API calls with -D ANNOTATE_EXPAND_NULL
# Switch off VTune Amplifier API calls with -D INTEL_NO_ITTNOTIFY_API
# Switch off ITAC API calls by linking libVTnull instead of libVT
#ICCOPT += -D ANNOTATE_EXPAND_NULL
#ICCOPT += -D INTEL_NO_ITTNOTIFY_API
endif

# Select the compiler
# -------------------------------------------------------------------------------------------------
ifeq ($(findstring icpc,$(CPP)),icpc)
OPT += $(ICCOPT)
else ifeq ($(findstring clang,$(CPP)),clang)
OPT += $(CLANGOPT)
else ifeq ($(findstring g++,$(CPP)),g++)
OPT += $(GCCOPT)
else ifeq ($(CPP),mpic++)
# mpic++ could be icpc, need to test for this somehow.
OPT += $(GCCOPT)
endif

ifeq ($(MPI),1)
CFLAGS += -DMPI_PARALLEL
ifeq ($(CPP),mpiicpc)
# libraries are added by mpiicpc
else ifneq ($(CPP), mpic++)
LIB += -lmpi -lmpi++
endif
endif


ifeq ($(CXX),clang++)
 OPT += -D__extern_always_inline=inline
endif

# Precision options
# -------------------------------------------------------------------------------------------------
ifeq ($(PRECISION),SINGLE)
CFLAGS += -DGANDALF_SINGLE_PRECISION
else ifeq ($(PRECISION),DOUBLE)
CFLAGS += -DGANDALF_DOUBLE_PRECISION
endif

# Snapshot precision
# -------------------------------------------------------------------------------------------------
ifeq ($(PYSNAP_PRECISION),SINGLE)
CFLAGS += -DGANDALF_SNAPSHOT_SINGLE_PRECISION
endif


# FFTW linkage
#--------------------------------------------------------------------------------------------------
ifeq ($(FFTW),1)
CFLAGS += -DFFTW_TURBULENCE
ifeq ($(FFTW_INCLUDE),)
GANDALF_INCLUDE += -I/sw/include
else
GANDALF_INCLUDE += $(FFTW_INCLUDE)
endif
ifeq ($(FFTW_LIBRARY),)
LIB += -lfftw3
else
LIB += $(FFTW_LIBRARY)
endif
endif


# GNU scientific library linkage
#--------------------------------------------------------------------------------------------------
ifeq ($(GSL),1)
CFLAGS += -DGANDALF_GSL
ifeq ($(GSL_INCLUDE),)
GANDALF_INCLUDE += -I/usr/local/include
else
GANDALF_INCLUDE += $(GSL_INCLUDE)
endif
ifeq ($(GSL_LIBRARY),)
ifeq ($(MKL),1)
# Use MKL BLAS
LIB += -L/usr/local/lib -lgsl
else
LIB += -L/usr/local/lib -lgsl -lgslcblas
endif
else
LIB += $(GSL_LIBRARY)
endif
endif


# MKL linkage
#--------------------------------------------------------------------------------------------------
ifeq ($(MKL),1)
CFLAGS += -DGANDALF_MKL
ifeq ($(MKL_INCLUDE),)
# Do nothing
else
GANDALF_INCLUDE += $(MKL_INCLUDE)
endif
ifeq ($(MKL_LIBRARY),)
# Do nothing
else
LIB += $(MKL_LIBRARY)
endif
endif


# Debug output flags
# -------------------------------------------------------------------------------------------------
ifeq ($(OUTPUT_LEVEL),1)
CFLAGS += -DDEBUG1
else ifeq ($(OUTPUT_LEVEL),2)
CFLAGS += -DDEBUG1 -DDEBUG2
else ifeq ($(OUTPUT_LEVEL),3)
CFLAGS += -DDEBUG1 -DDEBUG2 -DDEBUG3
endif


# Include expensive verification code
# -------------------------------------------------------------------------------------------------
ifeq ($(DEBUG_LEVEL),1)
CFLAGS +=
else ifeq ($(DEBUG_LEVEL),2)
CFLAGS += -DVERIFY_ALL
else
CFLAGS += -DNDEBUG
endif


# Reorder particles flag
#ifeq ($(REORDER_PARTICLES),1)
#CFLAGS += -DREORDER_PARTICLES
#endif



# Object files to be compiled
# -------------------------------------------------------------------------------------------------
SWIG_HEADERS = Simulation.i
WRAP_OBJ = Simulation_wrap.o
OBJ = Parameters.o SimUnits.o
OBJ += Simulation.o Hydrodynamics.o
OBJ += SphSimulation.o Sph.o
OBJ += GradhSphSimulation.o GradhSph.o
#OBJ += SM2012SphSimulation.o SM2012Sph.o
#OBJ += MeshlessFVSimulation.o FV.o MeshlessFV.o MfvCommon.o
#OBJ += MfvMusclSimulation.o MfvMuscl.o
#OBJ += MfvRungeKuttaSimulation.o MfvRungeKutta.o
OBJ += NbodySimulation.o
OBJ += M4Kernel.o QuinticKernel.o GaussianKernel.o TabulatedKernel.o
OBJ += Integration.o SphIntegration.o SphLeapfrogKDK.o SphLeapfrogDKD.o
#OBJ += MfvIntegration.o
OBJ += RiemannSolver.o
OBJ += SphNeighbourSearch.o
OBJ += HydroTree.o Tree.o KDTree.o OctTree.o BruteForceTree.o
#OBJ += MeshlessFVTree.o
OBJ += GradhSphTree.o
#OBJ += SM2012SphTree.o
OBJ += Ewald.o
OBJ += AdiabaticEOS.o BarotropicEOS.o Barotropic2EOS.o
OBJ += PolytropicEOS.o IsothermalEOS.o RadwsEOS.o LocallyIsothermal.o DiscLocallyIsothermal.o RadiativeFB.o
OBJ += BetaCooling.o
OBJ += IonisingRadiationEOS.o
OBJ += MultipleSourceIonisation.o
OBJ += KDRadiationTree.o
OBJ += chealpix.o TreeRay.o TreeRayOnTheSpot.o
OBJ += EnergyEquation.o EnergyRadws.o OpacityTable.o
OBJ += Nbody.o NbodyLeapfrogKDK.o NbodyLeapfrogDKD.o
OBJ += NbodyHermite4.o NbodyHermite4TS.o NbodyHermite6TS.o
OBJ += NbodySystemTree.o
OBJ += Sinks.o
OBJ += Ghosts.o
OBJ += SphSnapshot.o
OBJ += CodeTiming.o
#OBJ += Dust.o
OBJ += Particle.o RandomNumber.o
OBJ += Supernova.o SupernovaDriver.o
OBJ += Ic.o BasicIc.o BinaryAccretionIc.o BlobIc.o BondiAccretionIc.o
OBJ += BossBodenheimerIc.o ContactDiscontinuityIc.o DiscIc.o DustyBoxIc.o
OBJ += EvrardCollapseIc.o EwaldIc.o FilamentIc.o GaussianRingIc.o GreshoVortexIc.o
OBJ += HierarchicalSystemIc.o IsothermalSphereIc.o KelvinHelmholtzIc.o NohIc.o PlummerSphereIc.o
OBJ += PolytropeIc.o RayleighTaylorIc.o SedovBlastwaveIc.o ShearflowIc.o ShocktubeIc.o Shock2DIc.o
OBJ += SilccIc.o SoundwaveIc.o SpitzerExpansionIc.o TurbulentCoreIc.o UniformIc.o ICRegularization.o

ifeq ($(MPI),1)
OBJ += MpiNode.o MpiTree.o MpiControl.o MpiKDTreeDecomposition.o
endif


TEST_OBJ = TestScaling.o TestTree.o TestOctTree.o TestSinks.o

.SUFFIXES: .cpp .i .o


%_wrap.cxx: %.i
	swig -c++ -python -outdir $(PYFOLDER) $(CFLAGS) $(GANDALF_INCLUDE) $<

%_wrap.o: %_wrap.cxx
	$(CPP) $(OPT) $(CFLAGS) -fexceptions -frtti -c $< -I$(PYLIB) -I$(PYLIB)/config -I$(NUMPY) $(GANDALF_INCLUDE)

%.o: %.cxx
	$(CPP) $(CFLAGS) $(GANDALF_INCLUDE) -c $< $(OPT)

%.o: %.cpp
	$(CPP) $(OPT) $(CFLAGS) -I$(GTEST)/include $(GANDALF_INCLUDE) -c $<

%.o: %.c
	$(CPP) $(OPT) $(CFLAGS) -I$(GTEST)/include $(GANDALF_INCLUDE) -c $<

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
SHARED_OPTIONS = -bundle -flat_namespace -undefined suppress
else ifeq ($(UNAME_S),Linux)
SHARED_OPTIONS = -shared
endif
.SECONDARY:

#==================================================================================================
all : $(WRAP_OBJ) $(OBJ) Exception.o _SphSim.so executable
	cp _SphSim.so ../analysis/swig_generated

#==================================================================================================

executable : $(OBJ) gandalf.o Exception.o
	$(CPP) $(CFLAGS) $(GANDALF_INCLUDE) -o gandalf $(OBJ) Exception.o gandalf.o $(OPT) $(LIB)
	cp gandalf ../bin/gandalf

unittests : $(OBJ) $(TEST_OBJ) Exception.o
	$(CPP) $(CFLAGS) $(GANDALF_INCLUDE) $(OPT) $(LIB) -o testgandalf $(OBJ) $(TEST_OBJ) Exception.o $(GTEST)/make/gtest_main.a
	cp testgandalf ../bin/testgandalf

_SphSim.so : $(WRAP_OBJ) $(OBJ) Exception.o Render.o
	$(CPP) $(CFLAGS) $(GANDALF_INCLUDE) $(OPT) $(SHARED_OPTIONS) $(WRAP_OBJ) $(OBJ) Exception.o Render.o -o _SphSim.so $(LIB)

Exception.o : Exception.cpp
	$(CPP) $(OPT) $(GANDALF_INCLUDE) $(CFLAGS) -fexceptions -c $<

SphSimulation.o : SimulationIO.hpp SimulationIC.hpp SimAnalysis.hpp


force: ;

clean ::
	\rm -f *_wrap.cxx
	\rm -f */*_wrap.cxx
	\rm -f *.o
	\rm -f *.so
	\rm -f ../analysis/*.so
	\rm -f ../analysis/swig_generated/*.so
	\rm -f *.pyc
	\rm -f ../analysis/swig_generated/SphSim.py
	\rm -f *.d
	\rm -f *.optrpt
	\rm -f *.s
	\rm -f *.dbg
	\rm -f gandalf



ifeq ($(BUILD_DEPENDENCIES), 1)
# Handle the dependencies using black magic
get_dir_str = $(patsubst %/,%,$(dir $1))
WHICH_GPP   = which g++

ifeq ($(findstring g++,$(WHICH_GPP)),g++)
BUILD_DEPS = @set -e; g++ -MM -MG $(GANDALF_INCLUDE) $< |\
	 sed 's;\($*\)\.o[ :]*;$(call get_dir_str, $@)/\1.o $@ : ;g' > $@
else
BUILD_DEPS = @echo "Warning: Cannot build dependencies for $(notdir $<)"
endif

else
# Ignore dependencies
BUILD_DEPS =
endif

# Compute the list of .cpp files from the OBJ files
#    Note: chealpix is a .c file, so handle the special case
CPP_TGT := $(patsubst %.o,  %.cpp, $(OBJ))
CPP_TGT := $(patsubst chealpix.cpp, chealpix.c, $(CPP_TGT))
REL_DEPS := $(patsubst %.cpp,  %.d, $(CPP_TGT))
REL_DEPS := $(patsubst %.c,  %.d, $(REL_DEPS))



%.o : %.d

.PHONY : depends tidy

depends: $(REL_DEPS)
	@echo "Building Dependencies..."

%.d: %.cpp
	$(BUILD_DEPS)

%.d: %.c
	$(BUILD_DEPS)


-include $(REL_DEPS)
