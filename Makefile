#==================================================================================================
#  GANDALF v0.4.0 Makefile frontend
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

MPICPP             = mpiicpc
CPP                = icpc
PYTHON             = python
COMPILER_MODE      = FAST
PRECISION          = DOUBLE
# Switch off OpenMP while doing serial testing using Advisor and VTune
# - otherwise the analyses are confusing.  May need to switch off
# inlining as well for Advisor to see where things really are, but
# that will certainly prevent vectorisation.
#OPENMP             = 0
OPENMP             = 1
PYSNAP_PRECISION   = DOUBLE
OUTPUT_LEVEL       = 1
DEBUG_LEVEL        = 0
BUILD_DEPENDENCIES = 1

# FFTW libary flags and paths.  If paths are empty, tries standard
# default linux paths.  On DIaL the paths are set in CPATH and LIBRARY_PATH.
#--------------------------------------------------------------------------------------------------
FFTW               = 1
FFTW_INCLUDE       =
FFTW_LIBRARY       =


# GNU Scientific library flags and paths.  If paths are empty, tries standard default linux paths.  On DIaL the paths are set in CPATH and LIBRARY_PATH.
#--------------------------------------------------------------------------------------------------
GSL                = 1
GSL_INCLUDE        =
GSL_LIBRARY        =


# MKL flags and paths.  Use the link line advisor
# https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
# to set these.
# --------------------------------------------------------------------------------------------------
MKL                = 1
MKL_INCLUDE        = -I${MKLROOT}/include
MKL_LIBRARY        = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl

# Select location of python and numpy libraries.  If blank, make will try to
# find the location of the libraries automatically using installed python
# utilities.  If you have multiple versions of python installed on your
# computer, then select the prefered version with the PYTHON variable above.
#--------------------------------------------------------------------------------------------------
PYLIB =
NUMPY =
GTEST = $(GTEST_DIR)




# Don't delete this command! Makes sure that defined variables here are passed
# to the others make
export

gandalf:
	@+$(MAKE) -C src

all:
	@+$(MAKE) -C src

executable:
	@+$(MAKE) executable -C src

unittests:
	@+$(MAKE) unittests -C src

clean:
	@+$(MAKE) clean -C src

depends:
	@+$(MAKE) depends -C src

installpython:
	./setup.py install

install:
	cp bin/gandalf /usr/local/bin
