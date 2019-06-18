#include <aligned_new>
//=============================================================================
//  SIMD.h
//  Contains macro definitions for SIMD lengths in all routines.
//
//  This file is part of GANDALF :
//  Graphical Astrophysics code for N-body Dynamics And Lagrangian Fluids
//  https://github.com/gandalfcode/gandalf
//  Contact : gandalfcode@gmail.com
//
//  Copyright (C) 2013  D. A. Hubber, G. Rosotti
//
//  GANDALF is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 2 of the License, or
//  (at your option) any later version.
//
//  GANDALF is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  General Public License (http://www.gnu.org/licenses) for more details.
//=============================================================================


#ifndef _SIMD_H_
#define _SIMD_H_

// Could this be done by templates?

// For FLOAT as double (the usual case in Gandalf), SIMD lengths are 4 (32 bytes
// AVX2) or 8 (64 bytes AVX512).

// AVX512
#define SIMD_NBYTES 64
// AVX2
//#define SIMD_NBYTES 32

#if defined(GANDALF_SINGLE_PRECISION)
#define SIMD_LENGTH SIMD_NBYTES/4
#elif defined(GANDALF_DOUBLE_PRECISION)
#define SIMD_LENGTH SIMD_NBYTES/8
#else
#define SIMD_LENGTH 0 // error
#endif

// MAX_NPART is the maximum likely Nleafmax.  From the Gandalf paper, 16 is
// probably enough, 32 is more than sufficient. SIMD_LENGTH must divide
// MAX_NPART, so MAX_NPART >= 16 (16 floats in an AVX512 vector).  MAX_NPART
// should be <= 32, so that boolean masks can be held in an int.  So the choices
// are 16 and 32.  SIMD_BOOL should match MAX_NPART.
#define MAX_NPART 16
#define SIMD_BOOL short
//#define MAX_NPART 32
//#define SIMD_BOOL int

// Remember to alignas(SIMD_NBYTES) the MAX_NPART arrays.

#endif
