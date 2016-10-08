//=================================================================================================
//  MfvNeighbourSearch.h
//  Header file containing class definitions for all Meshless Finite-Volume neighbour searching
//  data structures and algorithms.
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
//=================================================================================================


#ifndef _MFV_NEIGHBOUR_SEARCH_H_
#define _MFV_NEIGHBOUR_SEARCH_H_


#include <assert.h>
#include <iostream>
#include <string>
#include <vector>
#include "Precision.h"
#include "Constants.h"
#include "CodeTiming.h"
#include "Hydrodynamics.h"
#include "InlineFuncs.h"
#include "Nbody.h"
#include "NeighbourSearch.h"
#include "SmoothingKernel.h"
#include "Particle.h"
#include "MeshlessFV.h"
#include "DomainBox.h"
#include "Ewald.h"
#include "Parameters.h"
#include "KDTree.h"
#include "OctTree.h"
#if defined MPI_PARALLEL
#include "MpiExport.h"
#include "MpiNode.h"
#endif
using namespace std;




//=================================================================================================
//  Class MeshlessFVNeighbourSearch
/// \brief   MeshlessFVNeighbourSearch class definition.
/// \details MeshlessFVNeighbourSearch class definition.
/// \author  D. A. Hubber, G. Rosotti
/// \date    03/04/2013
//=================================================================================================
template <int ndim>
class MeshlessFVNeighbourSearch : public virtual NeighbourSearch<ndim>
{
#if defined MPI_PARALLEL
protected:
#endif
 public:


  //-----------------------------------------------------------------------------------------------
  MeshlessFVNeighbourSearch(FLOAT kernrangeaux, DomainBox<ndim> *boxaux,
                            SmoothingKernel<ndim> *kernaux, CodeTiming *timingaux) :
    NeighbourSearch<ndim>(kernrangeaux, boxaux, kernaux, timingaux) {};
  virtual ~MeshlessFVNeighbourSearch() {};


  //-----------------------------------------------------------------------------------------------
  virtual void UpdateAllProperties(int, int, MeshlessFVParticle<ndim> *, MeshlessFV<ndim> *,
                                   Nbody<ndim> *, DomainBox<ndim> &) = 0;
  virtual void UpdateGradientMatrices(int, int, MeshlessFVParticle<ndim> *, MeshlessFV<ndim> *,
                                      Nbody<ndim> *, DomainBox<ndim> &) = 0;
  virtual void UpdateGodunovFluxes(int, int, FLOAT, MeshlessFVParticle<ndim> *,
                                   MeshlessFV<ndim> *, Nbody<ndim> *, DomainBox<ndim> &) = 0;
  virtual void UpdateAllGravForces(int, int, MeshlessFVParticle<ndim> *, MeshlessFV<ndim> *,
                                   Nbody<ndim> *, DomainBox<ndim> &, Ewald<ndim> *) = 0;

};


//=================================================================================================
//  Class MeshlessFVTree
/// \brief   MeshlessFVTree class definition.
/// \details MeshlessFVTree class definition.
/// \author  D. A. Hubber
/// \date    21/04/2015
//=================================================================================================
template <int ndim, template<int> class ParticleType>
class MeshlessFVTree: public MeshlessFVNeighbourSearch<ndim>, public HydroTree<ndim,ParticleType>
{
#if defined MPI_PARALLEL
protected:
#endif
 public:

  using HydroTree<ndim,ParticleType>::activelistbuf;
  using HydroTree<ndim,ParticleType>::activepartbuf;
  using HydroTree<ndim,ParticleType>::allocated_buffer;
  using HydroTree<ndim,ParticleType>::box;
  using HydroTree<ndim,ParticleType>::cellbuf;
  using HydroTree<ndim,ParticleType>::gravity_mac;
  using HydroTree<ndim,ParticleType>::kernp;
  using HydroTree<ndim,ParticleType>::kernrange;
  using HydroTree<ndim,ParticleType>::kernrangesqd;
  using HydroTree<ndim,ParticleType>::levelneibbuf;
  using HydroTree<ndim,ParticleType>::multipole;
  using HydroTree<ndim,ParticleType>::neibcheck;
  using HydroTree<ndim,ParticleType>::neibpartbuf;
  using HydroTree<ndim,ParticleType>::Ngravcellmaxbuf;
  using HydroTree<ndim,ParticleType>::Nleafmax;
  using HydroTree<ndim,ParticleType>::Nneibmaxbuf;
  using HydroTree<ndim,ParticleType>::Ntot;
  using HydroTree<ndim,ParticleType>::Ntotmax;
  using HydroTree<ndim,ParticleType>::Ntotmaxold;
  using HydroTree<ndim,ParticleType>::Ntotold;
  using HydroTree<ndim,ParticleType>::timing;
  using HydroTree<ndim,ParticleType>::tree;
  using HydroTree<ndim,ParticleType>::ghosttree;
#ifdef MPI_PARALLEL
  using HydroTree<ndim,ParticleType>::mpighosttree;
  using HydroTree<ndim,ParticleType>::Nmpi;
  using HydroTree<ndim,ParticleType>::prunedtree;
  using HydroTree<ndim,ParticleType>::sendprunedtree;
#endif


  //-----------------------------------------------------------------------------------------------
  MeshlessFVTree(string tree_type,
                 int _Nleafmax, int _Nmpi, int _pruning_level_min, int _pruning_level_max,
                 FLOAT _thetamaxsqd, FLOAT _kernrange, FLOAT _macerror,
                 string _gravity_mac, string _multipole, DomainBox<ndim> *_box,
                 SmoothingKernel<ndim> *_kern, CodeTiming *_timing, ParticleTypeRegister&); //:
    /*NeighbourSearch<ndim>(_kernrange, _box, _kern, _timing),
    MeshlessFVNeighbourSearch<ndim>(_kernrange, _box, _kern, _timing),
    HydroTree<ndim,ParticleType>(_Nleafmax, _Nmpi, _thetamaxsqd, _kernrange, _macerror,
                                          _gravity_mac, _multipole, _box, _kern, _timing) {};*/
  virtual ~MeshlessFVTree(); //{};


  //-----------------------------------------------------------------------------------------------
  virtual void UpdateAllProperties(int, int, MeshlessFVParticle<ndim> *, MeshlessFV<ndim> *,
                                   Nbody<ndim> *, DomainBox<ndim> &);
  virtual void UpdateGradientMatrices(int, int, MeshlessFVParticle<ndim> *, MeshlessFV<ndim> *,
                                      Nbody<ndim> *, DomainBox<ndim> &);
  virtual void UpdateGodunovFluxes(int, int, FLOAT, MeshlessFVParticle<ndim> *,
                                   MeshlessFV<ndim> *, Nbody<ndim> *, DomainBox<ndim> &);
  virtual void UpdateAllGravForces(int, int, MeshlessFVParticle<ndim> *, MeshlessFV<ndim> *,
                                   Nbody<ndim> *, DomainBox<ndim> &, Ewald<ndim> *);

};


#endif
