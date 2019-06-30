#include <aligned_new>
//=================================================================================================
//  GradhSph.cpp
//  Contains all functions for calculating conservative 'grad-h' SPH quantities
//  (See Springel & Hernquist (2002) and Price & Monaghan (2007)).
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


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <math.h>
#include <algorithm>
#include "SIMD.h"
#include "Precision.h"
#include "Sph.h"
#include "Particle.h"
#include "Parameters.h"
#include "SmoothingKernel.h"
#include "EOS.h"
#include "Debug.h"
#include "Exception.h"
#include "InlineFuncs.h"
#include "NeighbourManager.h"
using namespace std;



//=================================================================================================
//  GradhSph::GradhSph
/// GradhSph class constructor.  Calls main SPH class constructor and also
/// sets additional kernel-related quantities
//=================================================================================================
template <int ndim, template<int> class kernelclass>
GradhSph<ndim, kernelclass>::GradhSph(int hydro_forces_aux, int self_gravity_aux,
  FLOAT alpha_visc_aux, FLOAT beta_visc_aux, FLOAT h_fac_aux, FLOAT h_converge_aux,
  aviscenum avisc_aux, acondenum acond_aux, tdaviscenum tdavisc_aux,
  string gas_eos_aux, string KernelName, SimUnits &units, Parameters *params):
  GradhSphBase<ndim>(hydro_forces_aux, self_gravity_aux, alpha_visc_aux, beta_visc_aux,
            h_fac_aux, h_converge_aux, avisc_aux, acond_aux, tdavisc_aux,
            gas_eos_aux, KernelName, units, params),
  kern(kernelclass<ndim>(KernelName))
{
  this->kernp      = &kern;
  this->kernrange  = this->kernp->kernrange;
}



//=================================================================================================
//  GradhSph::~GradhSph
/// GradhSph class destructor
//=================================================================================================
template <int ndim, template<int> class kernelclass>
GradhSph<ndim, kernelclass>::~GradhSph()
{
  DeallocateMemory();
}



//=================================================================================================
//  GradhSph::AllocateMemory
/// Allocate main SPH particle array
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void GradhSph<ndim, kernelclass>::AllocateMemory(int N)
{
  debug2("[GradhSph::AllocateMemory]");

  if (N > Nhydromax || !allocated) {

    GradhSphParticle<ndim>* newsphdata =
        new struct GradhSphParticle<ndim>[N];

    // Swap so that sphdata points to the new memory
    std::swap(newsphdata, sphdata) ;
    if (allocated) {
      // Copy back particle data
      std::copy(newsphdata,newsphdata+Nhydromax,sphdata);
      delete[] newsphdata;
    }

    Nhydromax=N;
    allocated        = true;
    hydrodata_unsafe = sphdata;
    sphdata_unsafe   = sphdata;

  }

  assert(Nhydromax >= Nhydro);
  assert(sphdata);

  return;
}



//=================================================================================================
//  GradhSph::DeallocateMemory
/// Deallocate main array containing SPH particle data.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void GradhSph<ndim, kernelclass>::DeallocateMemory(void)
{
  debug2("[GradhSph::DeallocateMemory]");

  if (allocated) {
    delete[] sphdata;
  }
  allocated = false;

  return;
}


//=================================================================================================
//  GradhSph::ComputeHNonArray
/// Compute the value of the smoothing length of particle 'i' by iterating the relation :
/// h = h_fac*(m/rho)^(1/ndim).
/// Uses the previous value of h as a starting guess and then uses either a Newton-Rhapson solver,
/// or fixed-point iteration, to converge on the correct value of h.  The maximum tolerance used
/// for deciding whether the iteration has converged is given by the 'h_converge' parameter.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
int GradhSph<ndim, kernelclass>::ComputeHNonArray
 (SphParticle<ndim> &part,                                ///< [inout] Particle i data
  FLOAT hmax,                                             ///< [in] Maximum smoothing length
  const NeighbourList<DensityParticle> &ngbs,             ///< [in] Neighbour properties
  Nbody<ndim> *nbody)                                     ///< [in] Main N-body object
{
  int j;                               // Neighbour id
  int k;                               // Dimension counter
  int iteration = 0;                   // h-rho iteration counter
  int iteration_max = 30;              // Max. no of iterations
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT h_lower_bound = (FLOAT) 0.0;   // Lower bound on h
  FLOAT h_upper_bound = hmax;          // Upper bound on h
  FLOAT invh;                          // 1 / h
  FLOAT invhsqd;                       // (1 / h)^2
  FLOAT ssqd;                          // Kernel parameter squared, (r/h)^2
  GradhSphParticle<ndim>& parti = static_cast<GradhSphParticle<ndim>& > (part);


  // If there are sink particles present, check if the particle is inside one.
  // If so, then adjust the iteration bounds and ensure they are valid (i.e. hmax is large enough)
  if (sink_particles) {
    // This makes hmax far too large.
    // h_lower_bound = h_fac*pow(parti.m/this->rho_sink, Sph<ndim>::invndim);
    // The not-OK return value is 0.
    // if (hmax < h_lower_bound) return -1;
    // The above two lines mean that with sink particles, ComputeH is not
    // executed, and h is always too big (the initial value).
    //
    // This looks correct.  The check can be also done earlier, in the calling
    // routine.
    if (parti.flags.check(inside_sink)) {
      h_lower_bound = hmin_sink;
      if (hmax < h_lower_bound) return 0;
    }
  }

  int Nneib = ngbs.size();

  // Some basic sanity-checking in case of invalid input into routine
  assert(Nneib > 0);
  assert(hmax > (FLOAT) 0.0);
  assert(!parti.flags.is_dead());
  assert(parti.m > (FLOAT) 0.0);



  // Main smoothing length iteration loop
  //===============================================================================================
  do {

    // Initialise all variables for this value of h
    iteration++;
    invh           = (FLOAT) 1.0/parti.h;
    parti.rho      = (FLOAT) 0.0;
    parti.invomega = (FLOAT) 0.0;
    parti.zeta     = (FLOAT) 0.0;
    parti.hfactor  = pow(invh,ndim);
    invhsqd        = invh*invh;

    // Loop over all nearest neighbours in list to calculate density, omega and zeta.
    //---------------------------------------------------------------------------------------------
    for (j=0; j<Nneib; j++) {
      const DensityParticle &ngb = ngbs[j];
      for (k=0; k<ndim; k++) dr[k] = ngb.r[k] - parti.r[k];
      ssqd           = invhsqd*DotProduct(dr, dr, ndim);
      parti.rho      += ngb.m*kern.w0_s2(ssqd);
      parti.invomega += ngb.m*invh*kern.womega_s2(ssqd);
      parti.zeta     += ngb.m*kern.wzeta_s2(ssqd);
    }
    //---------------------------------------------------------------------------------------------

    parti.rho      *= parti.hfactor;
    parti.invomega *= parti.hfactor;
    parti.zeta     *= invhsqd;

    // Density must at least equal its self-contribution
    // (failure could indicate neighbour list problem)
    assert(parti.rho >= 0.99*parti.m*parti.hfactor*kern.w0_s2(0.0));

    // If h changes below some fixed tolerance, exit iteration loop
    if (parti.rho > (FLOAT) 0.0 && parti.h > h_lower_bound &&
        fabs(parti.h - h_rho_func(parti.m, parti.rho))*invh < h_converge) break;

    // Use fixed-point iteration, i.e. h_new = h_fac*(m/rho_old)^(1/ndim), for now.  If this does
    // not converge in a reasonable number of iterations (iteration_max), then assume something is
    // wrong and switch to a bisection method, which should be guaranteed to converge, albeit much
    // more slowly.  (N.B. will implement Newton-Raphson soon)
    //---------------------------------------------------------------------------------------------
    if (iteration < iteration_max) {
      parti.h = h_rho_func(parti.m, parti.rho);
    }
    else if (iteration == iteration_max) {
      parti.h = (FLOAT) 0.5*(h_lower_bound + h_upper_bound);
    }
    else if (iteration < 5*iteration_max) {
      if (parti.rho < small_number || parti.h > h_rho_func(parti.m, parti.rho)) {
        //(parti.rho + rho_hmin)*pow(parti.h,ndim) > pow(h_fac,ndim)*parti.m) {
        h_upper_bound = parti.h;
      }
      else {
        h_lower_bound = parti.h;
      }
      parti.h = (FLOAT) 0.5*(h_lower_bound + h_upper_bound);
    }
    else {
      cout << "H ITERATION : " << iteration << "    h : " << parti.h
           << "   rho : " << parti.rho << "   h_upper " << h_upper_bound << "    hmax :  " << hmax
           << "   h_lower : " << h_lower_bound << "    " << parti.hfactor << "    m : " << parti.m
           << "     " << parti.m*parti.hfactor*kern.w0(0.0) << "    " << Nneib << endl;

      string message = "Problem with convergence of h-rho iteration";
      ExceptionHandler::getIstance().raise(message);
    }

    // If the smoothing length is too large for the neighbour list, exit routine and flag neighbour
    // list error in order to generate a larger neighbour list (not properly implemented yet).
    if (parti.h > hmax) return 0;

  } while (parti.h > h_lower_bound && parti.h < h_upper_bound);
  //===============================================================================================


  // Normalise all SPH sums correctly
  parti.h         = max(h_rho_func(parti.m, parti.rho), h_lower_bound);
  invh            = 1/parti.h;
  parti.hfactor   = pow(invh, ndim+1);
  parti.hrangesqd = kern.kernrangesqd*parti.h*parti.h;
  parti.div_v     = (FLOAT) 0.0;
  assert(isnormal(part.h));
  assert(part.h >= h_lower_bound);

  // Calculate the minimum neighbour potential (used later to identify new sinks)
  if (create_sinks == 1) {
    parti.flags.set(potmin);
    for (j=0; j<Nneib; j++) {
      const DensityParticle &ngb = ngbs[j];
      // Surely these two lines should be swapped?
      FLOAT drsqd = DotProduct(dr,dr,ndim);
      for (k=0; k<ndim; k++) dr[k] = ngb.r[k] - parti.r[k];
      if (ngb.gpot > (FLOAT) 1.000000001*parti.gpot &&
          drsqd*invhsqd < kern.kernrangesqd) parti.flags.unset(potmin);
    }
  }

  // If there are star particles, compute N-body zeta correction term
  //-----------------------------------------------------------------------------------------------
  if (!part.flags.check(inside_sink)) {
    parti.invomega = (FLOAT) 1.0 - h_rho_deriv(parti.h, part.m, parti.rho)*parti.invomega;
    parti.invomega = (FLOAT) 1.0/parti.invomega;

    // Use Hubber et al. (2013) zeta SPH-star conservative-gravity correction terms
    if (this->conservative_sph_star_gravity == 1) {
      if (nbody->nbody_softening == 1) {
        for (j=0; j<nbody->Nstar; j++) {
          invhsqd = pow((FLOAT) 2.0 / (parti.h + nbody->stardata[j].h),2);
          for (k=0; k<ndim; k++) dr[k] = nbody->stardata[j].r[k] - parti.r[k];
          ssqd = DotProduct(dr,dr,ndim)*invhsqd;
          parti.zeta += nbody->stardata[j].m*invhsqd*kern.wzeta_s2(ssqd);
        }
      }
      else {
        invhsqd = (FLOAT) 4.0*invh*invh;
        for (j=0; j<nbody->Nstar; j++) {
          for (k=0; k<ndim; k++) dr[k] = nbody->stardata[j].r[k] - parti.r[k];
          ssqd = DotProduct(dr,dr,ndim)*invhsqd;
          parti.zeta += nbody->stardata[j].m*invhsqd*kern.wzeta_s2(ssqd);
        }
      }
    }
    parti.zeta = h_rho_deriv(parti.h, part.m, parti.rho)*parti.zeta*parti.invomega;

  }
  else {
    parti.invomega = (FLOAT) 1.0;
    parti.zeta     = (FLOAT) 0.0;
  }


  // Set important thermal variables here
  ComputeThermalProperties(parti);

  if (tdavisc == cd2010) {
    this->ComputeCullenAndDehnenViscosity(parti, ngbs, kern);
  }

  // If h is invalid (i.e. larger than maximum h), then return error code (0)
  if (parti.h <= hmax) return 1;
  // The not-OK return value is zero.
  // else return -1;
  else return 0;
}



//=================================================================================================
//  GradhSph::ComputeH
/// Compute the value of the smoothing length of particle 'i' by iterating the relation :
/// h = h_fac*(m/rho)^(1/ndim).
/// Uses the previous value of h as a starting guess and then uses either a Newton-Rhapson solver,
/// or fixed-point iteration, to converge on the correct value of h.  The maximum tolerance used
/// for deciding whether the iteration has converged is given by the 'h_converge' parameter.
//=================================================================================================

template <int ndim, template<int> class kernelclass>
int GradhSph<ndim, kernelclass>::ComputeH
 (SphParticle<ndim>* part,                                ///< [inout] Particle i data array
  const int npart,                                        ///< [in] Particle i data array length
  const FLOAT hrangesqd,
  const FLOAT hmax,                                       ///< [in] Maximum smoothing length
  // No const, otherwise we have to have const-qualified copies of several
  // NeighbourManager methods - which is really what we should do.
  NeighbourManager<ndim,DensityParticle>& neibmanager,    ///< [in] NeighbourManager object
  Nbody<ndim> *nbody)                                     ///< [in] Main N-body object
{
  int j;                               // Neighbour id
  int k;                               // Dimension counter
  int iteration[MAX_NPART];            // h-rho iteration counter
  int iteration_max = 30;              // Max. no of iterations
  FLOAT h_lower_bound[MAX_NPART];      // Lower bound on h
  FLOAT h_upper_bound[MAX_NPART];      // Upper bound on h
  FLOAT h[MAX_NPART];                  // h
  FLOAT invh[MAX_NPART];               // 1 / h
  FLOAT invhsqd[MAX_NPART];            // (1 / h)^2
  FLOAT ssqd[MAX_NPART];               // Kernel parameter squared, (r/h)^2
  FLOAT rho[MAX_NPART];                //
  FLOAT invomega[MAX_NPART];           //
  FLOAT zeta[MAX_NPART];               //
  FLOAT hfactor[MAX_NPART];            //
  FLOAT m[MAX_NPART];                  //
  FLOAT div_v[MAX_NPART];              //
  FLOAT gpot[MAX_NPART];               //
  bool _potmin[MAX_NPART];             //
  bool _inside_sink[MAX_NPART];        //

  bool any, all;

  int level[MAX_NPART];
  int iorig[MAX_NPART];
  FLOAT r[ndim][MAX_NPART];
  FLOAT dr[ndim][MAX_NPART];           // Relative position vector
  FLOAT drsqd[MAX_NPART];
  int neib;
  // Is it useful to have these as arrays of the size of the number of
  // neighbours and as arguments, so that the forces routines (and Cullen and
  // Dehnen viscosity ) don't each call MaskParticleNeib...?  No, because
  // MaskParticleNeib does not make aperisodic correction ot the neighbour r.
  bool culled[MAX_NPART];
  bool smoothed_grav[MAX_NPART];
  bool direct[MAX_NPART];

  bool converged[MAX_NPART];
  bool iterating[MAX_NPART];
  bool any_iterating;
  bool mask[MAX_NPART];

  assert(dynamic_cast<GradhSphParticle<ndim>*>(part));
  // This should be reinterpret_cast.
  // GradhSphParticle<ndim>* parti = static_cast<GradhSphParticle<ndim>*>(part);
  GradhSphParticle<ndim>* parti = reinterpret_cast<GradhSphParticle<ndim>*>(part);

  for (int p=0; p<npart; p++) h_lower_bound[p] = (FLOAT) 0.0;
  for (int p=0; p<npart; p++) _inside_sink[p] = parti[p].flags.check(inside_sink);

  // If there are sink particles present, check if a particle is inside one.  If
  // so, then adjust the iteration bounds and ensure they are valid (i.e. hmax
  // is large enough).  If hmax is too small, then no particle is processed:
  // hmax is the cell hmax, and if it is increased, then the tree has to be
  // traversed again.  Thus, there can't be a per-particle hmax.
  if (sink_particles) {
    for (int p=0; p<npart; p++) if (_inside_sink[p]) h_lower_bound[p] = hmin_sink;
    any = false;
    for (int p=0; p<npart; p++) any = any || hmax < h_lower_bound[p];
    if (any) return 0;
  }

  // MaskParticleNeib... works with an array of hrangesqd, corresponding to the
  // active particles.  For the h computation, the cell hrangesqd is used.
  FLOAT _hrangesqd[MAX_NPART];
  for (int p=0; p<npart; p++) _hrangesqd[p] = hrangesqd;

  Typemask densmask[MAX_NPART];
  for (int p=0; p<npart; p++) densmask[p] = types[parti[p].ptype].hmask;

  for (int p=0; p<npart; p++) h_upper_bound[p] = hmax;

  int Nneib = neibmanager.GetNumNeib();

  // Some basic sanity-checking in case of invalid input into routine
  assert(Nneib > 0);
  assert(hmax > (FLOAT) 0.0);
  for (int p=0; p<npart; p++) {
    assert(!parti[p].flags.is_dead());
    assert(parti[p].m > (FLOAT) 0.0);
  }

  for (int p=0; p<npart; p++) level[p] = parti[p].level;
  for (int p=0; p<npart; p++) iorig[p] = parti[p].iorig;
  for (int k=0; k<ndim; k++)
    for (int p=0; p<npart; p++) r[k][p] = parti[p].r[k];
  for (int p=0; p<npart; p++) m[p] = parti[p].m;
  // Set to true below for (int p=0; p<npart; p++) _potmin[p] = parti[p].flags.check(potmin);;

  // Shouldn't there be one convergence condition?
  for (int p=0; p<npart; p++) converged[p] = false;
  for (int p=0; p<npart; p++) iterating[p] = true;
  for (int p=0; p<npart; p++) iteration[p] = 0;
  for (int p=0; p<npart; p++) h[p] = parti[p].h;
  // Main smoothing length iteration loop
  //===============================================================================================
  do {
    // Everything is masked with iterating[p], but that may not make any
    // difference to the speed once SIMDed, unless a whole SIMD length is masked
    // out, but will make a difference to the results.

    // Initialise all variables for this value of h
    for (int p=0; p<npart; p++) if (iterating[p]) iteration[p]++;
    for (int p=0; p<npart; p++) if (iterating[p]) invh[p]     = (FLOAT) 1.0 / h[p];
    for (int p=0; p<npart; p++) if (iterating[p]) invhsqd[p]  = invh[p]*invh[p];
    for (int p=0; p<npart; p++) if (iterating[p]) rho[p]      = (FLOAT) 0.0;
    for (int p=0; p<npart; p++) if (iterating[p]) invomega[p] = (FLOAT) 0.0;
    for (int p=0; p<npart; p++) if (iterating[p]) zeta[p]     = (FLOAT) 0.0;
    for (int p=0; p<npart; p++) if (iterating[p]) hfactor[p]  = pow(invh[p],ndim);

    // Loop over all nearest neighbours in list to calculate density, omega and zeta.
    //---------------------------------------------------------------------------------------------
    for (int neib=0; neib<Nneib; neib++) {

      neibmanager.MaskParticleNeibGather(level, iorig, r, npart, neib, densmask, _hrangesqd, iterating,
					 dr, drsqd, culled, smoothed_grav, direct);

      // the following is only for culled
      const DensityParticle& ngb = neibmanager.GetNeib(neib); // Already in
							      // cache, from
							      // MaskParticleNeibGather
      for (int p=0; p<npart; p++) if (iterating[p] && culled[p]) ssqd[p]      = invhsqd[p]*drsqd[p];
      // The kernel functions need to be SIMDed (TabulatedKernel is easy, using
      // a gather; the others less so).  And the kernels are 15% of the total
      // ComputeHydroForces and ComputeH time.
      for (int p=0; p<npart; p++) if (iterating[p] && culled[p]) rho[p]      += ngb.m*kern.w0_s2(ssqd[p]);
      for (int p=0; p<npart; p++) if (iterating[p] && culled[p]) invomega[p] += ngb.m*invh[p]*kern.womega_s2(ssqd[p]);
      for (int p=0; p<npart; p++) if (iterating[p] && culled[p]) zeta[p]     += ngb.m*kern.wzeta_s2(ssqd[p]);
    }
    //---------------------------------------------------------------------------------------------

    for (int p=0; p<npart; p++) if (iterating[p]) rho[p]      *= hfactor[p];
    for (int p=0; p<npart; p++) if (iterating[p]) invomega[p] *= hfactor[p];
    for (int p=0; p<npart; p++) if (iterating[p]) zeta[p]     *= invhsqd[p];

    // Density must at least equal its self-contribution
    // (failure could indicate neighbour list problem)
    for (int p=0; p<npart; p++) if (iterating[p]) assert(rho[p] >= 0.99*m[p]*hfactor[p]*kern.w0_s2(0.0));

    // If h changes below some fixed tolerance for all particles, exit iteration loop
    // There is a SIMD h_rho_func, but here we are out of the neighbour loop, so
    // performance is not important(?).
    for (int p=0; p<npart; p++)
      if (iterating[p]) converged[p]
			  = rho[p] > (FLOAT) 0.0
			  && h[p] > h_lower_bound[p]
			  && fabs(h[p] - h_rho_func(m[p], rho[p]))*invh[p] < h_converge;

    // This is to match the scalar case, but we could test for no particle still
    // iterating.
    all = true;
    for (int p=0; p<npart; p++) all = all && converged[p];
    if (all) break;

    for (int p=0; p<npart; p++) if (iterating[p]) iterating[p] = !converged[p];

    // Use fixed-point iteration, i.e. h_new = h_fac*(m/rho_old)^(1/ndim), for now.  If this does
    // not converge in a reasonable number of iterations (iteration_max), then assume something is
    // wrong and switch to a bisection method, which should be guaranteed to converge, albeit much
    // more slowly.  (N.B. will implement Newton-Raphson soon)
    //---------------------------------------------------------------------------------------------
    for (int p=0; p<npart; p++) {
      if (iterating[p] && iteration[p] < iteration_max) {
	h[p] = h_rho_func(m[p], rho[p]);
      }
      else if (iterating[p] && iteration[p] == iteration_max) {
	h[p] = (FLOAT) 0.5*(h_lower_bound[p] + h_upper_bound[p]);
      }
      else if (iterating[p] && iteration[p] < 5*iteration_max) {
	if (rho[p] < small_number || h[p] > h_rho_func(m[p], rho[p])) {
	  //(rho[p] + rho_hmin)*pow(h[p],ndim) > pow(h_fac,ndim)*m[p]) {
	  h_upper_bound[p] = h[p];
	}
	else {
	  h_lower_bound[p] = h[p];
	}
	h[p] = (FLOAT) 0.5*(h_lower_bound[p] + h_upper_bound[p]);
      }
      else if (iterating[p]) {
	cout << "p : " << p << "   H ITERATION : " << iteration[p] << "   h[p] : " << h[p]
	     << "   rho[p] : " << rho[p] << "   h_upper[p] " << h_upper_bound[p] << "   hmax :  " << hmax
	     << "   h_lower : " << h_lower_bound[p] << "   hfactor :  " << hfactor[p] << "   m : " << m[p]
	     << "     " << m[p]*hfactor[p]*kern.w0(0.0) << "    " << Nneib << endl;
	
	string message = "Problem with convergence of h-rho iteration";
	ExceptionHandler::getIstance().raise(message);
      }
    }

    // If the smoothing length is too large for the neighbour list, exit routine
    // and flag neighbour list error in order to generate a larger neighbour
    // list (not properly implemented yet).
    any = false;
    for (int p=0; p<npart; p++) any = any || h[p] > hmax;
    if (any) {
      // To match the scalar version, copy changed SIMDs back to parti.  The
      // return value of zero should indicate to any calling function that the
      // active particle array is invalid.  The actual calling function
      // GradhSphTree does recreate the active particle working array again from
      // the main particle array.  So this code is not necessary except to match
      // the serial case exactly.
      for (int p=0; p<npart; p++) parti[p].h        = h[p];
      for (int p=0; p<npart; p++) parti[p].hfactor  = hfactor[p];
      for (int p=0; p<npart; p++) parti[p].rho      = rho[p];
      for (int p=0; p<npart; p++) parti[p].invomega = invomega[p];
      for (int p=0; p<npart; p++) parti[p].zeta     = zeta[p];

      return 0;
    }

    for (int p=0; p<npart; p++)
      if (iterating[p]) iterating[p] = h[p] > h_lower_bound[p] && h[p] < h_upper_bound[p];
    any = false;
    for (int p=0; p<npart; p++) any = any || iterating[p];
    any_iterating = any;

  } while (any_iterating);
  //===============================================================================================


  // Normalise all SPH sums correctly
  for (int p=0; p<npart; p++) h[p]         = max(h_rho_func(m[p], rho[p]), h_lower_bound[p]);
  for (int p=0; p<npart; p++) invh[p]      = (FLOAT) 1.0 / h[p]; // 1/h[p];
  // MJF invhsqd was not recalculated in the original code, why not?
  //for (int p=0; p<npart; p++) invhsqd[p]   = invh[p]*invh[p]; 

  // The particle hrangesqd is not set here, but only when copied back to parti,
  // because the hrangesqd used for MaskParticleNeigbourGather below is the cell
  // hrangesqsd broadcast to an array.

  // ndim+1 because this is being set for the kernel derivative functions used in
  // the forces calculation; this is should be calculated in those.
  for (int p=0; p<npart; p++) hfactor[p]   = pow(invh[p], ndim+1);
  for (int p=0; p<npart; p++) div_v[p]     = (FLOAT) 0.0;
  for (int p=0; p<npart; p++) assert(isnormal(h[p]));
  for (int p=0; p<npart; p++) assert(h[p] >= h_lower_bound[p]);

  // Calculate the minimum neighbour potential (used later to identify new sinks)
  if (create_sinks == 1) {
    for (int p=0; p<npart; p++) gpot[p] = parti[p].gpot;
    for (int p=0; p<npart; p++) _potmin[p] = true;
    for (int p=0; p<npart; p++) mask[p] = true;
    for (int neib=0; neib<Nneib; neib++) {
      neibmanager.MaskParticleNeibGather(level, iorig, r, npart, neib, densmask, _hrangesqd, mask,
					 dr, drsqd, culled, smoothed_grav, direct);

      // the following is only for culled
      const DensityParticle& ngb = neibmanager.GetNeib(neib);
      for (int p=0; p<npart; p++)
	if (mask[p] && culled[p]) _potmin[p] = !(ngb.gpot > (FLOAT) 1.000000001 * gpot[p]
				      && drsqd[p]*invhsqd[p] < kern.kernrangesqd
				      );
      // If further neighbours will make no difference to any of the active
      // particles, then exit the neighbour loop.
      any = false;
      for (int p=0; p<npart; p++) any = any || _potmin[p];
      if (!any) break;
    }
  }

  // If there are star particles, compute N-body zeta correction term
  //-----------------------------------------------------------------------------------------------
  for (int p=0; p<npart; p++) mask[p] = !_inside_sink[p];
  any = false;
  for (int p=0; p<npart; p++) any = any || mask[p];
  if (any) { // Some outside a sink
    for (int p=0; p<npart; p++)
      if (mask[p]) invomega[p] = (FLOAT) 1.0 / ((FLOAT) 1.0 - h_rho_deriv(h[p], m[p], rho[p])*invomega[p]);
      else         invomega[p] = (FLOAT) 1.0;

    // Use Hubber et al. (2013) zeta SPH-star conservative-gravity correction terms
    if (this->conservative_sph_star_gravity == 1) {
      if (nbody->nbody_softening == 1) {
        for (j=0; j<nbody->Nstar; j++) {
          for (int p=0; p<npart; p++) if (mask[p]) invhsqd[p] = pow((FLOAT) 2.0 / (h[p] + nbody->stardata[j].h),2);
          for (k=0; k<ndim; k++) for (int p=0; p<npart; p++) if (mask[p]) dr[k][p] = nbody->stardata[j].r[k] - r[k][p];
	  DotProduct<ndim>(dr, dr, mask, npart, drsqd);
	  for (int p=0; p<npart; p++) if (mask[p]) ssqd[p] = drsqd[p]*invhsqd[p];
          for (int p=0; p<npart; p++) if (mask[p]) zeta[p] += nbody->stardata[j].m*invhsqd[p]*kern.wzeta_s2(ssqd[p]);
        }
      }
      else {
        for (int p=0; p<npart; p++) if (mask[p]) invhsqd[p] = (FLOAT) 4.0 * invh[p]*invh[p];
        for (j=0; j<nbody->Nstar; j++) {
          for (k=0; k<ndim; k++) for (int p=0; p<npart; p++) if (mask[p]) dr[k][p] = nbody->stardata[j].r[k] - r[k][p];
	  DotProduct<ndim>(dr, dr, mask, npart, drsqd);
	  for (int p=0; p<npart; p++) if (mask[p]) ssqd[p] = drsqd[p]*invhsqd[p];
          for (int p=0; p<npart; p++) if (mask[p]) zeta[p] += nbody->stardata[j].m*invhsqd[p]*kern.wzeta_s2(ssqd[p]);
        }
      }
    }
    for (int p=0; p<npart; p++)
      if (mask[p]) zeta[p] = h_rho_deriv(h[p], m[p], rho[p])*zeta[p]*invomega[p];
      else         zeta[p] = (FLOAT) 0.0;

  }
  else { // None outside a sink
    for (int p=0; p<npart; p++) invomega[p] = (FLOAT) 1.0;
    for (int p=0; p<npart; p++) zeta[p]     = (FLOAT) 0.0;
  }


  // copy SIMDs back to parti
  for (int p=0; p<npart; p++) parti[p].h = h[p]; //
  for (int p=0; p<npart; p++) parti[p].hfactor  = hfactor[p]; //
  // The particle hrangesqd
  for (int p=0; p<npart; p++) parti[p].hrangesqd = kern.kernrangesqd*h[p]*h[p];
  for (int p=0; p<npart; p++) parti[p].div_v = div_v[p]; //
  for (int p=0; p<npart; p++) parti[p].rho = rho[p]; //
  for (int p=0; p<npart; p++) parti[p].invomega = invomega[p]; //
  for (int p=0; p<npart; p++) parti[p].zeta     = zeta[p]; //
  if (create_sinks == 1)
    for (int p=0; p<npart; p++) parti[p].flags.set(potmin, _potmin[p]);

  // Set important thermal variables here.  This will take a lot of changes to
  // make SIMD.
  for (int p=0; p<npart; p++) ComputeThermalProperties(parti[p]);

  // Cullen and Dehnen viscosity depends on sound speed, so has to come after
  // ComputeThermalProperties.  Cullen and Dehnen viscosity is not yet fully
  // SIMD.
  //
  // It would be safer to copy the information back from parti to the arrays.
  if (tdavisc == cd2010) {
    // Buffers for gradients of vectors
    FLOAT dv[MAX_NPART][ndim][ndim];
    FLOAT da[MAX_NPART][ndim][ndim];
    FLOAT rr[MAX_NPART][ndim][ndim];

    FLOAT v[ndim][MAX_NPART];
    FLOAT a[ndim][MAX_NPART];
    FLOAT sound[MAX_NPART];

    FLOAT hfac[MAX_NPART];
    FLOAT w[MAX_NPART];
    FLOAT alpha[MAX_NPART];
    FLOAT dalphadt[MAX_NPART];

    // These are from the current values:  h[p], invh[p], rho[p], hfactor[p],
    // level[p], iorig[p], r[k][p], m[p], densmask[p], _hrangesqd[p].
    for (int k=0; k<ndim; k++)
      for (int p=0; p<npart; p++) v[k][p]  = parti[p].v[k];
    for (int k=0; k<ndim; k++)
      for (int p=0; p<npart; p++) a[k][p]  = parti[p].a[k];
    for (int p=0; p<npart; p++) alpha[p]   = parti[p].alpha;
    // This has to be from parti (until ComputeThermalProperties is made SIMD to
    // use the current values).
    for (int p=0; p<npart; p++) sound[p]   = parti[p].sound;

    for (int p=0; p<npart; p++) hfac[p] = invh[p] * hfactor[p] / rho[p];

    for (int p=0; p<npart; p++)
      for (int i=0; i<ndim; ++i)
	for (int j=0; j<ndim; ++j)
	  rr[p][i][j] = da[p][i][j] = dv[p][i][j] = (FLOAT) 0.0;

    for (int p=0; p<npart; p++) mask[p] = true;
    for (int neib=0; neib<Nneib; neib++) {
      // The hrangesqd used for MaskParticleNeigbourGather is the cell
      // hrangesqsd broadcast to an array.
      neibmanager.MaskParticleNeibGather(level, iorig, r, npart, neib, densmask, _hrangesqd, mask,
					 dr, drsqd, culled, smoothed_grav, direct);

      const DensityParticle& neibpart = neibmanager.GetNeib(neib);
      for (int p=0; p<npart; p++)
	if (mask[p] && culled[p]) w[p] = neibpart.m * hfac[p] * kern.w1(invh[p] * sqrt(drsqd[p]));

      for (int p=0; p<npart; p++)
	if (mask[p] && culled[p])
	  for (int j=0; j < ndim; j++)
	    for (int k=0; k < ndim; k++) {
	      rr[p][j][k] += w[p] * dr[j][p] * dr[k][p];
	      dv[p][j][k] += w[p] * dr[j][p] * (neibpart.v[k] - v[k][p]);
	      da[p][j][k] += w[p] * dr[j][p] * (neibpart.a[k] - a[k][p]);
	    }
    }

    for (int p=0; p<npart; p++) {
      FLOAT dvdx[ndim][ndim];
      FLOAT dadx[ndim][ndim];
      for (int i=0; i<ndim; ++i)
	for (int j=0; j<ndim; ++j)
	  dadx[i][j] = dvdx[i][j] = (FLOAT) 0.0;

      // Invert the rr matrix and compute the gradients
      FLOAT T[ndim][ndim];
      InvertMatrix(rr[p], T);

      // Check the accuracy of the integral gradients (using the square of the condition number),
      // if it's bad, we'll set alpha_loc to alpha_max.
      double modR(0), modT(0) ;
      for (int i=0; i<ndim; i++)
	for (int j=0; j<ndim; j++){
	  modR += rr[p][i][j]*rr[p][i][j];
	  modT +=  T[i][j]* T[i][j];
	}
      double sqd_condition_number = modR*modT / (ndim*ndim) ;

      FLOAT alpha_loc = (FLOAT) 0.0;
      if (sqd_condition_number > (double) 1.0e4) {
	// Bad gradients
	alpha_loc = alpha_visc ;
      } else {
	// Ok gradients
	for (int i=0; i<ndim; i++)
	  for (int j=0; j<ndim; j++)
	    for (int k=0; k<ndim; k++) {
	      dvdx[i][j] += T[j][k] * dv[p][k][i];
	      dadx[i][j] += T[j][k] * da[p][k][i];
	    }

	// Now compute the components needed for the limiter:
	FLOAT ddivdt = 0;
	FLOAT divv2 = 0;
	for (int i=0; i<ndim; ++i) {
	  ddivdt += dadx[i][i] ;
	  for (int j=0; j<ndim; ++j)
	    ddivdt -= dvdx[i][j]*dvdx[j][i];

	  divv2 += dvdx[i][i] ;
	}

	divv2 *= divv2;
	FLOAT curlv2 = CurlVelSqd(dvdx) ;

	FLOAT f_balsara = 1 ;
	if (curlv2 > 0)
	  f_balsara = (divv2 / (divv2 + curlv2)) ;

	if (ddivdt < 0) {
	  alpha_loc = (10 * h[p]*h[p] / (sound[p]*sound[p])) * f_balsara * (-ddivdt) ;
	  alpha_loc = min(alpha_loc, alpha_visc) ;
	}
      }

      if (alpha_loc > alpha[p])
	alpha[p] = alpha_loc ;

      dalphadt[p] = (FLOAT) 0.1 * sound[p]*(max(alpha_visc_min, alpha_loc) - alpha[p])*invh[p];
    }

    // copy SIMDs back to parti
    for (int p=0; p<npart; p++) parti[p].alpha = alpha[p];
    for (int p=0; p<npart; p++) parti[p].dalphadt = dalphadt[p];
  }

  // If any h is invalid (i.e. larger than maximum h), then return error code
  // (0).
  all = true;
  for (int p=0; p<npart; p++) all = all && h[p] <= hmax;
  // This is safer - the previous two functions don't change parti[p].h, but
  // they could:
  //for (int p=0; p<npart; p++) all = all && parti[p].h <= hmax;
  if (all) return 1;
  else return 0;
}



//=================================================================================================
//  GradhSph::ComputeThermalProperties
/// Compute all thermal properties for grad-h SPH method for given particle.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void GradhSph<ndim, kernelclass>::ComputeThermalProperties
 (SphParticle<ndim> &part_gen)         ///< [inout] Particle i data
{
  GradhSphParticle<ndim>& part = static_cast<GradhSphParticle<ndim> &> (part_gen);

  if (part.ptype != dust_type) {
    part.u       = eos->SpecificInternalEnergy(part);
    part.sound   = eos->SoundSpeed(part);
    part.pressure = eos->Pressure(part);
  }

  return;
}



//=================================================================================================
//  GradhSph::ComputeSphHydroForces
/// Compute SPH neighbour force pairs for
/// (i) All neighbour interactions of particle i with i.d. j > i,
/// (ii) Active neighbour interactions of particle j with i.d. j > i
/// (iii) All inactive neighbour interactions of particle i with i.d. j < i.
/// This ensures that all particle-particle pair interactions are only
/// computed once only for efficiency.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void GradhSph<ndim, kernelclass>::ComputeSphHydroForces
 (GradhSphParticle<ndim> &parti,                                   ///< [inout] Particle i data
  NeighbourList<typename GradhSphBase<ndim>::HydroNeib>& neibpart) ///< [in] List of neighbours
{
  int k;                               // Dimension counter
  FLOAT alpha_mean;                    // Mean articial viscosity alpha value
  FLOAT draux[ndim];                   // Relative position vector
  FLOAT dvdr;                          // Dot product of dv and dr
  FLOAT wkerni;                        // Value of w1 kernel function for part i
  FLOAT wkernj;                        // Value of w1 kernel function for neighbour j
  FLOAT vsignal;                       // Signal velocity
  FLOAT paux;                          // Aux. pressure force variable
  FLOAT winvrho;                       // 0.5*(wkerni + wkernj)*invrhomean

  // Some basic sanity-checking in case of invalid input into routine
  assert(!parti.flags.is_dead());

  const FLOAT invh_i   = 1/parti.h;
  const FLOAT invrho_i = 1/parti.rho;

  // Loop over all potential neighbours in the list
  //-----------------------------------------------------------------------------------------------
  int Nneib = neibpart.size() ;
  for (int j=0; j<Nneib; j++) {
    assert(!neibpart[j].flags.is_dead());

    const FLOAT invh_j   = 1/neibpart[j].h;
    const FLOAT invrho_j = 1/neibpart[j].rho;

    for (k=0; k<ndim; k++) draux[k] = neibpart[j].r[k]-parti.r[k];
    const FLOAT drmag = sqrt(DotProduct(draux,draux,ndim));
    if (drmag>0) for (k=0; k<ndim; k++) draux[k] /= drmag;

    wkerni = parti.hfactor*kern.w1(drmag*invh_i);
    wkernj = neibpart[j].hfactor*kern.w1(drmag*invh_j);

    dvdr = DotProduct(neibpart[j].v,draux, ndim);
    dvdr -= DotProduct(parti.v,draux,ndim);

    // Add contribution to velocity divergence
    parti.div_v -= neibpart[j].m*dvdr*wkerni;

    // Main SPH pressure force term
    paux = ((parti.pressure*parti.invomega)/(parti.rho*parti.rho))*wkerni +
      ((neibpart[j].pressure*neibpart[j].invomega)/(neibpart[j].rho*neibpart[j].rho))*wkernj;

    // Add dissipation terms (for approaching particle pairs)
    //---------------------------------------------------------------------------------------------
    if (dvdr < (FLOAT) 0.0) {

      winvrho = (FLOAT) 0.25*(wkerni + wkernj)*(invrho_i + invrho_j);

      // Artificial viscosity term
      if (avisc == mon97) {
        vsignal    = parti.sound + neibpart[j].sound - beta_visc*alpha_visc*dvdr;
        paux       -= alpha_visc*vsignal*dvdr*winvrho;
        parti.dudt -= 0.5*neibpart[j].m*alpha_visc*vsignal*dvdr*dvdr*winvrho;
      }
      else if (avisc == mon97mm97 || avisc == mon97cd2010) {
        alpha_mean = (FLOAT) 0.5*(parti.alpha + neibpart[j].alpha);
        vsignal    = parti.sound + neibpart[j].sound - beta_visc*alpha_mean*dvdr;
        paux       -= alpha_mean*vsignal*dvdr*winvrho;
        parti.dudt -= 0.5*neibpart[j].m*alpha_mean*vsignal*dvdr*dvdr*winvrho;
      }

      // Artificial conductivity term
      if (acond == wadsley2008) {
        parti.dudt += neibpart[j].m*dvdr*(neibpart[j].u - parti.u)*
          (invrho_i*wkerni + invrho_j*wkernj);
      }
      else if (acond == price2008) {
        parti.dudt += (FLOAT) 0.5*neibpart[j].m*(parti.u - neibpart[j].u)*
          winvrho*(invrho_i + invrho_j)*
          sqrt(fabs(parti.pressure - neibpart[j].pressure));
      }

    }
    //---------------------------------------------------------------------------------------------

    // Add total hydro contribution to acceleration for particle i
    for (k=0; k<ndim; k++) {
    	parti.a[k] += neibpart[j].m*draux[k]*paux;
    	assert(parti.a[k]==parti.a[k]);
    }
    parti.levelneib = max(parti.levelneib,neibpart[j].level);
    neibpart[j].levelneib = max(neibpart[j].levelneib,parti.level);

  }
  //-----------------------------------------------------------------------------------------------

  // Set velocity divergence and compressional heating rate terms
  parti.div_v    *= invrho_i;
  parti.dudt     -= eos->Pressure(parti)*parti.div_v*invrho_i*parti.invomega;
  if (tdavisc == mm97) {
    parti.dalphadt = (FLOAT) 0.1*parti.sound*(alpha_visc_min - parti.alpha)*invh_i +
        max(-parti.div_v, (FLOAT) 0.0)*(alpha_visc - parti.alpha);
  }

  return;
}



//=================================================================================================
//  GradhSph::ComputeSphHydroGravForces
/// Compute SPH neighbour force pairs for
/// (i) All neighbour interactions of particle i with i.d. j > i,
/// (ii) Active neighbour interactions of particle j with i.d. j > i
/// (iii) All inactive neighbour interactions of particle i with i.d. j < i.
/// This ensures that all particle-particle pair interactions are only
/// computed once only for efficiency.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void GradhSph<ndim, kernelclass>::ComputeSphHydroGravForces
(GradhSphParticle<ndim> &parti,                                   ///< [inout] Particle i data
 NeighbourList<typename GradhSphBase<ndim>::HydroNeib>& neibpart) ///< [in] List of neighbours
{
  int k;                               // Dimension counter
  FLOAT alpha_mean;                    // Mean articial viscosity alpha value
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drmag;                         // Distance
  FLOAT dv[ndim];                      // Relative velocity vector
  FLOAT dvdr;                          // Dot product of dv and dr
  FLOAT invdrmag;                      // 1 / distance
  FLOAT wkerni;                        // Value of w1 kernel function
  FLOAT wkernj;                        // Value of w1 kernel function
  FLOAT vsignal;                       // Signal velocity
  FLOAT paux;                          // Aux. pressure force variable
  FLOAT winvrho;                       // 0.5*(wkerni + wkernj)*invrhomean

  FLOAT invh_i   = 1/parti.h;
  FLOAT invrho_i = 1/parti.rho;
  const FLOAT invhsqdi = invh_i*invh_i;

  // Loop over all potential neighbours in the list
  //-----------------------------------------------------------------------------------------------
  int Nneib = neibpart.size() ;
  for (int j=0; j<Nneib; j++) {
    assert(!neibpart[j].flags.is_dead());
    assert(neibpart[j].h > (FLOAT) 0.0);
    assert(neibpart[j].rho > (FLOAT) 0.0);
    assert(neibpart[j].m > (FLOAT) 0.0);

    for (k=0; k<ndim; k++) dr[k] = neibpart[j].r[k] - parti.r[k];
    for (k=0; k<ndim; k++) dv[k] = neibpart[j].v[k] - parti.v[k];
    drmag = sqrt(DotProduct(dr,dr,ndim) + small_number);
    invdrmag = (FLOAT) 1.0/drmag;
    for (k=0; k<ndim; k++) dr[k] *= invdrmag;
    dvdr = DotProduct(dv,dr,ndim);

    FLOAT invh_j   = 1/neibpart[j].h;
    FLOAT invrho_j = 1/neibpart[j].rho;

    wkerni = parti.hfactor*kern.w1(drmag*invh_i);
    wkernj = neibpart[j].hfactor*kern.w1(drmag*invh_j);

    // Main SPH pressure force term
    paux = ((parti.pressure*parti.invomega)/(parti.rho*parti.rho))*wkerni +
      ((neibpart[j].pressure*neibpart[j].invomega)/(neibpart[j].rho*neibpart[j].rho))*wkernj;

    // Add dissipation terms (for approaching particle pairs)
    //---------------------------------------------------------------------------------------------
    if (dvdr < (FLOAT) 0.0) {

      winvrho = (FLOAT) 0.25*(wkerni + wkernj)*(invrho_i + invrho_j) ;

      // Artificial viscosity term
      if (avisc == mon97) {
        vsignal     = parti.sound + neibpart[j].sound - beta_visc*alpha_visc*dvdr;
        paux       -= alpha_visc*vsignal*dvdr*winvrho;
        parti.dudt -= 0.5*neibpart[j].m*alpha_visc*vsignal*dvdr*dvdr*winvrho;
      }
      else if (avisc == mon97mm97 || avisc == mon97cd2010) {
        alpha_mean  = (FLOAT) 0.5*(parti.alpha + neibpart[j].alpha);
        vsignal     = parti.sound + neibpart[j].sound - beta_visc*alpha_mean*dvdr;
        paux       -= alpha_mean*vsignal*dvdr*winvrho;
        parti.dudt -= 0.5*neibpart[j].m*alpha_mean*vsignal*dvdr*dvdr*winvrho;
      }

      // Artificial conductivity term
      if (acond == wadsley2008) {
        parti.dudt += neibpart[j].m*dvdr*(neibpart[j].u - parti.u)*
          (invrho_i*wkerni + invrho_j*wkernj);
      }
      else if (acond == price2008) {
        parti.dudt += (FLOAT) 0.5*neibpart[j].m*(parti.u - neibpart[j].u)*
          winvrho*(invrho_i + invrho_j)*
          sqrt(fabs(parti.pressure - neibpart[j].pressure));
      }

    }
    //---------------------------------------------------------------------------------------------

    // Add total hydro contribution to acceleration for particle i
    for (k=0; k<ndim; k++) parti.a[k] += neibpart[j].m*dr[k]*paux;


    // Compute SPH gravity terms
    //---------------------------------------------------------------------------------------------
    paux = (FLOAT) 0.5*(invhsqdi*kern.wgrav(drmag*invh_i) + parti.zeta*wkerni +
                        invh_j*invh_j*kern.wgrav(drmag*invh_j) +
                        neibpart[j].zeta*wkernj);
    for (k=0; k<ndim; k++) parti.atree[k] += neibpart[j].m*dr[k]*paux;
    parti.gpot += (FLOAT) 0.5*neibpart[j].m*(invh_i*kern.wpot(drmag*invh_i) +
                                             invh_j*kern.wpot(drmag*invh_j));

    // Add velocity divergence and neighbour timestep level terms
    parti.div_v -= neibpart[j].m*dvdr*wkerni;
    parti.levelneib = max(parti.levelneib,neibpart[j].level);
    neibpart[j].levelneib = max(neibpart[j].levelneib,parti.level);

  }
  //-----------------------------------------------------------------------------------------------


  // Set velocity divergence and compressional heating rate terms
  parti.div_v   *= invrho_i;
  parti.dudt    -= eos->Pressure(parti)*parti.div_v*invrho_i*parti.invomega;
  if (tdavisc == mm97) {
    parti.dalphadt = (FLOAT) 0.1*parti.sound*(alpha_visc_min - parti.alpha)*invh_i +
        max(parti.div_v,(FLOAT) 0.0)*(alpha_visc - parti.alpha);
  }

  return;
}



//=================================================================================================
//  GradhSph::ComputeSphGravForces
/// Compute SPH neighbour force pairs for
/// (i) All neighbour interactions of particle i with i.d. j > i,
/// (ii) Active neighbour interactions of particle j with i.d. j > i
/// (iii) All inactive neighbour interactions of particle i with i.d. j < i.
/// This ensures that particle-particle pair interactions are computed once only for efficiency.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void GradhSph<ndim, kernelclass>::ComputeSphGravForces
(GradhSphParticle<ndim> &parti,                                   ///< [inout] Particle i data
 NeighbourList<typename GradhSphBase<ndim>::HydroNeib>& neibpart) ///< [in] List of neighbours
{
  int k;                               // Dimension counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drmag;                         // Distance
  //FLOAT dv[ndim];                      // Relative velocity vector
  //FLOAT dvdr;                          // Dot product of dv and dr
  FLOAT invdrmag;                      // 1 / distance
  FLOAT gaux;                          // Aux. grav. potential variable
  FLOAT paux;                          // Aux. pressure force variable

  FLOAT invh_i = 1/parti.h;

  // Loop over all potential neighbours in the list
  //-----------------------------------------------------------------------------------------------
  int Nneib = neibpart.size();
  for (int j=0; j<Nneib; j++) {
    assert(!neibpart[j].flags.is_dead());

    for (k=0; k<ndim; k++) dr[k] = neibpart[j].r[k] - parti.r[k];
    //for (k=0; k<ndim; k++) dv[k] = neibpart[j].v[k] - parti.v[k];
    //dvdr = DotProduct(dv,dr,ndim);
    drmag = sqrt(DotProduct(dr,dr,ndim) + small_number);
    invdrmag = (FLOAT) 1.0/drmag;
    for (k=0; k<ndim; k++) dr[k] *= invdrmag;

    FLOAT invh_j =  1/neibpart[j].h;

    // Main SPH gravity terms
    //---------------------------------------------------------------------------------------------
    paux = (FLOAT) 0.5*(invh_i*invh_i*kern.wgrav(drmag*invh_i) +
                        parti.zeta*parti.hfactor*kern.w1(drmag*invh_i) +
                        invh_j*invh_j*kern.wgrav(drmag*invh_j) +
                        neibpart[j].zeta*neibpart[j].hfactor*kern.w1(drmag*invh_j));
    gaux = (FLOAT) 0.5*(invh_i*kern.wpot(drmag*invh_i) +
                        invh_j*kern.wpot(drmag*invh_j));

    // Add total hydro contribution to acceleration for particle i
    for (k=0; k<ndim; k++) parti.atree[k] += neibpart[j].m*dr[k]*paux;
    parti.gpot += neibpart[j].m*gaux;

    parti.levelneib = max(parti.levelneib,sphdata[j].level);
  }

  //===============================================================================================

  return;
}



//=================================================================================================
//  GradhSph::ComputeDirectGravForces
/// Compute the contribution to the total gravitational force of particle 'i'
/// due to 'Nneib' neighbouring particles in the list 'neiblist'.
//=================================================================================================
template <int ndim>
void GradhSphBase<ndim>::ComputeDirectGravForces
 (GradhSphParticle<ndim>& parti,                                      ///< Particle i data
  NeighbourList<typename GradhSphBase<ndim>::DirectNeib>& gravdata)  ///< Neighbour particle data
{
  int k;                               // Dimension counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drsqd;                         // Distance squared
  FLOAT invdrmag;                      // 1 / distance
  FLOAT invdr3;                        // 1 / dist^3


  // Loop over all neighbouring particles in list
  //-----------------------------------------------------------------------------------------------
  int Ndirect = gravdata.size();
  for (int j=0; j<Ndirect; j++) {
    assert(!gravdata[j].flags.is_dead());

    for (k=0; k<ndim; k++) dr[k] = gravdata[j].r[k] - parti.r[k];
    drsqd    = DotProduct(dr,dr,ndim) + small_number;
    invdrmag = (FLOAT) 1.0/sqrt(drsqd);
    invdr3   = invdrmag*invdrmag*invdrmag;

    // Add contribution to current particle
    for (k=0; k<ndim; k++) parti.atree[k] += gravdata[j].m*dr[k]*invdr3;
    parti.gpot += gravdata[j].m*invdrmag;

    // Sanity-checkt to ensure particles are really un-softened direct-sum neighbours
    assert(drsqd >= parti.hrangesqd && drsqd >= gravdata[j].hrangesqd);

  }
  //-----------------------------------------------------------------------------------------------

  return;
}



//=================================================================================================
//  GradhSph::ComputeStarGravForces
/// Computes contribution of gravitational force and potential due to stars.
//=================================================================================================
template <int ndim, template<int> class kernelclass>
void GradhSph<ndim, kernelclass>::ComputeStarGravForces
 (const int N,                         ///< [in] No. of stars
  NbodyParticle<ndim> **nbodydata,     ///< [in] Array of star pointers
  SphParticle<ndim> &part)             ///< [inout] SPH particle reference
{
  int j;                               // Star counter
  int k;                               // Dimension counter
  FLOAT dr[ndim];                      // Relative position vector
  FLOAT drmag;                         // Distance
  FLOAT drsqd;                         // Distance squared
  //FLOAT drdt;                          // Rate of change of relative distance
  //FLOAT dv[ndim];                      // Relative velocity vector
  FLOAT invdrmag;                      // 1 / drmag
  FLOAT invhmean;                      // 1 / hmean
  FLOAT ms;                            // Star mass
  FLOAT paux;                          // Aux. force variable
  GradhSphParticle<ndim>& parti = static_cast<GradhSphParticle<ndim>& > (part);

  // Loop over all stars and add each contribution
  //-----------------------------------------------------------------------------------------------
  for (j=0; j<N; j++) {

    if (fixed_sink_mass) ms = msink_fixed;
    else ms = nbodydata[j]->m;

    for (k=0; k<ndim; k++) dr[k] = nbodydata[j]->r[k] - parti.r[k];
    //for (k=0; k<ndim; k++) dv[k] = nbodydata[j]->v[k] - parti.v[k];
    drsqd    = DotProduct(dr,dr,ndim) + small_number;
    drmag    = sqrt(drsqd);
    invdrmag = (FLOAT) 1.0/drmag;
    invhmean = (FLOAT) 2.0/(parti.h + nbodydata[j]->h);
    //drdt     = DotProduct(dv,dr,ndim)*invdrmag;
    paux     = ms*invhmean*invhmean*kern.wgrav(drmag*invhmean)*invdrmag;

    // Add total hydro contribution to acceleration for particle i
    for (k=0; k<ndim; k++) parti.atree[k] += paux*dr[k];
    parti.gpot += ms*invhmean*kern.wpot(drmag*invhmean);

    assert(invhmean > (FLOAT) 0.0);

  }
  //-----------------------------------------------------------------------------------------------

  return;
}

#if defined MPI_PARALLEL
template <int ndim, template<int> class kernelclass>
void GradhSph<ndim, kernelclass>::FinishReturnExport () {
  if (tdavisc == mm97) {
    for (int i=0; i<Nhydro; i++) {
      GradhSphParticle<ndim>& part = sphdata[i];
      part.dalphadt = (FLOAT) 0.1*part.sound*(alpha_visc_min - part.alpha)/part.h +
          max(-part.div_v, (FLOAT) 0.0)*(alpha_visc - part.alpha);
    }
  }
}
#endif


template class GradhSphBase<1>;
template class GradhSphBase<2>;
template class GradhSphBase<3>;
template class GradhSph<1, M4Kernel>;
template class GradhSph<2, M4Kernel>;
template class GradhSph<3, M4Kernel>;
template class GradhSph<1, QuinticKernel>;
template class GradhSph<2, QuinticKernel>;
template class GradhSph<3, QuinticKernel>;
template class GradhSph<1, GaussianKernel>;
template class GradhSph<2, GaussianKernel>;
template class GradhSph<3, GaussianKernel>;
template class GradhSph<1, TabulatedKernel>;
template class GradhSph<2, TabulatedKernel>;
template class GradhSph<3, TabulatedKernel>;
