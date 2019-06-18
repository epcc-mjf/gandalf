#include <aligned_new>
/*
 * NeighbourManager.h
 *
 *  Created on: 5 Dec 2016
 *      Author: rosotti
 */

#ifndef NEIGHBOURMANAGER_H_
#define NEIGHBOURMANAGER_H_


#ifdef INTEL_INTRINSICS
#include <immintrin.h>
#endif

#include <iterator>
#include <vector>
using namespace std;
#include "SIMD.h"
#include "TreeCell.h"
#include "GhostNeighbours.hpp"
#include "NeighbourManagerBase.h"
#include "Particle.h"


// Fwd declare
template<int ndim, class ParticleType> class NeighbourManager;


//=================================================================================================
//  Class NeighbourIterator
/// \brief   Template class for iterating of neighbour lists.
/// \details Random Access Iterator. Iterates over members of a neighbour list.
/// \author  R. A. Booth
/// \date    20/12/2016
//=================================================================================================
template <class ParticleType>
class NeighbourIterator :
    public std::iterator<std::random_access_iterator_tag, ParticleType>
{
  typedef std::iterator<std::random_access_iterator_tag, ParticleType> iterator_type ;
public:
  typedef typename iterator_type::reference reference ;
  typedef typename iterator_type::pointer   pointer ;

  NeighbourIterator(int* p_idx, ParticleType* p_part) :
    _p_idx(p_idx), _p_part(p_part) { };

  // Dereference
  reference operator*() const {
    return _p_part[*_p_idx];
  }
  reference operator[](std::size_t n) const {
    return _p_part[_p_idx[n]];
  }
  pointer operator->() const {
    return &(_p_part[*_p_idx]);
  }

  // Increment / Decrement
  NeighbourIterator<ParticleType>& operator++(int) {
    ++_p_idx;
    return *this;
  }
  NeighbourIterator<ParticleType>& operator--(int) {
    --_p_idx;
    return *this;
  }
  NeighbourIterator<ParticleType> operator++() {
    return NeighbourIterator<ParticleType>(_p_idx++, _p_part);
  }
  NeighbourIterator<ParticleType> operator--() {
    return NeighbourIterator<ParticleType>(_p_idx--, _p_part);
  }

  // Add / subtract
  NeighbourIterator<ParticleType>& operator+=(std::size_t n) {
    _p_idx += n;
    return *this;
  }
  NeighbourIterator<ParticleType>& operator-=(std::size_t n) {
    _p_idx -= n;
    return *this;
  }

  // Comparison
  bool operator==(const NeighbourIterator<ParticleType>& other) const {
    assert(_p_part == other._p_part);
    return _p_idx == other._p_idx;
  }
  bool operator!=(const NeighbourIterator<ParticleType>& other) const {
    return !(*this == other);
  }
  bool operator<(const NeighbourIterator<ParticleType>& other) const {
    assert(_p_part == other._p_part);
    return _p_idx < other._p_idx;
  }
  bool operator>(const NeighbourIterator<ParticleType>& other) const {
    assert(_p_part == other._p_part);
    return _p_idx < other._p_idx;
  }
  bool operator<=(const NeighbourIterator<ParticleType>& other) const {
    return !(*this > other);
  }
  bool operator>=(const NeighbourIterator<ParticleType>& other) const {
    return (*this < other);
  }


private:
  int* _p_idx;
  ParticleType* _p_part;
};

// Arithmetic operators should be free functions:
template<class ParticleType>
NeighbourIterator<ParticleType> operator+(NeighbourIterator<ParticleType>& p,
                                          std::size_t n) {
  return NeighbourIterator<ParticleType>(p) += n ;
}
template<class ParticleType>
NeighbourIterator<ParticleType> operator+(std::size_t n,
                                          NeighbourIterator<ParticleType>& p) {
  return NeighbourIterator<ParticleType>(p) += n ;
}
template<class ParticleType>
NeighbourIterator<ParticleType> operator-(NeighbourIterator<ParticleType>& p,
                                          std::size_t n) {
  return NeighbourIterator<ParticleType>(p) -= n ;
}
template<class ParticleType>
NeighbourIterator<ParticleType> operator-(std::size_t n,
                                          NeighbourIterator<ParticleType>& p) {
  return NeighbourIterator<ParticleType>(p) -= n ;
}


//=================================================================================================
//  Class NeighbourList
/// \brief   Template container class which holds a list of neighbour lists.
/// \details Allows the neighbour lists to be fully abstracted.
/// \author  R. A. Booth
/// \date    20/12/2016
//=================================================================================================
template<class ParticleType>
class NeighbourList {
public:
  NeighbourList(std::vector<int>& idx,
      std::vector<ParticleType>& neib) :
        _idx(idx), _neibpart(neib) { };

  NeighbourIterator<ParticleType> begin() {
    return NeighbourIterator<ParticleType>(&_idx.front(), &_neibpart.front());
  }
  NeighbourIterator<ParticleType> end() {
    return NeighbourIterator<ParticleType>(&_idx.back(), &_neibpart.front());
  }

  std::size_t size() const {
    return _idx.size();
  }
  const ParticleType& operator[](std::size_t i) const {
    return _neibpart[_idx[i]];
  }
  ParticleType& operator[](std::size_t i) {
    return _neibpart[_idx[i]];
  }

  template <class InParticleType>
  void VerifyNeighbourList(int i, int Ntot, const InParticleType& partdata,
                           const string& searchmode) ;

private:
  template<int ndim, class PartType> friend class NeighbourManager;

  std::vector<int>& _idx;
  std::vector<ParticleType>& _neibpart;
};

struct ListLength {
  int Nhydro;
  int Ndirect;
  int Ngrav;
};


template <class ParticleType>
struct GravityNeighbourLists {
  typedef ParticleType DirectType ;

  NeighbourList<ParticleType> neiblist;
  NeighbourList<ParticleType> smooth_gravlist;
  NeighbourList<DirectType>   directlist;
};

//=================================================================================================
//  Class NeighbourManager
/// \brief   Template class for neighbour searches.
/// \details Handles issues related to positioning of ghosts, along with culling the neighbour list
///          down to only those needed in a tree walk.
/// \author  G. Rosotti, R. A. Booth
/// \date    15/12/2016
//=================================================================================================
template <int ndim, class ParticleType>
class NeighbourManager : public NeighbourManagerDim<ndim> {
private:
  int _NPeriodicGhosts;
  int _NCellDirectNeib;
  int _NNonInteract;
  vector<int> neiblist;
  vector<int> directlist;
  vector<int> neib_idx;
  vector<ParticleType> neibdata;

  vector<int> culled_neiblist;
  vector<int> smoothgravlist;

  vector<FLOAT> dr;
  vector<FLOAT> drmag;

  const DomainBox<ndim>* _domain;
  const ParticleTypeRegister* _types;
  double _kernrange;

  template<bool __v>
  struct __boolean_constant {
    static const bool value = __v ;
    typedef bool value_type ;
    typedef __boolean_constant<__v> type ;
  };

  typedef __boolean_constant<true>  _true_type ;
  typedef __boolean_constant<false> _false_type ;

#ifdef INTEL_INTRINSICS
  static const int AHEAD_E = 0; /// Number of particles ahead to prefetch when
                                /// accessing the Particle array for _EndSearch.
                                /// 0 switches off prefetching.
  static const int AHEAD_T = 0; /// Number of particles ahead to prefetch when
                                /// accessing the DensityParticle or
                                /// HydroForcesParticle vectors for
                                /// TrimNeighbourLists.  0 switches off
                                /// prefetching.
#endif

public:
  using NeighbourManagerDim<ndim>::tempneib;
  using NeighbourManagerDim<ndim>::tempperneib;
  using NeighbourManagerDim<ndim>::tempdirectneib;
  using NeighbourManagerDim<ndim>::gravcell;

  NeighbourManager(const Hydrodynamics<ndim>* hydro, const DomainBox<ndim>& domain)
  : _NPeriodicGhosts(0), _NCellDirectNeib(0), _NNonInteract(0),
    _domain(&domain),
    _types(&(hydro->types)),
    _kernrange(hydro->kernp->kernrange)
  { };

  NeighbourManager(const ParticleTypeRegister& types, double kernrange,
      const DomainBox<ndim>& domain)
  : _NPeriodicGhosts(0), _NCellDirectNeib(0), _NNonInteract(0),
    _domain(&domain),
    _types(&types),
    _kernrange(kernrange)
  { };

  void set_target_cell(const TreeCellBase<ndim>& cell) {
    NeighbourManagerDim<ndim>::set_target_cell(cell);
    neiblist.clear();
    neibdata.clear();
    neib_idx.clear();
    directlist.clear();
  }

  //===============================================================================================
  // EndSearch
  /// \brief Collect particle data needed for the neighbours and cull distant particles that do
  ///        not interact with the cell hydrodynamically
  //===============================================================================================
  template<class InParticleType>
  void EndSearch(const TreeCellBase<ndim> &cell, const InParticleType* partdata, const int Npartmax) {
    _EndSearch<InParticleType,_false_type>(cell, partdata, Npartmax, false);
  }

  //===============================================================================================
  // EndSearchGather
  /// \brief Collect particle data needed for the neighbours and cull distant particles that do
  ///        not contribute to the density
  //===============================================================================================
  template<class InParticleType>
  void EndSearchGather(const TreeCellBase<ndim> &cell, const InParticleType* partdata, const int Npartmax) {
    _EndSearch<InParticleType,_true_type>(cell, partdata, Npartmax, false);
  }
  //===============================================================================================
  // EndSearchGravity
  /// \brief Collect particle data needed for the neighbours, and demote particles that do not
  ///        interact hydrodynamically to the directlist.
  //===============================================================================================
  template<class InParticleType>
  void EndSearchGravity(const TreeCellBase<ndim> &cell, const InParticleType* partdata, const int Npartmax) {
    _EndSearch<InParticleType,_false_type>(cell, partdata, Npartmax, true);
  }

  //===============================================================================================
  //  GetNumAllNeib
  /// \brief Return the total number of all particles found in the walk.
  //===============================================================================================
  int GetNumAllNeib() {
    return neib_idx.size();
  }
  //===============================================================================================
  //  GetNumNeib
  /// \brief Return the total number of neighbours.
  //===============================================================================================
  int GetNumNeib() {
    return neiblist.size();
  }
  //===============================================================================================
  //  GetNumPeriodicGhosts
  /// \brief Return the number of periodic ghosts.
  //===============================================================================================
  int GetNumPeriodicGhosts() {
    return _NPeriodicGhosts;
  }
  //===============================================================================================
  //  operator[]
  /// \brief Provide access to individual neighbours found in the tree walk.
  //===============================================================================================
  ParticleType& operator[](std::size_t i) {
    return neibdata[i] ;
  }
  const ParticleType& operator[](std::size_t i) const {
    return neibdata[i] ;
  }
  //===============================================================================================
  //  GetNeib
  /// \brief Get the reduced neighbour data for the i-th neighbour in the
  //  neighbour list.
  //===============================================================================================
  ParticleType& GetNeib(int i) {
    return neibdata[neiblist[i]];
  }
  //===============================================================================================
  //  GetNeibI
  /// \brief Get the index of thew i-th neighbour in the original particle array, and a pointer to
  ///        to the reduced neighbour data stored.
  //===============================================================================================
  std::pair<int,ParticleType*> GetNeibI(int i) {
    return make_pair(neib_idx[i],&neibdata[i]);
  }

  //===============================================================================================
  //  GetParticleNeib
  /// \brief    Get the list of particles that interact hydrodynamically with the Particle, p.
  /// \details  Gets a trimmed list of particles that interact with hydrodynamically with Particle,
  ///           p. This is determined by whether or not the particles smoothing spheres overlap,
  ///           and hydromask, which specifies the types needed (hydromask[ptype] == true for the
  ///           required particles). Also, if do_pair_once is included then the neighbour is only
  ///           included if this interaction will not have been already calculated in the update
  ///           for the neighbour itself.
  /// \returns  NeighbourList object, the list of neighbours found.
  //===============================================================================================
  template<class InParticleType>
  NeighbourList<ParticleType> GetParticleNeib
  (const InParticleType& p,                            ///< [in] Particle to collect the neibs for
   const Typemask& hydromask,                          ///< [in] Boolean flags listing types we need
   const bool do_pair_once)                            ///< [in]
  {
    if (do_pair_once)
      TrimNeighbourLists<InParticleType,_true_type,_false_type>(p, hydromask, p.hrangesqd, false);
    else
      TrimNeighbourLists<InParticleType,_false_type,_false_type>(p, hydromask, p.hrangesqd, false);

    return NeighbourList<ParticleType>(culled_neiblist, neibdata) ;
  }

  //===============================================================================================
  //  MaskParticleNeib
  /// \brief    Set the masks and data for the Particles for a single neighbour
  ///           that interacts hydrodynamically with them.
  /// \details  This is determined by whether or not the particles smoothing spheres overlap,
  ///           and hydromask, which specifies the types needed (hydromask[ptype] == true for the
  ///           required particles). Also, if do_pair_once is included then the neighbour is only
  ///           included if this interaction will not have been already calculated in the update
  ///           for the neighbour itself.
  /// \returns  Nothings.  The masks and data are returned in parameters
  //===============================================================================================
  template<class InParticleType>
  void MaskParticleNeib(const int level[MAX_NPART], const int iorig[MAX_NPART],
			const FLOAT r[ndim][MAX_NPART], const int npart,
			const int neib, const Typemask& hydromask, const bool do_pair_once,
			FLOAT draux[ndim][MAX_NPART], FLOAT drsqd[MAX_NPART],
			bool culled[MAX_NPART], bool smoothed_grav[MAX_NPART], bool direct[MAX_NPART])
  {
    if (do_pair_once)
      TrimNeighbour<InParticleType,_true_type,_false_type>
	(level, iorig, r, npart, neib, hydromask, hrangesqd, false,
	 draux, drsqd, culled, smoothed_grav, direct)

	p, hydromask, p.hrangesqd, false);
    else
      TrimNeighbour<InParticleType,_false_type,_false_type>(p, hydromask, p.hrangesqd, false);
  }

  //===============================================================================================
  //  GetParticleNeibGather
  /// \brief    Get the list of particles within a given smoothing range of the target particle.
  /// \details  Gets a trimmed list of particles that that are within hrangeqd of the particle,
  ///           for a density calculation.
  /// \returns  NeighbourList object, the list of neighbours found.
  //===============================================================================================
  template<class InParticleType>
  NeighbourList<ParticleType> GetParticleNeibGather
  (const InParticleType& p,                            ///< [in] Particle to collect the neibs for
   const Typemask& hydromask,                          ///< [in] Boolean flags listing types we need
   const FLOAT hrangesqd)                              ///< [in] Maximum smoothing range for gather.
  {
    TrimNeighbourLists<InParticleType,_false_type,_true_type>(p, hydromask, hrangesqd, false);

    return NeighbourList<ParticleType>(culled_neiblist, neibdata) ;
  }

  //===============================================================================================
  //  MaskParticleNeibGather
  /// \brief    Set the masks and data for the Particles for a single neighbour
  ///           within a given smoothing range.
  /// \details  The masks and data are set if the neighbour is within hrangeqd of a particle,
  ///           for a density calculation.
  /// \returns  Nothing.  The masks and data are returned in parameters.
  //===============================================================================================
  template<class InParticleType>
  void MaskParticleNeibGather(const int level[MAX_NPART], const int iorig[MAX_NPART],
			      const FLOAT r[ndim][MAX_NPART], const int npart,
			      const int neib, const Typemask hydromask[MAX_NPART], const FLOAT hrangesqd,
			      FLOAT dr[ndim][MAX_NPART], FLOAT drsqd[MAX_NPART],
			      bool culled[MAX_NPART], bool smoothed_grav[MAX_NPART], bool direct[MAX_NPART])
  {
    TrimNeighbour<InParticleType,_false_type,_true_type>
      (level, iorig, r, npart, neib, hydromask, hrangesqd, false,
       dr, drsqd, culled, smoothed_grav, direct);
  }

  //===============================================================================================
  //  GetParticleNeibGravity
  /// \brief    Get the list of particles that interact with the Particle, p.
  /// \details  As with GetParticleNeib, this function returns the list of neighbours that interact
  ///           hydrodynamically with the Particle, p (neiblist). Additionally this function also
  ///           generatesthe list of particles needed for gravitational interactions, including the
  ///            list of smoothed (smooth_gravlist) and unsmoothed (directlist) contributions.
  /// \returns  GravityNeighbourLists object. This struct contains the three wraps the three lists
  ///           of particles found.
  //===============================================================================================
  template<class InParticleType>
  GravityNeighbourLists<ParticleType> GetParticleNeibGravity
  (const InParticleType& p,                            ///< [in] Particle to collect the neibs for
   const Typemask& hydromask)                          ///< [in] Type flags for hydro interactions
   {
    TrimNeighbourLists<InParticleType,_false_type,_false_type>(p, hydromask, p.hrangesqd, true);

    assert((int) (culled_neiblist.size()+directlist.size()+smoothgravlist.size()+_NNonInteract) ==
                  GetNumAllNeib());


    typedef typename GravityNeighbourLists<ParticleType>::DirectType DirectType ;
    GravityNeighbourLists<ParticleType> neiblists =
      { NeighbourList<ParticleType>(culled_neiblist, neibdata),
        NeighbourList<ParticleType>(smoothgravlist,  neibdata),
        NeighbourList<DirectType>  (directlist,      neibdata) } ;

    return neiblists ;
   }

private:

  //===============================================================================================
  //  _EndSearch
  /// \details Gathers the particle data for the cell, and trims the list of hydrodynamic particles.
  ///          If keep_direct is true, then these become direct-list gravitational neighbours,
  ///          otherwise they are discarded. Periodic corrections are included.
  //===============================================================================================
  template<class InParticleType, class gather_only>
  void _EndSearch(const TreeCellBase<ndim> &cell, const InParticleType* partdata, const int Npartmax,
                  bool keep_direct=true) {

    assert(partdata != NULL);


    FLOAT dr[ndim];                      // Relative position vector
    FLOAT drsqd;                         // Distance squared
    FLOAT rc[ndim];                      // Position of cell
    const FLOAT hrangemaxsqd = pow(cell.rmax + _kernrange*cell.hmax,2);
    const FLOAT rmax = cell.rmax;
    //for (int k=0; k<ndim; k++) rc[k] = cell.rcell[k];
    cell.ComputeCellCentre(rc);

    const GhostNeighbourFinder<ndim> GhostFinder(*_domain, cell);
    int MaxGhosts = GhostFinder.MaxNumGhosts();

    Typemask gravmask = _types->gravmask;

    assert(!keep_direct || MaxGhosts==1 ||
        (gravcell.size() == 0 && tempdirectneib.size() == 0));

    assert(neibdata.size()==0);
    assert(directlist.size()==0);

    // Now load the particles

    // We could have 3 neibdata vectors:  direct, periodic, and ordinary neighbours.

    // Start from direct neighbours
    if (keep_direct) {
      for (int ii=0; ii<(int) tempdirectneib.size(); ii++) {
        NeighbourManagerBase::range rng = tempdirectneib[ii] ;
        for (int i=rng.begin; i < rng.end; ++i) {
#ifdef INTEL_INTRINSICS
	  if (AHEAD_E != 0) {
	    if (i+AHEAD_E < Npartmax) {
	      // Properly this should be
	      // _mm_prefetch(reinterpret_cast<char const*> &partdata[i+AHEAD_E].flags, _MM_HINT_T1);
	      // but Intel 18 and onwards allow this without casting.
	      _mm_prefetch(&partdata[i+AHEAD_E].flags, _MM_HINT_T1);
	      // use a[0] to consistently use &
	      _mm_prefetch(&partdata[i+AHEAD_E].a[0], _MM_HINT_T1);
	      // Only up to here is needed for DensityParticle, but we don't
	      // know that this is a DensityParticle, so prefetch the parts
	      // needed by HydroForcesParticle as well.
	      _mm_prefetch(&partdata[i+AHEAD_E].hrangesqd, _MM_HINT_T1);
	    }
	  }
#endif
          const InParticleType& part = partdata[i];
          // Forget immediately: direct particles and particles that do not interact gravitationally
          if (part.flags.is_dead()) continue;
          if (!gravmask[part.ptype]) continue;

          // Now create the particle / ghosts
          _AddParticleAndGhosts(GhostFinder, part, gather_only());
          directlist.push_back(neibdata.size()-1);
          neib_idx.push_back(i);
        }
      }
    }

    assert(directlist.size() == neibdata.size() &&
           neib_idx.size()   == neibdata.size());

    // Now look at the hydro candidate neighbours
    // First the ones that need ghosts to be created on the fly
    int Nneib = directlist.size();
    for (int ii=0; ii<(int) tempperneib.size(); ii++) {
      NeighbourManagerBase::range rng = tempperneib[ii] ;
      for (int i=rng.begin; i < rng.end; ++i) {
#ifdef INTEL_INTRINSICS
	if (AHEAD_E != 0) {
	  if (i+AHEAD_E < Npartmax) {
	    // Properly this should be
	    // _mm_prefetch(reinterpret_cast<char const*> &partdata[i+AHEAD_E].flags, _MM_HINT_T1);
	    // but Intel 18 and onwards allow this without casting.
	    _mm_prefetch(&partdata[i+AHEAD_E].flags, _MM_HINT_T1);
	    // use a[0] to consistently use &
	    _mm_prefetch(&partdata[i+AHEAD_E].a[0], _MM_HINT_T1);
	    // Only up to here is needed for DensityParticle, but we don't know
	    // that this is a DensityParticle, so prefetch the parts needed by
	    // HydroForcesParticle as well.
	    _mm_prefetch(&partdata[i+AHEAD_E].hrangesqd, _MM_HINT_T1);
	  }
	}
#endif
        if (partdata[i].flags.is_dead()) continue;

        _AddParticleAndGhosts(GhostFinder, partdata[i], gather_only());

        while (Nneib < (int) neibdata.size()) {
          int Nmax = neibdata.size();
          for (int k=0; k<ndim; k++) dr[k] = neibdata[Nneib].r[k] - rc[k];
          drsqd = DotProduct(dr, dr, ndim);
          if (drsqd < hrangemaxsqd || _scatter_overlap(neibdata[Nneib], drsqd, rmax, gather_only())) {
            neiblist.push_back(Nneib);
            neib_idx.push_back(i);
            Nneib++;
          }
          else if (keep_direct && gravmask[neibdata[Nneib].ptype]) {
            directlist.push_back(Nneib);
            neib_idx.push_back(i);
            Nneib++;
          }
          else if (Nmax > Nneib) {
            Nmax--;
            if (Nmax > Nneib) neibdata[Nneib] = neibdata[Nmax];
            neibdata.resize(neibdata.size()-1);
          }
        }// Loop over Ghosts
      }
    }

    // Store the number of periodic particles in the neiblist
    _NPeriodicGhosts = neiblist.size();

    // Find those particles that do not need ghosts on the fly
    for (int ii=0; ii<(int) tempneib.size(); ii++) {
      NeighbourManagerBase::range rng = tempneib[ii] ;
      for (int i=rng.begin; i < rng.end; ++i) {
#ifdef INTEL_INTRINSICS
	if (AHEAD_E != 0) {
	  if (i+AHEAD_E < Npartmax) {
	    // Properly this should be
	    // _mm_prefetch(reinterpret_cast<char const*> &partdata[i+AHEAD_E].flags, _MM_HINT_T1);
	    // but Intel 18 and onwards allow this without casting.
	    _mm_prefetch(&partdata[i+AHEAD_E].flags, _MM_HINT_T1);
	    // use a[0] to consistently use &
	    _mm_prefetch(&partdata[i+AHEAD_E].a[0], _MM_HINT_T1);
	    // Only up to here is needed for DensityParticle, but we don't know
	    // that this is a DensityParticle, so prefetch the parts needed by
	    // HydroForcesParticle as well.
	    _mm_prefetch(&partdata[i+AHEAD_E].hrangesqd, _MM_HINT_T1);
	  }
	}
#endif
        for (int k=0; k<ndim; k++) dr[k] = partdata[i].r[k] - rc[k];
        drsqd = DotProduct(dr,dr,ndim);
	// This looks wrong.  Nneib is neibdata.size(), and neibdata[Nneib] will
	// not have any data until partdata[i] is pushed.  (And a
	// DensityParticle does not have an h.)  But here gather_only() is true,
	// and _scatter_overlap does nothing.  Should be
	//        if (drsqd < hrangemaxsqd || _scatter_overlap(partdata[i], drsqd, rmax, gather_only())) {
	// or omit the test
	//        if (drsqd < hrangemaxsqd) {
        if (drsqd < hrangemaxsqd || _scatter_overlap(neibdata[Nneib], drsqd, rmax, gather_only())) {
          neibdata.emplace_back(partdata[i]);
          neiblist.push_back(Nneib);
          neib_idx.push_back(i);
          Nneib++;
	// This looks wrong also and should be
	//        } else if (keep_direct && gravmask[partdata[i].ptype]) {
        } else if (keep_direct && gravmask[neibdata[Nneib].ptype]) {
          // Hydro candidates that fail the test get demoted to direct neighbours
          neibdata.emplace_back(partdata[i]);
          directlist.push_back(Nneib);
          neib_idx.push_back(i);
          Nneib++;
        }
      }
    }

    _NCellDirectNeib = directlist.size();
    assert(neibdata.size() == (neiblist.size() + directlist.size()));
  }


  //===============================================================================================
  //  TrimNeighbourLists
  /// \detail This function trims the neighbour lists for a given particle, determining whether
  ///         neighbours are needed for hydro or gravity. Periodic corrections are applied.
  //===============================================================================================
  template<class InParticleType, class do_pair_once, class gather_only>
  void TrimNeighbourLists(const InParticleType& p, const Typemask& hydromask, FLOAT hrangesqd,
                          bool keep_grav)
  {
    FLOAT rp[ndim];
    FLOAT draux[ndim];

    for (int k=0; k<ndim; k++) rp[k] = p.r[k];

    const GhostNeighbourFinder<ndim> GhostFinder(*_domain);

    Typemask gravmask = _types->gravmask;

    culled_neiblist.clear();
    smoothgravlist.clear();
    // Particles that are already in the directlist stay there; we just add the ones that were demoted
    directlist.resize(_NCellDirectNeib);
    _NNonInteract = 0;

    // Go through the hydro neighbour candidates and check the distance. The ones that are not real neighbours
    // are demoted to the direct list
    int imax = neibdata.size(); // Or can the compiler work out that the
				// neibdata size is not changed in the loop?
    for (int j=0; j<(int) neiblist.size(); j++) {
      int i = neiblist[j];
#ifdef INTEL_INTRINSICS
      if (AHEAD_T != 0) {
	// This should properly look ahead only in neiblist, i.e.,
	// int jmax = neiblist.size()
	// ...
	// if (j+AHEAD_T < jmax) {
	// _mm_prefetch(&neibdata[neiblist[j+AHEAD_T]].flags, _MM_HINT_T1);
	// _mm_prefetch(&neibdata[neiblist[j+AHEAD_T]].a[0], _MM_HINT_T1);
	// }
	// but that means an indirection, and may be slower than the few
	// unnecessary prefetches looking ahead in neibdata.
	if (i+AHEAD_T < imax) {
	  // Properly this should be
	  // _mm_prefetch(reinterpret_cast<char const*> &xxx, _MM_HINT_T1);
	  // but Intel 18 and onwards allow this without casting.
	  _mm_prefetch(&neibdata[i+AHEAD_T].flags, _MM_HINT_T1);
	  // use a[0] to consistently use &
	  _mm_prefetch(&neibdata[i+AHEAD_T].a[0], _MM_HINT_T1);
	  // This will give all of DensityParticle, and what we need of
	  // HydroForcesParticle.
	}
      }
#endif
      ParticleType& neibpart = neibdata[i] ;

      // If do_pair_once is true then only get the neighbour for the first of the two times the
      if (not _first_appearance(p, neibpart, do_pair_once())) continue;

      // Compute relative position and distance quantities for pair
      for (int k=0; k<ndim; k++) draux[k] = neibpart.r[k] - rp[k];
      if (j < _NPeriodicGhosts)
	//  This changes neibpart.r - is that correct?
        GhostFinder.ApplyPeriodicDistanceCorrection(neibpart.r, draux);

      const FLOAT drsqd = DotProduct(draux,draux,ndim);

      //if (drsqd <= small_number) continue;

      // Record if neighbour is direct-sum or and SPH neighbour.
      // If SPH neighbour, also record max. timestep level for neighbour
      // How does this work, neibpart is ParticleType but _scatter_overlap needs 
      if (drsqd >= hrangesqd && !_scatter_overlap(neibpart, drsqd, 0, gather_only())) {
        if (keep_grav) {
          if(gravmask[neibpart.ptype]) {
            directlist.push_back(i);
          } else {
            _NNonInteract++;
          }
        }
      }
      else {
        if (hydromask[neibpart.ptype]){
          culled_neiblist.push_back(i);
        }
        else if (keep_grav) {
          if (gravmask[neibpart.ptype]) {
            smoothgravlist.push_back(i);
          } else {
            _NNonInteract++;
          }
        }
      }
    }
  }


  //===============================================================================================
  //  TrimNeighbour
  /// \brief Only for ComputeHArray
  /// \detail This function create masks for a single neighbour, determining whether
  ///         the neighbour is needed for hydro or gravity. Periodic corrections
  ///         are applied to dr only.  Better to pass in a ParticleType&, but
  //          the index neib is nneded to work out if this is could be a periodic ghost
  //          (maybe have a marked in the particle?)
  //===============================================================================================
  template<class InParticleType, class do_pair_once, class gather_only>
  void TrimNeighbour(const int level[MAX_NPART], const int iorig[MAX_NPART],
		     const FLOAT r[ndim][MAX_NPART], const int npart,
		     const int neib, const Typemask hydromask[MAX_NPART],
		     const FLOAT hrangesqd, const bool keep_grav,
		     FLOAT dr[ndim][MAX_NPART], FLOAT drsqd[MAX_NPART],
		     bool culled[MAX_NPART], bool smoothed_grav[MAX_NPART], bool direct[MAX_NPART])
  {
    bool do_mask[MAX_NPART]; // When using SIMD, mask1 should also mask off
                             // [npart+1,MAX_NPART] for the last SIMD_LENGTH
                             // chunk.  We also want a test for any(mask1[0:7]),
                             // etc. - or will the masked AVX instructions
                             // short-circuit and so this isn't needed?
    bool scatter_mask[MAX_NPART];
    bool hrange_mask[MAX_NPART];
    bool hydro_mask[MAX_NPART];
    bool grav_mask[MAX_NPART];
    bool do_grav;
    bool do_hydro[MAX_NPART];
    bool any, all;
    const GhostNeighbourFinder<ndim> GhostFinder(*_domain);
    Typemask gravmask = _types->gravmask;

    for (int p=0; p<npart; p++) culled[p] = false;
    for (int p=0; p<npart; p++) smoothed_grav[p] = false;
    for (int p=0; p<npart; p++) direct[p] = false;

    for (int p=0; p<npart; p++) do_mask[p] = false;

    // Check the distance. The ones that are not real neighbours
    // are demoted to the direct list
    ParticleType& neibpart = neibdata[neiblist[neib]];

    _first_appearance(level, iorig, npart,
		      neibpart, do_pair_once(),
		      do_mask);
    any = false;
    for (int p=0; p<npart; p++) any |= do_mask[p];
    if (!any) return;

    // Convert from maps of booleans to scalars and arrays.
    do_grav = keep_grav && gravmask[neibpart.ptype];
    for (int p=0; p<npart; p++) do_hydro[p] = hydromask[neibpart.ptype];

    // Compute relative position and distance quantities for pair
    for (int k=0; k<ndim; k++)
      for (int p=0; p<npart; p++) if (do_mask[p]) dr[k][p] = neibpart.r[k] - r[k][p];

    if (neib < _NPeriodicGhosts)
      //  This does not change neibpart.r
      GhostFinder.PartlyApplyPeriodicDistanceCorrection(neibpart.r, dr, do_mask, npart);

    // In the <ndim> needed?
    DotProduct<ndim>(dr, dr, do_mask, npart, drsqd);

    // Record if neighbour is direct-sum or and SPH neighbour.
    // If SPH neighbour, also record max. timestep level for neighbour

    for (int p=0; p<npart; p++) if (do_mask[p]) hrange_mask[p] = drsqd[p] < hrangesqd;
    _scatter_overlap(neibpart.h, drsqd, (FLOAT) 0.0, gather_only(), do_mask, npart, scatter_mask);
    for (int p=0; p<npart; p++) hydro_mask[p] = do_mask[p] &&  (hrange_mask[p] || scatter_mask[p]);
    for (int p=0; p<npart; p++) grav_mask[p]  = do_mask[p] && !(hrange_mask[p] || scatter_mask[p]);

    // Counting of _NNonInteract needs to be added?
    for (int p=0; p<npart; p++) culled[p] = hydro_mask[p] && do_hydro[p] && ;
    if (do_grav) {
      for (int p=0; p<npart; p++) smoothed_grav[p] = hydro_mask[p] && !do_hydro[p];
      for (int p=0; p<npart; p++) direct[p] = grav_mask[p];
    }
  }

  template<class InParticleType>
  bool _first_appearance(const InParticleType& p, const ParticleType& neibpart, _true_type) {
    // Get the neighbour only if this is it's first appearance, which will be true if any of
    // the following criteria are met.
    // (i)   the particle is a ghost
    // (ii)  same i.d. as current active particle
    // (iii) neighbour is on lower timestep level
    // (iv)  the neighbour is on same timestep level as current particle but has larger id value
    return
        (neibpart.flags.is_mirror()) ||
        (p.level > neibpart.level) ||
        (p.level == neibpart.level &&  p.iorig < neibpart.iorig);
  }

  template<class InParticleType>
  bool _first_appearance(const InParticleType& p, const ParticleType& neibpart, _false_type) {
    return true;
  }

  void _first_appearance(const int level[MAX_NPART], const int iorig[MAX_NPART], const int npart,
		         const ParticleType& neibpart, _true_type,
		         bool do_appearance[MAX_NPART]) {
    // Get the neighbour only if this is it's first appearance, which will be true if any of
    // the following criteria are met.
    // (i)   the particle is a ghost
    // (ii)  same i.d. as current active particle
    // (iii) neighbour is on lower timestep level
    // (iv)  the neighbour is on same timestep level as current particle but has larger id value
    for (int p=0; p<npart; p++) {
      do_appearance[p] = (neibpart.flags.is_mirror()) ||
        (level[p] > neibpart.level) ||
        (level[p] == neibpart.level &&  iorig[p] < neibpart.iorig);
    }
  }

  template<class InParticleType>
  void _first_appearance(const int level[MAX_NPART], const int iorig[MAX_NPART], const int npart,
			 const ParticleType& neibpart, _false_type,
			 bool do_appearance[MAX_NPART]) {
    for (int p=0; p<npart; p++) do_appearance[p] = true;
  }

  template<class InParticleType>
  bool _scatter_overlap(const InParticleType& p, double drsqd, double rmax, _false_type) {
    double h = rmax + _kernrange*p.h ;
    return drsqd < h*h;
  }
  template<class InParticleType>
  bool _scatter_overlap(const InParticleType& p, double drsqd, double rmax, _true_type) {
    return false;
  }

  inline void _scatter_overlap(const FLOAT h, const FLOAT drsqd[MAX_NPART], const FLOAT rmax,
			       _false_type,
			       bool mask[MAX_NPART], const int npart,
			       bool scatter_mask[MAX_NPART]) {
    FLOAT h1 = rmax + _kernrange*h;
    FLOAT h1sqd = h1*h1;
    for (int p=0; p<npart; p++) if (mask[p]) scatter_mask[p] = drsqd[p] < h1sqd;
  }
  inline void _scatter_overlap(const FLOAT h, const FLOAT drsqd[MAX_NPART], const FLOAT rmax,
			       _true_type,
			       bool mask[MAX_NPART], const int npart,
			       bool scatter_mask[MAX_NPART]) {
    for (int p=0; p<npart; p++) if (mask[p]) scatter_mask[p] = false;
  }

  template<class InParticleType>
  void _AddParticleAndGhosts
  (const GhostNeighbourFinder<ndim>& GhostFinder,
   InParticleType& part,
   _false_type)
  {
    GhostFinder.ConstructGhostsScatterGather(part, neibdata);
  }

  template<class InParticleType>
  void _AddParticleAndGhosts
  (const GhostNeighbourFinder<ndim>& GhostFinder,
   InParticleType& part,
   _true_type)
  {
    GhostFinder.ConstructGhostsGather(part, neibdata);
  }



public:

  //===============================================================================================
  //  VerifyNeighbourList
  /// \brief    Verify the list of neighbours is complete.
  /// \details  Check that the complete list of neighbours contains all the neighbours needed for
  ///           the particle (i). This check can be for either the "gather" or "all"
  ///           (scatter-gather) of the neighbours.
  //===============================================================================================
  template <class InParticleType>
  void VerifyNeighbourList(int i, int Ntot, const InParticleType& partdata,
                           const string& searchmode) {

    std::vector<int> reducedngb = neib_idx ;

    Typemask types ;
    for (int n=0; n<Ntypes; n++) types[n] = true ;

    __VerifyNeighbourList(i, reducedngb, Ntot, partdata, types, searchmode, "all") ;
  }

  //===============================================================================================
  //  VerifyReducedNeighbourList
  /// \brief    Verify the list of neighbours is complete.
  /// \details  Checks whether the list of neighbours held by the passed list object is complete.
  ///           Currently only works for full scatter/gather neigbours, does not take into account
  ///           the missing neighbours for flux calculations.
  //===============================================================================================
   template <class InParticleType>
   void VerifyReducedNeighbourList(int i, const NeighbourList<ParticleType>& ngbs,
                                   int Ntot, const InParticleType& partdata,
                                   const Typemask& types,
                                   const string& searchmode) {

     // Get the idx of the reduced neighbour list
     std::vector<int> reducedngb(ngbs.size());
     for (unsigned int k=0; k<ngbs.size(); k++)
       reducedngb[k] = neib_idx[ngbs._idx[k]];

     __VerifyNeighbourList(i, reducedngb, Ntot, partdata, types, searchmode, "reduced") ;
   }

private:
   //===============================================================================================
   //  __VerifyNeighbourList
   /// \brief    Verify the list of neighbours is complete.
   /// \details  Compares a list of neighbour indexes (locations in particle array) to the list of
   ///           neighbours found by a brute-force comparison to check whether any neighbours have
   ///           been missed.
   //===============================================================================================
   template <class InParticleType>
   void __VerifyNeighbourList(int i, std::vector<int>& reducedngb,
                              int Ntot, const InParticleType& partdata, const Typemask& types,
                              const string& searchmode, const string& listtype) {
     std::vector<int> truengb ;
     FLOAT drsqd ;
     FLOAT dr[ndim];

     const GhostNeighbourFinder<ndim> GhostFinder(*_domain);

     // Compute the complete (true) list of neighbours
     if (searchmode == "gather") {
       for (int j=0; j<Ntot; j++) {
         if (not types[partdata[j].ptype]) continue ;

         for (int k=0; k<ndim; k++) dr[k] = partdata[j].r[k] - partdata[i].r[k];
         ///GhostFinder.NearestPeriodicVector(dr);
         drsqd = DotProduct(dr,dr,ndim);
         if (partdata[j].flags.is_dead()) continue;
         if (drsqd <= partdata[i].hrangesqd) truengb.push_back(j);
       }
     }
     else if (searchmode == "all") {
       for (int j=0; j<Ntot; j++) {
         if (not types[partdata[j].ptype]) continue ;

         for (int k=0; k<ndim; k++) dr[k] = partdata[j].r[k] - partdata[i].r[k];
         GhostFinder.NearestPeriodicVector(dr);
         drsqd = DotProduct(dr,dr,ndim);
         if (partdata[j].flags.is_dead()) continue;
         if (drsqd <= partdata[i].hrangesqd ||
             drsqd <= partdata[j].hrangesqd) truengb.push_back(j);
       }
     }

     // Now compare each given neighbour with true neighbour list for validation
     int invalid_flag = 0;
     unsigned int j;
     for (j=0; j<truengb.size(); j++) {
       int count = 0;
       for (unsigned int k=0; k<reducedngb.size(); k++) {
         if (reducedngb[k] == truengb[j])
           count++;
       }

       // If the true neighbour is not in the list, or included multiple times,
       // then output to screen and terminate program
       if (count != 1) {
         for (int k=0; k<ndim; k++) dr[k] = partdata[truengb[j]].r[k] - partdata[i].r[k];
         GhostFinder.NearestPeriodicVector(dr);
         drsqd = DotProduct(dr,dr,ndim);
         cout << "Could not find neighbour " << j << "   " << truengb[j] << "     " << i
             << "      " << sqrt(drsqd)/sqrt(partdata[i].hrangesqd) << "     "
             << sqrt(drsqd)/sqrt(partdata[j].hrangesqd) << "    "
             << partdata[truengb[j]].r[0] << "   type : "
             << partdata[truengb[j]].ptype << endl;
         invalid_flag++;
       }

     }
     // If the true neighbour is not in the list, or included multiple times,
     // then output to screen and terminate program
     if (invalid_flag) {
       cout << "Problem with neighbour lists (" << listtype << ") : " << i << "  " << j << "   "
           <<  invalid_flag << "    " << partdata[i].r[0] << "   " << partdata[i].h << endl
           << "Nneib : " << neib_idx.size() << "   Ntrueneib : " << truengb.size()
           << "    searchmode : " << searchmode << endl;
       InsertionSort(reducedngb.size(), &reducedngb[0]);
       PrintArray("neiblist     : ",reducedngb.size(), &reducedngb[0]);
       PrintArray("trueneiblist : ",truengb.size() , &truengb[0]);
       string message = "Problem with neighbour lists in tree search";
       ExceptionHandler::getIstance().raise(message);
     }
   }
};



#endif /* NEIGHBOURMANAGER_H_ */
