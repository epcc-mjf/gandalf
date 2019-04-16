/*
 * TreeCell.h
 *
 *  Created on: 5 Dec 2016
 *      Author: rosotti
 */

#ifndef SRC_HEADERS_TREECELL_H_
#define SRC_HEADERS_TREECELL_H_

#include <type_traits>
#define CACHE_LINE 0x40
#define AVX512_LENGTH 0x40
#define AVX_LENGTH 0x20

//=================================================================================================
//  Struct TreeCellBase
/// Base tree cell data structure which contains all data elements common to all trees.
//=================================================================================================
template <int ndim>
struct TreeCellBase {
  // MJF All this may not be OK if int becomes long for larger numbers of cells.
  // MJF This is all only correct for FLOAT=double.
  alignas(CACHE_LINE) Box<ndim> bb ;   ///< Bounding box
  // For ndim=3 and FLOAT=double, this should be cacheline + 2*3*8 = 48 bytes
  // (no padding at the end of BoxOverlap).
  int cnext;                           ///< i.d. of next cell if not opened
  int copen;                           ///< i.d. of first child cell
  int N;                               ///< No. of particles in cell
  int id;                              ///< Cell id
  // 1 cache line.  id is put here to get alignment
  alignas(CACHE_LINE) Box<ndim> vbox ; ///< Velocity space bounding box
  int parent;                          ///< Id of the cell's parent
  int level;                           ///< Level of cell on tree
  int Nactive;                         ///< No. of active particles in cell
  float maxsound;                      ///< Maximum sound speed inside the cell
  // 1 cache line.
  alignas(CACHE_LINE) FLOAT r[ndim];   ///< Position of cell COM
  FLOAT cdistsqd;                      ///< Minimum distance to use COM values
  FLOAT mac;                           ///< Multipole-opening criterion value
  FLOAT m;                             ///< Mass contained in cell
  FLOAT rmax;                          ///< Radius of bounding sphere
  union {
    FLOAT amin;                        ///< Minimum grav accel of particles in the cell
    FLOAT macfactor;                   ///< Potential based accuracy factor.
  } ;
  // 1 cache line
  alignas(CACHE_LINE) Box<ndim> hbox;  ///< Bounding box for smoothing volume
  // For ndim=3 and FLOAT=double, this should be cacheline + 2*3*8 = 48 bytes
  int ifirst;                          ///< i.d. of first particle in cell
  int ilast;                           ///< i.d. of last particle in cell
  FLOAT hmax;                          ///< Maximum smoothing length inside cell
  // 1 cache line.
  alignas(CACHE_LINE) FLOAT q[5];      ///< Quadrupole moment tensor
#ifdef MPI_PARALLEL
  double worktot;                      ///< Total work in cell
#endif
  // Align c1 and c2 of KDTreeCell to a cache line (see KDTree.h) => 6 cache
  // lines per TreeCell.

  void ComputeCellCentre(FLOAT rc[ndim]) const {
    for (int k=0; k<ndim; k++) rc[k] = rcell(k);
  }

  FLOAT rcell(int k) const {
    return  (FLOAT) 0.5*(bb.min[k] + bb.max[k]);
  }


};

//=================================================================================================
//  Struct MultipoleMoment
/// Structure to hold the multipole moment data.
//=================================================================================================
template<int ndim>
struct MultipoleMoment {
  MultipoleMoment()
  {
    for (int k=0; k<ndim; k++) r[k] = 0 ;
    for (int k=0; k<5; k++) q[k] = 0 ;
    m = 0;
    id = 0;
  }

  explicit MultipoleMoment(const TreeCellBase<ndim>& cell)
  {
    for (int k=0; k<ndim; k++) r[k] = cell.r[k] ;
    for (int k=0; k<5; k++) q[k] = cell.q[k] ;
    m = cell.m;
    id = cell.id;
  }

  FLOAT r[ndim];                       ///< Position of cell COM
  FLOAT m;                             ///< Mass contained in cell
  FLOAT q[5];                          ///< Quadrupole moment tensor
  int id ;
};




#endif /* SRC_HEADERS_TREECELL_H_ */
