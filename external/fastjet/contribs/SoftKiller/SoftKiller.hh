// $Id$
//
// Copyright (c) 2014-, Matteo Cacciari, Gavin. P. Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------

#ifndef __FASTJET_CONTRIB_SOFTKILLER_HH__
#define __FASTJET_CONTRIB_SOFTKILLER_HH__

#include <fastjet/internal/base.hh>
#include <fastjet/Selector.hh>
#include "fastjet/config.h"

#if FASTJET_VERSION_NUMBER >= 30100
#define FJCONTRIB_SOFTKILLER_USEFJGRID
#endif

#ifdef FJCONTRIB_SOFTKILLER_USEFJGRID
#include "fastjet/RectangularGrid.hh"
#endif

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


namespace contrib{

//------------------------------------------------------------------------
/// \class SoftKiller
/// Progressively kills soft particles in order of increasing pt until
/// half of the event is empty.
///
/// More precisely, the event up to |rap|=rapmax is split into grid
/// cells of size "cell_size". Soft particles are removed from the
/// event in order of increasing pt until half the cells are empty.
/// By default, that same pt cut is applied to _all_ the particles in
/// the event, not just those in the grid region.
///
/// If a sifter is provided, then only particles that pass the selector
/// are put onto the grid, and the pt threshold is then applied just 
/// to those particles, while all other particles are left untouched.
///
/// This is convenient, for example, if you want to apply SK just to
/// neutral particles.
#ifdef FJCONTRIB_SOFTKILLER_USEFJGRID
class SoftKiller : public RectangularGrid {
#else 
class SoftKiller {
#endif
public:
  /// ctor with simple initialisation
  ///  \param rapmax     the maximal absolute rapidity extent of the grid
  ///  \param tile_size  the requested grid spacing (equivalently, tile size)
  ///  \param sifter     when provided, the soft killer is applied
  ///                    only to particles that pass the sifter (the
  ///                    others are kept untouched)
  SoftKiller(double rapmax, double tile_size,
	     Selector sifter = Selector());

  /// ctor with more control over initialisation
  ///  \param rapmin     the minimum rapidity extent of the grid
  ///  \param rapmax     the maximum rapidity extent of the grid
  ///  \param drap       the grid spacing in rapidity
  ///  \param dphi       the grid spacing in azimuth
  ///  \param sifter     when provided, the soft killer is applied
  ///                    only to particles that pass the sifter (the
  ///                    others are kept untouched)
  SoftKiller(double rapmin, double rapmax, double drap, double dphi,
	     Selector sifter = Selector());

#ifdef FJCONTRIB_SOFTKILLER_USEFJGRID
  /// constructor that takes RectangularGrid object to flexibly
  /// specify the details of the grid (works only with FJ3.1)
  SoftKiller(const RectangularGrid & grid, Selector sifter = Selector());
#endif 

  /// dummy ctor (will give an unusable SoftKiller)
  SoftKiller();

  /// returns description of the soft killer
  std::string description() const;
  
  /// applies the soft killer to a given event. The result is the
  /// "reduced" event.
  std::vector<PseudoJet> operator()(const std::vector<PseudoJet> & event) const{
    return result(event);
  }
  
  /// similarly to Transformers in FastJet, introduce a 'result'
  /// method equivalent to the () operator.
  //std::vector<PseudoJet> result(const std::vector<PseudoJet> & event) const;
  std::vector<PseudoJet> result(const std::vector<PseudoJet> & event) const {
    double pt_threshold;
    std::vector<PseudoJet> reduced_event;
    apply(event, reduced_event, pt_threshold);
    return reduced_event;
  }

  /// alternative invocation that puts the resulting SK "reduced" event
  /// into the reduced_event vector, and also provides information about the
  /// pt threshold that was applied.
  ///
  /// The event and reduced_event must _not_ be the same variable.
  void apply(const std::vector<PseudoJet> & event, 
             std::vector<PseudoJet> & reduced_event,
             double & pt_threshold) const;

private:

#ifndef FJCONTRIB_SOFTKILLER_USEFJGRID
  /// initial setup of the grid 
  void _setup_grid();

  // retrieve the grid cell index for a given PseudoJet
  inline int tile_index(const PseudoJet & p) const;

  inline int n_tiles() const {return _ntotal;}
  inline int n_good_tiles() const {return n_tiles();}
  inline bool tile_is_good(int itile) const {return true;}
  
  inline bool all_tiles_equal_area() const {return true;}

  // ctor arguments
  double _ymax, _ymin;  ///< maximal and minimal rapidity coverage of the grid
  //double _cell_size;    ///< grid cell size
  double _requested_drap; ///< requested rapidity spacing
  double _requested_dphi; ///< requested phi spacing

  // information about the grid
  double _dy, _dphi, _cell_area, _inverse_dy, _inverse_dphi;
  int _ny, _nphi, _ntotal;

#endif // FJCONTRIB_SOFTKILLER_USEFJGRID

  Selector _sifter;     ///< optional sifter

};

} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_SOFTKILLER_HH__
