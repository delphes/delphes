#ifndef __FASTJET_GRIDJETPLUGIN_HH__
#define __FASTJET_GRIDJETPLUGIN_HH__

//FJSTARTHEADER
// $Id: GridJetPlugin.hh 2267 2011-06-20 15:10:23Z salam $
//
// Copyright (c) 2005-2014, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of FastJet.
//
//  FastJet is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  The algorithms that underlie FastJet have required considerable
//  development. They are described in the original FastJet paper,
//  hep-ph/0512210 and in the manual, arXiv:1111.6097. If you use
//  FastJet as part of work towards a scientific publication, please
//  quote the version you use and include a citation to the manual and
//  optionally also to hep-ph/0512210.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//FJENDHEADER


#include "fastjet/JetDefinition.hh"

// makes it easy to switch back and forth between use of
// RectangularGrid or not; this got enabled in FJ3.1
#define FASTJET_GRIDJET_USEFJGRID

#ifdef FASTJET_GRIDJET_USEFJGRID
#include "fastjet/RectangularGrid.hh"
#endif

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

// forward declaration to reduce includes
class ClusterSequence;

//----------------------------------------------------------------------
//
/// @ingroup plugins
/// \class GridJetPlugin
/// plugin for fastjet (v3.0 upwards) that clusters particles such
/// that all particles in a given cell of a rectangular rapidity-phi
/// grid end up in a common "jet".
///
/// This is not intended for use as a regular jet clustering algorithm, 
/// but is rather provided for comparison purposes with the 
/// GridMedianBackgroundEstimator (which is even faster).
class GridJetPlugin : public JetDefinition::Plugin 
#ifdef FASTJET_GRIDJET_USEFJGRID
                                      , RectangularGrid
#endif 
{
public:
  /// Basic constructor for the GridJetPlugin Plugin class.
  ///
  /// \param ymax           The maximal rapidity extent of the grid
  /// \param requested_grid_spacing The requested grid spacing
  /// \param post_jet_def   if present, and not == JetDefinition()
  ///                       (which has undefined_jet_algorithm), then
  ///                       run the post_jet_def on the result of the grid
  ///                       clustering.
  GridJetPlugin (double ymax, double requested_grid_spacing,
		 const JetDefinition & post_jet_def = JetDefinition());

#ifdef FASTJET_GRIDJET_USEFJGRID
  /// Constructor for the GridJetPlugin Plugin class that allows
  /// full control over the underlying grid. New in FastJet 3.1.
  ///
  /// \param grid           The maximal rapidity extent of the grid
  /// \param post_jet_def   if present, and not == JetDefinition()
  ///                       (which has undefined_jet_algorithm), then
  ///                       run the post_jet_def on the result of the grid
  ///                       clustering.
  GridJetPlugin (const RectangularGrid & grid,
		 const JetDefinition & post_jet_def = JetDefinition());
#endif // FASTJET_GRIDJET_USEFJGRID

  

  // /// copy constructor
  // GridJetPlugin (const GridJetPlugin & plugin) {
  //   *this = plugin;
  // }

  // the things that are required by base class
  virtual std::string description () const;
  virtual void run_clustering(ClusterSequence &) const;

  /// This returns the sqrt(dphi*dy/pi) -- i.e. the radius that for a
  /// circular jet would give the same area.
  virtual double R() const;

  // As of FastJet 3.1 the following functions become available through
  // the underlying RectangularGrid class.
#ifndef FASTJET_GRIDJET_USEFJGRID
  /// returns the actual rapidity spacing of the grid
  double drap()   const {return _dy;}
  /// returns the actual phi spacing of the grid
  double dphi() const {return _dphi;}
  /// returns the minimum rapidity of the grid
  double rapmin() const {return _ymin;}
  /// returns the maximum rapidity of the grid
  double rapmax() const {return _ymax;}
#endif

private:

#ifndef FASTJET_GRIDJET_USEFJGRID
  void setup_grid();

  int n_tiles() const {return _ntotal;}
  int n_good_tiles() const {return _ntotal;}

  int tile_index(const PseudoJet & p) const;
  bool tile_is_good(int /* itile */) const {return true;}

  double _ymin, _ymax, _dy, _dphi, _requested_grid_spacing;
  int _ny, _nphi, _ntotal;
#endif

  JetDefinition _post_jet_def;

};

FASTJET_END_NAMESPACE        // defined in fastjet/internal/base.hh

#endif // __FASTJET_GRIDJETPLUGIN_HH__

