#ifndef __GRIDJETPLUGIN_HH__
#define __GRIDJETPLUGIN_HH__

//STARTHEADER
// $Id: GridJetPlugin.hh 2267 2011-06-20 15:10:23Z salam $
//
// Copyright (c) 2011, Matteo Cacciari, Gavin Salam and Gregory Soyez
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
//  development and are described in hep-ph/0512210. If you use
//  FastJet as part of work towards a scientific publication, please
//  include a citation to the FastJet paper.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//ENDHEADER


#include "fastjet/JetDefinition.hh"

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
class GridJetPlugin : public JetDefinition::Plugin {
public:
  /// Main constructor for the GridJetPlugin Plugin class.
  ///
  /// \param ymax           The maximal rapidity extent of the grid
  /// \param requested_grid_spacing The requested grid spacing
  /// \param post_jet_def   if present, and not == JetDefinition()
  ///                       (which has undefined_jet_algorithm), then
  ///                       run the post_jet_def on the result of the grid
  ///                       clustering.
  GridJetPlugin (double ymax, double requested_grid_spacing,
		 const JetDefinition & post_jet_def = JetDefinition());

  /// copy constructor
  GridJetPlugin (const GridJetPlugin & plugin) {
    *this = plugin;
  }

  // the things that are required by base class
  virtual std::string description () const;
  virtual void run_clustering(ClusterSequence &) const;

  /// This returns the sqrt(dphi*dy/pi) -- i.e. the radius that for a
  /// circular jet would give the same area.
  virtual double R() const;

private:

  void setup_grid();

  int igrid(const PseudoJet & p) const;

  double _ymin, _ymax, _dy, _dphi, _requested_grid_spacing;
  int _ny, _nphi, _ntotal;

  JetDefinition _post_jet_def;

};

FASTJET_END_NAMESPACE        // defined in fastjet/internal/base.hh

#endif // __GRIDJETPLUGIN_HH__

