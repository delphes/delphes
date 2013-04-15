//STARTHEADER
// $Id: TrackJetPlugin.hh 2758 2011-11-24 08:31:58Z soyez $
//
// Copyright (c) 2007-2011, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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

#ifndef __TRACKJETPLUGIN_HH__
#define __TRACKJETPLUGIN_HH__

#include "fastjet/JetDefinition.hh"

// questionable whether this should be in fastjet namespace or not...
FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

// another forward declaration to reduce includes
class PseudoJet;

//----------------------------------------------------------------------

/// @ingroup plugins
/// \class TrackJetPlugin
/// Implementation of the TrackJet algorithm (plugin for fastjet v2.4 upwards)
//
class TrackJetPlugin : public JetDefinition::Plugin {
public:
  /// Main constructor for the TrackJet Plugin class.  
  ///
  /// The argument is an initialised list of jet algorithms
  /// \param _radius  the distance at which point a particle is no longer
  ///                 recombied into the jet
  /// \param jet_recombination_scheme  the recombination scheme used to 
  ///                                  sum the 4-vecors inside the jet
  /// \param track_recombination_scheme  the recombination scheme used to 
  ///                                    sum the 4-vecors when accumulating
  ///                                    track into a the jet
  /// Both recombiners are defaulted to pt_scheme recomb as for the Rivet
  /// implementation.
  TrackJetPlugin (double radius, 
		  RecombinationScheme jet_recombination_scheme=pt_scheme, 
		  RecombinationScheme track_recombination_scheme=pt_scheme){
    _radius  = radius;
    _radius2 = radius*radius;
    _jet_recombiner = JetDefinition::DefaultRecombiner(jet_recombination_scheme);
    _track_recombiner = JetDefinition::DefaultRecombiner(track_recombination_scheme);
  }

  /// copy constructor
  TrackJetPlugin (const TrackJetPlugin & plugin) {
    *this = plugin;
  }

  // the things that are required by base class
  virtual std::string description () const;
  virtual void run_clustering(ClusterSequence &) const;

  /// the plugin mechanism's standard way of accessing the jet radius
  /// here we return the R of the last alg in the list
  virtual double R() const {return _radius;}

private:
  double _radius, _radius2;

  JetDefinition::DefaultRecombiner _jet_recombiner;
  JetDefinition::DefaultRecombiner _track_recombiner;

  static bool _first_time;

  /// print a banner for reference to the 3rd-party code
  void _print_banner(std::ostream *ostr) const;
};

FASTJET_END_NAMESPACE        // defined in fastjet/internal/base.hh

#endif // __TRACKJETPLUGIN_HH__

