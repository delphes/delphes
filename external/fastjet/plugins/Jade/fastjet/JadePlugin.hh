#ifndef __JADEPLUGIN_HH__
#define __JADEPLUGIN_HH__

//STARTHEADER
// $Id: JadePlugin.hh 2577 2011-09-13 15:11:38Z salam $
//
// Copyright (c) 2009, Matteo Cacciari, Gavin Salam and Gregory Soyez
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
/// \class JadePlugin
/// Implementation of the e+e- Jade algorithm (plugin for fastjet v2.4 upwards)
///
/// JadePlugin is a plugin for fastjet (v2.4 upwards)
/// It implements the JADE algorithm, which is an e+e- sequential
/// recombination algorithm with interparticle distance
///
///   dij = 2 E_i E_j (1 - cos theta_ij)
///
/// or equivalently
///
///   yij = dij/E_{vis}^2                
///
/// This corresponds to the distance measured used in 
///
///   "Experimental Investigation of the Energy Dependence of the Strong Coupling Strength."
///   JADE Collaboration (S. Bethke et al.)
///   Phys.Lett.B213:235,1988
///
/// The JADE article carries out particle recombinations in the
/// E-scheme (4-vector recombination), which is the default procedure for this 
/// plugin.
///
/// NOTE: other widely used schemes include E0, P, P0; however they also 
///       involve modifications to the distance measure. Be sure of
///       what you're doing before running a JADE type algorithm.
///
/// To access the jets with a given ycut value (clustering stops once
/// all yij > ycut), use
///
///   vector<PseudoJet> jets = cluster_sequence.exclusive_jets_ycut(ycut);
///
/// and related routines.
class JadePlugin : public JetDefinition::Plugin {
public:
  /// Main constructor for the Jade Plugin class.  
  JadePlugin (){}

  /// copy constructor
  JadePlugin (const JadePlugin & plugin) {
    *this = plugin;
  }

  // the things that are required by base class
  virtual std::string description () const;
  virtual void run_clustering(ClusterSequence &) const;

  /// the plugin mechanism's standard way of accessing the jet radius.
  /// This must be set to return something sensible, even if R
  /// does not make sense for this algorithm!
  virtual double R() const {return 1.0;}

  /// avoid the warning whenever the user requests "exclusive" jets
  /// from the cluster sequence
  virtual bool exclusive_sequence_meaningful() const {return true;}

private:

};

FASTJET_END_NAMESPACE        // defined in fastjet/internal/base.hh

#endif // __JADEPLUGIN_HH__

