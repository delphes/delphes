#ifndef __EECAMBRIDGEPLUGIN_HH__
#define __EECAMBRIDGEPLUGIN_HH__

//STARTHEADER
// $Id: EECambridgePlugin.hh 2692 2011-11-14 16:27:44Z soyez $
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
/// \class EECambridgePlugin
/// Implementation of the e+e- Cambridge algorithm (plugin for fastjet v2.4 upwards)
///
/// EECambridgePlugin is a plugin for fastjet (v2.4 upwards)
/// It implements the Cambridge algorithm, as defined in 
/// 
/// Better jet clustering algorithms
/// Yuri Dokshitzer, Garth Leder, Stefano Moretti,  Bryan Webber 
/// JHEP 9708 (1997) 001
/// http://www-spires.slac.stanford.edu/spires/find/hep/www?rawcmd=FIND+j+JHEPA%2C9708%2C001
///
/// On construction one must supply a ycut value.
///
/// To get the jets at the end call ClusterSequence::inclusive_jets();
class EECambridgePlugin : public JetDefinition::Plugin {
public:
  /// Main constructor for the EECambridge Plugin class.  
  /// It takes the dimensionless parameter ycut (the Q value for normalisation
  /// of the kt-distances is taken from the sum of all particle energies).
  EECambridgePlugin (double ycut_in) : _ycut(ycut_in) {}

  /// copy constructor
  EECambridgePlugin (const EECambridgePlugin & plugin) {
    *this = plugin;
  }

  // the things that are required by base class
  virtual std::string description () const;
  virtual void run_clustering(ClusterSequence &) const;

  double ycut() const {return _ycut;}

  /// the plugin mechanism's standard way of accessing the jet radius.
  /// This must be set to return something sensible, even if R
  /// does not make sense for this algorithm!
  virtual double R() const {return 1.0;}

  /// avoid the warning whenever the user requests "exclusive" jets
  /// from the cluster sequence
  virtual bool exclusive_sequence_meaningful() const {return true;}

private:
  double _ycut;
};

FASTJET_END_NAMESPACE        // defined in fastjet/internal/base.hh

#endif // __EECAMBRIDGEPLUGIN_HH__

