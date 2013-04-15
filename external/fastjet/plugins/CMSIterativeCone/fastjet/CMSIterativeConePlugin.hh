//STARTHEADER
// $Id: CMSIterativeConePlugin.hh 1508 2009-04-10 22:46:49Z soyez $
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

#ifndef __CMSITERATIVECONEPLUGIN_HH__
#define __CMSITERATIVECONEPLUGIN_HH__

#include "fastjet/JetDefinition.hh"

// questionable whether this should be in fastjet namespace or not...
FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

// forward declaration to reduce includes
class PseudoJet;

//----------------------------------------------------------------------
//

/// @ingroup plugins
/// \class CMSIterativeConePlugin
/// Implementation of the CMS Iterative Cone (plugin for fastjet v2.4 upwards)
class CMSIterativeConePlugin : public JetDefinition::Plugin {
public:
  /// Main constructor for the CMSIterativeCone Plugin class.  
  ///
  /// The arguments are
  ///   ConeRadius      the radius of the cone
  ///   SeedThreshold   a threshold for the seeds to iterate from
  ///
  /// NOTE: to be more coherent with all other fastjet plugins,
  ///       we've put the radius before the seed threshold. 
  ///       CMS does the opposite.
  ///       In this way, we also put a default value of 0 for the 
  ///       seed threshold.
  CMSIterativeConePlugin (double ConeRadius, double SeedThreshold=1.0) :
    theConeRadius(ConeRadius), theSeedThreshold(SeedThreshold){}

  /// copy constructor
  CMSIterativeConePlugin (const CMSIterativeConePlugin & plugin) {
    *this = plugin;
  }

  // the things that are required by base class
  virtual std::string description () const;
  virtual void run_clustering(ClusterSequence &) const;

  /// the plugin mechanism's standard way of accessing the jet radius
  /// here we return the R of the last alg in the list
  virtual double R() const {return theConeRadius;}

  /// get the seed threshold
  virtual double seed_threshold() const {return theSeedThreshold;}

private:
  double theConeRadius;     ///< cone radius
  double theSeedThreshold;  ///< seed threshold

  static bool _first_time;

  /// print a banner for reference to the 3rd-party code
  void _print_banner(std::ostream *ostr) const;
};

FASTJET_END_NAMESPACE        // defined in fastjet/internal/base.hh

#endif // __CMSITERATIVECONEPLUGIN_HH__

