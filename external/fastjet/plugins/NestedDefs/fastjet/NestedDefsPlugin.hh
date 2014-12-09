//FJSTARTHEADER
// $Id: NestedDefsPlugin.hh 3433 2014-07-23 08:17:03Z salam $
//
// Copyright (c) 2007-2014, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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

#ifndef __NESTEDALGSPLUGIN_HH__
#define __NESTEDALGSPLUGIN_HH__

#include "fastjet/JetDefinition.hh"
#include <list>
#include <memory>
#include <cmath>

// questionable whether this should be in fastjet namespace or not...
FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

// another forward declaration to reduce includes
class PseudoJet;

//----------------------------------------------------------------------
//
/// @ingroup plugins
/// \class NestedDefsPlugin
/// Plugin to run multiple jet definitions successively (plugin for fastjet v2.4 upwards)
///
/// NestedAglsPlugin is a plugin for fastjet (v2.4 upwards) that, given
/// a list of jet definitions, performs the clustering by feeding the 
/// particles to the first algorithm and then, successively feeding the 
/// output to the next algorithm in the list.
//
class NestedDefsPlugin : public JetDefinition::Plugin {
public:
  /// Main constructor for the NestedDefs Plugin class.  
  ///
  /// The argument is an initialised list of jet algorithms
  NestedDefsPlugin (std::list<JetDefinition> &defs) :
    _defs(defs){}

  /// copy constructor
  NestedDefsPlugin (const NestedDefsPlugin & plugin) {
    *this = plugin;
  }

  // the things that are required by base class
  virtual std::string description () const;
  virtual void run_clustering(ClusterSequence &) const;
  /// the plugin mechanism's standard way of accessing the jet radius
  /// here we return the R of the last alg in the list
  virtual double R() const {return _defs.rbegin()->R();}

private:
  std::list<JetDefinition> _defs;
};

FASTJET_END_NAMESPACE        // defined in fastjet/internal/base.hh

#endif // __SISCONEPLUGIN_HH__

