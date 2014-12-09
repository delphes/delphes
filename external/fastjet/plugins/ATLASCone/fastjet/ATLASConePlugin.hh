//FJSTARTHEADER
// $Id: ATLASConePlugin.hh 3433 2014-07-23 08:17:03Z salam $
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

// Note on the implementation:
//   this implementation of the ATLAS Cone is based on the SpartyJet
//   v2.20.0 implementation. See README for details

#ifndef __ATLASCONEPLUGIN_HH__
#define __ATLASCONEPLUGIN_HH__

#include "fastjet/JetDefinition.hh"

// questionable whether this should be in fastjet namespace or not...
FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

// forward declaration to reduce includes
class PseudoJet;

//----------------------------------------------------------------------
//
/// @ingroup plugins
/// \class ATLASConePlugin
/// Implementation of the ATLAS Cone (plugin for fastjet v2.4 upwards)
class ATLASConePlugin : public JetDefinition::Plugin {
public:
  /// Main constructor for the ATLASCone Plugin class.
  ///
  /// Apparently the default parameters in the ATLAS software are the
  /// ones used here. SpartyJet uses a radius of 0.7, a seed threshold
  /// of 1 GeV and an overlap threshold of 0.75
  /// For the ATLAS SW defaults, see
  ///   http://atlas-sw.cern.ch/cgi-bin/viewcvs-atlas.cgi/groups/JetRoutines/SpartyJet/atlas/
  /// in the JetdoneFinderTools.cxx (rev1.1) and JetSplitMergeTool.cxx (rev1.1)
  /// For SpartyJet, see atlas/ConeFinderTool.h
  ///
  /// Finally, to agree with FastJet standards, we do not specify a default R,
  /// that in the ATLAS code is 0.7
  ATLASConePlugin (double radius, double seedPt_in=2.0, double f_in=0.5)
    : _radius(radius), _seedPt(seedPt_in), _f(f_in){}

  /// copy constructor
  ATLASConePlugin (const ATLASConePlugin & plugin) {
    *this = plugin;
  }

  // the things that are required by base class
  virtual std::string description () const;
  virtual void run_clustering(ClusterSequence &) const;

  /// the plugin mechanism's standard way of accessing the jet radius
  /// here we return the R of the last alg in the list
  virtual double R() const {return _radius;}

  // access to the other parameters
  /// seed threshold
  double seedPt() const {return _seedPt;}

  /// split-merge overlap threshold
  double f() const {return _f;}

private:

  double _radius;   ///< the cone radius
  double _seedPt;   ///< the pt seed threshold used in stable-cone search
  double _f;        ///< the overlap thresholod used in the split-merge

  static bool _first_time;

  /// print a banner for reference to the 3rd-party code
  void _print_banner(std::ostream *ostr) const;
};

FASTJET_END_NAMESPACE        // defined in fastjet/internal/base.hh

#endif // __ATLASCONEPLUGIN_HH__

