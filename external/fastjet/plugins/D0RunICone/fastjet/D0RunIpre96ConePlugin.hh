#ifndef __D0RUNIPRE96CONEPLUGIN_HH__
#define __D0RUNIPRE96CONEPLUGIN_HH__

//FJSTARTHEADER
// $Id: D0RunIpre96ConePlugin.hh 1778 2010-10-25 10:02:58Z soyez $
//
// Copyright (c) 2009-2014, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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

#include "fastjet/internal/base.hh"     // namespace macros (include explicitly to help Doxygen)
#include "D0RunIBaseConePlugin.hh"

// questionable whether this should be in fastjet namespace or not...

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//----------------------------------------------------------------------
//
/// @ingroup plugins
/// \class D0RunIpre96ConePlugin
///
/// A plugin for FastJet (v3.0 or later) that provides an interface to
/// the pre 1996 D0 version of Run-I cone algorithm
///
/// The D0 code has been obtained from Lars Sonnenschein's web-space
/// http://www-d0.fnal.gov/~sonne/D0RunIcone.tgz
///
/// The version of the D0 Run I code distributed
/// here has been modified by the FastJet authors, so as to provide
/// access to the contents of the jets (as is necessary for the
/// plugin). This does not modify the results of the clustering.
///
/// The difference between this algorithm and the post-1996 version
/// relates to the way the final jet momenta are calculated. Details
/// are to be found in FERMILAB-PUB-97-242-E.
//
//----------------------------------------------------------------------
class D0RunIpre96ConePlugin : public D0RunIBaseConePlugin {
public:
  /// The D0RunIpre96ConePlugin constructor, which sets the "free" parameters of the
  /// algorithm:
  ///
  ///  \param CONErad is the cone radius
  ///
  ///  \param JETmne is a minimum ET requirement on every iteration
  ///   (jet dropped if Et < JETmne * Et_min_ratio ).
  ///   The value that has been used by D0 for JETmne: 8 GeV 
  ///   (and Et_min_ratio is 0.5)
  ///
  ///  \param SPlifr is the shared Et fraction splitting threshold, and
  ///   a value of 0.5 was usually used by D0
  ///
  /// The remaining parameters of the algorithm are not to be modified if the algorithm
  /// is to correspond to the one actually used by D0.
  ///
  ///
  D0RunIpre96ConePlugin (double CONErad_in, double JETmne_in , double SPLifr_in = _DEFAULT_SPLifr)
    : D0RunIBaseConePlugin(CONErad_in, JETmne_in , SPLifr_in){}

  // the things that are required by base class
  virtual std::string description () const;

  // the part that really does the clustering
  virtual void run_clustering(ClusterSequence &) const;

private:
  static bool _first_time;

  /// print a banner for reference to the 3rd-party code
  void _print_banner(std::ostream *ostr) const;
};


FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh

#endif // __D0RUNIPRE96CONEPLUGIN_HH__
