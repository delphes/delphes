//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//  $Id: ExtraRecombiners.hh 828 2015-07-20 14:52:06Z jthaler $
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

#ifndef __FASTJET_CONTRIB_WINNERTAKEALLRECOMBINER_HH__
#define __FASTJET_CONTRIB_WINNERTAKEALLRECOMBINER_HH__

#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"

#include <cmath>
#include <vector>
#include <list>
#include <limits>
#include <stdio.h>
#include <string.h>
#include <errno.h>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {

///------------------------------------------------------------------------
/// \class GeneralEtSchemeRecombiner
/// \brief Recombination scheme with generalized Et weighting
///
/// GeneralEtSchemeRecombiner defines a new recombination scheme by inheriting from JetDefinition::Recombiner.
/// This scheme compares the pT of two input particles, and then combines them into a particle with
/// a pT equal to the sum of the two particle pTs and a direction (in rapidity/phi) weighted by the respective momenta of the
/// particle. The weighting is dependent on the power delta. For delta = infinity, this should return the same result as the
/// WinnerTakeAllRecombiner.
///------------------------------------------------------------------------
class GeneralEtSchemeRecombiner : public fastjet::JetDefinition::Recombiner {
public:

   /// Constructor takes delta weighting
   /// (delta = 1.0 for Et-scheme, delta = infinity for winner-take-all scheme)
   GeneralEtSchemeRecombiner(double delta) : _delta(delta) {}
  
   /// Description
   virtual std::string description() const;
  
   /// Recombine pa and pb and put result into pab
   virtual void recombine(const fastjet::PseudoJet & pa,
                         const fastjet::PseudoJet & pb, 
                         fastjet::PseudoJet & pab) const;

private:
  double _delta;  ///< Weighting exponent
};

///------------------------------------------------------------------------
/// \class WinnerTakeAllRecombiner
/// \brief Recombination scheme with winner-take-all weighting
///
/// WinnerTakeAllRecombiner defines a new recombination scheme by inheriting from JetDefinition::Recombiner.
/// This scheme compares the pT of two input particles, and then combines them into a particle with
/// a pT equal to the sum of the two particle pTs and a direction (in rapidity/phi) identical to that of the harder
/// particle. This creates a jet with an axis guaranteed to align with a particle in the event.
///------------------------------------------------------------------------
class WinnerTakeAllRecombiner : public fastjet::JetDefinition::Recombiner {
public:
   
	/// Constructor to choose value of alpha (defaulted to 1 for normal pT sum)
   WinnerTakeAllRecombiner(double alpha = 1.0) : _alpha(alpha) {}

   /// Description
   virtual std::string description() const;
   
   /// recombine pa and pb and put result into pab
   virtual void recombine(const fastjet::PseudoJet & pa,
                          const fastjet::PseudoJet & pb,
                          fastjet::PseudoJet & pab) const;

private:
	double _alpha; //power of (pt/E) term when recombining particles
};

} //namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_WINNERTAKEALLRECOMBINER_HH__
