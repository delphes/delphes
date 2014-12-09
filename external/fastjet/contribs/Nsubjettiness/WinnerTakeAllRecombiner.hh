//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//  $Id: WinnerTakeAllRecombiner.hh 670 2014-06-06 01:24:42Z jthaler $
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

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {

//------------------------------------------------------------------------
/// \class WinnerTakeAllRecombiner
// WinnerTakeAllRecombiner defines a new recombination scheme by inheriting from JetDefinition::Recombiner.
// This scheme compares the energy of two input particles, and then combines them into a particle with
// an energy equal to the sum of the two particle energies and a direction identical to that of the harder 
// particle. This creates a jet with an axis guaranteed to align with a particle in the event.
class WinnerTakeAllRecombiner : public fastjet::JetDefinition::Recombiner {
public:
	// Constructor to choose value of alpha (defaulted to 1 for normal pT sum)
   WinnerTakeAllRecombiner(double alpha = 1.0) : _alpha(alpha) {}

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
