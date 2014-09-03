//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//  $Id: Nsubjettiness.cc 597 2014-04-16 23:07:55Z jthaler $
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

#include "Nsubjettiness.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {

//result returns tau_N with normalization dependent on what is specified in constructor
double Nsubjettiness::result(const PseudoJet& jet) const {
   std::vector<fastjet::PseudoJet> particles = jet.constituents();
   return _njettinessFinder.getTau(_N, particles);
}

TauComponents Nsubjettiness::component_result(const PseudoJet& jet) const {
   std::vector<fastjet::PseudoJet> particles = jet.constituents();
   return _njettinessFinder.getTauComponents(_N, particles);
}

//ratio result uses Nsubjettiness result to find the ratio tau_N/tau_M, where N and M are specified by user
double NsubjettinessRatio::result(const PseudoJet& jet) const {
   double numerator = _nsub_numerator.result(jet);
   double denominator = _nsub_denominator.result(jet);
   return numerator/denominator;
}

} // namespace contrib

FASTJET_END_NAMESPACE
