//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
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

#include "WinnerTakeAllRecombiner.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

std::string WinnerTakeAllRecombiner::description() const {
   return "Winner Take All scheme recombination";
}

// recombine pa and pb by creating pab with energy of the sum of particle energies in the direction of the harder particle

// updated recombiner to use more general form of a metric equal to E*(pT/E)^(alpha), which reduces to pT*cosh(rap)^(1-alpha)
// alpha is specified by the user. The default is alpha = 1, which is the typical behavior. alpha = 2 provides a metric which more 
// favors central jets
void WinnerTakeAllRecombiner::recombine(const fastjet::PseudoJet & pa, const fastjet::PseudoJet & pb, fastjet::PseudoJet & pab) const {
   double a_pt = pa.perp(), b_pt = pb.perp(), a_rap = pa.rap(), b_rap = pb.rap();
   
   // special case of alpha = 1, everything is just pt (made separate so that pow function isn't called)
   if (_alpha == 1.0) {
      if (a_pt >= b_pt) {
         pab.reset_PtYPhiM(a_pt + b_pt, a_rap, pa.phi());
      }
      else if (b_pt > a_pt) {
         pab.reset_PtYPhiM(a_pt + b_pt, b_rap, pb.phi());
      }
   }

   // every other case uses additional cosh(rap) term
   else {
      double a_metric = a_pt*pow(cosh(a_rap), 1.0-_alpha);
      double b_metric = b_pt*pow(cosh(b_rap), 1.0-_alpha);
      if (a_metric >= b_metric) {
   	  double new_pt = a_pt + b_pt*pow(cosh(b_rap)/cosh(a_rap), 1.0-_alpha);
   	  pab.reset_PtYPhiM(new_pt, a_rap, pa.phi());
      }
      if (b_metric > a_metric) {
   	  double new_pt = b_pt + a_pt*pow(cosh(a_rap)/cosh(b_rap), 1.0-_alpha);
   	  pab.reset_PtYPhiM(new_pt, b_rap, pb.phi());
      }
   }
}

} //namespace contrib

FASTJET_END_NAMESPACE
