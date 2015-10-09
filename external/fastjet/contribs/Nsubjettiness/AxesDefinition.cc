//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//  $Id: AxesDefinition.cc 833 2015-07-23 14:35:23Z jthaler $
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

#include "AxesDefinition.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {
  
// Repeatedly calls the one pass finder to try to find global minimum
std::vector<fastjet::PseudoJet> AxesDefinition::get_multi_pass_axes(int n_jets,
                                                                    const std::vector <fastjet::PseudoJet> & inputJets,
                                                                    const std::vector <fastjet::PseudoJet> & seedAxes,
                                                                    const MeasureDefinition* measure) const {
   
   assert(n_jets == (int)seedAxes.size()); //added int casting to get rid of compiler warning
   
   // first iteration
   std::vector<fastjet::PseudoJet> bestAxes = measure->get_one_pass_axes(n_jets, inputJets, seedAxes,_nAttempts,_accuracy);
   
   double bestTau = measure->result(inputJets,bestAxes);
   
   for (int l = 1; l < _Npass; l++) { // Do minimization procedure multiple times (l = 1 to start since first iteration is done already)
      
      // Add noise to current best axes
      std::vector< PseudoJet > noiseAxes(n_jets, PseudoJet(0,0,0,0));
      for (int k = 0; k < n_jets; k++) {
         noiseAxes[k] = jiggle(bestAxes[k]);
      }
      
      std::vector<fastjet::PseudoJet> testAxes = measure->get_one_pass_axes(n_jets, inputJets, noiseAxes,_nAttempts,_accuracy);
      double testTau = measure->result(inputJets,testAxes);
      
      if (testTau < bestTau) {
         bestTau = testTau;
         bestAxes = testAxes;
      }
   }
   
   return bestAxes;
}

   
// Used by multi-pass minimization to jiggle axes with noise range.
PseudoJet AxesDefinition::jiggle(const PseudoJet& axis) const {
   double phi_noise = ((double)rand()/(double)RAND_MAX) * _noise_range * 2.0 - _noise_range;
   double rap_noise = ((double)rand()/(double)RAND_MAX) * _noise_range * 2.0 - _noise_range;
   
   double new_phi = axis.phi() + phi_noise;
   if (new_phi >= 2.0*M_PI) new_phi -= 2.0*M_PI;
   if (new_phi <= -2.0*M_PI) new_phi += 2.0*M_PI;
   
   PseudoJet newAxis(0,0,0,0);
   newAxis.reset_PtYPhiM(axis.perp(),axis.rap() + rap_noise,new_phi);
   return newAxis;
}

   
LimitedWarning HardestJetAxes::_too_few_axes_warning;
LimitedWarning ExclusiveJetAxes::_too_few_axes_warning;
LimitedWarning ExclusiveCombinatorialJetAxes::_too_few_axes_warning;

std::vector<fastjet::PseudoJet> Manual_Axes::get_starting_axes(int,
                                                               const std::vector<fastjet::PseudoJet>&,
                                                               const MeasureDefinition *) const {
   // This is a function dummy and should never be called
   assert(false);
   std::vector<fastjet::PseudoJet> dummy;
   return dummy;
}

} // namespace contrib

FASTJET_END_NAMESPACE
