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

#include "AxesFinder.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

///////
//
// Functions for minimization.
//
///////

// Given starting axes, update to find better axes by using Kmeans clustering around the old axes
template <int N>
std::vector<LightLikeAxis> AxesFinderFromOnePassMinimization::UpdateAxesFast(const std::vector <LightLikeAxis> & old_axes, 
                                  const std::vector <fastjet::PseudoJet> & inputJets) {
   assert(old_axes.size() == N);
   
   // some storage, declared static to save allocation/re-allocation costs
   static LightLikeAxis new_axes[N];
   static fastjet::PseudoJet new_jets[N];
   for (int n = 0; n < N; ++n) {
      new_axes[n].reset(0.0,0.0,0.0,0.0);
#ifdef FASTJET2
      new_jets[n].reset(0.0,0.0,0.0,0.0);
#else
      // use cheaper reset if available
      new_jets[n].reset_momentum(0.0,0.0,0.0,0.0);
#endif
   }

   double precision = _precision;
   
   /////////////// Assignment Step //////////////////////////////////////////////////////////
   std::vector<int> assignment_index(inputJets.size()); 
   int k_assign = -1;
   
   for (unsigned i = 0; i < inputJets.size(); i++){
      double smallestDist = std::numeric_limits<double>::max();  //large number
      for (int k = 0; k < N; k++) {
         double thisDist = old_axes[k].DistanceSq(inputJets[i]);
         if (thisDist < smallestDist) {
            smallestDist = thisDist;
            k_assign = k;
         }
      }
      if (smallestDist > sq(_Rcutoff)) {k_assign = -1;}
      assignment_index[i] = k_assign;
   }
   
   //////////////// Update Step /////////////////////////////////////////////////////////////
   double distPhi, old_dist;
   for (unsigned i = 0; i < inputJets.size(); i++) {
      int old_jet_i = assignment_index[i];
      if (old_jet_i == -1) {continue;}

      const fastjet::PseudoJet& inputJet_i = inputJets[i];
      LightLikeAxis& new_axis_i = new_axes[old_jet_i];
      double inputPhi_i = inputJet_i.phi();
      double inputRap_i = inputJet_i.rap();
            
      // optimize pow() call
      // add noise (the precision term) to make sure we don't divide by zero
      if (_beta == 1.0) {
         double DR = std::sqrt(sq(precision) + old_axes[old_jet_i].DistanceSq(inputJet_i));
         old_dist = 1.0/DR;
      } else if (_beta == 2.0) {
         old_dist = 1.0;
      } else if (_beta == 0.0) {
         double DRSq = sq(precision) + old_axes[old_jet_i].DistanceSq(inputJet_i);
         old_dist = 1.0/DRSq;
      } else {
         old_dist = sq(precision) + old_axes[old_jet_i].DistanceSq(inputJet_i);
         old_dist = std::pow(old_dist, (0.5*_beta-1.0));
      }
      
      // TODO:  Put some of these addition functions into light-like axes
      // rapidity sum
      new_axis_i.set_rap(new_axis_i.rap() + inputJet_i.perp() * inputRap_i * old_dist);
      // phi sum
      distPhi = inputPhi_i - old_axes[old_jet_i].phi();
      if (fabs(distPhi) <= M_PI){
         new_axis_i.set_phi( new_axis_i.phi() + inputJet_i.perp() * inputPhi_i * old_dist );
      } else if (distPhi > M_PI) {
         new_axis_i.set_phi( new_axis_i.phi() + inputJet_i.perp() * (-2*M_PI + inputPhi_i) * old_dist );
      } else if (distPhi < -M_PI) {
         new_axis_i.set_phi( new_axis_i.phi() + inputJet_i.perp() * (+2*M_PI + inputPhi_i) * old_dist );
      }
      // weights sum
      new_axis_i.set_weight( new_axis_i.weight() + inputJet_i.perp() * old_dist );
      // momentum magnitude sum
      new_jets[old_jet_i] += inputJet_i;
   }
   // normalize sums
   for (int k = 0; k < N; k++) {
      if (new_axes[k].weight() == 0) {
         // no particles were closest to this axis!  Return to old axis instead of (0,0,0,0)
         new_axes[k] = old_axes[k];
      } else {
         new_axes[k].set_rap( new_axes[k].rap() / new_axes[k].weight() );
         new_axes[k].set_phi( new_axes[k].phi() / new_axes[k].weight() );
         new_axes[k].set_phi( std::fmod(new_axes[k].phi() + 2*M_PI, 2*M_PI) );
         new_axes[k].set_mom( std::sqrt(new_jets[k].modp2()) );
      }
   }
   std::vector<LightLikeAxis> new_axes_vec(N);
   for (unsigned k = 0; k < N; ++k) new_axes_vec[k] = new_axes[k];
   return new_axes_vec;
}

// Given N starting axes, this function updates all axes to find N better axes. 
// (This is just a wrapper for the templated version above.)
std::vector<LightLikeAxis> AxesFinderFromOnePassMinimization::UpdateAxes(const std::vector <LightLikeAxis> & old_axes,
                                      const std::vector <fastjet::PseudoJet> & inputJets) {
   int N = old_axes.size();
   switch (N) {
      case 1: return UpdateAxesFast<1>(old_axes, inputJets);
      case 2: return UpdateAxesFast<2>(old_axes, inputJets);
      case 3: return UpdateAxesFast<3>(old_axes, inputJets);
      case 4: return UpdateAxesFast<4>(old_axes, inputJets);
      case 5: return UpdateAxesFast<5>(old_axes, inputJets);
      case 6: return UpdateAxesFast<6>(old_axes, inputJets);
      case 7: return UpdateAxesFast<7>(old_axes, inputJets);
      case 8: return UpdateAxesFast<8>(old_axes, inputJets);
      case 9: return UpdateAxesFast<9>(old_axes, inputJets);
      case 10: return UpdateAxesFast<10>(old_axes, inputJets);
      case 11: return UpdateAxesFast<11>(old_axes, inputJets);
      case 12: return UpdateAxesFast<12>(old_axes, inputJets);
      case 13: return UpdateAxesFast<13>(old_axes, inputJets);
      case 14: return UpdateAxesFast<14>(old_axes, inputJets);
      case 15: return UpdateAxesFast<15>(old_axes, inputJets);
      case 16: return UpdateAxesFast<16>(old_axes, inputJets);
      case 17: return UpdateAxesFast<17>(old_axes, inputJets);
      case 18: return UpdateAxesFast<18>(old_axes, inputJets);
      case 19: return UpdateAxesFast<19>(old_axes, inputJets);
      case 20: return UpdateAxesFast<20>(old_axes, inputJets);
      default: std::cout << "N-jettiness is hard-coded to only allow up to 20 jets!" << std::endl;
         return std::vector<LightLikeAxis>();
   }

}

// uses minimization of N-jettiness to continually update axes until convergence.
// The function returns the axes found at the (local) minimum
std::vector<fastjet::PseudoJet> AxesFinderFromOnePassMinimization::getBetterAxes(int n_jets, const std::vector <fastjet::PseudoJet> & inputJets, const std::vector<fastjet::PseudoJet>& seedAxes) {
	  
   // convert from PseudoJets to LightLikeAxes
   std::vector< LightLikeAxis > old_axes(n_jets, LightLikeAxis(0,0,0,0));
   for (int k = 0; k < n_jets; k++) {
      old_axes[k].set_rap( seedAxes[k].rap() );
      old_axes[k].set_phi( seedAxes[k].phi() );
   }
   
   // Find new axes by iterating (only one pass here)
   std::vector< LightLikeAxis > new_axes(n_jets, LightLikeAxis(0,0,0,0));
   double cmp = std::numeric_limits<double>::max();  //large number
   int h = 0;
   while (cmp > _precision && h < _halt) { // Keep updating axes until near-convergence or too many update steps
      cmp = 0.0;
      h++;
      new_axes = UpdateAxes(old_axes, inputJets); // Update axes
      for (int k = 0; k < n_jets; k++) {
         cmp += old_axes[k].Distance(new_axes[k]);
      }
      cmp = cmp / ((double) n_jets);
      old_axes = new_axes;
   }
      
   // Convert from internal LightLikeAxes to PseudoJet
   std::vector<fastjet::PseudoJet> outputAxes;
   for (int k = 0; k < n_jets; k++) {
      fastjet::PseudoJet temp = old_axes[k].ConvertToPseudoJet();
      outputAxes.push_back(temp);
   }
   
   // this is used to debug the minimization routine to make sure that it works.
   bool do_debug = false;
   if (do_debug) {
      // get this information to make sure that minimization is working properly
      TauComponents seed_tau_components = _measureFunction.result(inputJets, seedAxes);
      double seed_tau = seed_tau_components.tau();
      TauComponents tau_components = _measureFunction.result(inputJets, outputAxes);
      double outputTau = tau_components.tau();
      assert(outputTau <= seed_tau);
   }
   
   return outputAxes;
}

PseudoJet AxesFinderFromKmeansMinimization::jiggle(const PseudoJet& axis) {
   double phi_noise = ((double)rand()/(double)RAND_MAX) * _noise_range * 2.0 - _noise_range;
   double rap_noise = ((double)rand()/(double)RAND_MAX) * _noise_range * 2.0 - _noise_range;
   
   double new_phi = axis.phi() + phi_noise;
   if (new_phi >= 2.0*M_PI) new_phi -= 2.0*M_PI;
   if (new_phi <= -2.0*M_PI) new_phi += 2.0*M_PI;

   PseudoJet newAxis(0,0,0,0);
   newAxis.reset_PtYPhiM(axis.perp(),axis.rap() + rap_noise,new_phi);
   return newAxis;
}
   
   
// Repeatedly calls the one pass finder to try to find global minimum
std::vector<fastjet::PseudoJet> AxesFinderFromKmeansMinimization::getBetterAxes(int n_jets, const std::vector <fastjet::PseudoJet> & inputJets, const std::vector<fastjet::PseudoJet>& seedAxes) {
   
   // first iteration
	std::vector<fastjet::PseudoJet> bestAxes = _onePassFinder.getAxes(n_jets, inputJets, seedAxes);
   
   double bestTau = (_measureFunction.result(inputJets,bestAxes)).tau();
   
   for (int l = 1; l < _n_iterations; l++) { // Do minimization procedure multiple times (l = 1 to start since first iteration is done already)
   
      // Add noise to current best axes
      std::vector< PseudoJet > noiseAxes(n_jets, PseudoJet(0,0,0,0));
      for (int k = 0; k < n_jets; k++) {
         noiseAxes[k] = jiggle(bestAxes[k]);
      }

      std::vector<fastjet::PseudoJet> testAxes = _onePassFinder.getAxes(n_jets, inputJets, noiseAxes);
      double testTau = (_measureFunction.result(inputJets,testAxes)).tau();
      
      if (testTau < bestTau) {
         bestTau = testTau;
         bestAxes = testAxes;
      }
   }
   
   return bestAxes;
}

// Uses minimization of the geometric distance in order to find the minimum axes.
// It continually updates until it reaches convergence or it reaches the maximum number of attempts.
// This is essentially the same as a stable cone finder.
std::vector<fastjet::PseudoJet> AxesFinderFromGeometricMinimization::getBetterAxes(int n_jets, const std::vector <fastjet::PseudoJet> & particles, const std::vector<fastjet::PseudoJet>& currentAxes) {

   std::vector<fastjet::PseudoJet> seedAxes = currentAxes;
   double seedTau = _function->tau(particles, seedAxes);
   
   for (int i = 0; i < _nAttempts; i++) {
      
      std::vector<fastjet::PseudoJet> newAxes(seedAxes.size(),fastjet::PseudoJet(0,0,0,0));
      
      // find closest axis and assign to that
      for (unsigned int i = 0; i < particles.size(); i++) {
         
         // start from unclustered beam measure
         int minJ = -1;
         double minDist = _function->beam_distance_squared(particles[i]);
         
         // which axis am I closest to?
         for (unsigned int j = 0; j < seedAxes.size(); j++) {
            double tempDist = _function->jet_distance_squared(particles[i],seedAxes[j]);
            if (tempDist < minDist) {
               minDist = tempDist;
               minJ = j;
            }
         }
         
         // if not unclustered, then cluster
         if (minJ != -1) newAxes[minJ] += particles[i];
      }
      
      // calculate tau on new axes
      seedAxes = newAxes;
      double tempTau = _function->tau(particles, newAxes);
      
      // close enough to stop?
      if (fabs(tempTau - seedTau) < _accuracy) break;
      seedTau = tempTau;
   }
   
   return seedAxes;
}

// Go from internal LightLikeAxis to PseudoJet
fastjet::PseudoJet LightLikeAxis::ConvertToPseudoJet() {
    double px, py, pz, E;
    E = _mom;
    pz = (std::exp(2.0*_rap) - 1.0) / (std::exp(2.0*_rap) + 1.0) * E;
    px = std::cos(_phi) * std::sqrt( std::pow(E,2) - std::pow(pz,2) );
    py = std::sin(_phi) * std::sqrt( std::pow(E,2) - std::pow(pz,2) );
    return fastjet::PseudoJet(px,py,pz,E);
}

} //namespace contrib

FASTJET_END_NAMESPACE
