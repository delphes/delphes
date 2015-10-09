//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//  $Id: MeasureDefinition.cc 819 2015-06-12 21:23:24Z jthaler $
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


// #include "AxesRefiner.hh"
#include "MeasureDefinition.hh"

#include <iomanip>


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {

///////
//
// Measure Function
//
///////


//descriptions updated to include measure type
std::string DefaultMeasure::description() const {
   std::stringstream stream;
   stream << std::fixed << std::setprecision(2)
   << "Default Measure (should not be used directly)";
   return stream.str();
};

std::string NormalizedMeasure::description() const {
   std::stringstream stream;
   stream << std::fixed << std::setprecision(2)
   << "Normalized Measure (beta = " << _beta << ", R0 = " << _R0 << ")";
   return stream.str();
};

std::string UnnormalizedMeasure::description() const {
   std::stringstream stream;
   stream << std::fixed << std::setprecision(2)
   << "Unnormalized Measure (beta = " << _beta << ", in GeV)";
   return stream.str();
};


std::string NormalizedCutoffMeasure::description() const {
   std::stringstream stream;
   stream << std::fixed << std::setprecision(2)
   << "Normalized Cutoff Measure (beta = " << _beta << ", R0 = " << _R0 << ", Rcut = " << _Rcutoff << ")";
   return stream.str();
};

std::string UnnormalizedCutoffMeasure::description() const {
   std::stringstream stream;
   stream << std::fixed << std::setprecision(2)
   << "Unnormalized Cutoff Measure (beta = " << _beta << ", Rcut = " << _Rcutoff << ", in GeV)";
   return stream.str();
};

//std::string DeprecatedGeometricMeasure::description() const {
//   std::stringstream stream;
//   stream << std::fixed << std::setprecision(2)
//   << "Deprecated Geometric Measure (beta = " << _jet_beta << ", in GeV)";
//   return stream.str();
//};
   
//std::string DeprecatedGeometricCutoffMeasure::description() const {
//   std::stringstream stream;
//   stream << std::fixed << std::setprecision(2)
//   << "Deprecated Geometric Cutoff Measure (beta = " << _jet_beta << ", Rcut = " << _Rcutoff << ", in GeV)";
//   return stream.str();
//};
   
std::string ConicalMeasure::description() const {
   std::stringstream stream;
   stream << std::fixed << std::setprecision(2)
   << "Conical Measure (beta = " << _beta << ", Rcut = " << _Rcutoff << ", in GeV)";
   return stream.str();
};

std::string OriginalGeometricMeasure::description() const {
   std::stringstream stream;
   stream << std::fixed << std::setprecision(2)
   << "Original Geometric Measure (Rcut = " << _Rcutoff << ", in GeV)";
   return stream.str();
}; 

std::string ModifiedGeometricMeasure::description() const {
   std::stringstream stream;
   stream << std::fixed << std::setprecision(2)
   << "Modified Geometric Measure (Rcut = " << _Rcutoff << ", in GeV)";
   return stream.str();
}; 

std::string ConicalGeometricMeasure::description() const {
   std::stringstream stream;
   stream << std::fixed << std::setprecision(2)
   << "Conical Geometric Measure (beta = " << _jet_beta << ", gamma = " << _beam_gamma << ", Rcut = " << _Rcutoff << ", in GeV)";
   return stream.str();
}; 
   

std::string XConeMeasure::description() const {
   std::stringstream stream;
   stream << std::fixed << std::setprecision(2)
   << "XCone Measure (beta = " << _jet_beta << ", Rcut = " << _Rcutoff << ", in GeV)";
   return stream.str();
}; 
    
// Return all of the necessary TauComponents for specific input particles and axes
TauComponents MeasureDefinition::component_result(const std::vector<fastjet::PseudoJet>& particles,
                                                  const std::vector<fastjet::PseudoJet>& axes) const {
   
   // first find partition
   TauPartition partition = get_partition(particles,axes);
   
   // then return result calculated from partition
   return component_result_from_partition(partition,axes);
}

TauPartition MeasureDefinition::get_partition(const std::vector<fastjet::PseudoJet>& particles,
                                            const std::vector<fastjet::PseudoJet>& axes) const {
   
   TauPartition myPartition(axes.size());
   
   // Figures out the partiting of the input particles into the various jet pieces
   // Based on which axis the parition is closest to
   for (unsigned i = 0; i < particles.size(); i++) {
      
      // find minimum distance; start with beam (-1) for reference
      int j_min = -1;
      double minRsq;
      if (has_beam()) minRsq = beam_distance_squared(particles[i]);
      else minRsq = std::numeric_limits<double>::max(); // make it large value
      
      
      // check to see which axis the particle is closest to
      for (unsigned j = 0; j < axes.size(); j++) {
         double tempRsq = jet_distance_squared(particles[i],axes[j]); // delta R distance
         
         if (tempRsq < minRsq) {
            minRsq = tempRsq;
            j_min = j;
         }
      }
      
      if (j_min == -1) {
         assert(has_beam());  // should have beam for this to make sense.
         myPartition.push_back_beam(particles[i],i);
      } else {
         myPartition.push_back_jet(j_min,particles[i],i);
      }
   }
   
   return myPartition;
}

// Uses existing partition and calculates result
// TODO:  Can we cache this for speed up when doing area subtraction?
TauComponents MeasureDefinition::component_result_from_partition(const TauPartition& partition,
                                                     const std::vector<fastjet::PseudoJet>& axes) const {
   
   std::vector<double> jetPieces(axes.size(), 0.0);
   double beamPiece = 0.0;
   
   double tauDen = 0.0;
   if (!has_denominator()) tauDen = 1.0;  // if no denominator, then 1.0 for no normalization factor
   
   // first find jet pieces
   for (unsigned j = 0; j < axes.size(); j++) {
      std::vector<PseudoJet> thisPartition = partition.jet(j).constituents();
      for (unsigned i = 0; i < thisPartition.size(); i++) {
         jetPieces[j] += jet_numerator(thisPartition[i],axes[j]); //numerator jet piece
         if (has_denominator()) tauDen += denominator(thisPartition[i]); // denominator
      }
   }
   
   // then find beam piece
   if (has_beam()) {
      std::vector<PseudoJet> beamPartition = partition.beam().constituents();

      for (unsigned i = 0; i < beamPartition.size(); i++) {
         beamPiece += beam_numerator(beamPartition[i]); //numerator beam piece
         if (has_denominator()) tauDen += denominator(beamPartition[i]); // denominator
      }
   }
   
   // create jets for storage in TauComponents
   std::vector<PseudoJet> jets = partition.jets();
   
   return TauComponents(_tau_mode, jetPieces, beamPiece, tauDen, jets, axes);
}

// new methods added to generalize energy and angle squared for different measure types
double DefaultMeasure::energy(const PseudoJet& jet) const {
   double energy;
   switch (_measure_type) {
      case pt_R : 
      case perp_lorentz_dot : 
         energy = jet.perp();
         break;
      case E_theta : 
      case lorentz_dot : 
         energy = jet.e();
         break;
      default : {
         assert(_measure_type == pt_R || _measure_type == E_theta || _measure_type == lorentz_dot || _measure_type == perp_lorentz_dot);
         energy = std::numeric_limits<double>::quiet_NaN();
         break;
      }
   }
   return energy;
}
   
double DefaultMeasure::angleSquared(const PseudoJet& jet1, const PseudoJet& jet2) const {
   double pseudoRsquared;
   switch(_measure_type) {
      case pt_R : {
         pseudoRsquared = jet1.squared_distance(jet2);
         break; 
      }
      case E_theta : {
         // doesn't seem to be a fastjet built in for this
         double dot = jet1.px()*jet2.px() + jet1.py()*jet2.py() + jet1.pz()*jet2.pz();
         double norm1 = sqrt(jet1.px()*jet1.px() + jet1.py()*jet1.py() + jet1.pz()*jet1.pz());
         double norm2 = sqrt(jet2.px()*jet2.px() + jet2.py()*jet2.py() + jet2.pz()*jet2.pz());
        
         double costheta = dot/(norm1 * norm2);
         if (costheta > 1.0) costheta = 1.0; // Need to handle case of numerical overflow
         double theta = acos(costheta);
         pseudoRsquared = theta*theta;   
         break;
      }
      case lorentz_dot : {
         double dotproduct = dot_product(jet1,jet2);
         pseudoRsquared = 2.0 * dotproduct / (jet1.e() * jet2.e());
         break;
      }
      case perp_lorentz_dot : {
         PseudoJet lightJet = lightFrom(jet2); // assuming jet2 is the axis
         double dotproduct = dot_product(jet1,lightJet);
         pseudoRsquared = 2.0 * dotproduct / (lightJet.pt() * jet1.pt());
         break;
      }
      default : {
         assert(_measure_type == pt_R || _measure_type == E_theta || _measure_type == lorentz_dot || _measure_type == perp_lorentz_dot);
         pseudoRsquared = std::numeric_limits<double>::quiet_NaN();
         break;
      }
   }

   return pseudoRsquared;

}


///////
//
// Axes Refining
//
///////

// uses minimization of N-jettiness to continually update axes until convergence.
// The function returns the axes found at the (local) minimum
// This is the general axes refiner that can be used for a generic measure (but is
// overwritten in the case of the conical measure and the deprecated geometric measure)
std::vector<fastjet::PseudoJet> MeasureDefinition::get_one_pass_axes(int n_jets,
                                                                     const std::vector <fastjet::PseudoJet> & particles,
                                                                     const std::vector<fastjet::PseudoJet>& currentAxes,
                                                                     int nAttempts,
                                                                     double accuracy) const {

   assert(n_jets == (int)currentAxes.size());
   
   std::vector<fastjet::PseudoJet> seedAxes = currentAxes;

   std::vector<fastjet::PseudoJet> temp_axes(seedAxes.size(),fastjet::PseudoJet(0,0,0,0));
   for (unsigned int k = 0; k < seedAxes.size(); k++) {
      seedAxes[k] = lightFrom(seedAxes[k]) * seedAxes[k].E(); // making light-like, but keeping energy
   }
   
   double seedTau = result(particles, seedAxes);

   std::vector<fastjet::PseudoJet> bestAxesSoFar = seedAxes;
   double bestTauSoFar = seedTau;
   
   for (int i_att = 0; i_att < nAttempts; i_att++) {
      
      std::vector<fastjet::PseudoJet> newAxes(seedAxes.size(),fastjet::PseudoJet(0,0,0,0));
      std::vector<fastjet::PseudoJet> summed_jets(seedAxes.size(), fastjet::PseudoJet(0,0,0,0));      
      
      // find closest axis and assign to that
      for (unsigned int i = 0; i < particles.size(); i++) {
         
         // start from unclustered beam measure
         int minJ = -1;
         double minDist = beam_distance_squared(particles[i]);
         
         // which axis am I closest to?
         for (unsigned int j = 0; j < seedAxes.size(); j++) {
            double tempDist = jet_distance_squared(particles[i],seedAxes[j]);
            if (tempDist < minDist) {
               minDist = tempDist;
               minJ = j;
            }
         }
         
         // if not unclustered, then cluster
         if (minJ != -1) {
            summed_jets[minJ] += particles[i]; // keep track of energy to use later.
            if (_useAxisScaling) {
               double pseudoMomentum = dot_product(lightFrom(seedAxes[minJ]),particles[i]) + accuracy; // need small offset to avoid potential divide by zero issues
               double axis_scaling = (double)jet_numerator(particles[i], seedAxes[minJ])/pseudoMomentum;

               newAxes[minJ] += particles[i]*axis_scaling;
            }
         }
      }
      if (!_useAxisScaling) newAxes = summed_jets;

      // convert the axes to LightLike and then back to PseudoJet
      for (unsigned int k = 0; k < newAxes.size(); k++) {
         if (newAxes[k].perp() > 0) {
            newAxes[k] = lightFrom(newAxes[k]);
            newAxes[k] *= summed_jets[k].E(); // scale by energy to get sensible result
         }
      }

      // calculate tau on new axes
      double newTau = result(particles, newAxes);
      
      // find the smallest value of tau (and the corresponding axes) so far
      if (newTau < bestTauSoFar) {
         bestAxesSoFar = newAxes;
         bestTauSoFar = newTau;
      }

      if (fabs(newTau - seedTau) < accuracy) {// close enough for jazz
         seedAxes = newAxes;
         seedTau = newTau;
         break;
      }

      seedAxes = newAxes;
      seedTau = newTau;

}

   // return the axes corresponding to the smallest tau found throughout all iterations
   // this is to prevent the minimization from returning a non-minimized of tau due to potential oscillations around the minimum
   return bestAxesSoFar;

}


// One pass minimization for the DefaultMeasure
   
// Given starting axes, update to find better axes by using Kmeans clustering around the old axes
template <int N>
std::vector<LightLikeAxis> DefaultMeasure::UpdateAxesFast(const std::vector <LightLikeAxis> & old_axes,
                                                          const std::vector <fastjet::PseudoJet> & inputJets,
                                                          double accuracy
                                                          ) const {
   assert(old_axes.size() == N);
   
   // some storage, declared static to save allocation/re-allocation costs
   static LightLikeAxis new_axes[N];
   static fastjet::PseudoJet new_jets[N];
   for (int n = 0; n < N; ++n) {
      new_axes[n].reset(0.0,0.0,0.0,0.0);
      new_jets[n].reset_momentum(0.0,0.0,0.0,0.0);
   }

   double precision = accuracy;  //TODO: actually cascade this in
   
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
//  TODO:  Consider removing this in a future version
std::vector<LightLikeAxis> DefaultMeasure::UpdateAxes(const std::vector <LightLikeAxis> & old_axes,
                                                      const std::vector <fastjet::PseudoJet> & inputJets,
                                                      double accuracy) const {
   int N = old_axes.size();
   switch (N) {
      case 1: return UpdateAxesFast<1>(old_axes, inputJets, accuracy);
      case 2: return UpdateAxesFast<2>(old_axes, inputJets, accuracy);
      case 3: return UpdateAxesFast<3>(old_axes, inputJets, accuracy);
      case 4: return UpdateAxesFast<4>(old_axes, inputJets, accuracy);
      case 5: return UpdateAxesFast<5>(old_axes, inputJets, accuracy);
      case 6: return UpdateAxesFast<6>(old_axes, inputJets, accuracy);
      case 7: return UpdateAxesFast<7>(old_axes, inputJets, accuracy);
      case 8: return UpdateAxesFast<8>(old_axes, inputJets, accuracy);
      case 9: return UpdateAxesFast<9>(old_axes, inputJets, accuracy);
      case 10: return UpdateAxesFast<10>(old_axes, inputJets, accuracy);
      case 11: return UpdateAxesFast<11>(old_axes, inputJets, accuracy);
      case 12: return UpdateAxesFast<12>(old_axes, inputJets, accuracy);
      case 13: return UpdateAxesFast<13>(old_axes, inputJets, accuracy);
      case 14: return UpdateAxesFast<14>(old_axes, inputJets, accuracy);
      case 15: return UpdateAxesFast<15>(old_axes, inputJets, accuracy);
      case 16: return UpdateAxesFast<16>(old_axes, inputJets, accuracy);
      case 17: return UpdateAxesFast<17>(old_axes, inputJets, accuracy);
      case 18: return UpdateAxesFast<18>(old_axes, inputJets, accuracy);
      case 19: return UpdateAxesFast<19>(old_axes, inputJets, accuracy);
      case 20: return UpdateAxesFast<20>(old_axes, inputJets, accuracy);
      default: std::cout << "N-jettiness is hard-coded to only allow up to 20 jets!" << std::endl;
         return std::vector<LightLikeAxis>();
   }

}

// uses minimization of N-jettiness to continually update axes until convergence.
// The function returns the axes found at the (local) minimum
std::vector<fastjet::PseudoJet> DefaultMeasure::get_one_pass_axes(int n_jets,
                                                                  const std::vector <fastjet::PseudoJet> & inputJets,
                                                                  const std::vector<fastjet::PseudoJet>& seedAxes,
                                                                  int nAttempts,
                                                                  double accuracy
                                                                  ) const {

   // if the measure type doesn't use the pt_R metric, then the standard minimization scheme should be used   
   if (_measure_type != pt_R) {
      return MeasureDefinition::get_one_pass_axes(n_jets, inputJets, seedAxes, nAttempts, accuracy);
   }

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

   while (cmp > accuracy && h < nAttempts) { // Keep updating axes until near-convergence or too many update steps
      cmp = 0.0;
      h++;
      new_axes = UpdateAxes(old_axes, inputJets,accuracy); // Update axes
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
      double seed_tau = result(inputJets, seedAxes);
      double outputTau = result(inputJets, outputAxes);
      assert(outputTau <= seed_tau);
   }
   
   return outputAxes;
}
   
//// One-pass minimization for the Deprecated Geometric Measure
//// Uses minimization of the geometric distance in order to find the minimum axes.
//// It continually updates until it reaches convergence or it reaches the maximum number of attempts.
//// This is essentially the same as a stable cone finder.
//std::vector<fastjet::PseudoJet> DeprecatedGeometricCutoffMeasure::get_one_pass_axes(int n_jets,
//                                                                                    const std::vector <fastjet::PseudoJet> & particles,
//                                                                                    const std::vector<fastjet::PseudoJet>& currentAxes,
//                                                                                    int nAttempts,
//                                                                                    double accuracy) const {
//
//   assert(n_jets == (int)currentAxes.size()); //added int casting to get rid of compiler warning
//   
//   std::vector<fastjet::PseudoJet> seedAxes = currentAxes;
//   double seedTau = result(particles, seedAxes);
//   
//   for (int i = 0; i < nAttempts; i++) {
//      
//      std::vector<fastjet::PseudoJet> newAxes(seedAxes.size(),fastjet::PseudoJet(0,0,0,0));
//      
//      // find closest axis and assign to that
//      for (unsigned int i = 0; i < particles.size(); i++) {
//         
//         // start from unclustered beam measure
//         int minJ = -1;
//         double minDist = beam_distance_squared(particles[i]);
//         
//         // which axis am I closest to?
//         for (unsigned int j = 0; j < seedAxes.size(); j++) {
//            double tempDist = jet_distance_squared(particles[i],seedAxes[j]);
//            if (tempDist < minDist) {
//               minDist = tempDist;
//               minJ = j;
//            }
//         }
//         
//         // if not unclustered, then cluster
//         if (minJ != -1) newAxes[minJ] += particles[i];
//      }
//
//      // calculate tau on new axes
//      seedAxes = newAxes;
//      double tempTau = result(particles, newAxes);
//      
//      // close enough to stop?
//      if (fabs(tempTau - seedTau) < accuracy) break;
//      seedTau = tempTau;
//   }
//   
//   return seedAxes;
//}


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
