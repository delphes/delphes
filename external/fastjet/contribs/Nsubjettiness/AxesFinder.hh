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


#ifndef __FASTJET_CONTRIB_AXESFINDER_HH__
#define __FASTJET_CONTRIB_AXESFINDER_HH__

#include "WinnerTakeAllRecombiner.hh"
#include "MeasureFunction.hh"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"

#include <cmath>
#include <vector>
#include <list>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

///////
//
// Axes Finder Options
//
///////

//------------------------------------------------------------------------
/// \class AxesFinder
// This is the base class for all axes finders. These axes are used along with the MeasureFunctions to calculate 
// tau_N. There are different implementations of axes finding that are defined in derived classes below.
class AxesFinder {
   
protected:
   AxesFinder* _startingFinder; // storing a possible starting finder if needed
   std::vector<fastjet::PseudoJet> _seedAxes;
   
   AxesFinder(AxesFinder* startingFinder = NULL) : _startingFinder(startingFinder) {}
   
public:
   virtual ~AxesFinder(){
      if (_startingFinder) delete _startingFinder;  //TODO: Convert to smart pointers to avoid this.
   }
   
   // Allow setting of seedAxes from a starting finder
   std::vector<fastjet::PseudoJet> getAxes(int n_jets, const std::vector<fastjet::PseudoJet> & inputs, const std::vector<fastjet::PseudoJet>& currentAxes) {
      if (_startingFinder) {
         _seedAxes = _startingFinder->getAxes(n_jets,inputs,currentAxes);
         return getBetterAxes(n_jets,inputs,_seedAxes);
      } else {
         _seedAxes = getBetterAxes(n_jets,inputs,currentAxes);
         return _seedAxes;
      }
   }
   
   // say what the current seed axes are
   std::vector<fastjet::PseudoJet> seedAxes() const {
      return _seedAxes;
   }
   
   // This function should be overloaded, and updates the seedAxes
   virtual std::vector<fastjet::PseudoJet> getBetterAxes(int n_jets, const std::vector<fastjet::PseudoJet> & inputs, const std::vector<fastjet::PseudoJet>& seedAxes) = 0;
   
};

//------------------------------------------------------------------------
/// \class AxesFinderFromExclusiveJetDefinition
// This class finds axes by clustering the particles and then finding the exclusive jets. This can be implemented
// with different jet algorithms.
class AxesFinderFromExclusiveJetDefinition : public AxesFinder {

   private:
      fastjet::JetDefinition _def;
   
   public:
      AxesFinderFromExclusiveJetDefinition(fastjet::JetDefinition def) : _def(def) {}
      
      virtual std::vector<fastjet::PseudoJet> getBetterAxes(int n_jets, const std::vector <fastjet::PseudoJet> & inputs, const std::vector<fastjet::PseudoJet>& currentAxes) {
         fastjet::ClusterSequence jet_clust_seq(inputs, _def);
         return jet_clust_seq.exclusive_jets(n_jets);
      }
};

//------------------------------------------------------------------------
/// \class AxesFinderFromWTA_KT
// This class finds axes by finding the exlusive jets after clustering according to a kT algorithm and a 
// winner take all recombination scheme.
class AxesFinderFromWTA_KT : public AxesFinderFromExclusiveJetDefinition { 
   private: 
      const WinnerTakeAllRecombiner *recomb;
   public:
      AxesFinderFromWTA_KT() : AxesFinderFromExclusiveJetDefinition(
         fastjet::JetDefinition(fastjet::kt_algorithm, 
         fastjet::JetDefinition::max_allowable_R, //maximum jet radius constant
         recomb = new WinnerTakeAllRecombiner(), 
         fastjet::Best)) {}
      ~AxesFinderFromWTA_KT() {delete recomb;}
   };
   
//------------------------------------------------------------------------
/// \class AxesFinderFromWTA_CA
// This class finds axes by finding the exlusive jets after clustering according to a CA algorithm and a 
// winner take all recombination scheme.
class AxesFinderFromWTA_CA : public AxesFinderFromExclusiveJetDefinition {
   private: 
      const WinnerTakeAllRecombiner *recomb;
   public:
      AxesFinderFromWTA_CA() : AxesFinderFromExclusiveJetDefinition(
         fastjet::JetDefinition(fastjet::cambridge_algorithm, 
         fastjet::JetDefinition::max_allowable_R, //maximum jet radius constant
         recomb = new WinnerTakeAllRecombiner(), 
         fastjet::Best)) {}
      ~AxesFinderFromWTA_CA() {delete recomb;}
};

//  The following classes are for testing, and are commented out for initial release
//
////------------------------------------------------------------------------
///// \class AxesFinderFromWTA2_KT
//// This class finds axes by finding the exlusive jets after clustering according to a kT algorithm and a 
//// winner take all recombination scheme with alpha = 2.
//class AxesFinderFromWTA2_KT : public AxesFinderFromExclusiveJetDefinition { 
//   private: 
//      const WinnerTakeAllRecombiner *recomb;
//   public:
//      AxesFinderFromWTA2_KT() : AxesFinderFromExclusiveJetDefinition(
//         fastjet::JetDefinition(fastjet::kt_algorithm, 
//         fastjet::JetDefinition::max_allowable_R, //maximum jet radius constant
//         recomb = new WinnerTakeAllRecombiner(2), // uses alpha = 2 here
//         fastjet::Best)) {}
//      ~AxesFinderFromWTA2_KT() {delete recomb;}
//   };
//   
////------------------------------------------------------------------------
///// \class AxesFinderFromWTA2_CA
//// This class finds axes by finding the exlusive jets after clustering according to a CA algorithm and a 
//// winner take all recombination scheme with alpha = 2.
//class AxesFinderFromWTA2_CA : public AxesFinderFromExclusiveJetDefinition {
//   private: 
//      const WinnerTakeAllRecombiner *recomb;
//   public:
//      AxesFinderFromWTA2_CA() : AxesFinderFromExclusiveJetDefinition(
//         fastjet::JetDefinition(fastjet::cambridge_algorithm, 
//         fastjet::JetDefinition::max_allowable_R, //maximum jet radius constant
//         recomb = new WinnerTakeAllRecombiner(2), //uses alpha = 2 here
//         fastjet::Best)) {}
//      ~AxesFinderFromWTA2_CA() {delete recomb;}
//};

//------------------------------------------------------------------------
/// \class AxesFinderFromKT
// This class finds axes by finding the exlusive jets after clustering according to a kT algorithm and a 
// E_scheme recombination.
class AxesFinderFromKT : public AxesFinderFromExclusiveJetDefinition {
   public:
      AxesFinderFromKT() : AxesFinderFromExclusiveJetDefinition(
         fastjet::JetDefinition(fastjet::kt_algorithm,
         fastjet::JetDefinition::max_allowable_R, //maximum jet radius constant
         fastjet::E_scheme,
         fastjet::Best)) {}
};

//------------------------------------------------------------------------
/// \class AxesFinderFromCA
// This class finds axes by finding the exlusive jets after clustering according to a CA algorithm and a 
// E_scheme recombination.
class AxesFinderFromCA : public AxesFinderFromExclusiveJetDefinition { 
   public:
      AxesFinderFromCA() : AxesFinderFromExclusiveJetDefinition(
         fastjet::JetDefinition(fastjet::cambridge_algorithm,
                                fastjet::JetDefinition::max_allowable_R,  //maximum jet radius constant
         fastjet::E_scheme,
         fastjet::Best)) {}
};


//------------------------------------------------------------------------
/// \class AxesFinderFromHardestJetDefinition
// This class finds axes by clustering the particles and then finding the n hardest inclusive jets. 
// This can be implemented with different jet algorithms.
class AxesFinderFromHardestJetDefinition : public AxesFinder {

   private:
      fastjet::JetDefinition _def;
   
   public:
      AxesFinderFromHardestJetDefinition(fastjet::JetDefinition def) : _def(def) {}
      
      virtual std::vector<fastjet::PseudoJet> getBetterAxes(int n_jets, const std::vector <fastjet::PseudoJet> & inputs, const std::vector<fastjet::PseudoJet>& currentAxes) {
         fastjet::ClusterSequence jet_clust_seq(inputs, _def);
         std::vector<fastjet::PseudoJet> myJets = sorted_by_pt(jet_clust_seq.inclusive_jets());
         myJets.resize(n_jets);  // only keep n hardest
         return myJets;
      }      
};

//------------------------------------------------------------------------
/// \class AxesFinderFromAntiKT
// This class finds axes by finding the n hardest jets after clustering the particles according 
// to an anti kT algorithm and E_scheme.
class AxesFinderFromAntiKT : public AxesFinderFromHardestJetDefinition {
   public:
      AxesFinderFromAntiKT(double R0) : AxesFinderFromHardestJetDefinition(fastjet::JetDefinition(fastjet::antikt_algorithm,R0,fastjet::E_scheme,fastjet::Best)) {}
};


//------------------------------------------------------------------------
/// \class AxesFinderFromUserInput
// This class allows the user to manually define the axes.
class AxesFinderFromUserInput : public AxesFinder {

   public:
      AxesFinderFromUserInput() {}
      
      virtual std::vector<fastjet::PseudoJet> getBetterAxes(int n_jets, const std::vector <fastjet::PseudoJet> & inputs, const std::vector<fastjet::PseudoJet>& currentAxes) {
         assert(currentAxes.size() == (unsigned int) n_jets);
         return currentAxes;
      }
};

//This is a helper class for the Minimum Axes Finders. It is defined later.
class LightLikeAxis;                                          


//------------------------------------------------------------------------
/// \class AxesFinderFromOnePassMinimization
// This class defines an AxesFinder that uses Kmeans minimization, but only on a single pass.
class AxesFinderFromOnePassMinimization : public AxesFinder {

   private:
      double _precision;  // Desired precision in axes alignment
      int _halt;  // maximum number of steps per iteration
      
      double _beta;
      double _Rcutoff;
      
      DefaultUnnormalizedMeasure _measureFunction;
   
   public:

      // From a startingFinder, try to minimize the unnormalized_measure
      AxesFinderFromOnePassMinimization(AxesFinder* startingFinder, double beta, double Rcutoff)
         : AxesFinder(startingFinder),
           _precision(0.0001), //hard coded for now
           _halt(1000), //hard coded for now
           _beta(beta),
           _Rcutoff(Rcutoff),
           _measureFunction(beta, Rcutoff)
           {}
   
      virtual std::vector<fastjet::PseudoJet> getBetterAxes(int n_jets, const std::vector <fastjet::PseudoJet> & inputJets, const std::vector<fastjet::PseudoJet>& currentAxes);

      template <int N> std::vector<LightLikeAxis> UpdateAxesFast(const std::vector <LightLikeAxis> & old_axes,
                                  const std::vector <fastjet::PseudoJet> & inputJets);
   
      std::vector<LightLikeAxis> UpdateAxes(const std::vector <LightLikeAxis> & old_axes, 
                                      const std::vector <fastjet::PseudoJet> & inputJets);

};


//------------------------------------------------------------------------
/// \class AxesFinderFromKmeansMinimization
// This class finds finds axes by using Kmeans clustering to minimizaiton N-jettiness. Given a first set of 
// starting axes, it updates n times to get as close to the global minimum as possible. This class calls OnePass many times,
// added noise to the axes.
class AxesFinderFromKmeansMinimization : public AxesFinder{

   private:
      int _n_iterations;   // Number of iterations to run  (0 for no minimization, 1 for one-pass, >>1 for global minimum)
      double _noise_range; // noise range for random initialization
   
      DefaultUnnormalizedMeasure _measureFunction; //function to test whether minimum is reached
   
      AxesFinderFromOnePassMinimization _onePassFinder;  //one pass finder for minimization

      PseudoJet jiggle(const PseudoJet& axis);
   
   public:
      AxesFinderFromKmeansMinimization(AxesFinder *startingFinder, double beta, double Rcutoff, int n_iterations) :
         AxesFinder(startingFinder),
         _n_iterations(n_iterations),
         _noise_range(1.0), // hard coded for the time being
         _measureFunction(beta, Rcutoff),
         _onePassFinder(NULL, beta, Rcutoff)
         {}

      virtual std::vector<fastjet::PseudoJet> getBetterAxes(int n_jets, const std::vector <fastjet::PseudoJet> & inputJets, const std::vector<fastjet::PseudoJet>& currentAxes);

};

//------------------------------------------------------------------------
/// \class AxesFinderFromGeometricMinimization
// This class finds axes by minimizing the Lorentz dot product distance between axes and particles. Given a first set of starting axes,
// it essentially does stable cone finxing.
class AxesFinderFromGeometricMinimization : public AxesFinder {

   private:
      MeasureFunction* _function;
      double _Rcutoff;
      double _nAttempts;
      double _accuracy;

   
   public:
      AxesFinderFromGeometricMinimization(AxesFinder* startingFinder, double beta, double Rcutoff) : AxesFinder(startingFinder), _Rcutoff(Rcutoff) {
         if (beta != 2.0) {
            std::cerr << "Geometric minimization is currently only defined for beta = 2.0." << std::endl;
            exit(1);
         }
         
         _nAttempts = 100;
         _accuracy = 0.000000001;
         _function = new GeometricMeasure(beta,_Rcutoff);
      }

      ~AxesFinderFromGeometricMinimization() {
         delete _function;
      }
   
      virtual std::vector<fastjet::PseudoJet> getBetterAxes(int n_jets, const std::vector <fastjet::PseudoJet> & particles, const std::vector<fastjet::PseudoJet>& currentAxes);
};

//------------------------------------------------------------------------
/// \class LightLikeAxis
// This is a helper class for the minimum Axes Finders classes above. It creates a convenient way of defining axes
// in order to better facilitate calculations.
class LightLikeAxis {
private:
   double _rap, _phi, _weight, _mom;
   
   double DistanceSq(double rap2, double phi2) const {
      double rap1 = _rap;
      double phi1 = _phi;
      
      double distRap = rap1-rap2;
      double distPhi = std::fabs(phi1-phi2);
      if (distPhi > M_PI) {distPhi = 2.0*M_PI - distPhi;}
      return sq(distRap) + sq(distPhi);
   }
   
   double Distance(double rap2, double phi2) const {
      return std::sqrt(DistanceSq(rap2,phi2));
   }
   
   
public:
   LightLikeAxis() : _rap(0.0), _phi(0.0), _weight(0.0), _mom(0.0) {}
   LightLikeAxis(double my_rap, double my_phi, double my_weight, double my_mom) :
   _rap(my_rap), _phi(my_phi), _weight(my_weight), _mom(my_mom) {}
   double rap() const {return _rap;}
   double phi() const {return _phi;}
   double weight() const {return _weight;}
   double mom() const {return _mom;}
   void set_rap(double my_set_rap) {_rap = my_set_rap;}
   void set_phi(double my_set_phi) {_phi = my_set_phi;}
   void set_weight(double my_set_weight) {_weight = my_set_weight;}
   void set_mom(double my_set_mom) {_mom = my_set_mom;}
   void reset(double my_rap, double my_phi, double my_weight, double my_mom) {_rap=my_rap; _phi=my_phi; _weight=my_weight; _mom=my_mom;}

   // return PseudoJet with information
   fastjet::PseudoJet ConvertToPseudoJet();
   
   double DistanceSq(const fastjet::PseudoJet& input) const {
      return DistanceSq(input.rap(),input.phi());
   }

   double Distance(const fastjet::PseudoJet& input) const {
      return std::sqrt(DistanceSq(input));
   }

   double DistanceSq(const LightLikeAxis& input) const {
      return DistanceSq(input.rap(),input.phi());
   }

   double Distance(const LightLikeAxis& input) const {
      return std::sqrt(DistanceSq(input));
   }

};

} //namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_AXESFINDER_HH__
