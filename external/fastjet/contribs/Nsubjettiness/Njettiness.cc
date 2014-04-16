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

#include "Njettiness.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {


///////
//
// Main Njettiness Class
//
///////

// Helper function to correlate one pass minimization with appropriate measure
void Njettiness::setOnePassAxesFinder(MeasureMode measure_mode, AxesFinder* startingFinder, double beta, double Rcutoff) {
   if (measure_mode == normalized_measure || measure_mode == unnormalized_measure || measure_mode == normalized_cutoff_measure || measure_mode == unnormalized_cutoff_measure) {
      _axesFinder = new AxesFinderFromOnePassMinimization(startingFinder, beta, Rcutoff);
   }
   else if (measure_mode == geometric_measure || measure_mode == geometric_cutoff_measure) {
      _axesFinder = new AxesFinderFromGeometricMinimization(startingFinder, beta, Rcutoff);
   }
   else {
      std::cerr << "Minimization only set up for normalized_measure, unnormalized_measure, normalized_cutoff_measure, unnormalized_cutoff_measure, geometric_measure, geometric_cutoff_measure" << std::endl;
      exit(1); }
}

// Parsing needed for constructor to set AxesFinder and MeasureFunction
// All of the parameter handling is here, and checking that number of parameters is correct.
void Njettiness::setMeasureFunctionandAxesFinder(AxesMode axes_mode, MeasureMode measure_mode, double para1, double para2, double para3, double para4) {

   // definition of maximum Rcutoff for non-cutoff measures, changed later by other measures
   double Rcutoff = std::numeric_limits<double>::max();  //large number
   // Most (but all measures have some kind of beta value)
   double beta = NAN;
   // The normalized measures have an R0 value.
   double R0 = NAN;

   // Find the MeasureFunction and set the parameters.
   switch (measure_mode) {
      case normalized_measure:
         beta = para1;
         R0 = para2;
         if(correctParameterCount(2, para1, para2, para3, para4)) 
            _measureFunction = new DefaultNormalizedMeasure(beta, R0, Rcutoff); //normalized_measure requires 2 parameters, beta and R0
         else { 
            std::cerr << "normalized_measure needs 2 parameters (beta and R0)" << std::endl;
            exit(1); }
         break;
      case unnormalized_measure:
         beta = para1;
         if(correctParameterCount(1, para1, para2, para3, para4)) 
            _measureFunction = new DefaultUnnormalizedMeasure(beta, Rcutoff); //unnormalized_measure requires 1 parameter, beta
         else {
            std::cerr << "unnormalized_measure needs 1 parameter (beta)" << std::endl;
            exit(1); }
         break;
      case geometric_measure:
         beta = para1;
         if(correctParameterCount(1, para1, para2, para3, para4))
            _measureFunction = new GeometricMeasure(beta,Rcutoff); //geometric_measure requires 1 parameter, beta
         else {
            std::cerr << "geometric_measure needs 1 parameter (beta)" << std::endl;
            exit(1); }
         break;
      case normalized_cutoff_measure:
         beta = para1;
         R0 = para2;
         Rcutoff = para3; //Rcutoff parameter is 3rd parameter in normalized_cutoff_measure
         if(correctParameterCount(3, para1, para2, para3, para4))
            _measureFunction = new DefaultNormalizedMeasure(beta, R0, Rcutoff); //normalized_cutoff_measure requires 3 parameters, beta, R0, and Rcutoff
         else { 
            std::cerr << "normalized_cutoff_measure has 3 parameters (beta, R0, Rcutoff)" << std::endl;
            exit(1); }
         break;
      case unnormalized_cutoff_measure:
         beta = para1;
         Rcutoff = para2; //Rcutoff parameter is 2nd parameter in normalized_cutoff_measure
         if (correctParameterCount(2, para1, para2, para3, para4))
            _measureFunction = new DefaultUnnormalizedMeasure(beta, Rcutoff); //unnormalized_cutoff_measure requires 2 parameters, beta and Rcutoff
         else {
            std::cerr << "unnormalized_cutoff_measure has 2 parameters (beta, Rcutoff)" << std::endl;
            exit(1); }
         break;
      case geometric_cutoff_measure:
         beta = para1;
         Rcutoff = para2; //Rcutoff parameter is 2nd parameter in geometric_cutoff_measure
         if(correctParameterCount(2, para1, para2, para3, para4))
            _measureFunction = new GeometricMeasure(beta,Rcutoff); //geometric_cutoff_measure requires 2 parameters, beta and Rcutoff
         else {
            std::cerr << "geometric_cutoff_measure has 2 parameters (beta,Rcutoff)" << std::endl;
            exit(1); }
         break;
      default:
         assert(false);
         break;
   }   

   // Choose which AxesFinder from user input.
   // Uses setOnePassAxesFinder helpful function to use beta and Rcutoff values about (if needed)
   switch (axes_mode) {
      case wta_kt_axes:
         _axesFinder = new AxesFinderFromWTA_KT(); 
         break;
      case wta_ca_axes:
         _axesFinder = new AxesFinderFromWTA_CA(); 
         break;
      case kt_axes:
         _axesFinder = new AxesFinderFromKT();
         break;
      case ca_axes:
         _axesFinder = new AxesFinderFromCA();
         break;
      case antikt_0p2_axes:
         _axesFinder = new AxesFinderFromAntiKT(0.2);     
         break;
      case onepass_wta_kt_axes:
         setOnePassAxesFinder(measure_mode, new AxesFinderFromWTA_KT(), beta, Rcutoff);
         break;
      case onepass_wta_ca_axes:
         setOnePassAxesFinder(measure_mode, new AxesFinderFromWTA_CA(), beta, Rcutoff);
         break;
      case onepass_kt_axes:
         setOnePassAxesFinder(measure_mode, new AxesFinderFromKT(), beta, Rcutoff);
         break;
      case onepass_ca_axes:
         setOnePassAxesFinder(measure_mode, new AxesFinderFromCA(), beta, Rcutoff);
         break;
      case onepass_antikt_0p2_axes:
         setOnePassAxesFinder(measure_mode, new AxesFinderFromAntiKT(0.2), beta, Rcutoff);
         break;
      case onepass_manual_axes:
         setOnePassAxesFinder(measure_mode, new AxesFinderFromUserInput(), beta, Rcutoff);
         break;
      case min_axes: //full minimization is not defined for geometric_measure.
         if (measure_mode == normalized_measure || measure_mode == unnormalized_measure || measure_mode == normalized_cutoff_measure || measure_mode == unnormalized_cutoff_measure)
            //Defaults to 100 iteration to find minimum
            _axesFinder = new AxesFinderFromKmeansMinimization(new AxesFinderFromKT(), beta, Rcutoff, 100);
         else {
            std::cerr << "Multi-pass minimization only set up for normalized_measure, unnormalized_measure, normalized_cutoff_measure, unnormalized_cutoff_measure." << std::endl;
            exit(1);
         }
         break;
      case manual_axes:
         _axesFinder = new AxesFinderFromUserInput();
         break;
// These options have been commented out because they have not been fully tested
//      case wta2_kt_axes: // option for alpha = 2 added
//         _axesFinder = new AxesFinderFromWTA2_KT();
//         break;
//      case wta2_ca_axes: // option for alpha = 2 added
//         _axesFinder = new AxesFinderFromWTA2_CA();
//         break;
//      case onepass_wta2_kt_axes: // option for alpha = 2 added
//         setOnePassAxesFinder(measure_mode, new AxesFinderFromWTA2_KT(), beta, Rcutoff);
//         break;
//      case onepass_wta2_ca_axes: // option for alpha = 2 added
//         setOnePassAxesFinder(measure_mode, new AxesFinderFromWTA2_CA(), beta, Rcutoff);
//         break;
      default:
         assert(false);
         break;
      }   

}

// setAxes for Manual mode
void Njettiness::setAxes(std::vector<fastjet::PseudoJet> myAxes) {
   if (_current_axes_mode == manual_axes || _current_axes_mode == onepass_manual_axes) {
      _currentAxes = myAxes;
   }
   else {
      std::cerr << "You can only use setAxes if using manual_axes or onepass_manual_axes measure mode" << std::endl;
      exit(1);
   }
}
   
// Calculates and returns all TauComponents that user would want.
// This information is stored in _current_tau_components for later access as well.
TauComponents Njettiness::getTauComponents(unsigned n_jets, const std::vector<fastjet::PseudoJet> & inputJets) {
   if (inputJets.size() <= n_jets) {  //if not enough particles, return zero
      _currentAxes = inputJets;
      _currentAxes.resize(n_jets,fastjet::PseudoJet(0.0,0.0,0.0,0.0));
      _current_tau_components = TauComponents();
      _seedAxes = _currentAxes;
   } else {
      _currentAxes = _axesFinder->getAxes(n_jets,inputJets,_currentAxes); // sets current Axes
      _seedAxes = _axesFinder->seedAxes(); // sets seed Axes (if one pass minimization was used)
      _current_tau_components = _measureFunction->result(inputJets, _currentAxes);  // sets current Tau Values
   }
   return _current_tau_components;
}
   
   
// Partition a list of particles according to which N-jettiness axis they are closest to.
// Return a vector of length _currentAxes.size() (which should be N).
// Each vector element is a list of ints corresponding to the indices in
// particles of the particles belonging to that jet.
// TODO:  Consider moving to MeasureFunction
std::vector<std::list<int> > Njettiness::getPartition(const std::vector<fastjet::PseudoJet> & particles) {
   std::vector<std::list<int> > partitions(_currentAxes.size());

   for (unsigned i = 0; i < particles.size(); i++) {
      
      int j_min = -1;
      // find minimum distance
      double minR = std::numeric_limits<double>::max();  //large number
      for (unsigned j = 0; j < _currentAxes.size(); j++) {
         double tempR = _measureFunction->jet_distance_squared(particles[i],_currentAxes[j]); // delta R distance
         if (tempR < minR) {
            minR = tempR;
            j_min = j;
         }
      }
      if (_measureFunction->do_cluster(particles[i],_currentAxes[j_min])) partitions[j_min].push_back(i);
   }
   return partitions;
}

// Having found axes, assign each particle in particles to an axis, and return a set of jets.
// Each jet is the sum of particles closest to an axis (Njet = Naxes).
// TODO:  Consider moving to MeasureFunction
std::vector<fastjet::PseudoJet> Njettiness::getJets(const std::vector<fastjet::PseudoJet> & particles) {
   
   std::vector<fastjet::PseudoJet> jets(_currentAxes.size());

   std::vector<std::list<int> > partition = getPartition(particles);
   for (unsigned j = 0; j < partition.size(); ++j) {
      std::list<int>::const_iterator it, itE;
      for (it = partition[j].begin(), itE = partition[j].end(); it != itE; ++it) {
         jets[j] += particles[*it];
      }
   }
   return jets;
}

} // namespace contrib

FASTJET_END_NAMESPACE
