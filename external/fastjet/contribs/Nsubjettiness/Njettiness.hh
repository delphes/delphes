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

#ifndef __FASTJET_CONTRIB_NJETTINESS_HH__
#define __FASTJET_CONTRIB_NJETTINESS_HH__

#include "MeasureFunction.hh"
#include "AxesFinder.hh"

#include "fastjet/PseudoJet.hh"
#include <cmath>
#include <vector>
#include <list>


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {

///////
//
// Main Njettiness Class
//
///////

//------------------------------------------------------------------------
/// \class Njettiness
// Njettiness uses AxesFinder and MeasureFunction together in order to find tau_N for the event. The user specifies
// which AxesFinder and which MeasureFunction to use in the calculation, and then Njettiness returns tau_N for the event.
// It also can return information about the axes and jets it used in the calculation, as well as information about 
// how the event was partitioned.
class Njettiness {
public:
   
   // The various axes choices available to the user
   enum AxesMode {
      kt_axes,             // exclusive kt axes
      ca_axes,             // exclusive ca axes
      antikt_0p2_axes,     // inclusive hardest axes with antikt-0.2
      wta_kt_axes,         // Winner Take All axes with kt
      wta_ca_axes,         // Winner Take All axes with CA
      onepass_kt_axes,     // one-pass minimization from kt starting point
      onepass_ca_axes,     // one-pass minimization from ca starting point
      onepass_antikt_0p2_axes,  // one-pass minimization from antikt-0.2 starting point
      onepass_wta_kt_axes, //one-pass minimization of WTA axes with kt 
      onepass_wta_ca_axes, //one-pass minimization of WTA axes with ca 
      min_axes,            // axes that minimize N-subjettiness (100 passes by default)
      manual_axes,         // set your own axes with setAxes()
      onepass_manual_axes  // one-pass minimization from manual starting point
      //  These options are commented out because they have not been fully tested
      //      wta2_kt_axes,        // Winner Take All (alpha = 2) with kt
      //      wta2_ca_axes,         // Winner Take All (alpha = 2) with CA
      //      onepass_wta2_kt_axes, //one-pass minimization of WTA (alpha = 2) axes with kt
      //      onepass_wta2_ca_axes, //one-pass minimization of WTA (alpha = 2) axes with ca
   };

   // The measures available to the user.
   // "normalized_cutoff_measure" was the default in v1.0 of Nsubjettiness
   // "unnormalized_measure" is now the recommended default usage
   enum MeasureMode {
      normalized_measure,           //default normalized measure
      unnormalized_measure,         //default unnormalized measure
      geometric_measure,            //geometric measure
      normalized_cutoff_measure,    //default normalized measure with explicit Rcutoff
      unnormalized_cutoff_measure,  //default unnormalized measure with explicit Rcutoff
      geometric_cutoff_measure      //geometric measure with explicit Rcutoff
   };

private:
   // The chosen axes/measure modes
   AxesFinder* _axesFinder;  // The chosen axes
   MeasureFunction* _measureFunction; // The chosen measure

   // Enum information so functions can specify output based on specific options, primarily for setAxes
   AxesMode _current_axes_mode;
   MeasureMode _current_measure_mode;
   
   // Information about the current information
   TauComponents _current_tau_components; //automatically set to have components of 0; these values will be set by the getTau function call
   std::vector<fastjet::PseudoJet> _currentAxes;
   std::vector<fastjet::PseudoJet> _seedAxes; // axes used prior to minimization (if applicable)
   
   // Needed for compilation of non C++11 users
   bool isnan(double para) { return para != para; }

   // Helpful function to check to make sure input has correct number of parameters
   bool correctParameterCount(int n, double para1, double para2, double para3, double para4){
      int numpara;
      if (!isnan(para1) && !isnan(para2) && !isnan(para3) && !isnan(para4)) numpara = 4;
      else if (!isnan(para1) && !isnan(para2) && !isnan(para3) && isnan(para4)) numpara = 3;
      else if (!isnan(para1) && !isnan(para2) && isnan(para3) && isnan(para4)) numpara = 2;
      else if (!isnan(para1) && isnan(para2) && isnan(para3) && isnan(para4)) numpara = 1;
      else numpara = 0;
      return n == numpara;
   }

   // Helper function to set onepass_axes depending on input measure_mode and startingFinder
   void setOnePassAxesFinder(MeasureMode measure_mode, AxesFinder* startingFinder, double para1, double Rcutoff);
 
   // created separate function to set MeasureFunction and AxesFinder in order to keep constructor cleaner.
   void setMeasureFunctionandAxesFinder(AxesMode axes_mode, MeasureMode measure_mode, double para1, double para2, double para3, double para4);

public:


   // Main constructor which takes axes/measure information, and possible parameters.
   // Unlike Nsubjettiness or NjettinessPlugin, the value N is not chosen
   Njettiness(AxesMode axes_mode,
              MeasureMode measure_mode,
              double para1 = NAN,
              double para2 = NAN,
              double para3 = NAN,
              double para4 = NAN)
   : _current_axes_mode(axes_mode),
   _current_measure_mode(measure_mode) {
      setMeasureFunctionandAxesFinder(axes_mode, measure_mode, para1, para2, para3, para4);  // call helper function to do the hard work
   }

   ~Njettiness() {
      // clean house
      delete _measureFunction;
      delete _axesFinder;
   }
   
   // setAxes for Manual mode
   void setAxes(std::vector<fastjet::PseudoJet> myAxes);
   
   // Calculates and returns all TauComponents that user would want.
   // This information is stored in _current_tau_components for later access as well.
   TauComponents getTauComponents(unsigned n_jets, const std::vector<fastjet::PseudoJet> & inputJets);

   // Calculates the value of N-subjettiness,
   // but only returns the tau value from _current_tau_components
   double getTau(unsigned n_jets, const std::vector<fastjet::PseudoJet> & inputJets) {
      return getTauComponents(n_jets, inputJets).tau();
   }
   
   // returns enum information
   MeasureMode currentMeasureMode() { return _current_measure_mode;}
   AxesMode currentAxesMode() { return _current_axes_mode;}

   // Return all relevant information about tau components
   TauComponents currentTauComponents() {return _current_tau_components;}
   
   // Return axes found by getTauComponents.
   std::vector<fastjet::PseudoJet> currentAxes() { return _currentAxes;}
   // Return seedAxes used if onepass minimization (otherwise, same as currentAxes)
   std::vector<fastjet::PseudoJet> seedAxes() { return _seedAxes;}
   
   // partition inputs by Voronoi (each vector stores indices corresponding to inputJets)
   std::vector<std::list<int> > getPartition(const std::vector<fastjet::PseudoJet> & inputJets);

   // partition inputs by Voronoi
   std::vector<fastjet::PseudoJet> getJets(const std::vector<fastjet::PseudoJet> & inputJets);

};

} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_NJETTINESS_HH__

