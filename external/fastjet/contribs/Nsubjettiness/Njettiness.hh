//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//  $Id: Njettiness.hh 670 2014-06-06 01:24:42Z jthaler $
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
#include "NjettinessDefinition.hh"

#include "fastjet/PseudoJet.hh"
#include "fastjet/SharedPtr.hh"
#include <fastjet/LimitedWarning.hh>

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
   // It is recommended to use AxesDefinition instead of these.
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
   // But it is recommended to use MeasureDefinition instead of these.
   enum MeasureMode {
      normalized_measure,           //default normalized measure
      unnormalized_measure,         //default unnormalized measure
      geometric_measure,            //geometric measure
      normalized_cutoff_measure,    //default normalized measure with explicit Rcutoff
      unnormalized_cutoff_measure,  //default unnormalized measure with explicit Rcutoff
      geometric_cutoff_measure      //geometric measure with explicit Rcutoff
   };

   // Main constructor that uses AxesMode and MeasureDefinition to specify measure
   // Unlike Nsubjettiness or NjettinessPlugin, the value N is not chosen
   Njettiness(const AxesDefinition & axes_def, const MeasureDefinition & measure_def);

   // Intermediate constructor (needed to enable v1.0.3 backwards compatibility?)
   Njettiness(AxesMode axes_mode, const MeasureDefinition & measure_def);

   // Alternative constructor which takes axes/measure information as enums with measure parameters
   // This version is not recommended
   Njettiness(AxesMode axes_mode,
              MeasureMode measure_mode,
              int num_para,
              double para1 = std::numeric_limits<double>::quiet_NaN(),
              double para2 = std::numeric_limits<double>::quiet_NaN(),
              double para3 = std::numeric_limits<double>::quiet_NaN())
   : _axes_def(createAxesDef(axes_mode)), _measure_def(createMeasureDef(measure_mode, num_para, para1, para2, para3)) {
      setMeasureFunctionAndAxesFinder();  // call helper function to do the hard work
   }

   // destructor
   ~Njettiness() {};
   
   // setAxes for Manual mode
   void setAxes(const std::vector<fastjet::PseudoJet> & myAxes);
   
   // Calculates and returns all TauComponents that user would want.
   // This information is stored in _current_tau_components for later access as well.
   TauComponents getTauComponents(unsigned n_jets, const std::vector<fastjet::PseudoJet> & inputJets) const;

   // Calculates the value of N-subjettiness,
   // but only returns the tau value from _current_tau_components
   double getTau(unsigned n_jets, const std::vector<fastjet::PseudoJet> & inputJets) const {
      return getTauComponents(n_jets, inputJets).tau();
   }

   // Return all relevant information about tau components
   TauComponents currentTauComponents() const {return _current_tau_components;}
   // Return axes found by getTauComponents.
   std::vector<fastjet::PseudoJet> currentAxes() const { return _currentAxes;}
   // Return seedAxes used if onepass minimization (otherwise, same as currentAxes)
   std::vector<fastjet::PseudoJet> seedAxes() const { return _seedAxes;}
   // Return jet partition found by getTauComponents.
   std::vector<fastjet::PseudoJet> currentJets() const {return _currentJets;}
   // Return beam partition found by getTauComponents.
   fastjet::PseudoJet currentBeam() const {return _currentBeam;}
   
   // partition inputs by Voronoi (each vector stores indices corresponding to inputJets)
   std::vector<std::list<int> > getPartitionList(const std::vector<fastjet::PseudoJet> & inputJets) const;

private:
   
   // Information about Axes and Measures to be Used
   // Implemented as SharedPtrs to avoid memory management headaches
   SharedPtr<const AxesDefinition> _axes_def;
   SharedPtr<const MeasureDefinition> _measure_def;
   
   // The chosen axes/measure mode workers
   // Implemented as SharedPtrs to avoid memory management headaches
   // TODO: make into a SharedPtr<const AxesFinder>?
   SharedPtr<MeasureFunction> _measureFunction;  // The chosen measure
   SharedPtr<AxesFinder> _startingAxesFinder;    // The initial axes finder
   SharedPtr<AxesFinder> _finishingAxesFinder;   // A possible minimization step
   
   // Information about the current information
   // Defined as mutables, so user should be aware that these change when getTau is called.
   mutable TauComponents _current_tau_components; //automatically set to have components of 0; these values will be set by the getTau function call
   mutable std::vector<fastjet::PseudoJet> _currentAxes; //axes found after minimization
   mutable std::vector<fastjet::PseudoJet> _seedAxes; // axes used prior to minimization (if applicable)
   mutable std::vector<fastjet::PseudoJet> _currentJets; //partitioning information
   mutable fastjet::PseudoJet _currentBeam; //return beam, if requested
   
   // created separate function to set MeasureFunction and AxesFinder in order to keep constructor cleaner.
   void setMeasureFunctionAndAxesFinder();
   
   // Convert old style enums into new style MeasureDefinition
   AxesDefinition* createAxesDef(AxesMode axes_mode) const;
   
   // Convert old style enums into new style MeasureDefinition
   MeasureDefinition* createMeasureDef(MeasureMode measure_mode, int num_para, double para1, double para2, double para3) const;

};
   
} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_NJETTINESS_HH__

