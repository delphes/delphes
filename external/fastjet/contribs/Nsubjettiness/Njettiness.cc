//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//  $Id: Njettiness.cc 677 2014-06-12 18:56:46Z jthaler $
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

Njettiness::Njettiness(const AxesDefinition & axes_def, const MeasureDefinition & measure_def)
: _axes_def(axes_def.create()), _measure_def(measure_def.create()) {
   setMeasureFunctionAndAxesFinder();  // call helper function to do the hard work
}

Njettiness::Njettiness(AxesMode axes_mode, const MeasureDefinition & measure_def)
: _axes_def(createAxesDef(axes_mode)), _measure_def(measure_def.create()) {
   setMeasureFunctionAndAxesFinder();  // call helper function to do the hard work
}
   
// Convert from MeasureMode enum to MeasureDefinition
// This returns a pointer that will be claimed by a SharedPtr
MeasureDefinition* Njettiness::createMeasureDef(MeasureMode measure_mode, int num_para, double para1, double para2, double para3) const {

   // definition of maximum Rcutoff for non-cutoff measures, changed later by other measures
   double Rcutoff = std::numeric_limits<double>::max();  //large number
   // Most (but not all) measures have some kind of beta value
   double beta = std::numeric_limits<double>::quiet_NaN();
   // The normalized measures have an R0 value.
   double R0 = std::numeric_limits<double>::quiet_NaN();
   
   // Find the MeasureFunction and set the parameters.
   switch (measure_mode) {
      case normalized_measure:
         beta = para1;
         R0 = para2;
         if(num_para == 2) {
            return new NormalizedMeasure(beta,R0);
         } else {
            throw Error("normalized_measure needs 2 parameters (beta and R0)");
         }
         break;
      case unnormalized_measure:
         beta = para1;
         if(num_para == 1) {
            return new UnnormalizedMeasure(beta);
         } else {
            throw Error("unnormalized_measure needs 1 parameter (beta)");
         }
         break;
      case geometric_measure:
         beta = para1;
         if (num_para == 1) {
            return new GeometricMeasure(beta);
         } else {
            throw Error("geometric_measure needs 1 parameter (beta)");
         }
         break;
      case normalized_cutoff_measure:
         beta = para1;
         R0 = para2;
         Rcutoff = para3; //Rcutoff parameter is 3rd parameter in normalized_cutoff_measure
         if (num_para == 3) {
            return new NormalizedCutoffMeasure(beta,R0,Rcutoff);
         } else {
            throw Error("normalized_cutoff_measure has 3 parameters (beta, R0, Rcutoff)");
         }
         break;
      case unnormalized_cutoff_measure:
         beta = para1;
         Rcutoff = para2; //Rcutoff parameter is 2nd parameter in normalized_cutoff_measure
         if (num_para == 2) {
            return new UnnormalizedCutoffMeasure(beta,Rcutoff);
         } else {
            throw Error("unnormalized_cutoff_measure has 2 parameters (beta, Rcutoff)");
         }
         break;
      case geometric_cutoff_measure:
         beta = para1;
         Rcutoff = para2; //Rcutoff parameter is 2nd parameter in geometric_cutoff_measure
         if(num_para == 2) {
           return new GeometricCutoffMeasure(beta,Rcutoff);
         } else {
            throw Error("geometric_cutoff_measure has 2 parameters (beta, Rcutoff)");
         }
         break;
      default:
         assert(false);
         break;
   }
   return NULL;
}

// Convert from AxesMode enum to AxesDefinition
// This returns a pointer that will be claimed by a SharedPtr
AxesDefinition* Njettiness::createAxesDef(Njettiness::AxesMode axes_mode) const {
   
   switch (axes_mode) {
      case wta_kt_axes:
         return new WTA_KT_Axes();
      case wta_ca_axes:
         return new WTA_CA_Axes();
      case kt_axes:
         return new KT_Axes();
      case ca_axes:
         return new CA_Axes();
      case antikt_0p2_axes:
         return new AntiKT_Axes(0.2);
      case onepass_wta_kt_axes:
         return new OnePass_WTA_KT_Axes();
      case onepass_wta_ca_axes:
         return new OnePass_WTA_CA_Axes();
      case onepass_kt_axes:
         return new OnePass_KT_Axes();
      case onepass_ca_axes:
         return new OnePass_CA_Axes();
      case onepass_antikt_0p2_axes:
         return new OnePass_AntiKT_Axes(0.2);
      case onepass_manual_axes:
         return new OnePass_Manual_Axes();
      case min_axes:
         return new MultiPass_Axes(100);
      case manual_axes:
         return new Manual_Axes();
      default:
         assert(false);
         return NULL;
   }
}

   
// Parsing needed for constructor to set AxesFinder and MeasureFunction
// All of the parameter handling is here, and checking that number of parameters is correct.
void Njettiness::setMeasureFunctionAndAxesFinder() {
   // Get the correct MeasureFunction and AxesFinders
   _measureFunction.reset(_measure_def->createMeasureFunction());
   _startingAxesFinder.reset(_axes_def->createStartingAxesFinder(*_measure_def));
   _finishingAxesFinder.reset(_axes_def->createFinishingAxesFinder(*_measure_def));
}

// setAxes for Manual mode
void Njettiness::setAxes(const std::vector<fastjet::PseudoJet> & myAxes) {
   if (_axes_def->supportsManualAxes()) {
      _currentAxes = myAxes;
   } else {
      throw Error("You can only use setAxes for manual AxesDefinitions");
   }
}
   
// Calculates and returns all TauComponents that user would want.
// This information is stored in _current_tau_components for later access as well.
TauComponents Njettiness::getTauComponents(unsigned n_jets, const std::vector<fastjet::PseudoJet> & inputJets) const {
   if (inputJets.size() <= n_jets) {  //if not enough particles, return zero
      _currentAxes = inputJets;
      _currentAxes.resize(n_jets,fastjet::PseudoJet(0.0,0.0,0.0,0.0));
      _current_tau_components = TauComponents();
      _seedAxes = _currentAxes;
      _currentJets = _currentAxes;
      _currentBeam = PseudoJet(0.0,0.0,0.0,0.0);
   } else {

      _seedAxes = _startingAxesFinder->getAxes(n_jets,inputJets,_currentAxes); //sets starting point for minimization
      if (_finishingAxesFinder) {
         _currentAxes = _finishingAxesFinder->getAxes(n_jets,inputJets,_seedAxes);
      } else {
         _currentAxes = _seedAxes;
      }
      
      // Find partition and store information
      // (jet information in _currentJets, beam in _currentBeam)
      _currentJets = _measureFunction->get_partition(inputJets,_currentAxes,&_currentBeam);
      
      // Find tau value and store information
      _current_tau_components = _measureFunction->result_from_partition(_currentJets, _currentAxes,&_currentBeam);  // sets current Tau Values
   }
   return _current_tau_components;
}
   
   
// Partition a list of particles according to which N-jettiness axis they are closest to.
// Return a vector of length _currentAxes.size() (which should be N).
// Each vector element is a list of ints corresponding to the indices in
// particles of the particles belonging to that jet.
std::vector<std::list<int> > Njettiness::getPartitionList(const std::vector<fastjet::PseudoJet> & particles) const {
   // core code is in MeasureFunction
   return _measureFunction->get_partition_list(particles,_currentAxes);
}

   
} // namespace contrib

FASTJET_END_NAMESPACE
