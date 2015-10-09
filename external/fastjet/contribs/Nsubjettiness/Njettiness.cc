//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//  $Id: Njettiness.cc 821 2015-06-15 18:50:53Z jthaler $
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

LimitedWarning Njettiness::_old_measure_warning;
LimitedWarning Njettiness::_old_axes_warning;

   
// Constructor
Njettiness::Njettiness(const AxesDefinition & axes_def, const MeasureDefinition & measure_def)
: _axes_def(axes_def.create()), _measure_def(measure_def.create()) {}
   
// setAxes for Manual mode
void Njettiness::setAxes(const std::vector<fastjet::PseudoJet> & myAxes) {
   if (_axes_def()->needsManualAxes()) {
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
      _currentPartition = TauPartition(n_jets); // empty partition
   } else {
      assert(_axes_def()); // this should never fail.
      
      if (_axes_def()->needsManualAxes()) { // if manual mode
         // take current axes as seeds
         _seedAxes = _currentAxes;
         
         // refine axes if requested
         _currentAxes = _axes_def->get_refined_axes(n_jets,inputJets,_seedAxes, _measure_def());
      } else { // non-manual axes
         
          //set starting point for minimization
         _seedAxes = _axes_def->get_starting_axes(n_jets,inputJets,_measure_def());
         
         // refine axes as needed
         _currentAxes = _axes_def->get_refined_axes(n_jets,inputJets,_seedAxes, _measure_def());
         
         // NOTE:  The above two function calls are combined in "AxesDefinition::get_axes"
         // but are separated here to allow seed axes to be stored.
      }
      
      // Find and store partition
      _currentPartition = _measure_def->get_partition(inputJets,_currentAxes);
      
      // Find and store tau value
      _current_tau_components = _measure_def->component_result_from_partition(_currentPartition, _currentAxes);  // sets current Tau Values
   }
   return _current_tau_components;
}
   

///////
//
// Below is code for backward compatibility to use the old interface.
// May be deleted in a future version
//
///////
   
Njettiness::Njettiness(AxesMode axes_mode, const MeasureDefinition & measure_def)
: _axes_def(createAxesDef(axes_mode)), _measure_def(measure_def.create()) {}

// Convert from MeasureMode enum to MeasureDefinition
// This returns a pointer that will be claimed by a SharedPtr
MeasureDefinition* Njettiness::createMeasureDef(MeasureMode measure_mode, int num_para, double para1, double para2, double para3) const {
   
   _old_measure_warning.warn("Njettiness::createMeasureDef:  You are using the old MeasureMode way of specifying N-subjettiness measures.  This is deprecated as of v2.1 and will be removed in v3.0.  Please use MeasureDefinition instead.");
   
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
         throw Error("This class has been removed. Please use OriginalGeometricMeasure, ModifiedGeometricMeasure, or ConicalGeometricMeasure with the new Njettiness constructor.");
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
         throw Error("This class has been removed. Please use OriginalGeometricMeasure, ModifiedGeometricMeasure, or ConicalGeometricMeasure with the new Njettiness constructor.");
      default:
         assert(false);
         break;
   }
   return NULL;
}

// Convert from AxesMode enum to AxesDefinition
// This returns a pointer that will be claimed by a SharedPtr
AxesDefinition* Njettiness::createAxesDef(Njettiness::AxesMode axes_mode) const {
   
   _old_axes_warning.warn("Njettiness::createAxesDef:  You are using the old AxesMode way of specifying N-subjettiness axes.  This is deprecated as of v2.1 and will be removed in v3.0.  Please use AxesDefinition instead.");
   
   
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


   
   
   
} // namespace contrib

FASTJET_END_NAMESPACE
