// Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//  $Id: Nsubjettiness.hh 822 2015-06-15 23:52:57Z jthaler $
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

#ifndef __FASTJET_CONTRIB_NSUBJETTINESS_HH__
#define __FASTJET_CONTRIB_NSUBJETTINESS_HH__

#include <fastjet/internal/base.hh>

#include "Njettiness.hh"

#include "fastjet/FunctionOfPseudoJet.hh"
#include <string>
#include <climits>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {

// Classes defined in this file.
class Nsubjettiness;
class NsubjettinessRatio;
   
///------------------------------------------------------------------------
/// \class Nsubjettiness
/// \brief Implements the N-subjettiness jet shape
///
/**
 * The N-jettiness jet shape.
 *
 * Nsubjettiness extends the concept of Njettiness to a jet shape, but other
 * than the set of particles considered, they are identical.  This class
 * wraps the core Njettiness code to provide the fastjet::FunctionOfPseudoJet
 * interface for convenience in larger analyses.
 *
 * The recommended AxesDefinitions are:
 *   KT_Axes             : exclusive kt axes
 *   WTA_KT_Axes         : exclusive kt with winner-take-all recombination
 *   OnePass_KT_Axes     : one-pass minimization from kt starting point
 *   OnePass_WTA_KT_Axes : one-pass min. from wta_kt starting point
 * More AxesDefinitions are listed in the README and defined in AxesDefinition.hh
 *
 * The recommended MeasureDefinitions are (with the corresponding parameters)
 *   NormalizedMeasure(beta,R0)
 *      :  This was the original N-subjettiness measure (dimensionless)
 *   UnnormalizedMeasure(beta)
 *      :  This is the new recommended default, same as above but without
 *      :  the normalization factor, and hence has units of GeV
 *   NormalizedCutoffMeasure(beta,R0,Rcutoff)
 *      :  Same as normalized_measure, but cuts off at Rcutoff
 *   UnnormalizedCutoffMeasure(beta,Rcutoff)
 *      :  Same as unnormalized_measure, but cuts off at Rcutoff
 * More MeasureDefinitions are listed in the README and defined in MeasureDefinition.hh
 *
 * For example, for the UnnormalizedMeasure(beta), N-subjettiness is defined as:
 *
 * tau_N = Sum_{i in jet} p_T^i min((DR_i1)^beta, (DR_i2)^beta, ...)
 *
 *   DR_ij is the distance sqrt(Delta_phi^2 + Delta_rap^2) between particle i
 *   and jet j.
 *
 * The NormalizedMeausure include an extra parameter R0, and the various cutoff
 * measures include an Rcutoff, which effectively defines an angular cutoff
 * similar in effect to a cone-jet radius.
 */
class Nsubjettiness : public FunctionOfPseudoJet<double> {

public:
   
   /// Main constructor, which takes N, the AxesDefinition, and the MeasureDefinition.
   /// The Definitions are given in AxesDefinition.hh and MeasureDefinition.hh
   Nsubjettiness(int N,
                 const AxesDefinition& axes_def,
                 const MeasureDefinition& measure_def)
   : _njettinessFinder(axes_def,measure_def), _N(N) {}
   
   /// Returns tau_N, measured on the constituents of this jet
   double result(const PseudoJet& jet) const;

   /// Returns components of tau_N, so that user can find individual tau values.
   TauComponents component_result(const PseudoJet& jet) const;
   
   /// To set axes in manual mode
   void setAxes(const std::vector<fastjet::PseudoJet> & myAxes) {
      // Note that cross check that manual mode has been set is in Njettiness
    	_njettinessFinder.setAxes(myAxes);
   }
   
   /// returns seed axes used for onepass minimization (otherwise same as currentAxes)
   std::vector<fastjet::PseudoJet> seedAxes() const {
      return _njettinessFinder.seedAxes();
   }
   
   /// returns current axes found by result() calculation
   std::vector<fastjet::PseudoJet> currentAxes() const {
      return _njettinessFinder.currentAxes();
   }

   /// returns subjet regions found by result() calculation (these have valid constituents)
   /// Note that the axes and the subjets are not the same
   std::vector<fastjet::PseudoJet> currentSubjets() const {
      return _njettinessFinder.currentJets();
   }

   /// returns components of tau_N without recalculating anything
   TauComponents currentTauComponents() const {
      return _njettinessFinder.currentTauComponents();
   }

   /// returns components of tau_N without recalculating anything
   TauPartition currentPartition() const {
      return _njettinessFinder.currentPartition();
   }
   
private:
   
   /// Core Njettiness code that is called
   Njettiness _njettinessFinder; // TODO:  should muck with this so result can be const without this mutable
   /// Number of subjets to find
   int _N;
   
   /// Warning if the user tries to use v1.0.3 constructor.
   static LimitedWarning _old_constructor_warning;

public:
   
   // The following interfaces are included for backwards compatibility, but no longer recommended.
   // They may be deleted in a future release
   
   /// \deprecated
   /// Alternative constructors that define the measure via enums and parameters
   /// These constructors will be removed in v3.0
   /// Zero parameter arguments
   /// (Currently, no measure uses this)
   Nsubjettiness(int N,
                 Njettiness::AxesMode axes_mode,
                 Njettiness::MeasureMode measure_mode)
   : _njettinessFinder(axes_mode, measure_mode, 0), _N(N) {
      _old_constructor_warning.warn("Nsubjettiness:  You are using the old style constructor.  This is deprecated as of v2.1 and will be removed in v3.0.  Please use the Nsubjettiness constructor based on AxesDefinition and MeasureDefinition instead.");
   }
   
   /// \deprecated
   /// Construcotr for one parameter argument
   /// (for unnormalized_measure, para1=beta)
   Nsubjettiness(int N,
                 Njettiness::AxesMode axes_mode,
                 Njettiness::MeasureMode measure_mode,
                 double para1)
   : _njettinessFinder(axes_mode, measure_mode, 1, para1), _N(N) {
      _old_constructor_warning.warn("Nsubjettiness:  You are using the old style constructor.  This is deprecated as of v2.1 and will be removed in v3.0.  Please use the Nsubjettiness constructor based on AxesDefinition and MeasureDefinition instead.");
   }
   
   /// \deprecated
   /// Constructor for two parameter arguments
   /// (for normalized_measure, para1=beta, para2=R0)
   /// (for unnormalized_cutoff_measure, para1=beta, para2=Rcutoff)
   Nsubjettiness(int N,
                 Njettiness::AxesMode axes_mode,
                 Njettiness::MeasureMode measure_mode,
                 double para1,
                 double para2)
   : _njettinessFinder(axes_mode, measure_mode, 2, para1, para2), _N(N) {
      _old_constructor_warning.warn("Nsubjettiness:  You are using the old style constructor.  This is deprecated as of v2.1 and will be removed in v3.0.  Please use the Nsubjettiness constructor based on AxesDefinition and MeasureDefinition instead.");
   }
   
   /// \deprecated
   /// Constructor for three parameter arguments
   /// (for unnormalized_cutoff_measure, para1=beta, para2=R0, para3=Rcutoff)
   Nsubjettiness(int N,
                 Njettiness::AxesMode axes_mode,
                 Njettiness::MeasureMode measure_mode,
                 double para1,
                 double para2,
                 double para3)
   : _njettinessFinder(axes_mode, measure_mode, 3, para1, para2, para3), _N(N) {
      _old_constructor_warning.warn("Nsubjettiness:  You are using the old style constructor.  This is deprecated as of v2.1 and will be removed in v3.0.  Please use the Nsubjettiness constructor based on AxesDefinition and MeasureDefinition instead.");
   }
   
   /// \deprecated
   /// Old constructor for backwards compatibility with v1.0,
   /// where normalized_cutoff_measure was the only option
   Nsubjettiness(int N,
                 Njettiness::AxesMode axes_mode,
                 double beta,
                 double R0,
                 double Rcutoff=std::numeric_limits<double>::max())
   : _njettinessFinder(axes_mode, NormalizedCutoffMeasure(beta,R0,Rcutoff)), _N(N) {
      _old_constructor_warning.warn("Nsubjettiness:  You are using the old style constructor.  This is deprecated as of v2.1 and will be removed in v3.0.  Please use the Nsubjettiness constructor based on AxesDefinition and MeasureDefinition instead.");
   }
   
};


///------------------------------------------------------------------------
/// \class NsubjettinessRatio
/// \brief Implements ratios of N-subjettiness jet shapes
///
/// NsubjettinessRatio uses the results from Nsubjettiness to calculate the ratio
/// tau_N/tau_M, where N and M are specified by the user. The ratio of different tau values
/// is often used in analyses, so this class is helpful to streamline code.  Note that
/// manual axis mode is not supported
class NsubjettinessRatio : public FunctionOfPseudoJet<double> {
public:

   /// Main constructor.  Apart from specifying both N and M, the same options as Nsubjettiness
   NsubjettinessRatio(int N,
                      int M,
                      const AxesDefinition & axes_def,
                      const MeasureDefinition & measure_def)
   : _nsub_numerator(N,axes_def,measure_def),
   _nsub_denominator(M,axes_def,measure_def) {
      if (axes_def.needsManualAxes()) {
         throw Error("NsubjettinessRatio does not support ManualAxes mode.");
      }
   }

   /// Returns tau_N/tau_M based off the input jet using result function from Nsubjettiness
   double result(const PseudoJet& jet) const;

private: 

   Nsubjettiness _nsub_numerator;   ///< Function for numerator
   Nsubjettiness _nsub_denominator; ///< Function for denominator

public:
   
   // The following interfaces are included for backwards compatibility, but no longer recommended.
   // They may be deprecated at some point.
   
   /// \deprecated
   /// Old-style constructor for zero arguments
   /// Alternative constructor with enums and parameters
   /// Again, likely to be removed
   NsubjettinessRatio(int N,
                      int M,
                      Njettiness::AxesMode axes_mode,
                      Njettiness::MeasureMode measure_mode)
   : _nsub_numerator(N, axes_mode, measure_mode),
   _nsub_denominator(M, axes_mode, measure_mode) {}
   
   /// \deprecated
   /// Old-style constructor for one argument
   NsubjettinessRatio(int N,
                      int M,
                      Njettiness::AxesMode axes_mode,
                      Njettiness::MeasureMode measure_mode,
                      double para1)
   : _nsub_numerator(N, axes_mode, measure_mode, para1),
   _nsub_denominator(M, axes_mode, measure_mode, para1) {}
   
   /// \deprecated
   /// Old-style constructor for 2 arguments
   NsubjettinessRatio(int N,
                      int M,
                      Njettiness::AxesMode axes_mode,
                      Njettiness::MeasureMode measure_mode,
                      double para1,
                      double para2)
   : _nsub_numerator(N, axes_mode, measure_mode, para1, para2),
   _nsub_denominator(M, axes_mode, measure_mode, para1, para2) {}
   
   /// \deprecated
   /// Old-style constructor for 3 arguments
   NsubjettinessRatio(int N,
                      int M,
                      Njettiness::AxesMode axes_mode,
                      Njettiness::MeasureMode measure_mode,
                      double para1,
                      double para2,
                      double para3)
   : _nsub_numerator(N, axes_mode, measure_mode, para1, para2, para3),
   _nsub_denominator(M, axes_mode, measure_mode, para1, para2, para3) {}

   
};

} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_NSUBJETTINESS_HH__
