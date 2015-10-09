//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//  $Id: NjettinessPlugin.hh 822 2015-06-15 23:52:57Z jthaler $
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

#ifndef __FASTJET_CONTRIB_NJETTINESSPLUGIN_HH__
#define __FASTJET_CONTRIB_NJETTINESSPLUGIN_HH__

#include <fastjet/config.h>

#include "Njettiness.hh"
#include "MeasureDefinition.hh"
#include "AxesDefinition.hh"
#include "TauComponents.hh"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"

#include <string>
#include <climits>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


namespace contrib {


/// \class NjettinessPlugin
/// \brief Implements the N-jettiness Jet Algorithm
/**
 * An exclusive jet finder that identifies N jets; first N axes are found, then
 * particles are assigned to the nearest (DeltaR) axis and for each axis the
 * corresponding jet is simply the four-momentum sum of these particles.
 *
 * As of version 2.2, it is recommended to use the XConePlugin, which has
 * sensible default values for jet finding.
 *
 * Axes can be found in several ways, specified by the AxesDefinition argument.
 * The recommended AxesDefinitions for jet finding (different than for jet shapes)
 *    OnePass_AntiKT(R0)  :  one-pass minimization from anti-kT starting point
 *    OnePass_GenET_GenKT_Axes(delta, p, R0) : one-pass min. from GenET/KT
 *    OnePass_WTA_GenKT_Axes(p, R0)          : one-pass min from WTA/GenKT
 * For recommendations on which axes to use, please see the README file.
 * 
 * Jet regions are determined by the MeasureDefinition.  The recommended choices
 * for jet finding are
 *    ConicalMeasure(beta,R0) : perfect cones in rapidity/azimuth plane
 *    XConeMeasure(beta,R0)   : approximate cones based on dot product distances.
 *
 * Other measures introduced in version 2.2 include OriginalGeometricMeasure,
 * ModifiedGeometricMeasure, and ConicalGeometricMeasure, which define N-jettiness
 * through dot products of particle momenta with light-like axes. OriginalGeometricMeasure
 * produces football-shaped jets due to its central weighting of the beam measure,
 * but ModifiedGeometric and ConicalGeometric both deform the original geometric measure
 * to allow for cone-shaped jets. The size of these cones can be controlled through Rcutoff
 * just as in the other measures. See the README file or MeasureDefinition.hh for information
 * on how to call these measures.
 * 
 */
class NjettinessPlugin : public JetDefinition::Plugin {
public:

   /// Constructor with same arguments as Nsubjettiness (N, AxesDefinition, MeasureDefinition)
   NjettinessPlugin(int N,
                    const AxesDefinition & axes_def,
                    const MeasureDefinition & measure_def)
   : _njettinessFinder(axes_def, measure_def), _N(N) {}
   

   /// Description
   virtual std::string description () const;
   /// Jet radius (this does not make sense yet)
   virtual double R() const {return -1.0;} // TODO: make this not stupid

   /// The actually clustering, which first called Njettiness and then creates a dummy ClusterSequence
   virtual void run_clustering(ClusterSequence&) const;

   /// For using manual axes with Njettiness Plugin
   void setAxes(const std::vector<fastjet::PseudoJet> & myAxes) {
      // Cross check that manual axes are being used is in Njettiness
      _njettinessFinder.setAxes(myAxes);
   }

   /// Destructor
   virtual ~NjettinessPlugin() {}

private:

   Njettiness _njettinessFinder;  ///< The core Njettiness that does the heavy lifting
   int _N;  ///< Number of exclusive jets to find.

   /// Warning if the user tries to use v1.0.3 constructor.
   static LimitedWarning _old_constructor_warning;
   
public:

   // Alternative constructors that define the measure via enums and parameters
   // These constructors are deprecated and will be removed in a future version.
   
   /// \deprecated
   /// Old-style constructor with 0 arguments (DEPRECATED)
   NjettinessPlugin(int N,
                    Njettiness::AxesMode axes_mode,
                    Njettiness::MeasureMode measure_mode)
   : _njettinessFinder(axes_mode, measure_mode, 0), _N(N) {
      _old_constructor_warning.warn("NjettinessPlugin:  You are using the old style constructor.  This is deprecated as of v2.1 and will be removed in v3.0.  Please use the NjettinessPlugin constructor based on AxesDefinition and MeasureDefinition instead.");
   }
   
   /// \deprecated
   /// Old-style constructor with 1 argument (DEPRECATED)
      NjettinessPlugin(int N,
                    Njettiness::AxesMode axes_mode,
                    Njettiness::MeasureMode measure_mode,
                    double para1)
   : _njettinessFinder(axes_mode, measure_mode, 1, para1), _N(N) {
      _old_constructor_warning.warn("NjettinessPlugin:  You are using the old style constructor.  This is deprecated as of v2.1 and will be removed in v3.0.  Please use the NjettinessPlugin constructor based on AxesDefinition and MeasureDefinition instead.");
   
   }
   
   /// \deprecated
   /// Old-style constructor with 2 arguments (DEPRECATED)
   NjettinessPlugin(int N,
                    Njettiness::AxesMode axes_mode,
                    Njettiness::MeasureMode measure_mode,
                    double para1,
                    double para2)
   : _njettinessFinder(axes_mode, measure_mode, 2, para1, para2), _N(N) {
      _old_constructor_warning.warn("NjettinessPlugin:  You are using the old style constructor.  This is deprecated as of v2.1 and will be removed in v3.0.  Please use the NjettinessPlugin constructor based on AxesDefinition and MeasureDefinition instead.");
   
   }
   
   /// \deprecated
   /// Old-style constructor with 3 arguments (DEPRECATED)
   NjettinessPlugin(int N,
                    Njettiness::AxesMode axes_mode,
                    Njettiness::MeasureMode measure_mode,
                    double para1,
                    double para2,
                    double para3)
   : _njettinessFinder(axes_mode, measure_mode, 3, para1, para2, para3), _N(N) {
      _old_constructor_warning.warn("NjettinessPlugin:  You are using the old style constructor.  This is deprecated as of v2.1 and will be removed in v3.0.  Please use the NjettinessPlugin constructor based on AxesDefinition and MeasureDefinition instead.");
   }
   
   /// \deprecated
   /// Old-style constructor for backwards compatibility with v1.0, when NormalizedCutoffMeasure was the only option
   NjettinessPlugin(int N,
                    Njettiness::AxesMode mode,
                    double beta,
                    double R0,
                    double Rcutoff=std::numeric_limits<double>::max())
   : _njettinessFinder(mode, NormalizedCutoffMeasure(beta, R0, Rcutoff)), _N(N) {
      _old_constructor_warning.warn("NjettinessPlugin:  You are using the old style constructor.  This is deprecated as of v2.1 and will be removed in v3.0.  Please use the NjettinessPlugin constructor based on AxesDefinition and MeasureDefinition instead.");
   }


};

} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_NJETTINESSPLUGIN_HH__