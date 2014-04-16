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


#include "MeasureFunction.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

///////
//
// Measure Function
//
///////

// Return all of the necessary TauComponents for specific input particles and axes
TauComponents MeasureFunction::result(const std::vector<fastjet::PseudoJet>& particles, const std::vector<fastjet::PseudoJet>& axes) {

   std::vector<double> jetPieces(axes.size(), 0.0);
   double beamPiece = 0.0;
   
   // Calculates the unnormalized sub-tau values, i.e. a std::vector of the contributions to tau_N of each Voronoi region (or region within R_0)
   for (unsigned i = 0; i < particles.size(); i++) {
      
      // find minimum distance; start with beam (-1) for reference
      int j_min = -1;
      double minRsq;
      if (_has_beam) minRsq = beam_distance_squared(particles[i]);
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
         if (_has_beam) beamPiece += beam_numerator(particles[i]);
         else assert(_has_beam);  // this should never happen.
      } else {
         jetPieces[j_min] += jet_numerator(particles[i],axes[j_min]);
      }
   }
   
   // Calculates normalization for tau and subTau if _has_denominator is true, otherwise returns 1.0 (i.e. no normalization)
   double tauDen = 0.0;
   if (_has_denominator) {
      for (unsigned i = 0; i < particles.size(); i++) {
         tauDen += denominator(particles[i]);
      }
   } else {
      tauDen = 1.0; // if no denominator, then 1.0 for no normalization factor
   }
   
   return TauComponents(jetPieces, beamPiece, tauDen, _has_denominator, _has_beam);
}
   
} //namespace contrib

FASTJET_END_NAMESPACE
