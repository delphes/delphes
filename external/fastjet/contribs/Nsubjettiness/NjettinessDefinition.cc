//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//  $Id: NjettinessDefinition.cc 704 2014-07-07 14:30:43Z jthaler $
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

#include "NjettinessDefinition.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib {

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
   
   std::string GeometricMeasure::description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "Geometric Measure (beta = " << _beta << ", in GeV)";
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
   
   std::string GeometricCutoffMeasure::description() const {
      std::stringstream stream;
      stream << std::fixed << std::setprecision(2)
      << "Geometric Cutoff Measure (beta = " << _beta << ", Rcut = " << _Rcutoff << ", in GeV)";
      return stream.str();
   };
   
   
} // namespace contrib

FASTJET_END_NAMESPACE
