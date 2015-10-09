//  Nsubjettiness Package
//  Questions/Comments?  jthaler@jthaler.net
//
//  Copyright (c) 2011-14
//  Jesse Thaler, Ken Van Tilburg, Christopher K. Vermilion, and TJ Wilkason
//
//  $Id: XConePlugin.cc 745 2014-08-26 23:51:48Z jthaler $
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

#include "XConePlugin.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

std::string XConePlugin::description() const {
   std::stringstream stream;
   stream << "XCone Jet Algorithm with N = " << _N << std::fixed << std::setprecision(2) << ", Rcut = " << _R0 << ", beta = " << _beta;
   return stream.str();
}

std::string PseudoXConePlugin::description() const {
   std::stringstream stream;
   stream
   << "PseudoXCone Jet Algorithm with N = " << _N << std::fixed << std::setprecision(2) << ", Rcut = " << _R0 << ", beta = " << _beta;
   return stream.str();
}
   
} // namespace contrib

FASTJET_END_NAMESPACE
