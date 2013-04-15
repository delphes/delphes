//STARTHEADER
// $Id$
//
// Copyright (c) 2005-2011, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of FastJet.
//
//  FastJet is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  The algorithms that underlie FastJet have required considerable
//  development and are described in hep-ph/0512210. If you use
//  FastJet as part of work towards a scientific publication, please
//  include a citation to the FastJet paper.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//ENDHEADER

#include "fastjet/tools/Subtractor.hh"
#include <cassert>
#include <sstream>
#include <limits>
using namespace std;

FASTJET_BEGIN_NAMESPACE     // defined in fastjet/internal/base.hh

const double Subtractor::_invalid_rho = -numeric_limits<double>::infinity();


Subtractor::Subtractor(double rho) : _bge(0), _rho(rho) {
  assert(_rho>0.0);
}

PseudoJet Subtractor::result(const PseudoJet & jet) const {
  if (!jet.has_area()){
    throw Error("Trying to subtract a jet without area support");
  }
  
  double rho;
  if (_bge != 0) {
    rho = _bge->rho(jet);
  } else if (_rho != _invalid_rho) {
    rho = _rho;
  } else {
    throw Error("default Subtractor does not have any information about the background, which is needed to perform the subtraction");
  }

  PseudoJet subtracted_jet = jet;
  PseudoJet area4vect = jet.area_4vector();
  // sanity check
  if (rho*area4vect.perp() < jet.perp() ) { 
    // this subtraction should retain the jet's structural
    // information
    subtracted_jet -= rho*area4vect;
  } else { 
    // this sets the jet's momentum to zero while
    // maintaining all of the jet's structural information
    subtracted_jet *= 0;
  }
  return subtracted_jet;
}

//----------------------------------------------------------------------
std::string Subtractor::description() const{
  if (_bge != 0) {
    return "Subtractor that uses the following background estimator to determine rho: "+_bge->description();
  } else if (_rho != _invalid_rho) {
    ostringstream ostr;
    ostr << "Subtractor that uses a fixed value of rho = " << _rho;
    return ostr.str();
  } else {
    return "Uninitialised subtractor";
  }
}

FASTJET_END_NAMESPACE
