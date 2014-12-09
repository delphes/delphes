//FJSTARTHEADER
// $Id: RangeDefinition.cc 3433 2014-07-23 08:17:03Z salam $
//
// Copyright (c) 2005-2014, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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
//  development. They are described in the original FastJet paper,
//  hep-ph/0512210 and in the manual, arXiv:1111.6097. If you use
//  FastJet as part of work towards a scientific publication, please
//  quote the version you use and include a citation to the manual and
//  optionally also to hep-ph/0512210.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//FJENDHEADER

#include "fastjet/RangeDefinition.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

using namespace std;
 
LimitedWarning RangeDefinition::_warnings_deprecated;

// calculate, and set in _total_area, the area with a numerical test
// takes a reasonable time with rapmax = 10, npoints = 100
void RangeDefinition::_numerical_total_area(double rapmax, int npoints) {

      int count = 0;
      double deltaphi = twopi/double(npoints);
      double deltarap = 2.0*rapmax/double(npoints);
      double phi = 0.0;
      for(int i = 0; i < npoints; i++) {
        double rap = -rapmax;
        for (int j = 0; j < npoints; j++) {
	  if ( is_in_range(rap,phi) ) { count++; }
          rap += deltarap;
        }
        phi += deltaphi;
      }

      _total_area = double(count)/double(npoints*npoints)*2.0*twopi*rapmax;
}


// the use of RangeDefinition is deprecated since FastJet version
// 3.0 onwards. Please use Selector instead.  
// RangeDefinition is only provided for backward compatibility
// reasons and is not guaranteed to work in future releases of
// FastJet.
void RangeDefinition::_warn_deprecated() const{
  _warnings_deprecated.warn("The use of RangeDefinition is deprecated since FastJet version 3.0 onwards. Please consider using Selector (defined in fastjet/Selector.hh) instead. There is no guarantee that support for RangeDefinition will be provided in future releases of FastJet.");
}

FASTJET_END_NAMESPACE
