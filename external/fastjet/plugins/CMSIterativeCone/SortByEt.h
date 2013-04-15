#ifndef __CMS_ITERATIVE_CONE__SORT_BY_ET_H__
#define __CMS_ITERATIVE_CONE__SORT_BY_ET_H__

//STARTHEADER
// $Id$
//
// Copyright (c) ????-????, CMS collaboration
// Copyright (c) 2009-2011, Matteo Cacciari, Gavin P. Salam and Gregory Soyez [for the changes listed below]
//
//----------------------------------------------------------------------
// This file distributed with FastJet has been obtained from the CMS
// collaboration, revision 1.2 of the EtComparator.h file in CMSSW,
// see
//   http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/PhysicsTools/Utilities/interface/EtComparator.h?hideattic=0&revision=1.2&view=markup
//
// Permission has been granted by the CMS collaboration to release it
// in FastJet under the terms of the GNU Public License(v2) (see the
// COPYING file in the main FastJet directory for details).
// Changes from the original file are listed below.
//
// FastJet is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// The algorithms that underlie FastJet have required considerable
// development and are described in hep-ph/0512210. If you use
// FastJet as part of work towards a scientific publication, please
// include a citation to the FastJet paper.
//
// FastJet is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//ENDHEADER

// List of changes compared to the original CMS code (revision 1.2 of
// EtComparator.h)
//
// 2009-01-06  Gregory Soyez  <soyez@fastjet.fr>
//
//        * Extracted (only) NumericSafeGreaterByEt from the CMS code
//          and adapted it to act on PseudoJet rather than CMS types
//          for 4-momenta
//        * Put the code in a fastjet::cms namespace

#include <limits>
#include <fastjet/internal/base.hh>

FASTJET_BEGIN_NAMESPACE

namespace cms{

template <class T>
struct NumericSafeGreaterByEt {
  typedef T first_argument_type;
  typedef T second_argument_type;
  bool operator()(const T& a1, const T& a2) {
    // FastJet::PseudoJet does not provide a direct access to Et2
    // Plus, we want it to be computed in the same way as in the CMS
    // code (actually the Root code that is used by CMS)
    double et1 = a1.Et();
    double et2 = a2.Et();

    // now we can come back to the CMS code
    return
      fabs (et1-et2) > std::numeric_limits<double>::epsilon() ? et1 > et2 :
      fabs (a1.px()-a2.px()) > std::numeric_limits<double>::epsilon() ? a1.px() > a2.px() :
      a1.pz() > a2.pz();
  }
};

}  // namespace cms

FASTJET_END_NAMESPACE


#endif   // __CMS_ITERATIVE_CONE__SORT_BY_ET_H__
