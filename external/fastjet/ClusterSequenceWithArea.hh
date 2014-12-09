//FJSTARTHEADER
// $Id: ClusterSequenceWithArea.hh 3433 2014-07-23 08:17:03Z salam $
//
// Copyright (c) 2006-2014, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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

#ifndef __FASTJET_CLUSTERSEQUENCEWITHAREA_HH__
#define __FASTJET_CLUSTERSEQUENCEWITHAREA_HH__

// for backwards compatibility fastjet/ClusterSequenceWithArea.hh
// provides the ClusterSequenceWithArea class that is equivalent to
// ClusterSequenceArea. The latter should now be used together with
// the fastjet/ClusterSequenceArea.hh. ClusterSequenceWithArea is not
// guaranteed to work in future release of FastJet
#warning You have included fastjet/ClusterSequenceWithArea.hh, \
a deprecated FastJet header provided only for backward compatibility. \
This is not guaranteed to work in future releases of FastJet. \
Please consider including fastjet/ClusterSequenceArea.hh directly. \
Similarily, if you use the (deprecated) class ClusterSequenceWithArea, \
please switch to the equivalent ClusterSequenceArea instead.

#include "fastjet/ClusterSequenceAreaBase.hh"

FASTJET_BEGIN_NAMESPACE

//----- backwards compatibility with version 2.1 ---------------
typedef ClusterSequenceAreaBase ClusterSequenceWithArea;

FASTJET_END_NAMESPACE

#endif // __FASTJET_CLUSTERSEQUENCEWITHAREA_HH__


