#ifndef __FASTJET_INTERNALLIMITEDWARNING_HH__
#define __FASTJET_INTERNALLIMITEDWARNING_HH__

//FJSTARTHEADER
// $Id: LimitedWarning.hh 3433 2014-07-23 08:17:03Z salam $
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



// we have moved LimitedWarning.hh one directory up; still allow
// old form of access
#include "fastjet/LimitedWarning.hh"

#warning *** You have included fastjet/internal/LimitedWarning.hh. \
Access to LimitedWarning through this header is deprecated as of FJ3.0. \
Please instead use fastjet/LimitedWarning.hh

#endif // __FASTJET_INTERNALLIMITEDWARNING_HH__
