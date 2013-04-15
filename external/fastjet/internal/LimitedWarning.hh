#ifndef __FASTJET_INTERNALLIMITEDWARNING_HH__
#define __FASTJET_INTERNALLIMITEDWARNING_HH__

//STARTHEADER
// $Id: LimitedWarning.hh 2577 2011-09-13 15:11:38Z salam $
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



// we have moved LimitedWarning.hh one directory up; still allow
// old form of access
#include "fastjet/LimitedWarning.hh"

#warning *** You have included fastjet/internal/LimitedWarning.hh. \
Access to LimitedWarning through this header is deprecated as of FJ3.0. \
Please instead use fastjet/LimitedWarning.hh

#endif // __FASTJET_INTERNALLIMITEDWARNING_HH__
