
//FJSTARTHEADER
// $Id: deprecated.hh 4098 2016-03-15 16:38:22Z salam $
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

#ifndef __FASTJET_FASTJET_DEPRECATED_HH__
#define __FASTJET_FASTJET_DEPRECATED_HH__

#include "fastjet/config.h"

// define a deprecation macro based on the capabilities of the compiler
// (as determined at configure time).
#if defined(FASTJET_HAVE_CXX14_DEPRECATED)
#define FASTJET_DEPRECATED               [[deprecated]]
#define FASTJET_DEPRECATED_MSG(message)  [[deprecated(message)]]
#elif defined(FASTJET_HAVE_GNUCXX_DEPRECATED)
#define FASTJET_DEPRECATED               __attribute__((__deprecated__))
#define FASTJET_DEPRECATED_MSG(message)  __attribute__((__deprecated__))
#else
#define FASTJET_DEPRECATED               
#define FASTJET_DEPRECATED_MSG(message) 
#endif

#endif // __FASTJET_FASTJET_DEPRECATED_HH__
