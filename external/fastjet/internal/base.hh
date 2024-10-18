
//FJSTARTHEADER
// $Id$
//
// Copyright (c) 2005-2024, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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

#ifndef __FASTJET_FASTJET_BASE_HH__
#define __FASTJET_FASTJET_BASE_HH__

#include "fastjet/config.h"

/// \namespace fastjet
/// the FastJet namespace
/// 
/// all the fastjet-related material is put under that namespace

// define this for easier readability (and obfuscation?) in
// a range of places
#define FASTJET_BEGIN_NAMESPACE namespace fastjet {
#define FASTJET_END_NAMESPACE   }

// define a macro to mark virtual function in derived classes as
// overriding the base-class definition
#ifdef FASTJET_HAVE_OVERRIDE
#define FASTJET_OVERRIDE  override
#else
#define FASTJET_OVERRIDE  
#endif

#endif // __FASTJET_FASTJET_BASE_HH__
