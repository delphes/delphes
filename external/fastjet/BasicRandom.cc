//FJSTARTHEADER
// $Id: BasicRandom.cc 3433 2014-07-23 08:17:03Z salam $
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

//   nlo includes
#include "fastjet/internal/BasicRandom.hh"


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


//
//                   random number generator
//  uses method of L'Ecuyer, (via F.James, comp. phys. comm. 60(1990)329)
//
int __default_random_generator(int *__iseed) 
{
  int __k = __iseed[0]/53668;
  __iseed[0] = (__iseed[0] - __k*53668)*40014 - __k*12211;
  if(__iseed[0] < 0) __iseed[0] += 2147483563;
  
  __k = __iseed[1]/52774;
  __iseed[1] = (__iseed[1] - __k*52774)*40692 - __k*3791;
  if(__iseed[1] < 0) __iseed[1] += 2147483399;
  
  int __iz = __iseed[0] - __iseed[1];
  if(__iz < 1) __iz += 2147483562;
  
  return __iz;
}

//   global defined random number generator
BasicRandom<int>     _G_random_int;
BasicRandom<double>  _G_random_double;


FASTJET_END_NAMESPACE

