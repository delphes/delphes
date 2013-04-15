//STARTHEADER
// simple random number generator class taken from nlojet++.
// $Id$
//
//  Copyright (C) 2002 Zoltan Nagy
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//ENDHEADER

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

