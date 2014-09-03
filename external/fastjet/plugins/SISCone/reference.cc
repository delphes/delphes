///////////////////////////////////////////////////////////////////////////////
// File: reference.cpp                                                       //
// Description: source file for checkxor management (Creference class)       //
// This file is part of the SISCone project.                                 //
// For more details, see http://projects.hepforge.org/siscone                //
//                                                                           //
// Copyright (c) 2006 Gavin Salam and Gregory Soyez                          //
//                                                                           //
// This program is free software; you can redistribute it and/or modify      //
// it under the terms of the GNU General Public License as published by      //
// the Free Software Foundation; either version 2 of the License, or         //
// (at your option) any later version.                                       //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program; if not, write to the Free Software               //
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA //
//                                                                           //
// $Revision:: 311                                                          $//
// $Date:: 2011-10-05 23:27:09 +0200 (Wed, 05 Oct 2011)                     $//
///////////////////////////////////////////////////////////////////////////////

#include "reference.h"
#include "ranlux.h"
#include <stdlib.h>

namespace siscone{

/*******************************************************
 * Creference implementation                           *
 * references used for checksums.                      *
 *                                                     *
 * This class implements some reference variable       *
 * that can be used for checksums. Those checksums     *
 * are useful to disentengle between contents of two   *
 * cones without looking into their explicit particle  *
 * contents.                                           *
 *******************************************************/

// default constructor
//////////////////////
Creference::Creference(){
  ref[0] = ref[1] = ref[2] = 0;
}

  //static unsigned int reference_bit = 1;

// create a random reference
//---------------------------
void Creference::randomize(){
//  ref[0] = reference_bit;
//  ref[1] = 0;
//  ref[2] = 0;
//  reference_bit <<= 1;

  unsigned int r1 = ranlux_get();
  unsigned int r2 = ranlux_get();
  unsigned int r3 = ranlux_get();
  unsigned int r4 = ranlux_get();
  // since ranlux only produces 24 bits, take r4 and add 8 bits
  // from it to each of r1,r2, r3 to get 3*32 bits.
  ref[0] = r1+((r4 & 0x00ff0000) <<  8);
  ref[1] = r2+((r4 & 0x0000ff00) << 16);
  ref[2] = r3+((r4 & 0x000000ff) << 24);

  if (is_empty()) randomize();
}

// test emptyness
//----------------
bool Creference::is_empty(){
  return (ref[0]==0) && (ref[1]==0) && (ref[2]==0);
}

// test non-emptyness
//--------------------
bool Creference::not_empty(){
  return (ref[0]!=0) || (ref[1]!=0) || (ref[2]!=0);
}

// assignment of reference
//-------------------------
Creference& Creference::operator = (const Creference &r){
  ref[0] = r.ref[0];  
  ref[1] = r.ref[1];
  ref[2] = r.ref[2];
  return *this;
}

// addition of reference
//-----------------------
Creference Creference::operator + (const Creference &r){
  Creference tmp = *this;
  return tmp+=r;
}

// incrementation of reference
//-----------------------------
Creference& Creference::operator += (const Creference &r){
  ref[0] ^= r.ref[0];  
  ref[1] ^= r.ref[1];
  ref[2] ^= r.ref[2];
  return *this; 
}

// decrementation of reference
//-----------------------------
Creference& Creference::operator -= (const Creference &r){
  ref[0] ^= r.ref[0];  
  ref[1] ^= r.ref[1];
  ref[2] ^= r.ref[2];
  return *this; 
}

}

