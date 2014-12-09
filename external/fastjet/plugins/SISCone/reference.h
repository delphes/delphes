// -*- C++ -*-
///////////////////////////////////////////////////////////////////////////////
// File: reference.h                                                         //
// Description: header file for checkxor management (Creference class)       //
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
// $Revision:: 123                                                          $//
// $Date:: 2007-03-01 02:52:16 +0100 (Thu, 01 Mar 2007)                     $//
///////////////////////////////////////////////////////////////////////////////

#ifndef __REFERENCE_H__
#define __REFERENCE_H__

namespace siscone{

/**
 * \class Creference
 * \brief references used for checksums.
 *
 * This class implements some reference variable
 * that can be used for checksums. Those checksums
 * are useful to disentengle between contents of two
 * cones without looking into their explicit particle
 * contents.
 */
class Creference{
 public:
  /// default constructor
  Creference();

  /// create a random reference
  void randomize();

  /// test emptyness
  bool is_empty();

  /// test non-emptyness
  bool not_empty();

  /// assignment of reference
  Creference& operator = (const Creference &r);

  /// addition of reference
  Creference operator + (const Creference &r);

  /// incrementation of reference
  Creference& operator += (const Creference &r);

  /// decrementation of reference
  Creference& operator -= (const Creference &r);

  /// accessing the reference
  inline unsigned int operator[] (int i) {return ref[i];}

  unsigned int ref[3];   ///< actual data for the reference
};

/// addition of two references
Creference operator + (Creference &r1, Creference &r2);

/// equality test of two references
bool operator == (const Creference &r1, const Creference &r2);

/// difference test of two references
bool operator != (const Creference &r1, const Creference &r2);

/// ordering of two references
bool operator < (const Creference &r1, const Creference &r2);


//=============== inline material ================

// equality test for two references
//----------------------------------
inline bool operator == (const Creference &r1, const Creference &r2){
  return (r1.ref[0]==r2.ref[0]) && (r1.ref[1]==r2.ref[1]) && (r1.ref[2]==r2.ref[2]);
}

// difference test for two references
//----------------------------------
inline bool operator != (const Creference &r1, const Creference &r2){
  return (r1.ref[0]!=r2.ref[0]) || (r1.ref[1]!=r2.ref[1]) || (r1.ref[2]!=r2.ref[2]);
}

// difference test for two references
//----------------------------------
inline bool operator < (const Creference &r1, const Creference &r2){
  return (r1.ref[0]<r2.ref[0]) || ((r1.ref[0]==r2.ref[0]) && 
				   ((r1.ref[1]<r2.ref[1]) || ((r1.ref[1]==r2.ref[1]) && (r1.ref[2]<r2.ref[2]))
				    ));
}

}
#endif
