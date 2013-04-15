//STARTHEADER
// $Id: ClosestPair2DBase.hh 2577 2011-09-13 15:11:38Z salam $
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

#ifndef __FASTJET_CLOSESTPAIR2DBASE__HH__
#define __FASTJET_CLOSESTPAIR2DBASE__HH__

#include<vector>
#include "fastjet/internal/base.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//----------------------------------------------------------------------
/// \if internal_doc
/// @ingroup internal
/// \class Coord2D
/// class for representing 2d coordinates and carrying out some basic 
/// operations on them
/// \endif
class Coord2D {
public:
  double x, y;

  Coord2D() {};

  Coord2D(double a, double b): x(a), y(b) {};

  /// return the vector difference between two coordinates
  Coord2D operator-(const Coord2D & other) const {
    return Coord2D(x - other.x,  y - other.y);};

  /// return the vector sum between two coordinates
  Coord2D operator+(const Coord2D & other) const {
    return Coord2D(x + other.x,  y + other.y);};

  /// return the product of the coordinate with the factor
  Coord2D operator*(double factor) const {return Coord2D(factor*x,factor*y);};
  friend Coord2D operator*(double factor, const Coord2D & coord) {
    return Coord2D(factor*coord.x,factor*coord.y);
  }

  /// division of each component of coordinate
  Coord2D operator/(double divisor) const {
    return Coord2D(x / divisor,  y / divisor);};

  /// return the squared distance between two coordinates
  friend double distance2(const Coord2D & a, const Coord2D & b) {
    double dx = a.x - b.x, dy = a.y-b.y;
    return dx*dx+dy*dy;
  };
  /// return the squared distance between two coordinates
  double distance2(const Coord2D & b) const {
    double dx = x - b.x, dy = y-b.y;
    return dx*dx+dy*dy;
  };
};


//----------------------------------------------------------------------
/// \if internal_doc
/// @ingroup internal
/// \class ClosestPair2DBase
/// abstract base class for finding closest pairs in 2D
/// \endif
class ClosestPair2DBase {
public:
  /// provides the IDs of the closest pair as well as the squared
  /// distance between them
  virtual void closest_pair(unsigned int & ID1, unsigned int & ID2, 
			    double & distance2) const = 0;
  
  /// removes the entry labelled by ID from the object;
  virtual void remove(unsigned int ID) = 0;
 
  /// inserts the position into the closest pair structure and returns the
  /// ID that has been allocated for the object.
  virtual unsigned int insert(const Coord2D & position) = 0;

  /// replaces the specified ID1 and ID2 with something at a new position
  /// assuming that ID1 and ID2 are in sequence wrt position; it returns
  /// the ID of the new object...
  virtual unsigned int replace(unsigned int ID1, unsigned int ID2, 
			       const Coord2D & position) {
    remove(ID1); 
    remove(ID2); 
    unsigned new_ID = insert(position);
    return(new_ID);
  };

  /// replaces IDs_to_remove with points at the new_positions
  /// indicating the IDs allocated to the new points in new_IDs
  virtual void replace_many(const std::vector<unsigned int> & IDs_to_remove,
		       const std::vector<Coord2D> & new_positions,
		       std::vector<unsigned int> & new_IDs) {
    for(unsigned i = 0; i < IDs_to_remove.size(); i++) {
      remove(IDs_to_remove[i]);}
    new_IDs.resize(0);
    for(unsigned i = 0; i < new_positions.size(); i++) {
      new_IDs.push_back(insert(new_positions[i]));}
  }

  virtual unsigned int size() = 0;

  virtual ~ClosestPair2DBase() {};
  
};


FASTJET_END_NAMESPACE

#endif // __FASTJET_CLOSESTPAIR2DBASE__HH__
