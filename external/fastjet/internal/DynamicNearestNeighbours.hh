//STARTHEADER
// $Id: DynamicNearestNeighbours.hh 2687 2011-11-14 11:17:51Z soyez $
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


#ifndef __FASTJET_DYNAMICNEARESTNEIGHBOURS_HH__
#define __FASTJET_DYNAMICNEARESTNEIGHBOURS_HH__

#include<vector>
#include<string>
#include<iostream>
#include<sstream>
#include<cassert>
#include "fastjet/internal/numconsts.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

/// Shortcut for dealing with eta-phi coordinates.
//typedef std::pair<double,double> EtaPhi;

/// \if internal_doc
/// @ingroup internal
/// \class EtaPhi
/// use a class instead of a pair so that phi can be sanitized
/// and put into proper range on initialization.
/// \endif
class EtaPhi {
public:
  double first, second;
  EtaPhi() {}
  EtaPhi(double a, double b) {first = a; second = b;}
  /// put things into the desired range.
  void sanitize() {    
    if (second <  0)     second += twopi; 
    if (second >= twopi) second -= twopi;
  }

};

/// \if internal_doc
/// @ingroup internal
/// \class DnnError
/// class corresponding to errors that will be thrown by Dynamic
/// Nearest Neighbours code
/// \endif
class DnnError {
public:
  // constructors
  DnnError() {;};
  DnnError(const std::string & message_in) {
    _message = message_in; std::cerr << message_in << std::endl;};

  std::string message() const {return _message;};

private:
  std::string _message;
};


/// \if internal_doc
/// @ingroup internal
/// \class DynamicNearestNeighbours
/// Abstract base class for quick location of nearest neighbours in a set of
/// points.
///
/// Abstract base class for quick location of nearest neighbours in a set of
/// points, with facilities for adding and removing points from the
/// set after initialisation. Derived classes will be
/// named according to the convention DnnSomeName (e.g. DnnPlane).
///
/// The main purpose of this abstract base class is to define the
/// general interface of a whole set of classes that deal with
/// nearest-neighbour location on different 2-d geometries and with
/// various underlying data structures and algorithms.
///
/// \endif
class DynamicNearestNeighbours {

public:
  /// Dummy initialiser --- does nothing!
  //virtual DynamicNearestNeighbours() {};
   
  /// Initialiser --- sets up the necessary structures to allow efficient
  /// nearest-neighbour finding on the std::vector<EtaPhi> of input points
  //virtual DynamicNearestNeighbours(const std::vector<EtaPhi> &, 
  //				   const bool & verbose = false ) = 0;

  /// Returns the index of the nearest neighbour of point labelled
  /// by ii (assumes ii is valid)
  virtual int NearestNeighbourIndex(const int & ii) const = 0;

  /// Returns the distance to the nearest neighbour of point labelled
  /// by index ii (assumes ii is valid)
  virtual double NearestNeighbourDistance(const int & ii) const = 0;

  /// Returns true iff the given index corresponds to a point that
  /// exists in the DNN structure (meaning that it has been added, and
  /// not removed in the meantime)
  virtual bool Valid(const int & index) const = 0;

  /// remove the points labelled by the std::vector indices_to_remove, and
  /// add the points specified by the std::vector points_to_add
  /// (corresponding indices will be calculated automatically); the
  /// idea behind this routine is that the points to be added will
  /// somehow be close to the one or other of the points being removed
  /// and this can be used by the implementation to provide hints for
  /// inserting the new points in whatever structure it is using.  In a
  /// kt-algorithm the points being added will be a result of a
  /// combination of the points to be removed -- hence the proximity
  /// is (more or less) guaranteed.
  virtual void RemoveAndAddPoints(const std::vector<int> & indices_to_remove,
			  const std::vector<EtaPhi> & points_to_add,
			  std::vector<int> & indices_added,
			  std::vector<int> & indices_of_updated_neighbours) = 0;


  /// Remove the point labelled by index and return the list of
  /// points whose nearest neighbours have changed in the process
  inline void RemovePoint (const int & index,
			   std::vector<int> & indices_of_updated_neighbours) {
    std::vector<int> indices_added;
    std::vector<EtaPhi> points_to_add;
    std::vector<int> indices_to_remove(1);
    indices_to_remove[0] = index;
    RemoveAndAddPoints(indices_to_remove, points_to_add, indices_added,
		       indices_of_updated_neighbours
		       );};


  /// Removes the two points labelled by index1, index2 and adds in the
  /// a point with coordinates newpoint; it returns an index for the new 
  /// point (index 3) and a std::vector of indices of neighbours whose
  /// nearest neighbour has changed (the list includes index3, i.e. the new
  /// point).
  inline void RemoveCombinedAddCombination(
			const int & index1, const int & index2,
			const EtaPhi & newpoint,
			int & index3,
			std::vector<int> & indices_of_updated_neighbours) {
    std::vector<int> indices_added(1);
    std::vector<EtaPhi> points_to_add(1);
    std::vector<int> indices_to_remove(2);
    indices_to_remove[0] = index1;
    indices_to_remove[1] = index2;
    points_to_add[0] = newpoint;
    RemoveAndAddPoints(indices_to_remove, points_to_add, indices_added,
		       indices_of_updated_neighbours
		       );
    index3 = indices_added[0];
  };

  /// destructor -- here it is now implemented
  virtual ~DynamicNearestNeighbours () {}
};
  

FASTJET_END_NAMESPACE

#endif // __FASTJET_DYNAMICNEARESTNEIGHBOURS_HH__
