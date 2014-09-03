//FJSTARTHEADER
// $Id: ClosestPair2D.hh 3433 2014-07-23 08:17:03Z salam $
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

#ifndef __FASTJET_CLOSESTPAIR2D__HH__
#define __FASTJET_CLOSESTPAIR2D__HH__

#include<vector>
#include<stack>
#include<iostream>
#include "fastjet/internal/ClosestPair2DBase.hh"
#include "fastjet/internal/SearchTree.hh"
#include "fastjet/internal/MinHeap.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//----------------------------------------------------------------------
/// \if internal_doc
/// @ingroup internal
/// \class ClosestPair2D
/// concrete implementation for finding closest pairs in 2D -- will
/// use Chan's (hopefully efficient) shuffle based structures
/// \endif
class ClosestPair2D : public ClosestPair2DBase {
public:
  /// constructor from a vector of 2D positions -- number of objects
  /// after insertion and deletion must never exceed positions.size();
  /// objects are given IDs that correspond to their index in the vector 
  /// of positions
  ClosestPair2D(const std::vector<Coord2D> & positions, 
		const Coord2D & left_corner, const Coord2D & right_corner) {
    _initialize(positions, left_corner, right_corner, positions.size());
  };

  /// constructor which allows structure to grow beyond positions.size(), up
  /// to max_size
  ClosestPair2D(const std::vector<Coord2D> & positions, 
		const Coord2D & left_corner, const Coord2D & right_corner,
		const unsigned int max_size) {
    _initialize(positions, left_corner, right_corner, max_size);
  };

  /// provides the IDs of the closest pair as well as the distance between
  /// them
  void closest_pair(unsigned int & ID1, unsigned int & ID2, 
		    double & distance2) const;

  /// removes the entry labelled by ID from the object;
  void remove(unsigned int ID);
 
  /// inserts the position into the closest pair structure and returns the
  /// ID that has been allocated for the object.
  unsigned int insert(const Coord2D &);

  /// removes ID1 and ID2 and inserts position, returning the ID 
  /// corresponding to position...
  virtual unsigned int replace(unsigned int ID1, unsigned int ID2, 
			       const Coord2D & position);

  /// replaces IDs_to_remove with points at the new_positions
  /// indicating the IDs allocated to the new points in new_IDs
  virtual void replace_many(const std::vector<unsigned int> & IDs_to_remove,
			    const std::vector<Coord2D> & new_positions,
			    std::vector<unsigned int> & new_IDs);

  // mostly for checking how things are working...
  inline void print_tree_depths(std::ostream & outdev) const {
    outdev    << _trees[0]->max_depth() << " "
	      << _trees[1]->max_depth() << " "
	      << _trees[2]->max_depth() << "\n";
  };

  unsigned int size();

private:
  
  void _initialize(const std::vector<Coord2D> & positions, 
	      const Coord2D & left_corner, const Coord2D & right_corner,
	      const unsigned int max_size);

  static const unsigned int _nshift = 3;

  class Point; // will be defined below

  /// since sets of three objects will crop up repeatedly, useful
  /// to have a triplet class?
  template<class T> class triplet {
  public:
    inline const T & operator[](unsigned int i) const {return _contents[i];};
    inline       T & operator[](unsigned int i)       {return _contents[i];};
  private:
    T _contents[_nshift];
  };


  /// class that will take care of ordering of shuffles for us
  class Shuffle {
  public:
    unsigned int x, y;
    Point * point;
    bool operator<(const Shuffle &) const;
    void operator+=(unsigned int shift) {x += shift; y+= shift;};
  };

  typedef SearchTree<Shuffle>     Tree;
  typedef Tree::circulator        circulator;
  typedef Tree::const_circulator  const_circulator;


  triplet<std::auto_ptr<Tree> >  _trees;
  std::auto_ptr<MinHeap> _heap;
  std::vector<Point>     _points;
  std::stack<Point *>    _available_points;

  /// points that are "under review" in some way
  std::vector<Point *>   _points_under_review;

  // different statuses for review
  static const unsigned int _remove_heap_entry = 1;
  static const unsigned int _review_heap_entry = 2;
  static const unsigned int _review_neighbour  = 4;

  /// add a label to a point as to the nature of review needed
  /// (includes adding it to list of points needing review) [doesn't
  /// affect other labels already set for the point]
  void _add_label(Point * point, unsigned int review_flag);

  /// sets the label for the point to be exclusively this 
  /// review flag (and adds it to list of points needing review
  /// if not already there)
  void _set_label(Point * point, unsigned int review_flag);

  /// for all entries of the _points_under_review[] vector, carry out
  /// the actions indicated by its review flag; the points are
  /// then removed from _points_under_review[] and their flags
  /// set to zero
  void _deal_with_points_to_review();

  /// carry out the search-tree related operations of point removal
  void _remove_from_search_tree(Point * point_to_remove);

  /// carry out the search-tree related operations of point insertion
  void _insert_into_search_tree(Point * new_point);

  /// takes a point and creates a shuffle with the given shift
  void _point2shuffle(Point & , Shuffle & , unsigned int shift);

  /// pieces needed for converting coordinates to integer
  Coord2D _left_corner;
  double _range;

  int _ID(const Point *) const;

  triplet<unsigned int> _shifts;     // absolute shifts
  triplet<unsigned int> _rel_shifts; // shifts relative to previous shift

  unsigned int _cp_search_range;
};


//----------------------------------------------------------------------
/// \if internal_doc
/// @ingroup internal
/// \class ClosestPair2D::Point
/// class for representing all info needed about a point
/// \endif
class ClosestPair2D::Point {
public:
  /// the point's coordinates
  Coord2D coord;
  /// a pointer to its closest neighbour in our structure
  Point * neighbour;
  /// the corresponding squared distance
  double  neighbour_dist2;
  /// circulators for each of the shifts of the shuffles
  triplet<circulator> circ;

  /// indicates that something special is currently happening to this point
  unsigned int review_flag;

  /// returns the distance between two of these objects
  double distance2(const Point & other) const {
    return coord.distance2(other.coord);
  };

  /// creates a shuffle for us with a given shift
  //void set_shuffle(Shuffle & shuffle);
};


//----------------------------------------------------------------------
/// returns true if floor(ln_base2(x)) < floor(ln_base2(y)), using
/// Chan's neat trick...
inline bool floor_ln2_less(unsigned x, unsigned y) {
  if (x>y) return false;
  return (x < (x^y)); // beware of operator precedence...
}


//----------------------------------------------------------------------
/// returns the ID for the specified point...
inline int ClosestPair2D::_ID(const Point * point) const {
  return point - &(_points[0]);
}


//
inline unsigned int ClosestPair2D::size() {
  return _points.size() - _available_points.size();
}



FASTJET_END_NAMESPACE

#endif // __FASTJET_CLOSESTPAIR2D__HH__
