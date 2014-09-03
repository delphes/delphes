//FJSTARTHEADER
// $Id: ClosestPair2D.cc 3433 2014-07-23 08:17:03Z salam $
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

#include "fastjet/internal/ClosestPair2D.hh"

#include<limits>
#include<iostream>
#include<iomanip>
#include<algorithm>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

const unsigned int twopow31      = 2147483648U;

using namespace std;

//----------------------------------------------------------------------
/// takes a point and sets a shuffle with the given shift...
void ClosestPair2D::_point2shuffle(Point & point, Shuffle & shuffle, 
				  unsigned int shift) {
  
  Coord2D renorm_point = (point.coord - _left_corner)/_range;
  // make sure the point is sensible
  //cerr << point.coord.x <<" "<<point.coord.y<<endl;
  assert(renorm_point.x >=0);
  assert(renorm_point.x <=1);
  assert(renorm_point.y >=0);
  assert(renorm_point.y <=1);
  
  shuffle.x = static_cast<unsigned int>(twopow31 * renorm_point.x) + shift;
  shuffle.y = static_cast<unsigned int>(twopow31 * renorm_point.y) + shift;
  shuffle.point = &point;
}

//----------------------------------------------------------------------
/// compares this shuffle with the other one
bool ClosestPair2D::Shuffle::operator<(const Shuffle & q) const {

  if (floor_ln2_less(x ^ q.x, y ^ q.y)) {
    // i = 2 in Chan's algorithm
    return (y < q.y);
  } else {
    // i = 1 in Chan's algorithm
    return (x < q.x);
  }
}



//----------------------------------------------------------------------
void ClosestPair2D::_initialize(const std::vector<Coord2D> & positions, 
			     const Coord2D & left_corner, 
			     const Coord2D & right_corner,
			     unsigned int max_size) {

  unsigned int n_positions = positions.size();
  assert(max_size >= n_positions);

  //_points(positions.size())

  // allow the points array to grow to the following size
  _points.resize(max_size);
  // currently unused points are immediately made available on the
  // stack
  for (unsigned int i = n_positions; i < max_size; i++) {
    _available_points.push(&(_points[i]));
  }

  _left_corner = left_corner;
  _range       = max((right_corner.x - left_corner.x),
		     (right_corner.y - left_corner.y));

  // initialise the coordinates for the points and create the zero-shifted
  // shuffle array
  vector<Shuffle> shuffles(n_positions);
  for (unsigned int i = 0; i < n_positions; i++) {
    // set up the points
    _points[i].coord = positions[i];
    _points[i].neighbour_dist2 = numeric_limits<double>::max();
    _points[i].review_flag = 0;

    // create shuffle with 0 shift.
    _point2shuffle(_points[i], shuffles[i], 0);
  }

  // establish what our shifts will be
  for (unsigned ishift = 0; ishift < _nshift; ishift++) {
    // make sure we use double-precision for calculating the shifts
    // since otherwise we will hit integer overflow.
   _shifts[ishift] = static_cast<unsigned int>(((twopow31*1.0)*ishift)/_nshift);
    if (ishift == 0) {_rel_shifts[ishift] = 0;}
    else {_rel_shifts[ishift] = _shifts[ishift] - _shifts[ishift-1];}
  }
  //_shifts[0] = 0;
  //_shifts[1] = static_cast<unsigned int>((twopow31*1.0)/3.0);
  //_shifts[2] = static_cast<unsigned int>((twopow31*2.0)/3.0);
  //_rel_shifts[0] = 0;
  //_rel_shifts[1] = _shifts[1];
  //_rel_shifts[2] = _shifts[2]-_shifts[1];

  // and how we will search...
  //_cp_search_range = 49;
  _cp_search_range = 30;
  _points_under_review.reserve(_nshift * _cp_search_range);

  // now initialise the three trees
  for (unsigned int ishift = 0; ishift < _nshift; ishift++) {

    // shift the shuffles if need be.
    if (ishift > 0) {
      unsigned rel_shift = _rel_shifts[ishift];
      for (unsigned int i = 0; i < shuffles.size(); i++) {
	shuffles[i] += rel_shift; }
    }

    // sort the shuffles
    sort(shuffles.begin(), shuffles.end());

    // and create the search tree
    _trees[ishift] = auto_ptr<Tree>(new Tree(shuffles, max_size));

    // now we look for the closest-pair candidates on this tree
    circulator circ = _trees[ishift]->somewhere(), start=circ;
    // the actual range in which we search 
    unsigned int CP_range = min(_cp_search_range, n_positions-1);
    do {
      Point * this_point = circ->point;
      //cout << _ID(this_point) << " ";
      this_point->circ[ishift] = circ;
      // run over all points within _cp_search_range of this_point on tree
      circulator other = circ;
      for (unsigned i=0; i < CP_range; i++) {
	++other;
	double dist2 = this_point->distance2(*other->point);
	if (dist2 < this_point->neighbour_dist2) {
	  this_point->neighbour_dist2 = dist2;
	  this_point->neighbour       = other->point;
	}
      }
    } while (++circ != start);
    //cout << endl<<endl;
  }

  // now initialise the heap object...
  vector<double> mindists2(n_positions);
  for (unsigned int i = 0; i < n_positions; i++) {
    mindists2[i] = _points[i].neighbour_dist2;}
  
  _heap = auto_ptr<MinHeap>(new MinHeap(mindists2, max_size));
}


//----------------------------------------------------------------------=
void ClosestPair2D::closest_pair(unsigned int & ID1, unsigned int & ID2, 
				 double & distance2) const {
  ID1 = _heap->minloc();
  ID2 = _ID(_points[ID1].neighbour);
  distance2 = _points[ID1].neighbour_dist2;
  // we make the swap explicitly in the std namespace to avoid
  // potential conflict with the fastjet::swap introduced by
  // SharedPtr.
  // This is only an issue because we are in the fastjet namespace
  if (ID1 > ID2) std::swap(ID1,ID2);
}


//----------------------------------------------------------------------
inline void ClosestPair2D::_add_label(Point * point, unsigned int review_flag) {

  // if it's not already under review, then put it on the list of
  // points needing review
  if (point->review_flag == 0) _points_under_review.push_back(point);

  // OR the point's current flag with the requested review flag
  point->review_flag |= review_flag;
}

//----------------------------------------------------------------------
inline void ClosestPair2D::_set_label(Point * point, unsigned int review_flag) {

  // if it's not already under review, then put it on the list of
  // points needing review
  if (point->review_flag == 0) _points_under_review.push_back(point);

  // SET the point's current flag to the requested review flag
  point->review_flag = review_flag;
}

//----------------------------------------------------------------------
void ClosestPair2D::remove(unsigned int ID) {

  //cout << "While removing " << ID <<endl;

  Point * point_to_remove = & (_points[ID]);

  // remove this point from the search tree
  _remove_from_search_tree(point_to_remove);

  // the above statement labels certain points as needing "review" --
  // deal with them...
  _deal_with_points_to_review();
}


//----------------------------------------------------------------------
void ClosestPair2D::_remove_from_search_tree(Point * point_to_remove) {

  // add this point to the list of "available" points (this is
  // relevant also from the point of view of determining the range
  // over which we circulate).
  _available_points.push(point_to_remove);

  // label it so that it goes from the heap
  _set_label(point_to_remove, _remove_heap_entry);

  // establish the range over which we search (a) for points that have
  // acquired a new neighbour and (b) for points which had ID as their
  // neighbour;
  
  unsigned int CP_range = min(_cp_search_range, size()-1);

  // then, for each shift
  for (unsigned int ishift = 0; ishift < _nshift; ishift++) {
    //cout << "   ishift = " << ishift <<endl;
    // get the circulator for the point we'll remove and its successor
    circulator removed_circ = point_to_remove->circ[ishift];
    circulator right_end = removed_circ.next();
    // remove the point
    _trees[ishift]->remove(removed_circ);
    
    // next find the point CP_range points to the left
    circulator left_end  = right_end, orig_right_end = right_end;
    for (unsigned int i = 0; i < CP_range; i++) {left_end--;}

    if (size()-1 < _cp_search_range) {
      // we have a smaller range now than before -- but when seeing who 
      // could have had ID as a neighbour, we still need to use the old
      // range for seeing how far back we search (but new separation between
      // points). [cf CCN28-42]
      left_end--; right_end--;
    }

    // and then for each left-end point: establish if the removed
    // point was its neighbour [in which case a new neighbour must be
    // found], otherwise see if the right-end point is a closer neighbour
    do {
      Point * left_point = left_end->point;

      //cout << "    visited " << setw(3)<<_ID(left_point)<<" (its neighbour was "<<	setw(3)<< _ID(left_point->neighbour)<<")"<<endl;

      if (left_point->neighbour == point_to_remove) {
	// we'll deal with it later...
	_add_label(left_point, _review_neighbour);
      } else {
	// check to see if right point has become its closest neighbour
	double dist2 = left_point->distance2(*right_end->point);
	if (dist2 < left_point->neighbour_dist2) {
	  left_point->neighbour = right_end->point;
	  left_point->neighbour_dist2 = dist2;
	  // NB: (LESSER) REVIEW NEEDED HERE TOO...
	  _add_label(left_point, _review_heap_entry);
	}
      }
      ++right_end;
    } while (++left_end != orig_right_end);
  } // ishift...

}


//----------------------------------------------------------------------
void ClosestPair2D::_deal_with_points_to_review() {

  // the range in which we carry out searches for new neighbours on
  // the search tree
  unsigned int CP_range = min(_cp_search_range, size()-1);

  // now deal with the points that are "under review" in some way
  // (have lost their neighbour, or need their heap entry updating)
  while(_points_under_review.size() > 0) {
    // get the point to be considered
    Point * this_point = _points_under_review.back();
    // remove it from the list
    _points_under_review.pop_back();  
    
    if (this_point->review_flag & _remove_heap_entry) {
      // make sure no other flags are on (it wouldn't be consistent?)
      assert(!(this_point->review_flag ^ _remove_heap_entry));
      _heap->remove(_ID(this_point));
    } 
    // check to see if the _review_neighbour flag is on
    else {
      if (this_point->review_flag & _review_neighbour) {
	this_point->neighbour_dist2 = numeric_limits<double>::max();
	// among all three shifts
	for (unsigned int ishift = 0; ishift < _nshift; ishift++) {
	  circulator other = this_point->circ[ishift];
	  // among points within CP_range
	  for (unsigned i=0; i < CP_range; i++) {
	    ++other;
	    double dist2 = this_point->distance2(*other->point);
	    if (dist2 < this_point->neighbour_dist2) {
	      this_point->neighbour_dist2 = dist2;
	      this_point->neighbour       = other->point;
	    }
	  }
	}
      }

      // for any non-zero review flag we'll have to update the heap
      _heap->update(_ID(this_point), this_point->neighbour_dist2);
    }

    // "delabel" the point
    this_point->review_flag = 0; 

  }

}

//----------------------------------------------------------------------
unsigned int ClosestPair2D::insert(const Coord2D & new_coord) {

  // get hold of a point
  assert(_available_points.size() > 0);
  Point * new_point = _available_points.top();
  _available_points.pop();

  // set the point's coordinate
  new_point->coord = new_coord;
  
  // now find it's neighbour in the search tree
  _insert_into_search_tree(new_point);

  // sort out other points that may have been affected by this, 
  // and/or for which the heap needs to be updated
  _deal_with_points_to_review();

  // 
  return _ID(new_point);
}

//----------------------------------------------------------------------
unsigned int ClosestPair2D::replace(unsigned int ID1, unsigned int ID2, 
				    const Coord2D & position) {
  
  // deletion from tree...
  Point * point_to_remove = & (_points[ID1]);
  _remove_from_search_tree(point_to_remove);

  point_to_remove = & (_points[ID2]);
  _remove_from_search_tree(point_to_remove);

  // insertion into tree
  // get hold of a point
  Point * new_point = _available_points.top();
  _available_points.pop();

  // set the point's coordinate
  new_point->coord = position;
  
  // now find it's neighbour in the search tree
  _insert_into_search_tree(new_point);

  // the above statement labels certain points as needing "review" --
  // deal with them...
  _deal_with_points_to_review();

  //
  return _ID(new_point);

}


//----------------------------------------------------------------------
void ClosestPair2D::replace_many(
                  const std::vector<unsigned int> & IDs_to_remove,
		  const std::vector<Coord2D> & new_positions,
		  std::vector<unsigned int> & new_IDs) {

  // deletion from tree...
  for (unsigned int i = 0; i < IDs_to_remove.size(); i++) {
    _remove_from_search_tree(& (_points[IDs_to_remove[i]]));
  }

  // insertion into tree
  new_IDs.resize(0);
  for (unsigned int i = 0; i < new_positions.size(); i++) {
    Point * new_point = _available_points.top();
    _available_points.pop();
    // set the point's coordinate
    new_point->coord = new_positions[i];
    // now find it's neighbour in the search tree
    _insert_into_search_tree(new_point);
    // record the ID
    new_IDs.push_back(_ID(new_point));
  }

  // the above statement labels certain points as needing "review" --
  // deal with them...
  _deal_with_points_to_review();

}


//----------------------------------------------------------------------
void ClosestPair2D::_insert_into_search_tree(Point * new_point) {

  // this point will have to have it's heap entry reviewed...
  _set_label(new_point, _review_heap_entry);

  // set the current distance to "infinity"
  new_point->neighbour_dist2 = numeric_limits<double>::max();
  
  // establish how far we will be searching;
  unsigned int CP_range = min(_cp_search_range, size()-1);

  for (unsigned ishift = 0; ishift < _nshift; ishift++) {
    // create the shuffle
    Shuffle new_shuffle;
    _point2shuffle(*new_point, new_shuffle, _shifts[ishift]);
    
    // insert it into the tree
    circulator new_circ = _trees[ishift]->insert(new_shuffle);
    new_point->circ[ishift] = new_circ;

    // now get hold of the right and left edges of the region we will be
    // looking at (cf CCN28-43)
    circulator right_edge = new_circ; right_edge++;
    circulator left_edge  = new_circ;
    for (unsigned int i = 0; i < CP_range; i++) {left_edge--;}

    // now 
    do {
      Point * left_point  = left_edge->point;
      Point * right_point = right_edge->point;

      // see if the new point is closer to the left-edge than the latter's
      // current neighbour
      double new_dist2 = left_point->distance2(*new_point);
      if (new_dist2 < left_point->neighbour_dist2) {
	left_point->neighbour_dist2 = new_dist2;
	left_point->neighbour       = new_point;
	_add_label(left_point, _review_heap_entry);
      }

      // see if the right-point is closer to the new point than it's current
      // neighbour
      new_dist2 = new_point->distance2(*right_point);
      if (new_dist2 < new_point->neighbour_dist2) {
	new_point->neighbour_dist2 = new_dist2;
	new_point->neighbour = right_point;
      }

      // if the right-edge point was the left-edge's neighbour, then
      // then it's just gone off-radar and the left-point will need to
      // have its neighbour recalculated [actually, this is overdoing
      // it a little, since right point may be an less "distant"
      // (circulator distance) in one of the other shifts -- but not
      // sure how to deal with this...]
      if (left_point->neighbour == right_point) {
	_add_label(left_point, _review_neighbour);
      }

      // shift the left and right edges until left edge hits new_circ
      right_edge++;
    } while (++left_edge != new_circ);
  }
}

FASTJET_END_NAMESPACE

