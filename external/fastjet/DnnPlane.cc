//STARTHEADER
// $Id$
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


#ifndef DROP_CGAL // in case we do not have the code for CGAL

#include<set>
#include<list>
#include "fastjet/internal/DnnPlane.hh"
using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


/// Initialiser from a set of points on an Eta-Phi plane, where both
/// eta and phi can have arbitrary ranges
DnnPlane::DnnPlane(const vector<EtaPhi> & input_points, 
		   const bool & verbose ) {

  _verbose = verbose;
  int n = input_points.size();
  
  // construct Voronoi diagram in such a way as to get the vertex handles
  // and remember to set CGAL info with the index of the vertex
  SuperVertex sv;
  for (int i = 0; i < n; i++) {
    sv.vertex = 
       _TR.insert(Point(input_points[i].first, input_points[i].second));

    // we are not up to dealing with coincident vertices, so make 
    // sure the user knows!
    _CrashIfVertexPresent(sv.vertex, i);
    
    // we need to assicate an index to each vertex -- thus when we get
    // a vertex (e.g. as a nearest neighbour) from CGAL, we will be
    // able to figure out which particle it corresponded to.
    sv.vertex->info() = i;
    _supervertex.push_back(sv);    
  }

  // label infinite vertex info with negative index 
  _TR.infinite_vertex()->info() = INFINITE_VERTEX;

  // set up the structure that holds nearest distances and neighbours
  for (int j = 0; j < n; j++) {_SetNearest(j);}

}


//----------------------------------------------------------------------
/// Crashes if the given vertex handle already exists. Otherwise
/// it does the bookkeeping for future such tests
void DnnPlane::_CrashIfVertexPresent(
	const Vertex_handle & vertex, const int & its_index) {
  if (!_crash_on_coincidence) return;

  // vertices that do not have the same geometric position as any
  // other vertex so far added have info().val() == NEW_VERTEX -- this
  // is ensured by the InitialisedInt class, which forms the "info"
  // part of our
  // CGAL::Triangulation_vertex_base_with_info_2<InitialisedInt,K>
  //
  // If the vertex coincides with one that already exists, then
  // info().val() it's info().val() will have been updated (in
  // DNN:DNN) to be equal to a vertex "index".
  if (vertex->info().val() != NEW_VERTEX) {
    ostringstream err;
    err << "ERROR in DnnPlane::_CrashIfVertexPresent"
	 <<endl << "Point "<<its_index<<" coincides with point "
	 <<vertex->info().val() << endl;
    throw DnnError(err.str());
  } 
}


//----------------------------------------------------------------------
/// remove the points labelled by the vector indices_to_remove, and
/// add the points specified by the vector points_to_add
/// (corresponding indices will be calculated automatically); the
/// idea behind this routine is that the points to be added will
/// somehow be close to the one or other of the points being removed
/// and this can be used by the implementation to provide hints for
/// inserting the new points in whatever structure it is using.  In a
/// kt-algorithm the points being added will be a result of a
/// combination of the points to be removed -- hence the proximity
/// is (more or less) guaranteed.
void DnnPlane::RemoveAndAddPoints(
			  const vector<int> & indices_to_remove,
			  const vector<EtaPhi> & points_to_add,
			  vector<int> & indices_added,
			  vector<int> & indices_of_updated_neighbours) {


  // build set of UNION of Voronoi neighbours of a pair of nearest
  // neighbours
  set<int> NeighbourUnion;
  // later on it will be convenient to have access to a set (rather
  // than vector) of indices being removed
  set<int> indices_removed;

  // for each of the indices to be removed add the voronoi neighbourhood to
  // the NeighbourUnion set.
  for (size_t ir = 0; ir < indices_to_remove.size(); ir++) {
    int index = indices_to_remove[ir];
    indices_removed.insert(index);
    if (_verbose) cout << " Starting  RemoveAndAddPoints" << endl; 
    if (_verbose) cout << " point " << index << endl;			 
    // have a circulators that will go round the Voronoi neighbours of
    // _supervertex[index1].vertex
    Vertex_circulator vc = _TR.incident_vertices(_supervertex[index].vertex);
    Vertex_circulator done = vc;
    do  {
      // if a neighbouring vertex not the infinite vertex, then add it
      // to our union of neighbouring vertices.
      if (_verbose) cout << "examining " << vc->info().val() << endl;
      if (vc->info().val() != INFINITE_VERTEX) {
	// NB: from it=1 onwards occasionally it might already have
	// been inserted -- but double insertion still leaves only one
	// copy in the set, so there's no problem
	NeighbourUnion.insert(vc->info().val());
	if (_verbose) cout << "inserted " << vc->info().val() << endl;
      } 
    } while (++vc != done);
  }
  
  if (_verbose) {
    set<int>::iterator it = NeighbourUnion.begin();
    cout << "Union of neighbours of combined points" << endl;
    for ( ; it != NeighbourUnion.end(); ++it ) {
      cout << *it << endl;
    }
  }

  // update set, triangulation and supervertex info
  for (size_t ir = 0; ir < indices_to_remove.size(); ir++) {
    int index = indices_to_remove[ir];

    // NeighbourUnion should not contain the points to be removed
    // (because later we will assume they still exist).
    NeighbourUnion.erase(indices_to_remove[ir]);
    
    // points to be removed should also be eliminated from the
    // triangulation and the supervertex structure should be updated
    // to reflect the fact that the points are no longer valid.
    _TR.remove(_supervertex[index].vertex);
    _supervertex[index].vertex = NULL;
  }

  // add new point: give a "hint" to the inserter that
  // the new point should be added close to old points -- the easiest way 
  // of getting this is to take a point from the NeighbourUnion AFTER we have
  // removed point1, point2, and to get one of its incident faces.
  // 
  // This hinting improves speed by c. 25% for 10^4 points because it
  // avoids the need for a costly (sqrt{N}) location step (at least
  // with a non-hierarchical triangulation -- with a hierarchical one,
  // this step could be done away with, though there will still be a
  // cost of O(ln N) to pay.
  // 
  // For some reason inserting the point before the two removals
  // slows things down by c. 25%. This importance of the order
  // is not understood.
  //
  // At some point it might be worth trying to select the "nearest"
  // of the various points in the neighbour union to avoid large 
  // steps in cases where we have 0..2pi periodicity and the first member
  // of the neighbour union happens to be on the wrong side.
  Face_handle face;
  if (indices_to_remove.size() > 0) {
    // face can only be found if there were points to remove in first place
    face = _TR.incident_faces(
   	                   _supervertex[*NeighbourUnion.begin()].vertex);}
  // make sure the output arrays are empty
  indices_added.clear();
  indices_of_updated_neighbours.clear();
  for (size_t ia = 0; ia < points_to_add.size(); ia++) {
    SuperVertex sv;
    _supervertex.push_back(sv);
    int index = _supervertex.size()-1;
    indices_added.push_back(index);

    if (indices_to_remove.size() > 0) {
      // be careful of using face (for location hinting) only when it exists
      _supervertex[index].vertex = _TR.insert(Point(points_to_add[ia].first, 
				  points_to_add[ia].second),face);}
    else { 
      _supervertex[index].vertex = _TR.insert(Point(points_to_add[ia].first, 
						    points_to_add[ia].second));
    }
    // we are not up to dealing with coincident vertices, so make 
    // sure the user knows!
    _CrashIfVertexPresent(_supervertex[index].vertex, index);
    _supervertex[index].vertex->info() = index;
    
    // first find nearest neighbour of "newpoint" (shorthand for
    // _supervertex[index].vertex); while we're at it, for each of the
    // voronoi neighbours, "D", of newpoint, examine whether newpoint is
    // closer to "D" than D's current nearest neighbour -- when this
    // occurs, put D into indices_of_updated_neighbours.
    // 
    // manually put newpoint on indices_of_updated_neighbours
    indices_of_updated_neighbours.push_back(index);
    _SetAndUpdateNearest(index, indices_of_updated_neighbours);
  }

  // for Voronoi neighbours j of any of the removed points for which
  // one of those removed points was the nearest neighbour,
  // redetermine the nearest neighbour of j and add j onto the vector
  // of indices_of_updated_neighbours.
  set<int>::iterator it2 = NeighbourUnion.begin();
  for ( ; it2 != NeighbourUnion.end(); ++it2 ) {
    int j = *it2;
    // the if avoids the vertex at infinity, which gets a negative index
    if( j != INFINITE_VERTEX ) {
      // this is where we check if the nearest neighbour of j was one
      // of the removed points
      if (indices_removed.count(_supervertex[j].NNindex)) {
	if (_verbose) cout << "j " << j << endl;
	_SetNearest(j);
	indices_of_updated_neighbours.push_back(j);
	if (_verbose) cout << "NN of " << j << " : " 
			  << _supervertex[j].NNindex
	                  << ", dist = " << _supervertex[j].NNdistance <<endl;
      }
    }
  }

}


//----------------------------------------------------------------------
/// Determines the index and distance of the nearest neighbour to 
/// point j and puts the information into the _supervertex entry for j.
void DnnPlane::_SetNearest (const int & j) {
  Vertex_handle current = _supervertex[j].vertex;
  Vertex_circulator vc = _TR.incident_vertices(current);
  Vertex_circulator done = vc;
  double dist;
  double mindist = HUGE_DOUBLE; // change this to "HUGE" or max_double?
  Vertex_handle nearest = _TR.infinite_vertex();

  // when there is only one finite point left in the triangulation, 
  // there are no triangles. Presumably this is why voronoi returns
  // NULL for the incident vertex circulator. Check if this is
  // happening before circulating over it... (Otherwise it crashes
  // when looking for neighbours of last point)
  if (vc != NULL) do { 
    if ( vc->info().val() != INFINITE_VERTEX) {
      // find distance between j and its Voronoi neighbour (vc)
      if (_verbose) cout << current->info().val() << " " << vc->info().val() << endl;
      dist = _euclid_distance(current->point(), vc->point());
      // check if j is closer to vc than vc's currently registered
      // nearest neighbour (and update things if it is)
      if (dist < mindist) {
	mindist = dist; nearest = vc; 
	if (_verbose) cout << "nearer ";
      } 
      if (_verbose) cout << vc->point() << "; "<< dist << endl;
    }
  } while (++vc != done); // move on to next Voronoi neighbour
  // set j's supervertex info about nearest neighbour
  _supervertex[j].NNindex = nearest->info().val();
  _supervertex[j].NNdistance = mindist;
}

//----------------------------------------------------------------------
/// Determines and stores the nearest neighbour of j, and where
/// necessary update the nearest-neighbour info of Voronoi neighbours
/// of j;
///
/// For each voronoi neighbour D of j if the distance between j and D
/// is less than D's own nearest neighbour, then update the
/// nearest-neighbour info in D; push D's index onto 
/// indices_of_updated_neighbours
///
/// Note that j is NOT pushed onto indices_of_updated_neighbours --
/// if you want it there, put it there yourself.
///
/// NB: note that we have _SetAndUpdateNearest as a completely
///     separate routine from _SetNearest because we want to
///     use one single ciruclation over voronoi neighbours to find the
///     nearest neighbour and to update the voronoi neighbours if need
///     be.
void DnnPlane::_SetAndUpdateNearest(
			  const int & j, 
			  vector<int> & indices_of_updated_neighbours) {

  Vertex_handle current = _supervertex[j].vertex;
  Vertex_circulator vc = _TR.incident_vertices(current);
  Vertex_circulator done = vc;
  double dist;
  double mindist = HUGE_DOUBLE; // change this to "HUGE" or max_double?
  Vertex_handle nearest = _TR.infinite_vertex();

  // when there is only one finite point left in the triangulation, 
  // there are no triangles. Presumably this is why voronoi returns
  // NULL for the incident vertex circulator. Check if this is
  // happening before circulating over it... (Otherwise it crashes
  // when looking for neighbours of last point)
  if (vc != NULL) do { 
    if (vc->info().val() != INFINITE_VERTEX) {
      if (_verbose) cout << current->info().val() << " " << vc->info().val() << endl;
      // find distance between j and its Voronoi neighbour (vc)
      dist = _euclid_distance(current->point(), vc->point());
      // update the mindist if we are closer than anything found so far
      if (dist < mindist) {
	mindist = dist; nearest = vc; 
	if (_verbose) cout << "nearer ";
      } 
      // find index corresponding to vc for easy manipulation
      int vcindx = vc->info().val();
      if (_verbose) cout << vc->point() << "; "<< dist << endl;
      // check if j is closer to vc than vc's currently registered
      // nearest neighbour (and update things if it is)
      if (dist < _supervertex[vcindx].NNdistance) {
	if (_verbose) cout << vcindx << "'s NN becomes " << current->info().val() << endl;
	_supervertex[vcindx].NNdistance = dist;
	_supervertex[vcindx].NNindex = j;
	indices_of_updated_neighbours.push_back(vcindx);
      }
    }
  } while (++vc != done); // move on to next Voronoi neighbour
  // set j's supervertex info about nearest neighbour
  _supervertex[j].NNindex = nearest->info().val();
  _supervertex[j].NNdistance = mindist;
}

FASTJET_END_NAMESPACE

#endif //  DROP_CGAL
