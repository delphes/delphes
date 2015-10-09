//FJSTARTHEADER
// $Id: DnnPlane.cc 3918 2015-07-03 14:19:13Z salam $
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


#ifndef DROP_CGAL // in case we do not have the code for CGAL

#include<set>
#include<list>
#include "fastjet/internal/DnnPlane.hh"

using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

const double DnnPlane::DISTANCE_FOR_CGAL_CHECKS=1.0e-12;  


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

    // check if we are dealing with coincident vertices
    int coinciding_index = _CheckIfVertexPresent(sv.vertex, i);
    if (coinciding_index == i){
      // we need to associate an index to each vertex -- thus when we get
      // a vertex (e.g. as a nearest neighbour) from CGAL, we will be
      // able to figure out which particle it corresponded to.
      sv.vertex->info() = sv.coincidence = i;
    } else {
      //cout << "  coincident with " << coinciding_index << endl;
      // the new vertex points to the already existing one and we
      // record the coincidence
      //
      // Note that we must not only set the coincidence of the
      // currently-added particle, the one it coincides with also
      // needs be updated (taking into account that it might already
      // coincide with another one)
      //
      // An example may help. Say coinciding_index = i1 and we're adding i2==i.
      // Then _sv[i2].coincidence = i1; _sv[i1].coincidence = i2. In both
      // cases sv.vertex->info() == i1;
      //
      // Later on we add i3; we find out that its coinciding index is i1;
      // so we set _sv[i3].coincidence = i2 and sv[i1].coincidence = i3. 
      //
      // This gives us the structure 
      //  _supervertex[i1].coincidence == in
      //  _supervertex[i2].coincidence == i1
      //  ...
      //  _supervertex[in].coincidence == in-1
      //
      sv.coincidence = _supervertex[coinciding_index].coincidence; // handles cases with previous coincidences
      _supervertex[coinciding_index].coincidence = i;
    }
      
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
int DnnPlane::_CheckIfVertexPresent(
	const Vertex_handle & vertex, const int its_index) {
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
    if (_crash_on_coincidence){
      ostringstream err;
      err << "Error: DnnPlane::_CheckIfVertexPresent"
	  << "Point "<<its_index<<" coincides with point "
	  <<vertex->info().val() << endl;
      throw DnnError(err.str());
    }
    return vertex->info().val();
  }

  return its_index;
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

  if (_verbose) cout << "Starting  DnnPlane::RemoveAndAddPoints" << endl;

  // build set of UNION of Voronoi neighbours of a pair of nearest
  // neighbours
  set<int> NeighbourUnion;
  // later on it will be convenient to have access to a set (rather
  // than vector) of indices being removed
  set<int> indices_removed;

  // for each of the indices to be removed add the voronoi
  // neighbourhood to the NeighbourUnion set as well as the coinciding
  // points that had the current point as coincidence before.
  for (size_t ir = 0; ir < indices_to_remove.size(); ir++) {
    int index = indices_to_remove[ir];
    indices_removed.insert(index);
    if (_verbose) cout << "  scheduling point " << index << " for removal" << endl;

    if (_supervertex[index].coincidence != index){
      // we have a coincidence
      //
      // The only one of the coincident points that has to be
      // inserted in the neighbourhood list (and thus updated) is the
      // one that has 'index' as coincidence.
      int new_index = _supervertex[index].coincidence;
      while (_supervertex[new_index].coincidence != index)
	new_index = _supervertex[new_index].coincidence;
      if (_verbose) cout << "  inserted coinciding " << new_index << " to neighbours union" << endl;
      NeighbourUnion.insert(new_index);

      // if this is the point among the coiciding ones that holds the
      // CGAL vertex, then also insert the CGAL neighbours, otherwise
      // just skip that step.
      if (index != _supervertex[index].vertex->info().val()) continue;
    } 

    // have a circulators that will go round the Voronoi neighbours of
    // _supervertex[index1].vertex
    Vertex_circulator vc = _TR.incident_vertices(_supervertex[index].vertex);
    Vertex_circulator done = vc;
    if (vc != NULL){ // a safety check in case there is no Voronoi
		     // neighbour (which may happen e.g. if we just
		     // have a bunch of coincident points)
      do  {
	// if a neighbouring vertex is not the infinite vertex, then add it
	// to our union of neighbouring vertices.
	if (_verbose) cout << "examining " << vc->info().val() << endl;
	if (vc->info().val() != INFINITE_VERTEX) {
	  // NB: from it=1 onwards occasionally it might already have
	  // been inserted -- but double insertion still leaves only one
	  // copy in the set, so there's no problem
	  NeighbourUnion.insert(vc->info().val());
	  if (_verbose) cout << "  inserted " << vc->info().val() << " to neighbours union" << endl;
	} 
      } while (++vc != done);
    }
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
    if (_verbose) cout << "  removing " << index << endl;

    // NeighbourUnion should not contain the points to be removed
    // (because later we will assume they still exist).
    NeighbourUnion.erase(indices_to_remove[ir]);

    // first deal with  coincidences
    if (_supervertex[index].coincidence != index){
      int new_index = _supervertex[index].coincidence;

      // if this is the point among the coiciding ones that "owns" the
      // CGAL vertex we need to re-label the CGAL vertex so that it
      // points to the coincident particle and set the current one to
      // NULL
      //
      // This can be done only on the first point as they all share
      // the same value
      //
      // Note that this has to be done before the following step since
      // it will alter the coincidence information
      if (index == _supervertex[index].vertex->info().val())
	_supervertex[new_index].vertex->info() = new_index;

      // we need to browse the coincidences until we end the loop, at
      // which point we reset the coincidence of the point that has
      // the current one as a coincidence
      while (_supervertex[new_index].coincidence != index)
	new_index = _supervertex[new_index].coincidence;
      _supervertex[new_index].coincidence = _supervertex[index].coincidence;

      // remove the coincidence on the point being removed and mark it
      // as removed
      _supervertex[index].coincidence = index;
      _supervertex[index].vertex = NULL;

      continue;
    }

    // points to be removed should also be eliminated from the
    // triangulation and the supervertex structure should be updated
    // to reflect the fact that the points are no longer valid.
    _TR.remove(_supervertex[index].vertex);
    if (_verbose) cout << "DnnPlane about to set _supervertex["<< index<<"].vertex to NULL" << endl;
    _supervertex[index].vertex = NULL;
    if (_verbose) cout << "                 value is " << (_is_not_null(_supervertex[index].vertex)) << endl;
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
  //if (indices_to_remove.size() > 0) { // GS: use NeighbourUnion instead
                                        //     (safe also in case of coincidences)
  if (NeighbourUnion.size() > 0) {
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
    if (_verbose) cout << "  adding " << index << " at "
                       << points_to_add[ia].first<< " " << points_to_add[ia].second << endl;

    //if (indices_to_remove.size() > 0) {
    if (NeighbourUnion.size() > 0) {
      // be careful of using face (for location hinting) only when it exists
      _supervertex[index].vertex = _TR.insert(Point(points_to_add[ia].first, 
				  points_to_add[ia].second),face);}
    else { 
      _supervertex[index].vertex = _TR.insert(Point(points_to_add[ia].first, 
						    points_to_add[ia].second));
    }

    // check if this leads to a coincidence
    int coinciding_index = _CheckIfVertexPresent(_supervertex[index].vertex, index);
    if (coinciding_index == index){
      // we need to associate an index to each vertex -- thus when we get
      // a vertex (e.g. as a nearest neighbour) from CGAL, we will be
      // able to figure out which particle it corresponded to.
      _supervertex[index].vertex->info() = _supervertex[index].coincidence = index;
    } else {
      if (_verbose) cout << "  coinciding with vertex " << coinciding_index << endl;
      // the new vertex points to an already existing one and we
      // record the coincidence
      //
      // we also update the NN of the coinciding particle (to avoid
      // having to loop over the list of coinciding neighbours later)
      // This is done first as it allows us to check if this is a new
      // coincidence or a coincidence added to a particle that was
      // previously "alone"
      _supervertex[coinciding_index].NNindex = index;
      _supervertex[coinciding_index].NNdistance = 0.0;
      indices_of_updated_neighbours.push_back(coinciding_index);

      // Note that we must not only set the coincidence of the
      // currently-added particle, the one it coincides with also
      // needs be updated (taking into account that it might already
      // coincide with another one)
      _supervertex[index].coincidence = _supervertex[coinciding_index].coincidence; // handles cases with previous coincidences
      _supervertex[coinciding_index].coincidence = index;

    }
    
    // first find nearest neighbour of "newpoint" (shorthand for
    // _supervertex[index].vertex); while we're at it, for each of the
    // voronoi neighbours, "D", of newpoint, examine whether newpoint is
    // closer to "D" than D's current nearest neighbour -- when this
    // occurs, put D into indices_of_updated_neighbours.
    // 
    // manually put newpoint on indices_of_updated_neighbours
    indices_of_updated_neighbours.push_back(index);
    _SetAndUpdateNearest(index, indices_of_updated_neighbours);

    //cout << "Added: " << setprecision(20) << " (" 
    //     << points_to_add[ia].first << "," << points_to_add[ia].second
    //     << ") with index " << index << endl;
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

  if (_verbose) cout << "Leaving  DnnPlane::RemoveAndAddPoints" << endl;
}

//----------------------------------------------------------------------
/// Determines the index and distance of the nearest neighbour to 
/// point j and puts the information into the _supervertex entry for j.
void DnnPlane::_SetNearest (const int j) {
  // first deal with the cases where we have a coincidence
  if (_supervertex[j].coincidence != j){
    _supervertex[j].NNindex = _supervertex[j].coincidence;
    _supervertex[j].NNdistance = 0.0;
    return;
  }

  // The code below entirely uses CGAL distance comparisons to compute
  // the nearest neighbour. It has the mais drawback to induice a
  // 10-20% time penalty so we switched to our own comparison (which
  // only turns to CGAL for dangerous situations)
  //
  //  Vertex_handle current = _supervertex[j].vertex;
  //  Vertex_circulator vc = _TR.incident_vertices(current);
  //  Vertex_circulator done = vc;
  //  Vertex_handle nearest = _TR.infinite_vertex();
  //  double mindist = HUGE_DOUBLE;
  //
  //   // when there is only one finite point left in the triangulation, 
  //   // there are no triangles. Presumably this is why voronoi returns
  //   // NULL for the incident vertex circulator. Check if this is
  //   // happening before circulating over it... (Otherwise it crashes
  //   // when looking for neighbours of last point)
  //   if (vc != NULL){
  //     // initialise the nearest vertex handle to the first incident
  //     // vertex that is not INFINITE_VERTEX
  //     while (vc->info().val() == INFINITE_VERTEX){
  //       vc++;
  //       if (vc==done) break; // if vc==done, then INFINITE_VERTEX is the
  // 			   // only element in the neighbourhood
  //     }
  // 
  //     // if there is just the infinite vertex, we have vc->info().val()
  //     // == INFINITE_VERTEX and nothing has to be done
  //     // otherwise, use the current vc as an initialisation
  //     if (vc->info().val() != INFINITE_VERTEX){
  //       nearest = vc; // initialisation to the first non-infinite vertex
  // 
  //       // and loop over the following ones
  //       while (++vc != done){
  // 	// we should not compare with the infinite vertex
  // 	if (vc->info().val() == INFINITE_VERTEX) continue;
  // 
  // 	if (_verbose) cout << current->info().val() << " " << vc->info().val() << endl; 
  // 	// use CGAL's distance comparison to check if 'vc' is closer to
  // 	// 'current' than the nearest so far (we include the == case for
  // 	// safety though it should not matter in this precise case)
  // 	if (CGAL::compare_distance_to_point(current->point(), vc->point(), nearest->point())!=CGAL::LARGER){
  // 	  nearest = vc;
  // 	  if (_verbose) cout << "nearer";
  // 	}
  //       }
  // 
  //       // now compute the distance
  //       //
  //       // Note that since we're always using CGAL to compare distances
  //       // (and never the distance computed using _euclid_distance) we
  //       // should not worry about rounding errors in mindist
  //       mindist = _euclid_distance(current->point(), nearest->point());
  //     }
  //   }
  //
  //  // set j's supervertex info about nearest neighbour
  //  _supervertex[j].NNindex = nearest->info().val();
  //  _supervertex[j].NNdistance = mindist;

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

      // check if j is closer to vc than vc's currently registered
      // nearest neighbour (and update things if it is)
      if (_is_closer_to(current->point(), vc->point(), nearest, dist, mindist)){
	nearest = vc; 
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
/// necessary updates the nearest-neighbour info of Voronoi neighbours
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
///     use one single circulation over voronoi neighbours to find the
///     nearest neighbour and to update the voronoi neighbours if need
///     be.
void DnnPlane::_SetAndUpdateNearest(
			  const int j, 
			  vector<int> & indices_of_updated_neighbours) {
  //cout << "SetAndUpdateNearest for point " << j << endl;
  // first deal with coincidences
  if (_supervertex[j].coincidence != j){
    _supervertex[j].NNindex = _supervertex[j].coincidence;
    _supervertex[j].NNdistance = 0.0;
    //cout << "  set to coinciding point " << _supervertex[j].coincidence << endl;
    return;
  }

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

      // update the mindist if we are closer than anything found so far
      if (_is_closer_to(current->point(), vc->point(), nearest, dist, mindist)){
	nearest = vc; 
      	if (_verbose) cout << "nearer ";
      } 

      // find index corresponding to vc for easy manipulation
      int vcindx = vc->info().val();
      if (_verbose) cout << vc->point() << "; "<< dist << endl;

      if (_is_closer_to_with_hint(vc->point(), current->point(), 
				  _supervertex[_supervertex[vcindx].NNindex].vertex,
				  dist, _supervertex[vcindx].NNdistance)){
	if (_verbose) cout << vcindx << "'s NN becomes " << current->info().val() << endl;
	_supervertex[vcindx].NNindex = j;
	indices_of_updated_neighbours.push_back(vcindx);
      }

      // original code without the use of CGAL distance in potentially
      // dangerous cases
      //
      // // check if j is closer to vc than vc's currently registered
      // // nearest neighbour (and update things if it is)
      // //
      // // GS: originally, the distance test below was a strict <. It
      // //     has to be <= because if the two distances are ==, it is
      // //     possible that the old NN is no longer connected to vc in
      // //     the triangulation, and we are sure that the newly
      // //     inserted point (j) is (since we loop over j's
      // //     neighbouring points in the triangulation).
      // if (dist <= _supervertex[vcindx].NNdistance) {
      // 	if (_verbose) cout << vcindx << "'s NN becomes " << current->info().val() << endl;
      // 	_supervertex[vcindx].NNdistance = dist;
      // 	_supervertex[vcindx].NNindex = j;
      // 	indices_of_updated_neighbours.push_back(vcindx);
      // }
    }
  } while (++vc != done); // move on to next Voronoi neighbour
  // set j's supervertex info about nearest neighbour
  //cout << "  set to point " << nearest->info().val() << endl;
  _supervertex[j].NNindex = nearest->info().val();
  _supervertex[j].NNdistance = mindist;
}

FASTJET_END_NAMESPACE

#endif //  DROP_CGAL
