//FJSTARTHEADER
// $Id: Dnn2piCylinder.cc 3433 2014-07-23 08:17:03Z salam $
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
#include <set>
#include "fastjet/internal/Dnn2piCylinder.hh"
using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//----------------------------------------------------------------------
/// initialiser...
Dnn2piCylinder::Dnn2piCylinder(
	const vector<EtaPhi> & input_points, 
	const bool & ignore_nearest_is_mirror,
	const bool & verbose) {
  
  _verbose = verbose;
  _ignore_nearest_is_mirror = ignore_nearest_is_mirror;
  vector<EtaPhi> plane_points;
  vector<int>    plane_point_indices(input_points.size());
  //plane_points.reserve(2*input_points.size());

  for (unsigned int i=0; i < input_points.size(); i++) {
    _RegisterCylinderPoint(input_points[i], plane_points);
    plane_point_indices[i] = i;
  }
  
  if (_verbose) cout << "============== Preparing _DNN" << endl;
  _DNN = new DnnPlane(plane_points, verbose);


  vector<int> updated_point_indices; // we'll not use information from this
  _CreateNecessaryMirrorPoints(plane_point_indices,updated_point_indices);
}


//----------------------------------------------------------------------
/// Actions here are similar to those in the
/// Dnn3piCylinder::_RegisterCylinderPoint case, however here we do
/// NOT create the mirror point -- instead we initialise the structure
/// as if there were no need for the mirror point.
///
/// ADDITIONALLY push the cylinder_point onto the vector plane_points.
void Dnn2piCylinder::_RegisterCylinderPoint (const EtaPhi & cylinder_point,
					     vector<EtaPhi> & plane_points) {
  double phi = cylinder_point.second;
  assert(phi >= 0.0 && phi < 2*pi);
  
  // do main point
  MirrorVertexInfo mvi;
  mvi.main_index = _cylinder_index_of_plane_vertex.size();
  _cylinder_index_of_plane_vertex.push_back(_mirror_info.size());
  plane_points.push_back(cylinder_point);
  mvi.mirror_index = INEXISTENT_VERTEX;
  
  // 
  _mirror_info.push_back(mvi);
}



//----------------------------------------------------------------------
/// For each plane point specified in the vector plane_indices,
/// establish whether there is a need to create a mirror point
/// according to the following criteria:
///
/// . phi < pi
/// . mirror does not already exist
/// . phi < NearestNeighbourDistance 
///   (if this is not true then there is no way that its mirror point
///   could have a nearer neighbour).
///
/// If conditions all hold, then create the mirror point, insert it
/// into the _DNN structure, adjusting any nearest neighbours, and
/// return the list of plane points whose nearest neighbours have
/// changed (this will include the new neighbours that have just been
/// added)
void Dnn2piCylinder::_CreateNecessaryMirrorPoints(
			  const vector<int> & plane_indices,
			  vector<int> & updated_plane_points) {

  vector<EtaPhi> new_plane_points;

  for (size_t i = 0; i < plane_indices.size(); i++) {

    int ip = plane_indices[i]; // plane index
    EtaPhi position = _DNN->etaphi(ip);
    double phi = position.second;

    //BAD // require phi < pi
    //BAD if (phi >= pi) {continue;}

    // require absence of mirror
    int ic = _cylinder_index_of_plane_vertex[ip];
    if (_mirror_info[ic].mirror_index != INEXISTENT_VERTEX) {continue;}

    //printf("%3d %3d %10.5f %10.5f %3d\n",ic, ip, phi, 
    //	   _DNN->NearestNeighbourDistance(ip),_DNN->NearestNeighbourIndex(ip));


    // check that we are sufficiently close to the border --
    // i.e. closer than nearest neighbour distance. But RECALL:
    // nearest neighbourDistance actually returns a squared distance
    // (this was really stupid on our part -- caused considerable loss
    // of time ... )
    double nndist = _DNN->NearestNeighbourDistance(ip); 
    if (phi*phi >= nndist && (twopi-phi)*(twopi-phi) >= nndist) {continue;}

    // now proceed to prepare the point for addition
    new_plane_points.push_back(_remap_phi(position));
    _mirror_info[ic].mirror_index = _cylinder_index_of_plane_vertex.size();
    _cylinder_index_of_plane_vertex.push_back(ic);
  }

  vector<int> indices_to_remove; // empty
  vector<int> indices_added;     // will be filled as result of call
  _DNN->RemoveAndAddPoints(indices_to_remove,new_plane_points,indices_added, 
			   updated_plane_points);

  // occasionally, addition of points might cause a change in the
  // nearest neighbour of a point in the 0--pi range? (But should be
  // impossible -- we add points beyond 2pi, so they can only be
  // nearest neighbours of points in the range pi--2pi; there is one
  // exception -- the nearest neighbour of one's self -- but in that
  // case this will already have been discovered, so there should be
  // nothing left to do). 

  // To be on the safe side, check to see if we have updated points
  // with phi<pi and no current mirror point. BUT: this check, as
  // written, only makes sense when the mirror image was created only
  // beyond 2pi, which is no longer the case. Only apparent
  // alternative is to run separate insertions for beyond 2pi and
  // below phi=0, with separate checks in the two cases. But, given
  // that across a large number of recombinations and events in the
  // old method (single mirror), we never ran into this problem, it's
  // probably safe to assume that the arguments given above are OK! So
  // comment out the check...
  //NOTNEEDED for (size_t i = 0; i < updated_plane_points.size(); i++) {
  //NOTNEEDED   int ip = updated_plane_points[i];
  //NOTNEEDED   double phi  = _DNN->phi(ip);
  //NOTNEEDED   int ic = _cylinder_index_of_plane_vertex[ip];
  //NOTNEEDED   assert (!(phi < pi && _mirror_info[ic].mirror_index == INEXISTENT_VERTEX));
  //NOTNEEDED }
  // alternative recursive code
  //vector<int> extra_updated_plane_points;
  //_CreateNecessaryMirrorPoints(updated_plane_points,extra_updated_plane_points)
  //updated_plane_points.push_back(extra_updated_plane_points);
}



//----------------------------------------------------------------------
/// insertion and removal of points
void Dnn2piCylinder::RemoveAndAddPoints(const vector<int> & indices_to_remove,
				const vector<EtaPhi> & points_to_add,
				vector<int> & indices_added,
				vector<int> & indices_of_updated_neighbours) {

  // translate from "cylinder" indices of points to remove to the
  // plane indices of points to remove, bearing in mind that sometimes
  // there are multple plane points to remove.
  vector<int> plane_indices_to_remove;
  for (unsigned int i=0; i < indices_to_remove.size(); i++) {
    MirrorVertexInfo * mvi;
    mvi = & _mirror_info[indices_to_remove[i]];
    plane_indices_to_remove.push_back(mvi->main_index);
    if (mvi->mirror_index != INEXISTENT_VERTEX) {
      plane_indices_to_remove.push_back(mvi->mirror_index);
    }
  }

  // given "cylinder" points to add get hold of the list of
  // plane-points to add.
  vector<EtaPhi> plane_points_to_add;
  indices_added.clear();
  for (unsigned int i=0; i < points_to_add.size(); i++) {
    indices_added.push_back(_mirror_info.size());
    _RegisterCylinderPoint(points_to_add[i], plane_points_to_add);
  }

  // now get the hard work done (note that we need to supply the
  // plane_indices_added vector but that we will not actually check
  // its contents in any way -- the indices_added that is actually
  // returned has been calculated above).
  vector<int> updated_plane_neighbours, plane_indices_added;
  _DNN->RemoveAndAddPoints(plane_indices_to_remove, plane_points_to_add,
			     plane_indices_added, updated_plane_neighbours);

  vector<int> extra_updated_plane_neighbours;
  _CreateNecessaryMirrorPoints(updated_plane_neighbours,
			       extra_updated_plane_neighbours);

  // extract, from the updated_plane_neighbours, and
  // extra_updated_plane_neighbours, the set of cylinder neighbours
  // that have changed
  set<int> index_set;
  unsigned int i;
  for (i=0; i < updated_plane_neighbours.size(); i++) {
    index_set.insert(
       _cylinder_index_of_plane_vertex[updated_plane_neighbours[i]]);}
  for (i=0; i < extra_updated_plane_neighbours.size(); i++) {
    index_set.insert(
       _cylinder_index_of_plane_vertex[extra_updated_plane_neighbours[i]]);}

  // decant the set into the vector that needs to be returned
  indices_of_updated_neighbours.clear();
  for (set<int>::iterator iter = index_set.begin(); 
       iter != index_set.end(); iter++) {
    indices_of_updated_neighbours.push_back(*iter);
  }
}

FASTJET_END_NAMESPACE

#endif //  DROP_CGAL 
