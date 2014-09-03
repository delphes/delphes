//FJSTARTHEADER
// $Id: Dnn3piCylinder.cc 3433 2014-07-23 08:17:03Z salam $
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
#include "fastjet/internal/Dnn3piCylinder.hh"
using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//----------------------------------------------------------------------
/// initialiser...
Dnn3piCylinder::Dnn3piCylinder(
	const vector<EtaPhi> & input_points, 
	const bool & ignore_nearest_is_mirror,
	const bool & verbose) {
  
  _verbose = verbose;
  _ignore_nearest_is_mirror = ignore_nearest_is_mirror;
  vector<EtaPhi> plane_points;
  //plane_points.reserve(2*input_points.size());

  for (unsigned int i=0; i < input_points.size(); i++) {
    _RegisterCylinderPoint(input_points[i], plane_points);
  }

  if (_verbose) cout << "============== Preparing _DNN" << endl;
  _DNN = new DnnPlane(plane_points, verbose);
}


//----------------------------------------------------------------------
/// What on earth does this do?
///
/// Example: last true "cylinder" index was 15
///          last plane index was 23
/// 
/// Then: _cylinder_index_of_plane_vertex.size() = 24 and 
///       _mirror_info.size() = 16
///
/// IF cylinder_point's phi < pi then
///   create:  _mirror_info[16] = (main_index = 24, mirror_index=25) 
///            _cylinder_index_of_plane_vertex[24] = 16
///            _cylinder_index_of_plane_vertex[25] = 16
/// ELSE
///   create:  _mirror_info[16] = (main_index = 24, mirror_index=INEXISTENT..) 
///            _cylinder_index_of_plane_vertex[24] = 16
///
/// ADDITIONALLY push the cylinder_point (and if it exists the mirror
/// copy) onto the vector plane_points.
void Dnn3piCylinder::_RegisterCylinderPoint (const EtaPhi & cylinder_point,
					     vector<EtaPhi> & plane_points) {
  double phi = cylinder_point.second;
  assert(phi >= 0.0 && phi < 2*pi);
  
  // do main point
  MirrorVertexInfo mvi;
  mvi.main_index = _cylinder_index_of_plane_vertex.size();
  _cylinder_index_of_plane_vertex.push_back(_mirror_info.size());
  plane_points.push_back(cylinder_point);
  
  // do mirror point if need be
  if (phi < pi) {
    mvi.mirror_index = _cylinder_index_of_plane_vertex.size();
    _cylinder_index_of_plane_vertex.push_back(_mirror_info.size());
    plane_points.push_back(_remap_phi(cylinder_point));
  } else {
    mvi.mirror_index = INEXISTENT_VERTEX;
  }
  
  // 
  _mirror_info.push_back(mvi);
}


//----------------------------------------------------------------------
/// insertion and removal of points
void Dnn3piCylinder::RemoveAndAddPoints(const vector<int> & indices_to_remove,
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

  // extract, from the updated_plane_neighbours, the set of cylinder
  // neighbours that have changed
  set<int> index_set;
  unsigned int i;
  for (i=0; i < updated_plane_neighbours.size(); i++) {
    index_set.insert(
       _cylinder_index_of_plane_vertex[updated_plane_neighbours[i]]);}

  // decant the set into the vector that needs to be returned
  indices_of_updated_neighbours.clear();
  for (set<int>::iterator iter = index_set.begin(); 
       iter != index_set.end(); iter++) {
    indices_of_updated_neighbours.push_back(*iter);
  }
}


FASTJET_END_NAMESPACE

#endif //  DROP_CGAL 

