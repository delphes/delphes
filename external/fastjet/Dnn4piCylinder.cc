//FJSTARTHEADER
// $Id: Dnn4piCylinder.cc 3433 2014-07-23 08:17:03Z salam $
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
#include "fastjet/internal/Dnn4piCylinder.hh"
using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//----------------------------------------------------------------------
/// initialiser...
Dnn4piCylinder::Dnn4piCylinder(
	const vector<EtaPhi> & input_points, const bool & verbose) {
  
  _verbose = verbose;
  vector<EtaPhi> copied_points(input_points.size());
  for (unsigned int i=0; i < input_points.size(); i++) {
    double phi = input_points[i].second;
    assert(phi >= 0.0 && phi < 2*pi);
    copied_points[i] = _remap_phi(input_points[i]);
  }

  if (_verbose) cout << "============== Preparing _DNN1" << endl;
  _DNN1 = new DnnPlane(input_points, verbose);
  if (_verbose) cout << "============== Preparing _DNN2" << endl;
  _DNN2 = new DnnPlane(copied_points, verbose);
}


//----------------------------------------------------------------------
/// insertion and removal of points
void Dnn4piCylinder::RemoveAndAddPoints(const vector<int> & indices_to_remove,
				const vector<EtaPhi> & points_to_add,
				vector<int> & indices_added,
				vector<int> & indices_of_updated_neighbours) {
  
  vector<int> indices1, indices2;
  
  _DNN1->RemoveAndAddPoints(indices_to_remove,points_to_add,
				      indices_added,indices1);

  // create a vector with the remapped points (pi..3pi)
  vector<EtaPhi> remapped_points(points_to_add.size());
  for (size_t i = 0; i < points_to_add.size(); i++) {
    remapped_points[i] = _remap_phi(points_to_add[i]);
  }
  _DNN2->RemoveAndAddPoints(indices_to_remove, remapped_points, 
				      indices_added,indices2);
  
  // merge the two sequences of updated vertices, avoiding double entries
  // of vertices with the same index
  set<int> index_set;
  unsigned int i;
  for (i=0; i < indices1.size(); i++) {index_set.insert(indices1[i]);}
  for (i=0; i < indices2.size(); i++) {index_set.insert(indices2[i]);}

  indices_of_updated_neighbours.clear();
  for (set<int>::iterator iter = index_set.begin(); 
       iter != index_set.end(); iter++) {
    indices_of_updated_neighbours.push_back(*iter);
  }
}

FASTJET_END_NAMESPACE

#endif //  DROP_CGAL 
