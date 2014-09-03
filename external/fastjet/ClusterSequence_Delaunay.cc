//FJSTARTHEADER
// $Id: ClusterSequence_Delaunay.cc 3475 2014-07-29 11:57:23Z salam $
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


#include "fastjet/Error.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/internal/DynamicNearestNeighbours.hh"
#include<iostream>
#include<sstream>
#include<cmath>
#include <cstdlib>
#include<cassert>
#include<memory>
//
#ifndef DROP_CGAL // in case we do not have the code for CGAL
#include "fastjet/internal/Dnn4piCylinder.hh"
#include "fastjet/internal/Dnn3piCylinder.hh"
#include "fastjet/internal/Dnn2piCylinder.hh"
#endif //  DROP_CGAL 

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

using namespace std;


//----------------------------------------------------------------------
/// Run the clustering using a Hierarchical Delaunay triangulation and
/// STL maps to achieve O(N*ln N) behaviour.
///
/// There may be internally asserted assumptions about absence of
/// points with coincident eta-phi coordinates.
void ClusterSequence::_delaunay_cluster () {

  int n = _jets.size();

  vector<EtaPhi> points(n); // recall EtaPhi is just a typedef'd pair<double>
  for (int i = 0; i < n; i++) {
    points[i] = EtaPhi(_jets[i].rap(),_jets[i].phi_02pi());
    points[i].sanitize(); // make sure things are in the right range
  }

  // initialise our DNN structure with the set of points
  auto_ptr<DynamicNearestNeighbours> DNN;
#ifndef DROP_CGAL // strategy = NlnN* are not supported if we drop CGAL...
  bool verbose = false;
  bool ignore_nearest_is_mirror = (_Rparam < twopi);
  if (_strategy == NlnN4pi) {
    DNN.reset(new Dnn4piCylinder(points,verbose));
  } else if (_strategy == NlnN3pi) {
    DNN.reset(new Dnn3piCylinder(points,ignore_nearest_is_mirror,verbose));
  } else if (_strategy == NlnN) {
    DNN.reset(new Dnn2piCylinder(points,ignore_nearest_is_mirror,verbose));
  } else 
#else
  if (_strategy == NlnN4pi || _strategy == NlnN3pi || _strategy == NlnN) {
    ostringstream err;
    err << "ERROR: Requested strategy "<<strategy_string()<<" but it is not"<<endl;
    err << "       supported because FastJet was compiled without CGAL"<<endl;
    throw Error(err.str());
    //assert(false);
  } else
#endif // DROP_CGAL
  {
    //ostringstream err;
    //err << "ERROR: Unrecognized value for strategy: "<<_strategy<<endl;
    //throw Error(err.str());
    //-----------------------------------------------------------------
    // The code should never reach this point, because the checks above
    // should always handle all _strategy values for which 
    // _delaunay_cluster() is called 
    assert(false);
  }

  // We will find nearest neighbour for each vertex, and include
  // distance in map (NB DistMap is a typedef given in the .h file)
  DistMap DijMap;

  // fill the map with the minimal (as far as we know) subset of Dij
  // distances (i.e. nearest neighbour ones).
  for (int ii = 0; ii < n; ii++) {
    _add_ktdistance_to_map(ii, DijMap, DNN.get());
  }

  // run the clustering (go up to i=n-1, but then will stop half-way down,
  // when we reach that point -- it will be the final beam jet and there
  // will be no nearest neighbours to find).
  for (int i=0;i<n;i++) {
    // find nearest vertices
    // NB: skip cases where the point is not there anymore!
    TwoVertices SmallestDijPair;
    int jet_i, jet_j;
    double SmallestDij;
    bool Valid2;
    bool recombine_with_beam;
    do { 
      SmallestDij = DijMap.begin()->first;
      SmallestDijPair = DijMap.begin()->second;
      jet_i = SmallestDijPair.first;
      jet_j = SmallestDijPair.second;
      // distance is immediately removed regardless of whether or not
      // it is used.
      // Some temporary testing code relating to problems with the gcc-3.2 compiler
      //cout << "got here and size is "<< DijMap.size()<< " and it is "<<SmallestDij <<"\n";
      //cout <<  jet_i << " "<< jet_j<<"\n";
      DijMap.erase(DijMap.begin());
      //cout << "got beyond here\n";

      // need to "prime" the validity of jet_j in such a way that 
      // if it corresponds to the beam then it is automatically valid.
      recombine_with_beam = (jet_j == BeamJet);
      if (!recombine_with_beam) {Valid2 = DNN->Valid(jet_j);} 
      else {Valid2 = true;}

    } while ( !DNN->Valid(jet_i) || !Valid2);


    // The following part acts just on jet momenta and on the history.
    // The action on the nearest-neighbour structures takes place
    // later (only if at least 2 jets are around).
    if (! recombine_with_beam) {
      int nn; // will be index of new jet
      _do_ij_recombination_step(jet_i, jet_j, SmallestDij, nn);
      //OBS // merge the two jets, add new jet, remove old ones
      //OBS _jets.push_back(_jets[jet_i] + _jets[jet_j]);
      //OBS 
      //OBS int nn = _jets.size()-1;
      //OBS _jets[nn].set_cluster_hist_index(n+i);
      //OBS 
      //OBS // get corresponding indices in history structure
      //OBS int hist_i = _jets[jet_i].cluster_hist_index();
      //OBS int hist_j = _jets[jet_j].cluster_hist_index();
      //OBS 
      //OBS 
      //OBS _add_step_to_history(n+i,min(hist_i,hist_j), max(hist_i,hist_j),
      //OBS 		      _jets.size()-1, SmallestDij);

      // add new point to points vector
      EtaPhi newpoint(_jets[nn].rap(), _jets[nn].phi_02pi());
      newpoint.sanitize(); // make sure it is in correct range
      points.push_back(newpoint);
    } else {
      // recombine the jet with the beam
      _do_iB_recombination_step(jet_i, SmallestDij);
      //OBS _add_step_to_history(n+i,_jets[jet_i].cluster_hist_index(),BeamJet,
      //OBS 			   Invalid, SmallestDij);
    }

    // exit the loop because we do not want to look for nearest neighbours
    // etc. of zero partons
    if (i == n-1) {break;}

    vector<int> updated_neighbours;
    if (! recombine_with_beam) {
      // update DNN
      int point3;
      DNN->RemoveCombinedAddCombination(jet_i, jet_j, 
				       points[points.size()-1], point3,
				       updated_neighbours);
      // C++ beginners' comment: static_cast to unsigned int is necessary
      // to do away with warnings about type mismatch between point3 (int) 
      // and points.size (unsigned int)
      if (static_cast<unsigned int> (point3) != points.size()-1) {
	throw Error("INTERNAL ERROR: point3 != points.size()-1");}
    } else {
      // update DNN
      DNN->RemovePoint(jet_i, updated_neighbours);
    }

    // update map
    vector<int>::iterator it = updated_neighbours.begin();
    for (; it != updated_neighbours.end(); ++it) {
      int ii = *it;
      _add_ktdistance_to_map(ii, DijMap, DNN.get());
    }
      
  } // end clustering loop 
  
}


//----------------------------------------------------------------------
/// Add the current kt distance for particle ii to the map (DijMap)
/// using information from the DNN object. Work as follows:
/// 
/// . if the kt is zero then it's nearest neighbour is taken to be the
///   the beam jet and the distance is zero.
///
/// . if cylinder distance to nearest neighbour > _Rparam then it is
///   yiB that is smallest and this is added to map.
///
/// . otherwise if the nearest neighbour jj has a larger kt then add
///   dij to the map.
///
/// . otherwise do nothing
///
void ClusterSequence::_add_ktdistance_to_map(
                          const int ii, 
			  DistMap & DijMap,
			  const DynamicNearestNeighbours * DNN) {
  
  double yiB = jet_scale_for_algorithm(_jets[ii]);
  if (yiB == 0.0) {
    // in this case convention is that we do not worry about distances
    // but directly state that nearest neighbour is beam
    DijMap.insert(DijEntry(yiB,  TwoVertices(ii,-1)));
  } else {
    double DeltaR2 = DNN->NearestNeighbourDistance(ii) * _invR2;
    // Logic of following bit is: only add point to map if it has
    // smaller kt2 than nearest neighbour j (if it has larger kt,
    // then: either it is j's nearest neighbour and then we will
    // include dij when we come to j; or it is not j's nearest
    // neighbour and j will recombine with someone else).
    
    // If DeltaR2 > 1.0 then in any case it will recombine with beam rather
    // than with any neighbours.
    // (put general normalisation here at some point)
    if (DeltaR2 > 1.0) {
      DijMap.insert(DijEntry(yiB,  TwoVertices(ii,-1)));
    } else {
      double kt2i = jet_scale_for_algorithm(_jets[ii]);
      int jj = DNN->NearestNeighbourIndex(ii);
      if (kt2i <= jet_scale_for_algorithm(_jets[jj])) {
	double dij = DeltaR2 * kt2i;
	DijMap.insert(DijEntry(dij, TwoVertices(ii,jj)));
      }
    }
  }
}


FASTJET_END_NAMESPACE

