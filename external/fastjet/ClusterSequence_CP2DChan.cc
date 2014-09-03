//FJSTARTHEADER
// $Id: ClusterSequence_CP2DChan.cc 3433 2014-07-23 08:17:03Z salam $
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

#include "fastjet/ClusterSequence.hh"
#include "fastjet/internal/ClosestPair2D.hh"
#include<limits>
#include<vector>
#include<cmath>
#include<iostream>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

using namespace std;

// place for things we don't want outside world to run into
namespace Private {
  /// class for helping us deal with mirror-image particles.
  class MirrorInfo{
  public:
    int orig, mirror;
    MirrorInfo(int a, int b) : orig(a), mirror(b) {};
    MirrorInfo() {};
  };

  /// if there is a need for a mirror when looking for closest pairs
  /// up to distance D, then return true and turn the supplied point
  /// into its mirror copy
  bool make_mirror(Coord2D & point, double Dlim) {
    if (point.y < Dlim)       {point.y += twopi; return true;}
    if (twopi-point.y < Dlim) {point.y -= twopi; return true;}
    return false;
  }
  
}

using namespace Private;


//----------------------------------------------------------------------
/// clusters only up to a distance Dlim -- does not deal with "inclusive" jets
/// -- these are left to some other part of the program
void ClusterSequence::_CP2DChan_limited_cluster (double Dlim) {
  
  unsigned int n = _initial_n;

  vector<MirrorInfo>   coordIDs(2*n); // coord IDs of a given jetID
  vector<int>          jetIDs(2*n);   // jet ID for a given coord ID
  vector<Coord2D>      coords(2*n);   // our coordinates (and copies)

  // particles within a distance Dlim of the phi edges (phi<Dlim ||
  // phi>2pi-Dli;) will be mirrored. For Dlim>pi, this could lead to
  // particles copies outside the fixed range in phi which is
  // [-pi:3pi] (see make_mirror above). Since in that case all
  // particles get copied anywaym we can just copy particles up to a
  // distance "pi" from the edges
  double Dlim4mirror = min(Dlim,pi);

  // start things off...
  double minrap = numeric_limits<double>::max();
  double maxrap = -minrap;
  int coord_index = -1;
  int n_active = 0;
  for (unsigned jet_i = 0; jet_i < _jets.size(); jet_i++) {

    // skip jets that already have children or that have infinite
    // rapidity
    if (_history[_jets[jet_i].cluster_hist_index()].child != Invalid ||
	(_jets[jet_i].E() == abs(_jets[jet_i].pz()) && 
	 _jets[jet_i].perp2() == 0.0)
	) {continue;}

    n_active++;

    coordIDs[jet_i].orig = ++coord_index;
    coords[coord_index]  = Coord2D(_jets[jet_i].rap(), _jets[jet_i].phi_02pi());
    jetIDs[coord_index]  = jet_i;
    minrap = min(coords[coord_index].x,minrap);
    maxrap = max(coords[coord_index].x,maxrap);

    Coord2D mirror_point(coords[coord_index]);
    if (make_mirror(mirror_point, Dlim4mirror)) {
      coordIDs[jet_i].mirror = ++coord_index;
      coords[coord_index] = mirror_point;
      jetIDs[coord_index] = jet_i;
    } else {
      coordIDs[jet_i].mirror = Invalid;
    }
  }

  // set them to their actual size...
  coords.resize(coord_index+1);

  // establish limits (with some leeway on rapidity)
  Coord2D left_edge(minrap-1.0, -3.15); // a security margin below  -pi
  Coord2D right_edge(maxrap+1.0, 9.45); // a security margin above 3*pi

  //cerr << "minrap, maxrap = " << minrap << " " << maxrap << endl;

  // now create the closest pair search object
  ClosestPair2D cp(coords, left_edge, right_edge);

  // cross check to see what's going on...
  //cerr << "Tree depths before: ";
  //cp.print_tree_depths(cerr);

  // and start iterating...
  vector<Coord2D> new_points(2);
  vector<unsigned int> cIDs_to_remove(4);
  vector<unsigned int> new_cIDs(2);

  do {
    // find out which pair we recombine
    unsigned int cID1, cID2;
    double distance2;
    cp.closest_pair(cID1,cID2,distance2);

    // if the closest distance > Dlim, we exit and go to "inclusive" stage
    if (distance2 > Dlim*Dlim) {break;}

    // normalise distance as we like it
    distance2 *= _invR2;

    // do the recombination...
    int jet_i = jetIDs[cID1];
    int jet_j = jetIDs[cID2];
    assert (jet_i != jet_j); // to catch issue of recombining with mirror point
    int newjet_k;
    _do_ij_recombination_step(jet_i, jet_j, distance2, newjet_k);

    // don't bother with any further action if only one active particle
    // is left (also avoid closest-pair error [cannot remove last particle]).
    if (--n_active == 1) {break;}

    // now prepare operations on CP structure
    cIDs_to_remove.resize(0);
    cIDs_to_remove.push_back(coordIDs[jet_i].orig);
    cIDs_to_remove.push_back(coordIDs[jet_j].orig);
    if (coordIDs[jet_i].mirror != Invalid) 
      cIDs_to_remove.push_back(coordIDs[jet_i].mirror);
    if (coordIDs[jet_j].mirror != Invalid) 
      cIDs_to_remove.push_back(coordIDs[jet_j].mirror);

    Coord2D new_point(_jets[newjet_k].rap(),_jets[newjet_k].phi_02pi());
    new_points.resize(0);
    new_points.push_back(new_point);
    if (make_mirror(new_point, Dlim4mirror)) new_points.push_back(new_point);  //< same warning as before concerning the mirroring
    
    // carry out actions on search tree
    cp.replace_many(cIDs_to_remove, new_points, new_cIDs);

    // now fill in info for new points...
    coordIDs[newjet_k].orig = new_cIDs[0];
    jetIDs[new_cIDs[0]]       = newjet_k;
    if (new_cIDs.size() == 2) {
      coordIDs[newjet_k].mirror = new_cIDs[1];
      jetIDs[new_cIDs[1]]         = newjet_k;
    } else {coordIDs[newjet_k].mirror = Invalid;}
    
    //// if we've reached one "active" jet we should exit...
    //n_active--;
    //if (n_active == 1) {break;}

  } while(true);
  
  // cross check to see what's going on...
  //cerr << "Max tree depths after: ";
  //cp.print_tree_depths(cerr);

}


//----------------------------------------------------------------------
/// a variant of the closest pair clustering which uses a region of
/// size 2pi+2R in phi.
void ClusterSequence::_CP2DChan_cluster_2pi2R () {

  if (_jet_algorithm != cambridge_algorithm) throw Error("CP2DChan clustering method called for a jet-finder that is not the cambridge algorithm");

  // run the clustering with mirror copies kept such that only things
  // within _Rparam of a border are mirrored
  _CP2DChan_limited_cluster(_Rparam);

  //
  _do_Cambridge_inclusive_jets();
}


//----------------------------------------------------------------------
/// a variant of the closest pair clustering which uses a region of
/// size 2pi + 2*0.3 and then carries on with 2pi+2*R
void ClusterSequence::_CP2DChan_cluster_2piMultD () {

  // do a first run of clustering up to a small distance parameter,
  if (_Rparam >= 0.39) {
    _CP2DChan_limited_cluster(min(_Rparam/2,0.3));
  }

  // and then the final round of clustering
  _CP2DChan_cluster_2pi2R ();
}


//----------------------------------------------------------------------
/// a 4pi variant of the closest pair clustering
void ClusterSequence::_CP2DChan_cluster () {

  if (_jet_algorithm != cambridge_algorithm) throw Error("_CP2DChan_cluster called for a jet-finder that is not the cambridge algorithm");

  unsigned int n = _jets.size();

  vector<MirrorInfo>   coordIDs(2*n);  // link from original to mirror indices
  vector<int>          jetIDs(2*n);     // link from mirror to original indices
  vector<Coord2D>      coords(2*n);   // our coordinates (and copies)

  // start things off...
  double minrap = numeric_limits<double>::max();
  double maxrap = -minrap;
  int coord_index = 0;
  for (unsigned i = 0; i < n; i++) {
    // separate out points with infinite rapidity
    if (_jets[i].E() == abs(_jets[i].pz()) && _jets[i].perp2() == 0.0) {
      coordIDs[i] = MirrorInfo(BeamJet,BeamJet);
    } else {
      coordIDs[i].orig   = coord_index;
      coordIDs[i].mirror = coord_index+1;
      coords[coord_index]   = Coord2D(_jets[i].rap(), _jets[i].phi_02pi());
      coords[coord_index+1] = Coord2D(_jets[i].rap(), _jets[i].phi_02pi()+twopi);
      jetIDs[coord_index]   = i;
      jetIDs[coord_index+1] = i;
      minrap = min(coords[coord_index].x,minrap);
      maxrap = max(coords[coord_index].x,maxrap);
      coord_index += 2;
    }
  }
  // label remaining "mirror info" as saying that there are no jets
  for (unsigned i = n; i < 2*n; i++) {coordIDs[i].orig = Invalid;}

  // set them to their actual size...
  coords.resize(coord_index);

  // establish limits (with some leeway on rapidity)
  Coord2D left_edge(minrap-1.0, 0.0);
  Coord2D right_edge(maxrap+1.0, 2*twopi);

  // now create the closest pair search object
  ClosestPair2D cp(coords, left_edge, right_edge);

  // and start iterating...
  vector<Coord2D> new_points(2);
  vector<unsigned int> cIDs_to_remove(4);
  vector<unsigned int> new_cIDs(2);
  do {
    // find out which pair we recombine
    unsigned int cID1, cID2;
    double distance2;
    cp.closest_pair(cID1,cID2,distance2);
    distance2 *= _invR2;

    // if the closest distance > 1, we exit and go to "inclusive" stage
    if (distance2 > 1.0) {break;}

    // do the recombination...
    int jet_i = jetIDs[cID1];
    int jet_j = jetIDs[cID2];
    assert (jet_i != jet_j); // to catch issue of recombining with mirror point
    int newjet_k;
    _do_ij_recombination_step(jet_i, jet_j, distance2, newjet_k);

    // now prepare operations on CP structure
    cIDs_to_remove[0] = coordIDs[jet_i].orig;
    cIDs_to_remove[1] = coordIDs[jet_i].mirror;
    cIDs_to_remove[2] = coordIDs[jet_j].orig;
    cIDs_to_remove[3] = coordIDs[jet_j].mirror;
    new_points[0] = Coord2D(_jets[newjet_k].rap(),_jets[newjet_k].phi_02pi());
    new_points[1] = Coord2D(_jets[newjet_k].rap(),_jets[newjet_k].phi_02pi()+twopi);
    // carry out the CP operation
    //cp.replace_many(cIDs_to_remove, new_points, new_cIDs);
    // remarkable the following is faster...
    new_cIDs[0] = cp.replace(cIDs_to_remove[0], cIDs_to_remove[2], new_points[0]);
    new_cIDs[1] = cp.replace(cIDs_to_remove[1], cIDs_to_remove[3], new_points[1]);
    // signal that the following jets no longer exist
    coordIDs[jet_i].orig = Invalid;
    coordIDs[jet_j].orig = Invalid;
    // and do bookkeeping for new jet
    coordIDs[newjet_k] = MirrorInfo(new_cIDs[0], new_cIDs[1]);
    jetIDs[new_cIDs[0]] = newjet_k;
    jetIDs[new_cIDs[1]] = newjet_k;

    // if we've reached one jet we should exit...
    n--;
    if (n == 1) {break;}

  } while(true);
  

  // now do "beam" recombinations 
  //for (unsigned int jet_i = 0; jet_i < coordIDs.size(); jet_i++) {
  //  // if we have a predeclared beam jet or a valid set of mirror
  //  // coordinates, recombine then recombine this jet with the beam
  //  if (coordIDs[jet_i].orig == BeamJet || coordIDs[jet_i].orig > 0) {
  //    // diB = 1 _always_ (for the cambridge alg.)
  //    _do_iB_recombination_step(jet_i, 1.0);
  //  }
  //}

  _do_Cambridge_inclusive_jets();

  //n = _history.size();
  //for (unsigned int hist_i = 0; hist_i < n; hist_i++) {
  //  if (_history[hist_i].child == Invalid) {
  //    _do_iB_recombination_step(_history[hist_i].jetp_index, 1.0);
  //  }
  //}

}


//----------------------------------------------------------------------
void ClusterSequence::_do_Cambridge_inclusive_jets () {
  unsigned int n = _history.size();
  for (unsigned int hist_i = 0; hist_i < n; hist_i++) {
    if (_history[hist_i].child == Invalid) {
      _do_iB_recombination_step(_history[hist_i].jetp_index, 1.0);
    }
  }
}

FASTJET_END_NAMESPACE

