#ifndef __FASTJET_LAZYTILING9ALT_HH__
#define __FASTJET_LAZYTILING9ALT_HH__

//FJSTARTHEADER
// $Id: LazyTiling9Alt.hh 3808 2015-02-20 11:24:53Z soyez $
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

//#include "fastjet/PseudoJet.hh"
#include "fastjet/internal/MinHeap.hh"
#include "fastjet/ClusterSequence.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

/// Rounding errors in the Lazy strategies may cause the following
/// problem: when browsing tiles in the vicinity of the particles
/// being clustered in order to decide which of these tiles may
/// contain particles that need to be updated (because theit NN is one
/// of the particles that are currently clustered), we discard tiles
/// that are deemed "too far from the cell" by the "max_NN_dist"
/// criterion. Because of rounding error, this condition can sometimes
/// miss cases where an update is needed.
///
/// An example of this happens if a particle '1' is, say, at the lower
/// edge of the rapidity of a given tile, with a particle '2' in the
/// tile directly on its left at the same rapidity. Assume also that
/// max_NN_dist in 2's tile corresponds to the distance between 2 and
/// teh tile of 1. If 2 is 1's NN then in case 2 gets clustered, 1's
/// NN needs to be updated. However, rounding errors in the
/// calculation of the distance between 1 and 2 may result is
/// something slightly larger than the max_NN_dist in 2's tile.
///
/// This situation corresponds to the bug reported by Jochen Olt on
/// February 12 2015 [see issue-tracker/2015-02-infinite-loop],
/// causing an infinite loop.
///
/// To prevent this, the simplest solution is, when looking at tiles
/// to browse for updateds, to add a margin of security close to the
/// edges of the cell, i.e. instead of updating only tiles for which
/// distance<=max_NN_dist, we will update tiles for which
/// distance<=max_NN_dist+tile_edge_security_margin.
///
/// Note that this does not need to be done when computing nearest
/// neighbours [rounding errors are tolerated there] but it is
/// critical when tracking points that have to be updated.
const double tile_edge_security_margin=1.0e-7;

/// structure analogous to BriefJet, but with the extra information
/// needed for dealing with tiles
class TiledJet {
public:
  double     eta, phi, kt2, NN_dist;
  TiledJet * NN, *previous, * next; 
  int        _jets_index, tile_index;
  bool _minheap_update_needed;

  // indicate whether jets need to have their minheap entries
  // updated).
  inline void label_minheap_update_needed() {_minheap_update_needed = true;}
  inline void label_minheap_update_done()   {_minheap_update_needed = false;}
  inline bool minheap_update_needed() const {return _minheap_update_needed;}
};

const int n_tile_neighbours = 9;

class Tile {
public:
  typedef double (Tile::*DistToTileFn)(const TiledJet*) const;
  typedef std::pair<Tile *, DistToTileFn> TileFnPair;
  /// pointers to neighbouring tiles, including self
  TileFnPair begin_tiles[n_tile_neighbours]; 
  /// neighbouring tiles, excluding self
  TileFnPair *  surrounding_tiles; 
  /// half of neighbouring tiles, no self
  TileFnPair *  RH_tiles;  
  /// just beyond end of tiles
  TileFnPair *  end_tiles; 
  /// start of list of BriefJets contained in this tile
  TiledJet * head;    
  /// sometimes useful to be able to tag a tile
  bool     tagged;    
  /// true for tiles where the delta phi calculation needs
  /// potentially to account for periodicity in phi
  bool     use_periodic_delta_phi;
  /// for all particles in the tile, this stores the largest of the
  /// (squared) nearest-neighbour distances.
  double max_NN_dist;
  double eta_min, eta_max, phi_min, phi_max;

  double distance_to_centre(const TiledJet *) const {return 0;}
  double distance_to_left(const TiledJet * jet) const {
    double deta = jet->eta - eta_min;
    return deta*deta;
  }
  double distance_to_right(const TiledJet * jet) const {
    double deta = jet->eta - eta_max;
    return deta*deta;
  }
  double distance_to_bottom(const TiledJet * jet) const {
    double dphi = jet->phi - phi_min;
    return dphi*dphi;
  }
  double distance_to_top(const TiledJet * jet) const {
    double dphi = jet->phi - phi_max;
    return dphi*dphi;
  }

  double distance_to_left_top(const TiledJet * jet) const {
    double deta = jet->eta - eta_min;
    double dphi = jet->phi - phi_max;
    return deta*deta + dphi*dphi;
  }
  double distance_to_left_bottom(const TiledJet * jet) const {
    double deta = jet->eta - eta_min;
    double dphi = jet->phi - phi_min;
    return deta*deta + dphi*dphi;
  }
  double distance_to_right_top(const TiledJet * jet) const {
    double deta = jet->eta - eta_max;
    double dphi = jet->phi - phi_max;
    return deta*deta + dphi*dphi;
  }
  double distance_to_right_bottom(const TiledJet * jet) const {
    double deta = jet->eta - eta_max;
    double dphi = jet->phi - phi_min;
    return deta*deta + dphi*dphi;
  }

  
};

//----------------------------------------------------------------------
class LazyTiling9Alt {
public:
  LazyTiling9Alt(ClusterSequence & cs);

  void run();

  //void get_next_clustering(int & jetA_index, int & jetB_index, double & dij);
  

protected:
  ClusterSequence & _cs;
  const std::vector<PseudoJet> & _jets;
  std::vector<Tile> _tiles;


  double _Rparam, _R2, _invR2;
  double _tiles_eta_min, _tiles_eta_max;
  double _tile_size_eta, _tile_size_phi;
  double _tile_half_size_eta, _tile_half_size_phi;
  int    _n_tiles_phi,_tiles_ieta_min,_tiles_ieta_max;

  std::vector<TiledJet *> _jets_for_minheap;
  
  //MinHeap _minheap;

  void _initialise_tiles();

  // reasonably robust return of tile index given ieta and iphi, in particular
  // it works even if iphi is negative
  inline int _tile_index (int ieta, int iphi) const {
    // note that (-1)%n = -1 so that we have to add _n_tiles_phi
    // before performing modulo operation
    return (ieta-_tiles_ieta_min)*_n_tiles_phi
                  + (iphi+_n_tiles_phi) % _n_tiles_phi;
  }

  void  _bj_remove_from_tiles(TiledJet * const jet);

  /// returns the tile index given the eta and phi values of a jet
  int _tile_index(const double eta, const double phi) const;

  // sets up information regarding the tiling of the given jet
  void _tj_set_jetinfo(TiledJet * const jet, const int _jets_index);

  void _print_tiles(TiledJet * briefjets ) const;
  void _add_neighbours_to_tile_union(const int tile_index, 
		 std::vector<int> & tile_union, int & n_near_tiles) const;
  void _add_untagged_neighbours_to_tile_union(const int tile_index, 
		 std::vector<int> & tile_union, int & n_near_tiles);
  void _add_untagged_neighbours_to_tile_union_using_max_info(const TiledJet * const jet, 
		 std::vector<int> & tile_union, int & n_near_tiles);
  //double _distance_to_tile(const TiledJet * bj, const Tile *) const;
  void _update_jetX_jetI_NN(TiledJet * jetX, TiledJet * jetI, std::vector<TiledJet *> & jets_for_minheap);

  void _set_NN(TiledJet * jetI, std::vector<TiledJet *> & jets_for_minheap);

  // return the diJ (multiplied by _R2) for this jet assuming its NN
  // info is correct
  template <class J> double _bj_diJ(const J * const jet) const {
    double kt2 = jet->kt2;
    if (jet->NN != NULL) {if (jet->NN->kt2 < kt2) {kt2 = jet->NN->kt2;}}
    return jet->NN_dist * kt2;
  }


  //----------------------------------------------------------------------
  template <class J> inline void _bj_set_jetinfo(
                            J * const jetA, const int _jets_index) const {
    jetA->eta  = _jets[_jets_index].rap();
    jetA->phi  = _jets[_jets_index].phi_02pi();
    jetA->kt2  = _cs.jet_scale_for_algorithm(_jets[_jets_index]);
    jetA->_jets_index = _jets_index;
    // initialise NN info as well
    jetA->NN_dist = _R2;
    jetA->NN      = NULL;
  }


  //----------------------------------------------------------------------
  template <class J> inline double _bj_dist(
                const J * const jetA, const J * const jetB) const {
    double dphi = std::abs(jetA->phi - jetB->phi);
    double deta = (jetA->eta - jetB->eta);
    if (dphi > pi) {dphi = twopi - dphi;}
    return dphi*dphi + deta*deta;
  }


  //----------------------------------------------------------------------
  template <class J> inline double _bj_dist_not_periodic(
                const J * const jetA, const J * const jetB) const {
    double dphi = jetA->phi - jetB->phi;
    double deta = (jetA->eta - jetB->eta);
    return dphi*dphi + deta*deta;
  }

};


FASTJET_END_NAMESPACE

#endif // __FASTJET_LAZYTILING9ALT_HH__
