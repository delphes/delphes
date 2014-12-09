#ifndef __FASTJET_LAZYTILING9_HH__
#define __FASTJET_LAZYTILING9_HH__

//#define INSTRUMENT2 1

//FJSTARTHEADER
// $Id: LazyTiling9.hh 3477 2014-07-29 14:34:39Z salam $
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
#include "fastjet/internal/LazyTiling9Alt.hh"


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh
template<int NN>
class Tile2Base {
public:
  /// pointers to neighbouring tiles, including self
  Tile2Base *   begin_tiles[NN]; 
  /// neighbouring tiles, excluding self
  Tile2Base **  surrounding_tiles; 
  /// half of neighbouring tiles, no self
  Tile2Base **  RH_tiles;  
  /// just beyond end of tiles
  Tile2Base **  end_tiles; 
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
  double eta_centre, phi_centre;

  /// returns the number of jets in the tile; useful principally for
  /// diagnostics
  int jet_count() const {
    int count = 0;
    const TiledJet * jet = head;
    while (jet != 0) {
      count++;
      jet = jet->next;
    }
    return count;
  }
};


typedef Tile2Base<9> Tile2;

// class Tile2 : public Tile2Base {
// public:
//   /// pointers to neighbouring tiles, including self
//   Tile2 *   begin_tiles[n_tile_neighbours]; 
//   bool is_near_zero_phi(double tile_size_phi) const {
//     return phi_centre < tile_size_phi || (twopi-phi_centre) < tile_size_phi;
//   }
// };

// class Tile2  {
// public:
//   /// pointers to neighbouring tiles, including self
//   Tile2 *   begin_tiles[n_tile_neighbours]; 
//   /// neighbouring tiles, excluding self
//   Tile2 **  surrounding_tiles; 
//   /// half of neighbouring tiles, no self
//   Tile2 **  RH_tiles;  
//   /// just beyond end of tiles
//   Tile2 **  end_tiles; 
//   /// start of list of BriefJets contained in this tile
//   TiledJet * head;    
//   /// sometimes useful to be able to tag a tile
//   bool     tagged;    
//   /// for all particles in the tile, this stores the largest of the
//   /// (squared) nearest-neighbour distances.
//   double max_NN_dist;
//   double eta_centre, phi_centre;
//   bool is_near_zero_phi(double tile_size_phi) const {
//     return phi_centre < tile_size_phi || (twopi-phi_centre) < tile_size_phi;
//   }
// };



//----------------------------------------------------------------------
class LazyTiling9 {
public:
  LazyTiling9(ClusterSequence & cs);

  void run();

  //void get_next_clustering(int & jetA_index, int & jetB_index, double & dij);
  

protected:
  ClusterSequence & _cs;
  const std::vector<PseudoJet> & _jets;
  std::vector<Tile2> _tiles;

#ifdef INSTRUMENT2
  int _ncall; // GPS tmp
  int _ncall_dtt; // GPS tmp
#endif // INSTRUMENT2

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
  double _distance_to_tile(const TiledJet * bj, const Tile2 *) 
#ifdef INSTRUMENT2
    ;
#else
    const;
#endif 
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
                const J * const jetA, const J * const jetB) 
#ifdef INSTRUMENT2
    {
    _ncall++; // GPS tmp
#else
    const {
#endif 
    double dphi = std::abs(jetA->phi - jetB->phi);
    double deta = (jetA->eta - jetB->eta);
    if (dphi > pi) {dphi = twopi - dphi;}
    return dphi*dphi + deta*deta;
  }


  //----------------------------------------------------------------------
  template <class J> inline double _bj_dist_not_periodic(
                const J * const jetA, const J * const jetB)
#ifdef INSTRUMENT2
    {
    _ncall++; // GPS tmp
#else
    const {
#endif 
    //_ncall++; // GPS tmp
    double dphi = jetA->phi - jetB->phi;
    double deta = (jetA->eta - jetB->eta);
    return dphi*dphi + deta*deta;
  }

};


FASTJET_END_NAMESPACE

#endif // __FASTJET_LAZYTILING9_HH__
