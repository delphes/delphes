//FJSTARTHEADER
// $Id: LazyTiling9Alt.cc 3477 2014-07-29 14:34:39Z salam $
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

#include <iomanip>
#include "fastjet/internal/LazyTiling9Alt.hh"
#include "fastjet/internal/TilingExtent.hh"
using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

LazyTiling9Alt::LazyTiling9Alt(ClusterSequence & cs) :
  _cs(cs), _jets(cs.jets())
  //, _minheap(_jets.size()) 
{
  _Rparam = cs.jet_def().R();
  _R2 = _Rparam * _Rparam;
  _invR2 = 1.0 / _R2;
  _initialise_tiles();
}


//----------------------------------------------------------------------
/// Set up the tiles:
///  - decide the range in eta
///  - allocate the tiles
///  - set up the cross-referencing info between tiles
///
/// The neighbourhood of a tile is set up as follows
///
/// 	      LRR
///           LXR
///           LLR
///
/// such that tiles is an array containing XLLLLRRRR with pointers
///                                         |   \ RH_tiles
///                                         \ surrounding_tiles
///
/// with appropriate precautions when close to the edge of the tiled
/// region.
///
void LazyTiling9Alt::_initialise_tiles() {

  // first decide tile sizes (with a lower bound to avoid huge memory use with
  // very small R)
  double default_size = max(0.1,_Rparam);
  _tile_size_eta = default_size;
  // it makes no sense to go below 3 tiles in phi -- 3 tiles is
  // sufficient to make sure all pair-wise combinations up to pi in
  // phi are possible
  _n_tiles_phi   = max(3,int(floor(twopi/default_size)));
  _tile_size_phi = twopi / _n_tiles_phi; // >= _Rparam and fits in 2pi

  // always include zero rapidity in the tiling region
  _tiles_eta_min = 0.0;
  _tiles_eta_max = 0.0;
  // but go no further than following
  const double maxrap = 7.0;

  // and find out how much further one should go
  for(unsigned int i = 0; i < _jets.size(); i++) {
    double eta = _jets[i].rap();
    // first check if eta is in range -- to avoid taking into account
    // very spurious rapidities due to particles with near-zero kt.
    if (abs(eta) < maxrap) {
      if (eta < _tiles_eta_min) {_tiles_eta_min = eta;}
      if (eta > _tiles_eta_max) {_tiles_eta_max = eta;}
    }
  }

  // now adjust the values
  _tiles_ieta_min = int(floor(_tiles_eta_min/_tile_size_eta));
  _tiles_ieta_max = int(floor( _tiles_eta_max/_tile_size_eta));
  _tiles_eta_min = _tiles_ieta_min * _tile_size_eta;
  _tiles_eta_max = _tiles_ieta_max * _tile_size_eta;

  _tile_half_size_eta = _tile_size_eta * 0.5;
  _tile_half_size_phi = _tile_size_phi * 0.5;

  // set up information about whether we need to allow for "periodic" 
  // wrapping tests in delta_phi calculations
  vector<bool> use_periodic_delta_phi(_n_tiles_phi, false);
  if (_n_tiles_phi <= 3) {
    fill(use_periodic_delta_phi.begin(), use_periodic_delta_phi.end(), true);
  } else {
    use_periodic_delta_phi[0] = true;
    use_periodic_delta_phi[_n_tiles_phi-1] = true;
  }

  // allocate the tiles
  _tiles.resize((_tiles_ieta_max-_tiles_ieta_min+1)*_n_tiles_phi);

  // now set up the cross-referencing between tiles
  for (int ieta = _tiles_ieta_min; ieta <= _tiles_ieta_max; ieta++) {
    for (int iphi = 0; iphi < _n_tiles_phi; iphi++) {
      Tile * tile = & _tiles[_tile_index(ieta,iphi)];
      // no jets in this tile yet
      tile->head = NULL; // first element of tiles points to itself
      tile->begin_tiles[0] =  Tile::TileFnPair(tile,&Tile::distance_to_centre);
      Tile::TileFnPair * pptile = & (tile->begin_tiles[0]);
      pptile++;
      //
      // set up L's in column to the left of X
      tile->surrounding_tiles = pptile;
      if (ieta > _tiles_ieta_min) {
	// with the itile subroutine, we can safely run tiles from
	// idphi=-1 to idphi=+1, because it takes care of
	// negative and positive boundaries
	//for (int idphi = -1; idphi <=+1; idphi++) {
        //  *pptile = & _tiles[_tile_index(ieta-1,iphi+idphi)];
        //  pptile++;
        //}	
        *pptile = Tile::TileFnPair(& _tiles[_tile_index(ieta-1,iphi-1)],
                                   &Tile::distance_to_left_bottom);
        pptile++;
        *pptile = Tile::TileFnPair(& _tiles[_tile_index(ieta-1,iphi)],
                                   &Tile::distance_to_left);
        pptile++;
        *pptile = Tile::TileFnPair(& _tiles[_tile_index(ieta-1,iphi+1)],
                                   &Tile::distance_to_left_top);
        pptile++;
      }
      // now set up last L (below X)
      *pptile = Tile::TileFnPair(& _tiles[_tile_index(ieta,iphi-1)], 
                                 &Tile::distance_to_bottom);
      pptile++;
      // set up first R (above X)
      tile->RH_tiles = pptile;
      *pptile = Tile::TileFnPair(& _tiles[_tile_index(ieta,iphi+1)], 
                                 &Tile::distance_to_top);
      pptile++;
      // set up remaining R's, to the right of X
      if (ieta < _tiles_ieta_max) {
	//for (int idphi = -1; idphi <= +1; idphi++) {
	//  *pptile = & _tiles[_tile_index(ieta+1,iphi+idphi)];
	//  pptile++;
	//}	
        *pptile = Tile::TileFnPair(& _tiles[_tile_index(ieta+1,iphi-1)],
                                   &Tile::distance_to_right_bottom);
        pptile++;
        *pptile = Tile::TileFnPair(& _tiles[_tile_index(ieta+1,iphi)],
                                   &Tile::distance_to_right);
        pptile++;
        *pptile = Tile::TileFnPair(& _tiles[_tile_index(ieta+1,iphi+1)],
                                   &Tile::distance_to_right_top);
        pptile++;
      }
      // now put semaphore for end tile
      tile->end_tiles = pptile;
      // finally make sure tiles are untagged
      tile->tagged = false;
      // and store the information about periodicity in phi
      tile->use_periodic_delta_phi = use_periodic_delta_phi[iphi];
      // and ensure max distance is sensibly initialised
      tile->max_NN_dist = 0;
      // and also position of centre of tile
      tile->eta_min = ieta*_tile_size_eta;
      tile->eta_max = (ieta+1)*_tile_size_eta;
      tile->phi_min = iphi*_tile_size_phi;
      tile->phi_max = (iphi+1)*_tile_size_phi;
    }
  }

}

//----------------------------------------------------------------------
/// return the tile index corresponding to the given eta,phi point
int LazyTiling9Alt::_tile_index(const double eta, const double phi) const {
  int ieta, iphi;
  if      (eta <= _tiles_eta_min) {ieta = 0;}
  else if (eta >= _tiles_eta_max) {ieta = _tiles_ieta_max-_tiles_ieta_min;}
  else {
    //ieta = int(floor((eta - _tiles_eta_min) / _tile_size_eta));
    ieta = int(((eta - _tiles_eta_min) / _tile_size_eta));
    // following needed in case of rare but nasty rounding errors
    if (ieta > _tiles_ieta_max-_tiles_ieta_min) {
      ieta = _tiles_ieta_max-_tiles_ieta_min;} 
  }
  // allow for some extent of being beyond range in calculation of phi
  // as well
  //iphi = (int(floor(phi/_tile_size_phi)) + _n_tiles_phi) % _n_tiles_phi;
  // with just int and no floor, things run faster but beware
  iphi = int((phi+twopi)/_tile_size_phi) % _n_tiles_phi;
  return (iphi + ieta * _n_tiles_phi);
}


//----------------------------------------------------------------------
// sets up information regarding the tiling of the given jet
inline void LazyTiling9Alt::_tj_set_jetinfo( TiledJet * const jet,
					      const int _jets_index) {
  // first call the generic setup
  _bj_set_jetinfo<>(jet, _jets_index);

  // Then do the setup specific to the tiled case.

  // Find out which tile it belonds to
  jet->tile_index = _tile_index(jet->eta, jet->phi);

  // Insert it into the tile's linked list of jets
  Tile * tile = &_tiles[jet->tile_index];
  jet->previous   = NULL;
  jet->next       = tile->head;
  if (jet->next != NULL) {jet->next->previous = jet;}
  tile->head      = jet;
}


//----------------------------------------------------------------------
void LazyTiling9Alt::_bj_remove_from_tiles(TiledJet * const jet) {
  Tile * tile = & _tiles[jet->tile_index];

  if (jet->previous == NULL) {
    // we are at head of the tile, so reset it.
    // If this was the only jet on the tile then tile->head will now be NULL
    tile->head = jet->next;
  } else {
    // adjust link from previous jet in this tile
    jet->previous->next = jet->next;
  }
  if (jet->next != NULL) {
    // adjust backwards-link from next jet in this tile
    jet->next->previous = jet->previous;
  }
}


//----------------------------------------------------------------------
/// output the contents of the tiles
void LazyTiling9Alt::_print_tiles(TiledJet * briefjets ) const {
  for (vector<Tile>::const_iterator tile = _tiles.begin(); 
       tile < _tiles.end(); tile++) {
    cout << "Tile " << tile - _tiles.begin()<<" = ";
    vector<int> list;
    for (TiledJet * jetI = tile->head; jetI != NULL; jetI = jetI->next) {
      list.push_back(jetI-briefjets);
      //cout <<" "<<jetI-briefjets;
    }
    sort(list.begin(),list.end());
    for (unsigned int i = 0; i < list.size(); i++) {cout <<" "<<list[i];}
    cout <<"\n";
  }
}


//----------------------------------------------------------------------
/// Add to the vector tile_union the tiles that are in the neighbourhood
/// of the specified tile_index, including itself -- start adding
/// from position n_near_tiles-1, and increase n_near_tiles as
/// you go along (could have done it more C++ like with vector with reserved
/// space, but fear is that it would have been slower, e.g. checking
/// for end of vector at each stage to decide whether to resize it)
void LazyTiling9Alt::_add_neighbours_to_tile_union(const int tile_index, 
	       vector<int> & tile_union, int & n_near_tiles) const {
  for (Tile::TileFnPair const * near_tile = _tiles[tile_index].begin_tiles; 
       near_tile != _tiles[tile_index].end_tiles; near_tile++){
    // get the tile number
    tile_union[n_near_tiles] = near_tile->first - & _tiles[0];
    n_near_tiles++;
  }
}


//----------------------------------------------------------------------
/// Like _add_neighbours_to_tile_union, but only adds neighbours if 
/// their "tagged" status is false; when a neighbour is added its
/// tagged status is set to true.
inline void LazyTiling9Alt::_add_untagged_neighbours_to_tile_union(
               const int tile_index, 
	       vector<int> & tile_union, int & n_near_tiles)  {
  for (Tile::TileFnPair * near_tile = _tiles[tile_index].begin_tiles; 
       near_tile != _tiles[tile_index].end_tiles; near_tile++){
    if (! (near_tile->first)->tagged) {
      (near_tile->first)->tagged = true;
      // get the tile number
      tile_union[n_near_tiles] = near_tile->first - & _tiles[0];
      n_near_tiles++;
    }
  }
}

//----------------------------------------------------------------------
/// Like _add_neighbours_to_tile_union, but adds tiles that are
/// "neighbours" of a jet (rather than a tile) and only if a
/// neighbouring tile's max_NN_dist is >= the distance between the jet
/// and the nearest point on the tile. It ignores tiles that have
/// already been tagged.
inline void LazyTiling9Alt::_add_untagged_neighbours_to_tile_union_using_max_info(
               const TiledJet * jet, 
	       vector<int> & tile_union, int & n_near_tiles)  {
  Tile & tile = _tiles[jet->tile_index];
  
  for (Tile::TileFnPair * near_tile = tile.begin_tiles; near_tile != tile.end_tiles; near_tile++){
    if ((near_tile->first)->tagged) continue;
    double dist = (tile.*(near_tile->second))(jet);
    // cout << "      max info looked at tile " << *near_tile - &_tiles[0] 
    // 	 << ", dist = " << dist << " " << (*near_tile)->max_NN_dist
    // 	 << endl;
    if (dist > (near_tile->first)->max_NN_dist) continue;

    // cout << "      max info tagged tile " << *near_tile - &_tiles[0] << endl;
    (near_tile->first)->tagged = true;
    // get the tile number
    tile_union[n_near_tiles] = near_tile->first - & _tiles[0];
    n_near_tiles++;
  }
}

//--------TMPTMPTMPTMPTMP-----GPS TEMP--------------------
ostream & operator<<(ostream & ostr, const TiledJet & jet) {
  ostr << "j" << setw(3) << jet._jets_index << ":pt2,rap,phi=" ; ostr.flush();
  ostr     << jet.kt2 << ","; ostr.flush();
  ostr     << jet.eta << ","; ostr.flush();
  ostr     << jet.phi; ostr.flush();
  ostr     << ", tile=" << jet.tile_index; ostr.flush();
  return ostr;
}


// //----------------------------------------------------------------------
// /// returns a particle's distance to the edge of the specified tile
// inline double LazyTiling9Alt::_distance_to_tile(const TiledJet * bj, const Tile * tile) const {
// 
// //   // Note the careful way of checking the minimum potential deta:
// //   // unlike the phi case below, we don't calculate the distance to the
// //   // centre and subtract spacing/2. This is because of issue of
// //   // boundary tiles, which can extend far beyond spacing/2 in eta. 
// //   // Using the positions of tile centers should instead be safe.
// //   double deta;
// //   if (_tiles[bj->tile_index].eta_centre == tile->eta_centre) deta = 0;
// //   //else   deta = std::abs(bj->eta - tile->eta_centre) - 0.5*_tile_size_eta;
// //   else   deta = std::abs(bj->eta - tile->eta_centre) - _tile_half_size_eta;
// //   // ------
// //   //   |
// //   // A | B
// //   // ------
// //   //   |
// //   // C | D
// //   // ------
// // 
// //   double dphi = std::abs(bj->phi - tile->phi_centre);
// //   if (dphi > pi) dphi = twopi-dphi;
// //   dphi -= _tile_half_size_phi;
// //   //dphi -= 0.5*_tile_size_phi;
// //   if (dphi < 0) dphi = 0;
// // 
// //   return dphi*dphi + deta*deta;
// 
//   return 0.0;
// }




//----------------------------------------------------------------------
/// looks at distance between jetX and jetI and updates the NN
/// information if relevant; also pushes identity of jetI onto
/// the vector of jets for minheap, to signal that it will have
/// to be handled later.
///
/// GPS TEMP GPS TMP: REMOVE THIS LATER: EVEN LABELLED AS INLINE, THE
/// CALL ADDS A SUBSTANTIAL PENALTY...
inline void LazyTiling9Alt::_update_jetX_jetI_NN(TiledJet * jetX, TiledJet * jetI, vector<TiledJet *> & jets_for_minheap) {
  double dist = _bj_dist(jetI,jetX);
  if (dist < jetI->NN_dist) {
    if (jetI != jetX) {
      jetI->NN_dist = dist;
      jetI->NN = jetX;
      // label jetI as needing heap action...
      if (!jetI->minheap_update_needed()) {
	jetI->label_minheap_update_needed();
	jets_for_minheap.push_back(jetI);
      }
    }
  }
  if (dist < jetX->NN_dist) {
    if (jetI != jetX) {
      jetX->NN_dist = dist;
      jetX->NN      = jetI;}
  }
}


inline void LazyTiling9Alt::_set_NN(TiledJet * jetI, 
                            vector<TiledJet *> & jets_for_minheap) {
  jetI->NN_dist = _R2;
  jetI->NN      = NULL;
  // label jetI as needing heap action...
  if (!jetI->minheap_update_needed()) {
    jetI->label_minheap_update_needed();
    jets_for_minheap.push_back(jetI);}
  // now go over tiles that are neighbours of I (include own tile)
  Tile * tile_ptr = &_tiles[jetI->tile_index];
  //if (tile_ptr->is_near_zero_phi(_tile_half_size_phi)) {
    for (Tile::TileFnPair * near_tile  = tile_ptr->begin_tiles; 
         near_tile != tile_ptr->end_tiles; near_tile++) {
      // for own tile, this will be zero automatically: should we be clever
      // and skip the test? (With some doubling of code?)
      if (jetI->NN_dist < (tile_ptr->*(near_tile->second))(jetI)) continue;
      // and then over the contents of that tile
      for (TiledJet * jetJ  = (near_tile->first)->head; 
           jetJ != NULL; jetJ = jetJ->next) {
        double dist = _bj_dist(jetI,jetJ);
        if (dist < jetI->NN_dist && jetJ != jetI) {
          jetI->NN_dist = dist; jetI->NN = jetJ;
        }
      }
    }
  // } else {
  //   // second copy that exploits the fact that for this tile we needn't worry
  //   // about periodicity
  //   for (Tile::TileFnPair * near_tile  = tile_ptr->begin_tiles; 
  //        near_tile != tile_ptr->end_tiles; near_tile++) {
  //     // for own tile, this will be zero automatically: should we be clever
  //     // and skip the test? (With some doubling of code?)
  //     if (jetI->NN_dist < tile_ptr->(*(near_tile->second)(jetI))) continue;
  //     // and then over the contents of that tile
  //     for (TiledJet * jetJ  = (*near_tile)->head; 
  //          jetJ != NULL; jetJ = jetJ->next) {
  //       double dist = _bj_dist_not_periodic(jetI,jetJ);
  //       if (dist < jetI->NN_dist && jetJ != jetI) {
  //         jetI->NN_dist = dist; jetI->NN = jetJ;
  //       }
  //     }
  //   }
  // }
}


void LazyTiling9Alt::run() {

  //_initialise_tiles();

  int n = _jets.size();
  TiledJet * briefjets = new TiledJet[n];
  TiledJet * jetA = briefjets, * jetB;
  TiledJet oldB;
  

  // will be used quite deep inside loops, but declare it here so that
  // memory (de)allocation gets done only once
  vector<int> tile_union(3*n_tile_neighbours);
  
  // initialise the basic jet info 
  for (int i = 0; i< n; i++) {
    _tj_set_jetinfo(jetA, i);
    //cout << i<<": "<<jetA->tile_index<<"\n";
    jetA++; // move on to next entry of briefjets
  }
  TiledJet * head = briefjets; // a nicer way of naming start

  // set up the initial nearest neighbour information
  vector<Tile>::iterator tile;
  for (tile = _tiles.begin(); tile != _tiles.end(); tile++) {
    // first do it on this tile
    for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
      for (jetB = tile->head; jetB != jetA; jetB = jetB->next) {
	double dist = _bj_dist_not_periodic(jetA,jetB);
	if (dist < jetA->NN_dist) {jetA->NN_dist = dist; jetA->NN = jetB;}
	if (dist < jetB->NN_dist) {jetB->NN_dist = dist; jetB->NN = jetA;}
      }
    }
    for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
      if (jetA->NN_dist > tile->max_NN_dist) tile->max_NN_dist = jetA->NN_dist;
    }
  }
  for (tile = _tiles.begin(); tile != _tiles.end(); tile++) {
    if (tile->use_periodic_delta_phi) {
      // then do it for RH tiles; 
      for (Tile::TileFnPair * RTileFnPair = tile->RH_tiles; 
           RTileFnPair != tile->end_tiles; RTileFnPair++) {
        Tile *RTile = RTileFnPair->first;
        for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
          double dist_to_tile = ((*tile).*(RTileFnPair->second))(jetA);
          // it only makes sense to do a tile if jetA is close enough to the Rtile
          // either for a jet in the Rtile to be closer to jetA than it's current NN
          // or if jetA could be closer to something in the Rtile than the largest
          // NN distance within the RTile.
          //
          // GPS note: also tried approach where we perform only the
          //           first test and run over all surrounding tiles
          //           (not just RH ones). The test is passed less
          //           frequently, but one is running over more tiles
          //           and on balance, for the trial event we used, it's
          //           a bit slower.
          bool relevant_for_jetA  = dist_to_tile <= jetA->NN_dist;
          bool relevant_for_RTile = dist_to_tile <= RTile->max_NN_dist;
          if (relevant_for_jetA || relevant_for_RTile) {
            for (jetB = RTile->head; jetB != NULL; jetB = jetB->next) {
              double dist = _bj_dist(jetA,jetB);
              if (dist < jetA->NN_dist) {jetA->NN_dist = dist; jetA->NN = jetB;}
              if (dist < jetB->NN_dist) {jetB->NN_dist = dist; jetB->NN = jetA;}
            }
          } 
        }
      }
    } else {
      // this second version of the code uses the faster
      // "not_periodic" version because it knows that the tile is
      // sufficiently far from the edge.
      for (Tile::TileFnPair* RTileFnPair = tile->RH_tiles;
           RTileFnPair != tile->end_tiles; RTileFnPair++) {
        Tile *RTile = RTileFnPair->first;
        for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
          double dist_to_tile = ((*tile).*(RTileFnPair->second))(jetA);
          bool relevant_for_jetA  = dist_to_tile <= jetA->NN_dist;
          bool relevant_for_RTile = dist_to_tile <= RTile->max_NN_dist;
          if (relevant_for_jetA || relevant_for_RTile) {
            for (jetB = RTile->head; jetB != NULL; jetB = jetB->next) {
              double dist = _bj_dist_not_periodic(jetA,jetB);
              if (dist < jetA->NN_dist) {jetA->NN_dist = dist; jetA->NN = jetB;}
              if (dist < jetB->NN_dist) {jetB->NN_dist = dist; jetB->NN = jetA;}
            }
          } 
        }
      }
    }
    // no need to do it for LH tiles, since they are implicitly done
    // when we set NN for both jetA and jetB on the RH tiles.
  }
  // Now update the max_NN_dist within each tile. Not strictly
  // necessary, because existing max_NN_dist is an upper bound.  but
  // costs little and may give some efficiency gain later.
  for (tile = _tiles.begin(); tile != _tiles.end(); tile++) {
    tile->max_NN_dist = 0;
    for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
      if (jetA->NN_dist > tile->max_NN_dist) tile->max_NN_dist = jetA->NN_dist;
    }
  }


  vector<double> diJs(n);
  for (int i = 0; i < n; i++) {
    diJs[i] = _bj_diJ(&briefjets[i]);
    briefjets[i].label_minheap_update_done();
  }
  MinHeap minheap(diJs);
  // have a stack telling us which jets we'll have to update on the heap
  vector<TiledJet *> jets_for_minheap;
  jets_for_minheap.reserve(n); 

  // now run the recombination loop
  int history_location = n-1;
  while (n > 0) {

    double diJ_min = minheap.minval() *_invR2;
    jetA = head + minheap.minloc();

    // do the recombination between A and B
    history_location++;
    jetB = jetA->NN;

    if (jetB != NULL) {
      // jet-jet recombination
      // If necessary relabel A & B to ensure jetB < jetA, that way if
      // the larger of them == newtail then that ends up being jetA and 
      // the new jet that is added as jetB is inserted in a position that
      // has a future!
      if (jetA < jetB) {std::swap(jetA,jetB);}

      int nn; // new jet index
      _cs.plugin_record_ij_recombination(jetA->_jets_index, jetB->_jets_index, diJ_min, nn);
      
      // what was jetB will now become the new jet
      _bj_remove_from_tiles(jetA);
      oldB = * jetB;  // take a copy because we will need it...
      _bj_remove_from_tiles(jetB);
      _tj_set_jetinfo(jetB, nn); // cause jetB to become _jets[nn]
                                 // (also registers the jet in the tiling)
    } else {
      // jet-beam recombination
      // get the hist_index
      _cs.plugin_record_iB_recombination(jetA->_jets_index, diJ_min);
      _bj_remove_from_tiles(jetA);
    }

    // remove the minheap entry for jetA
    minheap.remove(jetA-head);

    // first establish the set of tiles over which we are going to
    // have to run searches for updated and new nearest-neighbours --
    // basically a combination of vicinity of the tiles of the two old
    // and one new jet.
    int n_near_tiles = 0;
    _add_untagged_neighbours_to_tile_union_using_max_info(jetA, 
       					   tile_union, n_near_tiles);
    if (jetB != NULL) {
	_add_untagged_neighbours_to_tile_union_using_max_info(&oldB,
							      tile_union,n_near_tiles);
      jetB->label_minheap_update_needed();
      jets_for_minheap.push_back(jetB);
    }


    // Initialise jetB's NN distance as well as updating it for 
    // other particles.
    // Run over all tiles in our union 

    if (jetB != NULL) {
      Tile & jetB_tile = _tiles[jetB->tile_index];
      for (Tile::TileFnPair * near_tile_fn_pair  = jetB_tile.begin_tiles; 
	           near_tile_fn_pair != jetB_tile.end_tiles; near_tile_fn_pair++) {
        Tile * near_tile = near_tile_fn_pair->first;

    	double dist_to_tile = (jetB_tile.*(near_tile_fn_pair->second))(jetB);
        // use <= in next line so that on first tile, relevant_for_jetB is 
        // set to true
    	bool relevant_for_jetB  = dist_to_tile <= jetB->NN_dist;
    	bool relevant_for_near_tile = dist_to_tile <= near_tile->max_NN_dist;
        bool relevant = relevant_for_jetB || relevant_for_near_tile;
        // this first option decides exactly what loop to do based on whether 
        // the near tile was tagged. You'd think it's more efficient, but
        // not necessarily...
        if (relevant) {
          if (near_tile->tagged) {
            for (TiledJet * jetI = near_tile->head; jetI != NULL; jetI = jetI->next) {
              if (jetI->NN == jetA || jetI->NN == jetB) _set_NN(jetI, jets_for_minheap);
              _update_jetX_jetI_NN(jetB, jetI, jets_for_minheap);
            }
          near_tile->tagged = false;
          } else {
            for (TiledJet * jetI = near_tile->head; jetI != NULL; jetI = jetI->next) {
              _update_jetX_jetI_NN(jetB, jetI, jets_for_minheap);
            }
          }
        }

        // this second option does everything independently of whether the near tile
        // was tagged -- somehow you'd expect it to be slower, but it may actually be
        // marginally faster.
        // if (relevant_for_jetB || relevant_for_near_tile) {
        //   for (TiledJet * jetI = (*near_tile)->head; jetI != NULL; jetI = jetI->next) {
        // 
        //     if (jetI->NN == jetA || (jetI->NN == jetB && jetB != NULL)) {
        //       _set_NN(jetI, jets_for_minheap);
        //     }
        // 
        //     _update_jetX_jetI_NN(jetB, jetI, jets_for_minheap);
	//     // -- Keep this old inline code for later speed tests
        //     // double dist = _bj_dist(jetI,jetB);
        //     // if (dist < jetI->NN_dist) {
        //     //   if (jetI != jetB) {
        //     //     jetI->NN_dist = dist;
        //     //     jetI->NN = jetB;
        //     //     // label jetI as needing heap action...
        //     //     if (!jetI->minheap_update_needed()) {
        //     //       jetI->label_minheap_update_needed();
        //     //       jets_for_minheap.push_back(jetI);
        //     //     }
        //     //   }
        //     // }
        //     // if (dist < jetB->NN_dist) {
        //     //   if (jetI != jetB) {
        //     //     jetB->NN_dist = dist;
        //     //     jetB->NN      = jetI;}
        //     // }
        //   }
        //   (*near_tile)->tagged = false;
        // }
      }
    }

    // now run over the tiles that were tagged earlier and that we haven't yet
    // had a change to visit.
    for (int itile = 0; itile < n_near_tiles; itile++) {
      Tile * tile_ptr = &_tiles[tile_union[itile]];
      if (!tile_ptr->tagged) continue; // because earlier loop may have undone the tag
      tile_ptr->tagged = false;
      // run over all jets in the current tile
      for (TiledJet * jetI = tile_ptr->head; jetI != NULL; jetI = jetI->next) {
        // see if jetI had jetA or jetB as a NN -- if so recalculate the NN
        if (jetI->NN == jetA || (jetI->NN == jetB && jetB != NULL)) {
          _set_NN(jetI, jets_for_minheap);
        }
      }
    }

    // deal with jets whose minheap entry needs updating
    //if (verbose) cout << "  jets whose NN was modified: " << endl;
    while (jets_for_minheap.size() > 0) {
      TiledJet * jetI = jets_for_minheap.back(); 
      jets_for_minheap.pop_back();
      minheap.update(jetI-head, _bj_diJ(jetI));
      jetI->label_minheap_update_done();
      // handle max_NN_dist update for all jets that might have
      // seen a change (increase) of distance
      Tile & tile_I = _tiles[jetI->tile_index];
      if (tile_I.max_NN_dist < jetI->NN_dist) tile_I.max_NN_dist = jetI->NN_dist;
    }
    n--;
  }

  // final cleaning up;
  delete[] briefjets;
}


FASTJET_END_NAMESPACE
