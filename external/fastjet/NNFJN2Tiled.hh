#ifndef __FASTJET_NNFJN2TILED_HH__
#define __FASTJET_NNFJN2TILED_HH__

//FJSTARTHEADER
// $Id$
//
// Copyright (c) 2016-2024, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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

#include "fastjet/NNBase.hh"
#include "fastjet/internal/TilingExtent.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//----------------------------------------------------------------------
/// @ingroup advanced_usage
/// \class NNFJN2Tiled
///
/// Helps solve closest pair problems with factorised interparticle
/// and beam distances (ie satisfying the FastJet lemma) that are on
/// a cylindrical geometry and allow tiling.
///
/// (see NNBase.hh for an introductory description)
///
/// This variant provides an implementation based on the N2Tiled
/// clustering strategy in FastJet. As for the NNFJN2Plain case, the
/// interparticle and beam distances should be of the form
///
/// \code
///   dij = min(mom_factor(i), mom_factor(j)) * geometrical_distance(i,j)
///   diB = mom_factor(i) * geometrical_beam_distance(i)
/// \endcode
///
/// Additionally, the NNFJN2Tiled class takes a tile_size parameter
/// that controls the size of the tiles. It must be such that, for any
/// two points in non-neighbouring (and non-identical) tiles, the
/// geometrical distance between the 2 points is larger than the
/// geometrical beam distance of each of the 2 points.
///
/// It is templated with a BJ (brief jet) class and can be used with or
/// without an extra "Information" template, i.e. NNFJN2Tiled<BJ> or
/// NNFJN2Tiled<BJ,I>
///
/// For the NNFJN2Tiled<BJ> version of the class to function, BJ must provide 
/// three member functions
///  
/// \code
///   void   BJ::init(const PseudoJet & jet);                   // initialise with a PseudoJet
///   double BJ::geometrical_distance(const BJ * other_bj_jet); // distance between this and other_bj_jet (geometrical part)
///   double BJ::geometrical_beam_distance();                   // distance to the beam (geometrical part)
///   double BJ::momentum_factor();                             // extra momentum factor
/// \endcode
///
/// For the NNFJN2Tiled<BJ,I> version to function, the BJ::init(...) member
/// must accept an extra argument
///
/// \code
///   void BJ::init(const PseudoJet & jet, I * info);   // initialise with a PseudoJet + info
/// \endcode
///
/// NOTE: THE DISTANCE MUST BE SYMMETRIC I.E. SATISFY
/// \code
///     a.geometrical_distance(b) == b.geometrical_distance(a)
/// \endcode
///
/// Finally, the BJ class needs to provide access to the variables used
/// for the rectangular tiling:
///
/// \code
///   double BJ::rap(); // rapidity-like variable
///   double BJ::phi(); // azimutal-angle-like variable (should be > -2pi)
/// \endcode
///
/// Note that you are strongly advised to add the following lines
/// to your BJ class to allow it to be used also with NNH:
///
/// \code
///   /// make this BJ class compatible with the use of NNH
///   double BJ::distance(const BJ * other_bj_jet){
///     double mom1 = momentum_factor();
///     double mom2 = other_bj_jet->momentum_factor();
///     return (mom1<mom2 ? mom1 : mom2) * geometrical_distance(other_bj_jet);
///   }
///   double BJ::beam_distance(){
///     return momentum_factor() * geometrical_beam_distance();
///   }
/// \endcode
///
template<class BJ, class I = _NoInfo> class NNFJN2Tiled : public NNBase<I> {
public:

  /// constructor with an initial set of jets (which will be assigned indices
  /// `0...jets.size()-1`)
  NNFJN2Tiled(const std::vector<PseudoJet> & jets, double requested_tile_size)
    : NNBase<I>(),     _requested_tile_size(requested_tile_size) {start(jets);}
  NNFJN2Tiled(const std::vector<PseudoJet> & jets, double requested_tile_size, I * info)
    : NNBase<I>(info), _requested_tile_size(requested_tile_size) {start(jets);}
  
  void start(const std::vector<PseudoJet> & jets);

  /// return the dij_min and indices iA, iB, for the corresponding jets.
  /// If iB < 0 then iA recombines with the beam
  double dij_min(int & iA, int & iB);

  /// remove the jet pointed to by index iA
  void remove_jet(int iA);

  /// merge the jets pointed to by indices A and B and replace them with
  /// jet, assigning it an index jet_index.
  void merge_jets(int iA, int iB, const PseudoJet & jet, int jet_index);

  /// a destructor
  ~NNFJN2Tiled() {
    delete[] briefjets;
    delete[] diJ;
  }
    
private:
  class TiledJet; // forward declaration
  class Tile;
  class diJ_plus_link;

  // Set up the tiles:
  void _initialise_tiles(const std::vector<PseudoJet> & particles);

  // return the full distance of a particle to its NN
  inline double _compute_diJ(const TiledJet * const jet) const {
    double mom_fact = jet->momentum_factor();
    if (jet->NN != NULL) {
      double other_mom_fact = jet->NN->momentum_factor();
      if (other_mom_fact < mom_fact) {mom_fact = other_mom_fact;}
    }
    return jet->NN_dist * mom_fact;
  }

  // reasonably robust return of tile index given irap and iphi, in particular
  // it works even if iphi is negative
  inline int _tile_index (int irap, int iphi) const {
    // note that (-1)%n = -1 so that we have to add _n_tiles_phi
    // before performing modulo operation
    return (irap-_tiles_irap_min)*_n_tiles_phi
                  + (iphi+_n_tiles_phi) % _n_tiles_phi;
  }

  int  _tile_index(const double rap, const double phi) const;
  void _tiledjet_set_jetinfo ( TiledJet * const tiled_jet, const PseudoJet &jet, int index);
  void _bj_remove_from_tiles(TiledJet * const jet);
  void _initialise_tiles();
  void _print_tiles(TiledJet * briefjets ) const;
  void _add_neighbours_to_tile_union(const int tile_index, int & n_near_tiles) const;
  void _add_untagged_neighbours_to_tile_union(const int tile_index, int & n_near_tiles);


  /// contains the briefjets
  TiledJet * briefjets;

  /// semaphores for the current extent of our structure
  TiledJet * head;

  /// currently active number of jets
  int n;

  /// where_is[i] contains a pointer to the jet with index i
  std::vector<TiledJet *> where_is;

  /// helper to keep tracks of tiles to be checked for updates
  std::vector<int> tile_union;

  /// a table containing all the (full) distances to each object's NN
  diJ_plus_link * diJ;

  /// tiling information
  std::vector<Tile> _tiles;
  double _requested_tile_size;
  double _tiles_rap_min, _tiles_rap_max;
  double _tile_size_rap, _tile_size_phi;
  int    _n_tiles_phi,_tiles_irap_min,_tiles_irap_max;
  
  /// a class that wraps around the BJ, supplementing it with extra information
  /// such as pointers to neighbours, etc.
  class TiledJet : public BJ {
  public:
    void init(const PseudoJet & jet, int index_in) {
      BJ::init(jet);
      other_init(index_in);
    }
    void init(const PseudoJet & jet, int index_in, I * info) {
      BJ::init(jet, info);
      other_init(index_in);
    }
    void other_init(int index_in) {
      _index = index_in;
      NN_dist = BJ::geometrical_beam_distance();
      NN = NULL;
    }
    int jet_index() const {return _index;}
    
    double NN_dist;
    TiledJet * NN, *previous, * next;
    int tile_index, diJ_posn;
    // routines that are useful in the minheap version of tiled
    // clustering ("misuse" the otherwise unused diJ_posn, so as
    // to indicate whether jets need to have their minheap entries
    // updated).
    inline void label_minheap_update_needed() {diJ_posn = 1;}
    inline void label_minheap_update_done()   {diJ_posn = 0;}
    inline bool minheap_update_needed() const {return diJ_posn==1;}
    
  private:
    int _index;
  };

  /// number of neighbours that a tile will have (rectangular geometry
  /// gives 9 neighbours).
  static const int n_tile_neighbours = 9;
  //----------------------------------------------------------------------
  /// The fundamental structures to be used for the tiled N^2 algorithm
  /// (see CCN27-44 for some discussion of pattern of tiling)
  class Tile {
  public:
    /// pointers to neighbouring tiles, including self
    Tile *   begin_tiles[n_tile_neighbours]; 
    /// neighbouring tiles, excluding self
    Tile **  surrounding_tiles; 
    /// half of neighbouring tiles, no self
    Tile **  RH_tiles;  
    /// just beyond end of tiles
    Tile **  end_tiles; 
    /// start of list of BriefJets contained in this tile
    TiledJet * head;    
    /// sometimes useful to be able to tag a tile
    bool     tagged;    
  };

  // structure that holds the real, full, distance (as well as a pointer to the corresponding TiledJet)
  class diJ_plus_link {
  public:
    double     diJ; // the distance
    TiledJet * jet; // the jet (i) for which we've found this distance
                    // (whose NN will the J).
  };

};



//----------------------------------------------------------------------
template<class BJ, class I> void NNFJN2Tiled<BJ,I>::start(const std::vector<PseudoJet> & jets) {

  _initialise_tiles(jets);

  n = jets.size();

  briefjets = new TiledJet[n];
  where_is.resize(2*n);

  TiledJet * jetA = briefjets, * jetB;

  // will be used quite deep inside loops, but declare it here so that
  // memory (de)allocation gets done only once
  tile_union.resize(3*n_tile_neighbours);
  
  // initialise the basic jet info 
  for (int i = 0; i< n; i++) {
    _tiledjet_set_jetinfo(jetA, jets[i], i);
    where_is[i] = jetA;
    jetA++; // move on to next entry of briefjets
  }

  head = briefjets; // a nicer way of naming start

  // set up the initial nearest neighbour information
  typename std::vector<Tile>::const_iterator tile;
  for (tile = _tiles.begin(); tile != _tiles.end(); tile++) {
    // first do it on this tile
    for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
      for (jetB = tile->head; jetB != jetA; jetB = jetB->next) {
	double dist = jetA->geometrical_distance(jetB);
	if (dist < jetA->NN_dist) {jetA->NN_dist = dist; jetA->NN = jetB;}
	if (dist < jetB->NN_dist) {jetB->NN_dist = dist; jetB->NN = jetA;}
      }
    }
    // then do it for RH tiles
    for (Tile ** RTile = tile->RH_tiles; RTile != tile->end_tiles; RTile++) {
      for (jetA = tile->head; jetA != NULL; jetA = jetA->next) {
	for (jetB = (*RTile)->head; jetB != NULL; jetB = jetB->next) {
	  double dist = jetA->geometrical_distance(jetB);
	  if (dist < jetA->NN_dist) {jetA->NN_dist = dist; jetA->NN = jetB;}
	  if (dist < jetB->NN_dist) {jetB->NN_dist = dist; jetB->NN = jetA;}
	}
      }
    }
    // no need to do it for LH tiles, since they are implicitly done
    // when we set NN for both jetA and jetB on the RH tiles.
  }
  
  diJ = new diJ_plus_link[n];
  jetA = head;
  for (int i = 0; i < n; i++) {
    diJ[i].diJ = _compute_diJ(jetA); // kt distance * R^2
    diJ[i].jet = jetA;  // our compact diJ table will not be in	     
    jetA->diJ_posn = i; // one-to-one corresp. with non-compact jets,
                        // so set up bi-directional correspondence here.
    jetA++; // have jetA follow i 
  }
}


//----------------------------------------------------------------------
template<class BJ, class I> double NNFJN2Tiled<BJ,I>::dij_min(int & iA, int & iB) {
  // find the minimum of the diJ on this round
  diJ_plus_link * best, *stop; // pointers a bit faster than indices
                               // could use best to keep track of diJ
                               // min, but it turns out to be
                               // marginally faster to have a separate
                               // variable (avoids n dereferences at
                               // the expense of n/2 assignments).
  double diJ_min = diJ[0].diJ; // initialise the best one here.
  best = diJ;                  // and here
  stop = diJ+n;
  for (diJ_plus_link * here = diJ+1; here != stop; here++) {
    if (here->diJ < diJ_min) {best = here; diJ_min  = here->diJ;}
  }

  // return information to user about recombination
  TiledJet * jetA = best->jet;
  iA = jetA->jet_index();
  iB = jetA->NN ? jetA->NN->jet_index() : -1;
  return diJ_min;
}


//----------------------------------------------------------------------
// remove jetA from the list
template<class BJ, class I> void NNFJN2Tiled<BJ,I>::remove_jet(int iA) {
  TiledJet * jetA = where_is[iA];

  _bj_remove_from_tiles(jetA);

  // first establish the set of tiles over which we are going to
  // have to run searches for updated and new nearest-neighbours --
  // basically a combination of vicinity of the tiles of the two old
  // and one new jet.
  int n_near_tiles = 0;
  _add_untagged_neighbours_to_tile_union(jetA->tile_index, n_near_tiles);

  // now update our nearest neighbour info and diJ table
  // first reduce size of table
  n--;
  // then compactify the diJ by taking the last of the diJ and copying
  // it to the position occupied by the diJ for jetA
  diJ[n].jet->diJ_posn = jetA->diJ_posn;
  diJ[jetA->diJ_posn] = diJ[n];

  // updating other particles' NN.
  // Run over all tiles in our union 
  for (int itile = 0; itile < n_near_tiles; itile++) {
    Tile * tile_ptr = &_tiles[tile_union[itile]];
    tile_ptr->tagged = false; // reset tag, since we're done with unions
    // run over all jets in the current tile
    for (TiledJet * jetI = tile_ptr->head; jetI != NULL; jetI = jetI->next) {
      // see if jetI had jetA or jetB as a NN -- if so recalculate the NN
      if (jetI->NN == jetA) {
        jetI->NN_dist = jetI->geometrical_beam_distance();
        jetI->NN      = NULL;
        // now go over tiles that are neighbours of I (include own tile)
        for (Tile ** near_tile  = tile_ptr->begin_tiles; 
             near_tile != tile_ptr->end_tiles; near_tile++) {
          // and then over the contents of that tile
          for (TiledJet * jetJ  = (*near_tile)->head; jetJ != NULL; jetJ = jetJ->next) {
            double dist = jetI->geometrical_distance(jetJ);
            if (dist < jetI->NN_dist && jetJ != jetI) {
              jetI->NN_dist = dist; jetI->NN = jetJ;
            }
          }
        }
        diJ[jetI->diJ_posn].diJ = _compute_diJ(jetI); // update diJ kt-dist
      }
    }
  }

}


//----------------------------------------------------------------------
template<class BJ, class I> void NNFJN2Tiled<BJ,I>::merge_jets(int iA, int iB, 
					const PseudoJet & jet, int index) {

  TiledJet * jetA = where_is[iA];
  TiledJet * jetB = where_is[iB];

  // jet-jet recombination
  // If necessary relabel A & B to ensure jetB < jetA, that way if
  // the larger of them == newtail then that ends up being jetA and 
  // the new jet that is added as jetB is inserted in a position that
  // has a future!
  if (jetA < jetB) {std::swap(jetA,jetB);}

  // what was jetB will now become the new jet
  _bj_remove_from_tiles(jetA);
  TiledJet oldB = * jetB;  // take a copy because we will need it...
  _bj_remove_from_tiles(jetB);
  _tiledjet_set_jetinfo(jetB, jet, index); // cause jetB to become _jets[nn]
                                           // (also registers the jet in the tiling)
  where_is[index] = jetB;
  
  // first establish the set of tiles over which we are going to
  // have to run searches for updated and new nearest-neighbours --
  // basically a combination of vicinity of the tiles of the two old
  // and one new jet.
  int n_near_tiles = 0;
  _add_untagged_neighbours_to_tile_union(jetA->tile_index, n_near_tiles);
  if (jetB->tile_index != jetA->tile_index) {
    _add_untagged_neighbours_to_tile_union(jetB->tile_index, n_near_tiles);
  }
  if (oldB.tile_index != jetA->tile_index && 
      oldB.tile_index != jetB->tile_index) {
    _add_untagged_neighbours_to_tile_union(oldB.tile_index, n_near_tiles);
  }

  // now update our nearest neighbour info and diJ table
  // first reduce size of table
  n--;
  // then compactify the diJ by taking the last of the diJ and copying
  // it to the position occupied by the diJ for jetA
  diJ[n].jet->diJ_posn = jetA->diJ_posn;
  diJ[jetA->diJ_posn] = diJ[n];

  // Initialise jetB's NN distance as well as updating it for 
  // other particles.
  // Run over all tiles in our union 
  for (int itile = 0; itile < n_near_tiles; itile++) {
    Tile * tile_ptr = &_tiles[tile_union[itile]];
    tile_ptr->tagged = false; // reset tag, since we're done with unions
    // run over all jets in the current tile
    for (TiledJet * jetI = tile_ptr->head; jetI != NULL; jetI = jetI->next) {
      // see if jetI had jetA or jetB as a NN -- if so recalculate the NN
      if ((jetI->NN == jetA) || (jetI->NN == jetB)) {
        jetI->NN_dist = jetI->geometrical_beam_distance();
        jetI->NN      = NULL;
        // now go over tiles that are neighbours of I (include own tile)
        for (Tile ** near_tile  = tile_ptr->begin_tiles; near_tile != tile_ptr->end_tiles; near_tile++) {
          // and then over the contents of that tile
          for (TiledJet * jetJ  = (*near_tile)->head; jetJ != NULL; jetJ = jetJ->next) {
            double dist = jetI->geometrical_distance(jetJ);
            if (dist < jetI->NN_dist && jetJ != jetI) {
              jetI->NN_dist = dist; jetI->NN = jetJ;
            }
          }
        }
        diJ[jetI->diJ_posn].diJ = _compute_diJ(jetI); // update diJ kt-dist
      }
      // check whether new jetB is closer than jetI's current NN and
      // if jetI is closer than jetB's current (evolving) nearest
      // neighbour. Where relevant update things
      double dist = jetI->geometrical_distance(jetB);
      if (dist < jetI->NN_dist) {
        if (jetI != jetB) {
          jetI->NN_dist = dist;
          jetI->NN = jetB;
          diJ[jetI->diJ_posn].diJ = _compute_diJ(jetI); // update diJ...
        }
      }
      if (dist < jetB->NN_dist) {
        if (jetI != jetB) {
          jetB->NN_dist = dist;
          jetB->NN      = jetI;}
      }
    }
  }

  // finally, register the updated kt distance for B
  diJ[jetB->diJ_posn].diJ = _compute_diJ(jetB);
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
template <class BJ, class I>
void NNFJN2Tiled<BJ,I>::_initialise_tiles(const std::vector<PseudoJet> &particles) {

  // first decide tile sizes (with a lower bound to avoid huge memory use with
  // very small R)
  double default_size = _requested_tile_size>0.1 ? _requested_tile_size : 0.1;
  _tile_size_rap = default_size;
  // it makes no sense to go below 3 tiles in phi -- 3 tiles is
  // sufficient to make sure all pair-wise combinations up to pi in
  // phi are possible
  _n_tiles_phi   = int(floor(twopi/default_size));
  if (_n_tiles_phi<3) _n_tiles_phi = 3;
  _tile_size_phi = twopi / _n_tiles_phi; // >= _Rparam and fits in 2pi

  TilingExtent tiling_analysis(particles);
  _tiles_rap_min = tiling_analysis.minrap();
  _tiles_rap_max = tiling_analysis.maxrap();

  // now adjust the values
  _tiles_irap_min = int(floor(_tiles_rap_min/_tile_size_rap));
  _tiles_irap_max = int(floor( _tiles_rap_max/_tile_size_rap));
  _tiles_rap_min = _tiles_irap_min * _tile_size_rap;
  _tiles_rap_max = _tiles_irap_max * _tile_size_rap;

  // allocate the tiles
  _tiles.resize((_tiles_irap_max-_tiles_irap_min+1)*_n_tiles_phi);

  // now set up the cross-referencing between tiles
  for (int irap = _tiles_irap_min; irap <= _tiles_irap_max; irap++) {
    for (int iphi = 0; iphi < _n_tiles_phi; iphi++) {
      Tile * tile = & _tiles[_tile_index(irap,iphi)];
      // no jets in this tile yet
      tile->head = NULL; // first element of tiles points to itself
      tile->begin_tiles[0] =  tile;
      Tile ** pptile = & (tile->begin_tiles[0]);
      pptile++;
      //
      // set up L's in column to the left of X
      tile->surrounding_tiles = pptile;
      if (irap > _tiles_irap_min) {
	// with the itile subroutine, we can safely run tiles from
	// idphi=-1 to idphi=+1, because it takes care of
	// negative and positive boundaries
	for (int idphi = -1; idphi <=+1; idphi++) {
	  *pptile = & _tiles[_tile_index(irap-1,iphi+idphi)];
	  pptile++;
	}	
      }
      // now set up last L (below X)
      *pptile = & _tiles[_tile_index(irap,iphi-1)];
      pptile++;
      // set up first R (above X)
      tile->RH_tiles = pptile;
      *pptile = & _tiles[_tile_index(irap,iphi+1)];
      pptile++;
      // set up remaining R's, to the right of X
      if (irap < _tiles_irap_max) {
	for (int idphi = -1; idphi <= +1; idphi++) {
	  *pptile = & _tiles[_tile_index(irap+1,iphi+idphi)];
	  pptile++;
	}	
      }
      // now put semaphore for end tile
      tile->end_tiles = pptile;
      // finally make sure tiles are untagged
      tile->tagged = false;
    }
  }

}

//----------------------------------------------------------------------
/// return the tile index corresponding to the given rap,phi point
template <class BJ, class I>
int NNFJN2Tiled<BJ,I>::_tile_index(const double rap, const double phi) const {
  int irap, iphi;
  if      (rap <= _tiles_rap_min) {irap = 0;}
  else if (rap >= _tiles_rap_max) {irap = _tiles_irap_max-_tiles_irap_min;}
  else {
    //irap = int(floor((rap - _tiles_rap_min) / _tile_size_rap));
    irap = int(((rap - _tiles_rap_min) / _tile_size_rap));
    // following needed in case of rare but nasty rounding errors
    if (irap > _tiles_irap_max-_tiles_irap_min) {
      irap = _tiles_irap_max-_tiles_irap_min;} 
  }
  // allow for some extent of being beyond range in calculation of phi
  // as well
  //iphi = (int(floor(phi/_tile_size_phi)) + _n_tiles_phi) % _n_tiles_phi;
  // with just int and no floor, things run faster but beware
  iphi = int((phi+twopi)/_tile_size_phi) % _n_tiles_phi;
  return (iphi + irap * _n_tiles_phi);
}

//----------------------------------------------------------------------
template <class BJ, class I>
void NNFJN2Tiled<BJ,I>::_bj_remove_from_tiles(TiledJet * const jet) {
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
// overloaded version which additionally sets up information regarding the
// tiling
template <class BJ, class I>
inline void NNFJN2Tiled<BJ,I>::_tiledjet_set_jetinfo(TiledJet * const tile_jet,
                                                   const PseudoJet &jet, 
                                                   int index) {
  // the this-> in the next line is required by standard compiler
  // see e.g. http://stackoverflow.com/questions/10639053/name-lookups-in-c-templates
  this->init_jet(tile_jet, jet, index); 

  // Then do the setup specific to the tiled case.

  // Find out which tile it belonds to
  tile_jet->tile_index = _tile_index(tile_jet->rap(), tile_jet->phi());

  // Insert it into the tile's linked list of jets
  Tile * tile = &_tiles[tile_jet->tile_index];
  tile_jet->previous   = NULL;
  tile_jet->next       = tile->head;
  if (tile_jet->next != NULL) {tile_jet->next->previous = tile_jet;}
  tile->head      = tile_jet;
}

//----------------------------------------------------------------------
/// Add to the vector tile_union the tiles that are in the neighbourhood
/// of the specified tile_index, including itself -- start adding
/// from position n_near_tiles-1, and increase n_near_tiles as
/// you go along (could have done it more C++ like with vector with reserved
/// space, but fear is that it would have been slower, e.g. checking
/// for end of vector at each stage to decide whether to resize it)
template <class BJ, class I>
void NNFJN2Tiled<BJ,I>::_add_neighbours_to_tile_union(const int tile_index, 
                                                    int & n_near_tiles) const {
  for (Tile * const * near_tile = _tiles[tile_index].begin_tiles; 
       near_tile != _tiles[tile_index].end_tiles; near_tile++){
    // get the tile number
    tile_union[n_near_tiles] = *near_tile - & _tiles[0];
    n_near_tiles++;
  }
}

//----------------------------------------------------------------------
/// Like _add_neighbours_to_tile_union, but only adds neighbours if 
/// their "tagged" status is false; when a neighbour is added its
/// tagged status is set to true.
///
/// Note that with a high level of warnings (-pedantic -Wextra -ansi,
/// gcc complains about tile_index maybe being used uninitialised for
/// oldB in ClusterSequence::_minheap_faster_tiled_N2_cluster(). We
/// have explicitly checked that it was harmless so we could disable
/// the gcc warning by hand using the construct below
///
///  #pragma GCC diagnostic push
///  #pragma GCC diagnostic ignored "-Wpragmas"
///  #pragma GCC diagnostic ignored "-Wuninitialized"
///  #pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
///    ...
///  #pragma GCC diagnostic pop
///
/// the @GCC diagnostic push/pop directive was only introduced in
/// gcc-4.6, so for broader usage, we'd need to insert #pragma GCC
/// diagnostic ignored "-Wpragmas" at the top of this file
template <class BJ, class I>
inline void NNFJN2Tiled<BJ,I>::_add_untagged_neighbours_to_tile_union(
               const int tile_index, 
	       int & n_near_tiles)  {
  for (Tile ** near_tile = _tiles[tile_index].begin_tiles; 
       near_tile != _tiles[tile_index].end_tiles; near_tile++){
    if (! (*near_tile)->tagged) {
      (*near_tile)->tagged = true;
      // get the tile number
      tile_union[n_near_tiles] = *near_tile - & _tiles[0];
      n_near_tiles++;
    }
  }
}



FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh


#endif // __FASTJET_NNFJN2TILED_HH__
