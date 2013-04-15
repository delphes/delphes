#ifndef __NNH_HH__
#define __NNH_HH__

//STARTHEADER
// $Id: NNH.hh 2891 2012-06-15 12:47:59Z soyez $
//
// Copyright (c) 2009, Matteo Cacciari, Gavin Salam and Gregory Soyez
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
//  development and are described in hep-ph/0512210. If you use
//  FastJet as part of work towards a scientific publication, please
//  include a citation to the FastJet paper.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//ENDHEADER

#include<fastjet/ClusterSequence.hh>


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

/// @ingroup advanced_usage
/// \class _NoInfo
/// dummy class, used as a default template argument
class _NoInfo {};

/// @ingroup advanced_usage
/// \class NNHInfo
/// template that will help initialise a BJ with a PseudoJet and extra information
template<class I> class NNHInfo {
public:
  NNHInfo()         : _info(NULL) {}
  NNHInfo(I * info) : _info(info) {}
  template<class NNBJ> void init_jet(NNBJ * briefjet, const fastjet::PseudoJet & jet, int index) { briefjet->init(jet, index, _info);}
private:
  I * _info;
};

/// @ingroup advanced_usage
/// Specialisation of NNHInfo for cases where there is no extra info
template<> class NNHInfo<_NoInfo>  {
public:
  NNHInfo()           {}
  NNHInfo(_NoInfo * ) {}
  template<class NNBJ> void init_jet(NNBJ * briefjet, const fastjet::PseudoJet & jet, int index) { briefjet->init(jet, index);}
};


//----------------------------------------------------------------------
/// @ingroup advanced_usage
/// \class NNH
/// Help solve closest pair problems with generic interparticle and
/// beam distance.
///
/// Class to help solve closest pair problems with generic interparticle
/// distances and a beam distance, using Anderberg's Nearest Neighbour
/// Heuristic.
///
/// It is templated with a BJ (brief jet) class --- BJ should
/// basically cache the minimal amount of information that is needed
/// to efficiently calculate interparticle distances and particle-beam
/// distances.
///
/// This class can be used with or without an extra "Information" template, 
/// i.e. NNB<BJ> or NNH<BJ,I>
/// 
/// For the NNH<BJ> version of the class to function, BJ must provide 
/// three member functions
///  
///  - void   BJ::init(const PseudoJet & jet);       // initialise with a PseudoJet
///  - double BJ::distance(const BJ * other_bj_jet); // distance between this and other_bj_jet
///  - double BJ::beam_distance()                  ; // distance to the beam
///
/// For the NNH<BJ,I> version to function, the BJ::init(...) member
/// must accept an extra argument
///
///  - void   BJ::init(const PseudoJet & jet, I * info);   // initialise with a PseudoJet + info
///
/// where info might be a pointer to a class that contains, e.g., information
/// about R, or other parameters of the jet algorithm 
///
/// For an example of how the NNH<BJ> class is used, see the Jade (and
/// EECambridge) plugins
///
/// NB: the NNH algorithm is expected N^2, but has a worst case of
/// N^3. Many QCD problems tend to place one closer to the N^3 end of
/// the spectrum than one would like. There is scope for further
/// progress (cf Eppstein, Cardinal & Eppstein), nevertheless the
/// current class is already significantly faster than standard N^3
/// implementations.
///
///
/// Implementation note: this class derives from NNHInfo, which deals
/// with storing any global information that is needed during the clustering

template<class BJ, class I = _NoInfo> class NNH : public NNHInfo<I> {
public:

  /// constructor with an initial set of jets (which will be assigned indices
  /// 0 ... jets.size()-1
  NNH(const std::vector<PseudoJet> & jets) {start(jets);}
  NNH(const std::vector<PseudoJet> & jets, I * info) : NNHInfo<I>(info) {start(jets);}
  
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
  ~NNH() {
    delete[] briefjets;
  }
    
private:
  class NNBJ; // forward declaration
  
  /// establish the nearest neighbour for jet, and cross check constistency
  /// of distances for the other jets that are encountered. Assumes
  /// jet not contained within begin...end
  void set_NN_crosscheck(NNBJ * jet, NNBJ * begin, NNBJ * end);

  /// establish the nearest neighbour for jet; don't cross check other jets'
  /// distances and allow jet to be contained within begin...end
  void set_NN_nocross   (NNBJ * jet, NNBJ * begin, NNBJ * end);

  /// contains the briefjets
  NNBJ * briefjets;

  /// semaphores for the current extent of our structure
  NNBJ * head, * tail;

  /// currently active number of jets
  int n;

  /// where_is[i] contains a pointer to the jet with index i
  std::vector<NNBJ *> where_is;

  /// a class that wraps around the BJ, supplementing it with extra information
  /// such as pointers to neighbours, etc.
  class NNBJ : public BJ {
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
      NN_dist = BJ::beam_distance();
      NN = NULL;
    }
    int index() const {return _index;}
    
    double NN_dist;
    NNBJ * NN;
    
  private:
    int _index;
  };

};



//----------------------------------------------------------------------
template<class BJ, class I> void NNH<BJ,I>::start(const std::vector<PseudoJet> & jets) {
  n = jets.size();
  briefjets = new NNBJ[n];
  where_is.resize(2*n);

  NNBJ * jetA = briefjets;
  
  // initialise the basic jet info 
  for (int i = 0; i< n; i++) {
    //jetA->init(jets[i], i);
    this->init_jet(jetA, jets[i], i);
    where_is[i] = jetA;
    jetA++; // move on to next entry of briefjets
  }
  tail = jetA; // a semaphore for the end of briefjets
  head = briefjets; // a nicer way of naming start

  // now initialise the NN distances: jetA will run from 1..n-1; and
  // jetB from 0..jetA-1
  for (jetA = head + 1; jetA != tail; jetA++) {
    // set NN info for jetA based on jets running from head..jetA-1,
    // checking in the process whether jetA itself is an undiscovered
    // NN of one of those jets.
    set_NN_crosscheck(jetA, head, jetA);
  }
  //std::cout << "OOOO "  << briefjets[1].NN_dist << " " << briefjets[1].NN - briefjets << std::endl;
}


//----------------------------------------------------------------------
template<class BJ, class I> double NNH<BJ,I>::dij_min(int & iA, int & iB) {
  // find the minimum of the diJ on this round
  double diJ_min = briefjets[0].NN_dist;
  int diJ_min_jet = 0;
  for (int i = 1; i < n; i++) {
    if (briefjets[i].NN_dist < diJ_min) {
      diJ_min_jet = i; 
      diJ_min  = briefjets[i].NN_dist;
    }
  }
  
  // return information to user about recombination
  NNBJ * jetA = & briefjets[diJ_min_jet];
  //std::cout << jetA - briefjets << std::endl; 
  iA = jetA->index();
  iB = jetA->NN ? jetA->NN->index() : -1;
  return diJ_min;
}


//----------------------------------------------------------------------
// remove jetA from the list
template<class BJ, class I> void NNH<BJ,I>::remove_jet(int iA) {
  NNBJ * jetA = where_is[iA];
  // now update our nearest neighbour info and diJ table
  // first reduce size of table
  tail--; n--;
  // Copy last jet contents and diJ info into position of jetA
  *jetA = *tail;
  // update the info on where the given index is stored
  where_is[jetA->index()] = jetA;

  for (NNBJ * jetI = head; jetI != tail; jetI++) {
    // see if jetI had jetA or jetB as a NN -- if so recalculate the NN
    if (jetI->NN == jetA) set_NN_nocross(jetI, head, tail);

    // if jetI's NN is the new tail then relabel it so that it becomes jetA
    if (jetI->NN == tail) {jetI->NN = jetA;}
  }
}


//----------------------------------------------------------------------
template<class BJ, class I> void NNH<BJ,I>::merge_jets(int iA, int iB, 
					const PseudoJet & jet, int index) {

  NNBJ * jetA = where_is[iA];
  NNBJ * jetB = where_is[iB];

  // If necessary relabel A & B to ensure jetB < jetA, that way if
  // the larger of them == newtail then that ends up being jetA and 
  // the new jet that is added as jetB is inserted in a position that
  // has a future!
  if (jetA < jetB) std::swap(jetA,jetB);

  // initialise jetB based on the new jet
  //jetB->init(jet, index);
  this->init_jet(jetB, jet, index);
  // and record its position (making sure we have the space)
  if (index >= int(where_is.size())) where_is.resize(2*index);
  where_is[jetB->index()] = jetB;

  // now update our nearest neighbour info
  // first reduce size of table
  tail--; n--;
  // Copy last jet contents into position of jetA
  *jetA = *tail;
  // update the info on where the tail's index is stored
  where_is[jetA->index()] = jetA;

  for (NNBJ * jetI = head; jetI != tail; jetI++) {
    // see if jetI had jetA or jetB as a NN -- if so recalculate the NN
    if (jetI->NN == jetA || jetI->NN == jetB) {
      set_NN_nocross(jetI, head, tail);
    } 

    // check whether new jetB is closer than jetI's current NN and
    // if need be update things
    double dist = jetI->distance(jetB);
    if (dist < jetI->NN_dist) {
      if (jetI != jetB) {
	jetI->NN_dist = dist;
	jetI->NN = jetB;
      }
    }
    if (dist < jetB->NN_dist) {
      if (jetI != jetB) {
	jetB->NN_dist = dist;
	jetB->NN      = jetI;
      }
    }

    // if jetI's NN is the new tail then relabel it so that it becomes jetA
    if (jetI->NN == tail) {jetI->NN = jetA;}
  }
}


//----------------------------------------------------------------------
// this function assumes that jet is not contained within begin...end
template <class BJ, class I> void NNH<BJ,I>::set_NN_crosscheck(NNBJ * jet, 
		    NNBJ * begin, NNBJ * end) {
  double NN_dist = jet->beam_distance();
  NNBJ * NN      = NULL;
  for (NNBJ * jetB = begin; jetB != end; jetB++) {
    double dist = jet->distance(jetB);
    if (dist < NN_dist) {
      NN_dist = dist;
      NN = jetB;
    }
    if (dist < jetB->NN_dist) {
      jetB->NN_dist = dist;
      jetB->NN = jet;
    }
  }
  jet->NN = NN;
  jet->NN_dist = NN_dist;
}


//----------------------------------------------------------------------
// set the NN for jet without checking whether in the process you might
// have discovered a new nearest neighbour for another jet
template <class BJ, class I>  void NNH<BJ,I>::set_NN_nocross(
                 NNBJ * jet, NNBJ * begin, NNBJ * end) {
  double NN_dist = jet->beam_distance();
  NNBJ * NN      = NULL;
  // if (head < jet) {
  //   for (NNBJ * jetB = head; jetB != jet; jetB++) {
  if (begin < jet) {
    for (NNBJ * jetB = begin; jetB != jet; jetB++) {
      double dist = jet->distance(jetB);
      if (dist < NN_dist) {
	NN_dist = dist;
	NN = jetB;
      }
    }
  }
  // if (tail > jet) {
  //   for (NNBJ * jetB = jet+1; jetB != tail; jetB++) {
  if (end > jet) {
    for (NNBJ * jetB = jet+1; jetB != end; jetB++) {
      double dist = jet->distance (jetB);
      if (dist < NN_dist) {
	NN_dist = dist;
	NN = jetB;
      }
    }
  }
  jet->NN = NN;
  jet->NN_dist = NN_dist;
}




FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh


#endif // __NNH_HH__
