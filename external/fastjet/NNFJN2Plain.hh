#ifndef __FASTJET_NNFJN2PLAIN_HH__
#define __FASTJET_NNFJN2PLAIN_HH__

//FJSTARTHEADER
// $Id: NNFJN2Plain.hh 4442 2020-05-05 07:50:11Z soyez $
//
// Copyright (c) 2005-2020, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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

#include <fastjet/NNBase.hh>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//----------------------------------------------------------------------
/// @ingroup advanced_usage
/// \class NNFJN2Plain
///
/// Helps solve closest pair problems with factorised interparticle
/// and beam distances (ie satisfying the FastJet lemma)
///
/// (see NNBase.hh for an introductory description)
///
/// This variant provides an implementation based on the N2Plain
/// clustering strategy in FastJet. The interparticle and beam
/// distances should be of the form
///
/// \code
///   dij = min(mom_factor(i), mom_factor(j)) * geometrical_distance(i,j)
///   diB = mom_factor(i) * geometrical_beam_distance(i)
/// \endcode
///
/// The class is templated with a BJ (brief jet) class and can be used
/// with or without an extra "Information" template,
/// i.e. NNFJN2Plain<BJ> or NNFJN2Plain<BJ,I>
///
/// For the NNH_N2Plain<BJ> version of the class to function, BJ must
/// provide four member functions
///  
/// \code
///   void   BJ::init(const PseudoJet & jet);                   // initialise with a PseudoJet
///   double BJ::geometrical_distance(const BJ * other_bj_jet); // distance between this and other_bj_jet (geometrical part)
///   double BJ::geometrical_beam_distance();                   // distance to the beam (geometrical part)
///   double BJ::momentum_factor();                             // extra momentum factor
/// \endcode
///
/// For the NNH_N2Plain<BJ,I> version to function, the BJ::init(...)
/// member must accept an extra argument
///
/// \code
///   void  BJ::init(const PseudoJet & jet, I * info);   // initialise with a PseudoJet + info
/// \endcode
///
/// NOTE: THE DISTANCE MUST BE SYMMETRIC I.E. SATISFY
/// \code
///     a.geometrical_distance(b) == b.geometrical_distance(a)
/// \endcode
///
/// Note that you are strongly advised to add the following lines to
/// your BJ class to allow it to be used also with NNH:
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
template<class BJ, class I = _NoInfo> class NNFJN2Plain : public NNBase<I> {
public:

  /// constructor with an initial set of jets (which will be assigned indices
  /// `0...jets.size()-1`)
  NNFJN2Plain(const std::vector<PseudoJet> & jets)           : NNBase<I>()     {start(jets);}
  NNFJN2Plain(const std::vector<PseudoJet> & jets, I * info) : NNBase<I>(info) {start(jets);}
  
  /// initialisation from a given list of particles
  virtual void start(const std::vector<PseudoJet> & jets);

  /// returns the dij_min and indices iA, iB, for the corresponding jets.
  /// If iB < 0 then iA recombines with the beam
  double dij_min(int & iA, int & iB);

  /// removes the jet pointed to by index iA
  void remove_jet(int iA);

  /// merges the jets pointed to by indices A and B and replace them with
  /// jet, assigning it an index jet_index.
  void merge_jets(int iA, int iB, const PseudoJet & jet, int jet_index);

  /// a destructor
  ~NNFJN2Plain() {
    delete[] briefjets;
    delete[] diJ;
  }
    
private:
  class NNBJ; // forward declaration

  // return the full distance of a particle to its NN
  inline double compute_diJ(const NNBJ * const jet) const {
    double mom_fact = jet->momentum_factor();
    if (jet->NN != NULL) {
      double other_mom_fact = jet->NN->momentum_factor();
      if (other_mom_fact < mom_fact) {mom_fact = other_mom_fact;}
    }
    return jet->NN_dist * mom_fact;
  }

  /// establish the nearest neighbour for jet, and cross check consistency
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

  /// a table containing all the (full) distances to each object's NN
  double * diJ;

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
      NN_dist = BJ::geometrical_beam_distance();
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
template<class BJ, class I> void NNFJN2Plain<BJ,I>::start(const std::vector<PseudoJet> & jets) {
  n = jets.size();
  briefjets = new NNBJ[n];
  where_is.resize(2*n);

  NNBJ * jetA = briefjets;
  
  // initialise the basic jet info 
  for (int i = 0; i< n; i++) {
    // the this-> in the next line is required by standard compiler
    // see e.g. http://stackoverflow.com/questions/10639053/name-lookups-in-c-templates
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

  // now create the diJ (where J is i's NN) table -- remember that 
  // we differ from standard normalisation here by a factor of R2
  diJ = new double[n];
  jetA = head;
  for (int i = 0; i < n; i++) {
    diJ[i] = compute_diJ(jetA);
    jetA++; // have jetA follow i
  }
}


//----------------------------------------------------------------------
template<class BJ, class I> double NNFJN2Plain<BJ,I>::dij_min(int & iA, int & iB) {
  // find the minimum of the diJ on this round
  double diJ_min = diJ[0];
  int diJ_min_jet = 0;
  for (int i = 1; i < n; i++) {
    if (diJ[i] < diJ_min) {
      diJ_min_jet = i;
      diJ_min  = diJ[i];
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
template<class BJ, class I> void NNFJN2Plain<BJ,I>::remove_jet(int iA) {
  NNBJ * jetA = where_is[iA];
  // now update our nearest neighbour info and diJ table
  // first reduce size of table
  tail--; n--;
  // Copy last jet contents and diJ info into position of jetA
  *jetA = *tail;
  // update the info on where the given index is stored
  where_is[jetA->index()] = jetA;
  diJ[jetA - head] = diJ[tail-head];

  // updating NN infos. 2 cases to watch for (see below)
  for (NNBJ * jetI = head; jetI != tail; jetI++) {
    // see if jetI had jetA as a NN -- if so recalculate the NN
    if (jetI->NN == jetA) {
      set_NN_nocross(jetI, head, tail);
      diJ[jetI-head] = compute_diJ(jetI); // update diJ 
    } 
    // if jetI's NN is the new tail then relabel it so that it becomes jetA
    if (jetI->NN == tail) {jetI->NN = jetA;}
  }
}


//----------------------------------------------------------------------
template<class BJ, class I> void NNFJN2Plain<BJ,I>::merge_jets(int iA, int iB, 
					const PseudoJet & jet, int index) {

  NNBJ * jetA = where_is[iA];
  NNBJ * jetB = where_is[iB];

  // If necessary relabel A & B to ensure jetB < jetA, that way if
  // the larger of them == newtail then that ends up being jetA and 
  // the new jet that is added as jetB is inserted in a position that
  // has a future!
  if (jetA < jetB) std::swap(jetA,jetB);

  // initialise jetB based on the new jet
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
  diJ[jetA - head] = diJ[tail-head];

  // initialise jetB NN's distance and update NN infos
  for (NNBJ * jetI = head; jetI != tail; jetI++) {
    // see if jetI had jetA or jetB as a NN -- if so recalculate the NN
    if (jetI->NN == jetA || jetI->NN == jetB) {
      set_NN_nocross(jetI, head, tail);
      diJ[jetI-head] = compute_diJ(jetI); // update diJ 
    } 

    // check whether new jetB is closer than jetI's current NN and
    // if need be update things
    double dist = jetI->geometrical_distance(jetB);
    if (dist < jetI->NN_dist) {  // we need to update I
      if (jetI != jetB) {
        jetI->NN_dist = dist;
        jetI->NN = jetB;
        diJ[jetI-head] = compute_diJ(jetI); // update diJ...
      }
    }
    if (dist < jetB->NN_dist) { // we need to update B
      if (jetI != jetB) {
        jetB->NN_dist = dist;
        jetB->NN      = jetI;
      }
    }

    // if jetI's NN is the new tail then relabel it so that it becomes jetA
    if (jetI->NN == tail) {jetI->NN = jetA;}
  }

  // update info for jetB
  diJ[jetB-head] = compute_diJ(jetB);
}


//----------------------------------------------------------------------
// this function assumes that jet is not contained within begin...end
template <class BJ, class I> void NNFJN2Plain<BJ,I>::set_NN_crosscheck(NNBJ * jet, 
		    NNBJ * begin, NNBJ * end) {
  double NN_dist = jet->geometrical_beam_distance();
  NNBJ * NN      = NULL;
  for (NNBJ * jetB = begin; jetB != end; jetB++) {
    double dist = jet->geometrical_distance(jetB);
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
template <class BJ, class I>  void NNFJN2Plain<BJ,I>::set_NN_nocross(
                 NNBJ * jet, NNBJ * begin, NNBJ * end) {
  double NN_dist = jet->geometrical_beam_distance();
  NNBJ * NN      = NULL;
  // if (head < jet) {
  //   for (NNBJ * jetB = head; jetB != jet; jetB++) {
  if (begin < jet) {
    for (NNBJ * jetB = begin; jetB != jet; jetB++) {
      double dist = jet->geometrical_distance(jetB);
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
      double dist = jet->geometrical_distance(jetB);
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


#endif // __FASTJET_NNFJN2PLAIN_HH__
