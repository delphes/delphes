//FJSTARTHEADER
// $Id: ClusterSequenceStructure.cc 3433 2014-07-23 08:17:03Z salam $
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

#include "fastjet/ClusterSequenceStructure.hh"
#include "fastjet/Error.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#ifndef __FJCORE__
#include "fastjet/ClusterSequenceAreaBase.hh"
#endif  // __FJCORE__
#include <iostream>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

using namespace std;

ClusterSequenceStructure::~ClusterSequenceStructure(){
  if (_associated_cs != NULL 
      && _associated_cs->will_delete_self_when_unused()) {
    // automatically handle deletion of the cluster sequence;
    // execution should only ever reach this point if the user had
    // called CS::delete_self_when_unused, which resets the count of
    // the shared pointer to CSS (otherwise the CS's own destructor
    // will have zeroed the _associated_cs pointer before the shared
    // pointer count goes to zero [on destruction of the last of the
    // jets in the CS and the destruction of the CS's copy of the
    // shared pointer)
    _associated_cs->signal_imminent_self_deletion();
    delete _associated_cs;
  }
}


//----------------------------------------------------------------------
// Direct access to the associated ClusterSequence object.
//----------------------------------------------------------------------

// check whether this PseudoJet has an associated parent
// ClusterSequence
bool ClusterSequenceStructure::has_valid_cluster_sequence() const{
  return (_associated_cs != NULL);
}

// get a (const) pointer to the associated ClusterSequence (NULL if
// inexistent)
const ClusterSequence* ClusterSequenceStructure::associated_cluster_sequence() const{
  return _associated_cs;
}


// If there is a valid cluster sequence associated with this jet,
// returns a pointer to it; otherwise throws an Error.
//
// Open question: should these errors be upgraded to classes of their
// own so that they can be caught? [Maybe, but later]
const ClusterSequence * ClusterSequenceStructure::validated_cs() const {
  if (!_associated_cs) 
    throw Error("you requested information about the internal structure of a jet, but its associated ClusterSequence has gone out of scope.");
  return _associated_cs;
}


//----------------------------------------------------------------------
// Methods for access to information about jet structure
//----------------------------------------------------------------------

// check if it has been recombined with another PseudoJet in which
// case, return its partner through the argument. Otherwise,
// 'partner' is set to 0.
//
// false is also returned if this PseudoJet has no associated
// ClusterSequence
bool ClusterSequenceStructure::has_partner(const PseudoJet &reference, PseudoJet &partner) const{
  return validated_cs()->has_partner(reference, partner);
}

// check if it has been recombined with another PseudoJet in which
// case, return its child through the argument. Otherwise, 'child'
// is set to 0.
// 
// false is also returned if this PseudoJet has no associated
// ClusterSequence, with the child set to 0
bool ClusterSequenceStructure::has_child(const PseudoJet &reference, PseudoJet &child) const{
  return validated_cs()->has_child(reference, child);
}

// check if it is the product of a recombination, in which case
// return the 2 parents through the 'parent1' and 'parent2'
// arguments. Otherwise, set these to 0.
//
// false is also returned if this PseudoJet has no parent
// ClusterSequence
bool ClusterSequenceStructure::has_parents(const PseudoJet &reference, PseudoJet &parent1, PseudoJet &parent2) const{
  return validated_cs()->has_parents(reference, parent1, parent2);
}


// check if the reference PseudoJet is inside the "jet" passed as an argument
//
// an error is thrown if there is no CS associated with one of the 2 jets.
// fasle is returned if teh 2 jets do not belong to the same CS
bool ClusterSequenceStructure::object_in_jet(const PseudoJet &reference, const PseudoJet &jet) const{
  if ((!has_associated_cluster_sequence()) || (!jet.has_associated_cluster_sequence()))
    throw Error("you requested information about the internal structure of a jet, but it is not associated with a ClusterSequence or its associated ClusterSequence has gone out of scope."); 

  if (reference.associated_cluster_sequence() != jet.associated_cluster_sequence()) return false;

  return validated_cs()->object_in_jet(reference, jet);
}


// return true if the structure supports constituents. 
//
// an Error is thrown if this PseudoJet has no currently valid
// associated ClusterSequence
bool ClusterSequenceStructure::has_constituents() const{
  if (!has_associated_cluster_sequence())
    throw Error("you requested information about the internal structure of a jet, but it is not associated with a ClusterSequence or its associated ClusterSequence has gone out of scope."); 

  return true;
}


// retrieve the constituents. An empty vector is returned if there is
// no associated ClusterSequence
vector<PseudoJet> ClusterSequenceStructure::constituents(const PseudoJet &reference) const{
  return validated_cs()->constituents(reference);
}

// return true if the structure supports exclusive_subjets. 
//
// an Error is thrown if this PseudoJet has no currently valid
// associated ClusterSequence
bool ClusterSequenceStructure::has_exclusive_subjets() const{
  if (!has_associated_cluster_sequence())
    throw Error("you requested information about the internal structure of a jet, but it is not associated with a ClusterSequence or its associated ClusterSequence has gone out of scope."); 

  return true;
}

// return a vector of all subjets of the current jet (in the sense
// of the exclusive algorithm) that would be obtained when running
// the algorithm with the given dcut. 
//
// Time taken is O(m ln m), where m is the number of subjets that
// are found. If m gets to be of order of the total number of
// constituents in the jet, this could be substantially slower than
// just getting that list of constituents.
//
// an Error is thrown if this PseudoJet has no currently valid
// associated ClusterSequence
std::vector<PseudoJet> ClusterSequenceStructure::exclusive_subjets (const PseudoJet &reference, const double & dcut) const {
  return validated_cs()->exclusive_subjets(reference, dcut);
}

// return the size of exclusive_subjets(...); still n ln n with same
// coefficient, but marginally more efficient than manually taking
// exclusive_subjets.size()
//
// an Error is thrown if this PseudoJet has no currently valid
// associated ClusterSequence
int ClusterSequenceStructure::n_exclusive_subjets(const PseudoJet &reference, const double & dcut) const {
  return validated_cs()->n_exclusive_subjets(reference, dcut);
}

// return the list of subjets obtained by unclustering the supplied
// jet down to n subjets (or all constituents if there are fewer
// than n).
//
// requires n ln n time
//
// an Error is thrown if this PseudoJet has no currently valid
// associated ClusterSequence
std::vector<PseudoJet> ClusterSequenceStructure::exclusive_subjets_up_to (const PseudoJet &reference, int nsub) const {
  return validated_cs()->exclusive_subjets_up_to(reference, nsub);
}

// return the dij that was present in the merging nsub+1 -> nsub 
// subjets inside this jet.
//
// an Error is thrown if this PseudoJet has no currently valid
// associated ClusterSequence
double ClusterSequenceStructure::exclusive_subdmerge(const PseudoJet &reference, int nsub) const {
  return validated_cs()->exclusive_subdmerge(reference, nsub);
}

// return the maximum dij that occurred in the whole event at the
// stage that the nsub+1 -> nsub merge of subjets occurred inside 
// this jet.
//
// an Error is thrown if this PseudoJet has no currently valid
// associated ClusterSequence
double ClusterSequenceStructure::exclusive_subdmerge_max(const PseudoJet &reference, int nsub) const {
  return validated_cs()->exclusive_subdmerge_max(reference, nsub);
}


//----------------------------------------------------------------------
// information related to the pieces of the jet
//----------------------------------------------------------------------

// by convention, a jet associated with a ClusterSequence will have
// pieces if it has parents in the cluster sequence.
//
// an error is thrown if the ClusterSequence is out of scope (since
// the answer depends on information in the Cluster Sequence)
bool ClusterSequenceStructure::has_pieces(const PseudoJet &reference) const{
  PseudoJet dummy1, dummy2;
  return has_parents(reference, dummy1, dummy2);
}

// by convention, the pieces of a jet associated with a
// ClusterSequence are its parents in the Cluster Sequence. If it has
// no parents, an empty jet is returned.
//
// an error is thrown if the ClusterSequence is out of scope
vector<PseudoJet> ClusterSequenceStructure::pieces(const PseudoJet &reference) const{
  PseudoJet j1, j2;
  vector<PseudoJet> res;
  if (has_parents(reference, j1, j2)){
    res.push_back(j1);
    res.push_back(j2);
  }

  return res;
}


//----------------------------------------------------------------------
// the following ones require a computation of the area in the
// associated ClusterSequence (See ClusterSequenceAreaBase for details)
//----------------------------------------------------------------------

#ifndef __FJCORE__
// if possible, return a valid ClusterSequenceAreaBase pointer; otherwise
// throw an error
const ClusterSequenceAreaBase * ClusterSequenceStructure::validated_csab() const {
  const ClusterSequenceAreaBase *csab = dynamic_cast<const ClusterSequenceAreaBase*>(validated_cs());
  if (csab == NULL) throw Error("you requested jet-area related information, but the PseudoJet does not have associated area information.");
  return csab;
}

// check if it has a defined area
bool ClusterSequenceStructure::has_area() const{
  if (! has_associated_cluster_sequence()) return false;
  return (dynamic_cast<const ClusterSequenceAreaBase*>(_associated_cs) != NULL);
}

// return the jet (scalar) area.
// throw an Error if there is no support for area in the associated CS
double ClusterSequenceStructure::area(const PseudoJet &reference) const{
  return validated_csab()->area(reference);
}

// return the error (uncertainty) associated with the determination
// of the area of this jet.
// throws an Error if there is no support for area in the associated CS
double ClusterSequenceStructure::area_error(const PseudoJet &reference) const{
  return validated_csab()->area_error(reference);
}

// return the jet 4-vector area
// throws an Error if there is no support for area in the associated CS
PseudoJet ClusterSequenceStructure::area_4vector(const PseudoJet &reference) const{
  return validated_csab()->area_4vector(reference);
}

// true if this jet is made exclusively of ghosts
// throws an Error if there is no support for area in the associated CS
bool ClusterSequenceStructure::is_pure_ghost(const PseudoJet &reference) const{
  return validated_csab()->is_pure_ghost(reference);
}

#endif  // __FJCORE__



FASTJET_END_NAMESPACE
