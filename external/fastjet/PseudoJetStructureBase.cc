//FJSTARTHEADER
// $Id: PseudoJetStructureBase.cc 3433 2014-07-23 08:17:03Z salam $
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


#include "fastjet/PseudoJetStructureBase.hh"
#include "fastjet/Error.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#ifndef __FJCORE__
#include "fastjet/ClusterSequenceAreaBase.hh"
#endif  // __FJCORE__

using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

// PseudoJetStructureBase implementation
//
// Contains any information related to the clustering that should be
// directly accessible to PseudoJet.
//
// By default, this class implements basic access to the
// ClusterSequence related to a PseudoJet (like its constituents or
// its area). But it can be overloaded in order e.g. to give access
// to the jet substructure.
//
// Note that it accesses the underlying ClusterSequence through a
// ClusterSequenceWrapper object so it can check when the former goes
// out of scope.
//


//-------------------------------------------------------------
// Direct access to the associated ClusterSequence object.
//
// Get access to the associated ClusterSequence (if any)
//-------------------------------------------------------------

// get a (const) pointer to the parent ClusterSequence (NULL if
// inexistent)
const ClusterSequence* PseudoJetStructureBase::associated_cluster_sequence() const{
  return NULL;
}
  
// if the jet has a valid associated cluster sequence then return a
// pointer to it; otherwise throw an error
//
// by default, an Error is thrown
const ClusterSequence * PseudoJetStructureBase::validated_cs() const{
  throw Error("This PseudoJet structure is not associated with a valid ClusterSequence");
}

#ifndef __FJCORE__
// if the jet has valid area information then return a pointer to
// the associated ClusterSequenceAreaBase object; otherwise throw an error
//
// by default, an Error is thrown
const ClusterSequenceAreaBase * PseudoJetStructureBase::validated_csab() const{
  throw Error("This PseudoJet structure is not associated with a valid cluster sequence with area");
}
#endif


//-------------------------------------------------------------
// Methods for access to information about jet structure
//
// These allow access to jet constituents, and other jet
// subtructure information. They only work if the jet is associated
// with a ClusterSequence.
//-------------------------------------------------------------

// check if it has been recombined with another PseudoJet in which
// case, return its partner through the argument. Otherwise,
// 'partner' is set to 0.
//
// by default, an Error is thrown
bool PseudoJetStructureBase::has_partner(const PseudoJet & /*reference */, PseudoJet & /*partner*/) const{
  throw Error("This PseudoJet structure has no implementation for has_partner");
}

// check if it has been recombined with another PseudoJet in which
// case, return its child through the argument. Otherwise, 'child'
// is set to 0.
// 
// by default, an Error is thrown
bool PseudoJetStructureBase::has_child(const PseudoJet & /*reference*/, PseudoJet & /*child*/) const{
  throw Error("This PseudoJet structure has no implementation for has_child");
}

// check if it is the product of a recombination, in which case
// return the 2 parents through the 'parent1' and 'parent2'
// arguments. Otherwise, set these to 0.
//
// by default, an Error is thrown
bool PseudoJetStructureBase::has_parents(const PseudoJet & /*reference*/, PseudoJet &/*parent1*/, PseudoJet &/*parent2*/) const{
  throw Error("This PseudoJet structure has no implementation for has_parents");
}

// check if the reference PseudoJet is contained in the second one
// passed as argument.
//
// by default, an Error is thrown
bool PseudoJetStructureBase::object_in_jet(const PseudoJet & /*reference*/, const PseudoJet & /*jet*/) const{
  throw Error("This PseudoJet structure has no implementation for is_inside");
}

// retrieve the constituents. 
//
// by default, an Error is thrown
vector<PseudoJet> PseudoJetStructureBase::constituents(const PseudoJet &/*reference*/) const{
  throw Error("This PseudoJet structure has no implementation for constituents");
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
// by default, an Error is thrown
vector<PseudoJet> PseudoJetStructureBase::exclusive_subjets (const PseudoJet & /*reference*/, const double & /*dcut*/) const{
  throw Error("This PseudoJet structure has no implementation for exclusive_subjets");
}

// return the size of exclusive_subjets(...); still n ln n with same
// coefficient, but marginally more efficient than manually taking
// exclusive_subjets.size()
//
// by default, an Error is thrown
int PseudoJetStructureBase::n_exclusive_subjets(const PseudoJet & /*reference*/, const double & /*dcut*/) const{
  throw Error("This PseudoJet structure has no implementation for n_exclusive_subjets");
}

// return the list of subjets obtained by unclustering the supplied
// jet down to n subjets (or all constituents if there are fewer
// than n).
//
// by default, an Error is thrown
vector<PseudoJet> PseudoJetStructureBase::exclusive_subjets_up_to (const PseudoJet & /*reference*/, int /*nsub*/) const{
  throw Error("This PseudoJet structure has no implementation for exclusive_subjets");
}

// return the dij that was present in the merging nsub+1 -> nsub 
// subjets inside this jet.
//
// by default, an Error is thrown
double PseudoJetStructureBase::exclusive_subdmerge(const PseudoJet & /*reference*/, int /*nsub*/) const{
  throw Error("This PseudoJet structure has no implementation for exclusive_submerge");
}

// return the maximum dij that occurred in the whole event at the
// stage that the nsub+1 -> nsub merge of subjets occurred inside 
// this jet.
//
// by default, an Error is thrown
double PseudoJetStructureBase::exclusive_subdmerge_max(const PseudoJet & /*reference*/, int /*nsub*/) const{
  throw Error("This PseudoJet structure has no implementation for exclusive_submerge_max");
}


// retrieve the pieces building the jet. 
//
// by default, an Error is thrown
std::vector<PseudoJet> PseudoJetStructureBase::pieces(const PseudoJet & /*reference*/) const{
  throw Error("This PseudoJet structure has no implementation for pieces");  
}

// the following ones require a computation of the area in the
// parent ClusterSequence (See ClusterSequenceAreaBase for details)
//------------------------------------------------------------------
#ifndef __FJCORE__

// return the jet (scalar) area.
//
// by default, an Error is thrown
double PseudoJetStructureBase::area(const PseudoJet & /*reference*/) const{
  throw Error("This PseudoJet structure has no implementation for area");
}

// return the error (uncertainty) associated with the determination
// of the area of this jet.
//
// by default, an Error is thrown
double PseudoJetStructureBase::area_error(const PseudoJet & /*reference*/) const{
  throw Error("This PseudoJet structure has no implementation for area_error");
}

// return the jet 4-vector area.
//
// by default, an Error is thrown
PseudoJet PseudoJetStructureBase::area_4vector(const PseudoJet & /*reference*/) const{
  throw Error("This PseudoJet structure has no implementation for area_4vector");
}

// true if this jet is made exclusively of ghosts.
//
// by default, an Error is thrown
bool PseudoJetStructureBase::is_pure_ghost(const PseudoJet & /*reference*/) const{
  throw Error("This PseudoJet structure has no implementation for is_pure_ghost");
}
#endif  // __FJCORE__

FASTJET_END_NAMESPACE
