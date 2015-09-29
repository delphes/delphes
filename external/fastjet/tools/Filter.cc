//FJSTARTHEADER
// $Id: Filter.cc 3760 2014-12-19 10:05:10Z soyez $
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

#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Recluster.hh"
#include <fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh>
#include <cassert>
#include <algorithm>
#include <sstream>
#include <typeinfo>

using namespace std;


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//----------------------------------------------------------------------
// Filter class implementation
//----------------------------------------------------------------------

// class description
string Filter::description() const {
  if (!_initialised){
    return "uninitialised Filter";
  }

  ostringstream ostr;
  ostr << "Filter with subjet_def = ";
  if (_Rfiltfunc) {
    ostr << "Cambridge/Aachen algorithm with dynamic Rfilt"
         << " (recomb. scheme deduced from jet, or E-scheme if not unique)";
  } else if (_Rfilt > 0) {
    ostr << "Cambridge/Aachen algorithm with Rfilt = " 
         << _Rfilt 
         << " (recomb. scheme deduced from jet, or E-scheme if not unique)";
  } else {
    ostr << _subjet_def.description();
  }
  ostr<< ", selection " << _selector.description();
  if (_subtractor) {
    ostr << ", subtractor: " << _subtractor->description();
  } else if (_rho != 0) {
    ostr << ", subtracting with rho = " << _rho;
  }
  return ostr.str();
}


// core functions
//----------------------------------------------------------------------

// return a vector of subjets, which are the ones that would be kept
// by the filtering
PseudoJet Filter::result(const PseudoJet &jet) const {
  if (!_initialised){
    //Q: do we throw or do we return an empty PJ?
    throw Error("uninitialised Filter");
  }

  // start by getting the list of subjets (including a list of sanity
  // checks)
  // NB: subjets is empty to begin with (see the comment for
  //     _set_filtered_elements_cafilt)
  vector<PseudoJet> subjets; 
  //JetDefinition subjet_def;
  bool ca_optimised = _set_filtered_elements(jet, subjets);

  // apply subtraction if needed:
  if (_subtractor){
    subjets = (*_subtractor)(subjets);
  } else if (_rho!=0){
    if (subjets.size()>0){
      const ClusterSequenceAreaBase *csab = subjets[0].validated_csab();
      for (unsigned int i=0;i<subjets.size();i++){
        subjets[i]=csab->subtracted_jet(subjets[i], _rho);
      }
    }
  }

  // now build the vector of kept and rejected subjets
  vector<PseudoJet> kept, rejected;
  // Note that in the following line we make a copy of the _selector
  // to avoid issues with needing a mutable _selector
  Selector selector_copy = _selector;
  if (selector_copy.takes_reference()) selector_copy.set_reference(jet);
  selector_copy.sift(subjets, kept, rejected);

  // gather the info under the form of a PseudoJet
  return _finalise(jet, kept, rejected, ca_optimised);
}


// sets filtered_elements to be all the subjets on which filtering will work
//
// return true when the subjets have been optained using teh optimised
// method for C/A
bool Filter::_set_filtered_elements(const PseudoJet & jet,
                                    vector<PseudoJet> & filtered_elements) const {
  // create the recluster instance
  Recluster recluster;
  if ((_Rfilt>=0) || (_Rfiltfunc))
    recluster = Recluster(cambridge_algorithm, (_Rfiltfunc) ? (*_Rfiltfunc)(jet) : _Rfilt, Recluster::keep_all);
  else
    recluster = Recluster(_subjet_def, false, Recluster::keep_all);

  // get the subjets
  //JetDefinition subjet_def;
  return recluster.get_new_jets_and_def(jet, filtered_elements);
}

// gather the information about what is kept and rejected under the
// form of a PseudoJet with a special ClusterSequenceInfo
PseudoJet Filter::_finalise(const PseudoJet & /*jet*/, 
                            vector<PseudoJet> & kept, 
                            vector<PseudoJet> & rejected,
			    bool ca_optimisation_used) const {
  PseudoJet filtered_jet;

  if (kept.size()+rejected.size()>0){
    // figure out which recombiner to use
    const JetDefinition::Recombiner &rec = (kept.size()>0)
      ? *(kept[0].associated_cs()->jet_def().recombiner())
      : *(rejected[0].associated_cs()->jet_def().recombiner());

    // create an appropriate structure and transfer the info to it
    filtered_jet = join<StructureType>(kept, rec);
  } else {
    filtered_jet = join<StructureType>(kept);
  }
  StructureType *fs = (StructureType*) filtered_jet.structure_non_const_ptr();
  fs->_rejected = rejected;
  
  // if we've used C/A optimisation, we need to get rid of the area
  // information if it comes from a non-explicit-ghost clustering.
  // (because in that case it can be erroneous due the lack of
  // information about empty areas)
  if ((ca_optimisation_used) && (kept.size()+rejected.size()>0)){
    bool has_non_explicit_ghost_area = (kept.size()>0)
      ? (kept[0].has_area()     && (!(kept[0].validated_csab()->has_explicit_ghosts())))
      : (rejected[0].has_area() && (!(rejected[0].validated_csab()->has_explicit_ghosts())));
    if (has_non_explicit_ghost_area)
      fs->discard_area();
  }

  return filtered_jet;
}


FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh
