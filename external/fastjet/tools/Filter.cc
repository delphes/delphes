//STARTHEADER
// $Id$
//
// Copyright (c) 2005-2011, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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

#include "fastjet/tools/Filter.hh"
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
  // start by getting the list of subjets (including a list of sanity
  // checks)
  // NB: subjets is empty to begin with (see the comment for
  //     _set_filtered_elements_cafilt)
  vector<PseudoJet> subjets; 
  JetDefinition subjet_def;
  bool discard_area;
  // NB: on return, subjet_def is set to the jet definition actually
  //     used (so that we can make use of its recombination scheme 
  //     when joining the jets to be kept).
  _set_filtered_elements(jet, subjets, subjet_def, discard_area);

  // now build the vector of kept and rejected subjets
  vector<PseudoJet> kept, rejected;
  // Note that in the following line we make a copy of the _selector
  // to avoid issues with needing a mutable _selector
  Selector selector_copy = _selector;
  if (selector_copy.takes_reference()) selector_copy.set_reference(jet);
  selector_copy.sift(subjets, kept, rejected);

  // gather the info under the form of a PseudoJet
  return _finalise(jet, kept, rejected, subjet_def, discard_area);
}


// sets filtered_elements to be all the subjets on which filtering will work
void Filter::_set_filtered_elements(const PseudoJet & jet,
                                    vector<PseudoJet> & filtered_elements,
                                    JetDefinition & subjet_def,
                                    bool & discard_area) const {
  // sanity checks
  //-------------------------------------------------------------------
  // make sure that the jet has constituents
  if (! jet.has_constituents())
    throw Error("Filter can only be applied on jets having constituents");

  // for a whole variety of checks, we shall need the "recursive"
  // pieces of the jet (the jet itself or recursing down to its most
  // fundamental pieces).
  // So we do compute these once and for all
  vector<PseudoJet> all_pieces; //.clear();
  if ((!_get_all_pieces(jet, all_pieces)) || (all_pieces.size()==0))
    throw Error("Attempt to filter a jet that has no associated ClusterSequence or is not a superposition of jets associated with a ClusterSequence");
  
  // if the filter uses subtraction, make sure we have a CS that supports area and has
  // explicit ghosts 
  if (_uses_subtraction()) {
    if (!jet.has_area())   
      throw Error("Attempt to filter and subtract (non-zero rho or subtractor) without area info for the original jet");

    if (!_check_explicit_ghosts(all_pieces))
      throw Error("Attempt to filter and subtract (non-zero rho or subtractor) without explicit ghosts");
  }

  // if we're dealing with a dynamic determination of the filtering
  // radius, do it now
  if ((_Rfilt>=0) || (_Rfiltfunc)){
    double Rfilt = (_Rfiltfunc) ? (*_Rfiltfunc)(jet) : _Rfilt;
    const JetDefinition::Recombiner * common_recombiner = _get_common_recombiner(all_pieces);
    if (common_recombiner) {
      if (typeid(*common_recombiner) == typeid(JetDefinition::DefaultRecombiner)) {
        RecombinationScheme scheme = 
          static_cast<const JetDefinition::DefaultRecombiner *>(common_recombiner)->scheme();
        subjet_def = JetDefinition(cambridge_algorithm, Rfilt, scheme);
      } else {
        subjet_def = JetDefinition(cambridge_algorithm, Rfilt, common_recombiner);
      }
    } else {
      subjet_def = JetDefinition(cambridge_algorithm, Rfilt);
    }
  } else {
    subjet_def = _subjet_def;
  }

  // get the jet definition to be use and whether we can apply our
  // simplified C/A+C/A filter
  //
  // we apply C/A clustering iff
  //  - the request subjet_def is C/A
  //  - the jet is either directly coming from C/A or if it is a
  //    superposition of C/A jets
  //  - the pieces agree with the recombination scheme of subjet_def
  //------------------------------------------------------------------
  bool simple_cafilt = _check_ca(all_pieces);

  // extract the subjets
  //-------------------------------------------------------------------
  discard_area = false;
  if (simple_cafilt){
    // first make sure that 'filtered_elements' is empty
    filtered_elements.clear();
    _set_filtered_elements_cafilt(jet, filtered_elements, subjet_def.R());
    // in the following case, areas can be erroneous and will be discarded
    discard_area = (!_uses_subtraction()) && (jet.has_area()) && (!_check_explicit_ghosts(all_pieces));
  } else {
    // here we'll simply recluster the jets.
    //
    // To include an area support we need
    //  - the jet to have an area
    //  - subtraction requested or explicit ghosts
    bool do_areas = (jet.has_area()) && ((_uses_subtraction()) || (_check_explicit_ghosts(all_pieces))); 
    _set_filtered_elements_generic(jet, filtered_elements, subjet_def, do_areas);
  }

  // order the filtered elements in pt
  filtered_elements = sorted_by_pt(filtered_elements);
}

// set the filtered elements in the simple case of C/A+C/A
//
// WATCH OUT: this could be recursively called, so filtered elements
//            of 'jet' are APPENDED to 'filtered_elements'
void Filter::_set_filtered_elements_cafilt(const PseudoJet & jet, 
                                           vector<PseudoJet> & filtered_elements, 
                                           double Rfilt) const{
  // we know that the jet is either a C/A jet or a superposition of
  // such pieces
  if (jet.has_associated_cluster_sequence()){
    // just extract the exclusive subjets of 'jet'
    const ClusterSequence *cs = jet.associated_cluster_sequence(); 
    vector<PseudoJet> local_fe;

    double dcut = Rfilt / cs->jet_def().R();
    if (dcut>=1.0){
      local_fe.push_back(jet);
    } else {
      dcut *= dcut;
      local_fe = jet.exclusive_subjets(dcut);
    }

    // subtract the jets if needed
    // Note that this one would work on pieces!!
    //-----------------------------------------------------------------
    if (_uses_subtraction()){
      const ClusterSequenceAreaBase * csab = jet.validated_csab();
      for (unsigned int i=0;i<local_fe.size();i++) {
        if (_subtractor) {
          local_fe[i] = (*_subtractor)(local_fe[i]);
        } else {
          local_fe[i] = csab->subtracted_jet(local_fe[i], _rho);
        }
      }
    }

    copy(local_fe.begin(), local_fe.end(), back_inserter(filtered_elements));
    return;
  }

  // just recurse into the pieces
  const vector<PseudoJet> & pieces = jet.pieces();
  for (vector<PseudoJet>::const_iterator it = pieces.begin(); 
       it!=pieces.end(); it++)
    _set_filtered_elements_cafilt(*it, filtered_elements, Rfilt);
}


// set the filtered elements in the generic re-clustering case (wo
// subtraction)
void Filter::_set_filtered_elements_generic(const PseudoJet & jet, 
                                            vector<PseudoJet> & filtered_elements,
                                            const JetDefinition & subjet_def,
					    bool do_areas) const{
  // create a new, internal, ClusterSequence from the jet constituents
  // get the subjets directly from there
  //
  // If the jet has area support then we separate the ghosts from the
  // "regular" particles so the subjets will also have area
  // support. Note that we do this regardless of whether rho is zero
  // or not.
  //
  // Note that to be able to separate the ghosts, one needs explicit
  // ghosts!!
  // ---------------------------------------------------------------
  if (do_areas){
    vector<PseudoJet> all_constituents = jet.constituents();
    vector<PseudoJet> regular_constituents, ghosts;  

    for (vector<PseudoJet>::iterator it = all_constituents.begin(); 
         it != all_constituents.end(); it++){
      if (it->is_pure_ghost())
        ghosts.push_back(*it);
      else
        regular_constituents.push_back(*it);
    }

    // figure out the ghost area from the 1st ghost (if none, any value
    // would probably do as the area will be 0 and subtraction will have
    // no effect!)
    double ghost_area = (ghosts.size()) ? ghosts[0].area() : 0.01;
    ClusterSequenceActiveAreaExplicitGhosts * csa
      = new ClusterSequenceActiveAreaExplicitGhosts(regular_constituents, 
                                                    subjet_def, 
                                                    ghosts, ghost_area);

    // get the subjets: we use the subtracted or unsubtracted ones
    // depending on rho or _subtractor being non-zero
    if (_uses_subtraction()) {
      if (_subtractor) {
        filtered_elements = (*_subtractor)(csa->inclusive_jets());
      } else {
        filtered_elements = csa->subtracted_jets(_rho);
      }
    } else {
      filtered_elements = csa->inclusive_jets();
    }

    // allow the cs to be deleted when it's no longer used
    csa->delete_self_when_unused();
  } else {
    ClusterSequence * cs = new ClusterSequence(jet.constituents(), subjet_def);
    filtered_elements = cs->inclusive_jets();
    // allow the cs to be deleted when it's no longer used
    cs->delete_self_when_unused();
  }
}


// gather the information about what is kept and rejected under the
// form of a PseudoJet with a special ClusterSequenceInfo
PseudoJet Filter::_finalise(const PseudoJet & /*jet*/, 
                            vector<PseudoJet> & kept, 
                            vector<PseudoJet> & rejected,
                            const JetDefinition & subjet_def,
                            const bool discard_area) const {
  // figure out which recombiner to use
  const JetDefinition::Recombiner &rec = *(subjet_def.recombiner());

  // create an appropriate structure and transfer the info to it
  PseudoJet filtered_jet = join<StructureType>(kept, rec);
  StructureType *fs = (StructureType*) filtered_jet.structure_non_const_ptr();
//  fs->_original_jet = jet;
  fs->_rejected = rejected;

  if (discard_area){
    // safety check: make sure there is an area to discard!!!
    assert(fs->_area_4vector_ptr);
    delete fs->_area_4vector_ptr;
    fs->_area_4vector_ptr=0;
  }
  
  return filtered_jet;
}

// various checks
//----------------------------------------------------------------------

// get the pieces down to the fundamental pieces
// 
// Note that this just checks that there is an associated CS to the
// fundamental pieces, not that it is still valid
bool Filter::_get_all_pieces(const PseudoJet &jet, vector<PseudoJet> &all_pieces) const{
  if (jet.has_associated_cluster_sequence()){
    all_pieces.push_back(jet);
    return true;
  }

  if (jet.has_pieces()){
    const vector<PseudoJet> pieces = jet.pieces();
    for (vector<PseudoJet>::const_iterator it=pieces.begin(); it!=pieces.end(); it++)
      if (!_get_all_pieces(*it, all_pieces)) return false;
    return true;
  }

  return false;
}

// get the common recombiner to all pieces (NULL if none)
//
// Note that if the jet has an associated cluster sequence that is no
// longer valid, an error will be thrown (needed since it could be the
// 1st check called after the enumeration of the pieces)
const JetDefinition::Recombiner* Filter::_get_common_recombiner(const vector<PseudoJet> &all_pieces) const{
  const JetDefinition & jd_ref = all_pieces[0].validated_cs()->jet_def();
  for (unsigned int i=1; i<all_pieces.size(); i++)
    if (!all_pieces[i].validated_cs()->jet_def().has_same_recombiner(jd_ref)) return NULL;

  return jd_ref.recombiner();
}

// check if the jet (or all its pieces) have explicit ghosts
// (assuming the jet has area support
//
// Note that if the jet has an associated cluster sequence that is no
// longer valid, an error will be thrown (needed since it could be the
// 1st check called after the enumeration of the pieces)
bool Filter::_check_explicit_ghosts(const vector<PseudoJet> &all_pieces) const{
  for (vector<PseudoJet>::const_iterator it=all_pieces.begin(); it!=all_pieces.end(); it++)
    if (! it->validated_csab()->has_explicit_ghosts()) return false;
  return true;
}

// check if one can apply the simplification for C/A subjets
//
// This includes:
//  - the subjet definition asks for C/A subjets
//  - all the pieces share the same CS
//  - that CS is C/A with the same recombiner as the subjet def
//  - the filtering radius is not larger than any of the pairwise
//    distance between the pieces
//
// Note that if the jet has an associated cluster sequence that is no
// longer valid, an error will be thrown (needed since it could be the
// 1st check called after the enumeration of the pieces)
bool Filter::_check_ca(const vector<PseudoJet> &all_pieces) const{
  if (_subjet_def.jet_algorithm() != cambridge_algorithm) return false;

  // check that the 1st of all the pieces (we're sure there is at
  // least one) is coming from a C/A clustering. Then check that all
  // the following pieces share the same ClusterSequence
  const ClusterSequence * cs_ref = all_pieces[0].validated_cs();
  if (cs_ref->jet_def().jet_algorithm() != cambridge_algorithm) return false;
  for (unsigned int i=1; i<all_pieces.size(); i++)
    if (all_pieces[i].validated_cs() != cs_ref) return false;

  // check that the 1st peice has the same recombiner as the one used
  // for the subjet clustering
  // Note that since they share the same CS, checking the 2st one is enough
  if (!cs_ref->jet_def().has_same_recombiner(_subjet_def)) return false;

  // we also have to make sure that the filtering radius is not larger
  // than any of the inter-pieces distance
  double Rfilt2 = _subjet_def.R();
  Rfilt2 *= Rfilt2;
  for (unsigned int i=0; i<all_pieces.size()-1; i++){
    for (unsigned int j=i+1; j<all_pieces.size(); j++){
      if (all_pieces[i].squared_distance(all_pieces[j]) <  Rfilt2) return false;
    }
  }

  return true;
}

//----------------------------------------------------------------------
// FilterInterface implementation 
//----------------------------------------------------------------------


FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh
