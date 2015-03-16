// $Id: Recluster.cc 699 2014-07-07 09:58:12Z gsoyez $
//
// Copyright (c) 2014-, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------

#include "Recluster.hh"
#include <fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh>
#include <sstream>
#include <typeinfo>

using namespace std;

// Comments:
//
//  - If the jet comes from a C/A clustering (or is a composite jet
//    made of C/A clusterings) and we recluster it with a C/A
//    algorithm, we just use exclusive jets instead of doing the
//    clustering explicitly. In this specific case, we need to make
//    sure that all the pieces share the same cluster sequence.
//
//  - If the recombiner has to be taken from the original jet and that
//    jet is composite, we need to check that all the pieces were
//    obtained with the same recombiner.
//
// TODO:
//
//  - check this actually works!!!

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

LimitedWarning Recluster::_explicit_ghost_warning;

// class description
string Recluster::description() const {
  ostringstream ostr;
  ostr << "Recluster with subjet_def = ";
  if (_use_full_def) {
    ostr << _subjet_def.description();
  } else {
    if (_subjet_alg == kt_algorithm) {
      ostr << "Longitudinally invariant kt algorithm with R = " << _subjet_radius;
    } else if (_subjet_alg == cambridge_algorithm) {
      ostr << "Longitudinally invariant Cambridge/Aachen algorithm with R = " << _subjet_radius;
    } else if (_subjet_alg == antikt_algorithm) {
      ostr << "Longitudinally invariant anti-kt algorithm with R = " << _subjet_radius;
    } else if (_subjet_alg == genkt_algorithm) {
      ostr << "Longitudinally invariant generalised kt algorithm with R = " << _subjet_radius
           << ", p = " << _subjet_extra;
    } else if (_subjet_alg == cambridge_for_passive_algorithm) {
      ostr << "Longitudinally invariant Cambridge/Aachen algorithm with R = " << _subjet_radius
           << " and a special hack whereby particles with kt < " 
           << _subjet_extra << "are treated as passive ghosts";
    } else if (_subjet_alg == ee_kt_algorithm) {
      ostr << "e+e- kt (Durham) algorithm";
    } else if (_subjet_alg == ee_genkt_algorithm) {
      ostr << "e+e- generalised kt algorithm with R = " << _subjet_radius
           << ", p = " << _subjet_extra;
    } else if (_subjet_alg == undefined_jet_algorithm) {
      ostr << "uninitialised JetDefinition (jet_algorithm=undefined_jet_algorithm)" ;
    } else {
      ostr << "unrecognized jet_algorithm";
    }
    ostr << ", a recombiner obtained from the jet being reclustered";
  }

  if (_single)
    ostr << " and keeping the hardest subjet";
  else
    ostr << " and joining all subjets in a composite jet";

  return ostr.str();
}

// return a vector of subjets, which are the ones that would be kept
// by the filtering
PseudoJet Recluster::result(const PseudoJet &jet) const {
  // generic sanity checks
  //-------------------------------------------------------------------
  // make sure that the jet has constituents
  if (! jet.has_constituents())
    throw Error("Filter can only be applied on jets having constituents");

  // tests particular to certain configurations
  //-------------------------------------------------------------------

  // for a whole variety of checks, we shall need the "recursive"
  // pieces of the jet (the jet itself or recursing down to its most
  // fundamental pieces). So we do compute these once and for all.
  //
  // Note that the pieces are always needed (either for C/A or for the
  // area checks)
  vector<PseudoJet> all_pieces; //.clear();
  if ((!_get_all_pieces(jet, all_pieces)) || (all_pieces.size()==0)){
    throw Error("Recluster: failed to retrieve all the pieces composing the jet.");
  }

  // decide which jet definition to use
  //-------------------------------------------------------------------
  JetDefinition subjet_def;
  if (_use_full_def){
    subjet_def = _subjet_def;
  } else {
    _build_jet_def_with_recombiner(all_pieces, subjet_def);
  }


  // the vector that will ultimately hold the subjets
  vector<PseudoJet> subjets;

  // check if we can apply the simplification for C/A jets reclustered
  // with C/A
  //
  // we apply C/A clustering iff
  //  - the request subjet_def is C/A
  //  - the jet is either directly coming from C/A or if it is a
  //    superposition of C/A jets from the same cluster sequence
  //  - the pieces agree with the recombination scheme of subjet_def
  //
  // Note that in this case area support will be automatically
  // inherted so we can only worry about this later
  //-------------------------------------------------------------------
  if (_check_ca(all_pieces, subjet_def)){
    _recluster_cafilt(all_pieces, subjets, subjet_def.R());
    subjets = sorted_by_pt(subjets);
    return _single
      ? subjets[0]
      : join(subjets, *(subjet_def.recombiner()));
  }

  // decide if area support has to be kept
  //-------------------------------------------------------------------
  bool include_area_support = jet.has_area();
  if ((include_area_support) &&  (!_check_explicit_ghosts(all_pieces))){
    _explicit_ghost_warning.warn("Recluster: the original cluster sequence is lacking explicit ghosts; area support will no longer be available after re-clustering");
    include_area_support = false;
  }

  // extract the subjets
  //-------------------------------------------------------------------
  _recluster_generic(jet, subjets, subjet_def, include_area_support);
  subjets = sorted_by_pt(subjets);

  return _single
    ? subjets[0]
    : join(subjets, *(subjet_def.recombiner()));
}

//----------------------------------------------------------------------
// the parts that really do the reclustering
//----------------------------------------------------------------------

// get the subjets in the simple case of C/A+C/A
void Recluster::_recluster_cafilt(const vector<PseudoJet> & all_pieces, 
                                  vector<PseudoJet> & subjets, 
                                  double Rfilt) const{
  subjets.clear();

  // each individual piece should have a C/A cluster sequence
  // associated to it. So a simple loop would do the job
  for (vector<PseudoJet>::const_iterator piece_it = all_pieces.begin(); 
       piece_it!=all_pieces.end(); piece_it++){
    // just extract the exclusive subjets of 'jet'
    const ClusterSequence *cs = piece_it->associated_cluster_sequence(); 
    vector<PseudoJet> local_subjets;

    double dcut = Rfilt / cs->jet_def().R();
    if (dcut>=1.0){
      local_subjets.push_back(*piece_it);
    } else {
      local_subjets = piece_it->exclusive_subjets(dcut*dcut);
    }

    copy(local_subjets.begin(), local_subjets.end(), back_inserter(subjets));
  }
}


// set the filtered elements in the generic re-clustering case (w/o
// subtraction)
void Recluster::_recluster_generic(const PseudoJet & jet, 
                                   vector<PseudoJet> & subjets,
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

    subjets = csa->inclusive_jets();

    // allow the cs to be deleted when it's no longer used
    // 
    // Note that there is at least one constituent in the jet so there
    // is in principle at least one subjet But one may have used a
    // nasty recombiner that left an empty set of subjets, so we'd
    // rather play it safe
    if (subjets.size())
      csa->delete_self_when_unused();
    else
      delete csa;
  } else {
    ClusterSequence * cs = new ClusterSequence(jet.constituents(), subjet_def);
    subjets = cs->inclusive_jets();
    // allow the cs to be deleted when it's no longer used (again, we
    // add an extra safety check)
    if (subjets.size())
      cs->delete_self_when_unused();
    else 
      delete cs;
  }
}


//----------------------------------------------------------------------
// various checks and internal constructs
//----------------------------------------------------------------------

// fundamental info for CompositeJets
//----------------------------------------------------------------------
 
// get the pieces down to the fundamental pieces
// 
// Note that this just checks that there is an associated CS to the
// fundamental pieces, not that it is still valid
bool Recluster::_get_all_pieces(const PseudoJet &jet, vector<PseudoJet> &all_pieces) const{
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

// treatment of recombiners
//----------------------------------------------------------------------
// get the common recombiner to all pieces (NULL if none)
//
// Note that if the jet has an associated cluster sequence that is no
// longer valid, an error will be thrown (needed since it could be the
// 1st check called after the enumeration of the pieces)
const JetDefinition::Recombiner* Recluster::_get_common_recombiner(const vector<PseudoJet> &all_pieces) const{
  const JetDefinition & jd_ref = all_pieces[0].validated_cs()->jet_def();
  for (unsigned int i=1; i<all_pieces.size(); i++)
    if (!all_pieces[i].validated_cs()->jet_def().has_same_recombiner(jd_ref)) return NULL;

  return jd_ref.recombiner();
}
  
void Recluster::_build_jet_def_with_recombiner(const vector<PseudoJet> &all_pieces, 
                                               JetDefinition &subjet_def) const{
  // the recombiner has to be guessed from the pieces
  const JetDefinition::Recombiner * common_recombiner = _get_common_recombiner(all_pieces);
  if (common_recombiner) {
    if (typeid(*common_recombiner) == typeid(JetDefinition::DefaultRecombiner)) {
      RecombinationScheme scheme = 
        static_cast<const JetDefinition::DefaultRecombiner *>(common_recombiner)->scheme();
      if (_has_subjet_extra)
        subjet_def = JetDefinition(_subjet_alg, _subjet_radius, _subjet_extra, scheme);
      else if (_has_subjet_radius)
        subjet_def = JetDefinition(_subjet_alg, _subjet_radius, scheme);
      else 
        subjet_def = JetDefinition(_subjet_alg, scheme);
    } else {
      if (_has_subjet_extra)
        subjet_def = JetDefinition(_subjet_alg, _subjet_radius, _subjet_extra, common_recombiner);
      else if (_has_subjet_radius)
        subjet_def = JetDefinition(_subjet_alg, _subjet_radius, common_recombiner);
      else 
        subjet_def = JetDefinition(_subjet_alg, common_recombiner);
    }
  } else {
    throw Error("Recluster: requested to guess the recombination scheme (or recombiner) from the original jet but an inconsistency was found between the pieces constituing that jet.");
  }
}

// area support
//----------------------------------------------------------------------

// check if the jet (or all its pieces) have explicit ghosts
// (assuming the jet has area support).
//
// Note that if the jet has an associated cluster sequence that is no
// longer valid, an error will be thrown (needed since it could be the
// 1st check called after the enumeration of the pieces)
bool Recluster::_check_explicit_ghosts(const vector<PseudoJet> &all_pieces) const{
  for (vector<PseudoJet>::const_iterator it=all_pieces.begin(); it!=all_pieces.end(); it++)
    if (! it->validated_csab()->has_explicit_ghosts()) return false;
  return true;
}

// C/A specific tests
//----------------------------------------------------------------------

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
bool Recluster::_check_ca(const vector<PseudoJet> &all_pieces, 
                          const JetDefinition &subjet_def) const{
  if (subjet_def.jet_algorithm() != cambridge_algorithm) return false;

  // check that the 1st of all the pieces (we're sure there is at
  // least one) is coming from a C/A clustering. Then check that all
  // the following pieces share the same ClusterSequence
  const ClusterSequence * cs_ref = all_pieces[0].validated_cs();
  if (cs_ref->jet_def().jet_algorithm() != cambridge_algorithm) return false;
  for (unsigned int i=1; i<all_pieces.size(); i++)
    if (all_pieces[i].validated_cs() != cs_ref) return false;

  // check that the 1st peice has the same recombiner as the one used
  // for the subjet clustering
  // Note that since they share the same CS, checking the 1st one is enough
  if (!cs_ref->jet_def().has_same_recombiner(subjet_def)) return false;

  // we also have to make sure that the reclustering radius is not larger
  // than any of the inter-pieces distance
  double Rsub2 = subjet_def.R();
  Rsub2 *= Rsub2;
  for (unsigned int i=0; i<all_pieces.size()-1; i++){
    for (unsigned int j=i+1; j<all_pieces.size(); j++){
      if (all_pieces[i].squared_distance(all_pieces[j]) <  Rsub2) return false;
    }
  }

  return true;
}

} // contrib namespace

FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh
