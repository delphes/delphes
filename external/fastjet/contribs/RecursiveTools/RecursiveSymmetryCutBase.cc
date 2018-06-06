// $Id: RecursiveSymmetryCutBase.cc 1080 2017-09-28 07:51:37Z gsoyez $
//
// Copyright (c) 2014-, Gavin P. Salam, Gregory Soyez, Jesse Thaler
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

#include "RecursiveSymmetryCutBase.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequenceAreaBase.hh"
#include <cassert>
#include <algorithm> 
#include <cstdlib> 

using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

LimitedWarning RecursiveSymmetryCutBase::_negative_mass_warning;
LimitedWarning RecursiveSymmetryCutBase::_mu2_gt1_warning;
//LimitedWarning RecursiveSymmetryCutBase::_nonca_warning;
LimitedWarning RecursiveSymmetryCutBase::_explicit_ghost_warning;

bool RecursiveSymmetryCutBase::_verbose = false;

//----------------------------------------------------------------------
PseudoJet RecursiveSymmetryCutBase::result(const PseudoJet & jet) const {
  // construct the input jet (by default, recluster with C/A)
  if (! jet.has_constituents()){
    throw Error("RecursiveSymmetryCutBase can only be applied to jets with constituents");
  }

  PseudoJet j = _recluster_if_needed(jet);

  // sanity check: the jet must have a valid CS
  if (! j.has_valid_cluster_sequence()){
    throw Error("RecursiveSymmetryCutBase can only be applied to jets associated to a (valid) cluster sequence");
  }

  // check that area information is there in case we have a subtractor
  // GS: do we really need this since subtraction may not require areas?
  if (_subtractor) {
    const ClusterSequenceAreaBase * csab = 
      dynamic_cast<const ClusterSequenceAreaBase *>(j.associated_cs());
    if (csab == 0 || (!csab->has_explicit_ghosts()))
      _explicit_ghost_warning.warn("RecursiveSymmetryCutBase: there is no clustering sequence, or it lacks explicit ghosts: subtraction is not guaranteed to function properly");
  }

  // establish the first subjet and optionally subtract it
  PseudoJet subjet = j;
  if (_subtractor && (!_input_jet_is_subtracted)) {
    subjet = (*_subtractor)(subjet);
  }

  // variables for tracking what will happen
  PseudoJet piece1, piece2;
 
  // vectors for storing optional verbose structure
  // these hold the deltaR, symmetry, and mu values of dropped branches
  std::vector<double> dropped_delta_R;
  std::vector<double> dropped_symmetry;
  std::vector<double> dropped_mu;

  double sym, mu2;
  
  // now recurse into the jet's structure
  RecursionStatus status;
  while ((status=recurse_one_step(subjet, piece1, piece2, sym, mu2)) != recursion_success) {
    // start with sanity checks:
    if ((status == recursion_issue) || (status == recursion_no_parents)) {
      // we should return piece1 by our convention for recurse_one_step
      PseudoJet result;
      if (status == recursion_issue){
        result = piece1;
        if (_verbose) cout << "reached end; returning null jet " << endl;
      } else {
        result = _result_no_substructure(piece1);
        if (_verbose) cout << "no parents found; returning last PJ or empty jet" << endl;
      }
      
      if (result != 0) {
        // if in grooming mode, add dummy structure information
        StructureType * structure = new StructureType(result);
        // structure->_symmetry = 0.0;
        // structure->_mu       = 0.0;
        // structure->_delta_R  = 0.0;
        if (_verbose_structure) { // still want to store verbose information about dropped branches
          structure->_has_verbose = true;
          structure->_dropped_symmetry = dropped_symmetry;
          structure->_dropped_mu = dropped_mu;
          structure->_dropped_delta_R = dropped_delta_R;
        }
        result.set_structure_shared_ptr(SharedPtr<PseudoJetStructureBase>(structure));
      }
  
      return result;
    }

    assert(status == recursion_dropped);
    
    // if desired, store information about dropped branches before recursing
    if (_verbose_structure) {
      dropped_delta_R.push_back(piece1.delta_R(piece2));
      dropped_symmetry.push_back(sym);
      dropped_mu.push_back((mu2 >= 0) ? sqrt(mu2) : -sqrt(-mu2));
    }

    subjet = piece1;
  }
  

  // we've tagged the splitting, return the jet with its substructure
  StructureType * structure = new StructureType(subjet);
  structure->_symmetry = sym;
  structure->_mu       = (mu2 >= 0) ? sqrt(mu2) : -sqrt(-mu2);
  structure->_delta_R  = sqrt(squared_geometric_distance(piece1, piece2));
  if (_verbose_structure) {
    structure->_has_verbose = true;
    structure->_dropped_symmetry = dropped_symmetry;
    structure->_dropped_mu = dropped_mu;
    structure->_dropped_delta_R = dropped_delta_R;
  }
  subjet.set_structure_shared_ptr(SharedPtr<PseudoJetStructureBase>(structure));
  return subjet;
}


  
//----------------------------------------------------------------------
// the method below is the one actually performing one step of the
// recursion.
//
// It returns a status code (defined above)
//
// In case of success, all the information is filled
// In case of "no parents", piee1 is the same subjet
// In case of trouble, piece2 will be a 0 PJ and piece1 is the PJ we
//   should return (either 0 itself if the issue was critical, or
//   non-wero in case of a minor issue just causing the recursion to
//   stop)
RecursiveSymmetryCutBase::RecursionStatus
      RecursiveSymmetryCutBase::recurse_one_step(const PseudoJet & subjet,
                                                 PseudoJet &piece1, PseudoJet &piece2,
                                                 double &sym, double &mu2,
                                                 void *extra_parameters) const {
  if (!subjet.has_parents(piece1, piece2)){
    piece1 = subjet;    
    piece2 = PseudoJet();
    return recursion_no_parents;
  }

  // first sanity check: 
  // - zero or negative pts are not allowed for the input subjet
  // - zero or negative masses are not allowed for configurations
  //   in which the mass will effectively appear in a denominator
  //   (The masses will be checked later)
  if (subjet.pt2() <= 0){ // this is a critical problem, return an empty PJ
    piece1 = piece2 = PseudoJet();
    return recursion_issue;
  }

  if (_subtractor) {
    piece1 = (*_subtractor)(piece1);
    piece2 = (*_subtractor)(piece2);
  }
    
  // determine the symmetry parameter
  if        (_symmetry_measure == y) {
    // the original d_{ij}/m^2 choice from MDT
    // first make sure the mass is sensible
    if (subjet.m2() <= 0) { 
      _negative_mass_warning.warn("RecursiveSymmetryCutBase: cannot calculate y, because (sub)jet mass is negative; bailing out");
      // since rounding errors can give -ve masses, be a it more
      // tolerant and consider that no substructure has been found
      piece1 = _result_no_substructure(subjet);
      piece2 = PseudoJet();
      return recursion_issue;
    }
    sym = piece1.kt_distance(piece2) / subjet.m2();
    
  } else if (_symmetry_measure == vector_z) {
    // min(pt1, pt2)/(pt), where the denominator is a vector sum
    // of the two subjets
    sym = min(piece1.pt(), piece2.pt()) / subjet.pt();    
  } else if (_symmetry_measure == scalar_z) {
    // min(pt1, pt2)/(pt1+pt2), where the denominator is a scalar sum
    // of the two subjets
    double pt1 = piece1.pt();
    double pt2 = piece2.pt();
    // make sure denominator is non-zero
    sym = pt1 + pt2;
    if (sym == 0){ // this is a critical problem, return an empty PJ
      piece1 = piece2 = PseudoJet();
      return recursion_issue;
    }
    sym = min(pt1, pt2) / sym;
  } else if ((_symmetry_measure == theta_E) || (_symmetry_measure == cos_theta_E)){
    // min(E1, E2)/(E1+E2)
    double E1 = piece1.E();
    double E2 = piece2.E();
    // make sure denominator is non-zero
    sym = E1 + E2;
    if (sym == 0){ // this is a critical problem, return an empty PJ
      piece1 = piece2 = PseudoJet();
      return recursion_issue;
    }
    sym = min(E1, E2) / sym;
  } else {
    throw Error ("Unrecognized choice of symmetry_measure");
  }
  
  // determine the symmetry cut
  // (This function is specified in the derived classes)
  double this_symmetry_cut = symmetry_cut_fn(piece1, piece2, extra_parameters);
  
  // and make a first tagging decision based on symmetry cut
  bool tagged = (sym > this_symmetry_cut);

  // if tagged based on symmetry cut, then check the mu cut (if relevant)
  // and update the tagging decision. Calculate mu^2 regardless, for cases
  // of users not cutting on mu2, but still interested in its value.
  bool use_mu_cut = (_mu_cut != numeric_limits<double>::infinity());
  mu2 = max(piece1.m2(), piece2.m2())/subjet.m2();
  if (tagged && use_mu_cut) {
    // first a sanity check -- mu2 won't be sensible if the subjet mass 
    // is negative, so we can't then trust the mu cut - bail out
    if (subjet.m2() <= 0) {
      _negative_mass_warning.warn("RecursiveSymmetryCutBase: cannot trust mu, because (sub)jet mass is negative; bailing out");
      piece1 = piece2 = PseudoJet();
      return recursion_issue;
    }
    if (mu2 > 1) _mu2_gt1_warning.warn("RecursiveSymmetryCutBase encountered mu^2 value > 1");
    if (mu2 > pow(_mu_cut,2)) tagged = false;
  }

  // we'll continue unclustering, allowing for the different
  // ways of choosing which parent to look into
  if        (_recursion_choice == larger_pt) {
    if (piece1.pt2() < piece2.pt2()) std::swap(piece1, piece2);    
  } else if (_recursion_choice == larger_mt) {
    if (piece1.mt2() < piece2.mt2()) std::swap(piece1, piece2);    
  } else if (_recursion_choice == larger_m)  {
    if (piece1.m2() < piece2.m2()) std::swap(piece1, piece2);    
  } else if (_recursion_choice == larger_E)  {
    if (piece1.E()  < piece2.E())  std::swap(piece1, piece2);    
  } else {
    throw Error ("Unrecognized value for recursion_choice");
  }    

  return tagged ? recursion_success : recursion_dropped;
}

  
//----------------------------------------------------------------------
string RecursiveSymmetryCutBase::description() const {
  ostringstream ostr;
  ostr << "Recursive " << (_grooming_mode ? "Groomer" : "Tagger") << " with a symmetry cut ";

  switch(_symmetry_measure) {
  case y:
    ostr << "y"; break;
  case scalar_z:
    ostr << "scalar_z"; break;
  case vector_z:
    ostr << "vector_z"; break;
  case theta_E:
    ostr << "theta_E"; break;
  case cos_theta_E:
    ostr << "cos_theta_E"; break;
  default:
    cerr << "failed to interpret symmetry_measure" << endl; exit(-1);
  }
  ostr << " > " << symmetry_cut_description();

  if (_mu_cut != numeric_limits<double>::infinity()) {
    ostr << ", mass-drop cut mu=max(m1,m2)/m < " << _mu_cut;
  } else {
    ostr << ", no mass-drop requirement";
  }

  ostr << ", recursion into the subjet with larger ";
  switch(_recursion_choice) {
  case larger_pt:
    ostr << "pt"; break;
  case larger_mt:
    ostr << "mt(=sqrt(m^2+pt^2))"; break;
  case larger_m:
    ostr << "mass"; break;
  case larger_E:
    ostr << "energy"; break;
  default:
    cerr << "failed to interpret recursion_choice" << endl; exit(-1);
  }

  if (_subtractor) {
    ostr << ", subtractor: " << _subtractor->description();
    if (_input_jet_is_subtracted) {ostr << " (input jet is assumed already subtracted)";}
  }

  if (_recluster) {
    ostr << " and reclustering using " << _recluster->description();
  }
  
  return ostr.str();
}

//----------------------------------------------------------------------
// helper for handling the reclustering
PseudoJet RecursiveSymmetryCutBase::_recluster_if_needed(const PseudoJet &jet) const{
  if (! _do_reclustering) return jet;
  if (_recluster) return (*_recluster)(jet);
  if (is_ee()){
#if FASTJET_VERSION_NUMBER >= 30100
    return Recluster(JetDefinition(ee_genkt_algorithm, JetDefinition::max_allowable_R, 0.0), true)(jet);
#else
    return Recluster(JetDefinition(ee_genkt_algorithm, JetDefinition::max_allowable_R, 0.0))(jet);
#endif
  }

  return Recluster(cambridge_algorithm, JetDefinition::max_allowable_R)(jet);
}  
  
//----------------------------------------------------------------------
// decide what to return when no substructure has been found
double RecursiveSymmetryCutBase::squared_geometric_distance(const PseudoJet &j1,
                                                            const PseudoJet &j2) const{
  if (_symmetry_measure == theta_E){
    double dot_3d = j1.px()*j2.px() + j1.py()*j2.py() + j1.pz()*j2.pz();
    double cos_theta = max(-1.0,min(1.0, dot_3d/sqrt(j1.modp2()*j2.modp2())));
    double theta = acos(cos_theta);
    return theta*theta;
  } else if (_symmetry_measure == cos_theta_E){
    double dot_3d = j1.px()*j2.px() + j1.py()*j2.py() + j1.pz()*j2.pz();
    return max(0.0, 2*(1-dot_3d/sqrt(j1.modp2()*j2.modp2())));
  }

  return j1.squared_distance(j2);
}

//----------------------------------------------------------------------
PseudoJet RecursiveSymmetryCutBase::_result_no_substructure(const PseudoJet &last_parent) const{
  if (_grooming_mode){
    // in grooming mode, return the last parent
    return last_parent;
  } else {
    // in tagging mode, return an empty PseudoJet
    return PseudoJet();
  }
}


//========================================================================
// implementation of the details of the structure
  
// the number of dropped subjets
int RecursiveSymmetryCutBase::StructureType::dropped_count(bool global) const {
  check_verbose("dropped_count()");

  // if this jet has no substructure, just return an empty vector
  if (!has_substructure()) return _dropped_delta_R.size();

  // deal with the non-global case
  if (!global) return _dropped_delta_R.size();

  // for the global case, we've unfolded the recursion (likely more
  // efficient as it requires less copying)
  unsigned int count = 0;
  vector<const RecursiveSymmetryCutBase::StructureType*> to_parse;
  to_parse.push_back(this);

  unsigned int i_parse = 0;
  while (i_parse<to_parse.size()){
    const RecursiveSymmetryCutBase::StructureType *current = to_parse[i_parse];
    count += current->_dropped_delta_R.size();

    // check if we need to recurse deeper in the substructure
    //
    // we can have 2 situations here for the underlying structure (the
    // one we've wrapped around):
    //  - it's of the clustering type
    //  - it's a composite jet
    // only in the 2nd case do we have to recurse deeper
    const CompositeJetStructure *css = dynamic_cast<const CompositeJetStructure*>(current->_structure.get());
    if (css == 0){ ++i_parse; continue; }
    
    vector<PseudoJet> prongs = css->pieces(PseudoJet()); // argument irrelevant
    assert(prongs.size() == 2);
    for (unsigned int i_prong=0; i_prong<2; ++i_prong){
      if (prongs[i_prong].has_structure_of<RecursiveSymmetryCutBase>()){
        RecursiveSymmetryCutBase::StructureType* prong_structure
          = (RecursiveSymmetryCutBase::StructureType*) prongs[i_prong].structure_ptr();
        if (prong_structure->has_substructure())
          to_parse.push_back(prong_structure);
      }
    }

    ++i_parse;
  }
  return count;
}

// the delta_R of all the dropped subjets
vector<double> RecursiveSymmetryCutBase::StructureType::dropped_delta_R(bool global) const {
  check_verbose("dropped_delta_R()");

  // if this jet has no substructure, just return an empty vector
  if (!has_substructure()) return vector<double>();

  // deal with the non-global case
  if (!global) return _dropped_delta_R;

  // for the global case, we've unfolded the recursion (likely more
  // efficient as it requires less copying)
  vector<double> all_dropped;
  vector<const RecursiveSymmetryCutBase::StructureType*> to_parse;
  to_parse.push_back(this);

  unsigned int i_parse = 0;
  while (i_parse<to_parse.size()){
    const RecursiveSymmetryCutBase::StructureType *current = to_parse[i_parse];
    all_dropped.insert(all_dropped.end(), current->_dropped_delta_R.begin(), current->_dropped_delta_R.end());

    // check if we need to recurse deeper in the substructure
    //
    // we can have 2 situations here for the underlying structure (the
    // one we've wrapped around):
    //  - it's of the clustering type
    //  - it's a composite jet
    // only in the 2nd case do we have to recurse deeper
    const CompositeJetStructure *css = dynamic_cast<const CompositeJetStructure*>(current->_structure.get());
    if (css == 0){ ++i_parse; continue; }
    
    vector<PseudoJet> prongs = css->pieces(PseudoJet()); // argument irrelevant
    assert(prongs.size() == 2);
    for (unsigned int i_prong=0; i_prong<2; ++i_prong){
      if (prongs[i_prong].has_structure_of<RecursiveSymmetryCutBase>()){
        RecursiveSymmetryCutBase::StructureType* prong_structure
          = (RecursiveSymmetryCutBase::StructureType*) prongs[i_prong].structure_ptr();
        if (prong_structure->has_substructure())
          to_parse.push_back(prong_structure);
      }
    }

    ++i_parse;
  }
  return all_dropped;
}

// the symmetry of all the dropped subjets
vector<double> RecursiveSymmetryCutBase::StructureType::dropped_symmetry(bool global) const {
  check_verbose("dropped_symmetry()");

  // if this jet has no substructure, just return an empty vector
  if (!has_substructure()) return vector<double>();

  // deal with the non-global case
  if (!global) return _dropped_symmetry;

  // for the global case, we've unfolded the recursion (likely more
  // efficient as it requires less copying)
  vector<double> all_dropped;
  vector<const RecursiveSymmetryCutBase::StructureType*> to_parse;
  to_parse.push_back(this);

  unsigned int i_parse = 0;
  while (i_parse<to_parse.size()){
    const RecursiveSymmetryCutBase::StructureType *current = to_parse[i_parse];
    all_dropped.insert(all_dropped.end(), current->_dropped_symmetry.begin(), current->_dropped_symmetry.end());

    // check if we need to recurse deeper in the substructure
    //
    // we can have 2 situations here for the underlying structure (the
    // one we've wrapped around):
    //  - it's of the clustering type
    //  - it's a composite jet
    // only in the 2nd case do we have to recurse deeper
    const CompositeJetStructure *css = dynamic_cast<const CompositeJetStructure*>(current->_structure.get());
    if (css == 0){ ++i_parse; continue; }
    
    vector<PseudoJet> prongs = css->pieces(PseudoJet()); // argument irrelevant
    assert(prongs.size() == 2);
    for (unsigned int i_prong=0; i_prong<2; ++i_prong){
      if (prongs[i_prong].has_structure_of<RecursiveSymmetryCutBase>()){
        RecursiveSymmetryCutBase::StructureType* prong_structure
          = (RecursiveSymmetryCutBase::StructureType*) prongs[i_prong].structure_ptr();
        if (prong_structure->has_substructure())
          to_parse.push_back(prong_structure);
      }
    }

    ++i_parse;
  }
  return all_dropped;
}

// the mu of all the dropped subjets
vector<double> RecursiveSymmetryCutBase::StructureType::dropped_mu(bool global) const {
  check_verbose("dropped_mu()");

  // if this jet has no substructure, just return an empty vector
  if (!has_substructure()) return vector<double>();

  // deal with the non-global case
  if (!global) return _dropped_mu;

  // for the global case, we've unfolded the recursion (likely more
  // efficient as it requires less copying)
  vector<double> all_dropped;
  vector<const RecursiveSymmetryCutBase::StructureType*> to_parse;
  to_parse.push_back(this);

  unsigned int i_parse = 0;
  while (i_parse<to_parse.size()){
    const RecursiveSymmetryCutBase::StructureType *current = to_parse[i_parse];
    all_dropped.insert(all_dropped.end(), current->_dropped_mu.begin(), current->_dropped_mu.end());

    // check if we need to recurse deeper in the substructure
    //
    // we can have 2 situations here for the underlying structure (the
    // one we've wrapped around):
    //  - it's of the clustering type
    //  - it's a composite jet
    // only in the 2nd case do we have to recurse deeper
    const CompositeJetStructure *css = dynamic_cast<const CompositeJetStructure*>(current->_structure.get());
    if (css == 0){ ++i_parse; continue; }
    
    vector<PseudoJet> prongs = css->pieces(PseudoJet()); // argument irrelevant
    assert(prongs.size() == 2);
    for (unsigned int i_prong=0; i_prong<2; ++i_prong){
      if (prongs[i_prong].has_structure_of<RecursiveSymmetryCutBase>()){
        RecursiveSymmetryCutBase::StructureType* prong_structure
          = (RecursiveSymmetryCutBase::StructureType*) prongs[i_prong].structure_ptr();
        if (prong_structure->has_substructure())
          to_parse.push_back(prong_structure);
      }
    }

    ++i_parse;
  }
  return all_dropped;
}

// the maximum of the symmetry over the dropped subjets
double RecursiveSymmetryCutBase::StructureType::max_dropped_symmetry(bool global) const {
  check_verbose("max_dropped_symmetry()");

  // if there is no substructure, just exit
  if (!has_substructure()){ return 0.0; }

  // local value of the max_dropped_symmetry
  double local_max = (_dropped_symmetry.size() == 0)
    ? 0.0 : *max_element(_dropped_symmetry.begin(),_dropped_symmetry.end());

  // recurse down the structure if instructed to do so
  if (global){
    // we can have 2 situations here for the underlying structure (the
    // one we've wrapped around):
    //  - it's of the clustering type
    //  - it's a composite jet
    // only in the 2nd case do we have to recurse deeper
    const CompositeJetStructure *css = dynamic_cast<const CompositeJetStructure*>(_structure.get());
    if (css == 0) return local_max;
    
    vector<PseudoJet> prongs = css->pieces(PseudoJet()); // argument irrelevant
    assert(prongs.size() == 2);
    for (unsigned int i_prong=0; i_prong<2; ++i_prong){
      // check if the prong has further substructure
      if (prongs[i_prong].has_structure_of<RecursiveSymmetryCutBase>()){
        RecursiveSymmetryCutBase::StructureType* prong_structure
          = (RecursiveSymmetryCutBase::StructureType*) prongs[i_prong].structure_ptr();
        local_max = max(local_max, prong_structure->max_dropped_symmetry(true));
      }
    }
  }

  return local_max;
}

//------------------------------------------------------------------------
// helper class to sort by decreasing thetag
class SortRecursiveSoftDropStructureZgThetagPair{
public:
  bool operator()(const pair<double, double> &p1, const pair<double, double> &p2) const{
    return p1.second > p2.second;
  }
};
//------------------------------------------------------------------------

// the (zg,thetag) pairs of all the splitting that were found and passed the SD condition
vector<pair<double,double> > RecursiveSymmetryCutBase::StructureType::sorted_zg_and_thetag() const {
  //check_verbose("sorted_zg_and_thetag()");

  // if this jet has no substructure, just return an empty vector
  if (!has_substructure()) return vector<pair<double,double> >();
  
  // otherwise fill a vector with all the prongs (no specific ordering)
  vector<pair<double,double> > all;
  vector<const RecursiveSymmetryCutBase::StructureType*> to_parse;
  to_parse.push_back(this);

  unsigned int i_parse = 0;
  while (i_parse<to_parse.size()){
    const RecursiveSymmetryCutBase::StructureType *current = to_parse[i_parse];
    all.push_back(pair<double,double>(current->_symmetry, current->_delta_R));

    vector<PseudoJet> prongs = current->pieces(PseudoJet());
    assert(prongs.size() == 2);
    for (unsigned int i_prong=0; i_prong<2; ++i_prong){
      if (prongs[i_prong].has_structure_of<RecursiveSymmetryCutBase>()){
        RecursiveSymmetryCutBase::StructureType* prong_structure
          = (RecursiveSymmetryCutBase::StructureType*) prongs[i_prong].structure_ptr();
        if (prong_structure->has_substructure())
          to_parse.push_back(prong_structure);
      }
    }

    ++i_parse;
  }

  sort(all.begin(), all.end(), SortRecursiveSoftDropStructureZgThetagPair());
  return all;
}

} // namespace contrib

FASTJET_END_NAMESPACE
