// $Id: RecursiveSymmetryCutBase.cc 700 2014-07-07 12:50:05Z gsoyez $
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

  PseudoJet j = 
    _do_reclustering 
      ? _recluster ? (*_recluster)(jet)
                   : Recluster(cambridge_algorithm, JetDefinition::max_allowable_R)(jet)
      : jet;
    
  // issue a warning if the jet is not obtained through a C/A
  // clustering
  // if ((! j.has_associated_cluster_sequence()) ||
  //     (j.validated_cs()->jet_def().jet_algorithm() != cambridge_algorithm))
  //  _nonca_warning.warn("RecursiveSymmetryCutBase is designed to be applied on jets from a Cambridge/Aachen clustering; use it with other algorithms at your own risk.");

  if (! j.has_valid_cluster_sequence()){
    throw Error("RecursiveSymmetryCutBase can only be applied to jets associated to a (valid) cluster sequence");
  }

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

  bool use_mu_cut = (_mu_cut != numeric_limits<double>::infinity());

  // variables for tracking what will happen
  PseudoJet piece1, piece2;
 
  // vectors for storing optional verbose structure
  // these hold the deltaR, symmetry, and mu values of dropped branches
  std::vector<double> dropped_delta_R;
  std::vector<double> dropped_symmetry;
  std::vector<double> dropped_mu;
  
  // now recurse into the jet's structure
  while (subjet.has_parents(piece1, piece2)) {
    
    // first sanity check: 
    // - zero or negative pts are not allowed for the input subjet
    // - zero or negative masses are not allowed for configurations
    //   in which the mass will effectively appear in a denominator
    //   (The masses will be checked later)
    if (subjet.pt2() <= 0) return PseudoJet();

    if (_subtractor) {
      piece1 = (*_subtractor)(piece1);
      piece2 = (*_subtractor)(piece2);
    }
    
    // determine the symmetry parameter
    double sym;

    if        (_symmetry_measure == y) {
      // the original d_{ij}/m^2 choice from MDT
      // first make sure the mass is sensible
      if (subjet.m2() <= 0) {
        _negative_mass_warning.warn("RecursiveSymmetryCutBase: cannot calculate y, because (sub)jet mass is negative; bailing out");
        return _result_no_substructure(subjet); //TBC: do we return the hardest parent? A NULL PseudoJet?
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
      if (sym == 0) return PseudoJet(); 
      sym = min(pt1, pt2) / sym;

    } else {
      throw Error ("Unrecognized choice of symmetry_measure");
    }

    // determine the symmetry cut
    // (This function is specified in the derived classes)
    double this_symmetry_cut = symmetry_cut_fn(piece1, piece2);

    // and make a first tagging decision based on symmetry cut
    bool tagged = (sym > this_symmetry_cut);

    // if tagged based on symmetry cut, then check the mu cut (if relevant)
    // and update the tagging decision. Calculate mu^2 regardless, for cases
    // of users not cutting on mu2, but still interested in its value.
    double mu2 = max(piece1.m2(), piece2.m2())/subjet.m2();
    if (tagged && use_mu_cut) {
      // first a sanity check -- mu2 won't be sensible if the subjet mass 
      // is negative, so we can't then trust the mu cut - bail out
      if (subjet.m2() <= 0) {
        _negative_mass_warning.warn("RecursiveSymmetryCutBase: cannot trust mu, because (sub)jet mass is negative; bailing out");
        return PseudoJet();
      }
      if (mu2 > 1) _mu2_gt1_warning.warn("RecursiveSymmetryCutBase encountered mu^2 value > 1");
      if (mu2 > pow(_mu_cut,2)) tagged = false;
    }


    // if we've tagged the splitting, return the jet with its substructure
    if (tagged) {
      // record relevant information
      StructureType * structure = new StructureType(subjet);
      structure->_symmetry = sym;
      structure->_mu       = (mu2 >= 0) ? sqrt(mu2) : -sqrt(-mu2);
      structure->_delta_R  = piece1.delta_R(piece2);
      if (_verbose_structure) {
        structure->_has_verbose = true;
        structure->_dropped_symmetry = dropped_symmetry;
        structure->_dropped_mu = dropped_mu;
        structure->_dropped_delta_R = dropped_delta_R;
      }
      subjet.set_structure_shared_ptr(SharedPtr<PseudoJetStructureBase>(structure));
      return subjet;
    }
    
    // if desired, store information about dropped branches before recursing
    if (_verbose_structure) {
      dropped_delta_R.push_back(piece1.delta_R(piece2));
      dropped_symmetry.push_back(sym);
      dropped_mu.push_back((mu2 >= 0) ? sqrt(mu2) : -sqrt(-mu2));
    }
    
    // otherwise continue unclustering, allowing for the different
    // ways of choosing which parent to look into
    int choice;
    if        (_recursion_choice == larger_mt) {
      choice = piece1.mt2() > piece2.mt2() ? 1 : 2;

    } else if (_recursion_choice == larger_pt) {
      choice = piece1.pt2() > piece2.pt2() ? 1 : 2;

    } else if (_recursion_choice == larger_m)  {
      choice = piece1.m2()  > piece2.m2()  ? 1 : 2;

    } else {
      throw Error ("Unrecognized value for recursion_choice");
    }    
    if (_verbose) cout << "choice is " << choice << endl;;
    subjet = (choice == 1) ? piece1 : piece2;
  } // (subjet.has_parents(...))

  if (_verbose) cout << "reached end; returning null jet " << endl;
  
  // decide on tagging versus grooming mode here
  PseudoJet result = _result_no_substructure(subjet);
  
  if (result != 0) {
    // if in grooming mode, add dummy structure information
    StructureType * structure = new StructureType(result);
    structure->_symmetry = 0.0;
    structure->_mu       = 0.0;
    structure->_delta_R  = 0.0;
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
  default:
    cerr << "failed to interpret recursion_choice" << endl; exit(-1);
  }

  if (_subtractor) {
    ostr << " and subtractor: " << _subtractor->description();
    if (_input_jet_is_subtracted) {ostr << " (input jet is assumed already subtracted)";}
  }
  return ostr.str();
}

// decide what to return when no substructure has been found
PseudoJet RecursiveSymmetryCutBase::_result_no_substructure(const PseudoJet &last_parent) const{
  if (_grooming_mode){
    // in grooming mode, return the last parent
    return last_parent;
  } else {
    // in tagging mode, return an empty PseudoJet
    return PseudoJet();
  }
}


} // namespace contrib

FASTJET_END_NAMESPACE
