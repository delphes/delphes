// $Id: RecursiveSoftDrop.hh 1082 2017-10-10 12:00:13Z gsoyez $
//
// Copyright (c) 2014-, Gavin P. Salam, Gregory Soyez, Jesse Thaler,
// Kevin Zhou, Frederic Dreyer
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

#ifndef __RECURSIVESOFTDROP_HH__
#define __RECURSIVESOFTDROP_HH__

#include "Recluster.hh"
#include "SoftDrop.hh"
#include "fastjet/WrappedStructure.hh"

#include <iostream>
#include <queue>
#include <vector>

FASTJET_BEGIN_NAMESPACE

namespace contrib{

//------------------------------------------------------------------------
/// \class RecursiveSoftDrop
/// An implementation of the RecursiveSoftDrop.
///
/// Recursive Soft Drop will recursively groom a jet, removing
/// particles that fail the criterion
/// \f[
///     z > z_{\rm cut} (\theta/R0)^\beta
/// \f]
/// until n subjets have been found.
///
/// Several variants are supported:
///  - set_fixed_depth_mode() switches to fixed depth on all branches
///    of the clustering tree
///  - set_dynamical_R0() switches to dynamical R0 implementation of
///    RSD
///  - set_hardest_branch_only() switches to following only the
///    hardest branch (e.g. for Iterated Soft Drop)
///  - set_min_deltaR_square(val) sets a minimum angle considered for
///    substructure (e.g. for Iterated Soft Drop)
///
/// Notes:
///
///  - Even though the calls to "set_tagging_mode()" or
///    "set_grooming_mode(false)" are allowed, they should not be used
///    with n=-1, and the default grooming_mode has to remain
///    untouched (except for beta<0 and finite n).
///
//----------------------------------------------------------------------
class RecursiveSoftDrop : public SoftDrop {
public:
  /// Simplified constructor. This takes the value of the "beta"
  /// parameter and the symmetry cut (applied by default on the
  /// scalar_z variable, as for the mMDT). It also takes an optional
  /// subtractor.
  ///
  /// n is the number of times we require the SoftDrop condition to be
  /// satisfied. n=-1 means infinity, i.e. we recurse into the jet
  /// until individual constituents
  ///
  /// If the (optional) pileup subtractor can be supplied, then see
  /// also the documentation for the set_input_jet_is_subtracted() member
  /// function.
  ///
  /// \param beta               the value of the beta parameter
  /// \param symmetry_cut       the value of the cut on the symmetry measure
  /// \param n                  the requested number of iterations
  /// \param R0                 the angular distance normalisation [1 by default]
  RecursiveSoftDrop(double beta,
                    double symmetry_cut,
                    int n = -1,
                    double R0 = 1,
                    const FunctionOfPseudoJet<PseudoJet> * subtractor = 0) : 
    SoftDrop(beta, symmetry_cut, R0, subtractor), _n(n) { set_defaults(); }

  /// Full constructor, which takes the following parameters:
  ///
  /// \param beta               the value of the beta parameter
  /// \param symmetry_cut       the value of the cut on the symmetry measure
  /// \param symmetry_measure   the choice of measure to use to estimate the symmetry
  /// \param n                  the requested number of iterations
  /// \param R0                 the angular distance normalisation [1 by default]
  /// \param mu_cut             the maximal allowed value of mass drop variable mu = m_heavy/m_parent 
  /// \param recursion_choice   the strategy used to decide which subjet to recurse into
  /// \param subtractor         an optional pointer to a pileup subtractor (ignored if zero)
  RecursiveSoftDrop(double           beta,
                    double           symmetry_cut, 
                    SymmetryMeasure  symmetry_measure,
                    int              n = -1,
                    double           R0 = 1.0,
                    double           mu_cut = std::numeric_limits<double>::infinity(), 
                    RecursionChoice  recursion_choice = larger_pt,
                    const FunctionOfPseudoJet<PseudoJet> * subtractor = 0) : 
    SoftDrop(beta, symmetry_cut, symmetry_measure, R0, mu_cut, recursion_choice, subtractor),
    _n(n) { set_defaults(); }

  /// default destructor
  virtual ~RecursiveSoftDrop(){}

  //----------------------------------------------------------------------
  // access to class info
  int n() const { return _n; }
  
  //----------------------------------------------------------------------
  // on top of the tweaks that we inherit from SoftDrop (via
  // RecursiveSymmetryBase):
  //  - set_verbose_structure()
  //  - set_subtractor()
  //  - set_input_jet_is_subtracted()
  // we provide several other knobs, given below

  /// initialise all the flags below to their default value
  void set_defaults();
  
  /// switch to using the "same depth" variant where instead of
  /// recursing from large to small angles and requiring n SD
  /// conditions to be met (our default), we recurse simultaneously in
  /// all the branches found during the previous iteration, up to a
  /// maximum depth of n.
  /// default: false
  void set_fixed_depth_mode(bool value=true) { _fixed_depth = value; }
  bool fixed_depth_mode() const { return _fixed_depth; }
  
  /// switch to using a dynamical R0 (used for the normalisation of
  /// the symmetry measure) set by the last deltaR at which some
  /// substructure was found.
  /// default: false
  void set_dynamical_R0(bool value=true) { _dynamical_R0 = value; }
  bool use_dynamical_R0() const { return _dynamical_R0; }

  /// when finding some substructure, only follow the hardest branch
  /// for the recursion
  /// default: false (i.e. recurse in both branches)
  void set_hardest_branch_only(bool value=true) { _hardest_branch_only = value; }
  bool use_hardest_branch_only() const { return _hardest_branch_only; }

  /// set the minimum angle (squared) that we should consider for
  /// substructure
  /// default: -1.0 (i.e. no minimum)
  void set_min_deltaR_squared(double value=-1.0) { _min_dR2 = value; }
  double   min_deltaR_squared() const { return _min_dR2; }

  /// description of the tool
  virtual std::string description() const;

  //----------------------------------------------------------------------
  /// action on a single jet with RecursiveSoftDrop.
  ///
  /// uses "result_fixed_tags" by default (i.e. recurse from R0 to
  /// smaller angles until n SD conditions have been met), or
  /// "result_fixed_depth" where each of the previous SD branches are
  /// recirsed into down to a depth of n.
  virtual PseudoJet result(const PseudoJet &jet) const;

  /// this routine applies the Soft Drop criterion recursively on the
  /// CA tree until we find n subjets (or until it converges), and
  /// adds them together into a groomed PseudoJet
  PseudoJet result_fixed_tags(const PseudoJet &jet) const;

  /// this routine applies the Soft Drop criterion recursively on the
  /// CA tree, recursing into all the branches found during the previous iteration
  /// until n layers have been found (or until it converges)
  PseudoJet result_fixed_depth(const PseudoJet &jet) const;
    
protected:  
  /// return false if we reached desired layer of grooming _n
  bool continue_grooming(int current_n) const {
    return ((_n < 0) or (current_n < _n));
  }
  
private:
  int    _n;            ///< the value of n

  // behaviour tweaks
  bool _fixed_depth;         ///< look in parallel into each all branches until depth n
  bool _dynamical_R0;        ///< when true, use the last deltaR with substructure as D0
  bool _hardest_branch_only; ///< recurse only in the hardest branch
                             ///  when substructure is found
  double _min_dR2;           ///< the min allowed angle to search for substructure
};

// helper to get the (linear) list of prongs inside a jet resulting
// from RecursiveSoftDrop. This would avoid having amnually to go
// through the successive pairwise compositeness
std::vector<PseudoJet> recursive_soft_drop_prongs(const PseudoJet & rsd_jet);

}

FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh
#endif // __RECURSIVESOFTDROP_HH__
