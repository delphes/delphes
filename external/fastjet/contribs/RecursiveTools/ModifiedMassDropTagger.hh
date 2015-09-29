// $Id: ModifiedMassDropTagger.hh 688 2014-06-17 14:29:56Z jthaler $
//
// Copyright (c) 2014-, Gavin P. Salam
// based on arXiv:1307.007 by Mrinal Dasgupta, Simone Marzani and Gavin P. Salam
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

#ifndef __FASTJET_CONTRIB_MODIFIEDMASSDROPTAGGER_HH__
#define __FASTJET_CONTRIB_MODIFIEDMASSDROPTAGGER_HH__

#include "RecursiveSymmetryCutBase.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

//------------------------------------------------------------------------
/// \class ModifiedMassDropTagger
/// An implementation of the modified Mass-Drop Tagger from arXiv:1307.0007.
///
class ModifiedMassDropTagger : public RecursiveSymmetryCutBase {
public:
  
  /// Simplified constructor, which takes just a symmetry cut (applied
  /// on the scalar_z variable) and an optional subtractor.
  ///
  /// In this incarnation the ModifiedMassDropTagger is a bit of a
  /// misnomer, because there is no mass-drop condition
  /// applied. Recursion into the jet structure chooses the prong with
  /// largest pt. (Results from arXiv:1307.0007 were based on the
  /// largest mt, but this only makes a difference for values of the
  /// symmetry_cut close to 1/2).
  ///
  /// If the (optional) pileup subtractor can be supplied, then see
  /// also the documentation for the set_input_jet_is_subtracted() member
  /// function.
  ///
  /// NB: The configuration of MMDT provided by this constructor is
  /// probably the most robust for use with subtraction.
  ModifiedMassDropTagger(double symmetry_cut, 
                         const FunctionOfPseudoJet<PseudoJet> * subtractor = 0
                         ) :
    RecursiveSymmetryCutBase(scalar_z,  // the default SymmetryMeasure
                             std::numeric_limits<double>::infinity(), // the default is no mass drop
                             larger_pt, // the default RecursionChoice
                             subtractor),
    _symmetry_cut(symmetry_cut)
  {}

  /// Full constructor, which takes the following parameters:
  ///
  /// \param symmetry_cut       the value of the cut on the symmetry measure
  /// \param symmetry_measure   the choice of measure to use to estimate the symmetry
  /// \param mu_cut             the maximal allowed value of mass drop variable mu = m_heavy/m_parent 
  /// \param recursion_choice   the strategy used to decide which subjet to recurse into
  /// \param subtractor         an optional pointer to a pileup subtractor (ignored if zero)
  ///
  /// To obtain the mMDT as discussed in arXiv:1307.0007, use an
  /// symmetry_measure that's one of the following
  ///
  /// - RecursiveSymmetryCutBase::y         (for a cut on y)
  /// - RecursiveSymmetryCutBase::scalar_z  (for a cut on z)
  ///
  /// and use the default recursion choice of
  /// RecursiveSymmetryCutBase::larger_pt (larger_mt will give something
  /// very similar, while larger_m will give the behaviour of the
  /// original, but now deprecated MassDropTagger)
  ///
  /// Notes: 
  ///
  /// - By default the ModifiedMassDropTagger will relcuster the jets
  ///   with the C/A algorithm (if needed).
  ///
  /// - the mu_cut parameter is mostly irrelevant when it's taken
  ///   larger than about 1/2: the tagger is then one that cuts
  ///   essentially on the (a)symmetry of the jet's momentum
  ///   sharing. The default value of infinity turns off its use
  ///   entirely
  ModifiedMassDropTagger(double           symmetry_cut, 
                         SymmetryMeasure  symmetry_measure,
                         double           mu_cut = std::numeric_limits<double>::infinity(), 
                         RecursionChoice  recursion_choice = larger_pt,
                         const FunctionOfPseudoJet<PseudoJet> * subtractor = 0
                         ) : 
    RecursiveSymmetryCutBase(symmetry_measure, mu_cut, recursion_choice, subtractor),
    _symmetry_cut(symmetry_cut)
  {}

  /// default destructor
  virtual ~ModifiedMassDropTagger(){}

protected:

  /// The symmetry cut function for MMDT returns just a constant, since the cut value
  /// has no dependence on the subjet kinematics
  virtual double symmetry_cut_fn(const PseudoJet & /* p1 */, 
                                 const PseudoJet & /* p2 */
                                 ) const {return _symmetry_cut;}
  virtual std::string symmetry_cut_description() const;

  double           _symmetry_cut;
};

} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_MODIFIEDMASSDROPTAGGER_HH__
