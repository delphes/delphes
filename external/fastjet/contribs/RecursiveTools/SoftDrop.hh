// $Id: SoftDrop.hh 686 2014-06-14 03:25:09Z jthaler $
//
// Copyright (c) 2014-, Gregory Soyez, Jesse Thaler
// based on arXiv:1402.2657 by Andrew J. Larkoski, Simone Marzani,
// Gregory Soyez, Jesse Thaler
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

#ifndef __FASTJET_CONTRIB_SOFTDROP_HH__
#define __FASTJET_CONTRIB_SOFTDROP_HH__

#include "RecursiveSymmetryCutBase.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


namespace contrib{

//------------------------------------------------------------------------
/// \class SoftDrop
/// An implementation of the SoftDrop from arXiv:1402.2657.
///
/// For the basic functionalities, we refer the reader to the
/// documentation of the RecursiveSymmetryCutBase from which SoftDrop
/// inherits. Here, we mostly put the emphasis on things specific to
/// SoftDrop:
///
///  - the cut applied recursively is 
///     \f[
///        z > z_{\rm cut} (\theta/R0)^\beta
///     \f]
///    with z the asymmetry measure and \f$\theta\f$ the geometrical
///    distance between the two subjets. R0 is set to 1 by default.
///
///  - by default, we work in "grooming mode" i.s. if no substructure
///    is found, we return a jet made of a single parton. Note that
///    this behaviour differs from the mMDT (and can be a source of
///    differences when running SoftDrop with beta=0.)
///
class SoftDrop : public RecursiveSymmetryCutBase {
public:
  /// Simplified constructor. This takes the value of the "beta"
  /// parameter and the symmetry cut (applied by default on the
  /// scalar_z variable, as for the mMDT). It also takes an optional
  /// subtractor.
  ///
  /// If the (optional) pileup subtractor can be supplied, then see
  /// also the documentation for the set_input_jet_is_subtracted() member
  /// function.
  ///
  /// \param beta               the value of the beta parameter
  /// \param symmetry_cut       the value of the cut on the symmetry measure
  /// \param R0                 the angular distance normalisation [1 by default]
  SoftDrop(double beta,
           double symmetry_cut,
           double R0 = 1,
           const FunctionOfPseudoJet<PseudoJet> * subtractor = 0) : 
    RecursiveSymmetryCutBase(scalar_z,  // the default SymmetryMeasure
                             std::numeric_limits<double>::infinity(), // default is no mass drop
                             larger_pt, // the default RecursionChoice
                             subtractor),
    _beta(beta), _symmetry_cut(symmetry_cut), _R0sqr(R0*R0) {
    // change the default: use grooming mode
    set_grooming_mode();
  }

  /// Full constructor, which takes the following parameters:
  ///
  /// \param beta               the value of the beta parameter
  /// \param symmetry_cut       the value of the cut on the symmetry measure
  /// \param symmetry_measure   the choice of measure to use to estimate the symmetry
  /// \param R0                 the angular distance normalisation [1 by default]
  /// \param mu_cut             the maximal allowed value of mass drop variable mu = m_heavy/m_parent 
  /// \param recursion_choice   the strategy used to decide which subjet to recurse into
  /// \param subtractor         an optional pointer to a pileup subtractor (ignored if zero)
  ///
  /// The default values provided for this constructor are suited to
  /// obtain the SoftDrop as discussed in arXiv:1402.2657:
  ///  - no mass drop is requested
  ///  - recursion follows the branch with the largest pt
  /// The symmetry measure has to be specified (scalar_z is the recommended value)
  ///
  /// Notes: 
  ///
  /// - by default, SoftDrop will recluster the jet with the
  ///   Cambridge/Aachen algorithm if it is not already the case. This
  ///   behaviour can be changed using the "set_reclustering" method
  ///   defined below
  ///
  SoftDrop(double           beta,
           double           symmetry_cut, 
           SymmetryMeasure  symmetry_measure,
           double           R0 = 1.0,
           double           mu_cut = std::numeric_limits<double>::infinity(), 
           RecursionChoice  recursion_choice = larger_pt,
           const FunctionOfPseudoJet<PseudoJet> * subtractor = 0) : 
     RecursiveSymmetryCutBase(symmetry_measure, mu_cut, recursion_choice, subtractor),
     _beta(beta), _symmetry_cut(symmetry_cut), _R0sqr(R0*R0)
  {}

  /// default destructor
  virtual ~SoftDrop(){}

protected:

  // Unlike MMDT, the SoftDrop symmetry_cut_fn depends on the subjet kinematics
  // since the symmetry condition depends on the DeltaR between subjets.
  virtual double symmetry_cut_fn(const PseudoJet & /* p1 */, 
                                 const PseudoJet & /* p2 */) const;
  virtual std::string symmetry_cut_description() const;

private:
  double _beta;         ///< the power of the angular distance to be used
                        ///< in the symmetry condition
  double _symmetry_cut; ///< the value of zcut (the prefactor in the asymmetry cut)
  double _R0sqr;        ///< normalisation of the angular distance
                        ///< (typically set to the jet radius, 1 by default)
};

} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_SOFTDROP_HH__
