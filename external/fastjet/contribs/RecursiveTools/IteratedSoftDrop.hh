// $Id: IteratedSoftDrop.hh 1086 2017-10-11 08:07:26Z gsoyez $
//
// Copyright (c) 2017-, Jesse Thaler, Kevin Zhou, Gavin P. Salam,
// Gregory Soyez
//
// based on arXiv:1704.06266 by Christopher Frye, Andrew J. Larkoski,
// Jesse Thaler, Kevin Zhou
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

#ifndef __FASTJET_CONTRIB_ITERATEDSOFTDROP_HH__
#define __FASTJET_CONTRIB_ITERATEDSOFTDROP_HH__

#include "RecursiveSoftDrop.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

//------------------------------------------------------------------------
/// \class IteratedSoftDropInfo
/// helper class that carries all the relevant information one can get
/// from running IteratedSoftDrop on a given jet (or vector of jets)
///
class IteratedSoftDropInfo{
public:
  /// ctor without initialisation
  IteratedSoftDropInfo(){}

  /// ctor with initialisation
  IteratedSoftDropInfo(std::vector<std::pair<double,double> > zg_thetag_in)
    : _all_zg_thetag(zg_thetag_in){}

  /// get the raw list of (angular-ordered) zg and thetag
  const std::vector<std::pair<double,double> > &all_zg_thetag() const{
    return _all_zg_thetag;
  }

  /// overloadd the () operator so that it also returns the full (zg,thetag) list
  const std::vector<std::pair<double,double> > & operator()() const{
    return _all_zg_thetag;
  }

  /// overloadd the [] operator to access the ith (zg,thetag) pair
  const std::pair<double,double> & operator[](unsigned int i) const{
    return _all_zg_thetag[i];
  }

  /// returns the angularity with angular exponent alpha and z
  /// exponent kappa calculated on the zg's and thetag's found by
  /// iterated SoftDrop
  ///
  /// returns 0 if no substructure was found
  double angularity(double alpha, double kappa=1.0) const;

  /// returns the Iterated SoftDrop multiplicity
  unsigned int multiplicity() const{ return _all_zg_thetag.size(); }  

  /// returns the Iterated SoftDrop multiplicity (i.e. size)
  unsigned int size() const{ return _all_zg_thetag.size(); }  
  
protected:
  /// the real information: angular-ordered list of all the zg and
  /// thetag that passed the (recursive) SD conddition
  std::vector<std::pair<double,double> > _all_zg_thetag;
};
  

  
//------------------------------------------------------------------------
/// \class IteratedSoftDrop
/// implementation of the IteratedSoftDrop procedure
/// 
/// This class provides an implementation of the IteratedSoftDrop
/// procedure.  It is based on the SoftDrop procedure can be used to
/// define a 'groomed symmetry factor', equal to the symmetry factor
/// of the two subjets of the resulting groomed jet.  The Iterated
/// Soft Drop procedure recursively performs Soft Drop on the harder
/// branch of the groomed jet, halting at a specified angular cut
/// \f$\theta_{\rm cut}\f$, returning a list of symmetry factors which
/// can be used to define observables.
///
/// Like SoftDrop, the cut applied recursively is 
///   \f[
///      z > z_{\rm cut} (\theta/R_0)^\beta
///   \f]
/// with z the asymmetry measure and \f$\theta\f$ the geometrical
/// distance between the two subjets. The procedure halts when
/// \f$\theta < \theta_{\rm cut}\f$.
///
/// By default, this implementation returs the IteratedSoftDropInfo
/// obtained after running IteratedSoftDrop on a jet
///
/// Although all these quantities can be obtained from the returned
/// IteratedSoftDropInfo, we also provide helpers to directly get the
/// multiplicity, some (generalised) angularity, or the raw list of
/// (angular-ordered) (zg, thetag) pairs that passed the (recursive)
/// SoftDrop condition.
///
/// We stress the fact that IteratedSoftDrop is _not_ a Transformer
/// since it returns an IteratedSoftDropInfo and not a modified
/// PseudoJet
///
class IteratedSoftDrop : public FunctionOfPseudoJet<IteratedSoftDropInfo> {
public:
  /// Constructor. Takes in the standard Soft Drop parameters, an angular cut \f$\theta_{\rm cut}\f$, 
  /// and a choice of angular and symmetry measure. 
  ///
  /// \param beta               the Soft Drop beta parameter
  /// \param symmetry_cut       the Soft Drop symmetry cut
  /// \param angular_cut        the angular cutoff to halt Iterated Soft Drop
  /// \param R0                 the angular distance normalization
  /// \param subtractor         an optional pointer to a pileup subtractor (ignored if zero)
  IteratedSoftDrop(double beta, double symmetry_cut, double angular_cut, double R0 = 1.0,
                   const FunctionOfPseudoJet<PseudoJet> * subtractor = 0);
  
  /// Full constructor, which takes the following parameters:
  ///
  /// \param beta               the value of the beta parameter
  /// \param symmetry_cut       the value of the cut on the symmetry measure
  /// \param symmetry_measure   the choice of measure to use to estimate the symmetry
  /// \param angular_cut        the angular cutoff to halt Iterated Soft Drop
  /// \param R0                 the angular distance normalisation [1 by default]
  /// \param mu_cut             the maximal allowed value of mass drop variable mu = m_heavy/m_parent 
  /// \param recursion_choice   the strategy used to decide which subjet to recurse into
  /// \param subtractor         an optional pointer to a pileup subtractor (ignored if zero)
  ///
  /// Notes: 
  ///
  /// - by default, SoftDrop will recluster the jet with the
  ///   Cambridge/Aachen algorithm if it is not already the case. This
  ///   behaviour can be changed using the "set_reclustering" method
  ///   defined below
  ///
  IteratedSoftDrop(double  beta,
                   double  symmetry_cut, 
                   RecursiveSoftDrop::SymmetryMeasure  symmetry_measure,
                   double  angular_cut,
                   double  R0 = 1.0,
                   double  mu_cut = std::numeric_limits<double>::infinity(), 
                   RecursiveSoftDrop::RecursionChoice  recursion_choice = RecursiveSoftDrop::larger_pt,
                   const FunctionOfPseudoJet<PseudoJet> * subtractor = 0);

  /// default destructor
  virtual ~IteratedSoftDrop(){}

  //----------------------------------------------------------------------
  // behaviour tweaks (inherited from RecursiveSoftDrop and RecursiveSymmetryCutBase)

  /// switch to using a dynamical R0 (see RecursiveSoftDrop)
  void set_dynamical_R0(bool value=true) { _rsd.set_dynamical_R0(value); }
  bool use_dynamical_R0() const { return _rsd.use_dynamical_R0(); }

  /// an alternative way to set the subtractor (see RecursiveSymmetryCutBase)
  void set_subtractor(const FunctionOfPseudoJet<PseudoJet> * subtractor_) {_rsd.set_subtractor(subtractor_);}
  const FunctionOfPseudoJet<PseudoJet> * subtractor() const {return _rsd.subtractor();}
  
  /// returns the IteratedSoftDropInfo associated with the jet "jet"
  IteratedSoftDropInfo result(const PseudoJet& jet) const;

  /// Tells the tagger whether to assume that the input jet has
  /// already been subtracted (relevant only with a non-null
  /// subtractor, see RecursiveSymmetryCutBase)
  void set_input_jet_is_subtracted(bool is_subtracted) { _rsd.set_input_jet_is_subtracted(is_subtracted);}
  bool input_jet_is_subtracted() const {return _rsd.input_jet_is_subtracted();}
  
  /// configure the reclustering prior to the recursive de-clustering
  void set_reclustering(bool do_reclustering=true, const Recluster *recluster=0){
    _rsd.set_reclustering(do_reclustering, recluster);
  }
  
  //----------------------------------------------------------------------
  // actions on jets
  /// returns vector of ISD symmetry factors and splitting angles
  std::vector<std::pair<double,double> > all_zg_thetag(const PseudoJet& jet) const{
    return result(jet).all_zg_thetag();
  }

  /// returns the angularity with angular exponent alpha and z
  /// exponent kappa calculated on the zg's and thetag's found by
  /// iterated SoftDrop
  ///
  /// returns 0 if no substructure was found
  double angularity(const PseudoJet& jet, double alpha, double kappa=1.0) const{
    return result(jet).angularity(alpha, kappa);
  }
  
  /// returns the Iterated SoftDrop multiplicity
  double multiplicity(const PseudoJet& jet) const{ return result(jet).multiplicity(); }

  /// description of the class
  std::string description() const;

protected:
  RecursiveSoftDrop _rsd;
};

  

} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_ITERATEDSOFTDROP_HH__
