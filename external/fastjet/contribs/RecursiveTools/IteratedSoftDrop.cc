// $Id: IteratedSoftDrop.cc 1084 2017-10-10 20:36:50Z gsoyez $
//
// Copyright (c) 2017-, Jesse Thaler, Kevin Zhou, Gavin P. Salam
// andGregory Soyez
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

#include "IteratedSoftDrop.hh"
#include <sstream>
#include <cmath>
#include <fastjet/ClusterSequence.hh>

using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

//========================================================================
// implementation of IteratedSoftDropInfo
//========================================================================

// returns the angularity with angular exponent alpha and z
// exponent kappa calculated on the zg's and thetag's found by
// iterated SoftDrop
//
// returns 0 if no substructure was found
double IteratedSoftDropInfo::angularity(double alpha, double kappa) const{
  double sum = 0.0;
  for (unsigned int i=0; i< _all_zg_thetag.size(); ++i)
    sum += pow(_all_zg_thetag[i].first, kappa) * pow(_all_zg_thetag[i].second, alpha);
  return sum;    
}
  
//========================================================================
// implementation of IteratedSoftDrop
//========================================================================

// Constructor. Takes in the standard Soft Drop parameters, an angular cut \f$\theta_{\rm cut}\f$, 
// and a choice of angular and symmetry measure. 
//
//  - beta           the Soft Drop beta parameter
//  - symmetry_cut   the Soft Drop symmetry cut
//  - angular_cut    the angular cutoff to halt Iterated Soft Drop
//  - R0             the angular distance normalization
IteratedSoftDrop::IteratedSoftDrop(double beta, double symmetry_cut, double angular_cut, double R0,
                                   const FunctionOfPseudoJet<PseudoJet> * subtractor) :
  _rsd(beta, symmetry_cut, -1, R0, subtractor){
  _rsd.set_hardest_branch_only(true);
  if (angular_cut>0)
    _rsd.set_min_deltaR_squared(angular_cut*angular_cut);
}


// Full constructor, which takes the following parameters:
//
// \param beta               the value of the beta parameter
// \param symmetry_cut       the value of the cut on the symmetry measure
// \param symmetry_measure   the choice of measure to use to estimate the symmetry
// \param angular_cut        the angular cutoff to halt Iterated Soft Drop
// \param R0                 the angular distance normalisation [1 by default]
// \param mu_cut             the maximal allowed value of mass drop variable mu = m_heavy/m_parent 
// \param recursion_choice   the strategy used to decide which subjet to recurse into
// \param subtractor         an optional pointer to a pileup subtractor (ignored if zero)
IteratedSoftDrop::IteratedSoftDrop(double  beta,
                                   double  symmetry_cut, 
                                   RecursiveSoftDrop::SymmetryMeasure  symmetry_measure,
                                   double  angular_cut,
                                   double  R0,
                                   double  mu_cut, 
                                   RecursiveSoftDrop::RecursionChoice  recursion_choice,
                                   const FunctionOfPseudoJet<PseudoJet> * subtractor)
  : _rsd(beta, symmetry_cut, symmetry_measure, -1, R0, mu_cut, recursion_choice, subtractor){
  _rsd.set_hardest_branch_only(true);
  if (angular_cut>0)
    _rsd.set_min_deltaR_squared(angular_cut*angular_cut);
}

  
// returns vector of ISD symmetry factors and splitting angles
IteratedSoftDropInfo IteratedSoftDrop::result(const PseudoJet& jet) const{
  PseudoJet rsd_jet = _rsd(jet);
  if (! rsd_jet.has_structure_of<RecursiveSoftDrop>())
    return IteratedSoftDropInfo();
  return IteratedSoftDropInfo(rsd_jet.structure_of<RecursiveSoftDrop>().sorted_zg_and_thetag());
}


std::string IteratedSoftDrop::description() const{
  std::ostringstream oss;
  oss << "IteratedSoftDrop with beta =" << _rsd.beta()
      << ", symmetry_cut=" << _rsd.symmetry_cut()
      << ", R0=" << _rsd.R0();
  
  if (_rsd.min_deltaR_squared() >= 0){
    oss << " and angular_cut=" << sqrt(_rsd.min_deltaR_squared());
  } else {
    oss << " and no angular_cut";
  }
  
  if (_rsd.subtractor()){
    oss << ", and with internal subtraction using [" << _rsd.subtractor()->description() << "]";
  }
  return oss.str();
}


} // namespace contrib

FASTJET_END_NAMESPACE
