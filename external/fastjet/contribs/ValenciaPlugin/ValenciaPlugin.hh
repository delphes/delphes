// $Id: ValenciaPlugin.hh 771 2015-02-21 16:40:07Z vos $
//
// Copyright (c) 2014, Marcel Vos and Ignacio Garcia 
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

#ifndef __FASTJET_CONTRIB_VALENCIAJETALGORITHM_HH__
#define __FASTJET_CONTRIB_VALENCIAJETALGORITHM_HH__

#include <fastjet/internal/base.hh>
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

//----------------------------------------------------------------------
//
/// ValenciaPlugin is a plugin for fastjet (v2.4 upwards)
///
/// It implements the Valencia algorithm, as defined in 
/// Boronat, Garcia, Vos, 
/// A new jet reconstruction algorithm for lepton colliders
/// 
/// 
class ValenciaPlugin : public JetDefinition::Plugin {
public:

  /// Constructor for the Valencia Plugin class.  
  /// Three floating point arguments are specified to set the parameters
  /// the radius parameter R has the usual meaning,
  /// the clustering order beta (beta = 1 yields kt-style clustering,
  /// beta = 0 purely angular clustering a la C/A and beta = -1
  /// clusters hard, collinear radiation first, like anti-kt),
  /// and gamma, that governs the shrinking jet size in the forward region
  ValenciaPlugin (double R, double beta, double gamma) : _R(R), _beta(beta), _gamma(gamma){}

  /// Constructor for the Valencia Plugin class.  
  /// If two arguments are specified to set the parameters
  /// the gamma exponent is set equal to beta
  ValenciaPlugin (double R, double beta) : _R(R), _beta(beta), _gamma(beta){}

  /// copy constructor
  ValenciaPlugin (const ValenciaPlugin & plugin) {
    *this = plugin;
  }

  // the things that are required by base class
  virtual std::string description () const;
  virtual void run_clustering(ClusterSequence &) const;

  /// the plugin mechanism's standard way of accessing the jet radius.
  /// This must be set to return something sensible, even if R
  /// does not make sense for this algorithm!
  virtual double R() const {return _R;}
  
  // the Valencia algorithm has a second parameter beta that governs
  // the exponent of the energy in the inter-particle and beam distance 
  // criteria, and thus determines the clustering order
  virtual double beta() const {return _beta;}
  
  // the Valencia algorithm has a third parameter gamma that governs
  // the exponent of the sin(theta) in the beam distance
  // and thus the shrinking of the jet size in the forward region
   virtual double gamma() const {return _gamma;}
 
  

  /// avoid the warning whenever the user requests "exclusive" jets
  /// from the cluster sequence
  virtual bool exclusive_sequence_meaningful() const {return true;}

private:
  double _R;
  double _beta;
  double _gamma;
};




} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_VALENCIAJETALGORITHM_HH__
