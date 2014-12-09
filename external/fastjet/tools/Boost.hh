//FJSTARTHEADER
// $Id: Boost.hh 3433 2014-07-23 08:17:03Z salam $
//
// Copyright (c) 2005-2014, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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
//  development. They are described in the original FastJet paper,
//  hep-ph/0512210 and in the manual, arXiv:1111.6097. If you use
//  FastJet as part of work towards a scientific publication, please
//  quote the version you use and include a citation to the manual and
//  optionally also to hep-ph/0512210.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//FJENDHEADER

#ifndef __FASTJET_TOOL_BOOST_HH__
#define __FASTJET_TOOL_BOOST_HH__

#include <fastjet/PseudoJet.hh>
#include <fastjet/FunctionOfPseudoJet.hh>
#include <fastjet/PseudoJetStructureBase.hh>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

/// @ingroup tools_generic
/// \class Boost
/// Class to boost a PseudoJet
///
/// This is a FunctionOfPseudoJet with return type PseudoJet. Its
/// action if to boost the PseudoJet by a boost vector passed to its
/// constructor
class Boost : public FunctionOfPseudoJet<PseudoJet>{
public:
  /// default ctor
  Boost(const PseudoJet & jet_rest) : _jet_rest(jet_rest){}

  /// the action of the function: boost the PseudoJet by a boost
  /// vector _jet_rest
  PseudoJet result(const PseudoJet & original) const{
    PseudoJet res = original;
    return res.boost(_jet_rest);
  }

protected:
  PseudoJet _jet_rest;  ///< the boost vector
};

/// @ingroup tools_generic
/// \class Unboost
/// Class to un-boost a PseudoJet
///
/// This is a FunctionOfPseudoJet with return type PseudoJet. Its
/// action if to un-boost the PseudoJet back in the restframe of the
/// PseudoJet passed to its constructor
class Unboost : public FunctionOfPseudoJet<PseudoJet>{
public:
  /// default ctor
  Unboost(const PseudoJet & jet_rest) : _jet_rest(jet_rest){}

  /// the action of the function: boost the PseudoJet to the rest
  /// frame of _jet_rest
  PseudoJet result(const PseudoJet & original) const{
    PseudoJet res = original;
    return res.unboost(_jet_rest);
  }

protected:
  PseudoJet _jet_rest;  ///< the boost vector
};

FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh

#endif // __FASTJET_TRANSFORMER_HH__
