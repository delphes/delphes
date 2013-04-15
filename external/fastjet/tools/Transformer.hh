//STARTHEADER
// $Id: Transformer.hh 2577 2011-09-13 15:11:38Z salam $
//
// Copyright (c) 2005-2011, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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
//  development and are described in hep-ph/0512210. If you use
//  FastJet as part of work towards a scientific publication, please
//  include a citation to the FastJet paper.
//
//  FastJet is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with FastJet. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------
//ENDHEADER

#ifndef __FASTJET_TRANSFORMER_HH__
#define __FASTJET_TRANSFORMER_HH__

#include <fastjet/PseudoJet.hh>
#include <fastjet/FunctionOfPseudoJet.hh>
#include <fastjet/PseudoJetStructureBase.hh>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

// forward declarations of what we will have down here
class Transformer;

/// @ingroup tools_generic
/// \class Transformer
///
/// Base (abstract) class for a jet transformer.
///
/// A transformer, when it acts on a jet, returns a modified version
/// of that jet, one that may have a different momentum and/or
/// different internal structure.
///
/// The typical usage of a class derived from Transformer is
/// \code
///   SomeTransformer transformer(...);
///   PseudoJet transformed_jet = transformer(original_jet);
///   // or
///   vector<PseudoJet> transformed_jets = transformer(original_jets);
/// \endcode
///
/// For many transformers, the transformed jets have
/// transformer-specific information that can be accessed through the
///
/// \code
///   transformed_jet.structure_of<SomeTransformer>().transformer_specific_info();
/// \endcode
///
/// See the description of the Filter class for a more detailed usage
/// example. See the FastJet manual to find out how to implement
/// new transformers.
///
class Transformer : public FunctionOfPseudoJet<PseudoJet>{
public:
  /// default ctor
  Transformer(){}

  /// default dtor
  virtual ~Transformer(){}

  /// the result of the Transformer acting on the PseudoJet.
  /// this _has_ to be overloaded in derived classes
  /// \param original   the PseudoJet input to the Transformer
  virtual PseudoJet result(const PseudoJet & original) const = 0;

  /// This should be overloaded to return a description of the
  /// Transformer
  virtual std::string description() const = 0;

  /// A typedef that is needed to ensure that the
  /// PseudoJet::structure_of() template function works
  //
  // Make sure you reimplement this appropriately in any
  // derived classes
  typedef PseudoJetStructureBase StructureType;
};

FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh

#endif // __FASTJET_TRANSFORMER_HH__
