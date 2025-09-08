//FJSTARTHEADER
// $Id$
//
// Copyright (c) 2005-2025, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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


#ifndef __FASTJET_COMPOSITEJET_STRUCTURE_HH__
#define __FASTJET_COMPOSITEJET_STRUCTURE_HH__

#include "fastjet/PseudoJet.hh"
#include "fastjet/PseudoJetStructureBase.hh"

// to have access to the recombiner we need to include the JetDefinition header
#include "fastjet/JetDefinition.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

/// @ingroup extra_info
/// \class CompositeJetStructure
/// The structure for a jet made of pieces
///
/// This stores the vector of the pieces that make the jet and provide
/// the methods to access them
class CompositeJetStructure : public PseudoJetStructureBase{
public:
  // basic class info
  //-------------------------------------------------------------------
  /// default ctor
  CompositeJetStructure() : _area_4vector_ptr(0){};

  /// ctor with initialisation
  CompositeJetStructure(const std::vector<PseudoJet> & initial_pieces, 
			const JetDefinition::Recombiner * recombiner = 0);

  /// default dtor
  virtual ~CompositeJetStructure(){
    if (_area_4vector_ptr) delete _area_4vector_ptr;
  };

  /// description
  virtual std::string description() const FASTJET_OVERRIDE;

  // things reimplemented from the base structure
  //-------------------------------------------------------------------
  /// true unless the jet has no pieces (see also the description of
  /// constituents() below)
  virtual bool has_constituents() const FASTJET_OVERRIDE;

  /// return the constituents (i.e. the union of the constituents of each piece)
  /// 
  /// If any of the pieces has no constituent, the piece itself is
  /// considered as a constituent
  /// Note that as a consequence, a composite jet with no pieces will
  /// have an empty vector as constituents
  virtual std::vector<PseudoJet> constituents(const PseudoJet &jet) const FASTJET_OVERRIDE;

  //-------------------------------------------------------------------
  // information related to the pieces of the jet
  //-------------------------------------------------------------------
  /// true if it has pieces (always the case)
  virtual bool has_pieces(const PseudoJet & /*jet*/) const FASTJET_OVERRIDE {return true;}

  /// returns the pieces
  virtual std::vector<PseudoJet> pieces(const PseudoJet &jet) const FASTJET_OVERRIDE;

  // area-related material
#ifndef __FASTJET_ONLY_CORE__

  /// check if it has a well-defined area
  virtual bool has_area() const FASTJET_OVERRIDE;

  /// return the jet (scalar) area.
  virtual double area(const PseudoJet &reference) const FASTJET_OVERRIDE;

  /// return the error (uncertainty) associated with the determination
  /// of the area of this jet.
  ///
  /// Be conservative: return the sum of the errors
  virtual double area_error(const PseudoJet &reference) const FASTJET_OVERRIDE;

  /// return the jet 4-vector area.
  virtual PseudoJet area_4vector(const PseudoJet &reference) const FASTJET_OVERRIDE;

  /// true if this jet is made exclusively of ghosts.
  ///
  /// In this case, it will be true if all pieces are pure ghost
  virtual bool is_pure_ghost(const PseudoJet &reference) const FASTJET_OVERRIDE;

  //unused: // allows one to modify the area information
  //unused: // (for use in join())
  //unused: //
  //unused: // This member cannot be used by users who need to create a jet with
  //unused: // user-supplied area information, because it sets only the 4-vector
  //unused: // part of the area, but not all the other area information
  //unused: // (e.g. scalar area) -- that other information is always deduced
  //unused: // dynamically from the individual constituents.
  //unused: // ------------------------------------------------------------------------------
  //unused: void set_area_information(PseudoJet *area_4vector_ptr){
  //unused:   _area_4vector_ptr = area_4vector_ptr;
  //unused: }

  /// disable the area of the composite jet
  /// 
  /// this can be used e.g. to discard the area of a composite jet
  /// made of pieces with non-explicit-ghost area since the area may
  /// by erroneous in that case
  void discard_area(){
    if (_area_4vector_ptr) delete _area_4vector_ptr;
    _area_4vector_ptr = 0;
  }

#endif  // __FASTJET_ONLY_CORE__

protected:
  std::vector<PseudoJet> _pieces;  ///< the pieces building the jet
  PseudoJet * _area_4vector_ptr;   ///< pointer to the 4-vector jet area
};



// helpers to "join" jets and produce a structure derived from
// CompositeJetStructure
//
// The template structure T must have a constructor accepting as
// argument the pieces and of the composite jet
// ------------------------------------------------------------------------

/// build a "CompositeJet" from the vector of its pieces with an
/// extended structure of type T derived from CompositeJetStructure
template<typename T> PseudoJet join(const std::vector<PseudoJet> & pieces){
  PseudoJet result(0.0,0.0,0.0,0.0);
  for (unsigned int i=0; i<pieces.size(); i++){
    const PseudoJet it = pieces[i];
    result += it;
  }

  T *cj_struct = new T(pieces);
  result.set_structure_shared_ptr(SharedPtr<PseudoJetStructureBase>(cj_struct));

  return result;
}


/// build a "CompositeJet" from a single PseudoJet with an extended
/// structure of type T derived from CompositeJetStructure
template<typename T> PseudoJet join(const PseudoJet & j1){
  return join<T>(std::vector<PseudoJet>(1,j1));
}

/// build a "CompositeJet" from two PseudoJet with an extended
/// structure of type T derived from CompositeJetStructure
template<typename T> PseudoJet join(const PseudoJet & j1, const PseudoJet & j2){
  std::vector<PseudoJet> pieces;
  pieces.push_back(j1);
  pieces.push_back(j2);
  return join<T>(pieces);
}

/// build a "CompositeJet" from 3 PseudoJet with an extended structure
/// of type T derived from CompositeJetStructure
template<typename T> PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, 
				    const PseudoJet & j3){
  std::vector<PseudoJet> pieces;
  pieces.push_back(j1);
  pieces.push_back(j2);
  pieces.push_back(j3);
  return join<T>(pieces);
}

/// build a "CompositeJet" from 4 PseudoJet with an extended structure
/// of type T derived from CompositeJetStructure
template<typename T> PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, 
				    const PseudoJet & j3, const PseudoJet & j4){
  std::vector<PseudoJet> pieces;
  pieces.push_back(j1);
  pieces.push_back(j2);
  pieces.push_back(j3);
  pieces.push_back(j4);
  return join<T>(pieces);
}


// the same as above with an additional argument for a
// user-defined recombiner
//
// The template structure T must be derived from CompositeJetStructure
// and have a constructor accepting as arguments the pieces and a
// pointer to the recombination scheme
// ----------------------------------------------------------------------

/// build a "CompositeJet" from the vector of its pieces with an
/// extended structure of type T derived from CompositeJetStructure
template<typename T> PseudoJet join(const std::vector<PseudoJet> & pieces, 
				    const JetDefinition::Recombiner & recombiner){
  PseudoJet result;
  if (pieces.size()>0){
    result = pieces[0];
    for (unsigned int i=1; i<pieces.size(); i++){
      recombiner.plus_equal(result, pieces[i]);
    }
  }

  T *cj_struct = new T(pieces, &recombiner);
  result.set_structure_shared_ptr(SharedPtr<PseudoJetStructureBase>(cj_struct));

  return result;
}

/// build a "CompositeJet" from a single PseudoJet with an extended
/// structure of type T derived from CompositeJetStructure
template<typename T> PseudoJet join(const PseudoJet & j1, 
				    const JetDefinition::Recombiner & recombiner){
  return join<T>(std::vector<PseudoJet>(1,j1), recombiner);
}

/// build a "CompositeJet" from two PseudoJet with an extended
/// structure of type T derived from CompositeJetStructure
template<typename T> PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, 
				    const JetDefinition::Recombiner & recombiner){
  std::vector<PseudoJet> pieces;
  pieces.reserve(2);
  pieces.push_back(j1);
  pieces.push_back(j2);
  return join<T>(pieces, recombiner);
}

/// build a "CompositeJet" from 3 PseudoJet with an extended structure
/// of type T derived from CompositeJetStructure
template<typename T> PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, 
				    const PseudoJet & j3, 
				    const JetDefinition::Recombiner & recombiner){
  std::vector<PseudoJet> pieces;
  pieces.reserve(3);
  pieces.push_back(j1);
  pieces.push_back(j2);
  pieces.push_back(j3);
  return join<T>(pieces, recombiner);
}

/// build a "CompositeJet" from 4 PseudoJet with an extended structure
/// of type T derived from CompositeJetStructure
template<typename T> PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, 
				    const PseudoJet & j3, const PseudoJet & j4, 
				    const JetDefinition::Recombiner & recombiner){
  std::vector<PseudoJet> pieces;
  pieces.reserve(4);
  pieces.push_back(j1);
  pieces.push_back(j2);
  pieces.push_back(j3);
  pieces.push_back(j4);
  return join<T>(pieces, recombiner);
}


FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh

#endif // __FASTJET_MERGEDJET_STRUCTURE_HH__
