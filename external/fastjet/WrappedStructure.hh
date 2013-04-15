//STARTHEADER
// $Id: WrappedStructure.hh 2577 2011-09-13 15:11:38Z salam $
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


#ifndef __FASTJET_WRAPPED_STRUCTURE_HH__
#define __FASTJET_WRAPPED_STRUCTURE_HH__

#include "fastjet/PseudoJetStructureBase.hh"
#include "fastjet/Error.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

/// @ingroup extra_info
/// \class WrappedStructure
///
/// This wraps a (shared) pointer to an underlying structure
///
/// The typical use-case is when a PseusoJet needs to share its
/// structure with another PseudoJet but also include extra
/// information in its structure. For the memory management to be
/// handled properly, it should hold a shared pointer to the shared
/// structure. This is what this class ensures. Deriving a structure
/// from this class would then allow for the implementation of the
/// extra features.
///
class WrappedStructure : public PseudoJetStructureBase{
public:
  /// default ctor
  /// the argument is the structure we need to wrap
  WrappedStructure(const SharedPtr<PseudoJetStructureBase> & to_be_shared)
    : _structure(to_be_shared){
    if (!_structure())
      throw Error("Trying to construct a wrapped structure around an empty (NULL) structure");
  }

  /// default (virtual) dtor
  virtual ~WrappedStructure(){}

  /// description
  virtual std::string description() const{ 
    return "PseudoJet wrapping the structure ("+_structure->description()+")"; 
  }

  //-------------------------------------------------------------
  /// @name Direct access to the associated ClusterSequence object.
  ///
  /// Get access to the associated ClusterSequence (if any)
  //\{
  //-------------------------------------------------------------
  /// returns true if there is an associated ClusterSequence
  virtual bool has_associated_cluster_sequence() const {
    return _structure->has_associated_cluster_sequence();
  }

  /// get a (const) pointer to the parent ClusterSequence (NULL if
  /// inexistent)
  virtual const ClusterSequence* associated_cluster_sequence() const{
    return _structure->associated_cluster_sequence();
  }
  
  /// returns true if this PseudoJet has an associated and still
  /// valid ClusterSequence.
  virtual bool has_valid_cluster_sequence() const {
    return _structure->has_valid_cluster_sequence();
  }

  /// if the jet has a valid associated cluster sequence then return a
  /// pointer to it; otherwise throw an error
  virtual const ClusterSequence * validated_cs() const{
    return _structure->validated_cs();
  }

  /// if the jet has valid area information then return a pointer to
  /// the associated ClusterSequenceAreaBase object; otherwise throw an error
  virtual const ClusterSequenceAreaBase * validated_csab() const{
    return _structure->validated_csab();
  }

  //\}

  //-------------------------------------------------------------
  /// @name Methods for access to information about jet structure
  ///
  /// These allow access to jet constituents, and other jet
  /// subtructure information. They only work if the jet is associated
  /// with a ClusterSequence.
  //-------------------------------------------------------------
  //\{

  /// check if it has been recombined with another PseudoJet in which
  /// case, return its partner through the argument. Otherwise,
  /// 'partner' is set to 0.
  ///
  /// By default, throws an Error
  virtual bool has_partner(const PseudoJet &reference, PseudoJet &partner) const{
    return _structure->has_partner(reference, partner);
  }

  /// check if it has been recombined with another PseudoJet in which
  /// case, return its child through the argument. Otherwise, 'child'
  /// is set to 0.
  /// 
  /// By default, throws an Error
  virtual bool has_child(const PseudoJet &reference, PseudoJet &child) const{
    return _structure->has_child(reference, child);
  }

  /// check if it is the product of a recombination, in which case
  /// return the 2 parents through the 'parent1' and 'parent2'
  /// arguments. Otherwise, set these to 0.
  ///
  /// By default, throws an Error
  virtual bool has_parents(const PseudoJet &reference, PseudoJet &parent1, PseudoJet &parent2) const{
    return _structure->has_parents(reference, parent1, parent2);
  }

  /// check if the reference PseudoJet is contained the second one
  /// passed as argument.
  ///
  /// By default, throws an Error
  virtual bool object_in_jet(const PseudoJet &reference, const PseudoJet &jet) const{
    return _structure->object_in_jet(reference, jet);
  }


  /// return true if the structure supports constituents. 
  ///
  /// false by default
  virtual bool has_constituents() const {
    return _structure->has_constituents();
  }

  /// retrieve the constituents. 
  ///
  /// By default, throws an Error
  virtual std::vector<PseudoJet> constituents(const PseudoJet &reference) const{
    return _structure->constituents(reference);
  }

  /// return true if the structure supports exclusive_subjets. 
  virtual bool has_exclusive_subjets() const {
    return _structure->has_exclusive_subjets();
  }

  /// return a vector of all subjets of the current jet (in the sense
  /// of the exclusive algorithm) that would be obtained when running
  /// the algorithm with the given dcut. 
  ///
  /// Time taken is O(m ln m), where m is the number of subjets that
  /// are found. If m gets to be of order of the total number of
  /// constituents in the jet, this could be substantially slower than
  /// just getting that list of constituents.
  ///
  /// By default, throws an Error
  virtual std::vector<PseudoJet> exclusive_subjets(const PseudoJet &reference, const double & dcut) const{
    return _structure->exclusive_subjets(reference, dcut);
  }

  /// return the size of exclusive_subjets(...); still n ln n with same
  /// coefficient, but marginally more efficient than manually taking
  /// exclusive_subjets.size()
  ///
  /// By default, throws an Error
  virtual int n_exclusive_subjets(const PseudoJet &reference, const double & dcut) const{
    return _structure->n_exclusive_subjets(reference, dcut);
  }

  /// return the list of subjets obtained by unclustering the supplied
  /// jet down to n subjets (or all constituents if there are fewer
  /// than n).
  ///
  /// By default, throws an Error
  virtual std::vector<PseudoJet> exclusive_subjets_up_to (const PseudoJet &reference, int nsub) const{
    return _structure->exclusive_subjets_up_to (reference, nsub);
  }

  /// return the dij that was present in the merging nsub+1 -> nsub 
  /// subjets inside this jet.
  ///
  /// By default, throws an Error
  virtual double exclusive_subdmerge(const PseudoJet &reference, int nsub) const{
    return _structure->exclusive_subdmerge(reference, nsub);
  }

  /// return the maximum dij that occurred in the whole event at the
  /// stage that the nsub+1 -> nsub merge of subjets occurred inside 
  /// this jet.
  ///
  /// By default, throws an Error
  virtual double exclusive_subdmerge_max(const PseudoJet &reference, int nsub) const{
    return _structure->exclusive_subdmerge_max(reference, nsub);
  }


  //-------------------------------------------------------------------
  // information related to the pieces of the jet
  //-------------------------------------------------------------------
  /// return true if the structure supports pieces. 
  ///
  /// false by default
  virtual bool has_pieces(const PseudoJet &reference) const {
    return _structure->has_pieces(reference);
  }

  /// retrieve the pieces building the jet. 
  ///
  /// By default, throws an Error
  virtual std::vector<PseudoJet> pieces(const PseudoJet &reference) const{
    return _structure->pieces(reference);
  }

  // the following ones require a computation of the area in the
  // parent ClusterSequence (See ClusterSequenceAreaBase for details)
  //------------------------------------------------------------------

  /// check if it has a defined area
  ///
  /// false by default
  virtual bool has_area() const {
    return _structure->has_area();
  }

  /// return the jet (scalar) area.
  ///
  /// By default, throws an Error
  virtual double area(const PseudoJet &reference) const{
    return _structure->area(reference);
  }

  /// return the error (uncertainty) associated with the determination
  /// of the area of this jet.
  ///
  /// By default, throws an Error
  virtual double area_error(const PseudoJet &reference) const{
    return _structure->area_error(reference);
  }

  /// return the jet 4-vector area.
  ///
  /// By default, throws an Error
  virtual PseudoJet area_4vector(const PseudoJet &reference) const{
    return _structure->area_4vector(reference);
  }

  /// true if this jet is made exclusively of ghosts.
  ///
  /// By default, throws an Error
  virtual bool is_pure_ghost(const PseudoJet &reference) const{
    return _structure->is_pure_ghost(reference);
  }

  //\} --- end of jet structure -------------------------------------

protected:
  SharedPtr<PseudoJetStructureBase> _structure;  ///< the wrapped structure
};

FASTJET_END_NAMESPACE

#endif  //  __FASTJET_PSEUDOJET_STRUCTURE_BASE_HH__
