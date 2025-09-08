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


#ifndef __FASTJET_CLUSTER_SEQUENCE_STRUCTURE_HH__
#define __FASTJET_CLUSTER_SEQUENCE_STRUCTURE_HH__

#include "fastjet/internal/base.hh"
#include "fastjet/SharedPtr.hh"
#include "fastjet/PseudoJetStructureBase.hh"

#include <vector>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

/// @ingroup extra_info
/// \class ClusterSequenceStructure
///
/// Contains any information related to the clustering that should be
/// directly accessible to PseudoJet.
///
/// By default, this class implements basic access to the
/// ClusterSequence related to a PseudoJet (like its constituents or
/// its area). But it can be overloaded in order e.g. to give access
/// to the jet substructure.
///
// Design question: Do we only put the methods that can be overloaded
// or do we put everything a PJ can have access to? I think both cost
// the same number of indirections. The first option limits the amount
// of coding and maybe has a clearer structure. The second is more
// consistent (everything related to the same thing is at the same
// place) and gives better access for derived classes. We'll go for
// the second option.
class ClusterSequenceStructure : public PseudoJetStructureBase{
public:
  /// default ctor
  ClusterSequenceStructure() : _associated_cs(NULL){}

  /// ctor with initialisation to a given ClusterSequence
  /// 
  /// In principle, this is reserved for initialisation by the parent
  /// ClusterSequence
  ClusterSequenceStructure(const ClusterSequence *cs){
    set_associated_cs(cs);
  };

  /// default (virtual) dtor
  virtual ~ClusterSequenceStructure();

  /// description
  virtual std::string description() const FASTJET_OVERRIDE{
    return "PseudoJet with an associated ClusterSequence";
  }

  //-------------------------------------------------------------
  /// @name Direct access to the associated ClusterSequence object.
  ///
  /// Get access to the associated ClusterSequence (if any)
  //\{
  //-------------------------------------------------------------
  /// returns true if there is an associated ClusterSequence
  virtual bool has_associated_cluster_sequence() const FASTJET_OVERRIDE{ return true;}

  /// get a (const) pointer to the parent ClusterSequence (NULL if
  /// inexistent)
  virtual const ClusterSequence* associated_cluster_sequence() const FASTJET_OVERRIDE;
  
  /// returns true if there is a valid associated ClusterSequence
  virtual bool has_valid_cluster_sequence() const FASTJET_OVERRIDE;

  /// if the jet has a valid associated cluster sequence then return a
  /// pointer to it; otherwise throw an error
  virtual const ClusterSequence * validated_cs() const FASTJET_OVERRIDE;

#ifndef __FASTJET_ONLY_CORE__
  /// if the jet has valid area information then return a pointer to
  /// the associated ClusterSequenceAreaBase object; otherwise throw an error
  virtual const ClusterSequenceAreaBase * validated_csab() const FASTJET_OVERRIDE;
#endif  // __FASTJET_ONLY_CORE__

  /// set the associated csw
  virtual void set_associated_cs(const ClusterSequence * new_cs){
    _associated_cs = new_cs;
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
  /// an Error is thrown if this PseudoJet has no currently valid
  /// associated ClusterSequence
  virtual bool has_partner(const PseudoJet &reference, PseudoJet &partner) const FASTJET_OVERRIDE;

  /// check if it has been recombined with another PseudoJet in which
  /// case, return its child through the argument. Otherwise, 'child'
  /// is set to 0.
  /// 
  /// an Error is thrown if this PseudoJet has no currently valid
  /// associated ClusterSequence
  virtual bool has_child(const PseudoJet &reference, PseudoJet &child) const FASTJET_OVERRIDE;

  /// check if it is the product of a recombination, in which case
  /// return the 2 parents through the 'parent1' and 'parent2'
  /// arguments. Otherwise, set these to 0.
  ///
  /// an Error is thrown if this PseudoJet has no currently valid
  /// associated ClusterSequence
  virtual bool has_parents(const PseudoJet &reference, PseudoJet &parent1, PseudoJet &parent2) const FASTJET_OVERRIDE;

  /// check if the reference PseudoJet is contained in the second one
  /// passed as argument.
  ///
  /// an Error is thrown if this PseudoJet has no currently valid
  /// associated ClusterSequence
  ///
  /// false is returned if the 2 PseudoJet do not belong the same
  /// ClusterSequence
  virtual bool object_in_jet(const PseudoJet &reference, const PseudoJet &jet) const FASTJET_OVERRIDE;

  /// return true if the structure supports constituents. 
  ///
  /// an Error is thrown if this PseudoJet has no currently valid
  /// associated ClusterSequence
  virtual bool has_constituents() const FASTJET_OVERRIDE;

  /// retrieve the constituents. 
  ///
  /// an Error is thrown if this PseudoJet has no currently valid
  /// associated ClusterSequence
  virtual std::vector<PseudoJet> constituents(const PseudoJet &reference) const FASTJET_OVERRIDE;


  /// return true if the structure supports exclusive_subjets. 
  ///
  /// an Error is thrown if this PseudoJet has no currently valid
  /// associated ClusterSequence
  virtual bool has_exclusive_subjets() const FASTJET_OVERRIDE;

  /// return a vector of all subjets of the current jet (in the sense
  /// of the exclusive algorithm) that would be obtained when running
  /// the algorithm with the given dcut. 
  ///
  /// Time taken is O(m ln m), where m is the number of subjets that
  /// are found. If m gets to be of order of the total number of
  /// constituents in the jet, this could be substantially slower than
  /// just getting that list of constituents.
  ///
  /// an Error is thrown if this PseudoJet has no currently valid
  /// associated ClusterSequence
  virtual std::vector<PseudoJet> exclusive_subjets(const PseudoJet &reference, const double & dcut) const FASTJET_OVERRIDE;

  /// return the size of exclusive_subjets(...); still n ln n with same
  /// coefficient, but marginally more efficient than manually taking
  /// exclusive_subjets.size()
  ///
  /// an Error is thrown if this PseudoJet has no currently valid
  /// associated ClusterSequence
  virtual int n_exclusive_subjets(const PseudoJet &reference, const double & dcut) const FASTJET_OVERRIDE;

  /// return the list of subjets obtained by unclustering the supplied
  /// jet down to nsub subjets (or all constituents if there are fewer
  /// than nsub).
  ///
  /// requires nsub ln nsub time
  ///
  /// an Error is thrown if this PseudoJet has no currently valid
  /// associated ClusterSequence
  virtual std::vector<PseudoJet> exclusive_subjets_up_to (const PseudoJet &reference, int nsub) const FASTJET_OVERRIDE;

  /// return the dij that was present in the merging nsub+1 -> nsub 
  /// subjets inside this jet.
  ///
  /// an Error is thrown if this PseudoJet has no currently valid
  /// associated ClusterSequence
  virtual double exclusive_subdmerge(const PseudoJet &reference, int nsub) const FASTJET_OVERRIDE;

  /// return the maximum dij that occurred in the whole event at the
  /// stage that the nsub+1 -> nsub merge of subjets occurred inside 
  /// this jet.
  ///
  /// an Error is thrown if this PseudoJet has no currently valid
  /// associated ClusterSequence
  virtual double exclusive_subdmerge_max(const PseudoJet &reference, int nsub) const FASTJET_OVERRIDE;


  //-------------------------------------------------------------------
  // information related to the pieces of the jet
  //-------------------------------------------------------------------
  /// by convention, a jet associated with a ClusterSequence will have
  /// its parents as pieces
  virtual bool has_pieces(const PseudoJet &reference) const FASTJET_OVERRIDE;

  /// by convention, a jet associated with a ClusterSequence will have
  /// its parents as pieces
  ///
  /// if it has no parents, then there will only be a single piece:
  /// itself
  ///
  /// Note that to answer that question, we need to access the cluster
  /// sequence. If the cluster sequence has gone out of scope, an
  /// error will be thrown
  virtual std::vector<PseudoJet> pieces(const PseudoJet &reference) const FASTJET_OVERRIDE;


  // the following ones require a computation of the area in the
  // parent ClusterSequence (See ClusterSequenceAreaBase for details)
  //------------------------------------------------------------------
#ifndef __FASTJET_ONLY_CORE__

  /// check if it has a defined area
  virtual bool has_area() const FASTJET_OVERRIDE;

  /// return the jet (scalar) area.
  /// throws an Error if there is no support for area in the parent CS
  virtual double area(const PseudoJet &reference) const FASTJET_OVERRIDE;

  /// return the error (uncertainty) associated with the determination
  /// of the area of this jet.
  /// throws an Error if there is no support for area in the parent CS
  virtual double area_error(const PseudoJet &reference) const FASTJET_OVERRIDE;

  /// return the jet 4-vector area.
  /// throws an Error if there is no support for area in the parent CS
  virtual PseudoJet area_4vector(const PseudoJet &reference) const FASTJET_OVERRIDE;

  /// true if this jet is made exclusively of ghosts.
  /// throws an Error if there is no support for area in the parent CS
  virtual bool is_pure_ghost(const PseudoJet &reference) const FASTJET_OVERRIDE;

#endif  // __FASTJET_ONLY_CORE__
  //\} --- end of jet structure -------------------------------------

protected:
  const ClusterSequence *_associated_cs;
};

FASTJET_END_NAMESPACE

#endif  //  __FASTJET_CLUSTER_SEQUENCE_STRUCTURE_HH__
