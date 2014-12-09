#ifndef __FASTJET_PSEUDOJET_STRUCTURE_BASE_HH__
#define __FASTJET_PSEUDOJET_STRUCTURE_BASE_HH__

//FJSTARTHEADER
// $Id: PseudoJetStructureBase.hh 3433 2014-07-23 08:17:03Z salam $
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


#include "fastjet/internal/base.hh"

#include <vector>
#include <string>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

class PseudoJet;
class ClusterSequence;
#ifndef __FJCORE__
class ClusterSequenceAreaBase;
#endif  // __FJCORE__

/// @ingroup extra_info
/// \class PseudoJetStructureBase
///
/// Contains any information related to the clustering that should be
/// directly accessible to PseudoJet.
///
/// By default, this class implements basic access to the
/// ClusterSequence related to a PseudoJet (like its constituents or
/// its area). But it can be overloaded in order e.g. to give access
/// to the jet substructure.
///
class PseudoJetStructureBase{
public:
  /// default ctor
  PseudoJetStructureBase(){};

  /// default (virtual) dtor
  virtual ~PseudoJetStructureBase(){};

  /// description
  virtual std::string description() const{ return "PseudoJet with an unknown structure"; }

  //-------------------------------------------------------------
  /// @name Direct access to the associated ClusterSequence object.
  ///
  /// Get access to the associated ClusterSequence (if any)
  //\{
  //-------------------------------------------------------------
  /// returns true if there is an associated ClusterSequence
  virtual bool has_associated_cluster_sequence() const { return false;}

  /// get a (const) pointer to the parent ClusterSequence (NULL if
  /// inexistent)
  virtual const ClusterSequence* associated_cluster_sequence() const;
  
  /// returns true if this PseudoJet has an associated and still
  /// valid ClusterSequence.
  virtual bool has_valid_cluster_sequence() const {return false;}

  /// if the jet has a valid associated cluster sequence then return a
  /// pointer to it; otherwise throw an error
  virtual const ClusterSequence * validated_cs() const;

#ifndef __FJCORE__
  /// if the jet has valid area information then return a pointer to
  /// the associated ClusterSequenceAreaBase object; otherwise throw an error
  virtual const ClusterSequenceAreaBase * validated_csab() const;
#endif

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
  virtual bool has_partner(const PseudoJet &reference, PseudoJet &partner) const;

  /// check if it has been recombined with another PseudoJet in which
  /// case, return its child through the argument. Otherwise, 'child'
  /// is set to 0.
  /// 
  /// By default, throws an Error
  virtual bool has_child(const PseudoJet &reference, PseudoJet &child) const;

  /// check if it is the product of a recombination, in which case
  /// return the 2 parents through the 'parent1' and 'parent2'
  /// arguments. Otherwise, set these to 0.
  ///
  /// By default, throws an Error
  virtual bool has_parents(const PseudoJet &reference, PseudoJet &parent1, PseudoJet &parent2) const;

  /// check if the reference PseudoJet is contained the second one
  /// passed as argument.
  ///
  /// By default, throws an Error
  virtual bool object_in_jet(const PseudoJet &reference, const PseudoJet &jet) const;


  /// return true if the structure supports constituents. 
  ///
  /// false by default
  virtual bool has_constituents() const {return false;}

  /// retrieve the constituents. 
  ///
  /// By default, throws an Error
  virtual std::vector<PseudoJet> constituents(const PseudoJet &reference) const;


  /// return true if the structure supports exclusive_subjets. 
  virtual bool has_exclusive_subjets() const {return false;}

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
  ///
  /// Note: in a future major release of FastJet (4 or higher), 
  /// "const double & dcut" may be replaced with "const double dcut",
  /// requiring a modification of derived classes that overload
  /// this function.
  virtual std::vector<PseudoJet> exclusive_subjets(const PseudoJet &reference, const double & dcut) const;

  /// return the size of exclusive_subjets(...); still n ln n with same
  /// coefficient, but marginally more efficient than manually taking
  /// exclusive_subjets.size()
  ///
  /// By default, throws an Error
  ///
  /// Note: in a future major release of FastJet (4 or higher), 
  /// "const double & dcut" may be replaced with "const double dcut",
  /// requiring a modification of derived classes that overload
  /// this function.
  virtual int n_exclusive_subjets(const PseudoJet &reference, const double & dcut) const;

  /// return the list of subjets obtained by unclustering the supplied
  /// jet down to nsub subjets (or all constituents if there are fewer
  /// than nsub).
  ///
  /// By default, throws an Error
  virtual std::vector<PseudoJet> exclusive_subjets_up_to (const PseudoJet &reference, int nsub) const;

  /// return the dij that was present in the merging nsub+1 -> nsub 
  /// subjets inside this jet.
  ///
  /// By default, throws an Error
  virtual double exclusive_subdmerge(const PseudoJet &reference, int nsub) const;

  /// return the maximum dij that occurred in the whole event at the
  /// stage that the nsub+1 -> nsub merge of subjets occurred inside 
  /// this jet.
  ///
  /// By default, throws an Error
  virtual double exclusive_subdmerge_max(const PseudoJet &reference, int nsub) const;


  //-------------------------------------------------------------------
  // information related to the pieces of the jet
  //-------------------------------------------------------------------
  /// return true if the structure supports pieces. 
  ///
  /// false by default
  /// NB: "reference" is commented to avoid unused-variable compiler warnings
  virtual bool has_pieces(const PseudoJet & /* reference */) const {
    return false;}

  /// retrieve the pieces building the jet. 
  ///
  /// By default, throws an Error.
  /// NB: "reference" is commented to avoid unused-variable compiler warnings
  virtual std::vector<PseudoJet> pieces(const PseudoJet & /* reference */
                                        ) const;


  // the following ones require a computation of the area in the
  // parent ClusterSequence (See ClusterSequenceAreaBase for details)
  //------------------------------------------------------------------
#ifndef __FJCORE__

  /// check if it has a defined area
  ///
  /// false by default
  virtual bool has_area() const {return false;}

  /// return the jet (scalar) area.
  ///
  /// By default, throws an Error
  virtual double area(const PseudoJet &reference) const;

  /// return the error (uncertainty) associated with the determination
  /// of the area of this jet.
  ///
  /// By default, throws an Error
  virtual double area_error(const PseudoJet &reference) const;

  /// return the jet 4-vector area.
  ///
  /// By default, throws an Error
  virtual PseudoJet area_4vector(const PseudoJet &reference) const;

  /// true if this jet is made exclusively of ghosts.
  ///
  /// By default, throws an Error
  virtual bool is_pure_ghost(const PseudoJet &reference) const;

#endif  // __FJCORE__
  //\} --- end of jet structure -------------------------------------
};

FASTJET_END_NAMESPACE

#endif  //  __FASTJET_PSEUDOJET_STRUCTURE_BASE_HH__
