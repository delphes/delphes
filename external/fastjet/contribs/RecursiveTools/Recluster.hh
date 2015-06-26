#ifndef __FASTJET_TOOLS_RECLUSTER_HH__
#define __FASTJET_TOOLS_RECLUSTER_HH__

// $Id: Recluster.hh 700 2014-07-07 12:50:05Z gsoyez $
//
// Copyright (c) 2014-, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
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

#include <fastjet/JetDefinition.hh>
#include <fastjet/CompositeJetStructure.hh> // to derive the ReclusterStructure from CompositeJetStructure
#include <fastjet/tools/Transformer.hh>     // to derive Recluster from Transformer
#include <iostream>
#include <string>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

//----------------------------------------------------------------------
/// \class Recluster
/// Class that helps reclustering a jet with a new jet definition
///
/// The result of the reclustering is returned as a single PseudoJet
/// with a CompositeJet structure. The pieces of that PseudoJet will
/// be the individual subjets
///
/// When constructed from a JetDefinition, that definition will be
/// used to obtain the subjets.  When constructed from a JetAlgorithm
/// and parameters (0 parameters for e+e-, just R or R and an extra
/// parameter for others) the recombination scheme will be taken as
/// the same one used to initially cluster the original jet.
///
/// The result of this transformer depends on its usage. There are two
/// typical use-cases: either we recluster one fat jet into subjets,
/// OR, we recluster the jet with a different jet alg. When Recluster
/// is created from a full jet definition. The last parameter of the
/// constructors below dicatate that behaviour: if "single" is true
/// (the default), a single jet, issued from a regular clustering is
/// returned (if there are more than one, the hardest is taken);
/// otherwise (single==false), the result will be a composite jet with
/// each subjet as pieces
/// 
/// Open points for discussion:
///
///  - do we add an option to force area support? [could be useful
///    e.g. for the filter with a subtractor where area support is
///    mandatory]
///
class Recluster : public Transformer {
public:
  /// define a recluster that decomposes a jet into subjets using a
  /// generic JetDefinition
  ///
  ///  \param subjet_def   the jet definition applied to obtain the subjets
  ///  \param single       when true, cluster the jet in a single jet (the
  ///                      hardest one) with an associated ClusterSequence, 
  ///                      otherwise return a composite jet with subjets
  ///                      as pieces.
  Recluster(const JetDefinition & subjet_def, bool single=true)
    : _subjet_def(subjet_def), _use_full_def(true), _single(single) {}

  /// define a recluster that decomposes a jet into subjets using a
  /// JetAlgorithm and its parameters
  ///
  ///  \param subjet_alg    the jet algorithm applied to obtain the subjets
  ///  \param subjet_radius the jet radius if required
  ///  \param subjet_extra  optional extra parameters for the jet algorithm (only when needed)
  ///  \param single        when true, cluster the jet in a single jet (the
  ///                       hardest one) with an associated ClusterSequence, 
  ///                       otherwise return a composite jet with subjets
  ///                       as pieces.
  /// 
  /// Typically, for e+e- algoriothm you should use the third version
  /// below with no parameters, for "standard" pp algorithms, just the
  /// clustering radius has to be specified and for genkt-type of
  /// algorithms, both the radius and the extra parameter have to be
  /// specified.
  Recluster(JetAlgorithm subjet_alg, double subjet_radius, double subjet_extra,
            bool single=true)
    : _subjet_alg(subjet_alg), _use_full_def(false), 
      _subjet_radius(subjet_radius), _has_subjet_radius(true), 
      _subjet_extra(subjet_extra), _has_subjet_extra(true), _single(single) {}
  Recluster(JetAlgorithm subjet_alg, double subjet_radius, bool single=true)
    : _subjet_alg(subjet_alg), _use_full_def(false), 
      _subjet_radius(subjet_radius), _has_subjet_radius(true),
      _has_subjet_extra(false), _single(single) {}
  Recluster(JetAlgorithm subjet_alg, bool single=true)
    : _subjet_alg(subjet_alg), _use_full_def(false), 
      _has_subjet_radius(false), _has_subjet_extra(false), _single(single) {}

  /// default dtor
  virtual ~Recluster(){};

  //----------------------------------------------------------------------
  // standard Transformer behaviour inherited from the base class
  // (i.e. result(), description() and structural info)

  /// runs the reclustering and sets kept and rejected to be the jets of interest
  /// (with non-zero rho, they will have been subtracted).
  ///
  /// \param jet    the jet that gets reclustered
  /// \return the reclustered jet
  virtual PseudoJet result(const PseudoJet & jet) const;

  /// class description
  virtual std::string description() const;

  // the type of the associated structure
  typedef CompositeJetStructure StructureType;

private:
  /// set the reclustered elements in the simple case of C/A+C/A
  void _recluster_cafilt(const std::vector<PseudoJet> & all_pieces,
                         std::vector<PseudoJet> & subjets,
                         double Rfilt) const;

  /// set the reclustered elements in the generic re-clustering case
  void _recluster_generic(const PseudoJet & jet, 
                          std::vector<PseudoJet> & subjets,
                          const JetDefinition & subjet_def,
                          bool do_areas) const;
  
  // a series of checks
  //--------------------------------------------------------------------
  /// get the pieces down to the fundamental pieces
  bool _get_all_pieces(const PseudoJet &jet, std::vector<PseudoJet> &all_pieces) const;

  /// get the common recombiner to all pieces (NULL if none)
  const JetDefinition::Recombiner* _get_common_recombiner(const std::vector<PseudoJet> &all_pieces) const;

  /// construct the proper jet definition ensuring that the recombiner
  /// is taken from the underlying pieces (an error is thrown if the
  /// pieces do no share a common recombiner)
  void _build_jet_def_with_recombiner(const std::vector<PseudoJet> &all_pieces, 
                                      JetDefinition &subjet_def) const;

  /// check if one can apply the simplified trick for C/A subjets
  bool _check_ca(const std::vector<PseudoJet> &all_pieces, 
                 const JetDefinition &subjet_def) const;

  /// check if the jet (or all its pieces) have explicit ghosts
  /// (assuming the jet has area support
  ///
  /// Note that if the jet has an associated cluster sequence that is no
  /// longer valid, an error will be thrown
  bool _check_explicit_ghosts(const std::vector<PseudoJet> &all_pieces) const;

  JetDefinition _subjet_def;   ///< the jet definition to use to extract the subjets
  JetAlgorithm  _subjet_alg;   ///< the jet algorithm to be used
  bool _use_full_def;          ///< true when the full JetDefinition is supplied to the ctor
  double _subjet_radius;       ///< the jet radius (only if needed for the jet alg)
  bool _has_subjet_radius;     ///< the subjet radius has been specified
  double _subjet_extra;        ///< the jet alg extra param (only if needed)
  bool _has_subjet_extra;      ///< the extra param has been specified

  bool _single;                ///< (true) return a single jet with a
                               ///< regular clustering or (false) a
                               ///< composite jet with subjets as pieces

  static LimitedWarning   _explicit_ghost_warning;
};

} // namespace contrib

FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh

#endif   // __FASTJET_TOOLS_RECLUSTER_HH__
