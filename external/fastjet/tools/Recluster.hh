#ifndef __FASTJET_TOOLS_RECLUSTER_HH__
#define __FASTJET_TOOLS_RECLUSTER_HH__

// $Id: Recluster.hh 3714 2014-09-30 09:47:31Z soyez $
//
// Copyright (c) 2014, Matteo Cacciari, Gavin P. Salam and Gregory Soyez
//
//----------------------------------------------------------------------
// This file is part of FastJet
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
#include <fastjet/FunctionOfPseudoJet.hh>   // to derive Recluster from FOfPJ<PJ>
#include <iostream>
#include <string>

// TODO:
//
//  - maintain Voronoi areas? Requires CSAB:has_voronoi_area()    {->UNASSIGNED}
//  - make sure the description of the class is OK



FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//----------------------------------------------------------------------
/// @ingroup tools_generic
/// \class Recluster
/// Recluster a jet's constituents with a new jet definition.
///
/// When Recluster is constructed from a JetDefinition, it is that
/// definition that will be used to obtain the new jets. The user may
/// then decide if the recombiner should be the one from that jet
/// definition or if it should be acquired from the jet being
/// processed (the default).
///
/// Alternatively, Recluster can be constructed from a jet algorithm
/// and an optional radius. In that case the recombiner is
/// systematically obtained fromn the jet being processed (unless you
/// call set_acquire_recombiner(false)). If only the jet algorithm is
/// specified, a default radius of max_allowable_R will be assumed if
/// needed.
///
/// Recluster has two possible behaviours:
///
///  - if it is constructed with keep=keep_only_hardest the hardest
///    inclusive jet is returned as a "standard" jet with an
///    associated cluster sequence (unless there were no inclusive
///    jets, in which case a zero jet is returned, with no associated
///    cluster sequence)
///
///  - if it is constructed with keep=keep_all
///    all the inclusive jets are joined into a composite jet
///
/// [Note that since the structure of the resulting PseudoJet depends
/// on its usage, this class inherits from
/// FunctionOfPseudoJet<PseudoJet> (including a description) rather
/// than being a full-fledged Transformer]
///
class Recluster : public FunctionOfPseudoJet<PseudoJet> {
public:
  /// the various options for the output of Recluster
  enum Keep{
    /// keep only the hardest inclusive jet and return a "standard" jet with
    /// an associated ClusterSequence [this will be the default]
    keep_only_hardest,
    /// keep all the inclusive jets. result() will join them into a composite
    /// jet
    keep_all
  };

  /// default constructor (uses an undefined JetDefinition, and so cannot
  /// be used directly).
  Recluster() : _new_jet_def(), _acquire_recombiner(true),
                _keep(keep_only_hardest), _cambridge_optimisation_enabled(true){}

  /// Constructs a Recluster object that reclusters a jet into a new jet
  /// using a generic JetDefinition
  ///
  ///  \param new_jet_def   the jet definition applied to do the reclustering
  ///  \param acquire_recombiner
  ///                       when true, the reclustering will guess the
  ///                       recombiner from the input jet instead of
  ///                       the one in new_jet_def. An error is then
  ///                       thrown if no consistent recombiner is found
  ///  \param keep_in       Recluster::keep_only_hardest: the result is
  ///                       the hardest inclusive jet after reclustering,
  ///                       returned as a "standard" jet.
  ///                       Recluster::keep_all: the result is a
  ///                       composite jet with the inclusive jets as pieces.
  Recluster(const JetDefinition & new_jet_def, 
            bool acquire_recombiner_in = false, 
            Keep keep_in = keep_only_hardest)
    : _new_jet_def(new_jet_def), _acquire_recombiner(acquire_recombiner_in), 
      _keep(keep_in), _cambridge_optimisation_enabled(true) {}

  /// Constructs a Recluster object that reclusters a jet into a new jet
  /// using a JetAlgorithm and its parameters
  ///
  ///  \param new_jet_alg    the jet algorithm applied to obtain the new clustering
  ///  \param new_jet_radius the jet radius
  ///  \param keep_in        Recluster::keep_only_hardest: the result is
  ///                        the hardest inclusive jet after reclustering,
  ///                        returned as a "standard" jet.
  ///                        Recluster::keep_all: the result is a
  ///                        composite jet with the inclusive jets as pieces.
  /// 
  /// This ctor will always acquire the recombiner from the jet being
  /// reclustered (it will throw if none can be found).  If you wish
  /// to use Recluster with an algorithm that requires an extra
  /// parameter (like the genkt algorithm), please specify the jet
  /// definition fully using the constructor above.
  Recluster(JetAlgorithm snew_jet_alg, double new_jet_radius, Keep keep_in = keep_only_hardest);

  /// constructor with just a jet algorithm, but no jet radius. If the
  /// algorithm requires a jet radius, JetDefinition::max_allowable_R will be used. 
  ///
  Recluster(JetAlgorithm new_jet_alg, Keep keep_in = keep_only_hardest);

  /// default dtor
  virtual ~Recluster(){}

  //----------------------------------------------------------------------
  // tweaking the behaviour and corresponding enquiry functions

  /// set whether the reclustering should attempt to acquire a
  /// recombiner from the input jet
  void set_acquire_recombiner(bool acquire) {_acquire_recombiner = acquire;}

  /// returns true if this reclusterer is set to acquire the
  /// recombiner from the input jet
  bool acquire_recombiner() const{ return _acquire_recombiner;}


  /// sets whether to try to optimise reclustering with
  /// Cambridge/Aachen algorithms (by not reclustering if the
  /// requested C/A reclustering can be obtained by using subjets of
  /// an input C/A jet or one composed of multiple C/A pieces from the
  /// same clustering sequence). By default this is enabled, and
  /// _should_ always be correct; disable it to test this statement!
  void set_cambridge_optimisation(bool enabled){ _cambridge_optimisation_enabled = enabled;}
  /// sets whether to try to optimise reclustering with Cambridge/Aachen algorithms (US spelling!)
  void set_cambridge_optimization(bool enabled){ _cambridge_optimisation_enabled = enabled;}

  /// returns true if the reclusterer tries to optimise reclustering
  /// with Cambridge/Aachen algorithms
  bool cambridge_optimization(){return _cambridge_optimisation_enabled;}
  bool cambridge_optimisation(){return _cambridge_optimisation_enabled;}

  /// set the behaviour with regards to keeping all resulting jets or
  /// just the hardest.
  void set_keep(Keep keep_in) {_keep = keep_in;}

  /// returns the current "keep" mode i.e. whether only the hardest
  /// inclusive jet is returned or all of them (see the Keep enum above)
  Keep keep() const{ return _keep;}


  //----------------------------------------------------------------------
  // retrieving info about the behaviour

  /// class description
  virtual std::string description() const;


  //----------------------------------------------------------------------
  // core action of ths class

  /// runs the reclustering and sets kept and rejected to be the jets
  /// of interest (with non-zero rho, they will have been
  /// subtracted). Normally this will be accessed through the base
  /// class's operator().
  ///
  /// \param jet    the jet that gets reclustered
  /// \return the reclustered jet
  virtual PseudoJet result(const PseudoJet & jet) const;

  /// A lower-level method that does the actual work of reclustering
  /// the input jet. The resulting jets are stored in output_jets.
  /// The jet definition that has been used can be accessed from the
  /// output_jets' ClusterSequence.
  ///
  /// \param input_jet       the (input) jet that one wants to recluster
  /// \param output_jets     inclusive jets resulting from the new clustering
  ///
  /// Returns true if the C/A optimisation has been used (this means
  /// that generate_output_jet then has to watch out for non-explicit-ghost
  /// areas that might be leftover)
  bool get_new_jets_and_def(const PseudoJet & input_jet, 
                            std::vector<PseudoJet> & output_jets) const;

  /// given a set of inclusive jets and a jet definition used, create the
  /// resulting PseudoJet;
  /// 
  /// If ca_optimisation_used then special care will be taken in
  /// deciding whether the final jet can legitimately have an area.
  PseudoJet generate_output_jet(std::vector<PseudoJet> & incljets,
                                bool ca_optimisation_used) const;


private:
  /// set the reclustered elements in the simple case of C/A+C/A
  void _recluster_ca(const std::vector<PseudoJet> & all_pieces,
                     std::vector<PseudoJet> & incljets,
                     double Rfilt) const;

  /// set the reclustered elements in the generic re-clustering case
  void _recluster_generic(const PseudoJet & jet, 
                          std::vector<PseudoJet> & incljets,
                          const JetDefinition & new_jet_def,
                          bool do_areas) const;
  
  // a series of checks
  //--------------------------------------------------------------------
  /// get the pieces down to the fundamental pieces
  bool _get_all_pieces(const PseudoJet &jet, std::vector<PseudoJet> &all_pieces) const;

  /// associate the proper recombiner taken from the underlying pieces
  /// (an error is thrown if the pieces do no share a common
  /// recombiner)
  void _acquire_recombiner_from_pieces(const std::vector<PseudoJet> &all_pieces, 
                                       JetDefinition &new_jet_def) const;

  /// check if one can apply the simplified trick for C/A subjets
  bool _check_ca(const std::vector<PseudoJet> &all_pieces, 
                 const JetDefinition &new_jet_def) const;

  /// check if the jet (or all its pieces) have explicit ghosts
  /// (assuming the jet has area support
  ///
  /// Note that if the jet has an associated cluster sequence that is no
  /// longer valid, an error will be thrown
  bool _check_explicit_ghosts(const std::vector<PseudoJet> &all_pieces) const;

  JetDefinition _new_jet_def;  ///< the jet definition to use to extract the jets
  bool _acquire_recombiner;    ///< get the recombiner from the input
                               ///< jet rather than from _new_jet_def
  Keep _keep;                  ///< dictates which inclusive jets are kept and
                               ///< which are returned (see Keep above)

  bool _cambridge_optimisation_enabled; ///<enable the checks to
                                        ///< perform optimisation when
                                        ///< C/A reclustering is asked

  static LimitedWarning   _explicit_ghost_warning;
};

FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh

#endif   // __FASTJET_TOOLS_RECLUSTER_HH__
