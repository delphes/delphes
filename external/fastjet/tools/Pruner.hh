#ifndef __FASTJET_TOOLS_PRUNER_HH__
#define __FASTJET_TOOLS_PRUNER_HH__

//FJSTARTHEADER
// $Id: Pruner.hh 3481 2014-07-29 17:24:12Z soyez $
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

#include "fastjet/ClusterSequence.hh"
#include "fastjet/WrappedStructure.hh"
#include "fastjet/tools/Transformer.hh"
#include <iostream>
#include <string>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

// fwd declarations
class Pruner;
class PrunerStructure;
class PruningRecombiner;
class PruningPlugin;

// This tells third-party code that the pruner structure 
// stores Rcut info; the alternative is for the user to 
// get the information from the version number
#define FASTJET_PRUNER_STRUCTURE_STORES_RCUT

//----------------------------------------------------------------------
/// @ingroup tools_generic
/// \class Pruner
/// Transformer that prunes a jet
///
/// This transformer prunes a jet according to the ideas presented in 
/// arXiv:0903.5081 (S.D. Ellis, C.K. Vermilion and J.R. Walsh). 
/// 
/// The jet's constituents are reclustered with a user-specified jet
/// definition, with the modification that objects i and j are only
/// recombined if at least one of the following two criteria is
/// satisfied:
///
///  - the geometric distance between i and j is smaller than 'Rcut'
///    with Rcut = Rcut_factor*2m/pt (Rcut_factor is a parameter of
///    the Pruner and m and pt obtained from the jet being pruned)
///  - the transverse momenta of i and j are at least 'zcut' p_t(i+j)
///
/// If both these criteria fail, i and j are not recombined, the 
/// harder of i and j is kept, and the softer is rejected. 
///
/// Usage: 
/// \code
///    Pruner pruner(jet_def, zcut, Rcut_factor);
///    PseudoJet pruned_jet = pruner(jet);
/// \endcode
///
/// The pruned_jet has a valid associated cluster sequence. In addition
/// the subjets of the original jet that have been vetoed by pruning
/// (i.e. have been 'pruned away') can be accessed using
///
/// \code
///   vector<PseudoJet> rejected_subjets = pruned_jet.structure_of<Pruner>().rejected();
/// \endcode
///
/// If the re-clustering happens to find more than a single inclusive
/// jet (this should normally not happen if the radius of the jet
/// definition used for the reclustering was set large enough),
/// the hardest of these jets is retured as the result of the
/// Pruner. The other jets can be accessed through
///
/// \code
///   vector<PseudoJet> extra_jets = pruned_jet.structure_of<Pruner>().extra_jets();
/// \endcode
///
/// Instead of using Rcut_factor and zcut, one can alternatively
/// construct a Pruner by passing two (pointers to) functions of 
/// PseudoJet that dynamically compute the Rcut and zcut to 
/// be used for the jet being pruned.
///
/// When the jet being pruned has area support and explicit ghosts,
/// the resulting pruned jet will likewise have area.
///
//----------------------------------------------------------------------
class Pruner : public Transformer{
public:
  /// minimal constructor, which takes a jet algorithm, sets the radius
  /// to JetDefinition::max_allowable_R (practically equivalent to
  /// infinity) and also tries to use a recombiner based on the one in
  /// the jet definition of the particular jet being pruned.
  ///
  ///  \param jet_alg     the jet algorithm for the internal clustering
  ///  \param zcut        pt-fraction cut in the pruning
  ///  \param Rcut_factor the angular distance cut in the pruning will be
  ///                     Rcut_factor * 2m/pt
  Pruner(const JetAlgorithm jet_alg, double zcut, double Rcut_factor) 
    : _jet_def(jet_alg, JetDefinition::max_allowable_R),
      _zcut(zcut), _Rcut_factor(Rcut_factor),
      _zcut_dyn(0), _Rcut_dyn(0), _get_recombiner_from_jet(true) {}


  /// alternative ctor in which the full reclustering jet definition can
  /// be specified.
  ///
  ///  \param jet_def     the jet definition for the internal clustering
  ///  \param zcut        pt-fraction cut in the pruning
  ///  \param Rcut_factor the angular distance cut in the pruning will be
  ///                     Rcut_factor * 2m/pt
  Pruner(const JetDefinition &jet_def, double zcut, double Rcut_factor)
    : _jet_def(jet_def),
      _zcut(zcut), _Rcut_factor(Rcut_factor),
      _zcut_dyn(0), _Rcut_dyn(0), _get_recombiner_from_jet(false) {}


  /// alternative ctor in which the pt-fraction cut and angular distance
  /// cut are functions of the jet being pruned.
  ///
  ///  \param jet_def the jet definition for the internal clustering
  ///  \param zcut_dyn    dynamic pt-fraction cut in the pruning
  ///  \param Rcut_dyn    dynamic angular distance cut in the pruning
  Pruner(const JetDefinition &jet_def, 
         const FunctionOfPseudoJet<double> *zcut_dyn,
         const FunctionOfPseudoJet<double> *Rcut_dyn);

  /// action on a single jet
  virtual PseudoJet result(const PseudoJet &jet) const;

  /// description
  virtual std::string description() const;
  
  // the type of the associated structure
  typedef PrunerStructure StructureType;

private:
  /// check if the jet has explicit_ghosts (knowing that there is an
  /// area support)
  bool _check_explicit_ghosts(const PseudoJet &jet) const;

  /// see if there is a common recombiner among the pieces; if there
  /// is return true and set jet_def_for_recombiner so that the
  /// recombiner can be taken from that JetDefinition. Otherwise,
  /// return false. 'assigned' is initially false; when true, each
  /// time we meet a new jet definition, we'll check it shares the
  /// same recombiner as jet_def_for_recombiner.
  bool _check_common_recombiner(const PseudoJet &jet, 
				JetDefinition &jet_def_for_recombiner,
				bool assigned=false) const;

  JetDefinition _jet_def; ///< the internal jet definition
  double _zcut;        	  ///< the pt-fraction cut
  double _Rcut_factor;    ///< the angular separation cut factor
  const FunctionOfPseudoJet<double> *_zcut_dyn; ///< dynamic zcut
  const FunctionOfPseudoJet<double> *_Rcut_dyn; ///< dynamic Rcut
  bool   _get_recombiner_from_jet; ///< true for minimal constructor,
                                   ///< causes recombiner to be set equal 
                                   ///< to that already used in the jet 
                                   ///< (if it can be deduced)
};


//----------------------------------------------------------------------
/// @ingroup tools_generic
/// \class PrunerStructure
/// The structure associated with a PseudoJet thas has gone through a
/// Pruner transformer
//----------------------------------------------------------------------
class PrunerStructure : public WrappedStructure{
public:
  /// default ctor
  ///  \param result_jet  the jet for which we have to keep the structure
  PrunerStructure(const PseudoJet & result_jet)
    : WrappedStructure(result_jet.structure_shared_ptr()){}

  /// description
  virtual std::string description() const{ return "Pruned PseudoJet";}

  /// return the constituents that have been rejected
  std::vector<PseudoJet> rejected() const{ 
    return validated_cs()->childless_pseudojets();
  }

  /// return the other jets that may have been found along with the
  /// result of the pruning
  /// The resulting vector is sorted in pt
  std::vector<PseudoJet> extra_jets() const;

  /// return the value of Rcut that was used for this specific pruning.
  double Rcut() const {return _Rcut;}

  /// return the value of Rcut that was used for this specific pruning.
  double zcut() const {return _zcut;}

protected:
  friend class Pruner; ///< to allow setting the internal information

private:
  double _Rcut, _zcut;
};

//----------------------------------------------------------------------
/// \if internal_doc
/// @ingroup internal
/// \class PruningRecombiner
/// recombines the objects that are not vetoed by pruning
///
/// This recombiner only recombines, using the provided 'recombiner',
/// objects (i and j) that pass at least one of the following two criteria:
///
///  - the geometric distance between i and j is smaller than 'Rcut'
///  - the transverse momenta of i and j are at least 'zcut' p_t(i+j)
///
/// If both these criteria fail, the hardest of i and j is kept and
/// the softest is rejected.
///
/// Note that this in not meant for standalone use [in particular
/// because it could lead to memory issues due to the rejected indices
/// stored internally].
///
/// \endif
class PruningRecombiner : public JetDefinition::Recombiner{
public:
  /// ctor
  ///  \param zcut   transverse momentum fraction cut
  ///  \param Rcut   angular separation cut
  ///  \param recomb pointer to a recombiner to use to cluster pairs
  PruningRecombiner(double zcut, double Rcut, 
		    const JetDefinition::Recombiner *recombiner)
    : _zcut2(zcut*zcut), _Rcut2(Rcut*Rcut), 
      _recombiner(recombiner){}

  /// perform a recombination taking into account the pruning
  /// conditions
  virtual void recombine(const PseudoJet &pa, 
			 const PseudoJet &pb,
			 PseudoJet &pab) const;

  /// returns the description of the recombiner
  virtual std::string description() const;

  /// return the history indices that have been pruned away
  const std::vector<unsigned int> & rejected() const{ return _rejected;}

  /// clears the list of rejected indices
  ///
  /// If one decides to use this recombiner standalone, one has to
  /// call this after each clustering in order for the rejected() vector
  /// to remain sensible and not grow to infinite size.
  void clear_rejected(){ _rejected.clear();}

private:
  double _zcut2;  ///< transverse momentum fraction cut (squared)
  double _Rcut2;  ///< angular separation cut (squared)
  const JetDefinition::Recombiner *_recombiner; ///< the underlying recombiner to use
  mutable std::vector<unsigned int> _rejected;  ///< list of rejected history indices
};


//----------------------------------------------------------------------
/// \if internal_doc
/// @ingroup internal
/// \class PruningPlugin
/// FastJet internal plugin that clusters the particles using the
/// PruningRecombiner. 
///
/// See PruningRecombiner for a description of what pruning does.
///
/// Note that this is an internal FastJet class used by the Pruner
/// transformer and it is not meant to be used as a standalone clustering
/// tool.
///
/// \endif
//----------------------------------------------------------------------
class PruningPlugin : public JetDefinition::Plugin{
public:
  /// ctor
  ///  \param jet_def the jet definition to be used for the 
  ///                 internal clustering
  ///  \param zcut    transverse momentum fraction cut
  ///  \param Rcut    angular separation cut
  PruningPlugin(const JetDefinition &jet_def, double zcut, double Rcut)
    : _jet_def(jet_def), _zcut(zcut), _Rcut(Rcut){}

  /// the actual clustering work for the plugin
  virtual void run_clustering(ClusterSequence &input_cs) const;

  /// description of the plugin
  virtual std::string description() const;

  /// returns the radius
  virtual double R() const {return _jet_def.R();}

private:
  JetDefinition _jet_def; ///< the internal jet definition
  double _zcut;           ///< transverse momentum fraction cut 
  double _Rcut;           ///< angular separation cut
};



FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh

#endif   // __FASTJET_TOOLS_PRUNER_HH__
