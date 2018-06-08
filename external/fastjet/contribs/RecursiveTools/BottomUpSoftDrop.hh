// $Id: BottomUpSoftDrop.hh 1085 2017-10-11 02:16:59Z jthaler $
//
// Copyright (c) 2017-, Gavin P. Salam, Gregory Soyez, Jesse Thaler,
// Kevin Zhou, Frederic Dreyer
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

#ifndef __BOTTOMUPSOFTDROP_HH__
#define __BOTTOMUPSOFTDROP_HH__

#include "fastjet/ClusterSequence.hh"
#include "fastjet/WrappedStructure.hh"
#include "fastjet/tools/Transformer.hh"

#include <iostream>
#include <string>

// TODO
//
//  - missing class description
//
//  - check what to do when pta=ptb=0
//    for the moment, we recombine both for multiple reasons
//     . this avois breakingteh symemtry between pa and pb
//     . it would be groomed in later steps anyway
//    Note that this is slightly inconsistent with our use of
//    > (instead of >=) in the cdt 

FASTJET_BEGIN_NAMESPACE

namespace contrib{
  
// fwd declarations
class BottomUpSoftDrop;
class BottomUpSoftDropStructure;
class BottomUpSoftDropRecombiner;
class BottomUpSoftDropPlugin;

//----------------------------------------------------------------------
/// \class BottomUpSoftDrop
/// Implementation of the BottomUpSoftDrop transformer
///
/// Bottom-Up Soft drop grooms a jet by applying a modified
/// recombination scheme, where particles are recombined only if they
/// pass the Soft Drop condition
///
/// \f[
///   z < z_{\rm cut} (\theta/R0)^\beta
/// \f]
///
/// the groomed jet contains the particles remaining after this
/// pair-wise recombination
///
/// Note:
///  - one can use BottomUpSoftDrop on a full event with the
///    global_grooming(event) method.
///  - if two recombined particles a and b have momentum pta=ptb=0,
///    we recombine both.
/// 

  
class BottomUpSoftDrop : public Transformer {
public:
  /// minimal constructor, which the jet algorithm to CA, sets the radius
  /// to JetDefinition::max_allowable_R (practically equivalent to
  /// infinity) and also tries to use a recombiner based on the one in
  /// the jet definition of the particular jet being Soft Dropped.
  ///
  ///  \param beta           the value for beta
  ///  \param symmetry_cut   the value for symmetry_cut
  ///  \param R0             the value for R0
  BottomUpSoftDrop(double beta, double symmetry_cut, double R0 = 1.0) 
    : _jet_def(cambridge_algorithm, JetDefinition::max_allowable_R),
      _beta(beta),_symmetry_cut(symmetry_cut), _R0(R0),
      _get_recombiner_from_jet(true) {}

  /// alternative constructor which takes a specified jet algorithm
  ///
  ///  \param jet_alg        the jet algorithm for the internal clustering (uses R=infty)
  ///  \param symmetry_cut   the value of symmetry_cut
  ///  \param beta           the value for beta
  ///  \param R0             the value for R0
  BottomUpSoftDrop(const JetAlgorithm jet_alg, double beta, double symmetry_cut,
		   double R0 = 1.0) 
    : _jet_def(jet_alg, JetDefinition::max_allowable_R),
      _beta(beta), _symmetry_cut(symmetry_cut), _R0(R0),
      _get_recombiner_from_jet(true) {}


  /// alternative ctor in which the full reclustering jet definition can
  /// be specified.
  ///
  ///  \param jet_def        the jet definition for the internal clustering
  ///  \param symmetry_cut   the value of symmetry_cut
  ///  \param beta           the value for beta
  ///  \param R0             the value for R0
  BottomUpSoftDrop(const JetDefinition &jet_def, double beta, double symmetry_cut,
		   double R0 = 1.0)
    : _jet_def(jet_def), _beta(beta), _symmetry_cut(symmetry_cut), _R0(R0),
      _get_recombiner_from_jet(false) {}

  /// action on a single jet
  virtual PseudoJet result(const PseudoJet &jet) const;

  /// global grooming on a full event
  /// note: does not support jet areas
  virtual std::vector<PseudoJet> global_grooming(const std::vector<PseudoJet> & event) const;
  
  /// description
  virtual std::string description() const;
  
  // the type of the associated structure
  typedef BottomUpSoftDropStructure StructureType;

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
  double _beta;           ///< the value of beta
  double _symmetry_cut;   ///< the value of symmetry_cut
  double _R0;             ///< the value of R0
  bool   _get_recombiner_from_jet; ///< true for minimal constructor,
                                   ///< causes recombiner to be set equal 
                                   ///< to that already used in the jet 
                                   ///< (if it can be deduced)
};

//----------------------------------------------------------------------
/// The structure associated with a PseudoJet thas has gone through a
/// bottom/up SoftDrop transformer
class BottomUpSoftDropStructure : public WrappedStructure{
public:
  /// default ctor
  ///  \param result_jet  the jet for which we have to keep the structure
  BottomUpSoftDropStructure(const PseudoJet & result_jet)
    : WrappedStructure(result_jet.structure_shared_ptr()){}

  /// description
  virtual std::string description() const{
    return "Bottom/Up Soft Dropped PseudoJet";
  }

  /// return the constituents that have been rejected
  std::vector<PseudoJet> rejected() const{ 
    return validated_cs()->childless_pseudojets();
  }

  /// return the other jets that may have been found along with the
  /// result of the bottom/up Soft Drop
  /// The resulting vector is sorted in pt
  std::vector<PseudoJet> extra_jets() const {
    return sorted_by_pt((!SelectorNHardest(1))(validated_cs()->inclusive_jets()));
  }

  /// return the value of beta that was used for this specific Soft Drop.
  double beta() const {return _beta;}

  /// return the value of symmetry_cut that was used for this specific Soft Drop.
  double symmetry_cut() const {return _symmetry_cut;}

  /// return the value of R0 that was used for this specific Soft Drop.
  double R0() const {return _R0;}

protected:
  friend class BottomUpSoftDrop; ///< to allow setting the internal information

private:
  double _beta, _symmetry_cut, _R0;
};

//----------------------------------------------------------------------
/// Class for Soft Drop recombination
/// recombines the objects that are not vetoed by Bottom-Up SoftDrop
///
/// This recombiner only recombines, using the provided 'recombiner',
/// objects (i and j) that pass the following SoftDrop criterion:
///
///   min(pti, ptj) > zcut (pti+ptj) (theta_ij/R0)^beta
///
/// If the criterion fail, the hardest of i and j is kept and the
/// softest is rejected.
///
/// Note that this in not meant for standalone use [in particular
/// because it could lead to memory issues due to the rejected indices
/// stored internally].
///
/// This class is a direct adaptation of PruningRecombiner in Fastjet tools
class BottomUpSoftDropRecombiner : public JetDefinition::Recombiner {
public:
  /// ctor
  ///  \param symmetry_cut   value of cut on symmetry measure
  ///  \param beta   avalue of beta parameter
  ///  \param recomb pointer to a recombiner to use to cluster pairs
  BottomUpSoftDropRecombiner(double beta, double symmetry_cut, double R0,
                             const JetDefinition::Recombiner *recombiner)
    : _beta(beta), _symmetry_cut(symmetry_cut), _R0sqr(R0*R0),
      _recombiner(recombiner) {}

  /// perform a recombination taking into account the Soft Drop
  /// conditions
  virtual void recombine(const PseudoJet &pa, 
			 const PseudoJet &pb,
			 PseudoJet &pab) const;

  /// returns the description of the recombiner
  virtual std::string description() const {
    std::ostringstream oss;
    oss << "SoftDrop recombiner with symmetry_cut = " << _symmetry_cut
	<< ", beta = " << _beta
	<< ", and underlying recombiner = " << _recombiner->description();
    return oss.str();    
  }

  /// return the history indices that have been soft dropped away
  const std::vector<unsigned int> & rejected() const{ return _rejected;}

  /// clears the list of rejected indices
  ///
  /// If one decides to use this recombiner standalone, one has to
  /// call this after each clustering in order for the rejected() vector
  /// to remain sensible and not grow to infinite size.
  void clear_rejected(){ _rejected.clear();}

private:
  double _beta;          ///< beta parameter
  double _symmetry_cut;  ///< value of symmetry_cut
  double _R0sqr;         ///< normalisation of the angular distance
  const JetDefinition::Recombiner *_recombiner; ///< the underlying recombiner to use
  mutable std::vector<unsigned int> _rejected;  ///< list of rejected history indices
};

//----------------------------------------------------------------------
/// \class BottomUpSoftDropPlugin
/// Class for a bottom/up Soft Drop algorithm, based on the Pruner plugin
///
/// This is an internal plugin that clusters the particles using the
/// BottomUpRecombiner. 
///
/// See BottomUpRecombiner for a description of what bottom-up
/// SoftDrop does.
///
/// Note that this is an internal class used by the BottomUpSoftDrop
/// transformer and it is not meant to be used as a standalone
/// clustering tool.
class BottomUpSoftDropPlugin : public JetDefinition::Plugin {
public:
  /// ctor
  ///  \param jet_def the jet definition to be used for the 
  ///                 internal clustering
  ///  \param symmetry_cut    value of cut on symmetry measure
  ///  \param beta    value of beta parameter
  BottomUpSoftDropPlugin(const JetDefinition &jet_def, double beta, double symmetry_cut,
                         double R0 = 1.0)
    : _jet_def(jet_def), _beta(beta), _symmetry_cut(symmetry_cut), _R0(R0) {}

  /// the actual clustering work for the plugin
  virtual void run_clustering(ClusterSequence &input_cs) const;

  /// description of the plugin
  virtual std::string description() const;

  /// returns the radius
  virtual double R() const {return _jet_def.R();}

private:
  JetDefinition _jet_def; ///< the internal jet definition
  double _beta;           ///< beta parameter
  double _symmetry_cut;   ///< value of symmetry_cut 
  double _R0;             ///< normalisation of the angular distance
};

}

FASTJET_END_NAMESPACE      // defined in fastjet/internal/base.hh
#endif // __BOTTOMUPSOFTDROP_HH__
