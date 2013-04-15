//STARTHEADER
// $Id: CASubJetTagger.hh 2616 2011-09-30 18:03:40Z salam $
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

#ifndef __CASUBJET_TAGGER_HH__
#define __CASUBJET_TAGGER_HH__

#include <fastjet/PseudoJet.hh>
#include <fastjet/WrappedStructure.hh>
#include <fastjet/tools/Transformer.hh>
#include "fastjet/LimitedWarning.hh"

FASTJET_BEGIN_NAMESPACE

class CASubJetTagger;
class CASubJetTaggerStructure;

//----------------------------------------------------------------------
/// @ingroup tools_taggers
/// \class CASubJetTagger
/// clean (almost parameter-free) tagger searching for the element in
/// the clustering history that maximises a chosen distance
///
/// class to help us get a clean (almost parameter-free) handle on
/// substructure inside a C/A jet. It follows the logic described in
/// arXiv:0906.0728 (and is inspired by the original Cambridge
/// algorithm paper in its use of separate angular and dimensionful
/// distances), but provides some extra flexibility.
///
/// It searches for all splittings that pass a symmetry cut (zcut) and
/// then selects the one with the largest auxiliary scale choice
/// (e.g. jade distance of the splitting, kt distance of the
/// splitting, etc.)
///
/// By default, the zcut is calculated from the fraction of the child
/// pt carried by the parent jet. If one calls set_absolute_z_cut the
/// fraction of transverse momentum will be computed wrt the original
/// jet.
///
/// original code copyright (C) 2009 by Gavin Salam, released under
/// the GPL.
///
/// \section desc Options
/// 
///  - the distance choice: options are
///        kt2_distance        : usual min(kti^2,ktj^2)DeltaR_{ij}^2
///        jade_distance       : kti . ktj DeltaR_{ij}^2 (LI version of jade)
///        jade2_distance      : kti . ktj DeltaR_{ij}^4 (LI version of jade * DR^2)
///        plain_distance      :           DeltaR_{ij}^2
///        mass_drop_distance  : m_jet - max(m_parent1,m_parent2)
///        dot_product_distance: parent1.parent2
///    (kt2_distance by default)
///
///  - the z cut (0 by default)
///
///  - by calling set_absolute_z_cut(), one can ask that the pt
///    fraction if calculated wrt the original jet
///
///  - by calling set_dr_min(drmin), one can ask that only the
///    recombinations where the 2 objects are (geometrically) distant
///    by at least drmin are kept in the maximisation.
///
/// \section input Input conditions
/// 
///  - the jet must have been obtained from a Cambridge/Aachen cluster
///    sequence
///
/// \section output Output/structure
/// 
///  - the element of the cluster sequence maximising the requested
///    distance (and satisfying the zcut) is returned.
///
///  - if the original jet has no parents, it will be returned
///
///  - the value of the "z" and distance corresponding to that history
///    element are stored and accessible through
///      result.structure_of<CASubJetTagger>().z();
///      result.structure_of<CASubJetTagger>().distance();
///
class CASubJetTagger : public Transformer {
public:

  /// the available choices of auxiliary scale with respect to which
  /// to order the splittings
  enum ScaleChoice {
    kt2_distance,        //< usual min(kti^2,ktj^2)DeltaR_{ij}^2
    jade_distance,       //< kti . ktj DeltaR_{ij}^2 (LI version of jade)
    jade2_distance,      //< kti . ktj DeltaR_{ij}^4 (LI version of jade * DR^2)
    plain_distance,      //<           DeltaR_{ij}^2
    mass_drop_distance,  //<    m_jet - max(m_parent1,m_parent2)
    dot_product_distance //<  parent1.parent2
  };

  /// just constructs
  CASubJetTagger(ScaleChoice scale_choice = jade_distance,
                 double      z_threshold  = 0.1)
    : _scale_choice(scale_choice), _z_threshold(z_threshold),
      _dr2_min(0.0), _absolute_z_cut(false){};

  /// sets a minimum delta R below which spliting will be ignored
  /// (only relevant if set prior to calling run())
  void set_dr_min(double drmin) {_dr2_min = drmin*drmin;}

  /// returns a textual description of the tagger
  virtual std::string description() const;

  /// If (abs_z_cut) is set to false (the default) then for a
  /// splitting to be considered, each subjet must satisfy 
  /// 
  ///        p_{t,sub} > z_threshold * p_{t,parent}
  ///
  /// whereas if it is set to true, then each subject must satisfy
  ///
  ///        p_{t,sub} > z_threshold * p_{t,original-jet}
  ///
  /// where parent is the immediate parent of the splitting, and
  /// original jet is the one supplied to the run() function.
  ///
  /// Only relevant is called prior to run().
  void set_absolute_z_cut(bool abs_z_cut=true) {_absolute_z_cut = abs_z_cut;}
  
  /// runs the tagger on the given jet and
  /// returns the tagged PseudoJet if successful, or a PseudoJet==0 otherwise
  /// (standard access is through operator()).
  virtual PseudoJet result(const PseudoJet & jet) const;

  /// the type of Structure returned
  typedef CASubJetTaggerStructure StructureType;

protected:
  /// class that contains the result internally
  class JetAux {
  public:
    PseudoJet jet;          //< the subjet (immediate parent of splitting)
    double    aux_distance; //< the auxiliary distance between its two subjets
    double    delta_r;      //< the angular distance between its two subjets
    double    z;            //< the transverse momentum fraction
  };


  void _recurse_through_jet(const PseudoJet & current_jet, 
                            JetAux &aux_max,
                            const PseudoJet & original_jet) const;

  ScaleChoice _scale_choice;
  double      _z_threshold;
  double      _dr2_min;
  bool        _absolute_z_cut;

  static  LimitedWarning _non_ca_warnings;
};


//------------------------------------------------------------------------
/// @ingroup tools_taggers
/// the structure returned by a CASubJetTagger
///
/// Since this is directly an element of the ClusterSequence, we keep
/// basically the original ClusterSequenceStructure (wrapped for
/// memory-management reasons) and add information about the pt
/// fraction and distance of the subjet structure
class CASubJetTaggerStructure : public WrappedStructure{
public:
  /// default ctor
  ///  \param result_jet   the jet for which we have to keep the
  ///                      structure
  CASubJetTaggerStructure(const PseudoJet & result_jet)
    : WrappedStructure(result_jet.structure_shared_ptr()){}

  /// returns the scale choice asked for the maximisation
  CASubJetTagger::ScaleChoice scale_choice() const {return _scale_choice;}

  /// returns the value of the distance measure (corresponding to
  /// ScaleChoice) for this jet's splitting
  double distance() const {return _distance;}

  /// returns the pt fraction contained by the softer of the two component
  /// pieces of this jet (normalised relative to this jet)
  double z() const {return _z;}

  /// returns the pt fraction contained by the softer of the two component
  /// pieces of this jet (normalised relative to the original jet)
  bool absolute_z() const {return _absolute_z;}

//  /// returns the original jet (before tagging)
//  const PseudoJet & original() const {return _original_jet;}

protected:
  CASubJetTagger::ScaleChoice _scale_choice; ///< the user scale choice 
  double _distance;  ///< the maximal distance associated with the result
  bool _absolute_z;  ///< whether z is computed wrt to the original jet or not
  double _z;         ///< the transverse momentum fraction
//  PseudoJet _original_jet;  ///< the original jet (before tagging)

  friend class CASubJetTagger; ///< to allow setting the internal information
};

FASTJET_END_NAMESPACE

#endif // __CASUBJET_HH__
