// $Id: RecursiveSymmetryCutBase.hh 1074 2017-09-18 15:15:20Z gsoyez $
//
// Copyright (c) 2014-, Gavin P. Salam, Gregory Soyez, Jesse Thaler
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

#ifndef __FASTJET_CONTRIB_RECURSIVESYMMETRYCUTBASE_HH__
#define __FASTJET_CONTRIB_RECURSIVESYMMETRYCUTBASE_HH__

#include <limits>
#include <cassert>
#include <fastjet/internal/base.hh>
#include "fastjet/tools/Transformer.hh"
#include "fastjet/WrappedStructure.hh"
#include "fastjet/CompositeJetStructure.hh"

#include "fastjet/config.h"

// we'll use the native FJ class for reculstering if available
#if FASTJET_VERSION_NUMBER >= 30100
#include "fastjet/tools/Recluster.hh"
#else
#include "Recluster.hh"
#endif

/** \mainpage RecursiveTools contrib 

    The RecursiveTools contrib provides a set of tools for
    recursive investigation jet substructure. Currently it includes:
    - fastjet::contrib::ModifiedMassDropTagger
    - fastjet::contrib::SoftDrop
    - fastjet::contrib::RecursiveSymmetryCutBase (from which the above two classes derive)
    - fastjet::contrib::IteratedSoftDropSymmetryFactors (defines ISD procedure)
    - fastjet::contrib::IteratedSoftDropMultiplicity (defines a useful observable using ISD)  
    - example*.cc provides usage examples
 */


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

namespace contrib{

//------------------------------------------------------------------------
/// \class RecursiveSymmetryCutBase
/// A base class for all the tools that de-cluster a jet until a
/// sufficiently symmetric configuration if found.
///
/// Derived classes (so far, ModifiedMassDropTagger and SoftDrop) have to
/// implement the symmetry cut and its description
///
/// Note that by default, the jet will be reculstered with
/// Cambridge/Aachen before applying the de-clustering procedure. This
/// behaviour can be changed using set_clustering (see below).
///
/// By default, this class behaves as a tagger, i.e. returns an empty
/// PseudoJet if no substructure is found.  While the derived
/// ModifiedMassDropTagger is a tagger, the derived SoftDrop is a groomer
/// (i.e. it returns a non-zero jet even if no substructure is found).
///
/// This class provides support for
///  - an optional mass-drop cut (see ctor)
///  - options to select which symmetry measure should be used (see ctor)
///  - options to select how the recursion proceeds (see ctor)
///  - options for reclustering the jet before running the de-clustering
///    (see set_reclustering)
///  - an optional subtractor (see ctor and other methods below)
class RecursiveSymmetryCutBase : public Transformer {
public:
  // ctors and dtors
  //----------------------------------------------------------------------

  /// an enum of the different (a)symmetry measures that can be used
  enum SymmetryMeasure{scalar_z,   ///< \f$ \min(p_{ti}, p_{tj})/(p_{ti} + p_{tj}) \f$
                       vector_z,   ///< \f$ \min(p_{ti}, p_{tj})/p_{t(i+j)} \f$
                       y,          ///< \f$ \min(p_{ti}^2,p_{tj}^2) \Delta R_{ij}^2 / m_{ij}^2 \f$
                       theta_E,    ///< \f$ \min(E_i,E_j)/(E_i+E_j) \f$ with 3d angle (ee collisions)
                       cos_theta_E ///< \f$ \min(E_i,E_j)/(E_i+E_j) \f$ with
                                   ///  \f$ \sqrt{2[1-cos(theta)]}\f$ for angles (ee collisions)
  };

  /// an enum for the options of how to choose which of two subjets to recurse into
  enum RecursionChoice{larger_pt, ///< choose the subjet with larger \f$ p_t \f$
                       larger_mt, ///< choose the subjet with larger \f$ m_t \equiv (m^2+p_t^2)^{\frac12}] \f$
                       larger_m,  ///< choose the subjet with larger mass (deprecated)
                       larger_E   ///< choose the subjet with larger energy (meant for ee collisions)
  }; 

  /// Full constructor, which takes the following parameters:
  ///
  /// \param symmetry_measure   the choice of measure to use to estimate the symmetry
  /// \param mu_cut             the maximal allowed value of mass drop variable mu = m_heavy/m_parent 
  /// \param recursion_choice   the strategy used to decide which subjet to recurse into
  /// \param subtractor         an optional pointer to a pileup subtractor (ignored if zero)
  ///
  /// If the (optional) pileup subtractor is supplied, then, by
  /// default, the input jet is assumed unsubtracted and the
  /// RecursiveSymmetryCutBase returns a subtracted 4-vector. [see
  /// also the set_input_jet_is_subtracted() member function].
  ///
  RecursiveSymmetryCutBase(SymmetryMeasure  symmetry_measure = scalar_z,
			   double           mu_cut = std::numeric_limits<double>::infinity(), 
			   RecursionChoice  recursion_choice = larger_pt,
			   const FunctionOfPseudoJet<PseudoJet> * subtractor = 0
			   ) : 
    _symmetry_measure(symmetry_measure),
    _mu_cut(mu_cut),
    _recursion_choice(recursion_choice),
    _subtractor(subtractor), _input_jet_is_subtracted(false),
    _do_reclustering(true), _recluster(0),
    _grooming_mode(false),
    _verbose_structure(false) // by default, don't story verbose information
  {}
    
  /// default destructor
  virtual ~RecursiveSymmetryCutBase(){}

  // access to class info
  //----------------------------------------------------------------------
  SymmetryMeasure symmetry_measure() const { return _symmetry_measure; }
  double mu_cut() const { return _mu_cut; }
  RecursionChoice recursion_choice() const { return _recursion_choice; }

  // internal subtraction configuration
  //----------------------------------------------------------------------

  /// This tells the tagger whether to assume that the input jet has
  /// already been subtracted. This is relevant only if a non-null
  /// subtractor pointer has been supplied, and the default assymption
  /// is that the input jet is passed unsubtracted.
  /// 
  /// Note: given that subtractors usually change the momentum of the
  /// main jet, but not that of the subjets, subjets will continue to
  /// have subtraction applied to them.
  void set_input_jet_is_subtracted(bool is_subtracted) {_input_jet_is_subtracted = is_subtracted;}

  /// returns a bool to indicate if the input jet is assumed already subtracted
  bool input_jet_is_subtracted() const {return _input_jet_is_subtracted;}

  /// an alternative way to set the subtractor
  ///
  /// Note that when a subtractor is provided, the result of the
  /// RecursiveSymmetryCutBase will be a subtracted jet.
  void set_subtractor(const FunctionOfPseudoJet<PseudoJet> * subtractor_) {_subtractor = subtractor_;}

  /// returns a pointer to the subtractor
  const FunctionOfPseudoJet<PseudoJet> * subtractor() const {return _subtractor;}

  // reclustering behaviour
  //----------------------------------------------------------------------

  /// configure the reclustering prior to the recursive de-clustering
  ///  \param do_reclustering   recluster the jet or not?
  ///  \param recluster         how to recluster the jet 
  ///                           (only if do_recluster is true; 
  ///                           Cambridge/Aachen used if NULL)
  ///
  /// Note that the ModifiedMassDropTagger and SoftDrop are designed
  /// to work with a Cambridge/Aachen clustering. Use any other option
  /// at your own risk!
  void set_reclustering(bool do_reclustering=true, const Recluster *recluster=0){
    _do_reclustering = do_reclustering;
    _recluster = recluster;
  }

  // what to do when no substructure is found
  //----------------------------------------------------------------------
  /// specify the behaviour adopted when no substructure is found
  ///  - in tagging  mode, an empty PseudoJet will be returned
  ///  - in grooming mode, a single particle is returned
  /// for clarity, we provide both function although they are redundant
  void set_grooming_mode(bool enable=true){ _grooming_mode = enable;}
  void set_tagging_mode(bool enable=true){ _grooming_mode = !enable;}

  
  /// Allows access to verbose information about recursive declustering,
  /// in particular values of symmetry, delta_R, and mu of dropped branches
  void set_verbose_structure(bool enable=true) { _verbose_structure = enable; }
  bool has_verbose_structure() const { return _verbose_structure; }
  
  
  // inherited from the Transformer base
  //----------------------------------------------------------------------

  /// the function that carries out the tagging; if a subtractor is
  /// being used, then this function assumes that input jet is
  /// unsubtracted (unless set_input_jet_is_subtracted(true) has been
  /// explicitly called before) and the result of the MMDT will be a
  /// subtracted jet.
  virtual PseudoJet result(const PseudoJet & j) const;
  
  /// description of the tool
  virtual std::string description() const;

  /// returns the gepometrical distance between the two particles
  /// depending on the symmetry measure used
  double squared_geometric_distance(const PseudoJet &j1,
                                    const PseudoJet &j2) const;
  

  class StructureType;

  /// for testing 
  static bool _verbose;

protected:
  // the methods below have to be defined by deerived classes
  //----------------------------------------------------------------------
  /// the cut on the symmetry measure (typically z) that one need to
  /// apply for a given pair of subjets p1 and p2
  virtual double symmetry_cut_fn(const PseudoJet & /* p1 */, 
                                 const PseudoJet & /* p2 */,
                                 void *extra_parameters = 0) const = 0;
  /// the associated dwescription
  virtual std::string symmetry_cut_description() const = 0;

  //----------------------------------------------------------------------
  /// this defines status codes when checking for substructure
  enum RecursionStatus{
    recursion_success=0,   //< found some substructure
    recursion_dropped,     //< dropped softest prong; recursion continues
    recursion_no_parents,  //< down to constituents; bottom of recursion
    recursion_issue        //< something went wrong; recursion stops
  };
  
  //----------------------------------------------------------------------
  /// the method below is the one actually performing one step of the
  /// recursion.
  ///
  /// It returns a status code (defined above)
  ///
  /// In case of success, all the information is filled
  /// In case of "no parents", piee1 is the same subjet
  /// In case of trouble, piece2 will be a 0 PJ and piece1 is the PJ we
  ///   should return (either 0 itself if the issue was critical, or
  ///   non-wero in case of a minor issue just causing the recursion to
  ///   stop)
  ///
  /// The extra_parameter argument allows one to pass extra agruments
  /// to the symmetry condition
  RecursionStatus recurse_one_step(const PseudoJet & subjet,
                                   PseudoJet &piece1, PseudoJet &piece2,
                                   double &sym, double &mu2,
                                   void *extra_parameters = 0) const;

  //----------------------------------------------------------------------
  /// helper for handling the reclustering
  PseudoJet _recluster_if_needed(const PseudoJet &jet) const;
  
  //----------------------------------------------------------------------
  // helpers for selecting between ee and pp cases
  bool is_ee() const{
    return ((_symmetry_measure==theta_E) || (_symmetry_measure==cos_theta_E));
  }

private:
  SymmetryMeasure  _symmetry_measure;
  double           _mu_cut;
  RecursionChoice  _recursion_choice;
  const FunctionOfPseudoJet<PseudoJet> * _subtractor;
  bool             _input_jet_is_subtracted;

  bool _do_reclustering;       ///< start with a reclustering
  const Recluster *_recluster; ///< how to recluster the jet

  bool _grooming_mode;  ///< grooming or tagging mode

  static LimitedWarning   _negative_mass_warning;
  static LimitedWarning   _mu2_gt1_warning;
  //static LimitedWarning   _nonca_warning;
  static LimitedWarning   _explicit_ghost_warning;

  // additional verbose structure information
  bool _verbose_structure;
  
  /// decide what to return when no substructure has been found
  PseudoJet _result_no_substructure(const PseudoJet &last_parent) const;
};

  
  
//----------------------------------------------------------------------
/// class to hold the structure of a jet tagged by RecursiveSymmetryCutBase.
class RecursiveSymmetryCutBase::StructureType : public WrappedStructure {
public:
  StructureType(const PseudoJet & j) :
    WrappedStructure(j.structure_shared_ptr()), _delta_R(-1.0), _symmetry(-1.0), _mu(-1.0),
    _is_composite(false), _has_verbose(false) // by default, do not store verbose structure
  {}

  /// construct a structure with
  ///  - basic info inherited from the reference jet "j"
  ///  - a given deltaR       for substructure
  ///  - a given symmetry     for substructure
  ///  - a given mu parameter for substructure
  /// If j is a "copmposite jet", it means that it has further
  /// substructure to potentially recurse  into
  StructureType(const PseudoJet & j, double delta_R_in, double symmetry_in, double mu_in=-1.0) :
    WrappedStructure(j.structure_shared_ptr()), _delta_R(delta_R_in), _symmetry(symmetry_in), _mu(mu_in),
    _is_composite(dynamic_cast<const CompositeJetStructure*>(j.structure_ptr()) != NULL),
    _has_verbose(false) // by default, do not store verbose structure
  {}
  
  // information about kept branch
  double delta_R()  const {return _delta_R;}
  double thetag()   const {return _delta_R;}  // alternative name
  double symmetry() const {return _symmetry;}
  double zg()       const {return _symmetry;} // alternative name
  double mu()       const {return _mu;}
  
  // additional verbose information about dropped branches
  bool has_verbose() const { return _has_verbose;}
  void set_verbose(bool value) { _has_verbose = value;}

  /// returns true if the current jet has some substructure (i.e. has
  /// been tagged by the resursion) or not
  ///
  /// Note that this should include deltaR==0 (e.g. a perfectly
  /// collinear branching with SoftDrop)
  bool has_substructure() const { return _delta_R>=0; }
  
  /// number of dropped branches
  int dropped_count(bool global=true) const;
  
  /// delta_R of dropped branches
  /// when "global" is set, recurse into possile further substructure
  std::vector<double> dropped_delta_R(bool global=true) const;
  void set_dropped_delta_R(const std::vector<double> &v) { _dropped_delta_R = v; }

  /// symmetry values of dropped branches
  std::vector<double> dropped_symmetry(bool global=true) const;
  void set_dropped_symmetry(const std::vector<double> &v) { _dropped_symmetry = v; }

  /// mass drop values of dropped branches
  std::vector<double> dropped_mu(bool global=true) const;
  void set_dropped_mu(const std::vector<double> &v) { _dropped_mu = v; }
  
  /// maximum symmetry value dropped
  double max_dropped_symmetry(bool global=true) const;

  /// (global) list of groomed away elements as zg and thetag
  /// sorted from largest to smallest anlge
  std::vector<std::pair<double,double> > sorted_zg_and_thetag() const;

private:
  double _delta_R, _symmetry, _mu;
  friend class RecursiveSymmetryCutBase;

  bool _is_composite;
  
  // additional verbose information
  bool _has_verbose;
  // information about dropped values
  std::vector<double> _dropped_delta_R;
  std::vector<double> _dropped_symmetry;
  std::vector<double> _dropped_mu;

  bool check_verbose(const std::string &what) const{
    if (!_has_verbose){
      throw Error("RecursiveSymmetryCutBase::StructureType: Verbose structure must be turned on to get "+what+".");
      return false;
    }
    return true;
  }
    
  
};

} // namespace contrib

FASTJET_END_NAMESPACE

#endif  // __FASTJET_CONTRIB_RECURSIVESYMMETRYCUTBASE_HH__
