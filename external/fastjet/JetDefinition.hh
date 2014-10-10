#ifndef __FASTJET_JETDEFINITION_HH__
#define __FASTJET_JETDEFINITION_HH__

//FJSTARTHEADER
// $Id: JetDefinition.hh 3677 2014-09-09 22:45:25Z soyez $
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

#include<cassert>
#include "fastjet/internal/numconsts.hh"
#include "fastjet/PseudoJet.hh"
#include<string>
#include<memory>

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

/// return a string containing information about the release
//  NB: (implemented in ClusterSequence.cc but defined here because
//  this is a visible location)
std::string fastjet_version_string();

//======================================================================
/// the various options for the algorithmic strategy to adopt in
/// clustering events with kt and cambridge style algorithms.
enum Strategy {
  /// Like N2MHTLazy9 in a number of respects, but does not calculate
  /// ghost-ghost distances and so does not carry out ghost-ghost
  /// recombination. 
  ///
  /// If you want active ghosted areas, then this is only suitable for
  /// use with the anti-kt algorithm (or genkt with negative p), and
  /// does not produce any pure ghost jets. If used with active areas
  /// with Kt or Cam algorithms it will actually produce a passive
  /// area.
  /// 
  /// Particles are deemed to be ghosts if their pt is below a
  /// threshold (currently 1e-50, hard coded as ghost_limit in
  /// LazyTiling9SeparateGhosts).
  ///
  /// Currently for events with a couple of thousand normal particles
  /// and O(10k) ghosts, this can be quicker than N2MHTLazy9, which
  /// would otherwise be the best strategy. 
  ///
  /// New in FJ3.1
  N2MHTLazy9AntiKtSeparateGhosts   = -10, 
  /// only looks into a neighbouring tile for a particle's nearest
  /// neighbour (NN) if that particle's in-tile NN is further than the
  /// distance to the edge of the neighbouring tile. Uses tiles of
  /// size R and a 3x3 tile grid around the particle.
  /// New in FJ3.1
  N2MHTLazy9   = -7, 
  /// Similar to N2MHTLazy9, but uses tiles of size R/2 and a 5x5 tile
  /// grid around the particle.
  /// New in FJ3.1
  N2MHTLazy25   = -6, 
  /// Like to N2MHTLazy9 but uses slightly different optimizations,
  /// e.g. for calculations of distance to nearest tile; as of
  /// 2014-07-18 it is slightly slower and not recommended for
  /// production use. To considered deprecated.
  /// New in FJ3.1
  N2MHTLazy9Alt   = -5, 
  /// faster that N2Tiled above about 500 particles; differs from it
  /// by retainig the di(closest j) distances in a MinHeap (sort of
  /// priority queue) rather than a simple vector. 
  N2MinHeapTiled   = -4, 
  /// fastest from about 50..500
  N2Tiled     = -3, 
  /// legacy
  N2PoorTiled = -2, 
  /// fastest below 50
  N2Plain     = -1, 
  /// worse even than the usual N^3 algorithms
  N3Dumb      =  0, 
  /// automatic selection of the best (based on N), including 
  /// the LazyTiled strategies that are new to FJ3.1
  Best        =  1, 
  /// best of the NlnN variants -- best overall for N>10^4.
  /// (Does not work for R>=2pi)
  NlnN        =  2, 
  /// legacy N ln N using 3pi coverage of cylinder.
  /// (Does not work for R>=2pi)
  NlnN3pi     =  3, 
  /// legacy N ln N using 4pi coverage of cylinder
  NlnN4pi     =  4,
  /// Chan's closest pair method (in a variant with 4pi coverage),
  /// for use exclusively with the Cambridge algorithm.
  /// (Does not work for R>=2pi)
  NlnNCam4pi   = 14,
  /// Chan's closest pair method (in a variant with 2pi+2R coverage),
  /// for use exclusively with the Cambridge algorithm.
  /// (Does not work for R>=2pi)
  NlnNCam2pi2R = 13,
  /// Chan's closest pair method (in a variant with 2pi+minimal extra
  /// variant), for use exclusively with the Cambridge algorithm. 
  /// (Does not work for R>=2pi)
  NlnNCam      = 12, // 2piMultD
  /// the automatic strategy choice that was being made in FJ 3.0
  /// (restricted to strategies that were present in FJ 3.0)
  BestFJ30     =  21, 
  /// the plugin has been used...
  plugin_strategy = 999
};


//======================================================================
/// \enum JetAlgorithm
/// the various families of jet-clustering algorithm
//
// [Remember to update the "is_spherical()" routine if any further
// spherical algorithms are added to the list below]
enum JetAlgorithm {
  /// the longitudinally invariant kt algorithm
  kt_algorithm=0,
  /// the longitudinally invariant variant of the cambridge algorithm
  /// (aka Aachen algoithm).
  cambridge_algorithm=1,
  /// like the k_t but with distance measures 
  ///       dij = min(1/kti^2,1/ktj^2) Delta R_{ij}^2 / R^2
  ///       diB = 1/kti^2
  antikt_algorithm=2, 
  /// like the k_t but with distance measures 
  ///       dij = min(kti^{2p},ktj^{2p}) Delta R_{ij}^2 / R^2
  ///       diB = 1/kti^{2p}
  /// where p = extra_param()
  genkt_algorithm=3, 
  /// a version of cambridge with a special distance measure for
  /// particles whose pt is < extra_param(); this is not usually
  /// intended for end users, but is instead automatically selected
  /// when requesting a passive Cambridge area.
  cambridge_for_passive_algorithm=11,
  /// a version of genkt with a special distance measure for particles
  /// whose pt is < extra_param() [relevant for passive areas when p<=0]
  /// ***** NB: THERE IS CURRENTLY NO IMPLEMENTATION FOR THIS ALG *******
  genkt_for_passive_algorithm=13, 
  //.................................................................
  /// the e+e- kt algorithm
  ee_kt_algorithm=50,
  /// the e+e- genkt algorithm  (R > 2 and p=1 gives ee_kt)
  ee_genkt_algorithm=53,
  //.................................................................
  /// any plugin algorithm supplied by the user
  plugin_algorithm = 99,
  //.................................................................
  /// the value for the jet algorithm in a JetDefinition for which
  /// no algorithm has yet been defined
  undefined_jet_algorithm = 999
};

/// make standard Les Houches nomenclature JetAlgorithm (algorithm is general
/// recipe without the parameters) backward-compatible with old JetFinder
typedef JetAlgorithm JetFinder;

/// provide other possible names for the Cambridge/Aachen algorithm
const JetAlgorithm aachen_algorithm = cambridge_algorithm;
const JetAlgorithm cambridge_aachen_algorithm = cambridge_algorithm;

//======================================================================
/// The various recombination schemes
///
/// Note that the schemes that recombine with non-linear weighting of
/// the directions (e.g. pt2, winner-takes-all) are collinear safe
/// only for algorithms with a suitable ordering of the
/// recombinations: orderings in which, for particles of comparable
/// energies, small-angle clusterings take place before large-angle
/// clusterings. This property is satisfied by all gen-kt algorithms.
/// 
enum RecombinationScheme {
  /// summing the 4-momenta
  E_scheme=0,
  /// pt weighted recombination of y,phi (and summing of pt's)
  /// with preprocessing to make things massless by rescaling E=|\vec p|
  pt_scheme=1,
  /// pt^2 weighted recombination of y,phi (and summing of pt's)
  /// with preprocessing to make things massless by rescaling E=|\vec p|
  pt2_scheme=2,
  /// pt weighted recombination of y,phi (and summing of pt's)
  /// with preprocessing to make things massless by rescaling |\vec p|->=E
  Et_scheme=3,
  /// pt^2 weighted recombination of y,phi (and summing of pt's)
  /// with preprocessing to make things massless by rescaling |\vec p|->=E
  Et2_scheme=4,
  /// pt weighted recombination of y,phi (and summing of pt's), with 
  /// no preprocessing
  BIpt_scheme=5,
  /// pt^2 weighted recombination of y,phi (and summing of pt's)
  /// no preprocessing
  BIpt2_scheme=6,
  /// pt-based Winner-Takes-All (WTA) recombination: the
  /// result of the recombination has the rapidity, azimuth and mass
  /// of the the PseudoJet with the larger pt, and a pt equal to the
  /// sum of the two pt's
  WTA_pt_scheme=7,
  /// mod-p-based Winner-Takes-All (WTA) recombination: the result of
  /// the recombination gets the 3-vector direction and mass of the
  /// PseudoJet with the larger |3-momentum| (modp), and a
  /// |3-momentum| equal to the scalar sum of the two |3-momenta|.
  WTA_modp_scheme=8,
  // Energy-ordering can lead to dangerous situations with particles at
  // rest. We instead implement the WTA_modp_scheme
  //
  // // energy-based Winner-Takes-All (WTA) recombination: the result of
  // // the recombination gets the 3-vector direction and mass of the
  // // PseudoJet with the larger energy, and an energy equal to the
  // // to the sum of the two energies
  // WTA_E_scheme=8,
  /// for the user's external scheme
  external_scheme = 99
};



// forward declaration, needed in order to specify interface for the
// plugin.
class ClusterSequence;




//======================================================================
/// @ingroup basic_classes
/// \class JetDefinition
/// class that is intended to hold a full definition of the jet
/// clusterer
class JetDefinition {
  
public:

  /// forward declaration of a class that allows the user to introduce
  /// their own plugin 
  class Plugin;

  // forward declaration of a class that will provide the
  // recombination scheme facilities and/or allow a user to
  // extend these facilities
  class Recombiner;


  /// constructor with alternative ordering or arguments -- note that
  /// we have not provided a default jet finder, to avoid ambiguous
  /// JetDefinition() constructor.
  JetDefinition(JetAlgorithm jet_algorithm_in, 
                double R_in, 
                RecombinationScheme recomb_scheme_in = E_scheme,
                Strategy strategy_in = Best) {
    *this = JetDefinition(jet_algorithm_in, R_in, strategy_in, recomb_scheme_in, 1);
  }

  /// constructor for algorithms that have no free parameters
  /// (e.g. ee_kt_algorithm)
  JetDefinition(JetAlgorithm jet_algorithm_in, 
                RecombinationScheme recomb_scheme_in = E_scheme,
                Strategy strategy_in = Best) {
    double dummyR = 0.0;
    *this = JetDefinition(jet_algorithm_in, dummyR, strategy_in, recomb_scheme_in, 0);
  }

  /// constructor for algorithms that require R + one extra parameter to be set 
  /// (the gen-kt series for example)
  JetDefinition(JetAlgorithm jet_algorithm_in, 
                double R_in, 
                double xtra_param_in,
                RecombinationScheme recomb_scheme_in = E_scheme,
                Strategy strategy_in = Best) {
    *this = JetDefinition(jet_algorithm_in, R_in, strategy_in, recomb_scheme_in, 2);
    set_extra_param(xtra_param_in);
  }


  /// constructor in a form that allows the user to provide a pointer
  /// to an external recombiner class (which must remain valid for the
  /// life of the JetDefinition object).
  JetDefinition(JetAlgorithm jet_algorithm_in, 
                double R_in, 
                const Recombiner * recombiner_in,
                Strategy strategy_in = Best) {
    *this = JetDefinition(jet_algorithm_in, R_in, external_scheme, strategy_in);
    _recombiner = recombiner_in;
  }


  /// constructor for case with 0 parameters (ee_kt_algorithm) and
  /// and external recombiner
  JetDefinition(JetAlgorithm jet_algorithm_in, 
                const Recombiner * recombiner_in,
                Strategy strategy_in = Best) {
    *this = JetDefinition(jet_algorithm_in, external_scheme, strategy_in);
    _recombiner = recombiner_in;
  }

  /// constructor allowing the extra parameter to be set and a pointer to
  /// a recombiner
  JetDefinition(JetAlgorithm jet_algorithm_in, 
                double R_in, 
                double xtra_param_in,
                const Recombiner * recombiner_in,
                Strategy strategy_in = Best) {
    *this = JetDefinition(jet_algorithm_in, R_in, xtra_param_in, external_scheme, strategy_in);
    _recombiner = recombiner_in;
  }

  /// a default constructor which creates a jet definition that is in
  /// a well-defined internal state, but not actually usable for jet
  /// clustering.
  JetDefinition()  {
    *this = JetDefinition(undefined_jet_algorithm, 1.0);
  }
  

  // /// a default constructor
  // JetDefinition() {
  //   *this = JetDefinition(kt_algorithm, 1.0);
  // }

  /// constructor based on a pointer to a user's plugin; the object
  /// pointed to must remain valid for the whole duration of existence
  /// of the JetDefinition and any related ClusterSequences
  JetDefinition(const Plugin * plugin_in) {
    _plugin = plugin_in;
    _strategy = plugin_strategy;
    _Rparam = _plugin->R();
    _jet_algorithm = plugin_algorithm;
    set_recombination_scheme(E_scheme);
  }


  /// constructor to fully specify a jet-definition (together with
  /// information about how algorithically to run it).
  ///
  /// the ordering of arguments here is old and deprecated (except
  /// as the common constructor for internal use)
  JetDefinition(JetAlgorithm jet_algorithm_in, 
                double R_in, 
                Strategy strategy_in,
                RecombinationScheme recomb_scheme_in = E_scheme,
                int nparameters_in = 1);

  /// cluster the supplied particles and returns a vector of resulting
  /// jets, sorted by pt (or energy in the case of spherical,
  /// i.e. e+e-, algorithms). This routine currently only makes
  /// sense for "inclusive" type algorithms.
  template <class L> 
  std::vector<PseudoJet> operator()(const std::vector<L> & particles) const;
  
  /// R values larger than max_allowable_R are not allowed.
  ///
  /// We use a value of 1000, substantially smaller than
  /// numeric_limits<double>::max(), to leave room for the convention
  /// within PseudoJet of setting unphysical (infinite) rapidities to
  /// +-(MaxRap + abs(pz())), where MaxRap is 10^5.
  static const double max_allowable_R; //= 1000.0;

  /// set the recombination scheme to the one provided
  void set_recombination_scheme(RecombinationScheme);

  /// set the recombiner class to the one provided
  ///
  /// Note that in order to associate to a jet definition a recombiner
  /// from another jet definition, it is strongly recommended to use
  /// the set_recombiner(const JetDefinition &) method below. The
  /// latter correctly handles the situations where the jet definition
  /// owns the recombiner (i.e. where delete_recombiner_when_unused
  /// has been called). In such cases, using set_recombiner(const
  /// Recombiner *) may lead to memory corruption.
  void set_recombiner(const Recombiner * recomb) {
    if (_shared_recombiner()) _shared_recombiner.reset(recomb);
    _recombiner = recomb;
    _default_recombiner = DefaultRecombiner(external_scheme);
  }

  /// set the recombiner to be the same as the one of 'other_jet_def'
  ///
  /// Note that this is the recommended method to associate to a jet
  /// definition the recombiner from another jet definition. Compared
  /// to the set_recombiner(const Recombiner *) above, it correctly
  /// handles the case where the jet definition owns the recombiner
  /// (i.e. where delete_recombiner_when_unused has been called)
  void set_recombiner(const JetDefinition &other_jet_def);

  /// calling this tells the JetDefinition to handle the deletion of
  /// the recombiner when it is no longer used. (Should not be called
  /// if the recombiner was initialised from a JetDef whose recombiner
  /// was already scheduled to delete itself - memory handling will
  /// already be automatic across both JetDef's in that case).
  void delete_recombiner_when_unused();

  /// return a pointer to the plugin 
  const Plugin * plugin() const {return _plugin;};

  /// calling this causes the JetDefinition to handle the deletion of the
  /// plugin when it is no longer used
  void delete_plugin_when_unused();

  /// return information about the definition...
  JetAlgorithm jet_algorithm  () const {return _jet_algorithm  ;}
  /// same as above for backward compatibility
  JetAlgorithm jet_finder     () const {return _jet_algorithm  ;}
  double    R           () const {return _Rparam      ;}
  // a general purpose extra parameter, whose meaning depends on
  // the algorithm, and may often be unused.
  double    extra_param () const {return _extra_param ;}
  Strategy  strategy    () const {return _strategy    ;}
  RecombinationScheme recombination_scheme() const {
    return _default_recombiner.scheme();}

  /// (re)set the jet finder
  void set_jet_algorithm(JetAlgorithm njf) {_jet_algorithm = njf;}
  /// same as above for backward compatibility
  void set_jet_finder(JetAlgorithm njf)    {_jet_algorithm = njf;}
  /// (re)set the general purpose extra parameter
  void set_extra_param(double xtra_param) {_extra_param = xtra_param;}

  /// returns a pointer to the currently defined recombiner. 
  ///
  /// Warning: the pointer may be to an internal recombiner (for
  /// default recombination schemes), in which case if the
  /// JetDefinition becomes invalid (e.g. is deleted), the pointer
  /// will then point to an object that no longer exists.
  /// 
  /// Note also that if you copy a JetDefinition with a default
  /// recombination scheme, then the two copies will have distinct
  /// recombiners, and return different recombiner() pointers.
  const Recombiner * recombiner() const {
    return _recombiner == 0 ? & _default_recombiner : _recombiner;}

  /// returns true if the current jet definitions shares the same
  /// recombiner as the one passed as an argument
  bool has_same_recombiner(const JetDefinition &other_jd) const;

  /// returns true if the jet definition involves an algorithm
  /// intended for use on a spherical geometry (e.g. e+e- algorithms,
  /// as opposed to most pp algorithms, which use a cylindrical,
  /// rapidity-phi geometry).
  bool is_spherical() const;

  /// return a textual description of the current jet definition 
  std::string description() const;

  /// returns a description not including the recombiner information
  std::string description_no_recombiner() const;

  /// a short textual description of the algorithm jet_alg
  static std::string algorithm_description(const JetAlgorithm jet_alg);

  /// the number of parameters associated to a given jet algorithm
  static unsigned int n_parameters_for_algorithm(const JetAlgorithm jet_alg);

public:
  //======================================================================
  /// @ingroup advanced_usage
  /// \class Recombiner
  /// An abstract base class that will provide the recombination scheme
  /// facilities and/or allow a user to extend these facilities
  class Recombiner {
  public:
    /// return a textual description of the recombination scheme
    /// implemented here
    virtual std::string description() const = 0;
    
    /// recombine pa and pb and put result into pab
    virtual void recombine(const PseudoJet & pa, const PseudoJet & pb, 
                           PseudoJet & pab) const = 0;

    /// routine called to preprocess each input jet (to make all input
    /// jets compatible with the scheme requirements (e.g. massless).
    virtual void preprocess(PseudoJet & ) const {};
    
    /// a destructor to be replaced if necessary in derived classes...
    virtual ~Recombiner() {};

    /// pa += pb in the given recombination scheme. Not virtual -- the
    /// user should have no reason to want to redefine this!
    inline void plus_equal(PseudoJet & pa, const PseudoJet & pb) const {
      // put result in a temporary location in case the recombiner
      // does something funny (ours doesn't, but who knows about the
      // user's)
      PseudoJet pres; 
      recombine(pa,pb,pres);
      pa = pres;
    }

  };
  
  
  //======================================================================
  /// @ingroup advanced_usage
  /// \class DefaultRecombiner
  /// A class that will provide the recombination scheme facilities and/or
  /// allow a user to extend these facilities
  ///
  /// This class is derived from the (abstract) class Recombiner. It
  /// simply "sums" PseudoJets using a specified recombination scheme
  /// (E-scheme by default)
  class DefaultRecombiner : public Recombiner {
  public:
    DefaultRecombiner(RecombinationScheme recomb_scheme = E_scheme) : 
      _recomb_scheme(recomb_scheme) {}
    
    virtual std::string description() const;
    
    /// recombine pa and pb and put result into pab
    virtual void recombine(const PseudoJet & pa, const PseudoJet & pb, 
                           PseudoJet & pab) const;

    virtual void preprocess(PseudoJet & p) const;

    /// return the index of the recombination scheme
    RecombinationScheme scheme() const {return _recomb_scheme;}
    
  private:
    RecombinationScheme _recomb_scheme;
  };


  //======================================================================
  /// @ingroup advanced_usage
  /// \class Plugin
  /// a class that allows a user to introduce their own "plugin" jet
  /// finder
  ///
  /// Note that all the plugins provided with FastJet are derived from
  /// this class
  class Plugin{
  public:
    /// return a textual description of the jet-definition implemented
    /// in this plugin
    virtual std::string description() const = 0;
    
    /// given a ClusterSequence that has been filled up with initial
    /// particles, the following function should fill up the rest of the
    /// ClusterSequence, using the following member functions of
    /// ClusterSequence:
    ///   - plugin_do_ij_recombination(...)
    ///   - plugin_do_iB_recombination(...)
    virtual void run_clustering(ClusterSequence &) const = 0;
    
    virtual double R() const = 0;
    
    /// return true if there is specific support for the measurement
    /// of passive areas, in the sense that areas determined from all
    /// particles below the ghost separation scale will be a passive
    /// area. [If you don't understand this, ignore it!]
    virtual bool supports_ghosted_passive_areas() const {return false;}

    /// set the ghost separation scale for passive area determinations
    /// in future runs (strictly speaking that makes the routine
    /// a non const, so related internal info must be stored as a mutable)
    virtual void set_ghost_separation_scale(double scale) const;
    virtual double ghost_separation_scale() const {return 0.0;}

    /// if this returns false then a warning will be given
    /// whenever the user requests "exclusive" jets from the
    /// cluster sequence
    virtual bool exclusive_sequence_meaningful() const {return false;}

    /// returns true if the plugin implements an algorithm intended
    /// for use on a spherical geometry (e.g. e+e- algorithms, as
    /// opposed to most pp algorithms, which use a cylindrical,
    /// rapidity-phi geometry).
    virtual bool is_spherical() const {return false;}

    /// a destructor to be replaced if necessary in derived classes...
    virtual ~Plugin() {};
  };

private:


  JetAlgorithm _jet_algorithm;
  double    _Rparam;
  double    _extra_param ; ///< parameter whose meaning varies according to context
  Strategy  _strategy  ;

  const Plugin * _plugin;
  SharedPtr<const Plugin> _plugin_shared;

  // when we use our own recombiner it's useful to point to it here
  // so that we don't have to worry about deleting it etc...
  DefaultRecombiner _default_recombiner;
  const Recombiner * _recombiner;
  SharedPtr<const Recombiner> _shared_recombiner;

};


//-------------------------------------------------------------------------------
// helper functions to build a jet made of pieces
//
// These functions include an options recombiner used to compute the
// total composite jet momentum
// -------------------------------------------------------------------------------

/// build a "CompositeJet" from the vector of its pieces
///
/// In this case, E-scheme recombination is assumed to compute the
/// total momentum
PseudoJet join(const std::vector<PseudoJet> & pieces, const JetDefinition::Recombiner & recombiner);

/// build a MergedJet from a single PseudoJet
PseudoJet join(const PseudoJet & j1, 
	       const JetDefinition::Recombiner & recombiner);

/// build a MergedJet from 2 PseudoJet
PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, 
	       const JetDefinition::Recombiner & recombiner);

/// build a MergedJet from 3 PseudoJet
PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, const PseudoJet & j3, 
	       const JetDefinition::Recombiner & recombiner);

/// build a MergedJet from 4 PseudoJet
PseudoJet join(const PseudoJet & j1, const PseudoJet & j2, const PseudoJet & j3, const PseudoJet & j4, 
	       const JetDefinition::Recombiner & recombiner);


FASTJET_END_NAMESPACE

// include ClusterSequence which includes the implementation of the 
// templated JetDefinition::operator()(...) member
#include "fastjet/ClusterSequence.hh"


#endif // __FASTJET_JETDEFINITION_HH__
