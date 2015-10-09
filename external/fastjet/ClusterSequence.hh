#ifndef __FASTJET_CLUSTERSEQUENCE_HH__
#define __FASTJET_CLUSTERSEQUENCE_HH__

//FJSTARTHEADER
// $Id: ClusterSequence.hh 3911 2015-07-02 12:09:58Z salam $
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


#include<vector>
#include<map>
#include "fastjet/PseudoJet.hh"
#include<memory>
#include<cassert>
#include<iostream>
#include<string>
#include<set>
#include<cmath> // needed to get double std::abs(double)
#include "fastjet/Error.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/SharedPtr.hh"
#include "fastjet/LimitedWarning.hh"
#include "fastjet/FunctionOfPseudoJet.hh"
#include "fastjet/ClusterSequenceStructure.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


// forward declarations
class ClusterSequenceStructure;
class DynamicNearestNeighbours;

/// @ingroup basic_classes
/// \class ClusterSequence
/// deals with clustering
class ClusterSequence {


 public: 

  /// default constructor
  ClusterSequence () : _deletes_self_when_unused(false) {}

  /// create a ClusterSequence, starting from the supplied set
  /// of PseudoJets and clustering them with jet definition specified
  /// by jet_def (which also specifies the clustering strategy)
  template<class L> ClusterSequence (
			          const std::vector<L> & pseudojets,
				  const JetDefinition & jet_def,
				  const bool & writeout_combinations = false);
  
  /// copy constructor for a ClusterSequence
  ClusterSequence (const ClusterSequence & cs) : _deletes_self_when_unused(false) {
    transfer_from_sequence(cs);
  }

  // virtual ClusterSequence destructor, in case any derived class
  // thinks of needing a destructor at some point
  virtual ~ClusterSequence (); //{}

  // NB: in the routines that follow, for extracting lists of jets, a
  //     list structure might be more efficient, if sometimes a little
  //     more awkward to use (at least for old fortran hands).

  /// return a vector of all jets (in the sense of the inclusive
  /// algorithm) with pt >= ptmin. Time taken should be of the order
  /// of the number of jets returned.
  std::vector<PseudoJet> inclusive_jets (const double ptmin = 0.0) const;

  /// return the number of jets (in the sense of the exclusive
  /// algorithm) that would be obtained when running the algorithm
  /// with the given dcut.
  int n_exclusive_jets (const double dcut) const;

  /// return a vector of all jets (in the sense of the exclusive
  /// algorithm) that would be obtained when running the algorithm
  /// with the given dcut.
  std::vector<PseudoJet> exclusive_jets (const double dcut) const;

  /// return a vector of all jets when the event is clustered (in the
  /// exclusive sense) to exactly njets. 
  ///
  /// If there are fewer than njets particles in the ClusterSequence
  /// an error is thrown
  std::vector<PseudoJet> exclusive_jets (const int njets) const;

  /// return a vector of all jets when the event is clustered (in the
  /// exclusive sense) to exactly njets. 
  ///
  /// If there are fewer than njets particles in the ClusterSequence
  /// the function just returns however many particles there were.
  std::vector<PseudoJet> exclusive_jets_up_to (const int njets) const;

  /// return the dmin corresponding to the recombination that went
  /// from n+1 to n jets (sometimes known as d_{n n+1}). If the number
  /// of particles in the event is <= njets, the function returns 0.
  double exclusive_dmerge (const int njets) const;

  /// return the maximum of the dmin encountered during all recombinations 
  /// up to the one that led to an n-jet final state; identical to
  /// exclusive_dmerge, except in cases where the dmin do not increase
  /// monotonically.
  double exclusive_dmerge_max (const int njets) const;

  /// return the ymin corresponding to the recombination that went from
  /// n+1 to n jets (sometimes known as y_{n n+1}).
  double exclusive_ymerge (int njets) const {return exclusive_dmerge(njets) / Q2();}

  /// same as exclusive_dmerge_max, but normalised to squared total energy
  double exclusive_ymerge_max (int njets) const {return exclusive_dmerge_max(njets)/Q2();}

  /// the number of exclusive jets at the given ycut
  int n_exclusive_jets_ycut (double ycut) const {return n_exclusive_jets(ycut*Q2());}

  /// the exclusive jets obtained at the given ycut
  std::vector<PseudoJet> exclusive_jets_ycut (double ycut) const {
    int njets = n_exclusive_jets_ycut(ycut);
    return exclusive_jets(njets);
  }


  //int n_exclusive_jets (const PseudoJet & jet, const double dcut) const;

  /// return a vector of all subjets of the current jet (in the sense
  /// of the exclusive algorithm) that would be obtained when running
  /// the algorithm with the given dcut. 
  ///
  /// Time taken is O(m ln m), where m is the number of subjets that
  /// are found. If m gets to be of order of the total number of
  /// constituents in the jet, this could be substantially slower than
  /// just getting that list of constituents.
  std::vector<PseudoJet> exclusive_subjets (const PseudoJet & jet, 
                                            const double dcut) const;

  /// return the size of exclusive_subjets(...); still n ln n with same
  /// coefficient, but marginally more efficient than manually taking
  /// exclusive_subjets.size()
  int n_exclusive_subjets(const PseudoJet & jet, 
                          const double dcut) const;

  /// return the list of subjets obtained by unclustering the supplied
  /// jet down to nsub subjets. Throws an error if there are fewer than
  /// nsub particles in the jet.
  ///
  /// This requires nsub ln nsub time
  std::vector<PseudoJet> exclusive_subjets (const PseudoJet & jet, 
                                            int nsub) const;

  /// return the list of subjets obtained by unclustering the supplied
  /// jet down to nsub subjets (or all constituents if there are fewer
  /// than nsub).
  ///
  /// This requires nsub ln nsub time
  std::vector<PseudoJet> exclusive_subjets_up_to (const PseudoJet & jet, 
						  int nsub) const;

  /// returns the dij that was present in the merging nsub+1 -> nsub 
  /// subjets inside this jet.
  ///
  /// Returns 0 if there were nsub or fewer constituents in the jet.
  double exclusive_subdmerge(const PseudoJet & jet, int nsub) const;

  /// returns the maximum dij that occurred in the whole event at the
  /// stage that the nsub+1 -> nsub merge of subjets occurred inside 
  /// this jet.
  ///
  /// Returns 0 if there were nsub or fewer constituents in the jet.
  double exclusive_subdmerge_max(const PseudoJet & jet, int nsub) const;

  //std::vector<PseudoJet> exclusive_jets (const PseudoJet & jet, 
  //                                       const int njets) const;
  //double exclusive_dmerge (const PseudoJet & jet, const int njets) const;

  /// returns the sum of all energies in the event (relevant mainly for e+e-)
  double Q() const {return _Qtot;}
  /// return Q()^2
  double Q2() const {return _Qtot*_Qtot;}

  /// returns true iff the object is included in the jet. 
  ///
  /// NB: this is only sensible if the object is already registered
  /// within the cluster sequence, so you cannot use it with an input
  /// particle to the CS (since the particle won't have the history
  /// index set properly).
  ///
  /// For nice clustering structures it should run in O(ln(N)) time
  /// but in worst cases (certain cone plugins) it can take O(n) time,
  /// where n is the number of particles in the jet.
  bool object_in_jet(const PseudoJet & object, const PseudoJet & jet) const;

  /// if the jet has parents in the clustering, it returns true
  /// and sets parent1 and parent2 equal to them.
  ///
  /// if it has no parents it returns false and sets parent1 and
  /// parent2 to zero
  bool has_parents(const PseudoJet & jet, PseudoJet & parent1, 
               PseudoJet & parent2) const;

  /// if the jet has a child then return true and give the child jet
  /// otherwise return false and set the child to zero
  bool has_child(const PseudoJet & jet, PseudoJet & child) const;

  /// Version of has_child that sets a pointer to the child if the child
  /// exists;
  bool has_child(const PseudoJet & jet, const PseudoJet * & childp) const;

  /// if this jet has a child (and so a partner) return true
  /// and give the partner, otherwise return false and set the
  /// partner to zero
  bool has_partner(const PseudoJet & jet, PseudoJet & partner) const;

  
  /// return a vector of the particles that make up jet
  std::vector<PseudoJet> constituents (const PseudoJet & jet) const;


  /// output the supplied vector of jets in a format that can be read
  /// by an appropriate root script; the format is:
  /// jet-n jet-px jet-py jet-pz jet-E 
  ///   particle-n particle-rap particle-phi particle-pt
  ///   particle-n particle-rap particle-phi particle-pt
  ///   ...
  /// #END
  /// ... [i.e. above repeated]
  void print_jets_for_root(const std::vector<PseudoJet> & jets, 
                           std::ostream & ostr = std::cout) const;

  /// print jets for root to the file labelled filename, with an
  /// optional comment at the beginning
  void print_jets_for_root(const std::vector<PseudoJet> & jets, 
                           const std::string & filename,
			   const std::string & comment = "") const;

// Not yet. Perhaps in a future release.
//   /// print out all inclusive jets with pt > ptmin
//   virtual void print_jets (const double ptmin=0.0) const;

  /// add on to subjet_vector the constituents of jet (for internal use mainly)
  void add_constituents (const PseudoJet & jet, 
			 std::vector<PseudoJet> & subjet_vector) const;

  /// return the enum value of the strategy used to cluster the event
  inline Strategy strategy_used () const {return _strategy;}

  /// return the name of the strategy used to cluster the event
  std::string strategy_string () const {return strategy_string(_strategy);}

  /// return the name of the strategy associated with the enum strategy_in
  std::string strategy_string (Strategy strategy_in) const;


  /// return a reference to the jet definition
  const JetDefinition & jet_def() const {return _jet_def;}

  /// by calling this routine you tell the ClusterSequence to delete
  /// itself when all the Pseudojets associated with it have gone out
  /// of scope. 
  ///
  /// At the time you call this, there must be at least one jet or
  /// other object outside the CS that is associated with the CS
  /// (e.g. the result of inclusive_jets()).
  ///
  /// NB: after having made this call, the user is still allowed to
  /// delete the CS. Jets associated with it will then simply not be
  /// able to access their substructure after that point.
  void delete_self_when_unused();

  /// return true if the object has been told to delete itself
  /// when unused
  bool will_delete_self_when_unused() const {return _deletes_self_when_unused;}

  /// tell the ClusterSequence it's about to be self deleted (internal use only)
  void signal_imminent_self_deletion() const;

  /// returns the scale associated with a jet as required for this
  /// clustering algorithm (kt^2 for the kt-algorithm, 1 for the
  /// Cambridge algorithm). Intended mainly for internal use and not
  /// valid for plugin algorithms.
  double jet_scale_for_algorithm(const PseudoJet & jet) const;

  ///

  //----- next follow functions designed specifically for plugins, which
  //      may only be called when plugin_activated() returns true

  /// record the fact that there has been a recombination between
  /// jets()[jet_i] and jets()[jet_k], with the specified dij, and
  /// return the index (newjet_k) allocated to the new jet, whose
  /// momentum is assumed to be the 4-vector sum of that of jet_i and
  /// jet_j
  void plugin_record_ij_recombination(int jet_i, int jet_j, double dij, 
				      int & newjet_k) {
    assert(plugin_activated());
    _do_ij_recombination_step(jet_i, jet_j, dij, newjet_k);
  }

  /// as for the simpler variant of plugin_record_ij_recombination,
  /// except that the new jet is attributed the momentum and
  /// user_index of newjet
  void plugin_record_ij_recombination(int jet_i, int jet_j, double dij, 
				      const PseudoJet & newjet, 
				      int & newjet_k);

  /// record the fact that there has been a recombination between
  /// jets()[jet_i] and the beam, with the specified diB; when looking
  /// for inclusive jets, any iB recombination will returned to the user 
  /// as a jet.
  void plugin_record_iB_recombination(int jet_i, double diB) {
    assert(plugin_activated());
    _do_iB_recombination_step(jet_i, diB);
  }

  /// @ingroup extra_info
  /// \class Extras
  /// base class to store extra information that plugins may provide
  /// 
  /// a class intended to serve as a base in case a plugin needs to
  /// associate extra information with a ClusterSequence (see
  /// SISConePlugin.* for an example).
  class Extras {
  public:
    virtual ~Extras() {}
    virtual std::string description() const {return "This is a dummy extras class that contains no extra information! Derive from it if you want to use it to provide extra information from a plugin jet finder";}
  };

  /// the plugin can associate some extra information with the
  /// ClusterSequence object by calling this function. The
  /// ClusterSequence takes ownership of the pointer (and
  /// responsibility for deleting it when the CS gets deleted).
  inline void plugin_associate_extras(Extras * extras_in) {
    _extras.reset(extras_in);
  }

  /// the plugin can associate some extra information with the
  /// ClusterSequence object by calling this function
  /// 
  /// As of FJ v3.1, this is deprecated, in line with the deprecation
  /// of auto_ptr in C++11
  inline void plugin_associate_extras(std::auto_ptr<Extras> extras_in) {
    _extras.reset(extras_in.release());
  }

  /// returns true when the plugin is allowed to run the show.
  inline bool plugin_activated() const {return _plugin_activated;}

  /// returns a pointer to the extras object (may be null)
  const Extras * extras() const {return _extras.get();}

  /// allows a plugin to run a templated clustering (nearest-neighbour heuristic)
  ///
  /// This has N^2 behaviour on "good" distance, but a worst case behaviour
  /// of N^3 (and many algs trigger the worst case behaviour)
  ///
  /// 
  /// For more details on how this works, see GenBriefJet below
  template<class GBJ> void plugin_simple_N2_cluster () {
    assert(plugin_activated());
    _simple_N2_cluster<GBJ>();
  }


public:
//DEP   /// set the default (static) jet finder across all current and future
//DEP   /// ClusterSequence objects -- deprecated and obsolescent (i.e. may be
//DEP   /// suppressed in a future release).
//DEP   static void set_jet_algorithm (JetAlgorithm jet_algorithm) {_default_jet_algorithm = jet_algorithm;}
//DEP   /// same as above for backward compatibility
//DEP   static void set_jet_finder (JetAlgorithm jet_algorithm)    {_default_jet_algorithm = jet_algorithm;}


  /// \ingroup extra_info
  /// \struct history_element
  /// a single element in the clustering history
  /// 
  /// (see vector _history below).
  struct history_element{
    int parent1; /// index in _history where first parent of this jet
                 /// was created (InexistentParent if this jet is an
                 /// original particle)

    int parent2; /// index in _history where second parent of this jet
                 /// was created (InexistentParent if this jet is an
                 /// original particle); BeamJet if this history entry
                 /// just labels the fact that the jet has recombined
                 /// with the beam)

    int child;   /// index in _history where the current jet is
		 /// recombined with another jet to form its child. It
		 /// is Invalid if this jet does not further
		 /// recombine.

    int jetp_index; /// index in the _jets vector where we will find the
                 /// PseudoJet object corresponding to this jet
                 /// (i.e. the jet created at this entry of the
                 /// history). NB: if this element of the history
                 /// corresponds to a beam recombination, then
                 /// jetp_index=Invalid.

    double dij;  /// the distance corresponding to the recombination
		 /// at this stage of the clustering.

    double max_dij_so_far; /// the largest recombination distance seen
			   /// so far in the clustering history.
  };

  enum JetType {Invalid=-3, InexistentParent = -2, BeamJet = -1};

  /// allow the user to access the internally stored _jets() array,
  /// which contains both the initial particles and the various
  /// intermediate and final stages of recombination.
  ///
  /// The first n_particles() entries are the original particles,
  /// in the order in which they were supplied to the ClusterSequence
  /// constructor. It can be useful to access them for example when
  /// examining whether a given input object is part of a specific
  /// jet, via the objects_in_jet(...) member function (which only takes
  /// PseudoJets that are registered in the ClusterSequence).
  ///
  /// One of the other (internal uses) is related to the fact
  /// because we don't seem to be able to access protected elements of
  /// the class for an object that is not "this" (at least in case where
  /// "this" is of a slightly different kind from the object, both
  /// derived from ClusterSequence).
  const std::vector<PseudoJet> & jets()    const;

  /// allow the user to access the raw internal history.
  ///
  /// This is present (as for jets()) in part so that protected
  /// derived classes can access this information about other
  /// ClusterSequences.
  ///
  /// A user who wishes to follow the details of the ClusterSequence
  /// can also make use of this information (and should consult the
  /// history_element documentation for more information), but should
  /// be aware that these internal structures may evolve in future
  /// FastJet versions.
  const std::vector<history_element> & history() const;

  /// returns the number of particles that were provided to the
  /// clustering algorithm (helps the user find their way around the
  /// history and jets objects if they weren't paying attention
  /// beforehand).
  unsigned int n_particles() const;

  /// returns a vector of size n_particles() which indicates, for 
  /// each of the initial particles (in the order in which they were
  /// supplied), which of the supplied jets it belongs to; if it does
  /// not belong to any of the supplied jets, the index is set to -1;
  std::vector<int> particle_jet_indices(const std::vector<PseudoJet> &) const;

  /// routine that returns an order in which to read the history
  /// such that clusterings that lead to identical jet compositions
  /// but different histories (because of degeneracies in the
  /// clustering order) will have matching constituents for each
  /// matching entry in the unique_history_order.
  ///
  /// The order has the property that an entry's parents will always
  /// appear prior to that entry itself. 
  ///
  /// Roughly speaking the order is such that we first provide all
  /// steps that lead to the final jet containing particle 1; then we
  /// have the steps that lead to reconstruction of the jet containing
  /// the next-lowest-numbered unclustered particle, etc...
  /// [see GPS CCN28-12 for more info -- of course a full explanation
  /// here would be better...]
  std::vector<int> unique_history_order() const;

  /// return the set of particles that have not been clustered. For 
  /// kt and cam/aachen algorithms this should always be null, but for
  /// cone type algorithms it can be non-null;
  std::vector<PseudoJet> unclustered_particles() const;

  /// Return the list of pseudojets in the ClusterSequence that do not
  /// have children (and are not among the inclusive jets). They may
  /// result from a clustering step or may be one of the pseudojets
  /// returned by unclustered_particles().
  std::vector<PseudoJet> childless_pseudojets() const;

  /// returns true if the object (jet or particle) is contained by (ie
  /// belongs to) this cluster sequence.
  ///
  /// Tests performed: if thejet's interface is this cluster sequence
  /// and its cluster history index is in a consistent range.
  bool contains(const PseudoJet & object) const;

  /// transfer the sequence contained in other_seq into our own;
  /// any plugin "extras" contained in the from_seq will be lost
  /// from there.
  ///
  /// It also sets the ClusterSequence pointers of the PseudoJets in
  /// the history to point to this ClusterSequence
  ///
  /// When specified, the second argument is an action that will be
  /// applied on every jets in the resulting ClusterSequence
  void transfer_from_sequence(const ClusterSequence & from_seq,
			      const FunctionOfPseudoJet<PseudoJet> * action_on_jets = 0);

  /// retrieve a shared pointer to the wrapper to this ClusterSequence
  ///
  /// this may turn useful if you want to track when this
  /// ClusterSequence goes out of scope
  const SharedPtr<PseudoJetStructureBase> & structure_shared_ptr() const{
    return _structure_shared_ptr;
  }

  /// the structure type associated with a jet belonging to a ClusterSequence
  typedef ClusterSequenceStructure StructureType;

  /// This is the function that is automatically called during
  /// clustering to print the FastJet banner. Only the first call to
  /// this function will result in the printout of the banner. Users
  /// may wish to call this function themselves, during the
  /// initialization phase of their program, in order to ensure that
  /// the banner appears before other output. This call will not
  /// affect 3rd-party banners, e.g. those from plugins.
  static void print_banner();

  /// \cond internal_doc
  //  [this line must be left as is to hide the doxygen comment]
  /// A call to this function modifies the stream used to print
  /// banners (by default cout). If a null pointer is passed, banner
  /// printout is suppressed. This affects all banners, including
  /// those from plugins.
  ///
  /// Please note that if you distribute 3rd party code
  /// that links with FastJet, that 3rd party code must not
  /// use this call turn off the printing of FastJet banners
  /// by default. This requirement reflects the spirit of
  /// clause 2c of the GNU Public License (v2), under which
  /// FastJet and its plugins are distributed.
  static void set_fastjet_banner_stream(std::ostream * ostr) {_fastjet_banner_ostr = ostr;}
  //  [this line must be left as is to hide the doxygen comment]
  /// \endcond

  /// returns a pointer to the stream to be used to print banners
  /// (cout by default). This function is used by plugins to determine
  /// where to direct their banners. Plugins should properly handle
  /// the case where the pointer is null.
  static std::ostream * fastjet_banner_stream() {return _fastjet_banner_ostr;}

private:
  /// \cond internal_doc

  /// contains the actual stream to use for banners 
  static std::ostream * _fastjet_banner_ostr;

  /// \endcond

protected:
//DEP  static JetAlgorithm _default_jet_algorithm;
  JetDefinition _jet_def;

  /// transfer the vector<L> of input jets into our own vector<PseudoJet>
  /// _jets (with some reserved space for future growth).
  template<class L> void _transfer_input_jets(
                                     const std::vector<L> & pseudojets);

  /// This is what is called to do all the initialisation and
  /// then run the clustering (may be called by various constructors).
  /// It assumes _jets contains the momenta to be clustered.
  void _initialise_and_run (const JetDefinition & jet_def,
			    const bool & writeout_combinations);

  //// this performs the initialisation, minus the option-decanting
  //// stage; for low multiplicity, initialising a few things in the
  //// constructor, calling the decant_options_partial() and then this
  //// is faster than going through _initialise_and_run.
  void _initialise_and_run_no_decant();

//DEP   /// This is an alternative routine for initialising and running the
//DEP   /// clustering, provided for legacy purposes. The jet finder is that
//DEP   /// specified in the static member _default_jet_algorithm.
//DEP   void _initialise_and_run (const double R,
//DEP 			    const Strategy & strategy,
//DEP 			    const bool & writeout_combinations);

  /// fills in the various member variables with "decanted" options from
  /// the jet_definition and writeout_combinations variables
  void _decant_options(const JetDefinition & jet_def,
                       const bool & writeout_combinations);

  /// assuming that the jet definition, writeout_combinations and
  /// _structure_shared_ptr have been set (e.g. in an initialiser list
  /// in the constructor), it handles the remaining decanting of
  /// options.
  void _decant_options_partial();

  /// fill out the history (and jet cross refs) related to the initial
  /// set of jets (assumed already to have been "transferred"),
  /// without any clustering
  void _fill_initial_history();

  /// carry out the recombination between the jets numbered jet_i and
  /// jet_j, at distance scale dij; return the index newjet_k of the
  /// result of the recombination of i and j.
  void _do_ij_recombination_step(const int jet_i, const int jet_j, 
				 const double dij, int & newjet_k);

  /// carry out an recombination step in which _jets[jet_i] merges with
  /// the beam, 
  void _do_iB_recombination_step(const int jet_i, const double diB);

  /// every time a jet is added internally during clustering, this
  /// should be called to set the jet's structure shared ptr to point
  /// to the CS (and the count of internally associated objects is
  /// also updated). This should not be called outside construction of
  /// a CS object.
  void _set_structure_shared_ptr(PseudoJet & j);

  /// make sure that the CS's internal tally of the use count matches
  /// that of the _structure_shared_ptr
  void _update_structure_use_count();
  
  /// returns a suggestion for the best strategy to use on event
  /// multiplicity, algorithm, R, etc.
  Strategy _best_strategy() const;
  
  /// \if internal_doc
  /// \class _Parabola
  /// returns c*(a*R**2 + b*R + 1);
  /// Written as a class in case we want to give names to different
  /// parabolas
  /// \endif
  class _Parabola {
  public:
    _Parabola(double a, double b, double c) : _a(a), _b(b), _c(c) {}
    inline double operator()(const double R) const {return _c*(_a*R*R + _b*R + 1);}
  private:
    double _a, _b, _c;
  };

  /// \if internal_doc
  /// \class _Line
  /// operator()(R) returns a*R+b;
  /// \endif
  class _Line {
  public:
    _Line(double a, double b) : _a(a), _b(b) {}
    inline double operator()(const double R) const {return _a*R + _b;}
  private:
    double _a, _b;
  };

  /// This contains the physical PseudoJets; for each PseudoJet one
  /// can find the corresponding position in the _history by looking
  /// at _jets[i].cluster_hist_index().
  std::vector<PseudoJet> _jets;


  /// this vector will contain the branching history; for each stage,
  /// _history[i].jetp_index indicates where to look in the _jets
  /// vector to get the physical PseudoJet.
  std::vector<history_element> _history;

  /// set subhist to be a set pointers to history entries corresponding to the
  /// subjets of this jet; one stops going working down through the
  /// subjets either when 
  ///   - there is no further to go
  ///   - one has found maxjet entries
  ///   - max_dij_so_far <= dcut
  /// By setting maxjet=0 one can use just dcut; by setting dcut<0
  /// one can use jet maxjet
  void get_subhist_set(std::set<const history_element*> & subhist,
                       const  PseudoJet & jet, double dcut, int maxjet) const;

  bool _writeout_combinations;
  int  _initial_n;
  double _Rparam, _R2, _invR2;
  double _Qtot;
  Strategy    _strategy;
  JetAlgorithm  _jet_algorithm;

  SharedPtr<PseudoJetStructureBase> _structure_shared_ptr; //< will actually be of type ClusterSequenceStructure
  int _structure_use_count_after_construction; //< info of use when CS handles its own memory
  /// if true then the CS will delete itself when the last external
  /// object referring to it disappears. It is mutable so as to ensure
  /// that signal_imminent_self_deletion() [const] can make relevant
  /// changes.
  mutable bool _deletes_self_when_unused;

 private:

  bool _plugin_activated;
  SharedPtr<Extras> _extras; // things the plugin might want to add

  void _really_dumb_cluster ();
  void _delaunay_cluster ();
  //void _simple_N2_cluster ();
  template<class BJ> void _simple_N2_cluster ();
  void _tiled_N2_cluster ();
  void _faster_tiled_N2_cluster ();

  //
  void _minheap_faster_tiled_N2_cluster();

  // things needed specifically for Cambridge with Chan's 2D closest
  // pairs method
  void _CP2DChan_cluster();
  void _CP2DChan_cluster_2pi2R ();
  void _CP2DChan_cluster_2piMultD ();
  void _CP2DChan_limited_cluster(double D);
  void _do_Cambridge_inclusive_jets();

  // NSqrtN method for C/A
  void _fast_NsqrtN_cluster();

  void _add_step_to_history(const int step_number, const int parent1, 
			       const int parent2, const int jetp_index,
			       const double dij);

  /// internal routine associated with the construction of the unique
  /// history order (following children in the tree)
  void _extract_tree_children(int pos, std::valarray<bool> &, 
		const std::valarray<int> &, std::vector<int> &) const;

  /// internal routine associated with the construction of the unique
  /// history order (following parents in the tree)
  void _extract_tree_parents (int pos, std::valarray<bool> &, 
                const std::valarray<int> &,  std::vector<int> &) const;


  // these will be useful shorthands in the Voronoi-based code
  typedef std::pair<int,int> TwoVertices;
  typedef std::pair<double,TwoVertices> DijEntry;
  typedef std::multimap<double,TwoVertices> DistMap;

  /// currently used only in the Voronoi based code
  void _add_ktdistance_to_map(const int ii, 
			      DistMap & DijMap,
  			      const DynamicNearestNeighbours * DNN);


  /// will be set by default to be true for the first run
  static bool _first_time;

  /// manage warnings related to exclusive jets access
  static LimitedWarning _exclusive_warnings;

  /// the limited warning member for notification of user that 
  /// their requested strategy has been overridden (usually because
  /// they have R>2pi and not all strategies work then)
  static LimitedWarning _changed_strategy_warning;

  //----------------------------------------------------------------------
  /// the fundamental structure which contains the minimal info about
  /// a jet, as needed for our plain N^2 algorithm -- the idea is to
  /// put all info that will be accessed N^2 times into an array of
  /// BriefJets...
  struct BriefJet {
    double     eta, phi, kt2, NN_dist;
    BriefJet * NN;
    int        _jets_index;
  };

  /// structure analogous to BriefJet, but with the extra information
  /// needed for dealing with tiles
  class TiledJet {
  public:
    double     eta, phi, kt2, NN_dist;
    TiledJet * NN, *previous, * next; 
    int        _jets_index, tile_index, diJ_posn;
    // routines that are useful in the minheap version of tiled
    // clustering ("misuse" the otherwise unused diJ_posn, so as
    // to indicate whether jets need to have their minheap entries
    // updated).
    inline void label_minheap_update_needed() {diJ_posn = 1;}
    inline void label_minheap_update_done()   {diJ_posn = 0;}
    inline bool minheap_update_needed() const {return diJ_posn==1;}
  };

  //-- some of the functions that follow are templates and will work
  //as well for briefjet and tiled jets

  /// set the kinematic and labelling info for jeta so that it corresponds
  /// to _jets[_jets_index]
  template <class J> void _bj_set_jetinfo( J * const jet, 
						 const int _jets_index) const;

  /// "remove" this jet, which implies updating links of neighbours and
  /// perhaps modifying the tile structure
  void _bj_remove_from_tiles( TiledJet * const jet) const;

  /// return the distance between two BriefJet objects
  template <class J> double _bj_dist(const J * const jeta, 
			const J * const jetb) const;

  // return the diJ (multiplied by _R2) for this jet assuming its NN
  // info is correct
  template <class J> double _bj_diJ(const J * const jeta) const;

  /// for testing purposes only: if in the range head--tail-1 there is a
  /// a jet which corresponds to hist_index in the history, then
  /// return a pointer to that jet; otherwise return tail.
  template <class J> inline J * _bj_of_hindex(
                          const int hist_index, 
			  J * const head, J * const tail) 
    const {
    J * res;
    for(res = head; res<tail; res++) {
      if (_jets[res->_jets_index].cluster_hist_index() == hist_index) {break;}
    }
    return res;
  }


  //-- remaining functions are different in various cases, so we
  //   will use templates but are not sure if they're useful...

  /// updates (only towards smaller distances) the NN for jeta without checking
  /// whether in the process jeta itself might be a new NN of one of
  /// the jets being scanned -- span the range head to tail-1 with
  /// assumption that jeta is not contained in that range
  template <class J> void _bj_set_NN_nocross(J * const jeta, 
            J * const head, const J * const tail) const;

  /// reset the NN for jeta and DO check whether in the process jeta
  /// itself might be a new NN of one of the jets being scanned --
  /// span the range head to tail-1 with assumption that jeta is not
  /// contained in that range
  template <class J> void _bj_set_NN_crosscheck(J * const jeta, 
            J * const head, const J * const tail) const;
  


  /// number of neighbours that a tile will have (rectangular geometry
  /// gives 9 neighbours).
  static const int n_tile_neighbours = 9;
  //----------------------------------------------------------------------
  /// The fundamental structures to be used for the tiled N^2 algorithm
  /// (see CCN27-44 for some discussion of pattern of tiling)
  struct Tile {
    /// pointers to neighbouring tiles, including self
    Tile *   begin_tiles[n_tile_neighbours]; 
    /// neighbouring tiles, excluding self
    Tile **  surrounding_tiles; 
    /// half of neighbouring tiles, no self
    Tile **  RH_tiles;  
    /// just beyond end of tiles
    Tile **  end_tiles; 
    /// start of list of BriefJets contained in this tile
    TiledJet * head;    
    /// sometimes useful to be able to tag a tile
    bool     tagged;    
  };
  std::vector<Tile> _tiles;
  double _tiles_eta_min, _tiles_eta_max;
  double _tile_size_eta, _tile_size_phi;
  int    _n_tiles_phi,_tiles_ieta_min,_tiles_ieta_max;

  // reasonably robust return of tile index given ieta and iphi, in particular
  // it works even if iphi is negative
  inline int _tile_index (int ieta, int iphi) const {
    // note that (-1)%n = -1 so that we have to add _n_tiles_phi
    // before performing modulo operation
    return (ieta-_tiles_ieta_min)*_n_tiles_phi
                  + (iphi+_n_tiles_phi) % _n_tiles_phi;
  }

  // routines for tiled case, including some overloads of the plain
  // BriefJet cases
  int  _tile_index(const double eta, const double phi) const;
  void _tj_set_jetinfo ( TiledJet * const jet, const int _jets_index);
  void  _bj_remove_from_tiles(TiledJet * const jet);
  void _initialise_tiles();
  void _print_tiles(TiledJet * briefjets ) const;
  void _add_neighbours_to_tile_union(const int tile_index, 
		 std::vector<int> & tile_union, int & n_near_tiles) const;
  void _add_untagged_neighbours_to_tile_union(const int tile_index, 
		 std::vector<int> & tile_union, int & n_near_tiles);

  //----------------------------------------------------------------------
  /// fundamental structure for e+e- clustering
  struct EEBriefJet {
    double NN_dist;  // obligatorily present
    double kt2;      // obligatorily present == E^2 in general
    EEBriefJet * NN; // must be present too
    int    _jets_index; // must also be present!
    //...........................................................
    double nx, ny, nz;  // our internal storage for fast distance calcs
  };

  /// to help instantiation (fj 2.4.0; did not quite work on gcc 33 and os x 10.3?)
  //void _dummy_N2_cluster_instantiation();


  /// to avoid issues with template instantiation (OS X 10.3, gcc 3.3)
  void _simple_N2_cluster_BriefJet();
  /// to avoid issues with template instantiation (OS X 10.3, gcc 3.3)
  void _simple_N2_cluster_EEBriefJet();
};


//**********************************************************************
//**************    START   OF   INLINE   MATERIAL    ******************
//**********************************************************************


//----------------------------------------------------------------------
// Transfer the initial jets into our internal structure
template<class L> void ClusterSequence::_transfer_input_jets(
                                       const std::vector<L> & pseudojets) {

  // this will ensure that we can point to jets without difficulties
  // arising.
  _jets.reserve(pseudojets.size()*2);

  // insert initial jets this way so that any type L that can be
  // converted to a pseudojet will work fine (basically PseudoJet
  // and any type that has [] subscript access to the momentum
  // components, such as CLHEP HepLorentzVector).
  for (unsigned int i = 0; i < pseudojets.size(); i++) {
    _jets.push_back(pseudojets[i]);}
  
}

// //----------------------------------------------------------------------
// // initialise from some generic type... Has to be made available
// // here in order for it the template aspect of it to work...
// template<class L> ClusterSequence::ClusterSequence (
// 			          const std::vector<L> & pseudojets,
// 				  const double R,
// 				  const Strategy & strategy,
// 				  const bool & writeout_combinations) {
// 
//   // transfer the initial jets (type L) into our own array
//   _transfer_input_jets(pseudojets);
// 
//   // run the clustering
//   _initialise_and_run(R,strategy,writeout_combinations);
// }


//----------------------------------------------------------------------
/// constructor of a jet-clustering sequence from a vector of
/// four-momenta, with the jet definition specified by jet_def
template<class L> ClusterSequence::ClusterSequence (
			          const std::vector<L> & pseudojets,
				  const JetDefinition & jet_def_in,
				  const bool & writeout_combinations) :
  _jet_def(jet_def_in), _writeout_combinations(writeout_combinations),
  _structure_shared_ptr(new ClusterSequenceStructure(this))
{

  // transfer the initial jets (type L) into our own array
  _transfer_input_jets(pseudojets);

  // transfer the remaining options
  _decant_options_partial();

  // run the clustering
  _initialise_and_run_no_decant();
}


inline const std::vector<PseudoJet> & ClusterSequence::jets () const {
  return _jets;
}

inline const std::vector<ClusterSequence::history_element> & ClusterSequence::history () const {
  return _history;
}

inline unsigned int ClusterSequence::n_particles() const {return _initial_n;}

//----------------------------------------------------------------------
// implementation of JetDefinition::operator() is here to avoid nasty
// issues of order of implementations and includes
#ifndef __CINT__
template<class L>
std::vector<PseudoJet> JetDefinition::operator()(const std::vector<L> & particles) const {
  // create a new cluster sequence
  ClusterSequence * cs = new ClusterSequence(particles, *this);

  // get the jets, and sort them according to whether the algorithm
  // is spherical or not
  std::vector<PseudoJet> jets;
  if (is_spherical()) {
    jets = sorted_by_E(cs->inclusive_jets());
  } else {
    jets = sorted_by_pt(cs->inclusive_jets());
  }
  
  // make sure the ClusterSequence gets deleted once it's no longer
  // needed
  if (jets.size() != 0) {
    cs->delete_self_when_unused();
  } else {
    delete cs;
  }

  return jets;
}
#endif // __CINT__


//----------------------------------------------------------------------
template <class J> inline void ClusterSequence::_bj_set_jetinfo(
                            J * const jetA, const int _jets_index) const {
    jetA->eta  = _jets[_jets_index].rap();
    jetA->phi  = _jets[_jets_index].phi_02pi();
    jetA->kt2  = jet_scale_for_algorithm(_jets[_jets_index]);
    jetA->_jets_index = _jets_index;
    // initialise NN info as well
    jetA->NN_dist = _R2;
    jetA->NN      = NULL;
}




//----------------------------------------------------------------------
template <class J> inline double ClusterSequence::_bj_dist(
                const J * const jetA, const J * const jetB) const {
  double dphi = std::abs(jetA->phi - jetB->phi);
  double deta = (jetA->eta - jetB->eta);
  if (dphi > pi) {dphi = twopi - dphi;}
  return dphi*dphi + deta*deta;
}

//----------------------------------------------------------------------
template <class J> inline double ClusterSequence::_bj_diJ(const J * const jet) const {
  double kt2 = jet->kt2;
  if (jet->NN != NULL) {if (jet->NN->kt2 < kt2) {kt2 = jet->NN->kt2;}}
  return jet->NN_dist * kt2;
}


//----------------------------------------------------------------------
// set the NN for jet without checking whether in the process you might
// have discovered a new nearest neighbour for another jet
template <class J> inline void ClusterSequence::_bj_set_NN_nocross(
                 J * const jet, J * const head, const J * const tail) const {
  double NN_dist = _R2;
  J * NN  = NULL;
  if (head < jet) {
    for (J * jetB = head; jetB != jet; jetB++) {
      double dist = _bj_dist(jet,jetB);
      if (dist < NN_dist) {
	NN_dist = dist;
	NN = jetB;
      }
    }
  }
  if (tail > jet) {
    for (J * jetB = jet+1; jetB != tail; jetB++) {
      double dist = _bj_dist(jet,jetB);
      if (dist < NN_dist) {
	NN_dist = dist;
	NN = jetB;
      }
    }
  }
  jet->NN = NN;
  jet->NN_dist = NN_dist;
}


//----------------------------------------------------------------------
template <class J> inline void ClusterSequence::_bj_set_NN_crosscheck(J * const jet, 
		    J * const head, const J * const tail) const {
  double NN_dist = _R2;
  J * NN  = NULL;
  for (J * jetB = head; jetB != tail; jetB++) {
    double dist = _bj_dist(jet,jetB);
    if (dist < NN_dist) {
      NN_dist = dist;
      NN = jetB;
    }
    if (dist < jetB->NN_dist) {
      jetB->NN_dist = dist;
      jetB->NN = jet;
    }
  }
  jet->NN = NN;
  jet->NN_dist = NN_dist;
}

FASTJET_END_NAMESPACE

#endif // __FASTJET_CLUSTERSEQUENCE_HH__
