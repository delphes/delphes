#ifndef __SISCONEBASEPLUGIN_HH__
#define __SISCONEBASEPLUGIN_HH__

#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include <vector>
#include <memory>
#include <cmath>

#include <sstream>

// questionable whether this should be in fastjet namespace or not...
FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//----------------------------------------------------------------------
//
/// \if internal_doc
/// @ingroup internal
/// \class SISConeBasePlugin
/// Implementation of the SISCone algorithm, base class (plugin for fastjet v2.1 upwards)
///
/// SISConeBasePlugin is a plugin for fastjet (v2.1 upwards) that
/// provides a base interface to SISCone-type cone jet finder by
/// Gregory Soyez and Gavin Salam.
///
/// This is a purely virtual class that needs to be overloaded
/// for the specific implementations of SISCone (i.e. regular or
/// spherical as of July 16th 2008).
///
/// any derived plugin MUST overload the following methods:
///   description()
///   run_siscone_clustering()
///   reset_stored_plugin()
///
/// For further details, see the derived plugins or
/// http://projects.hepforge.com/siscone
///
/// \endif
//
class SISConeBasePlugin : public JetDefinition::Plugin {
public:
  /// default ctor
  SISConeBasePlugin (){
    _use_jet_def_recombiner = false;
    set_progressive_removal(false);
  }

  /// copy constructor
  SISConeBasePlugin (const SISConeBasePlugin & plugin) {
    *this = plugin;
  }

  /// set whether to use SISCone with progressive removal instead of
  /// the default split_merge step.
  ///
  /// If progressive removal is enabled, the following SISCone
  /// variables are not used:
  ///
  /// - overlap_threshold
  /// - caching
  /// - split_merge_stopping_scale
  ///
  /// The split_merge_scale choice is reinterpreted as the ordering
  /// variable for progressive removal. It is also possible for the
  /// user to supply his/her own function for the scale that orders
  /// progressive removal, with set_user_scale(...)
  void set_progressive_removal(bool progressive_removal_in=true){
    _progressive_removal = progressive_removal_in;
  }

  /// returns true if progressive_removal is enabled
  bool progressive_removal() const{ return _progressive_removal;}

  /// the cone radius
  double cone_radius        () const {return _cone_radius        ;}

  /// Fraction of overlap energy in a jet above which jets are merged
  /// and below which jets are split.
  double overlap_threshold  () const {return _overlap_threshold  ;}

  /// the maximum number of passes of stable-cone searching (<=0 is same
  /// as infinity).
  int n_pass_max  () const {return _n_pass_max  ;}

  /// set the "split_merge_stopping_scale": if the scale variable for
  /// all protojets is below this, then stop the split-merge procedure
  /// and keep only those jets found so far. This is useful in
  /// determination of areas of hard jets because it can be used to
  /// avoid running the split-merging on the pure ghost-part of the
  /// event.
  void set_split_merge_stopping_scale(double scale) {
    _split_merge_stopping_scale = scale;}

  /// return the value of the split_merge_stopping_scale (see
  /// set_split_merge_stopping_scale(...) for description)
  double split_merge_stopping_scale() {return _split_merge_stopping_scale;}

  /// allow the user to decide if one uses the jet_def's own recombination scheme
  void set_use_jet_def_recombiner(bool choice) {_use_jet_def_recombiner = choice;}

  /// indicate if the jet_def's recombination scheme is being used
  bool use_jet_def_recombiner() const {return _use_jet_def_recombiner;}

  /// indicates whether caching is turned on or not.
  bool caching() const {return _caching ;}

  /// the plugin mechanism's standard way of accessing the jet radius
  virtual double R() const {return cone_radius();}

  /// return true since there is specific support for the measurement
  /// of passive areas, in the sense that areas determined from all
  /// particles below the ghost separation scale will be a passive
  /// area. 
  virtual bool supports_ghosted_passive_areas() const {
    return true;
  }
  
  /// set the ghost separation scale for passive area determinations
  /// _just_ in the next run (strictly speaking that makes the routine
  /// a non const, so related internal info must be stored as a mutable)
  virtual void set_ghost_separation_scale(double scale) const {
    _ghost_sep_scale = scale;
  }

  virtual double ghost_separation_scale() const {
    return _ghost_sep_scale;
  }

  // user-defined scale for progressive removal
  //------------------------------------------------------------

  /// \class UserScaleBase
  /// base class for user-defined ordering of stable cones (used for
  /// prorgessive removal)
  ///
  /// derived classes have to implement the () operator that returns
  /// the scale associated with a given jet.
  /// 
  /// It is also highly recommended to implement the is_larger()
  /// method whenever possible, in order to avoid rounding issues
  /// known to lead to possible infrared unsafeties.
  ///
  /// The jets that are passed to this class will carry the structure
  /// of type SISConePlugin::StructureType which allows to retreive
  /// easily the following information:
  ///
  ///   vector<PseudoJet> constituents = jet.constituents();
  ///   unsigned int n_constituents = jet.structure_of<SISConePlugin::UserScaleBase>().size();
  ///   int index = jet.structure_of<SISConePlugin::UserScaleBase>().constituent_index(index i);
  ///   const PseudoJet & p = jet.structure_of<SISConePlugin::UserScaleBase>().constituent(index i);
  ///   double scalar_pt = jet.structure_of<SISConePlugin::UserScaleBase>().pt_tilde();
  ///
  /// see SISConePlugin::StructureType below for further details
  class UserScaleBase : public FunctionOfPseudoJet<double>{
  public:
    /// empty virtual dtor
    virtual ~UserScaleBase(){}

    /// returns the scale associated with a given jet
    ///
    /// "progressive removal" iteratively removes the stable cone with
    /// the largest scale
    virtual double result(const PseudoJet & jet) const = 0;

    /// returns true when the scale associated with jet a is larger than
    /// the scale associated with jet b
    ///
    /// By default this does a simple direct comparison but it can be
    /// overloaded for higher precision [recommended if possible]
    virtual bool is_larger(const PseudoJet & a, const PseudoJet & b) const;

    class StructureType; // defined below
  };

  // template class derived from UserScaleBase::StryctureType that
  // works for both SISCone jet classes
  // implemented below 
  template<class Tjet> 
  class UserScaleBaseStructureType;

  /// set a user-defined scale for stable-cone ordering in
  /// progressive removal
  void set_user_scale(const UserScaleBase *user_scale_in){ _user_scale = user_scale_in;}

  /// returns the user-defined scale in use (0 if none)
  const UserScaleBase * user_scale() const{ return _user_scale;}


  // the things that one MUST overload required by base class
  //---------------------------------------------------------

  /// plugin description
  virtual std::string description () const =0;

  /// really do the clustering work
  virtual void run_clustering(ClusterSequence &) const = 0;

protected:
  double _cone_radius, _overlap_threshold;
  int    _n_pass_max;
  bool   _caching;//, _split_merge_on_transverse_mass;
  double _split_merge_stopping_scale;
  bool   _use_jet_def_recombiner;
  bool   _progressive_removal;

  mutable double _ghost_sep_scale;

  // the part that HAS to be overloaded
  /// call the re-clustering itself 
  virtual void reset_stored_plugin() const =0;

  const UserScaleBase * _user_scale;

};


//======================================================================
/// @ingroup extra_info
/// \class SISConeBaseExtras
/// Class that provides extra information about a SISCone clustering
///
/// This is only the base class that the "regular" and "spherical"
/// implementations of SISCone will have to overload. The only thing
/// that needs to be done for the derived classes is to define
/// '_jet_def_plugin', implement
///   jet_def_plugin();
/// and add the corresponding plugin class as a friend
class SISConeBaseExtras : public ClusterSequence::Extras {
public:

  /// constructor
  //  it just initialises the pass information 
  SISConeBaseExtras(int nparticles) : _pass(nparticles*2,-1) {}

  /// purely virtual destructor
  inline virtual ~SISConeBaseExtras()=0;

  /// returns a reference to the vector of stable cones (aka protocones)
  const std::vector<PseudoJet> & stable_cones() const {return _protocones;}

  /// an old name for getting the vector of stable cones (aka protocones)
  const std::vector<PseudoJet> & protocones() const {return _protocones;}

  /// return the # of the pass at which a given jet was found; will
  /// return -1 if the pass is invalid
  int pass(const PseudoJet & jet) const {return _pass[jet.cluster_hist_index()];}

  /// return a brief summary of the contents of the extras object
  /// (specifically, the number of protocones.
  std::string description() const{
    std::ostringstream ostr;
    ostr << "This SISCone clustering found " << protocones().size()
	 << " stable protocones";
    return ostr.str();
  };

  /// return the smallest difference in squared distance encountered
  /// during splitting between a particle and two overlapping
  /// protojets.
  inline double most_ambiguous_split() const {return _most_ambiguous_split;}

protected:
  std::vector<PseudoJet> _protocones;
  std::vector<int>       _pass;
  double                _most_ambiguous_split;
  const SISConeBasePlugin * _jet_def_plugin;
};

/// give the destructor its required implementation
inline SISConeBaseExtras::~SISConeBaseExtras(){}

//----------------------------------------------------------------------
// implementation of the structure type associated with the UserScaleBase class

/// \class SISConeBasePlugin::UserScaleBase::StructureType
/// the structure that allows to store the information contained
/// into a siscone::Cjet (built internally in SISCone from a stable
/// cone) into a PseudoJet
class SISConeBasePlugin::UserScaleBase::StructureType : public PseudoJetStructureBase {
public:
  /// base ctor (constructed from a ClusterSequence tin order to have
  /// access to the initial particles
  StructureType(const ClusterSequence &cs)
    : _cs(cs){}

  /// empty virtual dtor 
  virtual ~StructureType(){}
  
  //--------------------------------------------------
  // members inherited from the base class
  /// the textual descripotion
  virtual std::string description() const{
    return "PseudoJet wrapping a siscone jet from a stable cone"; 
  }

  /// this structure has constituents
  virtual bool has_constituents() const {return true;}

  /// retrieve the constituents 
  ///
  /// if you simply need to iterate over the constituents, it will be
  /// faster to access them via constituent(i)
  virtual std::vector<PseudoJet> constituents(const PseudoJet & /*reference*/) const{ 
    std::vector<PseudoJet> constits;
    constits.reserve(size());
    for (unsigned int i=0; i<size();i++)
      constits.push_back(constituent(i));
    return constits;
  }
  
  //--------------------------------------------------
  // additional information relevant for this structure

  /// returns the number of constituents
  virtual unsigned int size() const = 0;

  /// returns the index (in the original particle list) of the ith
  /// constituent
  virtual int constituent_index(unsigned int i) const = 0;

  /// returns the ith constituent (as a PseusoJet)
  const PseudoJet & constituent(unsigned int i) const{
    return _cs.jets()[constituent_index(i)];
  }

  // /// returns the scalar pt of this stable cone
  // virtual double pt_tilde() const = 0;

  /// returns the sm_var2 (signed ordering variable squared) for this stable cone
  virtual double ordering_var2() const = 0;

protected:
  const ClusterSequence &_cs; ///< a reference to the CS (for access to the particles)
};


///@ingroup internal
/// template class derived from UserScaleBase::StryctureType that
/// works for both SISCone jet classes
/// implemented below 
template<class Tjet>
class SISConeBasePlugin::UserScaleBaseStructureType : public UserScaleBase::StructureType{
public:
  UserScaleBaseStructureType(const Tjet &jet, const ClusterSequence &cs)
    : UserScaleBase::StructureType(cs), _jet(jet){}

  /// empty virtual dtor 
  virtual ~UserScaleBaseStructureType(){}

  //--------------------------------------------------
  // additional information relevant for this structure

  /// returns the number of constituents
  virtual unsigned int size() const{
    return _jet.n;
  }

  /// returns the index (in the original particle list) of the ith
  /// constituent
  virtual int constituent_index(unsigned int i) const{
    return _jet.contents[i];
  }

  // /// returns the scalar pt of this stable cone
  // virtual double pt_tilde() const{
  //   return _jet.pt_tilde;
  // }

  /// returns the sm_var2 (signed ordering variable squared) for this stable cone
  virtual double ordering_var2() const{
    return _jet.sm_var2;
  }

protected:
  const Tjet &_jet; ///< a reference to the internal jet in SISCone
};


FASTJET_END_NAMESPACE        // defined in fastjet/internal/base.hh

#endif // __SISCONEBASEPLUGIN_HH__

