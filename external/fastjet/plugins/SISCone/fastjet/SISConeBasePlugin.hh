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
  }

  /// copy constructor
  SISConeBasePlugin (const SISConeBasePlugin & plugin) {
    *this = plugin;
  }

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

  mutable double _ghost_sep_scale;

  // the part that HAS to be overloaded
  /// call the re-clustering itself 
  virtual void reset_stored_plugin() const =0;

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


FASTJET_END_NAMESPACE        // defined in fastjet/internal/base.hh

#endif // __SISCONEBASEPLUGIN_HH__

