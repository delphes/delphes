#ifndef __SISCONESPHERICALPLUGIN_HH__
#define __SISCONESPHERICALPLUGIN_HH__

#include "SISConeBasePlugin.hh"

// forward declaration of the siscone classes we'll need
namespace siscone_spherical{
  class CSphsiscone;
}

// questionable whether this should be in fastjet namespace or not...
FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//----------------------------------------------------------------------
//
/// @ingroup plugins
/// \class SISConeSphericalPlugin
/// Implementation of the spherical version of the SISCone algorithm
/// (plugin for fastjet v2.1 upwards)
///
/// SISConeSphericalPlugin is a plugin for fastjet (v2.1 upwards) that
/// provides an interface to the seedless infrared safe cone jet
/// finder by Gregory Soyez and Gavin Salam.
///
/// This is the version of SISCone using spherical coordinates. Compared
/// to the original cylindrical version:
///
///  - Particles are within a cone if their opening angle relative to the
///    centre of the cone is less than R
///
///  - The split-merge step uses the total energy in the protojet as the
///    ordering and overlap-measure variable
///
///  - The IR safety of the split-merge step is _not_ guaranteed for
///    events consisting of two back-to-back identical heavy particles
///    that decay. This is because of potential degeneracies in the
///    ordering for the split-merge step. 
///
///    For moderate values of R the problem should not be too severe
///    (or may even be absent for some values of the overlap
///    parameter), however the user should be aware of the issue.
///
///    The default split-merge scale may change at a later date to
///    resolve this issue.
///
///
/// SISCone uses geometrical techniques to exhaustively consider all
/// possible distinct cones. It then finds out which ones are stable
/// and sends the result to the Tevatron Run-II type split-merge
/// procedure for overlapping cones.
///
/// Four parameters govern the "physics" of the algorithm:
///
///  - the cone_radius (this should be self-explanatory!)
///
///  - the overlap_threshold is the parameter which dictates how much
///    two jets must overlap (E_overlap/min(E1,E2)) if they are to be 
///    merged
///
///  - Not all particles are in stable cones in the first round of
///    searching for stable cones; one can therefore optionally have the
///    the jet finder carry out additional passes of searching for
///    stable cones among particles that were in no stable cone in
///    previous passes --- the maximum number of passes carried out is
///    n_pass_max. If this is zero then additional passes are carried
///    out until no new stable cones are found.
///
///  - Protojet Emin: protojets that are below this Emin
///    (default = 0) are discarded before each iteration of the
///    split-merge loop.
///
/// One parameter governs some internal algorithmic shortcuts: 
///
/// - if "caching" is turned on then the last event clustered by
///   siscone is stored -- if the current event is identical and the
///   cone_radius and n_pass_max are identical, then the only part of
///   the clustering that needs to be rerun is the split-merge part,
///   leading to significant speed gains; there is a small (O(N) storage
///   and speed) penalty for caching, so it should be kept off
///   (default) if only a single overlap_threshold is used.
///
/// The final jets can be accessed by requestion the
/// inclusive_jets(...) from the ClusterSequence object. Note that
/// these PseudoJets have their user_index() set to the index of the
/// pass in which they were found (first pass = 0). NB: This does not
/// currently work for jets that consist of a single particle.
///
/// For further information on the details of the algorithm see the
/// SISCone paper, arXiv:0704.0292 [JHEP 0705:086,2007].
///
/// For documentation about the implementation, see the
/// siscone/doc/html/index.html file.
//
class SISConeSphericalPlugin : public SISConeBasePlugin{
public:

  /// enum for the different split-merge scale choices;
  /// Note that order _must_ be the same as in siscone
  enum SplitMergeScale {SM_E,        ///< Energy (IR unsafe with momentum conservation)
			SM_Etilde    ///< sum_{i \in jet} E_i [1+sin^2(theta_iJ)]
  };


  /// Main constructor for the SISConeSpherical Plugin class.  
  ///
  ///
  SISConeSphericalPlugin (double cone_radius_in,
			  double overlap_threshold_in,
			  int    n_pass_max_in = 0,
			  double protojet_Emin_in = 0.0, 
			  bool   caching_in = false,
			  SplitMergeScale  split_merge_scale_in = SM_Etilde,
			  double split_merge_stopping_scale_in = 0.0){
    _cone_radius           =cone_radius_in;
    _overlap_threshold     =overlap_threshold_in;
    _n_pass_max            =n_pass_max_in;
    _protojet_Emin         =protojet_Emin_in;
    _caching               =caching_in;
    _split_merge_scale     =split_merge_scale_in;
    _split_merge_stopping_scale = split_merge_stopping_scale_in;
    _ghost_sep_scale       = 0.0;
    _use_E_weighted_splitting = false;
  }

  /// minimum energy for a protojet to be considered in the split-merge step
  /// of the algorithm
  double protojet_Emin  () const {return _protojet_Emin  ;}

  /// return the scale to be passed to SISCone as the protojet_Emin
  /// -- if we have a ghost separation scale that is above the
  /// protojet_ptmin, then the ghost_separation_scale becomes the
  /// relevant one to use here
  double protojet_or_ghost_Emin  () const {return std::max(_protojet_Emin,
                                                           _ghost_sep_scale);}

  /// indicates scale used in split-merge
  SplitMergeScale split_merge_scale() const {return _split_merge_scale;}
  /// sets scale used in split-merge
  void set_split_merge_scale(SplitMergeScale sms) {_split_merge_scale = sms;}

  /// indicate if the splittings are done using the anti-kt distance
  bool split_merge_use_E_weighted_splitting() const {return _use_E_weighted_splitting;}
  void set_split_merge_use_E_weighted_splitting(bool val) {
    _use_E_weighted_splitting = val;}

  /// overload the default as we don't provide support 
  /// for passive areas.
  virtual bool supports_ghosted_passive_areas() const {return true;}
  
  // the things that are required by base class
  virtual std::string description () const;
  virtual void run_clustering(ClusterSequence &) const ;

  /// returns true because this plugin is intended for spherical
  /// geometries (i.e. it's an e+e- algorithm).
  virtual bool is_spherical() const {return true;}

protected:
  virtual void reset_stored_plugin() const;

private:
  double _protojet_Emin;
  SplitMergeScale _split_merge_scale;
  bool _use_E_weighted_splitting;

  // part needed for the cache 
  // variables for caching the results and the input
  static std::auto_ptr<SISConeSphericalPlugin        > stored_plugin;
  static std::auto_ptr<std::vector<PseudoJet>        > stored_particles;
  static std::auto_ptr<siscone_spherical::CSphsiscone> stored_siscone;
};

//======================================================================
/// @ingroup extra_info
/// \class SISConeSphericalExtras
/// Class that provides extra information about a SISCone clustering
class SISConeSphericalExtras : public SISConeBaseExtras {
public:
  /// constructor
  //  it just initialises the pass information 
  SISConeSphericalExtras(int nparticles)
    : SISConeBaseExtras(nparticles){}

  /// access to the siscone jet def plugin (more convenient than
  /// getting it from the original jet definition, because here it's
  /// directly of the right type (rather than the base type)
  const SISConeSphericalPlugin* jet_def_plugin() const {
    return dynamic_cast<const SISConeSphericalPlugin*>(_jet_def_plugin);
  }

private:
  // let us be written to by SISConePlugin
  friend class SISConeSphericalPlugin;
};

FASTJET_END_NAMESPACE        // defined in fastjet/internal/base.hh

#endif // __SISCONEPLUGIN_HH__

