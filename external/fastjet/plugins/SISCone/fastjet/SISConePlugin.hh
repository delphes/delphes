#ifndef __SISCONEPLUGIN_HH__
#define __SISCONEPLUGIN_HH__

#include "SISConeBasePlugin.hh"

// forward declaration of the siscone classes we'll need
namespace siscone{
  class Csiscone;
  class Cjet;
}


FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//----------------------------------------------------------------------
//
/// @ingroup plugins
/// \class SISConePlugin
/// Implementation of the SISCone algorithm (plugin for fastjet v2.1 upwards)
///
/// SISConePlugin is a plugin for fastjet (v2.1 upwards) that provides
/// an interface to the seedless infrared safe cone jet finder by
/// Gregory Soyez and Gavin Salam.
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
///    two jets must overlap (pt_overlap/min(pt1,pt2)) if they are to be 
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
///  - Protojet ptmin: protojets that are below this ptmin
///    (default = 0) are discarded before each iteration of the
///    split-merge loop.
///
/// One parameter governs some internal algorithmic shortcuts: 
///
/// - if "caching" is turned on then the last event clustered by
///   siscone is stored -- if the current event is identical and the
///   cone_radius and n_pass_mass are identical, then the only part of
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
class SISConePlugin : public SISConeBasePlugin{
public:

  /// enum for the different split-merge scale choices;
  /// Note that order _must_ be the same as in siscone
  enum SplitMergeScale {
    SM_pt,     ///< transverse momentum (E-scheme), IR unsafe
    SM_Et,     ///< transverse energy (E-scheme), not long. boost invariant
               ///< original run-II choice [may not be implemented]
    SM_mt,     ///< transverse mass (E-scheme), IR safe except
               ///< in decays of two identical narrow heavy particles
    SM_pttilde ///< pt-scheme pt = \sum_{i in jet} |p_{ti}|, should
               ///< be IR safe in all cases
  };


  /// Main constructor for the SISCone Plugin class.  
  ///
  /// Note: wrt version prior to 2.4 this constructor differs in that a 
  /// the default value has been removed for overlap_threshold. The
  /// former has been removed because the old default of 0.5 was found
  /// to be unsuitable in high-noise environments; so the user should
  /// now explicitly think about the value for this -- we recommend
  /// 0.75.
  ///
  SISConePlugin (double cone_radius_in,
                 double overlap_threshold_in,
                 int    n_pass_max_in = 0,
                 double protojet_ptmin_in = 0.0, 
                 bool   caching_in = false,
                 SplitMergeScale  split_merge_scale_in = SM_pttilde,
                 double split_merge_stopping_scale_in = 0.0){
    _cone_radius           = cone_radius_in;
    _overlap_threshold     = overlap_threshold_in;
    _n_pass_max            = n_pass_max_in;
    _protojet_ptmin        = protojet_ptmin_in;
    _caching               = caching_in;   
    _split_merge_scale     = split_merge_scale_in;
    _split_merge_stopping_scale = split_merge_stopping_scale_in;
    _ghost_sep_scale       = 0.0;
    _use_pt_weighted_splitting = false;
    _user_scale = 0;}


  /// Backwards compatible constructor for the SISCone Plugin class
  SISConePlugin (double cone_radius_in,
                 double overlap_threshold_in,
                 int    n_pass_max_in,
                 double protojet_ptmin_in, 
                 bool   caching_in,
                 bool   split_merge_on_transverse_mass_in){
    _cone_radius           = cone_radius_in;
    _overlap_threshold     = overlap_threshold_in;
    _n_pass_max            = n_pass_max_in;
    _protojet_ptmin        = protojet_ptmin_in;
    _caching               = caching_in;
    _split_merge_stopping_scale = 0.0;
    _split_merge_scale     = split_merge_on_transverse_mass_in ? SM_mt : SM_pttilde;
    _ghost_sep_scale       = 0.0;
    _user_scale = 0;}
  
  /// backwards compatible constructor for the SISCone Plugin class
  /// (avoid using this in future).
  SISConePlugin (double cone_radius_in,
                 double overlap_threshold_in,
                 int    n_pass_max_in,
                 bool   caching_in) {
    _cone_radius           = cone_radius_in;
    _overlap_threshold     = overlap_threshold_in;
    _n_pass_max            = n_pass_max_in;
    _protojet_ptmin        = 0.0;
    _caching               = caching_in;   
    _split_merge_scale     = SM_mt;
    _split_merge_stopping_scale = 0.0;
    _ghost_sep_scale       = 0.0;
    _use_pt_weighted_splitting = false;
    _user_scale = 0;}

  /// minimum pt for a protojet to be considered in the split-merge step
  /// of the algorithm
  double protojet_ptmin  () const {return _protojet_ptmin  ;}

  /// return the scale to be passed to SISCone as the protojet_ptmin
  /// -- if we have a ghost separation scale that is above the
  /// protojet_ptmin, then the ghost_separation_scale becomes the
  /// relevant one to use here
  double protojet_or_ghost_ptmin  () const {return std::max(_protojet_ptmin,
                                                            _ghost_sep_scale);}

  /// indicates scale used in split-merge
  SplitMergeScale split_merge_scale() const {return _split_merge_scale;}
  /// sets scale used in split-merge
  void set_split_merge_scale(SplitMergeScale sms) {_split_merge_scale = sms;}

  /// indicates whether the split-merge orders on transverse mass or not.
  /// retained for backwards compatibility with 2.1.0b3
  bool split_merge_on_transverse_mass() const {return _split_merge_scale == SM_mt ;}
  void set_split_merge_on_transverse_mass(bool val) {
    _split_merge_scale = val  ? SM_mt : SM_pt;}

  /// indicates whether the split-merge orders on transverse mass or not.
  /// retained for backwards compatibility with 2.1.0b3
  bool split_merge_use_pt_weighted_splitting() const {return _use_pt_weighted_splitting;}
  void set_split_merge_use_pt_weighted_splitting(bool val) {
    _use_pt_weighted_splitting = val;}

  // the things that are required by base class
  virtual std::string description () const;
  virtual void run_clustering(ClusterSequence &) const ;

protected:
  virtual void reset_stored_plugin() const;

private:
  double _protojet_ptmin;
  SplitMergeScale _split_merge_scale;

  bool _use_pt_weighted_splitting;

  // part needed for the cache 
  // variables for caching the results and the input
  static std::auto_ptr<SISConePlugin          > stored_plugin;
  static std::auto_ptr<std::vector<PseudoJet> > stored_particles;
  static std::auto_ptr<siscone::Csiscone      > stored_siscone;
};


/////\class SISConePlugin::UserScaleBase::StructureType
///// the structure that allows to store the information contained
///// into a siscone::Cjet (built internally in SISCone from a stable
///// cone) into a PseudoJet
//class SISConePlugin::UserScaleBase::StructureType : public PseudoJetStructureBase {
//public:
//  StructureType(const siscone::Cjet & jet, const ClusterSequence &cs)
//    : _jet(jet), _cs(cs){}
//
//  //--------------------------------------------------
//  // members inherited from the base class
//  /// the textual descripotion
//  virtual std::string description() const;
//
//  /// this structure has constituents
//  virtual bool has_constituents() const {return true;}
//
//  /// retrieve the constituents 
//  ///
//  /// if you simply need to iterate over the constituents, it will be
//  /// faster to access them via constituent(i)
//  virtual std::vector<PseudoJet> constituents(const PseudoJet & /*reference*/) const;
//
//  //--------------------------------------------------
//  // additional information relevant for this structure
//
//  /// returns the number of constituents
//  unsigned int size() const;
//
//  /// returns the index (in the original particle list) of the ith
//  /// constituent
//  int constituent_index(unsigned int i) const;
//
//  /// returns the ith constituent (as a PseusoJet)
//  const PseudoJet & constituent(unsigned int i) const;
//
//  /// returns the scalar pt of this stable cone
//  double pt_tilde() const;
//
//  /// returns the sm_var2 (signed ordering variable squared) for this stable cone
//  double ordering_var2() const;
//
//protected:
//  const siscone::Cjet &_jet;  ///< a dreference to the internal jet in SISCone
//  const ClusterSequence &_cs; ///< a reference to the CS (for access to the particles)
//};

//======================================================================
/// @ingroup extra_info
/// \class SISConeExtras
/// Class that provides extra information about a SISCone clustering
class SISConeExtras : public SISConeBaseExtras {
public:
  /// constructor
  //  it just initialises the pass information 
  SISConeExtras(int nparticles)
    : SISConeBaseExtras(nparticles){}

  /// access to the siscone jet def plugin (more convenient than
  /// getting it from the original jet definition, because here it's
  /// directly of the right type (rather than the base type)
  const SISConePlugin* jet_def_plugin() const {
    return dynamic_cast<const SISConePlugin*>(_jet_def_plugin);
  }

private:
  // let us be written to by SISConePlugin
  friend class SISConePlugin;
};

FASTJET_END_NAMESPACE        // defined in fastjet/internal/base.hh

#endif // __SISCONEPLUGIN_HH__

